"""Sentinel-2 SAFE archive I/O: read JP2 bands, extract metadata, resample to 10m."""

from pathlib import Path
from xml.etree import ElementTree

import numpy as np


# Band ordering matching the 13-band Rust layout
# (B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12)
SENTINEL_BAND_ORDER = [
    "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08",
    "B8A", "B09", "B10", "B11", "B12",
]

# Native resolutions (meters)
BAND_RESOLUTION = {
    "B01": 60, "B02": 10, "B03": 10, "B04": 10,
    "B05": 20, "B06": 20, "B07": 20, "B08": 10,
    "B8A": 20, "B09": 60, "B10": 60, "B11": 20, "B12": 20,
}


def read_sentinel_safe(safe_dir: str | Path) -> dict:
    """Read a Sentinel-2 SAFE archive and return TOA bands + metadata.

    Parameters
    ----------
    safe_dir : path
        Path to the .SAFE directory.

    Returns
    -------
    dict with keys:
        toa_bands : list of 13 numpy arrays (float32, 10980x10980)
        angles : dict with solar_zenith, solar_azimuth, view_zenith, view_azimuth (scalars)
        profile : dict with epsg, ul_x, ul_y, pixel_size, nlines, nsamps
    """
    import rasterio

    safe_dir = Path(safe_dir)

    # Parse metadata
    mtd_msil1c = safe_dir / "MTD_MSIL1C.xml"
    mtd_tl = safe_dir / "MTD_TL.xml"
    angles = _parse_angles(mtd_tl if mtd_tl.exists() else mtd_msil1c)
    quantification_value, radiometric_offset = _parse_quantification(mtd_msil1c)

    # Read bands
    toa_bands = []
    nlines_10m = None
    nsamps_10m = None
    profile = {}

    # Track fill mask: DN == 0 for any band marks fill (bit 0 of QA)
    fill_mask = None

    for band_name in SENTINEL_BAND_ORDER:
        # Find JP2 file
        jp2_files = list(safe_dir.glob(f"**/*_{band_name}.jp2"))
        if not jp2_files:
            raise FileNotFoundError(f"No JP2 file found for {band_name} in {safe_dir}")

        with rasterio.open(jp2_files[0]) as src:
            dn = src.read(1).astype(np.float32)

            # Get 10m reference dimensions from B02
            if band_name == "B02":
                nlines_10m = src.height
                nsamps_10m = src.width
                transform = src.transform
                crs = src.crs
                profile = {
                    "epsg": crs.to_epsg(),
                    "ul_x": transform.c,
                    "ul_y": transform.f,
                    "pixel_size_x": transform.a,
                    "pixel_size_y": abs(transform.e),
                    "nlines": nlines_10m,
                    "nsamps": nsamps_10m,
                }

        # Detect fill at native resolution before resampling
        dn_is_zero = dn == 0.0

        # Resample fill mask to 10m if needed
        native_res = BAND_RESOLUTION[band_name]
        if native_res != 10:
            dn_is_zero = _resample_to_10m(
                dn_is_zero.astype(np.float32), native_res, nlines_10m, nsamps_10m
            ) > 0.5

        # Accumulate fill mask (any band with DN==0 → fill)
        if fill_mask is None:
            fill_mask = dn_is_zero
        else:
            fill_mask |= dn_is_zero

        # Unscale to TOA reflectance
        # Since Baseline 4.00: toa = (DN + offset) / quantification_value
        toa = (dn + radiometric_offset) / quantification_value

        # Resample to 10m if needed
        if native_res != 10:
            toa = _resample_to_10m(toa, native_res, nlines_10m, nsamps_10m)

        toa_bands.append(toa)

    # Build QA band: bit 0 = fill (matches C level1_qa_is_fill)
    qa_band = np.zeros((nlines_10m, nsamps_10m), dtype=np.uint16)
    qa_band[fill_mask] = 1

    return {
        "toa_bands": toa_bands,
        "qa_band": qa_band,
        "angles": angles,
        "profile": profile,
    }


def _resample_to_10m(
    data: np.ndarray, native_res: int, nlines_10m: int, nsamps_10m: int
) -> np.ndarray:
    """Nearest-neighbor resample from 20m or 60m to 10m.

    Matches C code's convert_to_10m: each native pixel maps to a
    (native_res/10) x (native_res/10) block of 10m pixels.
    """
    scale = native_res // 10
    # np.repeat along both axes = nearest-neighbor upsampling
    resampled = np.repeat(np.repeat(data, scale, axis=0), scale, axis=1)
    # Trim to exact 10m dimensions (handles edge cases)
    return resampled[:nlines_10m, :nsamps_10m]


def _parse_angles(xml_path: Path) -> dict:
    """Extract scene-center mean angles from MTD_TL.xml or MTD_MSIL1C.xml."""
    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    # Remove namespace prefix for easier searching
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"

    angles = {}

    # Sun angles
    sun_el = root.find(f".//{ns}Mean_Sun_Angle/{ns}ZENITH_ANGLE")
    if sun_el is None:
        sun_el = root.find(".//Mean_Sun_Angle/ZENITH_ANGLE")
    if sun_el is not None:
        angles["solar_zenith"] = float(sun_el.text)
    else:
        angles["solar_zenith"] = 30.0  # fallback

    sun_az = root.find(f".//{ns}Mean_Sun_Angle/{ns}AZIMUTH_ANGLE")
    if sun_az is None:
        sun_az = root.find(".//Mean_Sun_Angle/AZIMUTH_ANGLE")
    if sun_az is not None:
        angles["solar_azimuth"] = float(sun_az.text)
    else:
        angles["solar_azimuth"] = 150.0

    # View angles — average across all bands
    view_zen_els = root.findall(f".//{ns}Mean_Viewing_Incidence_Angle/{ns}ZENITH_ANGLE")
    if not view_zen_els:
        view_zen_els = root.findall(".//Mean_Viewing_Incidence_Angle/ZENITH_ANGLE")
    if view_zen_els:
        angles["view_zenith"] = np.mean([float(e.text) for e in view_zen_els])
    else:
        angles["view_zenith"] = 0.0

    view_az_els = root.findall(f".//{ns}Mean_Viewing_Incidence_Angle/{ns}AZIMUTH_ANGLE")
    if not view_az_els:
        view_az_els = root.findall(".//Mean_Viewing_Incidence_Angle/AZIMUTH_ANGLE")
    if view_az_els:
        angles["view_azimuth"] = np.mean([float(e.text) for e in view_az_els])
    else:
        angles["view_azimuth"] = 0.0

    return angles


def _parse_quantification(xml_path: Path) -> tuple[float, float]:
    """Extract quantification value and radiometric offset from MTD_MSIL1C.xml.

    Returns (quantification_value, radiometric_offset).
    """
    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    # Quantification value
    qv_el = root.find(".//QUANTIFICATION_VALUE")
    if qv_el is None:
        qv_el = root.find(".//{*}QUANTIFICATION_VALUE")
    quantification_value = float(qv_el.text) if qv_el is not None else 10000.0

    # Radiometric offset (Baseline 4.00+)
    offset_el = root.find(".//RADIO_ADD_OFFSET")
    if offset_el is None:
        offset_el = root.find(".//{*}RADIO_ADD_OFFSET")
    radiometric_offset = float(offset_el.text) if offset_el is not None else -1000.0

    return quantification_value, radiometric_offset
