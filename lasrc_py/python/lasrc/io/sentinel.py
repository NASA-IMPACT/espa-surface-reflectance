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
    angles = _parse_angles(_find_tile_metadata(safe_dir))
    quantification_value, radiometric_offsets = _parse_quantification(mtd_msil1c)

    # Read bands
    toa_bands = []
    nlines_10m = None
    nsamps_10m = None
    profile = {}

    # Track fill mask: DN == 0 for any band marks fill (bit 0 of QA)
    fill_mask = None

    for band_name, radiometric_offset in zip(SENTINEL_BAND_ORDER, radiometric_offsets):
        with rasterio.open(_find_band_file(safe_dir, band_name)) as src:
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

        toa = _dn_to_toa(dn, radiometric_offset, quantification_value)

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


def _dn_to_toa(dn: np.ndarray, offset: float, quantification_value: float) -> np.ndarray:
    """Convert L1C DN to TOA reflectance with C LaSRC's float arithmetic.

    C computes (DN + add_offset) * scale_factor in float, with scale_factor =
    1/QUANTIFICATION_VALUE stored as a float. Multiplying by that scale is not
    bit-identical to dividing by QUANTIFICATION_VALUE.
    """
    scale = np.float32(1.0 / quantification_value)
    return (dn.astype(np.float32, copy=False) + np.float32(offset)) * scale


def _find_band_file(safe_dir: Path, band_name: str) -> Path:
    """Locate the single IMG_DATA JP2 for ``band_name`` inside a SAFE archive.

    QI_DATA holds masks with the same ``*_<band>.jp2`` suffix (e.g.
    MSK_QUALIT_B01.jp2), so the search is restricted to IMG_DATA.
    """
    pattern = f"GRANULE/*/IMG_DATA/*_{band_name}.jp2"
    matches = sorted(Path(safe_dir).glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No {pattern} found in {safe_dir}")
    if len(matches) > 1:
        raise ValueError(f"Expected one {pattern} in {safe_dir}, found {len(matches)}")
    return matches[0]


def _find_tile_metadata(safe_dir: Path) -> Path:
    """Locate the single granule-level MTD_TL.xml inside a SAFE archive."""
    matches = sorted(Path(safe_dir).glob("GRANULE/*/MTD_TL.xml"))
    if not matches:
        raise FileNotFoundError(f"No GRANULE/*/MTD_TL.xml found in {safe_dir}")
    if len(matches) > 1:
        raise ValueError(
            f"Expected one GRANULE/*/MTD_TL.xml in {safe_dir}, found {len(matches)}"
        )
    return matches[0]


def _espa_angle(text: str) -> float:
    """Round an angle the way C LaSRC receives it through the ESPA XML.

    espa-product-formatter parses the value into a float, writes it with
    "%f", and LaSRC reads it back into a float. Each parse is atof (to
    double) then assignment to float, hence float() before np.float32().
    """
    return float(np.float32(float(f"{float(np.float32(float(text))):f}")))


def _parse_angles(mtd_tl: Path) -> dict:
    """Extract scene-center sun and view angles from a granule MTD_TL.xml.

    Matches the C LaSRC input chain: the view angle is the first
    Mean_Viewing_Incidence_Angle entry (espa-product-formatter ignores the
    rest), and every angle carries ESPA's float/"%f" rounding.
    """
    root = ElementTree.parse(mtd_tl).getroot()

    def _find_all(path: str) -> list[float]:
        return [_espa_angle(e.text) for e in root.iterfind(path)]

    sun_zen = _find_all(".//{*}Mean_Sun_Angle/{*}ZENITH_ANGLE")
    sun_az = _find_all(".//{*}Mean_Sun_Angle/{*}AZIMUTH_ANGLE")
    if len(sun_zen) != 1 or len(sun_az) != 1:
        raise ValueError(
            f"Expected one Mean_Sun_Angle ZENITH_ANGLE and AZIMUTH_ANGLE in {mtd_tl}"
        )

    view_zen = _find_all(".//{*}Mean_Viewing_Incidence_Angle/{*}ZENITH_ANGLE")
    view_az = _find_all(".//{*}Mean_Viewing_Incidence_Angle/{*}AZIMUTH_ANGLE")
    if not view_zen or len(view_zen) != len(view_az):
        raise ValueError(
            f"Missing or incomplete Mean_Viewing_Incidence_Angle entries in {mtd_tl}"
        )

    return {
        "solar_zenith": sun_zen[0],
        "solar_azimuth": sun_az[0],
        "view_zenith": view_zen[0],
        "view_azimuth": view_az[0],
    }


def _parse_quantification(mtd_msil1c: Path) -> tuple[float, list[float]]:
    """Extract the quantification value and per-band offsets from MTD_MSIL1C.xml.

    Returns (quantification_value, offsets) with offsets ordered as
    SENTINEL_BAND_ORDER, which matches the RADIO_ADD_OFFSET band_id order.
    """
    root = ElementTree.parse(mtd_msil1c).getroot()

    quant = [float(e.text) for e in root.iterfind(".//{*}QUANTIFICATION_VALUE")]
    if len(quant) != 1:
        raise ValueError(
            f"Expected one QUANTIFICATION_VALUE in {mtd_msil1c}, found {len(quant)}"
        )

    offsets: dict[int, float] = {}
    for e in root.iterfind(".//{*}RADIO_ADD_OFFSET"):
        band_id = int(e.get("band_id"))
        if band_id in offsets or not 0 <= band_id < len(SENTINEL_BAND_ORDER):
            raise ValueError(
                f"Duplicate or out-of-range RADIO_ADD_OFFSET band_id={band_id} in {mtd_msil1c}"
            )
        offsets[band_id] = float(e.text)
    if len(offsets) != len(SENTINEL_BAND_ORDER):
        raise ValueError(
            f"Expected RADIO_ADD_OFFSET for {len(SENTINEL_BAND_ORDER)} bands in "
            f"{mtd_msil1c}, found {len(offsets)}"
        )

    return quant[0], [offsets[i] for i in range(len(SENTINEL_BAND_ORDER))]
