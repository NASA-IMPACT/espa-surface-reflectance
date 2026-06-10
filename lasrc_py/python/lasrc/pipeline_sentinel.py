"""Sentinel-2 processing pipeline: SAFE archive → surface reflectance."""

from pathlib import Path

import numpy as np

from lasrc.aux import load_auxiliary_data, load_dem, load_ratio_file
from lasrc.io_sentinel import read_sentinel_safe
from lasrc.sensors import SENSORS

import lasrc as _lasrc


def process_sentinel_scene(
    safe_dir: str | Path,
    angle_hdf: str | Path,
    intref_hdf: str | Path,
    transm_hdf: str | Path,
    sphera_hdf: str | Path,
    viirs_aux_path: str | Path,
    dem_path: str | Path,
    ratio_path: str | Path,
    output_path: str | Path,
    sensor_name: str = "SENTINEL_2B",
    output_format: str = "espa",
) -> None:
    """Process a Sentinel-2 SAFE archive to surface reflectance."""
    sensor_config = SENSORS[sensor_name]
    nsr_bands = sensor_config["nsr_bands"]  # 13

    # Read SAFE archive
    scene = read_sentinel_safe(safe_dir)
    angles = scene["angles"]
    profile = scene["profile"]
    nlines = profile["nlines"]
    nsamps = profile["nsamps"]

    # Load LUTs (HDF4)
    lut_data = _load_sentinel_luts(
        str(angle_hdf), str(intref_hdf), str(transm_hdf), str(sphera_hdf), nsr_bands
    )
    lut = _lasrc.PyLookupTables(**lut_data)

    # Load auxiliary data
    aux_data = load_auxiliary_data(str(viirs_aux_path), aux_source="VIIRS")
    dem = load_dem(str(dem_path))
    ratios = load_ratio_file(str(ratio_path))
    aux = _lasrc.PyAuxiliaryData(dem=dem, **aux_data, **ratios)

    # Extract UTM zone from EPSG
    epsg = profile["epsg"]
    utm_zone = epsg % 100
    if epsg > 32700:
        utm_zone = -utm_zone

    # Run Rust correction
    result = _lasrc.process_sentinel_surface_reflectance(
        sensor_name=sensor_name,
        toa_bands=scene["toa_bands"],
        solar_zenith=angles["solar_zenith"],
        solar_azimuth=angles["solar_azimuth"],
        view_zenith=angles["view_zenith"],
        view_azimuth=angles["view_azimuth"],
        qa_band=scene["qa_band"],
        lut=lut,
        aux=aux,
        ul_corner_x=profile["ul_x"],
        ul_corner_y=profile["ul_y"],
        pixel_size_x=profile["pixel_size_x"],
        pixel_size_y=profile["pixel_size_y"],
        utm_zone=utm_zone,
    )

    # Reshape and write output
    for i in range(len(result["sr_bands"])):
        result["sr_bands"][i] = result["sr_bands"][i].reshape(nlines, nsamps)
    result["aerosol"] = result["aerosol"].reshape(nlines, nsamps)
    result["qa"] = result["qa"].reshape(nlines, nsamps)

    if output_format == "espa":
        _write_espa_sentinel(output_path, result, sensor_config)


def _write_espa_sentinel(output_dir, result, sensor_config):
    """Write Sentinel SR output in ESPA format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, name in enumerate(sensor_config["refl_band_names"]):
        result["sr_bands"][i].astype(np.uint16).tofile(output_dir / f"{name}.img")

    result["aerosol"].astype(np.int16).tofile(output_dir / "sr_aerosol.img")
    result["qa"].astype(np.uint8).tofile(output_dir / "sr_aerosol_qa.img")


def _load_sentinel_luts(
    angle_hdf_path: str,
    intref_hdf_path: str,
    transm_hdf_path: str,
    sphera_hdf_path: str,
    nsr_bands: int,
) -> dict:
    """Load Sentinel-2 LUTs from HDF4 files.

    The MSI LUTs use band names (NRLUT_BAND_1, NRLUT_BAND_8a, etc.)
    instead of sequential numbers.
    """
    from pyhdf.SD import SD, SDC

    NSOLAR = 8000
    NAOT = 22
    NPRES = 7

    # Band name mapping for RES_LUT dataset names
    band_ds_names = [
        "NRLUT_BAND_1", "NRLUT_BAND_2", "NRLUT_BAND_3", "NRLUT_BAND_4",
        "NRLUT_BAND_5", "NRLUT_BAND_6", "NRLUT_BAND_7", "NRLUT_BAND_8",
        "NRLUT_BAND_8a", "NRLUT_BAND_9", "NRLUT_BAND_10", "NRLUT_BAND_11",
        "NRLUT_BAND_12",
    ]

    # Angle LUT
    hdf = SD(angle_hdf_path, SDC.READ)
    angle_data = {}
    for name in ["TSMAX", "TSMIN", "NBFIC", "TTV", "TTS"]:
        angle_data[name.lower()] = hdf.select(name).get().astype(np.float64).ravel().tolist()
    angle_data["nbfi"] = hdf.select("NBFI").get().astype(np.int32).ravel().tolist()
    angle_data["indts"] = hdf.select("INDTS").get().astype(np.int32).ravel().tolist()
    hdf.end()

    # RES LUT (rolutt)
    hdf = SD(intref_hdf_path, SDC.READ)
    rolutt = np.zeros(nsr_bands * NPRES * NAOT * NSOLAR, dtype=np.float64)
    for ib, ds_name in enumerate(band_ds_names):
        data = hdf.select(ds_name).get().astype(np.float64)  # [NSOLAR, NAOT, NPRES]
        for ip in range(NPRES):
            for ia in range(NAOT):
                base = ib * NPRES * NAOT * NSOLAR + ip * NAOT * NSOLAR + ia * NSOLAR
                rolutt[base:base + NSOLAR] = data[:, ia, ip]
    hdf.end()

    # AERO LUT (sphalbt, normext) — ASCII format
    sphalbt, normext = _load_aero_ascii(sphera_hdf_path, nsr_bands)

    # TRANS LUT — ASCII format
    transt = _load_trans_ascii(transm_hdf_path, nsr_bands)

    return {
        "rolutt": rolutt.tolist(),
        "transt": transt,
        "sphalbt": sphalbt,
        "normext": normext,
        **angle_data,
        "nsr_bands": nsr_bands,
    }


def _load_aero_ascii(path: Path, nsr_bands: int) -> tuple[list, list]:
    """Load sphalbt and normext from ASCII file. Same format as Landsat."""
    NPRES = 7
    NAOT = 22

    sphalbt = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)
    normext = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        idx += 1  # band header
        for ip in range(NPRES):
            idx += 1  # pressure header
            for ia in range(NAOT):
                parts = lines[idx].split()
                base = ib * NPRES * NAOT + ip * NAOT + ia
                sphalbt[base] = float(parts[1])
                normext[base] = float(parts[2])
                idx += 1

    return sphalbt.tolist(), normext.tolist()


def _load_trans_ascii(path: Path, nsr_bands: int) -> list:
    """Load transmission LUT from ASCII file. Same format as Landsat."""
    NPRES = 7
    NSUNANGLE = 22
    NSUNANGLE_FILE = 21
    NAOT = 22

    transt = np.zeros(nsr_bands * NPRES * NAOT * NSUNANGLE, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        idx += 1  # band header
        for ip in range(NPRES):
            idx += 1  # pressure header
            for isun in range(NSUNANGLE_FILE):
                parts = lines[idx].split()
                for ia in range(NAOT):
                    base = (ib * NPRES * NAOT * NSUNANGLE
                            + ip * NAOT * NSUNANGLE
                            + isun + ia * NSUNANGLE)
                    transt[base] = float(parts[1 + ia])
                idx += 1

    return transt.tolist()
