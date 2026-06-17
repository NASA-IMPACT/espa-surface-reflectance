"""Sentinel-2 processing pipeline: SAFE archive → surface reflectance."""

from pathlib import Path

import numpy as np

from lasrc.aux import load_auxiliary_data, load_dem, load_lut, load_ratio_file
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
) -> dict:
    """Process a Sentinel-2 SAFE archive to surface reflectance."""
    sensor_config = SENSORS[sensor_name]
    nsr_bands = sensor_config["nsr_bands"]  # 13

    # Read SAFE archive
    scene = read_sentinel_safe(safe_dir)
    angles = scene["angles"]
    profile = scene["profile"]
    nlines = profile["nlines"]
    nsamps = profile["nsamps"]

    # Load LUTs (HDF4 angle + intref, ASCII transm + sphera)
    lut_data = load_lut(
        str(angle_hdf), str(intref_hdf), str(transm_hdf), str(sphera_hdf),
        band_names=sensor_config["lut_band_names"],
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

    return result


def _write_espa_sentinel(output_dir, result, sensor_config):
    """Write Sentinel SR output in ESPA format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, name in enumerate(sensor_config["refl_band_names"]):
        result["sr_bands"][i].astype(np.uint16).tofile(output_dir / f"{name}.img")

    result["aerosol"].astype(np.int16).tofile(output_dir / "sr_aerosol.img")
    result["qa"].astype(np.uint8).tofile(output_dir / "sr_aerosol_qa.img")

