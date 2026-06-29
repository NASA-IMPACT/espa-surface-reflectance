"""Top-level scene processing pipeline."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lasrc.aux import load_auxiliary_data, load_dem, load_lut, load_ratio_file
from lasrc.io import read_landsat_scene, write_cog_output, write_espa_output, write_numpy
from lasrc.sensors import SENSORS

import lasrc as _lasrc


@dataclass
class AuxFilePaths:
    """Paths to all auxiliary input files required for surface reflectance."""

    # LUT files: angle/intref are HDF4; transm/sphera are ASCII.
    angle_hdf: str | Path
    intref_hdf: str | Path
    transm_hdf: str | Path
    sphera_hdf: str | Path

    # Atmospheric grids
    wv_oz_hdf: str | Path
    dem_hdf: str | Path
    ratio_hdf: str | Path

    # Data source for water vapor / ozone
    aux_source: str = "VIIRS"


def process_scene(
    input_path: str | Path,
    aux_files: AuxFilePaths,
    output_path: str | Path,
    sensor_name: str = "LANDSAT_8",
    output_format: str = "cog",
    use_orig_aero: bool = False,
) -> dict:
    """Process a scene from TOA to surface reflectance.

    Parameters
    ----------
    input_path : path
        Landsat scene directory or file.
    aux_files : AuxFilePaths
        Resolved paths to all auxiliary data files.
    output_path : path
        Output file or directory.
    sensor_name : str
        One of LANDSAT_8, LANDSAT_9, SENTINEL_2A, SENTINEL_2B, SENTINEL_2C.
    output_format : str
        "cog", "espa", or "numpy".
    use_orig_aero : bool
        Use original aerosol algorithm (not yet implemented).
    """
    sensor_config = SENSORS[sensor_name]

    scene = read_landsat_scene(input_path)

    lut_data = load_lut(
        str(aux_files.angle_hdf),
        str(aux_files.intref_hdf),
        str(aux_files.transm_hdf),
        str(aux_files.sphera_hdf),
        band_names=sensor_config["lut_band_names"],
    )
    lut = _lasrc.PyLookupTables(**lut_data)

    aux_data = load_auxiliary_data(
        str(aux_files.wv_oz_hdf),
        aux_source=aux_files.aux_source,
    )
    dem = load_dem(str(aux_files.dem_hdf))
    ratios = load_ratio_file(str(aux_files.ratio_hdf))
    aux = _lasrc.PyAuxiliaryData(dem=dem, **aux_data, **ratios)

    profile = scene["profile"]
    transform = profile["transform"]
    nlines = profile["height"]
    nsamps = profile["width"]

    angles = scene["angles"]
    if "sza" in angles:
        sza = angles["sza"]
        saa = angles["saa"]
        vza = angles["vza"]
        vaa = angles["vaa"]
    else:
        sza = np.full((nlines, nsamps), 90.0 - angles["sun_elevation"], dtype=np.float32)
        saa = np.full((nlines, nsamps), 0.0, dtype=np.float32)
        vza = np.full((nlines, nsamps), 0.0, dtype=np.float32)
        vaa = np.full((nlines, nsamps), 0.0, dtype=np.float32)

    # Extract UTM geotransform for per-pixel coordinate conversion
    # transform.c = UL x, transform.f = UL y, transform.a = pixel width, transform.e = pixel height (neg)
    utm_zone = profile.get("utm_zone", profile["crs"].to_epsg() % 100)
    if profile["crs"].to_epsg() > 32700:
        utm_zone = -utm_zone  # Southern hemisphere

    result = _lasrc.process_surface_reflectance(
        sensor_name=sensor_name,
        toa_bands=scene["toa_bands"],
        bt_bands=scene["bt_bands"],
        solar_zenith=sza,
        solar_azimuth=saa,
        view_zenith=vza,
        view_azimuth=vaa,
        qa_band=scene["qa_band"],
        lut=lut,
        aux=aux,
        ul_corner_x=transform.c,
        ul_corner_y=transform.f,
        pixel_size_x=transform.a,
        pixel_size_y=abs(transform.e),
        utm_zone=utm_zone,
        use_orig_aero=use_orig_aero,
    )

    for i in range(len(result["sr_bands"])):
        result["sr_bands"][i] = result["sr_bands"][i].reshape(nlines, nsamps)
    if "bt_bands" in result:
        for i in range(len(result["bt_bands"])):
            result["bt_bands"][i] = result["bt_bands"][i].reshape(nlines, nsamps)
    result["aerosol"] = result["aerosol"].reshape(nlines, nsamps)
    result["qa"] = result["qa"].reshape(nlines, nsamps)

    if output_format == "espa":
        write_espa_output(output_path, result, sensor_config,
                          crs=profile["crs"], transform=transform)
    elif output_format == "numpy":
        write_numpy(output_path, result, sensor_config)
    else:
        write_cog_output(output_path, result, sensor_config, profile)

    return result
