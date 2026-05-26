"""Top-level scene processing pipeline."""

from pathlib import Path

import numpy as np

from lasrc.aux import load_auxiliary_data, load_dem, load_lut_from_hdf, load_ratio_file
from lasrc.io import read_landsat_scene, write_cog_output, write_espa_output
from lasrc.sensors import SENSORS

import lasrc as _lasrc


def process_scene(
    input_path: str | Path,
    aux_dir: str | Path,
    output_path: str | Path,
    sensor_name: str = "LANDSAT_8",
    output_format: str = "cog",
    use_orig_aero: bool = False,
    aux_source: str = "VIIRS",
) -> None:
    """Process a scene from TOA to surface reflectance."""
    sensor_config = SENSORS[sensor_name]
    aux_dir = Path(aux_dir)

    scene = read_landsat_scene(input_path)

    lut_data = load_lut_from_hdf(
        str(aux_dir / "ANGLE_NEW.hdf"),
        str(aux_dir / "RES_LUT_V3.0-LANDSAT.hdf"),
        str(aux_dir / "TRANS_LUT_V3.0-LANDSAT.hdf"),
        str(aux_dir / "AERO_LUT_V3.0-LANDSAT.hdf"),
        nsr_bands=sensor_config["nsr_bands"],
    )
    lut = _lasrc.PyLookupTables(**lut_data)

    aux_data = load_auxiliary_data(
        str(aux_dir / "VIIRS_CMG_DAILY.hdf"),
        aux_source=aux_source,
    )
    dem = load_dem(str(aux_dir / "CMGDEM.hdf"))
    ratios = load_ratio_file(str(aux_dir / "ratiomapndwiexp.hdf"))
    aux = _lasrc.PyAuxiliaryData(dem=dem, **aux_data, **ratios)

    profile = scene["profile"]
    transform = profile["transform"]
    nlines = profile["height"]
    nsamps = profile["width"]
    center_x = transform.c + nsamps / 2 * transform.a
    center_y = transform.f + nlines / 2 * transform.e

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
        scene_center_lat=center_y,
        scene_center_lon=center_x,
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
        write_espa_output(output_path, result, scene["metadata"], profile, sensor_config)
    else:
        write_cog_output(output_path, result, scene["metadata"], profile, sensor_config)
