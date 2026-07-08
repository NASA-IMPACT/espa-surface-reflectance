"""Output writers for surface reflectance products.

Three formats are supported, all driven by the same in-memory ``result`` dict
(keys: sr_bands, aerosol, qa) and the sensor configuration:

  - ESPA internal format : one ENVI file per band (flat binary ``.img`` plus a
    text ``.hdr`` sidecar), georeferenced. This is the real ESPA product format.
  - COG                  : a single multi-band Cloud-Optimized GeoTIFF.
  - NumPy raw            : headerless flat binary ``.img`` files (``ndarray.tofile``),
    kept only for compatibility with the analysis/test scripts.

Output dtypes match the C reference exactly: surface reflectance (and
brightness temperature) bands are uint16, the aerosol band is int16, and the
aerosol QA band is uint8.
"""

from pathlib import Path

import numpy as np
import rasterio

from lasrc.lasrc import AERO_FILL, SR_FILL_VALUE


def _write_band(path, data, crs, transform, driver, nodata=None, **creation_opts) -> None:
    """Write a single 2-D array as a georeferenced raster band."""
    height, width = data.shape
    with rasterio.open(
        path,
        "w",
        driver=driver,
        width=width,
        height=height,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
        **creation_opts,
    ) as dst:
        dst.write(data, 1)


def write_espa_output(output_dir: str | Path, result: dict,
                      sensor_config: dict, crs, transform,
                      product_id: str = "") -> None:
    """Write output in ESPA internal format: one ENVI band per file.

    Each surface reflectance band, the aerosol band, and the aerosol QA band
    are written as separate ENVI files (flat binary ``.img`` with a ``.hdr``
    header), georeferenced via ``crs`` and ``transform``.

    When ``product_id`` is given, each filename is prefixed with it, e.g.
    ``LC08_L1TP_230094_20250102_20250110_02_T1_sr_band5.img``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{product_id}_" if product_id else ""

    for i, name in enumerate(sensor_config["refl_band_names"]):
        band = result["sr_bands"][i].astype(np.uint16)
        _write_band(output_dir / f"{prefix}{name}.img", band, crs, transform,
                    driver="ENVI", nodata=SR_FILL_VALUE)

    _write_band(output_dir / f"{prefix}sr_aerosol.img",
               result["aerosol"].astype(np.int16), crs, transform,
               driver="ENVI", nodata=AERO_FILL)
    _write_band(output_dir / f"{prefix}sr_aerosol_qa.img",
               result["qa"].astype(np.uint8), crs, transform,
               driver="ENVI")


_COG_CREATION_OPTS = dict(compress="deflate", tiled=True, blockxsize=512, blockysize=512)


def write_cog_output(output_dir: str | Path, result: dict,
                     sensor_config: dict, crs, transform,
                     product_id: str = "") -> None:
    """Write output as one Cloud-Optimized GeoTIFF per band.

    Mirrors write_espa_output's one-file-per-band layout (and the C
    reference's per-band native dtypes: uint16 SR/BT, int16 aerosol, uint8
    QA), just using the GTiff/COG driver instead of ENVI.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{product_id}_" if product_id else ""

    for i, name in enumerate(sensor_config["refl_band_names"]):
        band = result["sr_bands"][i].astype(np.uint16)
        _write_band(output_dir / f"{prefix}{name}.tif", band, crs, transform,
                    driver="GTiff", nodata=SR_FILL_VALUE, **_COG_CREATION_OPTS)

    _write_band(output_dir / f"{prefix}sr_aerosol.tif",
               result["aerosol"].astype(np.int16), crs, transform,
               driver="GTiff", nodata=AERO_FILL, **_COG_CREATION_OPTS)
    _write_band(output_dir / f"{prefix}sr_aerosol_qa.tif",
               result["qa"].astype(np.uint8), crs, transform,
               driver="GTiff", **_COG_CREATION_OPTS)


def write_numpy(output_dir: str | Path, result: dict,
                sensor_config: dict, product_id: str = "") -> None:
    """Write raw headerless flat binary ``.bin`` files (ndarray.tofile)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{product_id}_" if product_id else ""

    for i, name in enumerate(sensor_config["refl_band_names"]):
        result["sr_bands"][i].astype(np.uint16).tofile(output_dir / f"{prefix}{name}.bin")

    result["aerosol"].astype(np.int16).tofile(output_dir / f"{prefix}sr_aerosol.bin")
    result["qa"].astype(np.uint8).tofile(output_dir / f"{prefix}sr_aerosol_qa.bin")
