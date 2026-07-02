"""Output writers for surface reflectance products.

Three formats are supported, all driven by the same in-memory ``result`` dict
(keys: sr_bands, aerosol, qa) and the sensor configuration:

  - ESPA internal format : one ENVI file per band (flat binary ``.img`` plus a
    text ``.hdr`` sidecar), georeferenced. This is the real ESPA product format.
  - COG                  : a single multi-band Cloud-Optimized GeoTIFF.
  - NumPy raw            : headerless flat binary ``.img`` files (``ndarray.tofile``),
    kept only for compatibility with the analysis/test scripts.

Output dtypes follow the ESPA convention: surface reflectance and aerosol are
int16, the aerosol QA band is uint8.
"""

from pathlib import Path

import numpy as np
import rasterio

from lasrc.lasrc import SR_FILL_VALUE


def _write_envi_band(path, data, crs, transform, nodata=None) -> None:
    """Write a single 2-D array as a georeferenced ENVI band (.img + .hdr)."""
    height, width = data.shape
    with rasterio.open(
        path,
        "w",
        driver="ENVI",
        width=width,
        height=height,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
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
        band = result["sr_bands"][i].astype(np.int16)
        _write_envi_band(output_dir / f"{prefix}{name}.img", band, crs, transform,
                         nodata=SR_FILL_VALUE)

    _write_envi_band(output_dir / f"{prefix}sr_aerosol.img",
                     result["aerosol"].astype(np.int16), crs, transform)
    _write_envi_band(output_dir / f"{prefix}sr_aerosol_qa.img",
                     result["qa"].astype(np.uint8), crs, transform)


def write_cog_output(output_path: str | Path, result: dict,
                     sensor_config: dict, profile: dict) -> None:
    """Write output as a single multi-band Cloud-Optimized GeoTIFF."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    nbands = len(result["sr_bands"])

    cog_profile = profile.copy()
    cog_profile.update(
        driver="GTiff",
        dtype="int16",
        count=nbands + 2,
        compress="deflate",
        tiled=True,
        blockxsize=512,
        blockysize=512,
        nodata=SR_FILL_VALUE,
    )

    with rasterio.open(output_path, "w", **cog_profile) as dst:
        for i, band in enumerate(result["sr_bands"]):
            dst.write(band.astype(np.int16), i + 1)
            dst.set_band_description(i + 1, sensor_config["refl_band_names"][i])
        dst.write(result["aerosol"].astype(np.int16), nbands + 1)
        dst.set_band_description(nbands + 1, "sr_aerosol")
        dst.write(result["qa"].astype(np.uint8), nbands + 2)
        dst.set_band_description(nbands + 2, "sr_aerosol_qa")


def write_numpy(output_dir: str | Path, result: dict,
                sensor_config: dict, product_id: str = "") -> None:
    """Write raw headerless flat binary ``.bin`` files (ndarray.tofile)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{product_id}_" if product_id else ""

    for i, name in enumerate(sensor_config["refl_band_names"]):
        result["sr_bands"][i].astype(np.int16).tofile(output_dir / f"{prefix}{name}.bin")

    result["aerosol"].astype(np.int16).tofile(output_dir / f"{prefix}sr_aerosol.bin")
    result["qa"].astype(np.uint8).tofile(output_dir / f"{prefix}sr_aerosol_qa.bin")
