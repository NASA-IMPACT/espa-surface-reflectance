"""Scene I/O: read input scenes and write output products."""

from pathlib import Path

import numpy as np
import rasterio


def read_landsat_scene(scene_dir: str | Path) -> dict:
    """Read a Landsat Level-1 scene from GeoTIFF files.

    Returns dict with keys: toa_bands, bt_bands, qa_band, angles, metadata, profile.
    """
    scene_dir = Path(scene_dir)
    mtl_files = list(scene_dir.glob("*_MTL.txt")) + list(scene_dir.glob("*_MTL.json"))
    if not mtl_files:
        raise FileNotFoundError(f"No MTL file found in {scene_dir}")

    metadata = _parse_landsat_mtl(mtl_files[0])
    band_files = sorted(scene_dir.glob("*_B*.TIF"))

    toa_bands = []
    bt_bands = []
    profile = None

    for band_num in [1, 2, 3, 4, 5, 6, 7, 9]:
        band_file = _find_band_file(band_files, band_num)
        with rasterio.open(band_file) as src:
            data = src.read(1).astype(np.float32)
            if profile is None:
                profile = src.profile.copy()
            gain = metadata["refl_mult"][band_num]
            bias = metadata["refl_add"][band_num]
            toa = data * gain + bias
            toa_bands.append(toa)

    for band_num in [10, 11]:
        band_file = _find_band_file(band_files, band_num)
        if band_file is not None:
            with rasterio.open(band_file) as src:
                data = src.read(1).astype(np.float32)
                radiance = data * metadata["rad_mult"][band_num] + metadata["rad_add"][band_num]
                k1 = metadata["k1"][band_num]
                k2 = metadata["k2"][band_num]
                bt = np.where(radiance > 0, k2 / np.log(k1 / radiance + 1.0), 0.0)
                bt_bands.append(bt.astype(np.float32))

    qa_files = list(scene_dir.glob("*_QA_PIXEL.TIF"))
    if not qa_files:
        raise FileNotFoundError(f"No *_QA_PIXEL.TIF file found in {scene_dir}")
    with rasterio.open(qa_files[0]) as src:
        qa_band = src.read(1).astype(np.uint16)

    angles = _read_landsat_angles(scene_dir, metadata)

    return {
        "toa_bands": toa_bands,
        "bt_bands": bt_bands,
        "qa_band": qa_band,
        "angles": angles,
        "metadata": metadata,
        "profile": profile,
    }


def _find_band_file(band_files, band_id):
    """Find the file matching a band number or name."""
    pattern = f"_B{band_id}." if isinstance(band_id, int) else f"_{band_id}."
    for f in band_files:
        if pattern in f.name:
            return f
    return None


def _parse_landsat_mtl(mtl_path):
    """Parse Landsat MTL metadata file. Returns dict of calibration params."""
    metadata = {"refl_mult": {}, "refl_add": {}, "rad_mult": {}, "rad_add": {},
                "k1": {}, "k2": {}}
    with open(mtl_path) as f:
        for line in f:
            line = line.strip()
            if "=" not in line:
                continue
            key, val = line.split("=", 1)
            key = key.strip()
            val = val.strip()
            for prefix, target in [
                ("REFLECTANCE_MULT_BAND_", "refl_mult"),
                ("REFLECTANCE_ADD_BAND_", "refl_add"),
                ("RADIANCE_MULT_BAND_", "rad_mult"),
                ("RADIANCE_ADD_BAND_", "rad_add"),
                ("K1_CONSTANT_BAND_", "k1"),
                ("K2_CONSTANT_BAND_", "k2"),
            ]:
                if key.startswith(prefix):
                    band_num = int(key[len(prefix):])
                    metadata[target][band_num] = float(val)
            if key == "SUN_ELEVATION":
                metadata["sun_elevation"] = float(val)
            if key == "SUN_AZIMUTH":
                metadata["sun_azimuth"] = float(val)
    return metadata


def _read_landsat_angles(scene_dir, metadata):
    """Read per-pixel angle bands or compute from scene-level metadata."""
    angle_files = list(scene_dir.glob("*_SZA.TIF"))
    if angle_files:
        angles = {}
        for name, suffix in [("sza", "SZA"), ("saa", "SAA"), ("vza", "VZA"), ("vaa", "VAA")]:
            f = list(scene_dir.glob(f"*_{suffix}.TIF"))
            if f:
                with rasterio.open(f[0]) as src:
                    angles[name] = src.read(1).astype(np.float32) * 0.01
        return angles
    else:
        return {"sun_elevation": metadata.get("sun_elevation", 45.0)}


def write_espa_output(output_dir: str | Path, result: dict,
                      metadata: dict, profile: dict,
                      sensor_config: dict) -> None:
    """Write output in ESPA internal format (flat binary + XML)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, name in enumerate(sensor_config["refl_band_names"]):
        band_data = result["sr_bands"][i]
        band_data.tofile(output_dir / f"{name}.img")

    result["aerosol"].tofile(output_dir / "sr_aerosol.img")
    result["qa"].tofile(output_dir / "sr_aerosol_qa.img")


def write_cog_output(output_path: str | Path, result: dict,
                     metadata: dict, profile: dict,
                     sensor_config: dict) -> None:
    """Write output as Cloud-Optimized GeoTIFF."""
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
    )

    with rasterio.open(output_path, "w", **cog_profile) as dst:
        for i, band in enumerate(result["sr_bands"]):
            dst.write(band.astype(np.int16), i + 1)
            dst.set_band_description(i + 1, sensor_config["refl_band_names"][i])
        dst.write(result["aerosol"].astype(np.int16), nbands + 1)
        dst.set_band_description(nbands + 1, "sr_aerosol")
        dst.write(result["qa"].astype(np.uint8), nbands + 2)
        dst.set_band_description(nbands + 2, "sr_aerosol_qa")
