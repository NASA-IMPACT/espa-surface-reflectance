"""Run LaSRC surface reflectance on the test Landsat 9 granule."""

import sys
import time
from pathlib import Path

import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────────
TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"
SCENE_DIR = TEST_DATA / "LC09_L1TP_001066_20260314_20260314_02_T1"
AUX_DIR = TEST_DATA / "aux_data"
LUT_DIR = AUX_DIR / "LDCMLUT"
OUTPUT_PATH = TEST_DATA / "output_sr.tif"

# Scene date: 2026-03-14 -> month 03
WV_IMG = AUX_DIR / "monthly_avgs" / "2026" / "monthly_avg_wv_2026_03.img"
OZ_IMG = AUX_DIR / "monthly_avgs" / "2026" / "monthly_avg_oz_2026_03.img"


# ── LUT loading ───────────────────────────────────────────────────────────────

def load_angle_lut(path: Path) -> dict:
    """Load angle LUT from HDF4 file."""
    from pyhdf.SD import SD, SDC

    hdf = SD(str(path), SDC.READ)
    result = {}
    for name in ["TSMAX", "TSMIN", "NBFIC", "TTV", "TTS"]:
        ds = hdf.select(name)
        result[name.lower()] = ds.get().astype(np.float64).ravel().tolist()
    # NBFI is integer
    result["nbfi"] = hdf.select("NBFI").get().astype(np.int32).ravel().tolist()
    hdf.end()
    return result


def load_rolutt(path: Path, nsr_bands: int) -> list:
    """Load intrinsic reflectance LUT from HDF4.

    HDF4 stores per-band as NRLUT_BAND_N with shape [NSOLAR, NAOT, NPRES].
    We need flat layout: band × NPRES × NAOT × NSOLAR.
    """
    from pyhdf.SD import SD, SDC

    NSOLAR = 8000
    NAOT = 22
    NPRES = 7

    hdf = SD(str(path), SDC.READ)
    rolutt = np.zeros(nsr_bands * NPRES * NAOT * NSOLAR, dtype=np.float64)

    for ib in range(nsr_bands):
        ds_name = f"NRLUT_BAND_{ib + 1}"
        data = hdf.select(ds_name).get().astype(np.float64)  # [NSOLAR, NAOT, NPRES]

        # Rearrange to [NPRES, NAOT, NSOLAR] then flatten
        for ip in range(NPRES):
            for ia in range(NAOT):
                base = ib * NPRES * NAOT * NSOLAR + ip * NAOT * NSOLAR + ia * NSOLAR
                rolutt[base:base + NSOLAR] = data[:, ia, ip]

    hdf.end()
    return rolutt.tolist()


def load_aero_ascii(path: Path, nsr_bands: int) -> tuple[list, list]:
    """Load spherical albedo and normalized extinction from ASCII file.

    Format: per band, 7 pressure blocks, 22 lines per block (aot sphalbt normext).
    Returns (sphalbt_flat, normext_flat) in band × NPRES × NAOT layout.
    """
    NPRES = 7
    NAOT = 22

    sphalbt = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)
    normext = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        # Skip band header line
        idx += 1
        for ip in range(NPRES):
            # Skip pressure header line
            idx += 1
            for ia in range(NAOT):
                parts = lines[idx].split()
                base = ib * NPRES * NAOT + ip * NAOT + ia
                sphalbt[base] = float(parts[1])
                normext[base] = float(parts[2])
                idx += 1

    return sphalbt.tolist(), normext.tolist()


def load_trans_ascii(path: Path, nsr_bands: int) -> list:
    """Load transmission LUT from ASCII file.

    Format: per band, 7 pressure blocks, 21 sun angle lines per block,
    each line has sun_angle followed by 22 AOT transmission values.
    The C code allocates NSUNANGLE_VALS=22 slots but only fills 21 from file.
    Layout: band × NPRES × NAOT × NSUNANGLE.
    """
    NPRES = 7
    NSUNANGLE = 22  # Array dimension (C code allocates 22)
    NSUNANGLE_FILE = 21  # Actual lines in the file per pressure block
    NAOT = 22

    transt = np.zeros(nsr_bands * NPRES * NAOT * NSUNANGLE, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        # Skip band header line
        idx += 1
        for ip in range(NPRES):
            # Skip pressure header line
            idx += 1
            for isun in range(NSUNANGLE_FILE):
                parts = lines[idx].split()
                # First value is sun angle, rest are transmission for each AOT
                for ia in range(NAOT):
                    base = (ib * NPRES * NAOT * NSUNANGLE
                            + ip * NAOT * NSUNANGLE
                            + isun + ia * NSUNANGLE)
                    transt[base] = float(parts[1 + ia])
                idx += 1

    return transt.tolist()


# ── Auxiliary data loading ─────────────────────────────────────────────────────

def load_monthly_wv_oz(wv_path: Path, oz_path: Path) -> dict:
    """Load monthly average water vapor and ozone from ENVI binary files.

    WV: uint16, 3600×7200, scale=200 (divide to get cm)
    OZ: uint8, 3600×7200, scale=400 (divide to get DU... but stored as uint8)
    """
    wv = np.fromfile(wv_path, dtype=np.uint16).reshape(3600, 7200).astype(np.int16)
    oz = np.fromfile(oz_path, dtype=np.uint8).reshape(3600, 7200).astype(np.int16)

    return {
        "wv": wv.ravel().tolist(),
        "oz": oz.ravel().tolist(),
        "wv_scale": 200.0,
        "oz_scale": 400.0,
        "wv_default": 2.5,
        "oz_default": 0.3,
    }


def load_dem(path: Path) -> list:
    """Load CMG DEM from HDF4."""
    from pyhdf.SD import SD, SDC

    hdf = SD(str(path), SDC.READ)
    dem = hdf.select("averaged elevation").get().astype(np.int16).ravel()
    hdf.end()
    return dem.tolist()


def load_ratios(path: Path) -> dict:
    """Load band ratios and NDWI from HDF4.

    The C code maps Landsat band names to internal ratio names:
      ratiob1 = "average ratio b9"   (Landsat band 9 = coastal aerosol proxy)
      ratiob2 = "average ratio b10"  (Landsat band 10 proxy)
      ratiob7 = "average ratio b7"
    """
    from pyhdf.SD import SD, SDC

    # Map from our internal names to HDF4 dataset names
    name_map = {
        "ratiob1": "average ratio b9",
        "ratiob2": "average ratio b10",
        "ratiob7": "average ratio b7",
        "intratiob1": "inter ratiob9",
        "intratiob2": "inter ratiob10",
        "intratiob7": "inter ratiob7",
        "slpratiob1": "slope ratiob9",
        "slpratiob2": "slope ratiob10",
        "slpratiob7": "slope ratiob7",
        "andwi": "average ndvi",
        "sndwi": "standard ndvi",
    }

    hdf = SD(str(path), SDC.READ)
    result = {}
    for key, sds_name in name_map.items():
        result[key] = np.array(hdf.select(sds_name).get(), dtype=np.int16).ravel().tolist()
    hdf.end()
    return result


# ── Scene reading ──────────────────────────────────────────────────────────────

def read_scene(scene_dir: Path) -> dict:
    """Read Landsat scene using rasterio."""
    import rasterio

    mtl_path = next(scene_dir.glob("*_MTL.txt"))
    metadata = parse_mtl(mtl_path)
    band_files = sorted(scene_dir.glob("*_B*.TIF"))

    toa_bands = []
    profile = None
    for band_num in [1, 2, 3, 4, 5, 6, 7, 9]:
        band_file = find_band_file(band_files, band_num)
        with rasterio.open(band_file) as src:
            data = src.read(1).astype(np.float32)
            if profile is None:
                profile = src.profile.copy()
            gain = metadata["refl_mult"][band_num]
            bias = metadata["refl_add"][band_num]
            toa_bands.append(data * gain + bias)

    bt_bands = []
    for band_num in [10, 11]:
        band_file = find_band_file(band_files, band_num)
        if band_file is not None:
            with rasterio.open(band_file) as src:
                data = src.read(1).astype(np.float32)
                radiance = data * metadata["rad_mult"][band_num] + metadata["rad_add"][band_num]
                k1 = metadata["k1"][band_num]
                k2 = metadata["k2"][band_num]
                bt = np.where(radiance > 0, k2 / np.log(k1 / radiance + 1.0), 0.0)
                bt_bands.append(bt.astype(np.float32))

    qa_file = next(scene_dir.glob("*_QA_PIXEL.TIF"))
    with rasterio.open(qa_file) as src:
        qa_band = src.read(1).astype(np.uint16)

    # Per-pixel angle bands (scaled by 0.01 in file)
    angles = {}
    for name, suffix in [("sza", "SZA"), ("saa", "SAA"), ("vza", "VZA"), ("vaa", "VAA")]:
        f = list(scene_dir.glob(f"*_{suffix}.TIF"))
        if f:
            with rasterio.open(f[0]) as src:
                angles[name] = src.read(1).astype(np.float32) * 0.01

    return {
        "toa_bands": toa_bands,
        "bt_bands": bt_bands,
        "qa_band": qa_band,
        "angles": angles,
        "profile": profile,
    }


def parse_mtl(mtl_path: Path) -> dict:
    """Parse Landsat MTL text file for calibration coefficients."""
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
    return metadata


def find_band_file(band_files, band_id):
    """Find file matching a band number or name."""
    pattern = f"_B{band_id}." if isinstance(band_id, int) else f"_{band_id}."
    for f in band_files:
        if pattern in f.name:
            return f
    return None


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    import lasrc as _lasrc

    nsr_bands = 8
    print(f"LaSRC version: {_lasrc.__version__}")

    # Load LUTs
    print("Loading LUTs...")
    t0 = time.time()
    angle_data = load_angle_lut(LUT_DIR / "ANGLE_NEW.hdf")
    rolutt = load_rolutt(LUT_DIR / "RES_LUT_V3.0-URBANCLEAN-V2.0.hdf", nsr_bands)
    sphalbt, normext = load_aero_ascii(LUT_DIR / "AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII", nsr_bands)
    transt = load_trans_ascii(LUT_DIR / "TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII", nsr_bands)
    print(f"  LUTs loaded in {time.time() - t0:.1f}s")

    lut = _lasrc.PyLookupTables(
        rolutt=rolutt,
        transt=transt,
        sphalbt=sphalbt,
        normext=normext,
        tsmax=angle_data["tsmax"],
        tsmin=angle_data["tsmin"],
        nbfic=angle_data["nbfic"],
        nbfi=angle_data["nbfi"],
        ttv=angle_data["ttv"],
        tts=angle_data["tts"],
        nsr_bands=nsr_bands,
    )

    # Load auxiliary data
    print("Loading auxiliary data...")
    t0 = time.time()
    wv_oz = load_monthly_wv_oz(WV_IMG, OZ_IMG)
    dem = load_dem(AUX_DIR / "CMGDEM.hdf")
    ratios = load_ratios(AUX_DIR / "ratiomapndwiexp.hdf")
    print(f"  Aux data loaded in {time.time() - t0:.1f}s")

    aux = _lasrc.PyAuxiliaryData(dem=dem, **wv_oz, **ratios)

    # Read scene
    print(f"Reading scene: {SCENE_DIR.name}")
    t0 = time.time()
    scene = read_scene(SCENE_DIR)
    profile = scene["profile"]
    nlines = profile["height"]
    nsamps = profile["width"]
    print(f"  Scene: {nlines} x {nsamps}, {len(scene['toa_bands'])} TOA bands, "
          f"{len(scene['bt_bands'])} BT bands")
    print(f"  Scene read in {time.time() - t0:.1f}s")

    # Extract UTM geotransform
    transform = profile["transform"]
    epsg = profile["crs"].to_epsg()
    utm_zone = epsg % 100
    if epsg > 32700:
        utm_zone = -utm_zone
    print(f"  EPSG: {epsg}, UTM zone: {utm_zone}")
    print(f"  UL corner: ({transform.c}, {transform.f}), pixel size: ({transform.a}, {transform.e})")

    # Run surface reflectance
    print("Running surface reflectance correction...")
    t0 = time.time()
    result = _lasrc.process_surface_reflectance(
        sensor_name="LANDSAT_9",
        toa_bands=scene["toa_bands"],
        bt_bands=scene["bt_bands"],
        solar_zenith=scene["angles"]["sza"],
        solar_azimuth=scene["angles"]["saa"],
        view_zenith=scene["angles"]["vza"],
        view_azimuth=scene["angles"]["vaa"],
        qa_band=scene["qa_band"],
        lut=lut,
        aux=aux,
        ul_corner_x=transform.c,
        ul_corner_y=transform.f,
        pixel_size_x=transform.a,
        pixel_size_y=abs(transform.e),
        utm_zone=utm_zone,
        use_orig_aero=False,
    )
    elapsed = time.time() - t0
    print(f"  Correction completed in {elapsed:.1f}s")

    # Reshape and write output
    import rasterio

    sr_bands = [b.reshape(nlines, nsamps) for b in result["sr_bands"]]
    aerosol = result["aerosol"].reshape(nlines, nsamps)
    qa = result["qa"].reshape(nlines, nsamps)

    print(f"  SR band 1 stats: min={sr_bands[0].min()}, max={sr_bands[0].max()}, "
          f"mean={sr_bands[0].mean():.1f}")
    print(f"  Aerosol stats: min={aerosol.min()}, max={aerosol.max()}, "
          f"mean={aerosol.mean():.1f}")

    # Write COG output
    nbands = len(sr_bands)
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

    print(f"Writing output to {OUTPUT_PATH}...")
    with rasterio.open(OUTPUT_PATH, "w", **cog_profile) as dst:
        for i, band in enumerate(sr_bands):
            dst.write(band.astype(np.int16), i + 1)
            dst.set_band_description(i + 1, f"sr_band{i + 1}")
        dst.write(aerosol.astype(np.int16), nbands + 1)
        dst.set_band_description(nbands + 1, "sr_aerosol")
        dst.write(qa.astype(np.uint8), nbands + 2)
        dst.set_band_description(nbands + 2, "sr_aerosol_qa")

    print("Done!")


if __name__ == "__main__":
    main()
