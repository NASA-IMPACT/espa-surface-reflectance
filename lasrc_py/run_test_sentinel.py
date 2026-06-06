"""Run LaSRC Sentinel-2 surface reflectance on the S2B test granule."""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"
SAFE_DIR = TEST_DATA / "S2B_MSIL1C_20260124T073109_N0511_R049_T38PNC_20260124T092241.SAFE"
AUX_DIR = TEST_DATA / "aux_data"
LUT_DIR = AUX_DIR / "MSILUT"
LADS_DIR = AUX_DIR / "LADS" / "2026"


def find_viirs_aux(lads_dir: Path, year: str, doy: str) -> Path:
    import glob as globmod
    pattern = str(lads_dir / f"V*04ANC.A{year}{doy}.*.h5")
    matches = sorted(globmod.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No VIIRS aux file found matching {pattern}")
    return Path(matches[0])


def main():
    parser = argparse.ArgumentParser(description="Run LaSRC on S2B test granule")
    parser.add_argument("--format", choices=["espa"], default="espa")
    parser.add_argument("--aux", type=Path, default=None)
    args = parser.parse_args()

    from lasrc.pipeline_sentinel import process_sentinel_scene

    # Scene date: 2026-01-24 -> DOY 024
    if args.aux:
        viirs_path = args.aux
    else:
        viirs_path = find_viirs_aux(LADS_DIR, "2026", "024")

    print(f"SAFE dir: {SAFE_DIR.name}")
    print(f"LUT dir: {LUT_DIR}")
    print(f"VIIRS aux: {viirs_path.name}")

    output_dir = TEST_DATA / "output_espa_s2"
    t0 = time.time()

    process_sentinel_scene(
        safe_dir=SAFE_DIR,
        lut_dir=LUT_DIR,
        viirs_aux_path=viirs_path,
        dem_path=AUX_DIR / "CMGDEM.hdf",
        ratio_path=AUX_DIR / "ratiomapndwiexp.hdf",
        output_path=output_dir,
        sensor_name="SENTINEL_2B",
        output_format="espa",
    )

    elapsed = time.time() - t0
    print(f"Completed in {elapsed:.1f}s")
    print(f"Output: {output_dir}")


if __name__ == "__main__":
    main()
