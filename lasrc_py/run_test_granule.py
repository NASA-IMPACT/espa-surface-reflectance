"""Run LaSRC surface reflectance on the test Landsat 9 granule.

Thin wrapper around the lasrc package pipeline (the single source of truth for
LUT/ancillary loading and correction). See lasrc.pipeline.process_scene.
"""

import argparse
import glob as globmod
import time
from pathlib import Path

# -- Paths -----------------------------------------------------------------------
TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"
SCENE_DIR = TEST_DATA / "LC09_L1TP_001066_20260314_20260314_02_T1"
AUX_DIR = TEST_DATA / "aux_data"
LUT_DIR = AUX_DIR / "LDCMLUT"

# Scene date: 2026-03-14 -> DOY 073
# VIIRS aux file: $LASRC_AUX_DIR/LADS/<year>/V*04ANC.A<year><doy>.*.h5
LADS_DIR = AUX_DIR / "LADS" / "2026"

SENSOR = "LANDSAT_9"


def find_viirs_aux(lads_dir: Path, year: str, doy: str) -> Path:
    """Find the VIIRS auxiliary file matching V*04ANC.A<year><doy>.*.h5."""
    pattern = str(lads_dir / f"V*04ANC.A{year}{doy}.*.h5")
    matches = sorted(globmod.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No VIIRS aux file found matching {pattern}. "
            f"Use --aux to specify the path explicitly."
        )
    return Path(matches[0])


def main():
    parser = argparse.ArgumentParser(description="Run LaSRC on test Landsat 9 granule")
    parser.add_argument(
        "--format", choices=["cog", "espa"], default="cog",
        help="Output format: 'cog' for GeoTIFF (default), 'espa' for flat binary .img files",
    )
    parser.add_argument(
        "--aux", type=Path, default=None,
        help="Path to VIIRS auxiliary HDF5 file (auto-detected from LADS dir if omitted)",
    )
    args = parser.parse_args()

    from lasrc.pipeline import AuxFilePaths, process_scene

    viirs_path = args.aux if args.aux is not None else find_viirs_aux(LADS_DIR, "2026", "073")
    print(f"Scene:     {SCENE_DIR.name}")
    print(f"VIIRS aux: {viirs_path.name}")
    print(f"Format:    {args.format}")

    aux_files = AuxFilePaths(
        angle_hdf=LUT_DIR / "ANGLE_NEW.hdf",
        intref_hdf=LUT_DIR / "RES_LUT_V3.0-URBANCLEAN-V2.0.hdf",
        transm_hdf=LUT_DIR / "TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
        sphera_hdf=LUT_DIR / "AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
        wv_oz_hdf=viirs_path,
        dem_hdf=AUX_DIR / "CMGDEM.hdf",
        ratio_hdf=AUX_DIR / "ratiomapndwiexp.hdf",
        aux_source="VIIRS",
    )

    if args.format == "espa":
        output_path = TEST_DATA / "output_espa"
    else:
        output_path = TEST_DATA / "output_sr.tif"

    print("Running surface reflectance correction...")
    t0 = time.time()
    result = process_scene(
        input_path=SCENE_DIR,
        aux_files=aux_files,
        output_path=output_path,
        sensor_name=SENSOR,
        output_format=args.format,
    )
    print(f"  Correction completed in {time.time() - t0:.1f}s")

    sr0 = result["sr_bands"][0]
    aerosol = result["aerosol"]
    print(f"  SR band 1 stats: min={sr0.min()}, max={sr0.max()}, mean={sr0.mean():.1f}")
    print(f"  Aerosol stats:   min={aerosol.min()}, max={aerosol.max()}, "
          f"mean={aerosol.mean():.1f}")
    print(f"Output written to {output_path}")
    print("Done!")


if __name__ == "__main__":
    main()
