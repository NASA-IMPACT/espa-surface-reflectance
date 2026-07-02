"""LUT IO cross-check, part 2: diff C readluts dumps vs Python lasrc.aux.load_lut.

The container-side `dump_c_readluts` test writes the arrays C `readluts` produces
from the four Landsat aux files. This script loads the SAME four files with the
Python loader and byte-compares, catching reshape/stride/band-remap divergence
in load_lut (the IO class the leaf-kernel differential tests cannot see).

Run on a host with the `lasrc` package importable (the repo .venv) and the aux
files available:

    python lasrc_cref/tools/compare_readluts.py \
        --dump   /path/to/dump_dir \
        --angle  ANGLE_NEW.hdf \
        --intref RES_LUT_*.hdf \
        --transm TRANS_LUT_*.ASCII \
        --sphera AERO_LUT_*.ASCII

Exit code 0 = all arrays match, 1 = at least one diverges.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from lasrc.aux import load_lut
from lasrc.sensors import SENSORS

# C readluts dump file -> load_lut dict key. Angle arrays + the four big LUTs.
FLOAT_ARRAYS = {
    "rolutt.bin": "rolutt",
    "transt.bin": "transt",
    "sphalbt.bin": "sphalbt",
    "normext.bin": "normext",
    "tsmax.bin": "tsmax",
    "tsmin.bin": "tsmin",
    "ttv.bin": "ttv",
    "tts.bin": "tts",
    "nbfic.bin": "nbfic",
    # C readluts stores nbfi as float; Python load_angle_lut casts it to int.
    # Values are integral, so compare as float (Python widened back).
    "nbfi.bin": "nbfi",
}
INT_ARRAYS = {
    # NOTE: C readluts reads INDTS with edges[0]=20 (first 20 only); Python
    # reads the full 22-element dataset, so indts[20:22] diverge (extra valid
    # tail values C never uses). Left here on purpose so the diff surfaces it.
    "indts.bin": "indts",
}


def compare(name, c_arr, py_arr):
    """Report exact-match status for one array; return True if identical."""
    if c_arr.shape != py_arr.shape:
        print(f"  FAIL {name}: size C={c_arr.size} Python={py_arr.size}")
        return False
    if np.array_equal(c_arr, py_arr):
        print(f"  ok   {name}: {c_arr.size} values identical")
        return True
    diff = np.flatnonzero(c_arr != py_arr)
    i = int(diff[0])
    print(f"  FAIL {name}: {diff.size}/{c_arr.size} differ; "
          f"first at [{i}] C={c_arr[i]} Python={py_arr[i]}")
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, type=Path)
    ap.add_argument("--angle", required=True)
    ap.add_argument("--intref", required=True)
    ap.add_argument("--transm", required=True)
    ap.add_argument("--sphera", required=True)
    ap.add_argument("--sensor", default="LANDSAT_8")
    args = ap.parse_args()

    band_names = SENSORS[args.sensor]["lut_band_names"]
    lut = load_lut(args.angle, args.intref, args.transm, args.sphera, band_names)

    ok = True
    for fname, key in FLOAT_ARRAYS.items():
        c = np.fromfile(args.dump / fname, dtype=np.float32)
        py = np.asarray(lut[key], dtype=np.float32)
        ok &= compare(key, c, py)
    for fname, key in INT_ARRAYS.items():
        c = np.fromfile(args.dump / fname, dtype=np.int32)
        py = np.asarray(lut[key], dtype=np.int32)
        ok &= compare(key, c, py)

    print("ALL MATCH" if ok else "DIVERGENCE FOUND")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
