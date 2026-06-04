"""Compare Rust and C ESPA output .img files."""
import argparse
import numpy as np
import os
import sys


DEFAULT_RUST_DIR = os.path.join(
    os.path.dirname(__file__), "..", "test_data", "output_espa"
)
DEFAULT_C_DIR = os.path.join(
    os.path.dirname(__file__), "..", "test_data", "c_output"
)
DEFAULT_PREFIX = "LC09_L1TP_001066_20260314_20260314_02_T1_"

BANDS = [
    ("sr_band1.img", np.uint16),
    ("sr_band2.img", np.uint16),
    ("sr_band3.img", np.uint16),
    ("sr_band4.img", np.uint16),
    ("sr_band5.img", np.uint16),
    ("sr_band6.img", np.uint16),
    ("sr_band7.img", np.uint16),
    ("sr_aerosol.img", np.int16),
    ("sr_aerosol_qa.img", np.uint8),
]

NROWS, NCOLS = 7771, 7641


def load_band(directory, fname, dtype, prefix=""):
    path = os.path.join(directory, prefix + fname)
    if not os.path.exists(path):
        return None
    return np.fromfile(path, dtype=dtype)


def compare(rust_dir, c_dir, c_prefix):
    # Load aerosol to build valid-pixel mask (non-fill)
    c_aero = load_band(c_dir, "sr_aerosol.img", np.int16, c_prefix)
    rust_aero = load_band(rust_dir, "sr_aerosol.img", np.int16)
    if c_aero is None or rust_aero is None:
        print("ERROR: Could not load aerosol bands for valid-pixel mask.")
        sys.exit(1)

    valid = (c_aero != -9999) & (rust_aero != -9999)
    n_valid = valid.sum()
    n_total = len(c_aero)
    print(f"Scene: {NROWS} x {NCOLS} = {n_total} pixels, {n_valid} valid ({100*n_valid/n_total:.1f}%)")
    print()

    header = f"{'Band':<18} {'% differ':>8} {'Median':>7} {'Mean':>7} {'P95':>7} {'P99':>7} {'Max':>7} {'Center':>7}"
    print(header)
    print("-" * len(header))

    center_idx = (NROWS // 2) * NCOLS + (NCOLS // 2)

    for fname, dtype in BANDS:
        rust = load_band(rust_dir, fname, dtype)
        c = load_band(c_dir, fname, dtype, c_prefix)

        if rust is None:
            print(f"{fname:<18} SKIP (no Rust output)")
            continue
        if c is None:
            print(f"{fname:<18} SKIP (no C output)")
            continue
        if rust.shape != c.shape:
            print(f"{fname:<18} SIZE MISMATCH rust={rust.shape} c={c.shape}")
            continue

        # All-pixel stats for % differ and center pixel
        all_diff = np.abs(rust.astype(np.int32) - c.astype(np.int32))
        n_nonzero = np.count_nonzero(all_diff)
        pct = 100.0 * n_nonzero / n_total
        center_diff = abs(int(rust[center_idx]) - int(c[center_idx]))

        # Valid-pixel stats for distribution
        vdiff = np.abs(rust[valid].astype(np.int32) - c[valid].astype(np.int32))
        nonzero_v = vdiff[vdiff > 0]

        if len(nonzero_v) == 0:
            print(f"{fname:<18} {'0.00%':>8} {'—':>7} {'—':>7} {'—':>7} {'—':>7} {'—':>7} {center_diff:>7}")
        else:
            median = int(np.median(nonzero_v))
            mean = nonzero_v.mean()
            p95 = int(np.percentile(nonzero_v, 95))
            p99 = int(np.percentile(nonzero_v, 99))
            mx = int(nonzero_v.max())
            print(f"{fname:<18} {pct:>7.2f}% {median:>7} {mean:>7.1f} {p95:>7} {p99:>7} {mx:>7} {center_diff:>7}")

    print()
    print("Stats are absolute differences on valid (non-fill) pixels.")
    print("Values are in scaled int16 units (÷10000 for reflectance or AOT).")


def main():
    parser = argparse.ArgumentParser(description="Compare Rust and C ESPA outputs")
    parser.add_argument("--rust-dir", default=DEFAULT_RUST_DIR,
                        help="Directory with Rust output .img files")
    parser.add_argument("--c-dir", default=DEFAULT_C_DIR,
                        help="Directory with C output .img files")
    parser.add_argument("--c-prefix", default=DEFAULT_PREFIX,
                        help="Filename prefix on C output files")
    args = parser.parse_args()

    compare(
        os.path.abspath(args.rust_dir),
        os.path.abspath(args.c_dir),
        args.c_prefix,
    )


if __name__ == "__main__":
    main()
