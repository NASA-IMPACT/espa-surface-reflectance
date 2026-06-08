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

SENTINEL_BANDS = [
    ("sr_band1.img", np.uint16),
    ("sr_band2.img", np.uint16),
    ("sr_band3.img", np.uint16),
    ("sr_band4.img", np.uint16),
    ("sr_band5.img", np.uint16),
    ("sr_band6.img", np.uint16),
    ("sr_band7.img", np.uint16),
    ("sr_band8.img", np.uint16),
    ("sr_band8a.img", np.uint16),
    ("sr_band9.img", np.uint16),
    ("sr_band10.img", np.uint16),
    ("sr_band11.img", np.uint16),
    ("sr_band12.img", np.uint16),
    ("sr_aerosol.img", np.int16),
    ("sr_aerosol_qa.img", np.uint8),
]

SENTINEL_NROWS, SENTINEL_NCOLS = 10980, 10980


def load_band(directory, fname, dtype, prefix=""):
    path = os.path.join(directory, prefix + fname)
    if not os.path.exists(path):
        return None
    return np.fromfile(path, dtype=dtype)


def compare(rust_dir, c_dir, c_prefix, bands=None, nrows=None, ncols=None):
    if bands is None:
        bands = BANDS
    if nrows is None:
        nrows = NROWS
    if ncols is None:
        ncols = NCOLS

    # Load aerosol to build valid-pixel mask (non-fill)
    c_aero = load_band(c_dir, "sr_aerosol.img", np.int16, c_prefix)
    rust_aero = load_band(rust_dir, "sr_aerosol.img", np.int16)
    if c_aero is None or rust_aero is None:
        print("ERROR: Could not load aerosol bands for valid-pixel mask.")
        sys.exit(1)

    valid = (c_aero != -9999) & (rust_aero != -9999)
    n_valid = valid.sum()
    n_total = len(c_aero)
    print(f"Scene: {nrows} x {ncols} = {n_total} pixels, {n_valid} valid ({100*n_valid/n_total:.1f}%)")
    print()

    header = f"| {'Band':<18} | {'% differ':>8} | {'Median':>7} | {'Mean':>7} | {'P95':>7} | {'P99':>7} | {'Max':>7} | {'Center':>7} |"
    sep = "|" + "-" * 20 + "|" + "-" * 10 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|"
    print(header)
    print(sep)

    center_idx = (nrows // 2) * ncols + (ncols // 2)

    for i, (fname, dtype) in enumerate(bands):
        # Print divider before aerosol bands
        if fname == "sr_aerosol.img":
            print(sep)

        rust = load_band(rust_dir, fname, dtype)
        c = load_band(c_dir, fname, dtype, c_prefix)

        if rust is None:
            print(f"| {fname:<18} | SKIP (no Rust output)")
            continue
        if c is None:
            print(f"| {fname:<18} | SKIP (no C output)")
            continue
        if rust.shape != c.shape:
            print(f"| {fname:<18} | SIZE MISMATCH rust={rust.shape} c={c.shape}")
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
            print(f"| {fname:<18} | {'0.00%':>8} | {'—':>7} | {'—':>7} | {'—':>7} | {'—':>7} | {'—':>7} | {center_diff:>7} |")
        else:
            median = int(np.median(nonzero_v))
            mean = nonzero_v.mean()
            p95 = int(np.percentile(nonzero_v, 95))
            p99 = int(np.percentile(nonzero_v, 99))
            mx = int(nonzero_v.max())
            print(f"| {fname:<18} | {pct:>7.2f}% | {median:>7} | {mean:>7.1f} | {p95:>7} | {p99:>7} | {mx:>7} | {center_diff:>7} |")

    print()
    print("Stats are absolute differences on valid (non-fill) pixels.")
    print("Values are in scaled int16 units (÷10000 for reflectance or AOT).")


def main():
    parser = argparse.ArgumentParser(description="Compare Rust and C ESPA outputs")
    parser.add_argument("--sensor", choices=["landsat", "sentinel"], default="landsat",
                        help="Sensor type: landsat (default) or sentinel")
    parser.add_argument("--rust-dir", default=None,
                        help="Directory with Rust output .img files")
    parser.add_argument("--c-dir", default=None,
                        help="Directory with C output .img files")
    parser.add_argument("--c-prefix", default=None,
                        help="Filename prefix on C output files")
    args = parser.parse_args()

    if args.sensor == "sentinel":
        rust_dir = args.rust_dir or os.path.join(
            os.path.dirname(__file__), "..", "test_data", "output_espa_s2"
        )
        c_dir = args.c_dir or os.path.join(
            os.path.dirname(__file__), "..", "test_data", "c_output"
        )
        c_prefix = args.c_prefix or "S2B_MSI_L1C_T38PNC_20260124_20260124_"
        bands = SENTINEL_BANDS
        nrows = SENTINEL_NROWS
        ncols = SENTINEL_NCOLS
    else:
        rust_dir = args.rust_dir or DEFAULT_RUST_DIR
        c_dir = args.c_dir or DEFAULT_C_DIR
        c_prefix = args.c_prefix or DEFAULT_PREFIX
        bands = BANDS
        nrows = NROWS
        ncols = NCOLS

    compare(
        os.path.abspath(rust_dir),
        os.path.abspath(c_dir),
        c_prefix,
        bands=bands,
        nrows=nrows,
        ncols=ncols,
    )


if __name__ == "__main__":
    main()
