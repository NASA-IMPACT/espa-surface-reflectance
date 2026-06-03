# LaSRC Rust Port — Next Steps

**Context:** The Landsat fast-path is functional with 45 passing unit tests. Remaining differences from C output are small (median 1-3 scaled int16 units for SR bands, median 2 for aerosol) but affect 55-65% of pixels due to f32/f64 precision cascading from aerosol retrieval.

**Current comparison results (branch `refactor`, commit `3cd9d5d`):**
- SR bands: 27-66% of pixels differ, median abs diff 1-3 (0.0001-0.0003 reflectance)
- Aerosol: 10.91% differ, median 2 but mean 20.1 (outliers from fill boundaries)
- QA: 0.04% differ

---

## Closing the Accuracy Gap

### 1. Match f32 precision in aerosol convergence loop

The ~11% aerosol pixel difference rate (and cascading ~55-65% SR differences) stems from `subaeroret_new` using f64 where C uses `float`. The convergence loop finds slightly different AOT optima at ~0.0002 AOT precision, which then produces different SR values via the polynomial fast path.

**Approach:** Cast the convergence arithmetic in `aerosol.rs` (`subaeroret_new`) to f32, matching the C code's float usage in `subaeroret.c`. This is the same pattern already applied to `atmcorlamb2_new`, `eval_cubic`, and `get_3rd_order_poly_coeff`.

**Files:** `lasrc_core/src/aerosol.rs`

**Expected impact:** Should reduce aerosol differences close to zero, which cascades to significantly reduce SR band differences.

### 2. Investigate fill-boundary artifacts

A small number of pixels show extreme diffs (65534/65535 on bands 1-5) from aerosol interpolation picking different fill-boundary neighbors. Band 7 (max diff=207) has no such artifacts.

**Approach:** Compare the aerosol fix/interpolation window logic at scene edges between C and Rust. Likely a difference in how fill neighbors are handled during the spatial interpolation step.

**Files:** `lasrc_core/src/aerosol.rs` (fix_aero / interp functions), `lasrc/c_version/src/compute_landsat_refl.c`

---

## Production Readiness

### 3. Add Rayon parallelism

Aerosol retrieval and final correction loops are single-threaded. C uses OpenMP. Per-pixel work is independent, making this straightforward.

**Files:** `lasrc_core/src/correction.rs` (main pixel loops), `lasrc_core/Cargo.toml` (add rayon dep)

### 4. Integration tests

Automate what `lasrc_py/compare_outputs.py` does manually: run the full pipeline against a known C output and assert all bands match within tolerance.

**Approach:** Either a Rust integration test that shells out to the Python pipeline, or a pytest that calls the pipeline and compares. Needs a small test scene committed or downloadable.

### 5. CI pipeline

Cargo test + clippy + format checks. Ideally includes a small end-to-end test scene.

---

## Scope Expansion

### 6. Sentinel-2 support

Several `todo!()` stubs remain: `tauray()`, `gas_coefficients()` on the Sentinel sensor, scene reader, different LUT paths. The trait-based sensor design was built to accommodate this.

**Files:** `lasrc_core/src/sensor.rs`, new `lasrc_py/python/lasrc/sentinel.py` reader, Sentinel-specific constants in `constants.rs`

### 7. `use_orig_aero` slow path

The full per-pixel `atmcorlamb2` atmospheric correction path (as opposed to the polynomial fast path). Not needed for current Landsat workflows but would complete the port.

**Files:** `lasrc_core/src/correction.rs`
