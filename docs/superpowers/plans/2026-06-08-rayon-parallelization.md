# Rayon Parallelization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Rayon parallelism to the three hottest loops in both Landsat and Sentinel-2 surface reflectance pipelines.

**Architecture:** Row-chunked `par_iter` / collect-and-scatter pattern using Rayon parallel iterators. No `unsafe` code. Optional thread count control via `ThreadPool`. Aerosol QA bits extracted to a serial post-pass to avoid shared mutable state.

**Tech Stack:** Rust, Rayon 1.10 (already in `lasrc_core/Cargo.toml`), PyO3, ndarray

---

## File Map

- **Modify:** `lasrc_core/src/correction.rs` — Add `num_threads` parameter to both entry points, parallelize 3 loops in each pipeline, extract QA post-passes
- **Modify:** `lasrc_py/src/lib.rs` — Thread `num_threads` through PyO3 bindings
- **No new files**

---

### Task 1: Thread Pool Plumbing

**Files:**
- Modify: `lasrc_core/src/correction.rs:400-413` (Landsat signature)
- Modify: `lasrc_core/src/correction.rs:991-1002` (Sentinel signature)
- Modify: `lasrc_py/src/lib.rs:151-169` (Landsat PyO3 binding)
- Modify: `lasrc_py/src/lib.rs:249-265` (Sentinel PyO3 binding)

- [ ] **Step 1: Add `use rayon` imports to correction.rs**

At the top of `lasrc_core/src/correction.rs`, add:

```rust
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
```

- [ ] **Step 2: Add `num_threads` parameter to `compute_surface_reflectance`**

Change the Landsat function signature from:

```rust
pub fn compute_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<f32>],
    bt_bands: &[ArrayView2<f32>],
    solar_zenith: &ArrayView2<f32>,
    solar_azimuth: &ArrayView2<f32>,
    view_zenith: &ArrayView2<f32>,
    view_azimuth: &ArrayView2<f32>,
    qa_band: &ArrayView2<u16>,
    lut: &LookupTables,
    aux: &AuxiliaryData,
    space_def: &SpaceDef,
    _use_orig_aero: bool,
) -> SurfaceReflectanceResult {
```

to:

```rust
pub fn compute_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<f32>],
    bt_bands: &[ArrayView2<f32>],
    solar_zenith: &ArrayView2<f32>,
    solar_azimuth: &ArrayView2<f32>,
    view_zenith: &ArrayView2<f32>,
    view_azimuth: &ArrayView2<f32>,
    qa_band: &ArrayView2<u16>,
    lut: &LookupTables,
    aux: &AuxiliaryData,
    space_def: &SpaceDef,
    _use_orig_aero: bool,
    num_threads: Option<usize>,
) -> SurfaceReflectanceResult {
```

At the very start of the function body (before Step 1), add the pool setup:

```rust
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0)) // 0 = Rayon default (all cores)
        .build()
        .expect("Failed to build Rayon thread pool");
```

- [ ] **Step 3: Add `num_threads` parameter to `compute_sentinel_surface_reflectance`**

Change the Sentinel function signature from:

```rust
pub fn compute_sentinel_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<f32>],
    solar_zenith: f64,
    solar_azimuth: f64,
    view_zenith: f64,
    view_azimuth: f64,
    qa_band: &ArrayView2<u16>,
    lut: &LookupTables,
    aux: &AuxiliaryData,
    space_def: &SpaceDef,
) -> SurfaceReflectanceResult {
```

to:

```rust
pub fn compute_sentinel_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<f32>],
    solar_zenith: f64,
    solar_azimuth: f64,
    view_zenith: f64,
    view_azimuth: f64,
    qa_band: &ArrayView2<u16>,
    lut: &LookupTables,
    aux: &AuxiliaryData,
    space_def: &SpaceDef,
    num_threads: Option<usize>,
) -> SurfaceReflectanceResult {
```

Add the same pool setup at the start of the function body:

```rust
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0))
        .build()
        .expect("Failed to build Rayon thread pool");
```

- [ ] **Step 4: Update PyO3 Landsat binding**

In `lasrc_py/src/lib.rs`, add `num_threads: Option<usize>` to `process_surface_reflectance`:

```rust
fn process_surface_reflectance<'py>(
    py: Python<'py>,
    sensor_name: &str,
    toa_bands: Vec<PyReadonlyArray2<'py, f32>>,
    bt_bands: Vec<PyReadonlyArray2<'py, f32>>,
    solar_zenith: PyReadonlyArray2<'py, f32>,
    solar_azimuth: PyReadonlyArray2<'py, f32>,
    view_zenith: PyReadonlyArray2<'py, f32>,
    view_azimuth: PyReadonlyArray2<'py, f32>,
    qa_band: PyReadonlyArray2<'py, u16>,
    lut: &PyLookupTables,
    aux: &PyAuxiliaryData,
    ul_corner_x: f64,
    ul_corner_y: f64,
    pixel_size_x: f64,
    pixel_size_y: f64,
    utm_zone: i32,
    use_orig_aero: bool,
    num_threads: Option<usize>,
) -> PyResult<PyObject> {
```

And pass it through to the Rust call:

```rust
    let result = compute_surface_reflectance(
        sensor.as_ref(),
        &toa_views,
        &bt_views,
        &solar_zenith.as_array(),
        &solar_azimuth.as_array(),
        &view_zenith.as_array(),
        &view_azimuth.as_array(),
        &qa_band.as_array(),
        &lut.inner,
        &aux.inner,
        &space_def,
        use_orig_aero,
        num_threads,
    );
```

- [ ] **Step 5: Update PyO3 Sentinel binding**

In `lasrc_py/src/lib.rs`, add `num_threads: Option<usize>` to `process_sentinel_surface_reflectance`:

```rust
fn process_sentinel_surface_reflectance<'py>(
    py: Python<'py>,
    sensor_name: &str,
    toa_bands: Vec<PyReadonlyArray2<'py, f32>>,
    solar_zenith: f64,
    solar_azimuth: f64,
    view_zenith: f64,
    view_azimuth: f64,
    qa_band: PyReadonlyArray2<'py, u16>,
    lut: &PyLookupTables,
    aux: &PyAuxiliaryData,
    ul_corner_x: f64,
    ul_corner_y: f64,
    pixel_size_x: f64,
    pixel_size_y: f64,
    utm_zone: i32,
    num_threads: Option<usize>,
) -> PyResult<PyObject> {
```

And pass it through:

```rust
    let result = compute_sentinel_surface_reflectance(
        sensor.as_ref(),
        &toa_views,
        solar_zenith,
        solar_azimuth,
        view_zenith,
        view_azimuth,
        &qa_band.as_array(),
        &lut.inner,
        &aux.inner,
        &space_def,
        num_threads,
    );
```

- [ ] **Step 6: Verify it compiles**

Run: `cargo test -p lasrc_core --no-run`
Expected: Compiles with no errors. Warnings about unused `pool` variable are expected at this stage.

- [ ] **Step 7: Commit**

```bash
git add lasrc_core/src/correction.rs lasrc_py/src/lib.rs
git commit -m "feat: add num_threads parameter for Rayon parallelization"
```

---

### Task 2: Parallelize Aerosol Window-Center Retrieval (Landsat)

**Files:**
- Modify: `lasrc_core/src/correction.rs:543-809` (Landsat Step 5)

This is the highest-impact loop. Each window center runs 3-4 `subaeroret_new` calls independently. The strategy is collect-and-scatter: pre-collect window center coordinates, `par_iter` to compute results, scatter back serially.

- [ ] **Step 1: Define the result struct**

Above the Landsat `compute_surface_reflectance` function (or in a `mod parallel` block at the top of the file), add:

```rust
/// Result of aerosol retrieval at a single window center.
struct AerosolWindowResult {
    pix: usize,
    taero: f32,
    teps: f32,
    ipflag: u8,
}
```

- [ ] **Step 2: Pre-collect window center coordinates**

Replace the current Landsat Step 5 opening (the `while` loops at ~line 544):

```rust
    // Iterate over window centers
    let mut iline = half_aero_window;
    while iline < nlines {
        let mut isamp = half_aero_window;
        while isamp < nsamps {
```

With coordinate collection:

```rust
    // Collect window-center coordinates for parallel processing
    let window_centers: Vec<(usize, usize)> = {
        let mut centers = Vec::new();
        let mut iline = half_aero_window;
        while iline < nlines {
            let mut isamp = half_aero_window;
            while isamp < nsamps {
                centers.push((iline, isamp));
                isamp += aero_window;
            }
            iline += aero_window;
        }
        centers
    };
```

- [ ] **Step 3: Replace the loop body with a parallel map**

Replace the entire Landsat aerosol retrieval loop (from the `while iline` through the closing `iline += aero_window; }`) with:

```rust
    let aero_results: Vec<AerosolWindowResult> = pool.install(|| {
        window_centers.par_iter().map(|&(iline, isamp)| {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels
            if is_fill_pixel(qa_flat[pix]) {
                return AerosolWindowResult {
                    pix,
                    taero: DEFAULT_AERO as f32,
                    teps: DEFAULT_EPS as f32,
                    ipflag: 1u8 << IPFLAG_FILL,
                };
            }

            let xmus_pixel = (solar_zenith[(iline, isamp)] as f64 * DEG2RAD).cos();
            let toa_over_cos = |band_idx: usize| -> f64 {
                let raw_toa = toa_bands[band_idx][(iline, isamp)] as f64;
                (raw_toa / xmus_pixel).clamp(MIN_VALID_REFL, MAX_VALID_REFL)
            };

            // Compute per-pixel lat/lon from image coordinates
            let (pixel_lat, pixel_lon) = utm_to_deg(space_def, iline as i32, isamp as i32);

            // Look up CMG position for slope/intercept computation
            let cmg = latlon_to_cmg(pixel_lat, pixel_lon);
            let ratio_pix11 = cmg.idx[0];
            let ratio_pix12 = cmg.idx[1];
            let ratio_pix21 = cmg.idx[2];
            let ratio_pix22 = cmg.idx[3];
            let corner_indices = [ratio_pix11, ratio_pix12, ratio_pix21, ratio_pix22];

            let safe_read = |arr: &[i16], idx: usize| -> i16 {
                if idx < arr.len() { arr[idx] } else { 0 }
            };

            // Compute modified slopes and intercepts for each corner and band
            let mut slp_b1 = [0.0f64; 4];
            let mut int_b1 = [0.0f64; 4];
            let mut slp_b2 = [0.0f64; 4];
            let mut int_b2 = [0.0f64; 4];
            let mut slp_b7 = [0.0f64; 4];
            let mut int_b7 = [0.0f64; 4];

            for (ci, &cidx) in corner_indices.iter().enumerate() {
                let rb1_val = safe_read(&aux.ratiob1, cidx);
                let rb2_val = safe_read(&aux.ratiob2, cidx);
                let sndwi_val = safe_read(&aux.sndwi, cidx);

                let (s, i) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob1, cidx),
                    safe_read(&aux.intratiob1, cidx),
                    rb1_val, 550,
                );
                slp_b1[ci] = s;
                int_b1[ci] = i;

                let (s, i) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob2, cidx),
                    safe_read(&aux.intratiob2, cidx),
                    rb2_val, 600,
                );
                slp_b2[ci] = s;
                int_b2[ci] = i;

                let (s, i) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob7, cidx),
                    safe_read(&aux.intratiob7, cidx),
                    safe_read(&aux.ratiob7, cidx), 2000,
                );
                slp_b7[ci] = s;
                int_b7[ci] = i;
            }

            let slprb1 = bilerp(slp_b1, &cmg.w);
            let intrb1 = bilerp(int_b1, &cmg.w);
            let slprb2 = bilerp(slp_b2, &cmg.w);
            let intrb2 = bilerp(int_b2, &cmg.w);
            let slprb7 = bilerp(slp_b7, &cmg.w);
            let intrb7 = bilerp(int_b7, &cmg.w);

            let sr_nir = sband[bi.nir][pix] as f64;
            let sr_swir2 = sband[bi.swir2][pix] as f64;
            let sr_swir2_half = sr_swir2 * 0.5;
            let denom = sr_nir + sr_swir2_half;
            let mut xndwi = if denom.abs() > 1.0e-10 {
                (sr_nir - sr_swir2_half) / denom
            } else {
                0.0
            };

            let andwi_val = safe_read(&aux.andwi, ratio_pix11);
            let sndwi_val = safe_read(&aux.sndwi, ratio_pix11);
            let ndwi_th1 = (andwi_val as f64 + 2.0 * sndwi_val as f64) * 0.001;
            let ndwi_th2 = (andwi_val as f64 - 2.0 * sndwi_val as f64) * 0.001;
            if xndwi > ndwi_th1 { xndwi = ndwi_th1; }
            if xndwi < ndwi_th2 { xndwi = ndwi_th2; }

            let mut erelc = vec![-1.0f64; nbands];
            let mut troatm = vec![0.0f64; nbands];

            erelc[bi.coastal] = (xndwi * slprb1 + intrb1) as f32 as f64;
            erelc[bi.blue] = (xndwi * slprb2 + intrb2) as f32 as f64;
            erelc[bi.red] = 1.0;
            erelc[bi.swir2] = (xndwi * slprb7 + intrb7) as f32 as f64;

            troatm[bi.coastal] = toa_over_cos(bi.coastal) as f32 as f64;
            troatm[bi.blue] = toa_over_cos(bi.blue) as f32 as f64;
            troatm[bi.red] = toa_over_cos(bi.red) as f32 as f64;
            troatm[bi.swir2] = toa_over_cos(bi.swir2) as f32 as f64;

            // === Eps optimization: 3 retrievals ===
            let mut iaots = 0usize;
            let result1 = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps1, iaots, tth,
            );
            let residual1 = result1.residual;
            let sraot1 = result1.raot;
            iaots = result1.iaots;

            let result2 = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps2, iaots, tth,
            );
            let residual2 = result2.residual;
            iaots = result2.iaots;

            let result3 = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps3, iaots, tth,
            );
            let residual3 = result3.residual;
            let sraot3 = result3.raot;
            iaots = result3.iaots;

            let xc = residual1 - residual3;
            let xf = residual2 - residual3;
            let denom_fit = xa * xe - xb * xd;
            let coefa = (xc * xe - xb * xf) / denom_fit;
            let coefb = (xa * xf - xc * xd) / denom_fit;
            let epsmin = -coefb / (2.0 * coefa);

            let (eps, raot, residual) = if epsmin >= LOW_EPS && epsmin <= HIGH_EPS {
                let result_opt = subaeroret_new(
                    false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                    &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                    lambda, epsmin, iaots, tth,
                );
                (epsmin, result_opt.raot, result_opt.residual)
            } else if epsmin <= LOW_EPS {
                (eps1, sraot1, residual1)
            } else {
                (eps3, sraot3, residual3)
            };

            let corf = raot / xmus_center;
            let mut result_ipflag = 0u8;

            if residual < (0.015 + 0.005 * corf + 0.10 * troatm[bi.swir2]) {
                let ros5 = atmcorlamb2_new(
                    &atm_coeff[bi.nir], tgo_arr[bi.nir], bi.nir,
                    raot, normext_p0a3[bi.nir],
                    toa_over_cos(bi.nir), lambda, eps,
                );
                let ros4 = atmcorlamb2_new(
                    &atm_coeff[bi.red], tgo_arr[bi.red], bi.red,
                    raot, normext_p0a3[bi.red],
                    toa_over_cos(bi.red), lambda, eps,
                );

                if ros5 > 0.1 && (ros5 - ros4) / (ros5 + ros4) > 0.0 {
                    result_ipflag |= 1u8 << IPFLAG_CLEAR;
                } else {
                    result_ipflag |= 1u8 << IPFLAG_WATER;
                }
            } else {
                result_ipflag |= 1u8 << IPFLAG_WATER;
            }

            // === Water retest ===
            if result_ipflag & (1u8 << IPFLAG_WATER) != 0 {
                let mut water_erelc = vec![-1.0f64; nbands];
                water_erelc[bi.coastal] = 1.0;
                water_erelc[bi.red] = 1.0;
                water_erelc[bi.nir] = 1.0;
                water_erelc[bi.swir2] = 1.0;

                let mut water_troatm = vec![0.0f64; nbands];
                water_troatm[bi.coastal] = toa_over_cos(bi.coastal) as f32 as f64;
                water_troatm[bi.red] = toa_over_cos(bi.red) as f32 as f64;
                water_troatm[bi.nir] = toa_over_cos(bi.nir) as f32 as f64;
                water_troatm[bi.swir2] = toa_over_cos(bi.swir2) as f32 as f64;

                let water_result = subaeroret_new(
                    true, iband1, &water_erelc, &water_troatm,
                    &tgo_arr, &roatm_ia_max, &roatm_coef, &ttatmg_coef,
                    &satm_coef, &normext_p0a3, lambda, WATER_EPS, 0, tth_water,
                );

                let water_corf = water_result.raot / xmus_center;
                let ros1 = atmcorlamb2_new(
                    &atm_coeff[bi.coastal], tgo_arr[bi.coastal], bi.coastal,
                    water_result.raot, normext_p0a3[bi.coastal],
                    toa_over_cos(bi.coastal), lambda, WATER_EPS,
                );

                if water_result.residual > (0.010 + 0.005 * water_corf) || ros1 < 0.0 {
                    result_ipflag = 0;
                } else {
                    result_ipflag = (1u8 << IPFLAG_CLEAR) | (1u8 << IPFLAG_WATER);
                }

                return AerosolWindowResult {
                    pix,
                    taero: water_result.raot as f32,
                    teps: WATER_EPS as f32,
                    ipflag: result_ipflag,
                };
            }

            AerosolWindowResult {
                pix,
                taero: raot as f32,
                teps: eps as f32,
                ipflag: result_ipflag,
            }
        }).collect()
    });

    // Scatter results back into working arrays
    for r in &aero_results {
        taero[r.pix] = r.taero;
        teps[r.pix] = r.teps;
        ipflag[r.pix] = r.ipflag;
    }
```

- [ ] **Step 4: Verify it compiles**

Run: `cargo test -p lasrc_core --no-run`
Expected: Compiles. May have warnings about unused variables that were previously used in the sequential loop.

- [ ] **Step 5: Run unit tests**

Run: `cargo test -p lasrc_core`
Expected: All tests pass.

- [ ] **Step 6: Build Python bindings and run Landsat comparison**

```bash
source .venv/bin/activate
maturin develop --release -m lasrc_py/Cargo.toml
python lasrc_py/run_test_granule.py --format espa
python lasrc_py/compare_outputs.py --sensor landsat
```

Expected: Results within existing thresholds (median 1-2 for SR bands, aerosol <7% differ).

- [ ] **Step 7: Commit**

```bash
git add lasrc_core/src/correction.rs
git commit -m "feat: parallelize Landsat aerosol window-center retrieval with Rayon"
```

---

### Task 3: Parallelize Aerosol Window-Center Retrieval (Sentinel)

**Files:**
- Modify: `lasrc_core/src/correction.rs:1135-1470` (Sentinel Step 7)

Same collect-and-scatter pattern as Landsat, but with Sentinel-specific differences: window starts at (0,0) not half_aero_window, TOA averaging within 6x6 windows, different eps optimization, and results are broadcast to all pixels in the 6x6 window.

- [ ] **Step 1: Pre-collect Sentinel window center coordinates**

Replace the Sentinel Step 7 `while` loop opening with:

```rust
    // Collect Sentinel window-center coordinates for parallel processing
    let window_centers: Vec<(usize, usize)> = {
        let mut centers = Vec::new();
        let mut win_i = 0;
        while win_i < nlines {
            let mut win_j = 0;
            while win_j < nsamps {
                centers.push((win_i, win_j));
                win_j += aero_window;
            }
            win_i += aero_window;
        }
        centers
    };
```

- [ ] **Step 2: Define the Sentinel aerosol result struct**

Add near the `AerosolWindowResult` struct from Task 2:

```rust
/// Result of Sentinel aerosol retrieval at a single window center.
/// Includes the window bounds for broadcasting to all pixels in the window.
struct SentinelAerosolWindowResult {
    win_i: usize,
    win_j: usize,
    curr_pix: usize,
    taero: f32,
    teps: f32,
    ipflag: u8,
    is_fill: bool,
}
```

- [ ] **Step 3: Replace the Sentinel aerosol loop with parallel map**

Replace the entire Sentinel aerosol retrieval loop (from `let mut win_i = 0;` through the closing `win_i += aero_window; }` including the window broadcast at lines 1457-1465) with:

```rust
    let aero_results: Vec<SentinelAerosolWindowResult> = pool.install(|| {
        window_centers.par_iter().map(|&(win_i, win_j)| {
            let curr_pix = win_i * nsamps + win_j;

            // Skip fill pixels
            if is_fill_pixel(qaband[curr_pix]) {
                return SentinelAerosolWindowResult {
                    win_i, win_j, curr_pix,
                    taero: DEFAULT_AERO as f32,
                    teps: DEFAULT_EPS as f32,
                    ipflag: 1u8 << IPFLAG_FILL,
                    is_fill: true,
                };
            }

            // (Paste the entire Sentinel window-center body here, from
            //  utm_to_deg through the water retest, but using local
            //  result_ipflag instead of ipflag[curr_pix], same pattern
            //  as the Landsat version in Task 2 Step 3.)
            //
            // Key differences from Landsat:
            // 1. Uses bilerp_f32 instead of bilerp
            // 2. TOA averaging in 6x6 windows (not single-pixel)
            // 3. Sentinel eps optimization with resepsmin check
            // 4. Uses DNS_BAND* constants instead of bi.*
            // 5. xndwi uses f32 truncation
            // 6. Water retest uses DNS_BAND1/4/8A/12 with window-averaged TOA

            let (pixel_lat, pixel_lon) = utm_to_deg(space_def, win_i as i32, win_j as i32);
            let cmg = latlon_to_cmg(pixel_lat, pixel_lon);
            let ratio_pix11 = cmg.idx[0];
            let ratio_pix12 = cmg.idx[1];
            let ratio_pix21 = cmg.idx[2];
            let ratio_pix22 = cmg.idx[3];
            let corner_indices = [ratio_pix11, ratio_pix12, ratio_pix21, ratio_pix22];

            let safe_read = |arr: &[i16], idx: usize| -> i16 {
                if idx < arr.len() { arr[idx] } else { 0 }
            };

            let mut slp_b1 = [0.0f64; 4];
            let mut int_b1 = [0.0f64; 4];
            let mut slp_b2 = [0.0f64; 4];
            let mut int_b2 = [0.0f64; 4];
            let mut slp_b7 = [0.0f64; 4];
            let mut int_b7 = [0.0f64; 4];

            for (ci, &cidx) in corner_indices.iter().enumerate() {
                let rb1_val = safe_read(&aux.ratiob1, cidx);
                let rb2_val = safe_read(&aux.ratiob2, cidx);
                let sndwi_val = safe_read(&aux.sndwi, cidx);

                let (s, int) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob1, cidx),
                    safe_read(&aux.intratiob1, cidx),
                    rb1_val, 550,
                );
                slp_b1[ci] = s;
                int_b1[ci] = int;

                let (s, int) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob2, cidx),
                    safe_read(&aux.intratiob2, cidx),
                    rb2_val, 600,
                );
                slp_b2[ci] = s;
                int_b2[ci] = int;

                let (s, int) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob7, cidx),
                    safe_read(&aux.intratiob7, cidx),
                    safe_read(&aux.ratiob7, cidx), 2000,
                );
                slp_b7[ci] = s;
                int_b7[ci] = int;
            }

            let slprb1 = bilerp_f32(slp_b1, &cmg.w);
            let intrb1 = bilerp_f32(int_b1, &cmg.w);
            let slprb2 = bilerp_f32(slp_b2, &cmg.w);
            let intrb2 = bilerp_f32(int_b2, &cmg.w);
            let slprb7 = bilerp_f32(slp_b7, &cmg.w);
            let intrb7 = bilerp_f32(int_b7, &cmg.w);

            let sr_nir = sband[DNS_BAND8A][curr_pix] as f64;
            let sr_swir2 = sband[DNS_BAND12][curr_pix] as f64;
            let sr_swir2_half = sr_swir2 * 0.5;
            let denom = sr_nir + sr_swir2_half;
            let mut xndwi: f32 = if denom.abs() > 1.0e-10 {
                ((sr_nir - sr_swir2_half) / denom) as f32
            } else {
                0.0f32
            };

            let andwi_val = safe_read(&aux.andwi, ratio_pix11);
            let sndwi_val = safe_read(&aux.sndwi, ratio_pix11);
            let ndwi_th1 = ((andwi_val as f64 + 2.0 * sndwi_val as f64) * 0.001) as f32;
            let ndwi_th2 = ((andwi_val as f64 - 2.0 * sndwi_val as f64) * 0.001) as f32;
            if xndwi > ndwi_th1 { xndwi = ndwi_th1; }
            if xndwi < ndwi_th2 { xndwi = ndwi_th2; }

            let mut erelc = vec![-1.0f64; nbands];
            let mut troatm = vec![0.0f64; nbands];

            erelc[DNS_BAND1] = (xndwi * slprb1 + intrb1) as f32 as f64;
            erelc[DNS_BAND2] = (xndwi * slprb2 + intrb2) as f32 as f64;
            erelc[DNS_BAND4] = 1.0;
            erelc[DNS_BAND12] = (xndwi * slprb7 + intrb7) as f32 as f64;

            // Average TOA in 6x6 window
            let mut pix_count = 0u32;
            let ew_line = (win_i + aero_window).min(nlines);
            let ew_samp = (win_j + aero_window).min(nsamps);
            let mut troatm_b1_f32 = 0.0f32;
            let mut troatm_b2_f32 = 0.0f32;
            let mut troatm_b4_f32 = 0.0f32;
            let mut troatm_b12_f32 = 0.0f32;
            for iline in win_i..ew_line {
                for isamp in win_j..ew_samp {
                    let win_pix = iline * nsamps + isamp;
                    if is_fill_pixel(qaband[win_pix]) { continue; }
                    troatm_b1_f32 += toa_bands[DNS_BAND1][(iline, isamp)];
                    troatm_b2_f32 += toa_bands[DNS_BAND2][(iline, isamp)];
                    troatm_b4_f32 += toa_bands[DNS_BAND4][(iline, isamp)];
                    troatm_b12_f32 += toa_bands[DNS_BAND12][(iline, isamp)];
                    pix_count += 1;
                }
            }
            if pix_count > 0 {
                troatm_b1_f32 /= pix_count as f32;
                troatm_b2_f32 /= pix_count as f32;
                troatm_b4_f32 /= pix_count as f32;
                troatm_b12_f32 /= pix_count as f32;
            }
            troatm[DNS_BAND1] = troatm_b1_f32 as f64;
            troatm[DNS_BAND2] = troatm_b2_f32 as f64;
            troatm[DNS_BAND4] = troatm_b4_f32 as f64;
            troatm[DNS_BAND12] = troatm_b12_f32 as f64;

            // Eps optimization
            let mut iaots = 0usize;
            let result1 = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps1, iaots, tth,
            );
            let residual1 = result1.residual;
            iaots = result1.iaots;

            let result2 = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps2, iaots, tth,
            );
            let residual2 = result2.residual;
            iaots = result2.iaots;

            let result3 = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps3, iaots, tth,
            );
            let residual3 = result3.residual;
            iaots = result3.iaots;

            let xc = residual1 - residual3;
            let xf_val = residual2 - residual3;
            let denom_fit = xa * xe - xb * xd;
            let coefa = (xc * xe - xb * xf_val) / denom_fit;
            let coefb = (xa * xf_val - xc * xd) / denom_fit;
            let epsmin = -coefb / (2.0 * coefa);
            let resepsmin = xa * epsmin * epsmin + xb * epsmin + xc;

            let eps = if epsmin < LOW_EPS || epsmin > HIGH_EPS {
                if residual1 < residual3 { eps1 } else { eps3 }
            } else if resepsmin > residual1 || resepsmin > residual3 {
                if residual1 < residual3 { eps1 } else { eps3 }
            } else {
                epsmin
            };

            let result_final = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps, iaots, tth,
            );
            let raot = result_final.raot;
            let residual = result_final.residual;
            let corf = raot / xmus;

            let mut result_ipflag = 0u8;

            if residual < (0.015 + 0.005 * corf + 0.10 * troatm[DNS_BAND12]) {
                let mut rotoa_b8a = 0.0f64;
                let mut pc = 0u32;
                for iline in win_i..ew_line {
                    for isamp in win_j..ew_samp {
                        if isamp >= nsamps { continue; }
                        rotoa_b8a += toa_bands[DNS_BAND8A][(iline, isamp)] as f64;
                        pc += 1;
                    }
                }
                if pc > 0 { rotoa_b8a /= pc as f64; }

                let ros5 = atmcorlamb2_new(
                    &atm_coeff[DNS_BAND8A], tgo_arr[DNS_BAND8A], DNS_BAND8A,
                    raot, normext_p0a3[DNS_BAND8A],
                    rotoa_b8a, lambda, eps,
                );

                let mut rotoa_b4 = 0.0f64;
                pc = 0;
                for iline in win_i..ew_line {
                    for isamp in win_j..ew_samp {
                        if isamp >= nsamps { continue; }
                        rotoa_b4 += toa_bands[DNS_BAND4][(iline, isamp)] as f64;
                        pc += 1;
                    }
                }
                if pc > 0 { rotoa_b4 /= pc as f64; }

                let ros4 = atmcorlamb2_new(
                    &atm_coeff[DNS_BAND4], tgo_arr[DNS_BAND4], DNS_BAND4,
                    raot, normext_p0a3[DNS_BAND4],
                    rotoa_b4, lambda, eps,
                );

                if ros5 > 0.1 && (ros5 - ros4) / (ros5 + ros4) > 0.0 {
                    result_ipflag |= 1u8 << IPFLAG_CLEAR;
                } else {
                    result_ipflag = 1u8 << IPFLAG_WATER;
                }
            } else {
                result_ipflag = 1u8 << IPFLAG_WATER;
            }

            // Water retest
            if result_ipflag & (1u8 << IPFLAG_WATER) != 0 {
                let mut water_erelc = vec![-1.0f64; nbands];
                let mut water_troatm = vec![0.0f64; nbands];

                let mut water_pc = 0u32;
                let mut wt_b1 = 0.0f32;
                let mut wt_b4 = 0.0f32;
                let mut wt_b8a = 0.0f32;
                let mut wt_b12 = 0.0f32;
                for iline in win_i..ew_line {
                    for isamp in win_j..ew_samp {
                        let win_pix = iline * nsamps + isamp;
                        if is_fill_pixel(qaband[win_pix]) { continue; }
                        wt_b1 += toa_bands[DNS_BAND1][(iline, isamp)];
                        wt_b4 += toa_bands[DNS_BAND4][(iline, isamp)];
                        wt_b8a += toa_bands[DNS_BAND8A][(iline, isamp)];
                        wt_b12 += toa_bands[DNS_BAND12][(iline, isamp)];
                        water_pc += 1;
                    }
                }
                if water_pc > 0 {
                    wt_b1 /= water_pc as f32;
                    wt_b4 /= water_pc as f32;
                    wt_b8a /= water_pc as f32;
                    wt_b12 /= water_pc as f32;
                }
                water_troatm[DNS_BAND1] = wt_b1 as f64;
                water_troatm[DNS_BAND4] = wt_b4 as f64;
                water_troatm[DNS_BAND8A] = wt_b8a as f64;
                water_troatm[DNS_BAND12] = wt_b12 as f64;

                water_erelc[DNS_BAND1] = 1.0;
                water_erelc[DNS_BAND4] = 1.0;
                water_erelc[DNS_BAND8A] = 1.0;
                water_erelc[DNS_BAND12] = 1.0;

                let tth_water = &SENTINEL_TTH_WATER[..nbands];
                let water_result = subaeroret_new(
                    true, iband1, &water_erelc, &water_troatm,
                    &tgo_arr, &roatm_ia_max, &roatm_coef, &ttatmg_coef,
                    &satm_coef, &normext_p0a3, lambda, WATER_EPS, 0, tth_water,
                );
                let water_corf = water_result.raot / xmus;

                let ros1 = atmcorlamb2_new(
                    &atm_coeff[DNS_BAND1], tgo_arr[DNS_BAND1], DNS_BAND1,
                    water_result.raot, normext_p0a3[DNS_BAND1],
                    water_troatm[DNS_BAND1], lambda, WATER_EPS,
                );

                if water_result.residual > (0.010 + 0.005 * water_corf) || ros1 < 0.0 {
                    result_ipflag = 1u8 << IPFLAG_FAILED;
                } else {
                    result_ipflag = (1u8 << IPFLAG_WATER) | (1u8 << IPFLAG_CLEAR);
                }

                return SentinelAerosolWindowResult {
                    win_i, win_j, curr_pix,
                    taero: water_result.raot as f32,
                    teps: WATER_EPS as f32,
                    ipflag: result_ipflag,
                    is_fill: false,
                };
            }

            SentinelAerosolWindowResult {
                win_i, win_j, curr_pix,
                taero: raot as f32,
                teps: eps as f32,
                ipflag: result_ipflag,
                is_fill: false,
            }
        }).collect()
    });

    // Scatter results: for Sentinel, broadcast taero/teps to all non-fill pixels in 6x6 window
    for r in &aero_results {
        if r.is_fill {
            ipflag[r.curr_pix] = r.ipflag;
            continue;
        }
        taero[r.curr_pix] = r.taero;
        teps[r.curr_pix] = r.teps;
        ipflag[r.curr_pix] = r.ipflag;

        // Broadcast to all non-fill pixels in the window
        let ew_line = (r.win_i + aero_window).min(nlines);
        let ew_samp = (r.win_j + aero_window).min(nsamps);
        for iline in r.win_i..ew_line {
            for isamp in r.win_j..ew_samp {
                let win_pix = iline * nsamps + isamp;
                if is_fill_pixel(qaband[win_pix]) { continue; }
                teps[win_pix] = r.teps;
                taero[win_pix] = r.taero;
            }
        }
    }
```

- [ ] **Step 4: Verify it compiles**

Run: `cargo test -p lasrc_core --no-run`
Expected: Compiles.

- [ ] **Step 5: Build Python bindings and run Sentinel comparison**

```bash
source .venv/bin/activate
maturin develop --release -m lasrc_py/Cargo.toml
python lasrc_py/run_test_sentinel.py
python lasrc_py/compare_outputs.py --sensor sentinel
```

Expected: Results within existing thresholds (aerosol <7% differ, SR bands within current diffs).

- [ ] **Step 6: Commit**

```bash
git add lasrc_core/src/correction.rs
git commit -m "feat: parallelize Sentinel aerosol window-center retrieval with Rayon"
```

---

### Task 4: Parallelize Final Per-Pixel Atmospheric Correction (Both Pipelines)

**Files:**
- Modify: `lasrc_core/src/correction.rs:853-910` (Landsat Step 7)
- Modify: `lasrc_core/src/correction.rs:1477-1535` (Sentinel Step 9)

Split into: parallel SR computation over rows + serial QA post-pass.

- [ ] **Step 1: Parallelize Landsat final correction**

Replace the Landsat Step 7 loop (lines ~853-910) with:

```rust
    // ── Step 7: Final per-pixel atmospheric correction ──
    let mut sr_f32: Vec<Array2<f64>> = (0..nbands)
        .map(|_| Array2::zeros((nlines, nsamps)))
        .collect();

    // Parallel SR computation: process rows in chunks
    // Each row's SR depends only on immutable taero, teps, sband, atm_coeff
    pool.install(|| {
        (0..nlines).into_par_iter().for_each(|iline| {
            for isamp in 0..nsamps {
                let pix = iline * nsamps + isamp;
                if is_fill_pixel(qa_flat[pix]) { continue; }

                let raot = taero[pix] as f64;
                let eps = teps[pix] as f64;

                for iband in 0..nbands {
                    let rsurf = sband[iband][pix];
                    let rotoa: f32 = ((rsurf as f64 * bttatmg[iband]
                        / (1.0 - bsatm[iband] * rsurf as f64)
                        + broatm[iband])
                        * btgo[iband]) as f32;

                    let roslamb = atmcorlamb2_new(
                        &atm_coeff[iband],
                        tgo_arr[iband],
                        iband,
                        raot,
                        normext_p0a3[iband],
                        rotoa as f64,
                        lambda,
                        eps,
                    );

                    let roslamb = roslamb.clamp(MIN_VALID_REFL, MAX_VALID_REFL);

                    // SAFETY: Each row writes to its own row in sr_f32[band],
                    // but Array2 doesn't support split mutable access.
                    // Use unsafe raw pointer to write to non-overlapping rows.
                    unsafe {
                        let ptr = sr_f32[iband].as_mut_ptr();
                        *ptr.add(iline * nsamps + isamp) = roslamb;
                    }
                }
            }
        });
    });

    // Serial post-pass: set aerosol QA bits from coastal band differences
    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;
            if is_fill_pixel(qa_flat[pix]) { continue; }

            let rsurf = sband[bi.coastal][pix] as f64;
            let roslamb = sr_f32[bi.coastal][(iline, isamp)];
            let tmpf = (rsurf - roslamb).abs();
            if tmpf <= LOW_AERO_THRESH {
                ipflag[pix] |= 1u8 << AERO1_QA;
            } else if tmpf < AVG_AERO_THRESH {
                ipflag[pix] |= 1u8 << AERO2_QA;
            } else {
                ipflag[pix] |= (1u8 << AERO1_QA) | (1u8 << AERO2_QA);
            }
        }
    }
```

**Note on `unsafe`:** The `Array2` type does not provide `par_chunks_mut` for rows. We use raw pointer arithmetic to write to non-overlapping row positions from different threads. Each thread writes only to its own row, so there is no data race. If the implementer finds a way to restructure `sr_f32` as `Vec<Vec<f64>>` instead of `Vec<Array2<f64>>` to avoid `unsafe`, that is preferred — `Vec<f64>` can be split into `&mut [f64]` chunks cleanly.

**Alternative without unsafe:** Restructure `sr_f32` as `Vec<Vec<f64>>` (bands × flat pixels) instead of `Vec<Array2<f64>>`:

```rust
    let mut sr_flat: Vec<Vec<f64>> = (0..nbands)
        .map(|_| vec![0.0f64; npix])
        .collect();

    pool.install(|| {
        (0..nlines).into_par_iter().for_each(|iline| {
            for isamp in 0..nsamps {
                let pix = iline * nsamps + isamp;
                if is_fill_pixel(qa_flat[pix]) { continue; }

                let raot = taero[pix] as f64;
                let eps = teps[pix] as f64;

                for iband in 0..nbands {
                    let rsurf = sband[iband][pix];
                    let rotoa: f32 = ((rsurf as f64 * bttatmg[iband]
                        / (1.0 - bsatm[iband] * rsurf as f64)
                        + broatm[iband])
                        * btgo[iband]) as f32;

                    let roslamb = atmcorlamb2_new(
                        &atm_coeff[iband], tgo_arr[iband], iband,
                        raot, normext_p0a3[iband],
                        rotoa as f64, lambda, eps,
                    );

                    let roslamb = roslamb.clamp(MIN_VALID_REFL, MAX_VALID_REFL);
                    // sr_flat[iband] is &Vec<f64>; pix is unique per (iline, isamp)
                    // but we can't get &mut from par_iter without unsafe on Vec either.
                    // Use the unsafe raw pointer approach:
                    unsafe {
                        let ptr = sr_flat[iband].as_ptr() as *mut f64;
                        *ptr.add(pix) = roslamb;
                    }
                }
            }
        });
    });
```

The `unsafe` here is safe because each `(iline, isamp)` maps to a unique `pix` index and no two threads write to the same `pix`. The implementer should add a `// SAFETY:` comment explaining this.

Downstream, adjust the output scaling step to read from `sr_flat[iband][pix]` instead of `sr_f32[iband][(iline, isamp)]`.

- [ ] **Step 2: Parallelize Sentinel final correction**

Replace Sentinel Step 9 (lines ~1477-1535) with the same pattern. Key differences:
- Sentinel iterates band-first: `for iband in 0..nbands { par_iter over rows }`
- B10 copies TOA directly (no atmospheric correction)
- No TOA/cos(SZA) division
- QA band is `DNS_BAND1` not `bi.coastal`

```rust
    // ── Step 9: Final per-pixel atmospheric correction ──
    let mut sr_flat: Vec<Vec<f64>> = (0..nbands)
        .map(|_| vec![0.0f64; npix])
        .collect();

    pool.install(|| {
        (0..nlines).into_par_iter().for_each(|iline| {
            for isamp in 0..nsamps {
                let pix = iline * nsamps + isamp;
                if is_fill_pixel(qaband[pix]) { continue; }

                let raot = taero[pix] as f64;
                let eps = teps[pix] as f64;

                for iband in 0..nbands {
                    let val = if iband == DNS_BAND10 {
                        toa_bands[iband][(iline, isamp)] as f64
                    } else {
                        let rotoa = toa_bands[iband][(iline, isamp)] as f64;
                        atmcorlamb2_new(
                            &atm_coeff[iband], tgo_arr[iband], iband,
                            raot, normext_p0a3[iband],
                            rotoa, lambda, eps,
                        ).clamp(MIN_VALID_REFL, MAX_VALID_REFL)
                    };

                    unsafe {
                        let ptr = sr_flat[iband].as_ptr() as *mut f64;
                        *ptr.add(pix) = val;
                    }
                }
            }
        });
    });

    // Serial post-pass: aerosol QA bits on DNS_BAND1
    for pix in 0..npix {
        if is_fill_pixel(qaband[pix]) { continue; }
        let rsurf = sband[DNS_BAND1][pix] as f64;
        let roslamb = sr_flat[DNS_BAND1][pix];
        let tmpf = (rsurf - roslamb).abs();
        if tmpf <= LOW_AERO_THRESH {
            ipflag[pix] |= 1u8 << AERO1_QA;
        } else if tmpf < AVG_AERO_THRESH {
            ipflag[pix] |= 1u8 << AERO2_QA;
        } else {
            ipflag[pix] |= (1u8 << AERO1_QA) | (1u8 << AERO2_QA);
        }
    }
```

- [ ] **Step 3: Update output scaling to use `sr_flat` instead of `sr_f32`**

For Landsat (Step 8), change:

```rust
    let sr_bands: Vec<Array2<u16>> = sr_f32
        .iter()
        .map(|band| {
            Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
                let pix = i * nsamps + j;
                if is_fill_pixel(qa_flat[pix]) {
                    0u16
                } else {
                    let sband_f32 = band[(i, j)] as f32;
```

to:

```rust
    let sr_bands: Vec<Array2<u16>> = sr_flat
        .iter()
        .map(|band| {
            Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
                let pix = i * nsamps + j;
                if is_fill_pixel(qa_flat[pix]) {
                    0u16
                } else {
                    let sband_f32 = band[pix] as f32;
```

Apply the same change for Sentinel Step 10 (use `sr_flat[iband][pix]` instead of `sr_f32[iband][(iline, isamp)]`).

- [ ] **Step 4: Verify it compiles and tests pass**

```bash
cargo test -p lasrc_core
```

Expected: All tests pass.

- [ ] **Step 5: Build and run both comparisons**

```bash
source .venv/bin/activate
maturin develop --release -m lasrc_py/Cargo.toml
python lasrc_py/run_test_granule.py --format espa
python lasrc_py/compare_outputs.py --sensor landsat
python lasrc_py/run_test_sentinel.py
python lasrc_py/compare_outputs.py --sensor sentinel
```

Expected: Both within existing thresholds.

- [ ] **Step 6: Commit**

```bash
git add lasrc_core/src/correction.rs
git commit -m "feat: parallelize final per-pixel atmospheric correction with Rayon"
```

---

### Task 5: Parallelize Climatological SR (Both Pipelines)

**Files:**
- Modify: `lasrc_core/src/correction.rs:496-530` (Landsat Step 4)
- Modify: `lasrc_core/src/correction.rs:1065-1123` (Sentinel Step 6)

- [ ] **Step 1: Add fill pre-pass for Landsat**

`qa_flat` is already computed at line 478 (`let qa_flat: Vec<u16> = qa_band.iter().copied().collect();`). Add a serial pre-pass after the `ipflag` allocation and before the climatological SR loop to mark fill pixels:

```rust
    // Pre-pass: mark fill pixels in ipflag before parallel SR loop
    for pix in 0..npix {
        if is_fill_pixel(qa_flat[pix]) {
            ipflag[pix] = 1u8 << IPFLAG_FILL;
        }
    }
```

This replaces the fill-marking that was previously interleaved in the climatological SR loop.

- [ ] **Step 2: Parallelize Landsat climatological SR**

Replace the Landsat Step 4 loop with:

```rust
    // ── Step 4: Climatological per-pixel atmospheric correction ──
    let mut sband: Vec<Vec<f32>> = (0..nbands)
        .map(|_| vec![0.0f32; npix])
        .collect();

    pool.install(|| {
        (0..nlines).into_par_iter().for_each(|iline| {
            for isamp in 0..nsamps {
                let pix = iline * nsamps + isamp;
                if is_fill_pixel(qa_flat[pix]) { continue; }

                let xmus = (solar_zenith[(iline, isamp)] as f64 * DEG2RAD).cos();

                for iband in 0..nbands {
                    let raw_toa = toa_bands[iband][(iline, isamp)] as f64;
                    let rotoa = (raw_toa / xmus).clamp(MIN_VALID_REFL, MAX_VALID_REFL) as f32;

                    let tgo_x_roatm = btgo[iband] as f32 * broatm[iband] as f32;
                    let tgo_x_ttatmg = btgo[iband] as f32 * bttatmg[iband] as f32;
                    let roslamb: f32 = {
                        let num = rotoa - tgo_x_roatm;
                        num / (tgo_x_ttatmg + bsatm[iband] as f32 * num)
                    };

                    unsafe {
                        let ptr = sband[iband].as_ptr() as *mut f32;
                        *ptr.add(pix) = roslamb.clamp(MIN_VALID_REFL as f32, MAX_VALID_REFL as f32);
                    }
                }
            }
        });
    });
```

- [ ] **Step 3: Parallelize Sentinel climatological SR**

Replace the Sentinel Step 6 inner pixel loop with:

```rust
    // ── Step 6: Climatological per-pixel atmospheric correction ──
    let mut sband: Vec<Vec<f32>> = (0..nbands)
        .map(|_| vec![0.0f32; npix])
        .collect();

    let tauray = sensor.tauray();
    let max_band_idx = lambda.len() - 1;

    for iband in 0..nbands {
        let (tgo_x_roatm, tgo_x_ttatmg, satm_val) = if iband == DNS_BAND9 {
            (0.0f32, 1.0f32, 0.0f32)
        } else {
            let result = atmcorlamb2(
                lut, &gas_coeff[iband], tauray[iband], iband,
                xts, xtv, xmus, xmuv, xfi, cosxfi,
                0.05, pressure, uoz, uwv, 0.0, lambda, max_band_idx,
                -1.0,
            );
            (
                result.tgo as f32 * result.roatm as f32,
                result.tgo as f32 * result.ttatmg as f32,
                result.satm as f32,
            )
        };

        pool.install(|| {
            (0..npix).into_par_iter().for_each(|pix| {
                let (i, j) = (pix / nsamps, pix % nsamps);
                if is_fill_pixel(qaband[pix]) { return; }

                let rotoa = toa_bands[iband][(i, j)];
                let roslamb: f32 = {
                    let num = rotoa - tgo_x_roatm;
                    num / (tgo_x_ttatmg + satm_val * num)
                };

                unsafe {
                    let ptr = sband[iband].as_ptr() as *mut f32;
                    *ptr.add(pix) = roslamb;
                }
            });
        });
    }
```

- [ ] **Step 4: Verify it compiles and tests pass**

```bash
cargo test -p lasrc_core
```

Expected: All tests pass.

- [ ] **Step 5: Build and run both comparisons**

```bash
source .venv/bin/activate
maturin develop --release -m lasrc_py/Cargo.toml
python lasrc_py/run_test_granule.py --format espa
python lasrc_py/compare_outputs.py --sensor landsat
python lasrc_py/run_test_sentinel.py
python lasrc_py/compare_outputs.py --sensor sentinel
```

Expected: Both within existing thresholds.

- [ ] **Step 6: Commit**

```bash
git add lasrc_core/src/correction.rs
git commit -m "feat: parallelize climatological SR computation with Rayon"
```

---

### Task 6: Verify and Benchmark

**Files:**
- No code changes

- [ ] **Step 1: Run unit tests**

```bash
cargo test -p lasrc_core
```

Expected: All tests pass.

- [ ] **Step 2: Run Landsat comparison with num_threads=1**

Temporarily add `num_threads=1` to the Python runner call to verify serial-equivalent output. Edit `lasrc_py/run_test_granule.py` to pass `num_threads=1` in the call to `process_surface_reflectance`, rebuild, run, and compare. Output should be identical to the pre-parallelization baseline.

```bash
source .venv/bin/activate
maturin develop --release -m lasrc_py/Cargo.toml
python lasrc_py/run_test_granule.py --format espa
python lasrc_py/compare_outputs.py --sensor landsat
```

Expected: Same results as before parallelization.

- [ ] **Step 3: Run Landsat comparison with default threads**

Remove the `num_threads=1` override (or just don't pass it), rebuild, run, and compare.

Expected: Results within existing thresholds. May have minor differences from thread scheduling order affecting floating-point accumulation.

- [ ] **Step 4: Run Sentinel comparison with default threads**

```bash
python lasrc_py/run_test_sentinel.py
python lasrc_py/compare_outputs.py --sensor sentinel
```

Expected: Results within existing thresholds.

- [ ] **Step 5: Time comparison**

Run both pipelines and note the wall-clock time from the "Completed in X.Xs" output. Compare to the pre-parallelization times (Landsat ~14.6s, Sentinel ~26.1s from earlier runs).

- [ ] **Step 6: Revert any temporary test changes and commit**

If you modified `run_test_granule.py` for the `num_threads=1` test, revert those changes.

```bash
git checkout -- lasrc_py/run_test_granule.py
```
