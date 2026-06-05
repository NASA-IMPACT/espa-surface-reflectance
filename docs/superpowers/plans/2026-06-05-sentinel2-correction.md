# Sentinel-2 Surface Reflectance Correction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Sentinel-2 atmospheric correction to the Rust+Python LaSRC port, processing all 13 bands from SAFE archive input, matching the C reference output.

**Architecture:** The Sentinel-2 path reuses all existing Rust core modules (LUT, atmospheric, aerosol retrieval, gas transmission, Rayleigh, geometry). Changes are: (1) new 13-band constants and sensor config, (2) a separate Sentinel correction orchestrator with Sentinel-specific fill detection, aerosol window averaging, NDWI formula, and three new aerosol post-processing functions, (3) Python SAFE archive I/O and pipeline, (4) updated comparison tooling.

**Tech Stack:** Rust (lasrc_core), PyO3 (lasrc_py), Python (pipeline, I/O), pyhdf (HDF4 LUTs), h5py (VIIRS aux), rasterio/GDAL (JP2 reading), numpy

**Spec:** `docs/superpowers/specs/2026-06-05-sentinel2-correction-design.md`

**C reference:** `lasrc/c_version/src/compute_sentinel_refl.c` and `lasrc/c_version/src/aero_interp.c`

---

## File Structure

### Files to modify (Rust)

| File | Responsibility |
|------|---------------|
| `lasrc_core/src/constants.rs` | Add 13-band Sentinel tauray, lambda, gas coefficients, IPFLAG_FAILED_TMP, Sentinel TTH arrays |
| `lasrc_core/src/sensor.rs` | Expand SENTINEL_REFL_BANDS to 13, update SENTINEL_BAND_INDICES, fill todo!() stubs |
| `lasrc_core/src/aerosol.rs` | Add `aerosol_interp_sentinel`, `ipflag_expand_failed_sentinel`, `aero_avg_failed_sentinel` |
| `lasrc_core/src/correction.rs` | Add `compute_sentinel_surface_reflectance` orchestrator |

### Files to create (Python)

| File | Responsibility |
|------|---------------|
| `lasrc_py/python/lasrc/io_sentinel.py` | Read SAFE archive: JP2 bands, XML metadata, resample to 10m |
| `lasrc_py/python/lasrc/pipeline_sentinel.py` | Sentinel-2 pipeline: read SAFE, load LUTs, call Rust, write output |
| `lasrc_py/run_test_sentinel.py` | Test runner for the S2B test granule |

### Files to modify (Python)

| File | Responsibility |
|------|---------------|
| `lasrc_py/python/lasrc/sensors.py` | Add SENTINEL_2B, SENTINEL_2C configs with 13-band names |
| `lasrc_py/compare_outputs.py` | Add `--sensor sentinel` mode with 13-band S2 names and 10980x10980 dimensions |

---

### Task 1: Add Sentinel-2 constants (13 bands)

**Files:**
- Modify: `lasrc_core/src/constants.rs`
- Test: `cargo test -p lasrc_core`

- [ ] **Step 1: Add TAURAY, LAMBDA, and gas coefficient arrays for 13 Sentinel bands**

Add after the existing Landsat constants (after line ~152):

```rust
// Sentinel-2 Rayleigh optical depth per band (13 bands, from tauray-msi.ASC)
// Order: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12
pub const TAURAY_SENTINEL: [f64; 13] = [
    0.23432, 0.15106, 0.09102, 0.04535, 0.03584, 0.02924, 0.02338, 0.01847,
    0.01560, 0.01092, 0.00243, 0.00128, 0.00037,
];

// Sentinel-2 band wavelengths (micrometers), all 13 bands
pub const LAMBDA_SENTINEL_ALL: [f64; 13] = [
    0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842,
    0.865, 0.945, 1.375, 1.610, 2.190,
];

// Sentinel-2 band count (all 13 bands including B09/B10)
pub const NSRS_BANDS: usize = 13;

// Gas transmission coefficients per Sentinel band (from gascoef-msi.ASC)
// Order: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12
pub const OZTRANSA_SENTINEL: [f64; 13] = [
    -0.00264691, -0.0272572, -0.0986512, -0.0500348, -0.0204295,
    -0.0108641, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
];
pub const WVTRANSA_SENTINEL: [f64; 13] = [
    2.29849e-27, 2.29849e-27, 0.000777307, 0.00361051, 0.0141249,
    0.0137067, 0.00410217, 0.0285871, 0.000390755, 0.00001, 0.01,
    0.000640155, 0.018006,
];
pub const WVTRANSB_SENTINEL: [f64; 13] = [
    0.999742, 0.999742, 0.891099, 0.754895, 0.75596, 0.763497, 0.74117,
    0.578722, 0.900899, 0.45818, 1.0, 0.943712, 0.647517,
];
pub const OGTRANSA1_SENTINEL: [f64; 13] = [
    4.91586e-20, 4.91586e-20, 4.91586e-20, 4.91586e-20, 5.3367e-06,
    4.91586e-20, 9.03583e-05, 1.64109e-09, 1.90458e-05, 4.91586e-20,
    7.62429e-06, 0.0212751, 0.0243065,
];
pub const OGTRANSB0_SENTINEL: [f64; 13] = [
    0.000197019, 0.000197019, 0.000197019, 0.000197019, -0.980313,
    0.000197019, 0.0265393, 1.0e-10, 0.0322844, 0.000197019, 0.000197019,
    0.000197019, 0.000197019,
];
pub const OGTRANSB1_SENTINEL: [f64; 13] = [
    9.57011e-16, 9.57011e-16, 9.57011e-16, 9.57011e-16, 1.33639,
    9.57011e-16, 0.0532256, 1.0e-10, -0.0219907, 9.57011e-16, -0.216849,
    0.0116062, 0.0604312,
];
```

- [ ] **Step 2: Add IPFLAG_FAILED_TMP and Sentinel band index constants**

Add alongside the existing IPFLAG constants:

```rust
// Sentinel-specific QA flag (shares bit 5 with IPFLAG_INTERP_WINDOW, different sensor)
pub const IPFLAG_FAILED_TMP: u8 = 5;

// Sentinel-2 band indices (for the 13-band array)
pub const DNS_BAND1: usize = 0;   // B01 - Coastal aerosol
pub const DNS_BAND2: usize = 1;   // B02 - Blue
pub const DNS_BAND3: usize = 2;   // B03 - Green
pub const DNS_BAND4: usize = 3;   // B04 - Red (aerosol reference)
pub const DNS_BAND5: usize = 4;   // B05
pub const DNS_BAND6: usize = 5;   // B06
pub const DNS_BAND7: usize = 6;   // B07
pub const DNS_BAND8: usize = 7;   // B08 - NIR
pub const DNS_BAND8A: usize = 8;  // B8A - NIR narrow (NDWI/NDVI)
pub const DNS_BAND9: usize = 9;   // B09 - Water vapor (tgo=1 bypass)
pub const DNS_BAND10: usize = 10; // B10 - Cirrus (copy TOA)
pub const DNS_BAND11: usize = 11; // B11 - SWIR1
pub const DNS_BAND12: usize = 12; // B12 - SWIR2 (NDWI, aerosol)

// Sentinel-2 surface reflectance threshold arrays (13 bands)
// All zeros — Sentinel does not use the tth negative-check in subaeroret_new
pub const SENTINEL_TTH: [f64; 13] = [0.0; 13];
```

- [ ] **Step 3: Run tests**

Run: `cargo test -p lasrc_core`
Expected: All existing tests pass. New constants compile without error.

- [ ] **Step 4: Commit**

```bash
git add lasrc_core/src/constants.rs
git commit -m "feat: add Sentinel-2 13-band constants (tauray, lambda, gas coefficients)"
```

---

### Task 2: Update sensor.rs for 13-band Sentinel config

**Files:**
- Modify: `lasrc_core/src/sensor.rs`
- Test: `cargo test -p lasrc_core`

- [ ] **Step 1: Add imports for new Sentinel constants**

Update the imports at the top of `sensor.rs`:

```rust
use crate::constants::{
    LAERO_WINDOW, LHALF_AERO_WINDOW, LFIX_AERO_WINDOW, LHALF_FIX_AERO_WINDOW, LMIN_CLEAR_PIX,
    SAERO_WINDOW, SFIX_AERO_WINDOW, SHALF_FIX_AERO_WINDOW, SMIN_CLEAR_PIX,
    TAURAY_LANDSAT, TAURAY_SENTINEL, LAMBDA_LANDSAT, LAMBDA_SENTINEL_ALL,
    OZTRANSA_LANDSAT, WVTRANSA_LANDSAT, WVTRANSB_LANDSAT,
    OGTRANSA1_LANDSAT, OGTRANSB0_LANDSAT, OGTRANSB1_LANDSAT,
    OZTRANSA_SENTINEL, WVTRANSA_SENTINEL, WVTRANSB_SENTINEL,
    OGTRANSA1_SENTINEL, OGTRANSB0_SENTINEL, OGTRANSB1_SENTINEL,
};
```

- [ ] **Step 2: Expand SENTINEL_REFL_BANDS to 13 bands**

Replace the existing 11-band `SENTINEL_REFL_BANDS` with:

```rust
const SENTINEL_REFL_BANDS: [BandConfig; 13] = [
    BandConfig { name: "B01", wavelength_um: 0.443, native_resolution_m: 60.0 },
    BandConfig { name: "B02", wavelength_um: 0.490, native_resolution_m: 10.0 },
    BandConfig { name: "B03", wavelength_um: 0.560, native_resolution_m: 10.0 },
    BandConfig { name: "B04", wavelength_um: 0.665, native_resolution_m: 10.0 },
    BandConfig { name: "B05", wavelength_um: 0.705, native_resolution_m: 20.0 },
    BandConfig { name: "B06", wavelength_um: 0.740, native_resolution_m: 20.0 },
    BandConfig { name: "B07", wavelength_um: 0.783, native_resolution_m: 20.0 },
    BandConfig { name: "B08", wavelength_um: 0.842, native_resolution_m: 10.0 },
    BandConfig { name: "B8A", wavelength_um: 0.865, native_resolution_m: 20.0 },
    BandConfig { name: "B09", wavelength_um: 0.945, native_resolution_m: 60.0 },
    BandConfig { name: "B10", wavelength_um: 1.375, native_resolution_m: 60.0 },
    BandConfig { name: "B11", wavelength_um: 1.610, native_resolution_m: 20.0 },
    BandConfig { name: "B12", wavelength_um: 2.190, native_resolution_m: 20.0 },
];
```

- [ ] **Step 3: Update SENTINEL_BAND_INDICES for 13-band layout**

```rust
const SENTINEL_BAND_INDICES: BandIndices = BandIndices {
    coastal: 0,   // B01
    blue: 1,      // B02
    green: 2,     // B03
    red: 3,       // B04
    nir: 8,       // B8A (used in NDVI/NDWI)
    swir1: 11,    // B11
    swir2: 12,    // B12
};
```

- [ ] **Step 4: Add sentinel_gas_coefficients helper and fill the todo!() stubs**

Add the helper function:

```rust
fn sentinel_gas_coefficients() -> Vec<GasCoefficients> {
    (0..OZTRANSA_SENTINEL.len())
        .map(|i| GasCoefficients {
            oztransa: OZTRANSA_SENTINEL[i],
            wvtransa: WVTRANSA_SENTINEL[i],
            wvtransb: WVTRANSB_SENTINEL[i],
            ogtransa1: OGTRANSA1_SENTINEL[i],
            ogtransb0: OGTRANSB0_SENTINEL[i],
            ogtransb1: OGTRANSB1_SENTINEL[i],
        })
        .collect()
}
```

Update the `impl_sentinel_sensor!` macro to replace the two `todo!()` stubs:

```rust
fn tauray(&self) -> &[f64] { &TAURAY_SENTINEL }
fn lambda(&self) -> &[f64] { &LAMBDA_SENTINEL_ALL }
fn gas_coefficients(&self) -> Vec<GasCoefficients> {
    sentinel_gas_coefficients()
}
```

- [ ] **Step 5: Update tests for 13 bands**

Update the existing Sentinel test assertions:

```rust
#[test]
fn test_sentinel2a_band_count() {
    let sensor = Sentinel2A;
    assert_eq!(sensor.reflectance_bands().len(), 13);
    assert_eq!(sensor.thermal_bands().len(), 0);
}

#[test]
fn test_sentinel2a_band_indices() {
    let sensor = Sentinel2A;
    let idx = sensor.band_indices();
    assert_eq!(idx.coastal, 0);
    assert_eq!(idx.blue, 1);
    assert_eq!(idx.red, 3);
    assert_eq!(idx.nir, 8);    // B8A
    assert_eq!(idx.swir1, 11); // B11
    assert_eq!(idx.swir2, 12); // B12
}

#[test]
fn test_sentinel2a_lambda() {
    let sensor = Sentinel2A;
    assert_eq!(sensor.lambda().len(), 13);
    assert!((sensor.lambda()[0] - 0.443).abs() < 1e-5);
    assert!((sensor.lambda()[12] - 2.190).abs() < 1e-5);
}

#[test]
fn test_sentinel2a_tauray() {
    let sensor = Sentinel2A;
    assert_eq!(sensor.tauray().len(), 13);
    assert!((sensor.tauray()[0] - 0.23432).abs() < 1e-5);
}

#[test]
fn test_sentinel2a_gas_coefficients() {
    let sensor = Sentinel2A;
    let gc = sensor.gas_coefficients();
    assert_eq!(gc.len(), 13);
    assert!((gc[0].oztransa - (-0.00264691)).abs() < 1e-8);
}

#[test]
fn test_num_refl_bands_default_impl() {
    assert_eq!(Landsat8.num_refl_bands(), 8);
    assert_eq!(Sentinel2A.num_refl_bands(), 13);
}
```

- [ ] **Step 6: Run tests**

Run: `cargo test -p lasrc_core`
Expected: All tests pass including updated Sentinel tests.

- [ ] **Step 7: Commit**

```bash
git add lasrc_core/src/sensor.rs
git commit -m "feat: expand Sentinel-2 sensor config to 13 bands, fill tauray/gas_coefficients stubs"
```

---

### Task 3: Add three Sentinel aerosol post-processing functions

**Files:**
- Modify: `lasrc_core/src/aerosol.rs`
- Test: `cargo test -p lasrc_core`

**C reference:**
- `aerosol_interp_sentinel`: `lasrc/c_version/src/aero_interp.c:527-637`
- `ipflag_expand_failed_sentinel`: `lasrc/c_version/src/aero_interp.c:658-722`
- `aero_avg_failed_sentinel`: `lasrc/c_version/src/aero_interp.c:741-927`

- [ ] **Step 1: Add `aerosol_interp_sentinel`**

This function interpolates aerosol values within 6x6 windows from the UL corner pixel to all other pixels. It differs from the Landsat `aerosol_interp` — it uses the 4 surrounding UL corners of adjacent windows and handles edge cases differently.

Add to `aerosol.rs`:

```rust
/// Bilinear interpolation of aerosol from window UL corners for Sentinel-2.
///
/// Ported from `aerosol_interp_sentinel()` in C `aero_interp.c:527-637`.
///
/// Unlike the Landsat version, this interpolates from UL corner pixels
/// (not window centers) using a bilinear weighting across the 4 surrounding
/// UL corner pixels. Handles edge cases where awline >= nlines or
/// awsamp >= nsamps by falling back to the current UL corner value.
pub fn aerosol_interp_sentinel(
    aero_window: usize,
    qaband: &[u16],
    ipflag: &mut [u8],
    taero: &mut [f32],
    nlines: usize,
    nsamps: usize,
) {
    let sq_aero_win = (aero_window * aero_window) as f32;

    for line in 0..nlines {
        let awline = line + aero_window;

        for samp in 0..nsamps {
            let curr_pix = line * nsamps + samp;

            // Skip fill
            if crate::constants::is_fill_pixel(qaband[curr_pix]) {
                continue;
            }

            let awsamp = samp + aero_window;

            // Pixel indices for the 4 UL corners
            let next_samp_pix = line * nsamps + awsamp;
            let next_line_pix = awline * nsamps + samp;
            let next_line_samp_pix = awline * nsamps + awsamp;

            // Loop through NxN window with current pixel as UL
            for iline in line..awline {
                if iline >= nlines {
                    continue;
                }
                let awline_iline = (awline - iline) as f32;
                let iline_line = (iline - line) as f32;

                for isamp in samp..awsamp {
                    if isamp >= nsamps {
                        continue;
                    }

                    let curr_win_pix = iline * nsamps + isamp;
                    if crate::constants::is_fill_pixel(qaband[curr_win_pix]) {
                        continue;
                    }

                    let awsamp_isamp = (awsamp - isamp) as f32;
                    let isamp_samp = (isamp - samp) as f32;

                    // Start with contribution from current UL corner
                    let mut val = taero[curr_pix] * awline_iline * awsamp_isamp;

                    // Add contributions from surrounding UL corners based on edge cases
                    if awline < nlines && awsamp < nsamps {
                        val += isamp_samp * awline_iline * taero[next_samp_pix]
                            + awsamp_isamp * iline_line * taero[next_line_pix]
                            + isamp_samp * iline_line * taero[next_line_samp_pix];
                    } else if awline >= nlines && awsamp < nsamps {
                        val += isamp_samp * awline_iline * taero[next_samp_pix]
                            + awsamp_isamp * iline_line * taero[curr_pix]
                            + isamp_samp * iline_line * taero[next_samp_pix];
                    } else if awline < nlines && awsamp >= nsamps {
                        val += isamp_samp * awline_iline * taero[curr_pix]
                            + awsamp_isamp * iline_line * taero[next_line_pix]
                            + isamp_samp * iline_line * taero[next_line_pix];
                    } else {
                        // Both awline >= nlines and awsamp >= nsamps
                        val += isamp_samp * awline_iline * taero[curr_pix]
                            + awsamp_isamp * iline_line * taero[curr_pix]
                            + isamp_samp * iline_line * taero[curr_pix];
                    }

                    taero[curr_win_pix] = val / sq_aero_win;
                }
            }
        }
    }
}
```

- [ ] **Step 2: Add `ipflag_expand_failed_sentinel`**

Ported from `aero_interp.c:658-722`. Expands `IPFLAG_FAILED` pixels to a ±12 pixel radius using `IPFLAG_FAILED_TMP` as an intermediate, skipping fill and water pixels.

```rust
/// Expand failed aerosol pixels to surrounding area for Sentinel-2.
///
/// Ported from `ipflag_expand_failed_sentinel()` in C `aero_interp.c:658-722`.
///
/// For each pixel with IPFLAG_FAILED set, mark all non-fill, non-water
/// pixels within HALF_EXPAND_WIN radius as IPFLAG_FAILED_TMP. Then
/// convert all FAILED_TMP to FAILED and clear the TMP bit.
pub fn ipflag_expand_failed_sentinel(
    ipflag: &mut [u8],
    nlines: usize,
    nsamps: usize,
) {
    let half_win = HALF_EXPAND_WIN as isize;

    // Pass 1: Mark surrounding pixels with FAILED_TMP
    for line in 0..nlines {
        for samp in 0..nsamps {
            let curr_pix = line * nsamps + samp;

            // Only expand from FAILED pixels
            if ipflag[curr_pix] & (1u8 << IPFLAG_FAILED) == 0 {
                continue;
            }

            for iline in -half_win..=half_win {
                let win_line = line as isize + iline;
                if win_line < 0 || win_line >= nlines as isize {
                    continue;
                }
                for isamp in -half_win..=half_win {
                    let win_samp = samp as isize + isamp;
                    if win_samp < 0 || win_samp >= nsamps as isize {
                        continue;
                    }

                    let curr_win_pix = win_line as usize * nsamps + win_samp as usize;
                    // Skip fill, water, and already-temp-failed pixels
                    if ipflag[curr_win_pix] & (1u8 << IPFLAG_FILL) != 0 {
                        continue;
                    }
                    if ipflag[curr_win_pix] & (1u8 << IPFLAG_WATER) != 0 {
                        continue;
                    }
                    if ipflag[curr_win_pix] & (1u8 << IPFLAG_FAILED_TMP) != 0 {
                        continue;
                    }
                    ipflag[curr_win_pix] |= 1u8 << IPFLAG_FAILED_TMP;
                }
            }
        }
    }

    // Pass 2: Convert FAILED_TMP to FAILED and clear TMP bit
    let npixels = nlines * nsamps;
    for pix in 0..npixels {
        if ipflag[pix] & (1u8 << IPFLAG_FAILED_TMP) != 0 {
            ipflag[pix] |= 1u8 << IPFLAG_FAILED;
            ipflag[pix] &= !(1u8 << IPFLAG_FAILED_TMP);
        }
    }
}
```

- [ ] **Step 3: Add `aero_avg_failed_sentinel`**

Ported from `aero_interp.c:741-927`. Two-pass averaging for failed pixels using separate taeros/tepss/smflag arrays.

```rust
/// Average aerosol values for failed Sentinel-2 pixels.
///
/// Ported from `aero_avg_failed_sentinel()` in C `aero_interp.c:741-927`.
///
/// Pass 1: For each pixel, compute average of non-fill, non-failed pixels
/// within ±HALF_FAILED_WIN. Require MIN_VALID_WINDOW_PIX valid neighbors.
/// Results stored in temporary arrays (taeros, tepss).
///
/// Pass 2+: For unfilled pixels, expand using already-filled neighbors.
/// Repeat until all pixels filled or stuck, then use defaults.
///
/// Finally, copy averaged values back to taero/teps for failed pixels.
pub fn aero_avg_failed_sentinel(
    qaband: &[u16],
    ipflag: &mut [u8],
    taero: &mut [f32],
    teps: &mut [f32],
    nlines: usize,
    nsamps: usize,
) {
    let npixels = nlines * nsamps;
    let half_win = HALF_FAILED_WIN as isize;

    let mut taeros = vec![0.0f32; npixels];
    let mut tepss = vec![0.0f32; npixels];
    let mut smflag = vec![false; npixels];

    // Pass 1: Average from non-fill, non-failed neighbors
    let mut one_filled = false;
    let mut nbpixnf = 0usize;

    for line in 0..nlines {
        for samp in 0..nsamps {
            let curr_pix = line * nsamps + samp;
            smflag[curr_pix] = false;

            if crate::constants::is_fill_pixel(qaband[curr_pix]) {
                continue;
            }

            let mut taerosum: f32 = 0.0;
            let mut tepssum: f32 = 0.0;
            let mut nbaeroavg = 0usize;

            for iline in -half_win..=half_win {
                let wl = line as isize + iline;
                if wl < 0 || wl >= nlines as isize {
                    continue;
                }
                for isamp in -half_win..=half_win {
                    let ws = samp as isize + isamp;
                    if ws < 0 || ws >= nsamps as isize {
                        continue;
                    }

                    let curr_win_pix = wl as usize * nsamps + ws as usize;
                    // Include non-fill, non-failed pixels
                    if ipflag[curr_win_pix] & (1u8 << IPFLAG_FILL) == 0
                        && ipflag[curr_win_pix] & (1u8 << IPFLAG_FAILED) == 0
                    {
                        nbaeroavg += 1;
                        taerosum += taero[curr_win_pix];
                        tepssum += teps[curr_win_pix];
                    }
                }
            }

            if nbaeroavg > MIN_VALID_WINDOW_PIX {
                taeros[curr_pix] = taerosum / nbaeroavg as f32;
                tepss[curr_pix] = tepssum / nbaeroavg as f32;
                smflag[curr_pix] = true;
                one_filled = true;
            } else {
                nbpixnf += 1;
            }
        }
    }

    // If nothing filled, use defaults for everything
    if !one_filled {
        for pix in 0..npixels {
            if ipflag[pix] & (1u8 << IPFLAG_FILL) == 0 {
                taero[pix] = DEFAULT_AERO as f32;
                teps[pix] = DEFAULT_EPS as f32;
                smflag[pix] = true;
            }
        }
        // Copy back and return
        for pix in 0..npixels {
            if smflag[pix] {
                taero[pix] = taeros[pix].max(taero[pix]); // defaults already set
                teps[pix] = tepss[pix].max(teps[pix]);
            }
        }
        return;
    }

    // Pass 2+: Fill remaining pixels using already-filled neighbors
    while nbpixnf > 0 {
        let prev_nbpixnf = nbpixnf;
        nbpixnf = 0;

        for line in 0..nlines {
            for samp in 0..nsamps {
                let curr_pix = line * nsamps + samp;

                if crate::constants::is_fill_pixel(qaband[curr_pix]) || smflag[curr_pix] {
                    continue;
                }

                let mut taerosum: f32 = 0.0;
                let mut tepssum: f32 = 0.0;
                let mut nbaeroavg = 0usize;

                for iline in -half_win..=half_win {
                    let wl = line as isize + iline;
                    if wl < 0 || wl >= nlines as isize {
                        continue;
                    }
                    for isamp in -half_win..=half_win {
                        let ws = samp as isize + isamp;
                        if ws < 0 || ws >= nsamps as isize {
                            continue;
                        }

                        let curr_win_pix = wl as usize * nsamps + ws as usize;
                        if smflag[curr_win_pix] {
                            nbaeroavg += 1;
                            taerosum += taeros[curr_win_pix];
                            tepssum += tepss[curr_win_pix];
                        }
                    }
                }

                if nbaeroavg > 0 {
                    taeros[curr_pix] = taerosum / nbaeroavg as f32;
                    tepss[curr_pix] = tepssum / nbaeroavg as f32;
                    smflag[curr_pix] = true;
                } else {
                    nbpixnf += 1;
                }
            }
        }

        // If no progress, fill remaining with defaults
        if nbpixnf >= prev_nbpixnf {
            for pix in 0..npixels {
                if !smflag[pix] && ipflag[pix] & (1u8 << IPFLAG_FILL) == 0 {
                    taeros[pix] = DEFAULT_AERO as f32;
                    tepss[pix] = DEFAULT_EPS as f32;
                    smflag[pix] = true;
                }
            }
            break;
        }
    }

    // Copy averaged values back to taero/teps for filled pixels
    for pix in 0..npixels {
        if smflag[pix] {
            taero[pix] = taeros[pix];
            teps[pix] = tepss[pix];
        }
    }
}
```

- [ ] **Step 4: Write tests for the three functions**

```rust
#[test]
fn test_aerosol_interp_sentinel_basic() {
    // 12x12 image, aero_window=6
    let nlines = 12;
    let nsamps = 12;
    let npix = nlines * nsamps;
    let mut taero = vec![0.0f32; npix];
    let mut ipflag = vec![0u8; npix];
    let qaband = vec![0u16; npix]; // no fill

    // Set UL corners (0,0), (0,6), (6,0), (6,6) with known values
    taero[0] = 0.1;
    taero[6] = 0.2;
    taero[6 * nsamps] = 0.3;
    taero[6 * nsamps + 6] = 0.4;

    aerosol_interp_sentinel(6, &qaband, &mut ipflag, &mut taero, nlines, nsamps);

    // Center of first window (3,3) should be interpolated
    let pix33 = 3 * nsamps + 3;
    assert!(taero[pix33] > 0.0, "Interior pixel should be interpolated");
    // UL corner should keep its value (after self-interpolation)
    assert!((taero[0] - 0.1).abs() < 0.01);
}

#[test]
fn test_ipflag_expand_failed_sentinel() {
    let nlines = 30;
    let nsamps = 30;
    let npix = nlines * nsamps;
    let mut ipflag = vec![0u8; npix];

    // Set center pixel as failed
    let center = 15 * nsamps + 15;
    ipflag[center] = 1u8 << IPFLAG_FAILED;

    ipflag_expand_failed_sentinel(&mut ipflag, nlines, nsamps);

    // Pixel 12 away should be marked as failed
    let edge_pix = 3 * nsamps + 15; // 12 lines away
    assert!(ipflag[edge_pix] & (1u8 << IPFLAG_FAILED) != 0);

    // Pixel 13 away should NOT be marked
    let far_pix = 2 * nsamps + 15; // 13 lines away
    assert!(ipflag[far_pix] & (1u8 << IPFLAG_FAILED) == 0);
}

#[test]
fn test_aero_avg_failed_sentinel() {
    let nlines = 10;
    let nsamps = 10;
    let npix = nlines * nsamps;
    let qaband = vec![0u16; npix]; // no fill
    let mut ipflag = vec![0u8; npix];
    let mut taero = vec![0.1f32; npix];
    let mut teps = vec![1.5f32; npix];

    // Mark center pixel as failed
    let center = 5 * nsamps + 5;
    ipflag[center] = 1u8 << IPFLAG_FAILED;
    taero[center] = 0.0;
    teps[center] = 0.0;

    aero_avg_failed_sentinel(&qaband, &mut ipflag, &mut taero, &mut teps, nlines, nsamps);

    // Failed pixel should now have averaged values from neighbors
    assert!(taero[center] > 0.05, "Failed pixel should be filled: {}", taero[center]);
    assert!(teps[center] > 0.5, "Failed pixel eps should be filled: {}", teps[center]);
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p lasrc_core`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add lasrc_core/src/aerosol.rs
git commit -m "feat: add Sentinel-2 aerosol post-processing (interp, expand_failed, avg_failed)"
```

---

### Task 4: Add Sentinel correction orchestrator

**Files:**
- Modify: `lasrc_core/src/correction.rs`
- Test: `cargo test -p lasrc_core`

**C reference:** `lasrc/c_version/src/compute_sentinel_refl.c:277-1915`

This is the largest task. The Sentinel orchestrator differs from Landsat in:
1. Fill detection: ANY band == 0.0 (not QA_PIXEL bit)
2. Aerosol window: starts at (0,0) stepping by 6 (not center-offset)
3. TOA averaging within 6x6 windows for aerosol retrieval
4. NDWI formula uses B8A and B12*0.5
5. Eps optimization uses quadratic on residuals (not subaeroret_new refit)
6. NDVI test uses NxN-averaged TOA for B8A and B04
7. Water test uses B01/B04/B8A/B12 with erelc=1.0
8. Water failure sets IPFLAG_FAILED (not just clearing bits)
9. Post-processing: aerosol_interp_sentinel + ipflag_expand_failed + aero_avg_failed
10. Final correction: B09 skip (already handled in precompute), B10 copy TOA
11. Aerosol QA on B01 (not B00)

- [ ] **Step 1: Add the `compute_sentinel_surface_reflectance` function signature**

Add to `correction.rs`, after the existing `compute_surface_reflectance` function:

```rust
/// Compute surface reflectance for a Sentinel-2 scene.
///
/// This is the Sentinel-specific orchestrator. Key differences from Landsat:
/// - Fill detection: any band with value 0.0 marks the pixel as fill
/// - 13 bands (including B09 water vapor bypass and B10 TOA copy)
/// - Aerosol retrieval on 6x6 windows starting from UL corner (0,0)
/// - TOA averaging within windows for aerosol bands
/// - Sentinel-specific aerosol post-processing (interp, expand_failed, avg_failed)
/// - Eps optimization via quadratic fit on 3 residuals
/// - No QA_PIXEL band — fill mask derived from band values
#[allow(clippy::too_many_arguments)]
pub fn compute_sentinel_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<f32>],
    solar_zenith: f64,    // scene-center mean, degrees
    solar_azimuth: f64,   // scene-center mean, degrees
    view_zenith: f64,     // scene-center mean, degrees
    view_azimuth: f64,    // scene-center mean, degrees
    lut: &LookupTables,
    aux: &AuxiliaryData,
    space_def: &SpaceDef,
) -> SurfaceReflectanceResult {
```

The key difference from the Landsat signature: Sentinel uses scene-center angles (scalars) instead of per-pixel angle grids, and there's no qa_band or bt_bands.

- [ ] **Step 2: Implement fill detection**

```rust
    let (nlines, nsamps) = toa_bands[0].dim();
    let nbands = sensor.num_refl_bands(); // 13
    let lambda = sensor.lambda();
    let aero_window = sensor.aerosol_window(); // 6

    let npix = nlines * nsamps;

    // Fill detection: any band with value == 0.0 marks pixel as fill
    // Build a QA-like array for fill detection (matches C: qaband with fill bit)
    let mut qaband = vec![0u16; npix];
    for pix in 0..npix {
        let (i, j) = (pix / nsamps, pix % nsamps);
        for ib in 0..nbands {
            if toa_bands[ib][(i, j)] == 0.0 {
                qaband[pix] = 1; // bit 0 = fill
                break;
            }
        }
    }
```

- [ ] **Step 3: Implement scene-center geometry and precompute coefficients**

```rust
    // Scene-center geometry
    let xts = solar_zenith;
    let xtv = view_zenith;
    let xmus = (xts * DEG2RAD).cos();
    let xmuv = (xtv * DEG2RAD).cos();
    let xfi = (view_azimuth - solar_azimuth).abs();
    let xfi = if xfi > 180.0 { 360.0 - xfi } else { xfi };
    let cosxfi = (xfi * DEG2RAD).cos();

    // Scene-center atmospheric params
    let center_line = (nlines / 2) as i32;
    let center_samp = (nsamps / 2) as i32;
    let (scene_center_lat, scene_center_lon) =
        utm_to_deg(space_def, center_line, center_samp);
    let (pressure, uoz, uwv) = extract_atm_params(aux, scene_center_lat, scene_center_lon);

    let gas_coeff: Vec<GasCoefficients> = sensor.gas_coefficients();
    let bi = sensor.band_indices();

    // Precompute polynomial coefficients (same as Landsat)
    let (atm_coeff, tgo_arr, normext_p0a3, btgo, broatm, bttatmg, bsatm) =
        precompute_coefficients(
            sensor, lut, &gas_coeff,
            xts, xtv, xmus, xmuv, xfi, cosxfi, pressure, uoz, uwv,
        );

    // Extract polynomial arrays for subaeroret_new
    let roatm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.roatm_coef).collect();
    let ttatmg_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.ttatmg_coef).collect();
    let satm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.satm_coef).collect();
    let roatm_ia_max: Vec<f64> = atm_coeff.iter().map(|c| c.roatm_upper).collect();
```

- [ ] **Step 4: Implement initial atmospheric correction (climatological)**

This matches C lines 682-752. B09 gets tgo=1, roatm=0, ttatmg=1, satm=0.

```rust
    // Allocate working arrays
    let mut ipflag = vec![0u8; npix];
    let mut taero = vec![DEFAULT_AERO as f32; npix];
    let mut teps = vec![DEFAULT_EPS as f32; npix];
    let mut sband: Vec<Vec<f32>> = (0..nbands).map(|_| vec![0.0f32; npix]).collect();

    // Set fill flags
    for pix in 0..npix {
        if is_fill_pixel(qaband[pix]) {
            ipflag[pix] = 1u8 << IPFLAG_FILL;
            for ib in 0..nbands {
                sband[ib][pix] = -9999.0; // fill value in sband
            }
        }
    }

    // Initial atmospheric correction for each band
    for ib in 0..nbands {
        let (tgo_x_roatm, tgo_x_ttatmg, satm_val) = if ib == DNS_BAND9 {
            // B09 water vapor band: bypass atmospheric correction
            (1.0f32 * 0.0f32, 1.0f32 * 1.0f32, 0.0f32)
        } else {
            (btgo[ib] as f32 * broatm[ib] as f32,
             btgo[ib] as f32 * bttatmg[ib] as f32,
             bsatm[ib] as f32)
        };

        for pix in 0..npix {
            if is_fill_pixel(qaband[pix]) {
                continue;
            }
            let toa = toa_bands[ib][(pix / nsamps, pix % nsamps)];
            // C: roslamb = toaband[ib][i] - tgo_x_roatm;
            //    roslamb /= tgo_x_ttatmg + satm * roslamb;
            let mut roslamb = toa - tgo_x_roatm;
            roslamb /= tgo_x_ttatmg + satm_val * roslamb;
            sband[ib][pix] = roslamb;
        }
    }
```

Note: Sentinel TOA bands are already unscaled reflectance (not divided by cos(SZA)). The C code uses `toaband[ib][i]` directly, not `toaband[ib][i] / cos(sza)`.

- [ ] **Step 5: Implement aerosol retrieval loop**

The Sentinel aerosol loop starts at (0,0) and steps by SAERO_WINDOW=6. For each window, it averages TOA within the 6x6 window, computes NDWI from climatological sband, and runs eps optimization.

Key differences from Landsat:
- Loops `i in (0..nlines).step_by(6)`, `j in (0..nsamps).step_by(6)`
- Averages TOA within window for bands B01, B02, B04, B12
- NDWI: `(sband[B8A] - sband[B12]*0.5) / (sband[B8A] + sband[B12]*0.5)` at curr_pix
- Eps optimization: quadratic fit on residual1/2/3, NOT via subaeroret_new refit
- NDVI validation: averages TOA in window for B8A and B04, then runs atmcorlamb2_new
- Water retest: averages TOA in window for B01/B04/B8A/B12
- Water failure: sets `IPFLAG_FAILED` (Sentinel), not just clearing bits (Landsat)
- Copies taero/teps to all pixels in window

This is a large function — implement it matching the C code's exact logic from `compute_sentinel_refl.c:1051-1726`.

```rust
    // Eps quadratic fit constants
    let eps1 = LOW_EPS;
    let eps2 = MOD_EPS;
    let eps3 = HIGH_EPS;
    let xa = eps1 * eps1 - eps3 * eps3;
    let xd = ((eps2 * eps2 - eps3 * eps3) as i32) as f64;
    let xb = eps1 - eps3;
    let xe = eps2 - eps3;

    let tth = &SENTINEL_TTH[..nbands];

    let iband1 = DNS_BAND4; // red band is reference

    // Iterate over window UL corners
    for i in (0..nlines).step_by(aero_window) {
        for j in (0..nsamps).step_by(aero_window) {
            let curr_pix = i * nsamps + j;

            if is_fill_pixel(qaband[curr_pix]) {
                ipflag[curr_pix] = 1u8 << IPFLAG_FILL;
                continue;
            }

            // Per-pixel geolocation for CMG lookup
            let (pixel_lat, pixel_lon) = utm_to_deg(space_def, i as i32, j as i32);
            let cmg = latlon_to_cmg(pixel_lat, pixel_lon);

            // ... (CMG slope/intercept computation — same pattern as Landsat)
            // ... (NDWI from sband[DNS_BAND8A] and sband[DNS_BAND12])
            // ... (erelc setup: B01, B02, B04=1.0, B12)

            // Average TOA within 6x6 window for bands B01, B02, B04, B12
            let mut troatm = vec![0.0f64; nbands];
            let mut pix_count = 0usize;
            for iline in i..(i + aero_window).min(nlines) {
                for isamp in j..(j + aero_window).min(nsamps) {
                    let wp = iline * nsamps + isamp;
                    if is_fill_pixel(qaband[wp]) { continue; }
                    troatm[DNS_BAND1] += toa_bands[DNS_BAND1][(iline, isamp)] as f64;
                    troatm[DNS_BAND2] += toa_bands[DNS_BAND2][(iline, isamp)] as f64;
                    troatm[DNS_BAND4] += toa_bands[DNS_BAND4][(iline, isamp)] as f64;
                    troatm[DNS_BAND12] += toa_bands[DNS_BAND12][(iline, isamp)] as f64;
                    pix_count += 1;
                }
            }
            if pix_count == 0 { continue; }
            troatm[DNS_BAND1] /= pix_count as f64;
            troatm[DNS_BAND2] /= pix_count as f64;
            troatm[DNS_BAND4] /= pix_count as f64;
            troatm[DNS_BAND12] /= pix_count as f64;

            // Run subaeroret_new at 3 eps values
            // ... (same pattern as Landsat but with Sentinel bands)

            // Eps quadratic optimization (C lines 1430-1454)
            // NOTE: This is different from Landsat — Sentinel uses quadratic
            // on residuals directly, with additional resepsmin check
            let xc = residual1 - residual3;
            let xf_val = residual2 - residual3;
            let coefa = (xc * xe - xb * xf_val) / (xa * xe - xb * xd);
            let coefb = (xa * xf_val - xc * xd) / (xa * xe - xb * xd);
            let epsmin = -coefb / (2.0 * coefa);
            let resepsmin = xa * epsmin * epsmin + xb * epsmin + xc;

            let eps = if epsmin < LOW_EPS || epsmin > HIGH_EPS {
                if residual1 < residual3 { eps1 } else { eps3 }
            } else if resepsmin > residual1 || resepsmin > residual3 {
                if residual1 < residual3 { eps1 } else { eps3 }
            } else {
                epsmin
            };

            // Run final subaeroret_new at chosen eps
            // ... (get raot and residual)

            teps[curr_pix] = eps as f32;
            taero[curr_pix] = raot as f32;
            let corf = raot / xmus;

            // Residual test (C line 1484)
            if residual < (0.015 + 0.005 * corf + 0.10 * troatm[DNS_BAND12]) {
                // NDVI validation using NxN-averaged TOA for B8A and B04
                // ... (average TOA in window, run atmcorlamb2_new)
                // if ros5 > 0.1 && (ros5-ros4)/(ros5+ros4) > 0: VALID
                // else: WATER
            } else {
                ipflag[curr_pix] = 1u8 << IPFLAG_WATER;
            }

            // Water retest (C lines 1595-1708)
            if ipflag[curr_pix] & (1u8 << IPFLAG_WATER) != 0 {
                // Average TOA for B01/B04/B8A/B12, erelc all 1.0, eps=1.5
                // ... (run subaeroret_new as water)
                // if residual > (0.010 + 0.005*corf) || ros1 < 0: FAILED
                // else: WATER | VALID
            }

            // Copy taero/teps to all pixels in window (C lines 1712-1724)
            for iline in i..(i + aero_window).min(nlines) {
                for isamp in j..(j + aero_window).min(nsamps) {
                    let wp = iline * nsamps + isamp;
                    if is_fill_pixel(qaband[wp]) { continue; }
                    teps[wp] = teps[curr_pix];
                    taero[wp] = taero[curr_pix];
                }
            }
        }
    }
```

The implementer should fill in the `...` sections following the exact C code patterns. The code snippets above show the structure; every `...` is a well-defined block from the C reference (lines noted in comments).

- [ ] **Step 6: Implement aerosol post-processing (3 calls)**

```rust
    // Sentinel aerosol post-processing (replaces Landsat fix_invalid + interp)
    // C lines 1770-1806
    aerosol_interp_sentinel(
        aero_window, &qaband, &mut ipflag, &mut taero, nlines, nsamps,
    );
    ipflag_expand_failed_sentinel(&mut ipflag, nlines, nsamps);
    aero_avg_failed_sentinel(
        &qaband, &mut ipflag, &mut taero, &mut teps, nlines, nsamps,
    );
```

- [ ] **Step 7: Implement final per-pixel correction**

```rust
    // Final atmospheric correction (C lines 1824-1915)
    let mut sr_f32: Vec<Array2<f64>> = (0..nbands)
        .map(|_| Array2::zeros((nlines, nsamps)))
        .collect();

    for ib in 0..nbands {
        // B10: just copy TOA values
        if ib == DNS_BAND10 {
            for pix in 0..npix {
                let (i, j) = (pix / nsamps, pix % nsamps);
                sr_f32[ib][(i, j)] = toa_bands[ib][(i, j)] as f64;
            }
            continue;
        }

        for pix in 0..npix {
            if is_fill_pixel(qaband[pix]) { continue; }
            let (i, j) = (pix / nsamps, pix % nsamps);

            let rotoa = toa_bands[ib][(i, j)] as f64;
            let raot = taero[pix] as f64;
            let eps = teps[pix] as f64;

            let roslamb = atmcorlamb2_new(
                &atm_coeff[ib], tgo_arr[ib], ib,
                raot, normext_p0a3[ib], rotoa, lambda, eps,
            );

            // Aerosol QA on B01 (C line 1884)
            if ib == DNS_BAND1 {
                let rsurf = sband[ib][pix] as f64;
                let tmpf = (rsurf - roslamb).abs();
                if tmpf <= LOW_AERO_THRESH {
                    ipflag[pix] |= 1u8 << AERO1_QA;
                } else if tmpf < AVG_AERO_THRESH {
                    ipflag[pix] |= 1u8 << AERO2_QA;
                } else {
                    ipflag[pix] |= (1u8 << AERO1_QA) | (1u8 << AERO2_QA);
                }
            }

            let roslamb = roslamb.clamp(MIN_VALID_REFL, MAX_VALID_REFL);
            sr_f32[ib][(i, j)] = roslamb;
        }
    }
```

- [ ] **Step 8: Implement output scaling**

Same as Landsat — scale SR to uint16, aerosol to int16, QA as uint8. Reuse the same scaling constants.

```rust
    // Scale to output (same as Landsat)
    let offset_f32 = BAND_OFFSET_REFL as f32;
    let mult_f32 = MULT_FACTOR_REFL as f32;
    let sr_bands: Vec<Array2<u16>> = sr_f32.iter().map(|band| {
        Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
            let pix = i * nsamps + j;
            if is_fill_pixel(qaband[pix]) { 0u16 }
            else {
                let sband_f32 = band[(i, j)] as f32;
                let tmpf = (sband_f32 + offset_f32) * mult_f32;
                tmpf.round().clamp(0.0, u16::MAX as f32) as u16
            }
        })
    }).collect();

    let aero_mult_f32 = MULT_FACTOR_AERO as f32;
    let aerosol = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        let pix = i * nsamps + j;
        if is_fill_pixel(qaband[pix]) { AERO_FILL }
        else {
            let tmpf = taero[pix] * aero_mult_f32;
            tmpf.round().clamp(0.0, 5000.0) as i16
        }
    });

    let qa = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        ipflag[i * nsamps + j]
    });

    SurfaceReflectanceResult {
        sr_bands,
        bt_bands: vec![], // No thermal bands for Sentinel
        aerosol,
        qa,
    }
}
```

- [ ] **Step 9: Add import for new aerosol functions**

Update the imports at the top of `correction.rs`:

```rust
use crate::aerosol::{
    aerosol_interp, fix_invalid_aerosols, subaeroret_new,
    aerosol_interp_sentinel, ipflag_expand_failed_sentinel, aero_avg_failed_sentinel,
};
```

- [ ] **Step 10: Run tests**

Run: `cargo test -p lasrc_core`
Expected: All tests pass. The new function compiles without errors.

- [ ] **Step 11: Commit**

```bash
git add lasrc_core/src/correction.rs
git commit -m "feat: add compute_sentinel_surface_reflectance orchestrator"
```

---

### Task 5: Update PyO3 bindings for Sentinel

**Files:**
- Modify: `lasrc_py/src/lib.rs`
- Test: `maturin develop --release -m lasrc_py/Cargo.toml`

- [ ] **Step 1: Add a Sentinel-specific processing function**

Add a new `process_sentinel_surface_reflectance` PyO3 function that takes scalar angles instead of per-pixel grids, and no qa_band or bt_bands:

```rust
use lasrc_core::correction::{
    AuxiliaryData, compute_surface_reflectance, compute_sentinel_surface_reflectance,
};

/// Compute surface reflectance for a Sentinel-2 scene.
///
/// Parameters
/// ----------
/// sensor_name : str
///     One of "SENTINEL_2A", "SENTINEL_2B", "SENTINEL_2C"
/// toa_bands : list of 2-D float32 arrays
///     TOA reflectance per band (13 bands, all resampled to 10m).
/// solar_zenith, solar_azimuth, view_zenith, view_azimuth : float
///     Scene-center mean angles in degrees.
/// lut : PyLookupTables
/// aux : PyAuxiliaryData
/// ul_corner_x, ul_corner_y, pixel_size_x, pixel_size_y : float
/// utm_zone : int
///
/// Returns dict with sr_bands, bt_bands (empty), aerosol, qa.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn process_sentinel_surface_reflectance<'py>(
    py: Python<'py>,
    sensor_name: &str,
    toa_bands: Vec<PyReadonlyArray2<'py, f32>>,
    solar_zenith: f64,
    solar_azimuth: f64,
    view_zenith: f64,
    view_azimuth: f64,
    lut: &PyLookupTables,
    aux: &PyAuxiliaryData,
    ul_corner_x: f64,
    ul_corner_y: f64,
    pixel_size_x: f64,
    pixel_size_y: f64,
    utm_zone: i32,
) -> PyResult<PyObject> {
    let sensor: Box<dyn Sensor> = match sensor_name {
        "SENTINEL_2A" => Box::new(Sentinel2A),
        "SENTINEL_2B" => Box::new(Sentinel2B),
        "SENTINEL_2C" => Box::new(Sentinel2C),
        _ => return Err(PyValueError::new_err(format!(
            "Unknown Sentinel sensor: {sensor_name}"
        ))),
    };

    let toa_views: Vec<_> = toa_bands.iter().map(|a| a.as_array()).collect();

    let space_def = SpaceDef {
        ul_corner_x,
        ul_corner_y,
        pixel_size: [pixel_size_x, pixel_size_y],
        zone: utm_zone,
    };

    let result = compute_sentinel_surface_reflectance(
        sensor.as_ref(),
        &toa_views,
        solar_zenith,
        solar_azimuth,
        view_zenith,
        view_azimuth,
        &lut.inner,
        &aux.inner,
        &space_def,
    );

    // Build output dict (same pattern as Landsat)
    let dict = pyo3::types::PyDict::new(py);

    let sr_list: Vec<Bound<'py, PyArray1<u16>>> = result.sr_bands.into_iter().map(|a| {
        let (v, _) = a.into_raw_vec_and_offset();
        v.into_pyarray(py)
    }).collect();
    dict.set_item("sr_bands", sr_list)?;
    dict.set_item("bt_bands", Vec::<Bound<'py, PyArray1<u16>>>::new())?;

    let (aerosol_vec, _) = result.aerosol.into_raw_vec_and_offset();
    dict.set_item("aerosol", aerosol_vec.into_pyarray(py))?;
    let (qa_vec, _) = result.qa.into_raw_vec_and_offset();
    dict.set_item("qa", qa_vec.into_pyarray(py))?;

    Ok(dict.into())
}
```

- [ ] **Step 2: Register the new function in the module**

```rust
#[pymodule]
fn lasrc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;
    m.add_class::<PyLookupTables>()?;
    m.add_class::<PyAuxiliaryData>()?;
    m.add_function(wrap_pyfunction!(process_surface_reflectance, m)?)?;
    m.add_function(wrap_pyfunction!(process_sentinel_surface_reflectance, m)?)?;
    Ok(())
}
```

- [ ] **Step 3: Build**

Run: `source .venv/bin/activate && maturin develop --release -m lasrc_py/Cargo.toml`
Expected: Build succeeds.

- [ ] **Step 4: Commit**

```bash
git add lasrc_py/src/lib.rs
git commit -m "feat: add PyO3 bindings for Sentinel-2 surface reflectance"
```

---

### Task 6: Create Python SAFE I/O module

**Files:**
- Create: `lasrc_py/python/lasrc/io_sentinel.py`
- Test: Manual — `python -c "from lasrc.io_sentinel import read_sentinel_safe; ..."`

- [ ] **Step 1: Create `io_sentinel.py`**

This module reads a Sentinel-2 SAFE archive: JP2 bands at native resolution, resamples 20m/60m to 10m using nearest-neighbor, extracts metadata from XML files.

```python
"""Sentinel-2 SAFE archive I/O: read JP2 bands, extract metadata, resample to 10m."""

from pathlib import Path
from xml.etree import ElementTree

import numpy as np


# Band ordering matching the 13-band Rust layout
# (B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12)
SENTINEL_BAND_ORDER = [
    "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08",
    "B8A", "B09", "B10", "B11", "B12",
]

# Native resolutions (meters)
BAND_RESOLUTION = {
    "B01": 60, "B02": 10, "B03": 10, "B04": 10,
    "B05": 20, "B06": 20, "B07": 20, "B08": 10,
    "B8A": 20, "B09": 60, "B10": 60, "B11": 20, "B12": 20,
}


def read_sentinel_safe(safe_dir: str | Path) -> dict:
    """Read a Sentinel-2 SAFE archive and return TOA bands + metadata.

    Parameters
    ----------
    safe_dir : path
        Path to the .SAFE directory.

    Returns
    -------
    dict with keys:
        toa_bands : list of 13 numpy arrays (float32, 10980x10980)
        angles : dict with solar_zenith, solar_azimuth, view_zenith, view_azimuth (scalars)
        profile : dict with epsg, ul_x, ul_y, pixel_size, nlines, nsamps
    """
    import rasterio

    safe_dir = Path(safe_dir)

    # Parse metadata
    mtd_msil1c = safe_dir / "MTD_MSIL1C.xml"
    mtd_tl = safe_dir / "MTD_TL.xml"
    angles = _parse_angles(mtd_tl if mtd_tl.exists() else mtd_msil1c)
    quantification_value, radiometric_offset = _parse_quantification(mtd_msil1c)

    # Read bands
    toa_bands = []
    nlines_10m = None
    nsamps_10m = None
    profile = {}

    for band_name in SENTINEL_BAND_ORDER:
        # Find JP2 file
        jp2_files = list(safe_dir.glob(f"*_{band_name}.jp2"))
        if not jp2_files:
            raise FileNotFoundError(f"No JP2 file found for {band_name} in {safe_dir}")

        with rasterio.open(jp2_files[0]) as src:
            dn = src.read(1).astype(np.float32)

            # Get 10m reference dimensions from B02
            if band_name == "B02":
                nlines_10m = src.height
                nsamps_10m = src.width
                transform = src.transform
                crs = src.crs
                profile = {
                    "epsg": crs.to_epsg(),
                    "ul_x": transform.c,
                    "ul_y": transform.f,
                    "pixel_size_x": transform.a,
                    "pixel_size_y": abs(transform.e),
                    "nlines": nlines_10m,
                    "nsamps": nsamps_10m,
                }

        # Unscale to TOA reflectance
        # Since Baseline 4.00: toa = (DN + offset) / quantification_value
        toa = (dn + radiometric_offset) / quantification_value

        # Resample to 10m if needed
        native_res = BAND_RESOLUTION[band_name]
        if native_res != 10:
            toa = _resample_to_10m(toa, native_res, nlines_10m, nsamps_10m)

        toa_bands.append(toa)

    return {
        "toa_bands": toa_bands,
        "angles": angles,
        "profile": profile,
    }


def _resample_to_10m(
    data: np.ndarray, native_res: int, nlines_10m: int, nsamps_10m: int
) -> np.ndarray:
    """Nearest-neighbor resample from 20m or 60m to 10m.

    Matches C code's convert_to_10m: each native pixel maps to a
    (native_res/10) x (native_res/10) block of 10m pixels.
    """
    scale = native_res // 10
    # np.repeat along both axes = nearest-neighbor upsampling
    resampled = np.repeat(np.repeat(data, scale, axis=0), scale, axis=1)
    # Trim to exact 10m dimensions (handles edge cases)
    return resampled[:nlines_10m, :nsamps_10m]


def _parse_angles(xml_path: Path) -> dict:
    """Extract scene-center mean angles from MTD_TL.xml or MTD_MSIL1C.xml."""
    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    # Remove namespace prefix for easier searching
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"

    # Try MTD_TL.xml format first (Tile-level metadata)
    # Mean_Sun_Angle and Mean_Viewing_Incidence_Angle
    angles = {}

    # Sun angles
    sun_el = root.find(f".//{ns}Mean_Sun_Angle/{ns}ZENITH_ANGLE")
    if sun_el is None:
        sun_el = root.find(".//Mean_Sun_Angle/ZENITH_ANGLE")
    if sun_el is not None:
        angles["solar_zenith"] = float(sun_el.text)
    else:
        angles["solar_zenith"] = 30.0  # fallback

    sun_az = root.find(f".//{ns}Mean_Sun_Angle/{ns}AZIMUTH_ANGLE")
    if sun_az is None:
        sun_az = root.find(".//Mean_Sun_Angle/AZIMUTH_ANGLE")
    if sun_az is not None:
        angles["solar_azimuth"] = float(sun_az.text)
    else:
        angles["solar_azimuth"] = 150.0

    # View angles — average across all bands
    view_zen_els = root.findall(f".//{ns}Mean_Viewing_Incidence_Angle/{ns}ZENITH_ANGLE")
    if not view_zen_els:
        view_zen_els = root.findall(".//Mean_Viewing_Incidence_Angle/ZENITH_ANGLE")
    if view_zen_els:
        angles["view_zenith"] = np.mean([float(e.text) for e in view_zen_els])
    else:
        angles["view_zenith"] = 0.0

    view_az_els = root.findall(f".//{ns}Mean_Viewing_Incidence_Angle/{ns}AZIMUTH_ANGLE")
    if not view_az_els:
        view_az_els = root.findall(".//Mean_Viewing_Incidence_Angle/AZIMUTH_ANGLE")
    if view_az_els:
        angles["view_azimuth"] = np.mean([float(e.text) for e in view_az_els])
    else:
        angles["view_azimuth"] = 0.0

    return angles


def _parse_quantification(xml_path: Path) -> tuple[float, float]:
    """Extract quantification value and radiometric offset from MTD_MSIL1C.xml.

    Returns (quantification_value, radiometric_offset).
    """
    tree = ElementTree.parse(xml_path)
    root = tree.getroot()

    # Quantification value
    qv_el = root.find(".//QUANTIFICATION_VALUE")
    if qv_el is None:
        qv_el = root.find(".//{*}QUANTIFICATION_VALUE")
    quantification_value = float(qv_el.text) if qv_el is not None else 10000.0

    # Radiometric offset (Baseline 4.00+)
    offset_el = root.find(".//RADIO_ADD_OFFSET")
    if offset_el is None:
        offset_el = root.find(".//{*}RADIO_ADD_OFFSET")
    radiometric_offset = float(offset_el.text) if offset_el is not None else -1000.0

    return quantification_value, radiometric_offset
```

- [ ] **Step 2: Verify it can read the test SAFE archive**

Run:
```bash
source .venv/bin/activate
python -c "
from lasrc.io_sentinel import read_sentinel_safe
data = read_sentinel_safe('test_data/S2B_MSIL1C_20260124T073109_N0511_R049_T38PNC_20260124T092241.SAFE')
print(f'Bands: {len(data[\"toa_bands\"])}')
for i, b in enumerate(data['toa_bands']):
    print(f'  Band {i}: shape={b.shape}, min={b.min():.4f}, max={b.max():.4f}')
print(f'Angles: {data[\"angles\"]}')
print(f'Profile: {data[\"profile\"]}')
"
```
Expected: 13 bands, all 10980x10980, reasonable TOA values (0-1.6), valid angles.

- [ ] **Step 3: Commit**

```bash
git add lasrc_py/python/lasrc/io_sentinel.py
git commit -m "feat: add Sentinel-2 SAFE archive I/O with JP2 reading and 10m resampling"
```

---

### Task 7: Update sensor configs and create Sentinel pipeline

**Files:**
- Modify: `lasrc_py/python/lasrc/sensors.py`
- Create: `lasrc_py/python/lasrc/pipeline_sentinel.py`
- Test: Integration test via `run_test_sentinel.py`

- [ ] **Step 1: Update sensors.py with 13-band Sentinel configs**

Add SENTINEL_2B and SENTINEL_2C configs, and update SENTINEL_2A to 13 bands:

```python
"SENTINEL_2A": {
    "name": "SENTINEL_2A",
    "satellite": "SENTINEL-2A",
    "instrument": "MSI",
    "refl_band_names": [
        "sr_band1", "sr_band2", "sr_band3", "sr_band4",
        "sr_band5", "sr_band6", "sr_band7", "sr_band8",
        "sr_band8a", "sr_band9", "sr_band10", "sr_band11", "sr_band12",
    ],
    "thm_band_names": [],
    "input_band_names": [
        "B01", "B02", "B03", "B04", "B05", "B06",
        "B07", "B08", "B8A", "B09", "B10", "B11", "B12",
    ],
    "resolution_m": 10.0,
    "nsr_bands": 13,
},
"SENTINEL_2B": {
    "name": "SENTINEL_2B",
    "satellite": "SENTINEL-2B",
    "instrument": "MSI",
    "refl_band_names": [
        "sr_band1", "sr_band2", "sr_band3", "sr_band4",
        "sr_band5", "sr_band6", "sr_band7", "sr_band8",
        "sr_band8a", "sr_band9", "sr_band10", "sr_band11", "sr_band12",
    ],
    "thm_band_names": [],
    "input_band_names": [
        "B01", "B02", "B03", "B04", "B05", "B06",
        "B07", "B08", "B8A", "B09", "B10", "B11", "B12",
    ],
    "resolution_m": 10.0,
    "nsr_bands": 13,
},
"SENTINEL_2C": {
    "name": "SENTINEL_2C",
    "satellite": "SENTINEL-2C",
    "instrument": "MSI",
    "refl_band_names": [
        "sr_band1", "sr_band2", "sr_band3", "sr_band4",
        "sr_band5", "sr_band6", "sr_band7", "sr_band8",
        "sr_band8a", "sr_band9", "sr_band10", "sr_band11", "sr_band12",
    ],
    "thm_band_names": [],
    "input_band_names": [
        "B01", "B02", "B03", "B04", "B05", "B06",
        "B07", "B08", "B8A", "B09", "B10", "B11", "B12",
    ],
    "resolution_m": 10.0,
    "nsr_bands": 13,
},
```

- [ ] **Step 2: Create `pipeline_sentinel.py`**

```python
"""Sentinel-2 processing pipeline: SAFE archive → surface reflectance."""

from pathlib import Path

import numpy as np

from lasrc.aux import load_auxiliary_data, load_dem, load_ratio_file
from lasrc.io_sentinel import read_sentinel_safe
from lasrc.sensors import SENSORS

import lasrc as _lasrc


def process_sentinel_scene(
    safe_dir: str | Path,
    lut_dir: str | Path,
    viirs_aux_path: str | Path,
    dem_path: str | Path,
    ratio_path: str | Path,
    output_path: str | Path,
    sensor_name: str = "SENTINEL_2B",
    output_format: str = "espa",
) -> None:
    """Process a Sentinel-2 SAFE archive to surface reflectance.

    Parameters
    ----------
    safe_dir : path
        Path to the .SAFE directory.
    lut_dir : path
        Path to directory containing MSI LUT files (ANGLE_NEW.hdf, etc.)
    viirs_aux_path : path
        Path to VIIRS auxiliary HDF5 file.
    dem_path : path
        Path to CMGDEM HDF4 file.
    ratio_path : path
        Path to band ratio HDF4 file.
    output_path : path
        Output directory (ESPA) or file (COG).
    sensor_name : str
        One of SENTINEL_2A, SENTINEL_2B, SENTINEL_2C.
    output_format : str
        "espa" for flat binary .img files.
    """
    sensor_config = SENSORS[sensor_name]
    nsr_bands = sensor_config["nsr_bands"]  # 13

    # Read SAFE archive
    scene = read_sentinel_safe(safe_dir)
    angles = scene["angles"]
    profile = scene["profile"]
    nlines = profile["nlines"]
    nsamps = profile["nsamps"]

    # Load LUTs (HDF4)
    lut_dir = Path(lut_dir)
    lut_data = _load_sentinel_luts(lut_dir, nsr_bands)
    lut = _lasrc.PyLookupTables(**lut_data)

    # Load auxiliary data
    aux_data = load_auxiliary_data(str(viirs_aux_path), aux_source="VIIRS")
    dem = load_dem(str(dem_path))
    ratios = load_ratio_file(str(ratio_path))
    aux = _lasrc.PyAuxiliaryData(dem=dem, **aux_data, **ratios)

    # Extract UTM zone from EPSG
    epsg = profile["epsg"]
    utm_zone = epsg % 100
    if epsg > 32700:
        utm_zone = -utm_zone

    # Run Rust correction
    result = _lasrc.process_sentinel_surface_reflectance(
        sensor_name=sensor_name,
        toa_bands=scene["toa_bands"],
        solar_zenith=angles["solar_zenith"],
        solar_azimuth=angles["solar_azimuth"],
        view_zenith=angles["view_zenith"],
        view_azimuth=angles["view_azimuth"],
        lut=lut,
        aux=aux,
        ul_corner_x=profile["ul_x"],
        ul_corner_y=profile["ul_y"],
        pixel_size_x=profile["pixel_size_x"],
        pixel_size_y=profile["pixel_size_y"],
        utm_zone=utm_zone,
    )

    # Reshape and write output
    for i in range(len(result["sr_bands"])):
        result["sr_bands"][i] = result["sr_bands"][i].reshape(nlines, nsamps)
    result["aerosol"] = result["aerosol"].reshape(nlines, nsamps)
    result["qa"] = result["qa"].reshape(nlines, nsamps)

    if output_format == "espa":
        _write_espa_sentinel(output_path, result, sensor_config)


def _write_espa_sentinel(output_dir, result, sensor_config):
    """Write Sentinel SR output in ESPA format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, name in enumerate(sensor_config["refl_band_names"]):
        result["sr_bands"][i].astype(np.uint16).tofile(output_dir / f"{name}.img")

    result["aerosol"].astype(np.int16).tofile(output_dir / "sr_aerosol.img")
    result["qa"].astype(np.uint8).tofile(output_dir / "sr_aerosol_qa.img")


def _load_sentinel_luts(lut_dir: Path, nsr_bands: int) -> dict:
    """Load Sentinel-2 LUTs from HDF4 files.

    The MSI LUTs use band names (NRLUT_BAND_1, NRLUT_BAND_8a, etc.)
    instead of sequential numbers.
    """
    from pyhdf.SD import SD, SDC

    NSOLAR = 8000
    NAOT = 22
    NPRES = 7

    # Band name mapping for RES_LUT dataset names
    # Order must match the 13-band Rust layout
    band_ds_names = [
        "NRLUT_BAND_1", "NRLUT_BAND_2", "NRLUT_BAND_3", "NRLUT_BAND_4",
        "NRLUT_BAND_5", "NRLUT_BAND_6", "NRLUT_BAND_7", "NRLUT_BAND_8",
        "NRLUT_BAND_8a", "NRLUT_BAND_9", "NRLUT_BAND_10", "NRLUT_BAND_11",
        "NRLUT_BAND_12",
    ]

    # Angle LUT
    hdf = SD(str(lut_dir / "ANGLE_NEW.hdf"), SDC.READ)
    angle_data = {}
    for name in ["TSMAX", "TSMIN", "NBFIC", "TTV", "TTS"]:
        angle_data[name.lower()] = hdf.select(name).get().astype(np.float64).ravel().tolist()
    angle_data["nbfi"] = hdf.select("NBFI").get().astype(np.int32).ravel().tolist()
    angle_data["indts"] = hdf.select("INDTS").get().astype(np.int32).ravel().tolist()
    hdf.end()

    # RES LUT (rolutt)
    hdf = SD(str(lut_dir / "RES_LUT_V3.0-URBANCLEAN-V3.0.hdf"), SDC.READ)
    rolutt = np.zeros(nsr_bands * NPRES * NAOT * NSOLAR, dtype=np.float64)
    for ib, ds_name in enumerate(band_ds_names):
        data = hdf.select(ds_name).get().astype(np.float64)  # [NSOLAR, NAOT, NPRES]
        for ip in range(NPRES):
            for ia in range(NAOT):
                base = ib * NPRES * NAOT * NSOLAR + ip * NAOT * NSOLAR + ia * NSOLAR
                rolutt[base:base + NSOLAR] = data[:, ia, ip]
    hdf.end()

    # AERO LUT (sphalbt, normext) — ASCII format
    sphalbt, normext = _load_aero_ascii(
        lut_dir / "AERO_LUT_V3.0-URBANCLEAN-V3.0.ASCII", nsr_bands
    )

    # TRANS LUT — ASCII format
    transt = _load_trans_ascii(
        lut_dir / "TRANS_LUT_V3.0-URBANCLEAN-V3.0.ASCII", nsr_bands
    )

    return {
        "rolutt": rolutt.tolist(),
        "transt": transt,
        "sphalbt": sphalbt,
        "normext": normext,
        **angle_data,
        "nsr_bands": nsr_bands,
    }


def _load_aero_ascii(path: Path, nsr_bands: int) -> tuple[list, list]:
    """Load sphalbt and normext from ASCII file. Same format as Landsat."""
    NPRES = 7
    NAOT = 22

    sphalbt = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)
    normext = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        idx += 1  # band header
        for ip in range(NPRES):
            idx += 1  # pressure header
            for ia in range(NAOT):
                parts = lines[idx].split()
                base = ib * NPRES * NAOT + ip * NAOT + ia
                sphalbt[base] = float(parts[1])
                normext[base] = float(parts[2])
                idx += 1

    return sphalbt.tolist(), normext.tolist()


def _load_trans_ascii(path: Path, nsr_bands: int) -> list:
    """Load transmission LUT from ASCII file. Same format as Landsat."""
    NPRES = 7
    NSUNANGLE = 22
    NSUNANGLE_FILE = 21
    NAOT = 22

    transt = np.zeros(nsr_bands * NPRES * NAOT * NSUNANGLE, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        idx += 1  # band header
        for ip in range(NPRES):
            idx += 1  # pressure header
            for isun in range(NSUNANGLE_FILE):
                parts = lines[idx].split()
                for ia in range(NAOT):
                    base = (ib * NPRES * NAOT * NSUNANGLE
                            + ip * NAOT * NSUNANGLE
                            + isun + ia * NSUNANGLE)
                    transt[base] = float(parts[1 + ia])
                idx += 1

    return transt.tolist()
```

- [ ] **Step 3: Commit**

```bash
git add lasrc_py/python/lasrc/sensors.py lasrc_py/python/lasrc/pipeline_sentinel.py
git commit -m "feat: add Sentinel-2 pipeline and update sensor configs to 13 bands"
```

---

### Task 8: Create test runner and extend comparison script

**Files:**
- Create: `lasrc_py/run_test_sentinel.py`
- Modify: `lasrc_py/compare_outputs.py`
- Test: Run on test granule and compare against C reference

- [ ] **Step 1: Create `run_test_sentinel.py`**

```python
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
```

- [ ] **Step 2: Extend `compare_outputs.py` with `--sensor sentinel` mode**

Add Sentinel-2 band definitions and dimensions. The comparison should support:
- 13 SR bands (sr_band1..sr_band8, sr_band8a, sr_band9, sr_band10, sr_band11, sr_band12) as uint16
- sr_aerosol as int16
- sr_aerosol_qa as uint8
- Dimensions: 10980x10980

Add a `--sensor` argument that selects between Landsat and Sentinel band lists/dimensions:

```python
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
```

Update `main()` to accept `--sensor landsat|sentinel` and select the appropriate band list and dimensions.

- [ ] **Step 3: Run the pipeline on the test granule**

Run:
```bash
source .venv/bin/activate
maturin develop --release -m lasrc_py/Cargo.toml
python lasrc_py/run_test_sentinel.py --format espa
```
Expected: Pipeline completes, output written to `test_data/output_espa_s2/`.

- [ ] **Step 4: Compare against C reference**

Run:
```bash
python lasrc_py/compare_outputs.py --sensor sentinel \
    --rust-dir test_data/output_espa_s2 \
    --c-dir test_data/c_output \
    --c-prefix S2B_MSI_L1C_T38PNC_20260124_20260124_
```
Expected: Comparison table for all 15 files. Initial run may have diffs — iterate on fixes.

- [ ] **Step 5: Commit**

```bash
git add lasrc_py/run_test_sentinel.py lasrc_py/compare_outputs.py
git commit -m "feat: add Sentinel-2 test runner and extend comparison script"
```

---

## Verification and Iteration

After all 8 tasks, the comparison should show results matching the acceptance criteria from the spec:
- **Median diff**: 1-2 scaled integer units
- **P95 diff**: under 5 units
- **P99 diff**: under 30 units
- **Aerosol**: under 7% of pixels differing
- **Aerosol QA**: under 1% of pixels differing

If diffs are larger than expected, investigate using the same approach as the Landsat port:
1. Check f32/f64 precision mismatches
2. Verify band index mapping (the 13-band layout shifts many indices)
3. Compare TOA values (input reading correctness)
4. Check fill detection (ANY band vs specific band)
5. Verify the eps optimization logic (quadratic fit differs from Landsat)
6. Check aerosol post-processing function outputs against C intermediate files

The C code supports `WRITE_TAERO` for dumping intermediate aerosol files (`ipflag.img`, `aerosols.img`, etc.) — use these for debugging if available.
