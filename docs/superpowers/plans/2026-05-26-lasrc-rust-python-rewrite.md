# LaSRC Rust+Python Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite the LaSRC atmospheric correction C code as a Rust core library with Python orchestration, matching C output for validation.

**Architecture:** Rust `lasrc_core` crate implements compute kernels (LUT interpolation, aerosol retrieval, atmospheric correction). A `lasrc_py` crate wraps it via PyO3. Python handles I/O (rasterio, h5py), auxiliary data loading, and CLI orchestration.

**Tech Stack:** Rust (ndarray, rayon, PyO3), Python (rasterio, numpy, h5py, pyhdf, click, maturin)

---

## File Structure

### Rust: `lasrc_core/`
- `Cargo.toml` - Crate config, depends on ndarray and rayon
- `src/lib.rs` - Public API re-exports
- `src/constants.rs` - All physical constants, LUT dimensions, thresholds
- `src/sensor.rs` - `Sensor` trait + `Landsat8`, `Landsat9`, `Sentinel2A/B/C` impls
- `src/geometry.rs` - UTM-to-degree conversion, scattering angle computation
- `src/utils.rs` - Quick select (median), polynomial fitting
- `src/lut.rs` - LookupTable struct, multi-dimensional interpolation
- `src/gas_transmission.rs` - Ozone, water vapor, other gas transmittance
- `src/rayleigh.rs` - Rayleigh (molecular) scattering reflectance
- `src/atmospheric.rs` - `atmcorlamb2` atmospheric correction function
- `src/aerosol.rs` - Aerosol retrieval (`subaeroret`) and interpolation
- `src/correction.rs` - Top-level `compute_surface_reflectance` orchestrator

### Rust: `lasrc_py/`
- `Cargo.toml` - Depends on lasrc_core, pyo3, numpy
- `pyproject.toml` - maturin build config
- `src/lib.rs` - PyO3 module definition, NumPy array conversion

### Python: `lasrc_py/python/lasrc/`
- `__init__.py` - Package init
- `sensors.py` - Sensor configuration dicts (band names, wavelengths, resolutions)
- `aux.py` - Auxiliary data loading (HDF4/5 → NumPy arrays)
- `io.py` - Scene I/O (read Landsat/Sentinel input, write ESPA or COG output)
- `pipeline.py` - Top-level `process_scene()` orchestration
- `cli.py` - Click-based CLI entry point

### Tests
- `lasrc_core/src/*.rs` - Inline `#[cfg(test)]` modules in each Rust file
- `tests/test_geometry.py` - Python-side geometry validation
- `tests/test_lut.py` - LUT interpolation validation
- `tests/test_aerosol.py` - Aerosol retrieval validation
- `tests/test_pipeline.py` - End-to-end pipeline validation

---

## Task 1: Project Scaffolding

**Files:**
- Create: `lasrc_core/Cargo.toml`
- Create: `lasrc_core/src/lib.rs`
- Create: `lasrc_core/src/constants.rs`
- Create: `lasrc_py/Cargo.toml`
- Create: `lasrc_py/pyproject.toml`
- Create: `lasrc_py/src/lib.rs`
- Create: `lasrc_py/python/lasrc/__init__.py`

- [ ] **Step 1: Create lasrc_core Cargo.toml**

```toml
[package]
name = "lasrc_core"
version = "0.1.0"
edition = "2021"
description = "LaSRC atmospheric correction core algorithms"

[dependencies]
ndarray = "0.16"
rayon = "1.10"
```

- [ ] **Step 2: Create lasrc_core/src/lib.rs**

```rust
pub mod constants;
```

- [ ] **Step 3: Create lasrc_core/src/constants.rs with all C constants**

```rust
//! Physical constants, LUT dimensions, and thresholds from the C LaSRC code.

// Version
pub const SR_VERSION: &str = "3.5.1.0 (Collection 2)";

// Angle conversions
pub const DEG2RAD: f64 = 0.017453293;
pub const RAD2DEG: f64 = 57.29577951;

// Atmospheric pressure
pub const ATMOS_PRES_0: f64 = 1013.0;
pub const ONE_DIV_ATMOS_PRES_0: f64 = 0.000987166;
pub const ONE_DIV_8500: f64 = 0.000117647;

// Default aerosol values
pub const DEFAULT_AERO: f64 = 0.05;
pub const DEFAULT_EPS: f64 = 1.5;

// Angstrom coefficient categories
pub const LOW_EPS: f64 = 1.0;
pub const WATER_EPS: f64 = 1.5;
pub const MOD_EPS: f64 = 1.75;
pub const HIGH_EPS: f64 = 2.5;

// Aerosol thresholds
pub const LOW_AERO_THRESH: f64 = 0.015;
pub const AVG_AERO_THRESH: f64 = 0.03;

// Aerosol window sizes
pub const LAERO_WINDOW: usize = 3;
pub const LHALF_AERO_WINDOW: usize = 1;
pub const SAERO_WINDOW: usize = 6;

pub const LFIX_AERO_WINDOW: usize = 15;
pub const LHALF_FIX_AERO_WINDOW: usize = 7;
pub const LMIN_CLEAR_PIX: usize = 4;

pub const SFIX_AERO_WINDOW: usize = 45;
pub const SHALF_FIX_AERO_WINDOW: usize = 22;
pub const SMIN_CLEAR_PIX: usize = 8;

// LUT dimensions
pub const NPRES_VALS: usize = 7;
pub const NAOT_VALS: usize = 22;
pub const NSOLAR_VALS: usize = 8000;
pub const NSUNANGLE_VALS: usize = 22;
pub const NVIEW_ZEN_VALS: usize = 20;
pub const NSOLAR_ZEN_VALS: usize = 22;
pub const NCOEF: usize = 4;

// Derived LUT dimensions
pub const NAOT_X_NSOLAR: usize = NAOT_VALS * NSOLAR_VALS; // 176000
pub const NAOT_X_NSUNANGLE: usize = NAOT_VALS * NSUNANGLE_VALS; // 484

// CMG/DEM/Ratio grid dimensions (0.05 degree resolution)
pub const CMG_NBLAT: usize = 3600;
pub const CMG_NBLON: usize = 7200;
pub const DEM_NBLAT: usize = 3600;
pub const DEM_NBLON: usize = 7200;
pub const RATIO_NBLAT: usize = 3600;
pub const RATIO_NBLON: usize = 7200;

// Band counts
pub const NREFLL_BANDS: usize = 7;
pub const NSRL_BANDS: usize = 8;
pub const NBANDL_THM: usize = 2;

// Reflectance valid range
pub const MIN_VALID_REFL: f64 = -0.2;
pub const MAX_VALID_REFL: f64 = 1.60;
pub const MIN_VALID_TH: f64 = 150.0;
pub const MAX_VALID_TH: f64 = 350.0;

// Output scale factors
pub const SCALE_FACTOR_REFL: f64 = 0.0000275;
pub const OFFSET_REFL: f64 = -0.20;
pub const MULT_FACTOR_REFL: f64 = 36363.636363636;
pub const BAND_OFFSET_REFL: f64 = 0.20;
pub const SCALE_FACTOR_TH: f64 = 0.00341802;
pub const OFFSET_TH: f64 = 149.0;
pub const MULT_FACTOR_TH: f64 = 292.668;
pub const BAND_OFFSET_TH: f64 = 149.0;
pub const SCALE_FACTOR_AERO: f64 = 0.001;
pub const MULT_FACTOR_AERO: f64 = 1000.0;
pub const AERO_FILL: i16 = -9999;

// AOT values at 550nm (22 values)
pub const AOT550NM: [f64; NAOT_VALS] = [
    0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60,
    0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.30,
    2.60, 3.00, 3.50, 4.00, 4.50, 5.00,
];

// Log of AOT values (for atmospheric reflectance interpolation)
pub const LOG_AOT550NM: [f64; NAOT_VALS] = [
    -4.605170186, -2.995732274, -2.302585093, -1.897119985,
    -1.609437912, -1.203972804, -0.916290732, -0.510825624,
    -0.223143551, 0.000000000, 0.182321557, 0.336472237,
    0.470003629, 0.587786665, 0.693157181, 0.832909123,
    0.955511445, 1.098612289, 1.252762969, 1.386294361,
    1.504077397, 1.609437912,
];

// Pressure table (millibars)
pub const TPRES: [f64; NPRES_VALS] = [1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0];

// Rayleigh optical depth per Landsat band
pub const TAURAY_LANDSAT: [f64; NSRL_BANDS] = [
    0.23638, 0.16933, 0.09070, 0.04827, 0.01563, 0.00129, 0.00037, 0.07984,
];

// Band wavelengths (micrometers)
pub const LAMBDA_LANDSAT: [f64; NREFLL_BANDS] = [0.443, 0.480, 0.585, 0.655, 0.865, 1.61, 2.2];
pub const LAMBDA_SENTINEL: [f64; 11] = [
    0.443, 0.490, 0.560, 0.665, 0.705, 0.740, 0.783, 0.842, 0.865, 1.61, 2.19,
];

// LUT angle parameters
pub const XTS_MIN: f64 = 0.0;
pub const XTS_STEP: f64 = 4.0;
pub const XTV_MIN: f64 = 2.84090;
pub const XTV_STEP: f64 = 3.68017;

// QA flag bit positions
pub const IPFLAG_FILL: u8 = 0;
pub const IPFLAG_CLEAR: u8 = 1;
pub const IPFLAG_WATER: u8 = 2;
pub const IPFLAG_FAILED: u8 = 3;
pub const IPFLAG_FIXED: u8 = 4;
pub const IPFLAG_INTERP_WINDOW: u8 = 5;
pub const AERO1_QA: u8 = 6;
pub const AERO2_QA: u8 = 7;

// Sentinel-specific aerosol window constants
pub const HALF_EXPAND_WIN: usize = 12;
pub const HALF_FAILED_WIN: usize = 30;
pub const MIN_VALID_WINDOW_PIX: usize = 20;

// Fill values
pub const INPUT_FILL: u16 = 0;
pub const ANGLE_FILL: f64 = -999.0;
```

- [ ] **Step 4: Create lasrc_py scaffolding**

`lasrc_py/Cargo.toml`:
```toml
[package]
name = "lasrc_py"
version = "0.1.0"
edition = "2021"

[lib]
name = "lasrc"
crate-type = ["cdylib"]

[dependencies]
lasrc_core = { path = "../lasrc_core" }
pyo3 = { version = "0.23", features = ["extension-module"] }
numpy = "0.23"
```

`lasrc_py/pyproject.toml`:
```toml
[build-system]
requires = ["maturin>=1.0,<2.0"]
build-backend = "maturin"

[project]
name = "lasrc"
version = "0.1.0"
requires-python = ">=3.10"
dependencies = [
    "numpy>=1.24",
    "rasterio>=1.3",
    "h5py>=3.8",
    "pyhdf>=0.10",
    "lxml>=4.9",
    "click>=8.1",
]

[tool.maturin]
python-source = "python"
features = ["pyo3/extension-module"]

[project.scripts]
lasrc = "lasrc.cli:main"
```

`lasrc_py/src/lib.rs`:
```rust
use pyo3::prelude::*;

#[pymodule]
fn lasrc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;
    Ok(())
}
```

`lasrc_py/python/lasrc/__init__.py`:
```python
from lasrc.lasrc import __version__

__all__ = ["__version__"]
```

- [ ] **Step 5: Verify build**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo check`
Expected: Compiles successfully

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_py && pip install -e . 2>&1 | tail -5`
Expected: Builds and installs successfully

- [ ] **Step 6: Commit**

```bash
git add lasrc_core/ lasrc_py/
git commit -m "feat: scaffold lasrc_core and lasrc_py crate structure"
```

---

## Task 2: Geometry Utilities

**Files:**
- Create: `lasrc_core/src/geometry.rs`
- Modify: `lasrc_core/src/lib.rs`

Porting `utmtodeg.c` and scattering angle computation from `comproatm` in `lut_subr.c`.

- [ ] **Step 1: Write failing tests**

Add to `lasrc_core/src/geometry.rs`:

```rust
//! UTM-to-degree conversion and scattering angle computation.

use crate::constants::*;

/// Spatial definition for UTM projection.
pub struct SpaceDef {
    pub ul_corner_x: f64,
    pub ul_corner_y: f64,
    pub pixel_size: [f64; 2],
    pub zone: i32,
}

/// Convert UTM image coordinates to latitude/longitude (WGS84).
pub fn utm_to_deg(space_def: &SpaceDef, line: i32, samp: i32) -> (f64, f64) {
    todo!()
}

/// Compute the scattering angle (degrees) from solar and view geometry.
pub fn scattering_angle(xmus: f64, xmuv: f64, cosxfi: f64) -> f64 {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_utm_to_deg_zone_10n() {
        // UTM Zone 10N, a point near San Francisco
        let space_def = SpaceDef {
            ul_corner_x: 545_000.0,
            ul_corner_y: 4_185_000.0,
            pixel_size: [30.0, 30.0],
            zone: 10,
        };
        let (lat, lon) = utm_to_deg(&space_def, 0, 0);
        // Should be approximately 37.79N, 122.40W
        assert!((lat - 37.79).abs() < 0.1, "lat={lat}");
        assert!((lon - (-122.40)).abs() < 0.1, "lon={lon}");
    }

    #[test]
    fn test_utm_to_deg_southern_hemisphere() {
        // Negative zone for southern hemisphere
        let space_def = SpaceDef {
            ul_corner_x: 500_000.0,
            ul_corner_y: 6_200_000.0,
            pixel_size: [30.0, 30.0],
            zone: -23,
        };
        let (lat, lon) = utm_to_deg(&space_def, 0, 0);
        assert!(lat < 0.0, "Expected southern hemisphere, got lat={lat}");
    }

    #[test]
    fn test_utm_to_deg_pixel_offset() {
        let space_def = SpaceDef {
            ul_corner_x: 500_000.0,
            ul_corner_y: 4_000_000.0,
            pixel_size: [30.0, 30.0],
            zone: 11,
        };
        let (lat0, lon0) = utm_to_deg(&space_def, 0, 0);
        let (lat1, lon1) = utm_to_deg(&space_def, 100, 100);
        // Moving down and right: lat should decrease, lon should increase
        assert!(lat1 < lat0, "lat should decrease: {lat1} vs {lat0}");
        assert!(lon1 > lon0, "lon should increase: {lon1} vs {lon0}");
    }

    #[test]
    fn test_scattering_angle_backscatter() {
        // Direct backscatter: sun and view from same direction
        let xmus = 0.5_f64.cos(); // ~60 deg solar zenith
        let xmuv = 0.0_f64.cos(); // nadir view
        let cosxfi = 1.0; // same azimuth
        let sca = scattering_angle(xmus, xmuv, cosxfi);
        // Scattering angle should be close to 180 for backscatter geometry
        assert!(sca > 90.0, "Expected backscatter, got sca={sca}");
    }

    #[test]
    fn test_scattering_angle_range() {
        let sca = scattering_angle(0.7, 0.95, 0.5);
        assert!(sca >= 0.0 && sca <= 180.0, "sca={sca}");
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test geometry`
Expected: FAIL with "not yet implemented"

- [ ] **Step 3: Implement utm_to_deg**

Replace the `todo!()` in `utm_to_deg` with:

```rust
pub fn utm_to_deg(space_def: &SpaceDef, line: i32, samp: i32) -> (f64, f64) {
    let sa: f64 = 6378137.0; // WGS84 semi-major axis
    let inv_flattening: f64 = 298.257223563;
    let false_easting: f64 = 500000.0;
    let false_northing: f64 = 10000000.0;
    let scale_fact: f64 = 0.9996;

    let sb = sa - (sa / inv_flattening);
    let e2 = ((sa * sa - sb * sb) as f64).sqrt() / sb;
    let e2sq = e2 * e2;
    let c = (sa * sa) / sb;

    let x_utm = space_def.ul_corner_x + (samp as f64 * space_def.pixel_size[0]);
    let y_utm = space_def.ul_corner_y - (line as f64 * space_def.pixel_size[1]);

    let mut x = x_utm - false_easting;
    let mut y = y_utm;
    if space_def.zone < 0 {
        y -= false_northing;
    }

    let zone = space_def.zone.unsigned_abs() as f64;
    let central_meridian = zone * 6.0 - 183.0;

    let mut lat = y / (6366197.724 * scale_fact);
    let cos_lat = lat.cos();
    let sqr_cos_lat = cos_lat * cos_lat;

    let v = (c / (1.0 + e2sq * sqr_cos_lat).sqrt()) * scale_fact;
    let a = x / v;

    let a1 = (2.0 * lat).sin();
    let a2 = a1 * sqr_cos_lat;
    let j2 = lat + a1 / 2.0;
    let j4 = (3.0 * j2 + a2) / 4.0;
    let j6 = (5.0 * j4 + a2 * sqr_cos_lat) / 3.0;

    let alpha = 0.75 * e2sq;
    let beta = (5.0 / 3.0) * alpha * alpha;
    let gama = (35.0 / 27.0) * alpha * alpha * alpha;

    let bm = scale_fact * c * (lat - alpha * j2 + beta * j4 - gama * j6);
    let b = (y - bm) / v;

    let epsi = e2sq * a * a / 2.0 * sqr_cos_lat;
    let eps = a * (1.0 - epsi / 3.0);
    let nab = b * (1.0 - epsi) + lat;
    let senoheps = (eps.exp() - (-eps).exp()) / 2.0;
    let delta = (senoheps / nab.cos()).atan();
    let ta0 = (delta.cos() * nab.tan()).atan();

    let lon_deg = delta * RAD2DEG + central_meridian;
    let lat_deg = (lat
        + (1.0 + e2sq * sqr_cos_lat
            - 1.5 * e2sq * lat.sin() * cos_lat * (ta0 - lat))
            * (ta0 - lat))
        * RAD2DEG;

    (lat_deg, lon_deg)
}
```

- [ ] **Step 4: Implement scattering_angle**

Replace the `todo!()` in `scattering_angle` with:

```rust
pub fn scattering_angle(xmus: f64, xmuv: f64, cosxfi: f64) -> f64 {
    let cscaa = -xmus * xmuv
        - cosxfi * (1.0 - xmus * xmus).sqrt() * (1.0 - xmuv * xmuv).sqrt();
    cscaa.clamp(-1.0, 1.0).acos() * RAD2DEG
}
```

- [ ] **Step 5: Update lib.rs**

```rust
pub mod constants;
pub mod geometry;
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test geometry`
Expected: All 5 tests PASS

- [ ] **Step 7: Commit**

```bash
git add lasrc_core/src/geometry.rs lasrc_core/src/lib.rs
git commit -m "feat: implement UTM-to-degree and scattering angle geometry"
```

---

## Task 3: Utility Functions

**Files:**
- Create: `lasrc_core/src/utils.rs`
- Modify: `lasrc_core/src/lib.rs`

Porting `quick_select.c` and `poly_coeff.c`.

- [ ] **Step 1: Write failing tests**

Add to `lasrc_core/src/utils.rs`:

```rust
//! Quick select median and polynomial fitting utilities.

use crate::constants::NCOEF;

/// Find the median of a mutable slice using the quick select algorithm.
/// Modifies the input slice order. Returns the median value.
pub fn quick_select(arr: &mut [f32]) -> f32 {
    todo!()
}

/// Fit a 3rd-order polynomial to (aot, atm) data pairs.
/// Returns coefficients [a3, a2, a1, a0] for: a3*x^3 + a2*x^2 + a1*x + a0.
pub fn get_3rd_order_poly_coeff(aot: &[f64], atm: &[f64]) -> [f64; NCOEF] {
    todo!()
}

/// Evaluate a cubic polynomial: coeff[0]*x^3 + coeff[1]*x^2 + coeff[2]*x + coeff[3].
pub fn eval_polynomial(coeff: &[f64; NCOEF], x: f64) -> f64 {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_select_odd() {
        let mut arr = vec![9.0, 1.0, 5.0, 3.0, 7.0];
        assert_eq!(quick_select(&mut arr), 5.0);
    }

    #[test]
    fn test_quick_select_even() {
        let mut arr = vec![4.0, 2.0, 6.0, 8.0];
        // Median of even-length array: middle element at index n/2
        let med = quick_select(&mut arr);
        assert_eq!(med, 6.0);
    }

    #[test]
    fn test_quick_select_single() {
        let mut arr = vec![42.0];
        assert_eq!(quick_select(&mut arr), 42.0);
    }

    #[test]
    fn test_quick_select_two() {
        let mut arr = vec![10.0, 5.0];
        let med = quick_select(&mut arr);
        assert_eq!(med, 10.0); // median index = 1 for n=2: (0+1)/2 = 0... actually the larger
    }

    #[test]
    fn test_poly_coeff_linear() {
        // If atm = 2*aot + 1, cubic and quadratic terms should be ~0
        let aot: Vec<f64> = (0..10).map(|i| i as f64 * 0.5).collect();
        let atm: Vec<f64> = aot.iter().map(|&x| 2.0 * x + 1.0).collect();
        let coeff = get_3rd_order_poly_coeff(&aot, &atm);
        assert!(coeff[0].abs() < 1e-6, "a3 should be ~0: {}", coeff[0]);
        assert!(coeff[1].abs() < 1e-6, "a2 should be ~0: {}", coeff[1]);
        assert!((coeff[2] - 2.0).abs() < 1e-6, "a1 should be ~2: {}", coeff[2]);
        assert!((coeff[3] - 1.0).abs() < 1e-6, "a0 should be ~1: {}", coeff[3]);
    }

    #[test]
    fn test_eval_polynomial() {
        let coeff = [1.0, -2.0, 3.0, 4.0]; // x^3 - 2x^2 + 3x + 4
        assert!((eval_polynomial(&coeff, 0.0) - 4.0).abs() < 1e-10);
        assert!((eval_polynomial(&coeff, 1.0) - 6.0).abs() < 1e-10); // 1 - 2 + 3 + 4
        assert!((eval_polynomial(&coeff, 2.0) - 14.0).abs() < 1e-10); // 8 - 8 + 6 + 4
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test utils`
Expected: FAIL with "not yet implemented"

- [ ] **Step 3: Implement quick_select**

```rust
pub fn quick_select(arr: &mut [f32]) -> f32 {
    let n = arr.len();
    if n == 0 {
        return 0.0;
    }
    let median = n / 2;
    let mut low: usize = 0;
    let mut high: usize = n - 1;

    loop {
        if high <= low {
            return arr[median];
        }
        if high == low + 1 {
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        let middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }
        arr.swap(middle, low + 1);

        let mut ll = low + 1;
        let mut hh = high;
        loop {
            ll += 1;
            while arr[low] > arr[ll] {
                ll += 1;
            }
            hh -= 1;
            while arr[hh] > arr[low] {
                hh -= 1;
            }
            if hh < ll {
                break;
            }
            arr.swap(ll, hh);
        }
        arr.swap(low, hh);

        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}
```

- [ ] **Step 4: Implement polynomial fitting and evaluation**

```rust
pub fn eval_polynomial(coeff: &[f64; NCOEF], x: f64) -> f64 {
    coeff[0] * x * x * x + coeff[1] * x * x + coeff[2] * x + coeff[3]
}

pub fn get_3rd_order_poly_coeff(aot: &[f64], atm: &[f64]) -> [f64; NCOEF] {
    let n = atm.len();

    // Build design matrix columns: x^3, x^2, x, 1
    // Compute normal equations: Z'Z * coeff = Z'y
    let mut z = [[0.0f64; NCOEF]; NCOEF];
    let mut y1 = [0.0f64; NCOEF];

    for i in 0..NCOEF {
        for j in 0..NCOEF {
            let mut sum = 0.0;
            for k in 0..n {
                let xi = aot[k].powi((3 - i) as i32);
                let xj = aot[k].powi((3 - j) as i32);
                sum += xi * xj;
            }
            z[i][j] = sum;
        }
        let mut sum = 0.0;
        for k in 0..n {
            let xi = aot[k].powi((3 - i) as i32);
            sum += xi * atm[k];
        }
        y1[i] = sum;
    }

    // LU decomposition with partial pivoting
    let mut a = z;
    let mut p = [0usize; NCOEF];
    for i in 0..NCOEF {
        p[i] = i;
    }

    for i in 0..NCOEF {
        let mut max_val = 0.0f64;
        let mut max_idx = i;
        for k in i..NCOEF {
            let abs_val = a[k][i].abs();
            if abs_val > max_val {
                max_val = abs_val;
                max_idx = k;
            }
        }
        if max_val < 1e-10 {
            return [0.0; NCOEF]; // Degenerate matrix
        }
        if max_idx != i {
            p.swap(i, max_idx);
            a.swap(i, max_idx);
        }
        for j in (i + 1)..NCOEF {
            a[j][i] /= a[i][i];
            for k in (i + 1)..NCOEF {
                a[j][k] -= a[j][i] * a[i][k];
            }
        }
    }

    // Forward substitution
    let mut x = [0.0f64; NCOEF];
    for i in 0..NCOEF {
        x[i] = y1[p[i]];
        for j in 0..i {
            x[i] -= a[i][j] * x[j];
        }
    }

    // Back substitution
    for i in (0..NCOEF).rev() {
        for j in (i + 1)..NCOEF {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }

    x
}
```

- [ ] **Step 5: Update lib.rs**

```rust
pub mod constants;
pub mod geometry;
pub mod utils;
```

- [ ] **Step 6: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test utils`
Expected: All 6 tests PASS

- [ ] **Step 7: Commit**

```bash
git add lasrc_core/src/utils.rs lasrc_core/src/lib.rs
git commit -m "feat: implement quick select median and polynomial fitting"
```

---

## Task 4: Sensor Trait and Implementations

**Files:**
- Create: `lasrc_core/src/sensor.rs`
- Modify: `lasrc_core/src/lib.rs`

- [ ] **Step 1: Write failing tests**

Add to `lasrc_core/src/sensor.rs`:

```rust
//! Sensor trait and concrete implementations for Landsat 8/9 and Sentinel-2.

use crate::constants::*;

/// Configuration for a spectral band.
#[derive(Debug, Clone)]
pub struct BandConfig {
    pub name: &'static str,
    pub wavelength_um: f64,
    pub native_resolution_m: f64,
}

/// Maps logical band roles to indices in the sensor's band array.
#[derive(Debug, Clone)]
pub struct BandIndices {
    pub coastal: usize,
    pub blue: usize,
    pub green: usize,
    pub red: usize,
    pub nir: usize,
    pub swir1: usize,
    pub swir2: usize,
}

/// Defines a satellite sensor's spectral and spatial configuration.
pub trait Sensor: Send + Sync {
    fn name(&self) -> &str;
    fn reflectance_bands(&self) -> &[BandConfig];
    fn thermal_bands(&self) -> &[BandConfig];
    fn output_resolution_m(&self) -> f64;
    fn aerosol_window(&self) -> usize;
    fn half_aerosol_window(&self) -> usize;
    fn fix_aerosol_window(&self) -> usize;
    fn half_fix_aerosol_window(&self) -> usize;
    fn min_clear_pix(&self) -> usize;
    fn band_indices(&self) -> &BandIndices;
    fn tauray(&self) -> &[f64];
    fn lambda(&self) -> &[f64];
    fn num_refl_bands(&self) -> usize {
        self.reflectance_bands().len()
    }
}

pub struct Landsat8;
pub struct Landsat9;
pub struct Sentinel2A;
pub struct Sentinel2B;
pub struct Sentinel2C;

impl Sensor for Landsat8 {
    fn name(&self) -> &str { todo!() }
    fn reflectance_bands(&self) -> &[BandConfig] { todo!() }
    fn thermal_bands(&self) -> &[BandConfig] { todo!() }
    fn output_resolution_m(&self) -> f64 { todo!() }
    fn aerosol_window(&self) -> usize { todo!() }
    fn half_aerosol_window(&self) -> usize { todo!() }
    fn fix_aerosol_window(&self) -> usize { todo!() }
    fn half_fix_aerosol_window(&self) -> usize { todo!() }
    fn min_clear_pix(&self) -> usize { todo!() }
    fn band_indices(&self) -> &BandIndices { todo!() }
    fn tauray(&self) -> &[f64] { todo!() }
    fn lambda(&self) -> &[f64] { todo!() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_landsat8_band_count() {
        let sensor = Landsat8;
        assert_eq!(sensor.reflectance_bands().len(), 8);
        assert_eq!(sensor.thermal_bands().len(), 2);
    }

    #[test]
    fn test_landsat8_resolution() {
        let sensor = Landsat8;
        assert_eq!(sensor.output_resolution_m(), 30.0);
    }

    #[test]
    fn test_landsat8_aerosol_window() {
        let sensor = Landsat8;
        assert_eq!(sensor.aerosol_window(), LAERO_WINDOW);
        assert_eq!(sensor.half_aerosol_window(), LHALF_AERO_WINDOW);
    }

    #[test]
    fn test_landsat8_band_indices() {
        let sensor = Landsat8;
        let idx = sensor.band_indices();
        assert_eq!(idx.coastal, 0);
        assert_eq!(idx.blue, 1);
        assert_eq!(idx.red, 3);
        assert_eq!(idx.nir, 4);
        assert_eq!(idx.swir2, 6);
    }

    #[test]
    fn test_landsat8_tauray() {
        let sensor = Landsat8;
        assert_eq!(sensor.tauray().len(), 8);
        assert!((sensor.tauray()[0] - 0.23638).abs() < 1e-5);
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test sensor`
Expected: FAIL with "not yet implemented"

- [ ] **Step 3: Implement Landsat8**

Replace the `Landsat8` impl:

```rust
const LANDSAT_REFL_BANDS: [BandConfig; 8] = [
    BandConfig { name: "band1", wavelength_um: 0.443, native_resolution_m: 30.0 },
    BandConfig { name: "band2", wavelength_um: 0.480, native_resolution_m: 30.0 },
    BandConfig { name: "band3", wavelength_um: 0.585, native_resolution_m: 30.0 },
    BandConfig { name: "band4", wavelength_um: 0.655, native_resolution_m: 30.0 },
    BandConfig { name: "band5", wavelength_um: 0.865, native_resolution_m: 30.0 },
    BandConfig { name: "band6", wavelength_um: 1.610, native_resolution_m: 30.0 },
    BandConfig { name: "band7", wavelength_um: 2.200, native_resolution_m: 30.0 },
    BandConfig { name: "band9", wavelength_um: 1.370, native_resolution_m: 30.0 },
];

const LANDSAT_THM_BANDS: [BandConfig; 2] = [
    BandConfig { name: "band10", wavelength_um: 10.9, native_resolution_m: 100.0 },
    BandConfig { name: "band11", wavelength_um: 12.0, native_resolution_m: 100.0 },
];

const LANDSAT_BAND_INDICES: BandIndices = BandIndices {
    coastal: 0,
    blue: 1,
    green: 2,
    red: 3,
    nir: 4,
    swir1: 5,
    swir2: 6,
};

impl Sensor for Landsat8 {
    fn name(&self) -> &str { "LANDSAT_8" }
    fn reflectance_bands(&self) -> &[BandConfig] { &LANDSAT_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &LANDSAT_THM_BANDS }
    fn output_resolution_m(&self) -> f64 { 30.0 }
    fn aerosol_window(&self) -> usize { LAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { LHALF_AERO_WINDOW }
    fn fix_aerosol_window(&self) -> usize { LFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { LHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { LMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &LANDSAT_BAND_INDICES }
    fn tauray(&self) -> &[f64] { &TAURAY_LANDSAT }
    fn lambda(&self) -> &[f64] { &LAMBDA_LANDSAT }
}
```

- [ ] **Step 4: Implement Landsat9 (delegates to same constants as Landsat8)**

```rust
impl Sensor for Landsat9 {
    fn name(&self) -> &str { "LANDSAT_9" }
    fn reflectance_bands(&self) -> &[BandConfig] { &LANDSAT_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &LANDSAT_THM_BANDS }
    fn output_resolution_m(&self) -> f64 { 30.0 }
    fn aerosol_window(&self) -> usize { LAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { LHALF_AERO_WINDOW }
    fn fix_aerosol_window(&self) -> usize { LFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { LHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { LMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &LANDSAT_BAND_INDICES }
    fn tauray(&self) -> &[f64] { &TAURAY_LANDSAT }
    fn lambda(&self) -> &[f64] { &LAMBDA_LANDSAT }
}
```

- [ ] **Step 5: Implement Sentinel2A/B/C**

```rust
const SENTINEL_REFL_BANDS: [BandConfig; 11] = [
    BandConfig { name: "B01", wavelength_um: 0.443, native_resolution_m: 60.0 },
    BandConfig { name: "B02", wavelength_um: 0.490, native_resolution_m: 10.0 },
    BandConfig { name: "B03", wavelength_um: 0.560, native_resolution_m: 10.0 },
    BandConfig { name: "B04", wavelength_um: 0.665, native_resolution_m: 10.0 },
    BandConfig { name: "B05", wavelength_um: 0.705, native_resolution_m: 20.0 },
    BandConfig { name: "B06", wavelength_um: 0.740, native_resolution_m: 20.0 },
    BandConfig { name: "B07", wavelength_um: 0.783, native_resolution_m: 20.0 },
    BandConfig { name: "B08", wavelength_um: 0.842, native_resolution_m: 10.0 },
    BandConfig { name: "B8A", wavelength_um: 0.865, native_resolution_m: 20.0 },
    BandConfig { name: "B11", wavelength_um: 1.610, native_resolution_m: 20.0 },
    BandConfig { name: "B12", wavelength_um: 2.190, native_resolution_m: 20.0 },
];

const SENTINEL_BAND_INDICES: BandIndices = BandIndices {
    coastal: 0,
    blue: 1,
    green: 2,
    red: 3,
    nir: 8, // B8A (narrow NIR at 865nm)
    swir1: 9,
    swir2: 10,
};

impl Sensor for Sentinel2A {
    fn name(&self) -> &str { "SENTINEL_2A" }
    fn reflectance_bands(&self) -> &[BandConfig] { &SENTINEL_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &[] }
    fn output_resolution_m(&self) -> f64 { 10.0 }
    fn aerosol_window(&self) -> usize { SAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { SAERO_WINDOW / 2 }
    fn fix_aerosol_window(&self) -> usize { SFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { SHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { SMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &SENTINEL_BAND_INDICES }
    fn tauray(&self) -> &[f64] { todo!("Sentinel Rayleigh values loaded from LUT") }
    fn lambda(&self) -> &[f64] { &LAMBDA_SENTINEL }
}

impl Sensor for Sentinel2B {
    fn name(&self) -> &str { "SENTINEL_2B" }
    fn reflectance_bands(&self) -> &[BandConfig] { &SENTINEL_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &[] }
    fn output_resolution_m(&self) -> f64 { 10.0 }
    fn aerosol_window(&self) -> usize { SAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { SAERO_WINDOW / 2 }
    fn fix_aerosol_window(&self) -> usize { SFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { SHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { SMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &SENTINEL_BAND_INDICES }
    fn tauray(&self) -> &[f64] { todo!("Sentinel Rayleigh values loaded from LUT") }
    fn lambda(&self) -> &[f64] { &LAMBDA_SENTINEL }
}

impl Sensor for Sentinel2C {
    fn name(&self) -> &str { "SENTINEL_2C" }
    fn reflectance_bands(&self) -> &[BandConfig] { &SENTINEL_REFL_BANDS }
    fn thermal_bands(&self) -> &[BandConfig] { &[] }
    fn output_resolution_m(&self) -> f64 { 10.0 }
    fn aerosol_window(&self) -> usize { SAERO_WINDOW }
    fn half_aerosol_window(&self) -> usize { SAERO_WINDOW / 2 }
    fn fix_aerosol_window(&self) -> usize { SFIX_AERO_WINDOW }
    fn half_fix_aerosol_window(&self) -> usize { SHALF_FIX_AERO_WINDOW }
    fn min_clear_pix(&self) -> usize { SMIN_CLEAR_PIX }
    fn band_indices(&self) -> &BandIndices { &SENTINEL_BAND_INDICES }
    fn tauray(&self) -> &[f64] { todo!("Sentinel Rayleigh values loaded from LUT") }
    fn lambda(&self) -> &[f64] { &LAMBDA_SENTINEL }
}
```

- [ ] **Step 6: Update lib.rs**

```rust
pub mod constants;
pub mod geometry;
pub mod sensor;
pub mod utils;
```

- [ ] **Step 7: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test sensor`
Expected: All 5 tests PASS

- [ ] **Step 8: Commit**

```bash
git add lasrc_core/src/sensor.rs lasrc_core/src/lib.rs
git commit -m "feat: implement Sensor trait with Landsat and Sentinel configs"
```

---

## Task 5: Gas Transmission and Rayleigh Scattering

**Files:**
- Create: `lasrc_core/src/gas_transmission.rs`
- Create: `lasrc_core/src/rayleigh.rs`
- Modify: `lasrc_core/src/lib.rs`

Porting `comptg()` and `local_chand()` from `lut_subr.c`.

- [ ] **Step 1: Write failing tests for gas transmission**

Add to `lasrc_core/src/gas_transmission.rs`:

```rust
//! Gas transmission computations for ozone, water vapor, and other gases.

/// Gas transmission coefficients for a single band, loaded from LUT.
pub struct GasCoefficients {
    pub ogtransa1: f64,
    pub ogtransb0: f64,
    pub ogtransb1: f64,
    pub wvtransa: f64,
    pub wvtransb: f64,
    pub oztransa: f64,
}

/// Computed gas transmissions for a pixel.
pub struct GasTransmission {
    pub tgoz: f64,   // ozone transmission
    pub tgwv: f64,   // water vapor transmission
    pub tgwvhalf: f64, // water vapor half-path transmission
    pub tgog: f64,   // other gases transmission
    pub tgo: f64,    // total other-gas transmission (tgog * tgoz)
}

/// Compute gas transmissions for a given band and atmospheric state.
pub fn compute_gas_transmission(
    coeff: &GasCoefficients,
    xmus: f64,     // cosine of solar zenith angle
    xmuv: f64,     // cosine of view zenith angle
    uoz: f64,      // ozone amount (cm-atm)
    uwv: f64,      // water vapor (g/cm^2)
    atm_pres: f64, // normalized atmospheric pressure (pres/1013)
) -> GasTransmission {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_coeff() -> GasCoefficients {
        // Approximate values for Landsat Band 1 (coastal aerosol)
        GasCoefficients {
            ogtransa1: -0.00014,
            ogtransb0: 0.00527,
            ogtransb1: 0.00015,
            wvtransa: 2.29849e-27,
            wvtransb: 14.77305,
            oztransa: -0.00160,
        }
    }

    #[test]
    fn test_gas_transmission_sea_level() {
        let coeff = sample_coeff();
        let gt = compute_gas_transmission(&coeff, 0.5, 0.95, 0.3, 2.5, 1.0);
        assert!(gt.tgoz > 0.0 && gt.tgoz <= 1.0, "tgoz={}", gt.tgoz);
        assert!(gt.tgwv > 0.0 && gt.tgwv <= 1.0, "tgwv={}", gt.tgwv);
        assert!(gt.tgog > 0.0 && gt.tgog <= 1.0, "tgog={}", gt.tgog);
        assert!(gt.tgo > 0.0 && gt.tgo <= 1.0, "tgo={}", gt.tgo);
    }

    #[test]
    fn test_gas_transmission_zero_wv() {
        let coeff = sample_coeff();
        let gt = compute_gas_transmission(&coeff, 0.5, 0.95, 0.3, 0.0, 1.0);
        assert!((gt.tgwv - 1.0).abs() < 1e-6, "No water vapor: tgwv should be ~1.0");
    }
}
```

- [ ] **Step 2: Write failing tests for Rayleigh**

Add to `lasrc_core/src/rayleigh.rs`:

```rust
//! Rayleigh (molecular) scattering reflectance computation.

/// Compute Rayleigh scattering reflectance.
///
/// Implements the `local_chand` function from C code. Uses depolarization
/// factor and 10-coefficient Legendre polynomial expansion.
pub fn rayleigh_reflectance(
    xphi: f64,   // relative azimuth angle (degrees)
    xmuv: f64,   // cosine of view zenith angle
    xmus: f64,   // cosine of solar zenith angle
    xtau: f64,   // Rayleigh optical depth
) -> f64 {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rayleigh_positive() {
        let rray = rayleigh_reflectance(90.0, 0.95, 0.5, 0.17);
        assert!(rray > 0.0, "Rayleigh reflectance should be positive: {rray}");
    }

    #[test]
    fn test_rayleigh_increases_with_tau() {
        let r1 = rayleigh_reflectance(90.0, 0.95, 0.5, 0.05);
        let r2 = rayleigh_reflectance(90.0, 0.95, 0.5, 0.25);
        assert!(r2 > r1, "More optical depth = more scattering: r1={r1}, r2={r2}");
    }

    #[test]
    fn test_rayleigh_range() {
        // Rayleigh reflectance should be small (< 0.5 for reasonable conditions)
        let rray = rayleigh_reflectance(45.0, 0.7, 0.6, 0.24);
        assert!(rray > 0.0 && rray < 0.5, "rray={rray}");
    }
}
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test gas_transmission && cargo test rayleigh`
Expected: FAIL with "not yet implemented"

- [ ] **Step 4: Implement compute_gas_transmission**

```rust
pub fn compute_gas_transmission(
    coeff: &GasCoefficients,
    xmus: f64,
    xmuv: f64,
    uoz: f64,
    uwv: f64,
    atm_pres: f64,
) -> GasTransmission {
    let m = 1.0 / xmus + 1.0 / xmuv;

    // Ozone transmission
    let tgoz = (coeff.oztransa * m * uoz).exp();

    // Water vapor transmission
    let x = m * uwv;
    let tgwv = if x > 1.0e-06 {
        (-coeff.wvtransa * x.powf(coeff.wvtransb)).exp()
    } else {
        1.0
    };

    // Water vapor half-path (solar path only)
    let xhalf = uwv / xmus;
    let tgwvhalf = if xhalf > 1.0e-06 {
        (-coeff.wvtransa * xhalf.powf(coeff.wvtransb)).exp()
    } else {
        1.0
    };

    // Other gases transmission
    let tgog = (-(coeff.ogtransa1 * atm_pres)
        * m.powf((-coeff.ogtransb0 + coeff.ogtransb1 * atm_pres).exp().neg()))
    .exp();

    let tgo = tgog * tgoz;

    GasTransmission {
        tgoz,
        tgwv,
        tgwvhalf,
        tgog,
        tgo,
    }
}

// Need this for the neg() call
use std::ops::Neg;
```

Wait — the C code for `tgog` is:
```c
tgog = exp(-(ogtransa1 * atm_pres) * pow(m, exp(-(ogtransb0 + ogtransb1 * atm_pres))))
```

Let me correct:

```rust
pub fn compute_gas_transmission(
    coeff: &GasCoefficients,
    xmus: f64,
    xmuv: f64,
    uoz: f64,
    uwv: f64,
    atm_pres: f64,
) -> GasTransmission {
    let m = 1.0 / xmus + 1.0 / xmuv;

    // Ozone transmission
    let tgoz = (coeff.oztransa * m * uoz).exp();

    // Water vapor transmission
    let x = m * uwv;
    let tgwv = if x > 1.0e-06 {
        (-coeff.wvtransa * x.powf(coeff.wvtransb)).exp()
    } else {
        1.0
    };

    // Water vapor half-path (solar path only)
    let xhalf = uwv / xmus;
    let tgwvhalf = if xhalf > 1.0e-06 {
        (-coeff.wvtransa * xhalf.powf(coeff.wvtransb)).exp()
    } else {
        1.0
    };

    // Other gases transmission
    let exponent = (-(coeff.ogtransb0 + coeff.ogtransb1 * atm_pres)).exp();
    let tgog = (-(coeff.ogtransa1 * atm_pres) * m.powf(exponent)).exp();

    let tgo = tgog * tgoz;

    GasTransmission {
        tgoz,
        tgwv,
        tgwvhalf,
        tgog,
        tgo,
    }
}
```

- [ ] **Step 5: Implement rayleigh_reflectance**

```rust
pub fn rayleigh_reflectance(
    xphi: f64,
    xmuv: f64,
    xmus: f64,
    xtau: f64,
) -> f64 {
    use std::f64::consts::PI;

    let xfd = 0.958725777; // derived from depolarization factor 0.0279
    let xmus2 = xmus * xmus;
    let xmuv2 = xmuv * xmuv;

    // Phase function components
    let xph1 = 1.0 + (3.0 * xmus2 - 1.0) * (3.0 * xmuv2 - 1.0) * xfd * 0.125;
    let xph3 = (1.0 - xmus2) * (1.0 - xmuv2);
    let xph2 = -xmus * xmuv * xph3.sqrt() * xfd * 0.75;
    let xph3 = xph3 * xfd * 0.1875;

    // Azimuth terms
    let phi_rad = xphi * PI / 180.0;
    let xcosf1 = 1.0;
    let xcosf2 = (2.0 * phi_rad).cos();
    let xcosf3 = (4.0 * phi_rad).cos(); // cos(2*2*phi)

    // Legendre coefficients
    let as0: [f64; 10] = [
        0.33243832, -6.777104e-02, 0.16285370, 1.577425e-03,
        -0.30924818, -1.240906e-02, -0.10324388, 3.241678e-02,
        0.11493334, -3.503695e-02,
    ];
    let as1: [f64; 2] = [0.19666292, -5.439061e-02];
    let as2: [f64; 2] = [0.14545937, -2.910845e-02];

    // Build polynomial terms
    let log_tau = xtau.ln();
    let pl = [
        1.0,
        log_tau,
        xmus + xmuv,
        log_tau * (xmus + xmuv),
        xmus * xmuv,
        log_tau * xmus * xmuv,
        xmus2 + xmuv2,
        log_tau * (xmus2 + xmuv2),
        xmus2 * xmuv2,
        log_tau * xmus2 * xmuv2,
    ];

    let mut xitot1 = 0.0;
    for i in 0..10 {
        xitot1 += as0[i] * pl[i];
    }
    let mut xitot2 = 0.0;
    for i in 0..2 {
        xitot2 += as1[i] * pl[i];
    }
    let mut xitot3 = 0.0;
    for i in 0..2 {
        xitot3 += as2[i] * pl[i];
    }

    let xrray = xph1 * xitot1 + 2.0 * (xph2 * xcosf2 * xitot2 + xph3 * xcosf3 * xitot3);
    xrray
}
```

- [ ] **Step 6: Update lib.rs**

```rust
pub mod constants;
pub mod gas_transmission;
pub mod geometry;
pub mod rayleigh;
pub mod sensor;
pub mod utils;
```

- [ ] **Step 7: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test gas_transmission && cargo test rayleigh`
Expected: All 5 tests PASS

- [ ] **Step 8: Commit**

```bash
git add lasrc_core/src/gas_transmission.rs lasrc_core/src/rayleigh.rs lasrc_core/src/lib.rs
git commit -m "feat: implement gas transmission and Rayleigh scattering"
```

---

## Task 6: LUT Interpolation

**Files:**
- Create: `lasrc_core/src/lut.rs`
- Modify: `lasrc_core/src/lib.rs`

Porting `compsalb`, `comptrans`, `comproatm`, and the scattering angle interpolation from `lut_subr.c`.

- [ ] **Step 1: Write LUT struct and failing tests**

Add to `lasrc_core/src/lut.rs`:

```rust
//! Lookup table storage and multi-dimensional interpolation.
//!
//! The LUT data is pre-computed by the 6S radiative transfer model and stored
//! in HDF files. Python loads these into flat arrays and passes them to Rust.

use crate::constants::*;

/// Pre-computed atmospheric lookup tables.
pub struct LookupTables {
    /// Intrinsic atmospheric reflectance.
    /// Dimensions: [NSR_BANDS][NPRES_VALS][NAOT_VALS][NSOLAR_VALS]
    pub rolutt: Vec<f64>,
    /// Atmospheric transmission.
    /// Dimensions: [NSR_BANDS][NPRES_VALS][NAOT_VALS][NSUNANGLE_VALS]
    pub transt: Vec<f64>,
    /// Spherical albedo.
    /// Dimensions: [NSR_BANDS][NPRES_VALS][NAOT_VALS]
    pub sphalbt: Vec<f64>,
    /// Aerosol extinction normalization.
    /// Dimensions: [NSR_BANDS][NPRES_VALS][NAOT_VALS]
    pub normext: Vec<f64>,

    // Angle tables: [NVIEW_ZEN_VALS][NSOLAR_ZEN_VALS]
    pub tsmax: Vec<f64>,
    pub tsmin: Vec<f64>,
    pub nbfic: Vec<f64>,
    pub nbfi: Vec<i32>,
    pub ttv: Vec<f64>,
    pub tts: [f64; NSOLAR_ZEN_VALS],

    /// Number of SR bands in this LUT (8 for Landsat, 11 for Sentinel)
    pub nsr_bands: usize,
}

/// Indices into the pressure and AOT tables for a given pixel state.
#[derive(Debug, Clone, Copy)]
pub struct LutIndices {
    pub ip1: usize,
    pub ip2: usize,
    pub iaot1: usize,
    pub iaot2: usize,
}

impl LookupTables {
    /// Find pressure and AOT interpolation indices.
    pub fn find_indices(&self, pressure: f64, raot550nm: f64) -> LutIndices {
        todo!()
    }

    /// Interpolate spherical albedo for given band, pressure, and AOT.
    pub fn interp_spherical_albedo(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
    ) -> f64 {
        todo!()
    }

    /// Interpolate atmospheric transmission for a zenith angle.
    pub fn interp_transmission(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64, // zenith angle in degrees
    ) -> f64 {
        todo!()
    }

    /// Interpolate atmospheric intrinsic reflectance.
    pub fn interp_atmospheric_reflectance(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64,
        xtv: f64,
        xmus: f64,
        xmuv: f64,
        cosxfi: f64,
    ) -> f64 {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_lut() -> LookupTables {
        // Create a small LUT with known values for testing.
        // Use 1 band to keep it manageable.
        let nsr = 1;

        // sphalbt: linearly increasing with AOT, constant across pressure
        let mut sphalbt = vec![0.0; nsr * NPRES_VALS * NAOT_VALS];
        for ip in 0..NPRES_VALS {
            for ia in 0..NAOT_VALS {
                sphalbt[ip * NAOT_VALS + ia] = AOT550NM[ia] * 0.1;
            }
        }

        // transt: simple decreasing with AOT
        let mut transt = vec![0.0; nsr * NPRES_VALS * NAOT_X_NSUNANGLE];
        for ip in 0..NPRES_VALS {
            for ia in 0..NAOT_VALS {
                for is_ in 0..NSUNANGLE_VALS {
                    let idx = ip * NAOT_X_NSUNANGLE + ia * NSUNANGLE_VALS + is_;
                    transt[idx] = 1.0 - AOT550NM[ia] * 0.1;
                }
            }
        }

        let normext = vec![1.0; nsr * NPRES_VALS * NAOT_VALS];
        let rolutt = vec![0.01; nsr * NPRES_VALS * NAOT_X_NSOLAR];
        let tsmax = vec![180.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let tsmin = vec![0.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let nbfic = vec![45.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let nbfi = vec![45; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let ttv = vec![3.0; NVIEW_ZEN_VALS * NSOLAR_ZEN_VALS];
        let mut tts = [0.0f64; NSOLAR_ZEN_VALS];
        for i in 0..NSOLAR_ZEN_VALS {
            tts[i] = i as f64 * 4.0;
        }

        LookupTables {
            rolutt,
            transt,
            sphalbt,
            normext,
            tsmax,
            tsmin,
            nbfic,
            nbfi,
            ttv,
            tts,
            nsr_bands: nsr,
        }
    }

    #[test]
    fn test_find_indices_sea_level() {
        let lut = make_test_lut();
        let idx = lut.find_indices(1013.0, 0.15);
        // 1013.0 is tpres[1], so ip1=0, ip2=1
        assert_eq!(idx.ip1, 0);
        assert_eq!(idx.ip2, 1);
        // 0.15 is aot550nm[3], so iaot1=2, iaot2=3
        assert_eq!(idx.iaot1, 2);
        assert_eq!(idx.iaot2, 3);
    }

    #[test]
    fn test_interp_spherical_albedo_monotonic() {
        let lut = make_test_lut();
        let idx1 = lut.find_indices(1013.0, 0.1);
        let sa1 = lut.interp_spherical_albedo(&idx1, 0, 1013.0, 0.1);
        let idx2 = lut.find_indices(1013.0, 0.5);
        let sa2 = lut.interp_spherical_albedo(&idx2, 0, 1013.0, 0.5);
        assert!(sa2 > sa1, "Higher AOT = higher albedo: {sa1} vs {sa2}");
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test lut`
Expected: FAIL with "not yet implemented"

- [ ] **Step 3: Implement find_indices**

```rust
impl LookupTables {
    pub fn find_indices(&self, pressure: f64, raot550nm: f64) -> LutIndices {
        // Find pressure indices
        let mut ip1 = 0usize;
        for i in 0..(NPRES_VALS - 1) {
            if pressure < TPRES[i] {
                ip1 = i;
            }
        }
        let ip2 = (ip1 + 1).min(NPRES_VALS - 1);

        // Find AOT indices
        let mut iaot1 = 0usize;
        for i in 0..(NAOT_VALS - 1) {
            if raot550nm > AOT550NM[i] {
                iaot1 = i;
            }
        }
        let iaot2 = (iaot1 + 1).min(NAOT_VALS - 1);

        LutIndices { ip1, ip2, iaot1, iaot2 }
    }
}
```

- [ ] **Step 4: Implement interp_spherical_albedo (compsalb)**

```rust
    pub fn interp_spherical_albedo(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;
        let iband_offset = iband * NPRES_VALS * NAOT_VALS;

        let deltaaot = if iaot1 != iaot2 {
            (raot550nm - AOT550NM[iaot1]) / (AOT550NM[iaot2] - AOT550NM[iaot1])
        } else {
            0.0
        };

        // Interpolate along AOT at pressure level 1
        let v11 = self.sphalbt[iband_offset + ip1 * NAOT_VALS + iaot1];
        let v12 = self.sphalbt[iband_offset + ip1 * NAOT_VALS + iaot2];
        let satm1 = v11 + (v12 - v11) * deltaaot;

        // Interpolate along AOT at pressure level 2
        let v21 = self.sphalbt[iband_offset + ip2 * NAOT_VALS + iaot1];
        let v22 = self.sphalbt[iband_offset + ip2 * NAOT_VALS + iaot2];
        let satm2 = v21 + (v22 - v21) * deltaaot;

        // Interpolate along pressure
        let dpres = if ip1 != ip2 {
            (pressure - TPRES[ip1]) / (TPRES[ip2] - TPRES[ip1])
        } else {
            0.0
        };

        satm1 + (satm2 - satm1) * dpres
    }
```

- [ ] **Step 5: Implement interp_transmission (comptrans)**

```rust
    pub fn interp_transmission(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64,
    ) -> f64 {
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;
        let iband_offset = iband * NPRES_VALS * NAOT_X_NSUNANGLE;

        // Find sun angle index
        let its = if xts <= self.tts[0] {
            0usize
        } else {
            let idx = ((xts - self.tts[0]) / XTS_STEP) as usize;
            idx.min(NSUNANGLE_VALS - 2)
        };

        let xmts = (xts - self.tts[its]) * 0.25; // fractional position (step=4.0)

        let deltaaot = if iaot1 != iaot2 {
            (raot550nm - AOT550NM[iaot1]) / (AOT550NM[iaot2] - AOT550NM[iaot1])
        } else {
            0.0
        };

        // Helper to read transt with angle interpolation
        let interp_angle = |ip: usize, iaot: usize| -> f64 {
            let base = iband_offset + ip * NAOT_X_NSUNANGLE + iaot * NSUNANGLE_VALS;
            let t0 = self.transt[base + its];
            let t1 = self.transt[base + its + 1];
            t0 + (t1 - t0) * xmts
        };

        // 4-point interpolation: (ip1,iaot1), (ip1,iaot2), (ip2,iaot1), (ip2,iaot2)
        let t_ip1_ia1 = interp_angle(ip1, iaot1);
        let t_ip1_ia2 = interp_angle(ip1, iaot2);
        let trans1 = t_ip1_ia1 + (t_ip1_ia2 - t_ip1_ia1) * deltaaot;

        let t_ip2_ia1 = interp_angle(ip2, iaot1);
        let t_ip2_ia2 = interp_angle(ip2, iaot2);
        let trans2 = t_ip2_ia1 + (t_ip2_ia2 - t_ip2_ia1) * deltaaot;

        let dpres = if ip1 != ip2 {
            (pressure - TPRES[ip1]) / (TPRES[ip2] - TPRES[ip1])
        } else {
            0.0
        };

        trans1 + (trans2 - trans1) * dpres
    }
```

- [ ] **Step 6: Implement interp_atmospheric_reflectance (comproatm)**

This is the most complex interpolation - uses scattering angle interpolation across 4 pressure/AOT combinations.

```rust
    pub fn interp_atmospheric_reflectance(
        &self,
        indices: &LutIndices,
        iband: usize,
        pressure: f64,
        raot550nm: f64,
        xts: f64,
        xtv: f64,
        xmus: f64,
        xmuv: f64,
        cosxfi: f64,
    ) -> f64 {
        use crate::geometry::scattering_angle;
        let LutIndices { ip1, ip2, iaot1, iaot2 } = *indices;

        let scaa = scattering_angle(xmus, xmuv, cosxfi);

        // Find view/sun angle indices
        let itv = if xtv <= XTV_MIN {
            0usize
        } else {
            ((xtv - XTV_MIN) / XTV_STEP + 1.0) as usize
        };
        let itv = itv.min(NVIEW_ZEN_VALS - 2);

        let its = if xts <= XTS_MIN {
            0usize
        } else {
            ((xts - XTS_MIN) / XTS_STEP) as usize
        };
        let its = its.min(NSOLAR_ZEN_VALS - 2);

        // Interpolation parameters for sun and view angles
        let t = if self.tts[its + 1] != self.tts[its] {
            (self.tts[its + 1] - xts) / (self.tts[its + 1] - self.tts[its])
        } else {
            0.0
        };

        let itv_its = itv * NSOLAR_ZEN_VALS + its;
        let itv1_its = (itv + 1) * NSOLAR_ZEN_VALS + its;
        let u = if self.ttv[itv1_its] != self.ttv[itv_its] {
            (self.ttv[itv1_its] - xtv) / (self.ttv[itv1_its] - self.ttv[itv_its])
        } else {
            0.0
        };

        // AOT interpolation in log space
        let log_raot = raot550nm.max(0.0001).ln();
        let deltaaot = if iaot1 != iaot2 {
            (log_raot - LOG_AOT550NM[iaot1]) / (LOG_AOT550NM[iaot2] - LOG_AOT550NM[iaot1])
        } else {
            0.0
        };

        // Interpolate reflectance at 4 (pressure, AOT) combinations
        let interp_scat = |ip: usize, iaot: usize| -> f64 {
            self.interp_refl_at_scat_angle(iband, ip, iaot, scaa, its, itv, t, u)
        };

        let ro_ip1_ia1 = interp_scat(ip1, iaot1);
        let ro_ip1_ia2 = interp_scat(ip1, iaot2);
        let rop1 = ro_ip1_ia1 + (ro_ip1_ia2 - ro_ip1_ia1) * deltaaot;

        let ro_ip2_ia1 = interp_scat(ip2, iaot1);
        let ro_ip2_ia2 = interp_scat(ip2, iaot2);
        let rop2 = ro_ip2_ia1 + (ro_ip2_ia2 - ro_ip2_ia1) * deltaaot;

        let dpres = if ip1 != ip2 {
            (pressure - TPRES[ip1]) / (TPRES[ip2] - TPRES[ip1])
        } else {
            0.0
        };

        rop1 + (rop2 - rop1) * dpres
    }

    /// Interpolate reflectance at a specific scattering angle.
    /// Implements interp_refl_using_scat_angle from C code.
    fn interp_refl_at_scat_angle(
        &self,
        iband: usize,
        ip: usize,
        iaot: usize,
        scaa: f64,
        its: usize,
        itv: usize,
        t: f64,
        u: f64,
    ) -> f64 {
        let iband_offset = iband * NPRES_VALS * NAOT_X_NSOLAR;
        let ip_offset = ip * NAOT_X_NSOLAR;
        let iaot_offset = iaot * NSOLAR_VALS;

        // Interpolate at 4 corners of (view_zen, solar_zen) grid
        let mut ro = [0.0f64; 4];
        let corners = [
            (itv, its),
            (itv, its + 1),
            (itv + 1, its),
            (itv + 1, its + 1),
        ];

        for (c, &(iv, is)) in corners.iter().enumerate() {
            let angle_idx = iv * NSOLAR_ZEN_VALS + is;
            let xtsmax_val = self.tsmax[angle_idx];
            let xtsmin_val = self.tsmin[angle_idx];
            let nbfi_val = self.nbfi[angle_idx] as usize;

            // Find scattering angle index
            let mut isca = ((xtsmax_val - scaa) * 0.25 + 1.0) as usize;
            if isca < 1 {
                isca = 1;
            }

            let (sca1, sca2);
            if isca + 1 < nbfi_val {
                sca1 = xtsmax_val - (isca as f64 - 1.0) * 4.0;
                sca2 = sca1 - 4.0;
            } else {
                isca = if nbfi_val > 1 { nbfi_val - 1 } else { 1 };
                sca1 = xtsmax_val - (isca as f64 - 1.0) * 4.0;
                sca2 = xtsmin_val;
            }

            let base = iband_offset + ip_offset + iaot_offset;
            let nbfic_val = self.nbfic[angle_idx] as usize;
            let rolutt_idx = base + nbfic_val + isca - 1;

            if rolutt_idx + 1 < self.rolutt.len() {
                let roinf = self.rolutt[rolutt_idx];
                let rosup = self.rolutt[rolutt_idx + 1];
                ro[c] = if (sca2 - sca1).abs() > 1e-10 {
                    roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1)
                } else {
                    roinf
                };
            }
        }

        // Bilinear interpolation of the 4 corner values
        ro[3] + u * (ro[1] - ro[3]) + t * (ro[2] - ro[3]) + u * t * (ro[0] - ro[1] - ro[2] + ro[3])
    }
```

- [ ] **Step 7: Update lib.rs**

```rust
pub mod constants;
pub mod gas_transmission;
pub mod geometry;
pub mod lut;
pub mod rayleigh;
pub mod sensor;
pub mod utils;
```

- [ ] **Step 8: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test lut`
Expected: All tests PASS

- [ ] **Step 9: Commit**

```bash
git add lasrc_core/src/lut.rs lasrc_core/src/lib.rs
git commit -m "feat: implement LUT interpolation for spherical albedo, transmission, and atmospheric reflectance"
```

---

## Task 7: Atmospheric Correction Function

**Files:**
- Create: `lasrc_core/src/atmospheric.rs`
- Modify: `lasrc_core/src/lib.rs`

Porting `atmcorlamb2` and `atmcorlamb2_new` from `lut_subr.c`.

- [ ] **Step 1: Write struct and failing tests**

Add to `lasrc_core/src/atmospheric.rs`:

```rust
//! Atmospheric correction: convert TOA reflectance to surface reflectance.

use crate::constants::*;
use crate::gas_transmission::{GasCoefficients, GasTransmission, compute_gas_transmission};
use crate::lut::{LookupTables, LutIndices};
use crate::rayleigh::rayleigh_reflectance;

/// Result of atmospheric correction for a single pixel and band.
#[derive(Debug, Clone)]
pub struct AtmCorrResult {
    pub roslamb: f64,  // surface reflectance
    pub tgo: f64,      // other gas transmittance
    pub roatm: f64,    // atmospheric intrinsic reflectance
    pub ttatmg: f64,   // total atmospheric transmittance
    pub satm: f64,     // spherical albedo
    pub xrorayp: f64,  // Rayleigh reflectance
}

/// Pre-computed polynomial coefficients for the "new" fast atmospheric correction.
#[derive(Debug, Clone)]
pub struct AtmCorrCoefficients {
    pub roatm_upper: f64,
    pub roatm_coef: [f64; NCOEF],
    pub ttatmg_coef: [f64; NCOEF],
    pub satm_coef: [f64; NCOEF],
}

/// Full atmospheric correction using LUT interpolation (atmcorlamb2).
/// Used for original aerosol algorithm and per-pixel angle correction.
pub fn atmcorlamb2(
    lut: &LookupTables,
    gas_coeff: &GasCoefficients,
    tauray_band: f64,
    iband: usize,
    xts: f64,       // solar zenith angle (degrees)
    xtv: f64,       // view zenith angle (degrees)
    xmus: f64,      // cosine solar zenith
    xmuv: f64,      // cosine view zenith
    xfi: f64,       // relative azimuth (degrees)
    cosxfi: f64,    // cosine of relative azimuth
    raot550nm: f64, // AOT at 550nm
    pressure: f64,  // surface pressure (mbar)
    uoz: f64,       // ozone (cm-atm)
    uwv: f64,       // water vapor (g/cm^2)
    rotoa: f64,     // TOA reflectance
    eps: f64,       // Angstrom coefficient
) -> AtmCorrResult {
    todo!()
}

/// Fast atmospheric correction using pre-computed polynomial coefficients (atmcorlamb2_new).
/// Used for semi-empirical aerosol algorithm with scene-center geometry.
pub fn atmcorlamb2_new(
    coeff: &AtmCorrCoefficients,
    tgo: f64,
    iband: usize,
    raot550nm: f64,
    normext_ib_0_3: f64, // normext at band 0, pressure index 3
    rotoa: f64,
    lambda: &[f64],
    eps: f64,
) -> f64 {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atmcorlamb2_new_low_aot() {
        let coeff = AtmCorrCoefficients {
            roatm_upper: 0.4,
            roatm_coef: [0.0, 0.0, 0.1, 0.01],  // ~linear: 0.1*aot + 0.01
            ttatmg_coef: [0.0, 0.0, -0.1, 0.9],  // ~0.9 - 0.1*aot
            satm_coef: [0.0, 0.0, 0.05, 0.01],   // ~0.05*aot + 0.01
        };

        let roslamb = atmcorlamb2_new(
            &coeff,
            0.95,          // tgo
            0,             // band
            0.05,          // low AOT
            1.0,           // normext
            0.15,          // rotoa
            &LAMBDA_LANDSAT,
            1.5,           // eps
        );
        // With low AOT, surface reflectance should be less than TOA
        assert!(roslamb < 0.15, "roslamb={roslamb}");
        assert!(roslamb > -0.2, "roslamb should be > -0.2: {roslamb}");
    }

    #[test]
    fn test_atmcorlamb2_new_higher_aot_lower_sr() {
        let coeff = AtmCorrCoefficients {
            roatm_upper: 0.5,
            roatm_coef: [0.0, 0.0, 0.2, 0.01],
            ttatmg_coef: [0.0, 0.0, -0.15, 0.9],
            satm_coef: [0.0, 0.0, 0.08, 0.01],
        };

        let sr1 = atmcorlamb2_new(&coeff, 0.95, 0, 0.05, 1.0, 0.15, &LAMBDA_LANDSAT, 1.5);
        let sr2 = atmcorlamb2_new(&coeff, 0.95, 0, 0.50, 1.0, 0.15, &LAMBDA_LANDSAT, 1.5);
        // Higher AOT means more atmospheric path radiance removed, so SR should differ
        assert!((sr1 - sr2).abs() > 0.001, "Different AOT should give different SR");
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test atmospheric`
Expected: FAIL with "not yet implemented"

- [ ] **Step 3: Implement atmcorlamb2_new**

```rust
pub fn atmcorlamb2_new(
    coeff: &AtmCorrCoefficients,
    tgo: f64,
    iband: usize,
    raot550nm: f64,
    normext_ib_0_3: f64,
    rotoa: f64,
    lambda: &[f64],
    eps: f64,
) -> f64 {
    let lambda_sf = 1.0 / 0.55;
    let mut mraot550nm = (raot550nm / normext_ib_0_3)
        * (lambda[iband] * lambda_sf).powf(-eps);

    if mraot550nm >= coeff.roatm_upper {
        mraot550nm = coeff.roatm_upper;
    }

    let mraot_sq = mraot550nm * mraot550nm;
    let mraot_cu = mraot_sq * mraot550nm;

    let roatm = coeff.roatm_coef[0] * mraot_cu
        + coeff.roatm_coef[1] * mraot_sq
        + coeff.roatm_coef[2] * mraot550nm
        + coeff.roatm_coef[3];

    let ttatmg = coeff.ttatmg_coef[0] * mraot_cu
        + coeff.ttatmg_coef[1] * mraot_sq
        + coeff.ttatmg_coef[2] * mraot550nm
        + coeff.ttatmg_coef[3];

    let satm = coeff.satm_coef[0] * mraot_cu
        + coeff.satm_coef[1] * mraot_sq
        + coeff.satm_coef[2] * mraot550nm
        + coeff.satm_coef[3];

    // Solve for surface reflectance
    let xroslamb = rotoa / tgo - roatm;
    xroslamb / (ttatmg + satm * xroslamb)
}
```

- [ ] **Step 4: Implement atmcorlamb2**

```rust
pub fn atmcorlamb2(
    lut: &LookupTables,
    gas_coeff: &GasCoefficients,
    tauray_band: f64,
    iband: usize,
    xts: f64,
    xtv: f64,
    xmus: f64,
    xmuv: f64,
    xfi: f64,
    cosxfi: f64,
    raot550nm: f64,
    pressure: f64,
    uoz: f64,
    uwv: f64,
    rotoa: f64,
    eps: f64,
) -> AtmCorrResult {
    let atm_pres = pressure * ONE_DIV_ATMOS_PRES_0;
    let xtaur = tauray_band * atm_pres;

    // Compute Rayleigh scattering
    let xrorayp = rayleigh_reflectance(xfi, xmuv, xmus, xtaur);

    // Find LUT indices
    let indices = lut.find_indices(pressure, raot550nm);

    // Interpolate atmospheric parameters from LUTs
    let satm = lut.interp_spherical_albedo(&indices, iband, pressure, raot550nm);
    let roatm = lut.interp_atmospheric_reflectance(
        &indices, iband, pressure, raot550nm, xts, xtv, xmus, xmuv, cosxfi,
    );

    // Transmission: downward (solar) and upward (view)
    let xtts = lut.interp_transmission(&indices, iband, pressure, raot550nm, xts);
    let xttv = lut.interp_transmission(&indices, iband, pressure, raot550nm, xtv);
    let ttatm = xtts * xttv;

    // Gas transmissions
    let gt = compute_gas_transmission(gas_coeff, xmus, xmuv, uoz, uwv, atm_pres);

    // Combine atmospheric and Rayleigh reflectance
    let roatm_corrected = (roatm - xrorayp) * gt.tgwvhalf + xrorayp;
    let ttatmg = ttatm * gt.tgwv;

    // Solve for surface reflectance
    let xroslamb = rotoa / gt.tgo - roatm_corrected;
    let roslamb = xroslamb / (ttatmg + satm * xroslamb);

    AtmCorrResult {
        roslamb,
        tgo: gt.tgo,
        roatm: roatm_corrected,
        ttatmg,
        satm,
        xrorayp,
    }
}
```

- [ ] **Step 5: Update lib.rs**

```rust
pub mod atmospheric;
pub mod constants;
pub mod gas_transmission;
pub mod geometry;
pub mod lut;
pub mod rayleigh;
pub mod sensor;
pub mod utils;
```

- [ ] **Step 6: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test atmospheric`
Expected: All tests PASS

- [ ] **Step 7: Commit**

```bash
git add lasrc_core/src/atmospheric.rs lasrc_core/src/lib.rs
git commit -m "feat: implement atmospheric correction (atmcorlamb2 and atmcorlamb2_new)"
```

---

## Task 8: Aerosol Retrieval

**Files:**
- Create: `lasrc_core/src/aerosol.rs`
- Modify: `lasrc_core/src/lib.rs`

Porting `subaeroret_new` and `subaeroret` from `subaeroret.c`, and the aerosol interpolation from `aero_interp.c`.

- [ ] **Step 1: Write structs and failing tests**

Add to `lasrc_core/src/aerosol.rs`:

```rust
//! Aerosol optical thickness retrieval and spatial interpolation.

use crate::atmospheric::{AtmCorrCoefficients, atmcorlamb2_new};
use crate::constants::*;

/// Result of aerosol retrieval for a single window.
#[derive(Debug, Clone)]
pub struct AerosolResult {
    pub raot: f64,      // retrieved AOT at 550nm
    pub residual: f64,  // model fit residual
    pub eps: f64,       // Angstrom coefficient
    pub iaots: usize,   // AOT index (warm start for next pixel)
}

/// Retrieve aerosol optical thickness using the semi-empirical method.
///
/// Iteratively searches AOT values to minimize residual between observed
/// and modeled band ratios.
pub fn subaeroret_new(
    is_water: bool,
    iband1: usize,           // primary band index (coastal aerosol)
    erelc: &[f64],           // expected band ratio coefficients
    troatm: &[f64],          // TOA reflectance per band
    tgo_arr: &[f64],         // gas transmission per band
    roatm_ia_max: &[f64],    // atmospheric reflectance at max AOT
    roatm_coef: &[[f64; NCOEF]], // polynomial coeff per band
    ttatmg_coef: &[[f64; NCOEF]],
    satm_coef: &[[f64; NCOEF]],
    normext_p0a3: &[f64],    // normext at pressure=0, AOT index=3
    lambda: &[f64],
    eps: f64,
    iaots: usize,            // starting AOT index
    tth: &[f64],             // threshold for negative SR test
) -> AerosolResult {
    todo!()
}

/// Bilinear interpolation of aerosol values from window centers to all pixels.
/// For Landsat: interpolates from center of 3x3 window grid.
pub fn aerosol_interp(
    taero: &mut [f64],       // AOT array (nlines x nsamps), centers filled
    teps: &mut [f64],        // eps array (nlines x nsamps), centers filled
    ipflag: &mut [u8],       // QA flags
    qa_band: &[u16],         // input QA for fill detection
    nlines: usize,
    nsamps: usize,
    aero_window: usize,
    half_aero_window: usize,
) {
    todo!()
}

/// Fix invalid aerosol retrievals using local averaging.
/// Multi-pass: forward, forward with relaxed threshold, reverse.
pub fn fix_invalid_aerosols(
    taero: &mut [f64],
    teps: &mut [f64],
    ipflag: &mut [u8],
    nlines: usize,
    nsamps: usize,
    aero_window: usize,
    half_aero_window: usize,
    fix_aero_window: usize,
    half_fix_aero_window: usize,
    min_clear_pix: usize,
) {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subaeroret_new_clear_pixel() {
        // Simulate a clear land pixel with moderate vegetation
        let nband = 7;
        let erelc = vec![1.0, 0.8, 0.0, 0.5, 0.0, 0.0, 0.2];
        let troatm = vec![0.15, 0.14, 0.12, 0.11, 0.08, 0.04, 0.03];
        let tgo_arr = vec![0.95; nband];
        let roatm_ia_max = vec![0.4; nband];

        let roatm_coef: Vec<[f64; NCOEF]> =
            (0..nband).map(|_| [0.0, 0.0, 0.15, 0.01]).collect();
        let ttatmg_coef: Vec<[f64; NCOEF]> =
            (0..nband).map(|_| [0.0, 0.0, -0.1, 0.85]).collect();
        let satm_coef: Vec<[f64; NCOEF]> =
            (0..nband).map(|_| [0.0, 0.0, 0.05, 0.005]).collect();
        let normext_p0a3 = vec![1.0; nband];
        let tth = vec![1.0e-3, 1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-4];

        let result = subaeroret_new(
            false, 0, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
            &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
            &LAMBDA_LANDSAT, 1.5, 0, &tth,
        );

        assert!(result.raot >= 0.01, "AOT should be >= 0.01: {}", result.raot);
        assert!(result.raot <= 5.0, "AOT should be <= 5.0: {}", result.raot);
        assert!(result.residual >= 0.0, "Residual should be non-negative");
    }

    #[test]
    fn test_aerosol_interp_fills_gaps() {
        let nlines = 9;
        let nsamps = 9;
        let npix = nlines * nsamps;
        let mut taero = vec![0.0f64; npix];
        let mut teps = vec![0.0f64; npix];
        let mut ipflag = vec![0u8; npix];
        let qa_band = vec![1u16; npix]; // all valid (non-fill)

        // Set center pixels of 3x3 windows
        for i in (1..nlines).step_by(3) {
            for j in (1..nsamps).step_by(3) {
                let pix = i * nsamps + j;
                taero[pix] = 0.1;
                teps[pix] = 1.5;
                ipflag[pix] = 1 << IPFLAG_CLEAR;
            }
        }

        aerosol_interp(&mut taero, &mut teps, &mut ipflag, &qa_band, nlines, nsamps, 3, 1);

        // All non-fill pixels should now have aerosol values
        for i in 0..npix {
            if qa_band[i] != 0 {
                assert!(taero[i] > 0.0, "Pixel {i} should have interpolated aerosol");
            }
        }
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test aerosol`
Expected: FAIL with "not yet implemented"

- [ ] **Step 3: Implement subaeroret_new**

```rust
pub fn subaeroret_new(
    is_water: bool,
    iband1: usize,
    erelc: &[f64],
    troatm: &[f64],
    tgo_arr: &[f64],
    roatm_ia_max: &[f64],
    roatm_coef: &[[f64; NCOEF]],
    ttatmg_coef: &[[f64; NCOEF]],
    satm_coef: &[[f64; NCOEF]],
    normext_p0a3: &[f64],
    lambda: &[f64],
    eps: f64,
    iaots: usize,
    tth: &[f64],
) -> AerosolResult {
    let nbands = erelc.len();
    let mut iaot = iaots.min(NAOT_VALS - 1);
    let mut residual1 = 2000.0f64;
    let mut residual2 = 1000.0f64;
    let mut raot1 = AOT550NM[0];
    let mut raot2 = AOT550NM[0];
    let mut raot3;
    let mut residual3;
    let mut best_raot = AOT550NM[0];
    let mut best_residual = f64::MAX;

    // Iterate through AOT values to find minimum residual
    loop {
        if iaot >= NAOT_VALS {
            break;
        }

        let aot = AOT550NM[iaot];
        let mut testth = false;

        // Correct band1 first
        let coeff1 = AtmCorrCoefficients {
            roatm_upper: roatm_ia_max[iband1],
            roatm_coef: roatm_coef[iband1],
            ttatmg_coef: ttatmg_coef[iband1],
            satm_coef: satm_coef[iband1],
        };
        let ros1 = atmcorlamb2_new(&coeff1, tgo_arr[iband1], iband1, aot, normext_p0a3[iband1], troatm[iband1], lambda, eps);

        if ros1 - tth[iband1] < 0.0 {
            testth = true;
        }

        // Calculate residual across all active bands
        let mut residual = 0.0;
        let mut nbval = 0;
        for ib in 0..nbands {
            if erelc[ib] == 0.0 && !is_water {
                continue;
            }
            let coeff_ib = AtmCorrCoefficients {
                roatm_upper: roatm_ia_max[ib],
                roatm_coef: roatm_coef[ib],
                ttatmg_coef: ttatmg_coef[ib],
                satm_coef: satm_coef[ib],
            };
            let roslamb = atmcorlamb2_new(&coeff_ib, tgo_arr[ib], ib, aot, normext_p0a3[ib], troatm[ib], lambda, eps);

            if is_water {
                residual += roslamb * roslamb;
            } else {
                let err = roslamb - erelc[ib] * ros1;
                residual += err * err;
            }
            nbval += 1;
        }
        if nbval > 0 {
            residual = (residual).sqrt() / nbval as f64;
        }

        if residual < best_residual {
            best_residual = residual;
            best_raot = aot;
        }

        if testth {
            break;
        }

        // Shift history
        residual3 = residual2;
        raot3 = raot2;
        residual2 = residual1;
        raot2 = raot1;
        residual1 = residual;
        raot1 = aot;

        if residual >= residual2 {
            break; // Residual increasing: passed minimum
        }

        iaot += 1;
    }

    // Parabolic refinement using 3 most recent points
    if iaot > 1 && residual2 < residual1 {
        // We have 3 points: (raot3, residual3), (raot2, residual2), (raot1, residual1)
        // Fit parabola: res = a*raot^2 + b*raot + c
        let x1 = raot3;
        let x2 = raot2;
        let x3 = raot1;
        let y1 = residual3;
        let y2 = residual2;
        let y3 = residual1;

        let denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
        if denom.abs() > 1e-20 {
            let a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
            let b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;

            if a > 0.0 {
                let raotmin = -b / (2.0 * a);
                if raotmin >= 0.01 && raotmin <= 4.0 {
                    // Evaluate at refined point
                    let coeff1 = AtmCorrCoefficients {
                        roatm_upper: roatm_ia_max[iband1],
                        roatm_coef: roatm_coef[iband1],
                        ttatmg_coef: ttatmg_coef[iband1],
                        satm_coef: satm_coef[iband1],
                    };
                    let ros1 = atmcorlamb2_new(&coeff1, tgo_arr[iband1], iband1, raotmin, normext_p0a3[iband1], troatm[iband1], lambda, eps);
                    let mut res = 0.0;
                    let mut nbval = 0;
                    for ib in 0..nbands {
                        if erelc[ib] == 0.0 && !is_water { continue; }
                        let coeff_ib = AtmCorrCoefficients {
                            roatm_upper: roatm_ia_max[ib],
                            roatm_coef: roatm_coef[ib],
                            ttatmg_coef: ttatmg_coef[ib],
                            satm_coef: satm_coef[ib],
                        };
                        let roslamb = atmcorlamb2_new(&coeff_ib, tgo_arr[ib], ib, raotmin, normext_p0a3[ib], troatm[ib], lambda, eps);
                        if is_water {
                            res += roslamb * roslamb;
                        } else {
                            let err = roslamb - erelc[ib] * ros1;
                            res += err * err;
                        }
                        nbval += 1;
                    }
                    if nbval > 0 {
                        res = res.sqrt() / nbval as f64;
                    }
                    if res < best_residual {
                        best_raot = raotmin;
                        best_residual = res;
                    }
                }
            }
        }
    }

    AerosolResult {
        raot: best_raot,
        residual: best_residual,
        eps,
        iaots: iaot.min(NAOT_VALS - 1),
    }
}
```

- [ ] **Step 4: Implement aerosol_interp (bilinear from window centers)**

```rust
pub fn aerosol_interp(
    taero: &mut [f64],
    teps: &mut [f64],
    ipflag: &mut [u8],
    qa_band: &[u16],
    nlines: usize,
    nsamps: usize,
    aero_window: usize,
    half_aero_window: usize,
) {
    let aero_step = 1.0 / aero_window as f64;

    for line in 0..nlines {
        for samp in 0..nsamps {
            let curr_pix = line * nsamps + samp;

            // Skip fill pixels
            if qa_band[curr_pix] == 0 {
                ipflag[curr_pix] = 1 << IPFLAG_FILL;
                continue;
            }

            // Find the surrounding window centers
            let center_line = ((line as f64 * aero_step) as usize) * aero_window + half_aero_window;
            let center_samp = ((samp as f64 * aero_step) as usize) * aero_window + half_aero_window;

            // Skip pixels that are already window centers
            if line == center_line && samp == center_samp {
                continue;
            }

            // Determine 4 surrounding centers for bilinear interp
            let yaero = (line as f64 - center_line as f64) * aero_step;
            let xaero = (samp as f64 - center_samp as f64) * aero_step;

            let l0 = center_line;
            let l1 = if yaero >= 0.0 { (l0 + aero_window).min(nlines - 1) } else { l0.saturating_sub(aero_window) };
            let s0 = center_samp;
            let s1 = if xaero >= 0.0 { (s0 + aero_window).min(nsamps - 1) } else { s0.saturating_sub(aero_window) };

            let u = yaero.abs();
            let v = xaero.abs();

            let p00 = l0 * nsamps + s0;
            let p01 = l0 * nsamps + s1;
            let p10 = l1 * nsamps + s0;
            let p11 = l1 * nsamps + s1;

            let a00 = taero[p00.min(taero.len() - 1)];
            let a01 = taero[p01.min(taero.len() - 1)];
            let a10 = taero[p10.min(taero.len() - 1)];
            let a11 = taero[p11.min(taero.len() - 1)];

            taero[curr_pix] = a00 + u * (a10 - a00) + v * (a01 - a00)
                + u * v * (a00 - a01 - a10 + a11);

            let e00 = teps[p00.min(teps.len() - 1)];
            let e01 = teps[p01.min(teps.len() - 1)];
            let e10 = teps[p10.min(teps.len() - 1)];
            let e11 = teps[p11.min(teps.len() - 1)];

            teps[curr_pix] = e00 + u * (e10 - e00) + v * (e01 - e00)
                + u * v * (e00 - e01 - e10 + e11);

            ipflag[curr_pix] |= 1 << IPFLAG_INTERP_WINDOW;
        }
    }
}
```

- [ ] **Step 5: Implement fix_invalid_aerosols (3-pass)**

```rust
pub fn fix_invalid_aerosols(
    taero: &mut [f64],
    teps: &mut [f64],
    ipflag: &mut [u8],
    nlines: usize,
    nsamps: usize,
    aero_window: usize,
    half_aero_window: usize,
    fix_aero_window: usize,
    half_fix_aero_window: usize,
    min_clear_pix: usize,
) {
    let fill_with_avg = |taero: &mut [f64], teps: &mut [f64], ipflag: &mut [u8],
                          nlines: usize, nsamps: usize, aero_window: usize,
                          half_aero_window: usize, half_fix: usize,
                          required_clear: usize, use_filled: bool, reverse: bool|
                          -> usize {
        let mut unfilled = 0usize;

        let lines: Vec<usize> = if reverse {
            (0..nlines).rev().step_by(aero_window).collect()
        } else {
            (half_aero_window..nlines).step_by(aero_window).collect()
        };

        for &i in &lines {
            let samps: Vec<usize> = if reverse {
                (0..nsamps).rev().step_by(aero_window).collect()
            } else {
                (half_aero_window..nsamps).step_by(aero_window).collect()
            };
            for &j in &samps {
                let curr_pix = i * nsamps + j;
                if ipflag[curr_pix] & (1 << IPFLAG_FILL) != 0 {
                    continue;
                }
                if ipflag[curr_pix] & (1 << IPFLAG_CLEAR) != 0 {
                    continue;
                }
                if use_filled && ipflag[curr_pix] & (1 << IPFLAG_FIXED) != 0 {
                    continue;
                }

                // Search in fix window
                let mut sum_aero = 0.0;
                let mut sum_eps = 0.0;
                let mut count = 0usize;
                let li_start = i.saturating_sub(half_fix * aero_window);
                let li_end = (i + half_fix * aero_window + 1).min(nlines);
                let si_start = j.saturating_sub(half_fix * aero_window);
                let si_end = (j + half_fix * aero_window + 1).min(nsamps);

                let mut li = li_start + half_aero_window;
                while li < li_end {
                    let mut si = si_start + half_aero_window;
                    while si < si_end {
                        let pix = li * nsamps + si;
                        if pix < ipflag.len() {
                            let is_valid = ipflag[pix] & (1 << IPFLAG_CLEAR) != 0;
                            let is_filled = use_filled && ipflag[pix] & (1 << IPFLAG_FIXED) != 0;
                            if is_valid || is_filled {
                                sum_aero += taero[pix];
                                sum_eps += teps[pix];
                                count += 1;
                            }
                        }
                        si += aero_window;
                    }
                    li += aero_window;
                }

                if count >= required_clear {
                    taero[curr_pix] = sum_aero / count as f64;
                    teps[curr_pix] = sum_eps / count as f64;
                    ipflag[curr_pix] |= 1 << IPFLAG_FIXED;
                } else {
                    unfilled += 1;
                }
            }
        }
        unfilled
    };

    // Pass 1: Forward, require min_clear_pix, no filled pixels
    let unfilled = fill_with_avg(
        taero, teps, ipflag, nlines, nsamps, aero_window,
        half_aero_window, half_fix_aero_window, min_clear_pix, false, false,
    );

    if unfilled > 0 {
        // Pass 2: Forward, require 1 pixel, allow filled
        let unfilled = fill_with_avg(
            taero, teps, ipflag, nlines, nsamps, aero_window,
            half_aero_window, half_fix_aero_window, 1, true, false,
        );

        if unfilled > 0 {
            // Pass 3: Reverse direction
            fill_with_avg(
                taero, teps, ipflag, nlines, nsamps, aero_window,
                half_aero_window, half_fix_aero_window, 1, true, true,
            );
        }
    }
}
```

- [ ] **Step 6: Update lib.rs**

```rust
pub mod aerosol;
pub mod atmospheric;
pub mod constants;
pub mod gas_transmission;
pub mod geometry;
pub mod lut;
pub mod rayleigh;
pub mod sensor;
pub mod utils;
```

- [ ] **Step 7: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test aerosol`
Expected: All tests PASS

- [ ] **Step 8: Commit**

```bash
git add lasrc_core/src/aerosol.rs lasrc_core/src/lib.rs
git commit -m "feat: implement aerosol retrieval and spatial interpolation"
```

---

## Task 9: Top-Level Correction Orchestrator

**Files:**
- Create: `lasrc_core/src/correction.rs`
- Modify: `lasrc_core/src/lib.rs`

This ties together all the compute modules into the main `compute_surface_reflectance` function.

- [ ] **Step 1: Write struct definitions and failing test**

Add to `lasrc_core/src/correction.rs`:

```rust
//! Top-level surface reflectance computation orchestrating all sub-modules.

use ndarray::{Array2, ArrayView2};
use rayon::prelude::*;

use crate::aerosol::{AerosolResult, aerosol_interp, fix_invalid_aerosols, subaeroret_new};
use crate::atmospheric::{AtmCorrCoefficients, AtmCorrResult, atmcorlamb2, atmcorlamb2_new};
use crate::constants::*;
use crate::gas_transmission::{GasCoefficients, compute_gas_transmission};
use crate::geometry::utm_to_deg;
use crate::lut::LookupTables;
use crate::sensor::Sensor;
use crate::utils::get_3rd_order_poly_coeff;

/// Auxiliary atmospheric data for the scene, loaded by Python.
pub struct AuxiliaryData {
    /// DEM elevation (meters). Dimensions: [DEM_NBLAT x DEM_NBLON].
    pub dem: Vec<i16>,
    /// Water vapor (scaled). Dimensions: [CMG_NBLAT x CMG_NBLON].
    pub wv: Vec<i16>,
    /// Ozone (scaled). Dimensions: [CMG_NBLAT x CMG_NBLON].
    pub oz: Vec<i16>,
    /// Band ratio arrays for aerosol retrieval.
    pub ratiob1: Vec<i16>,
    pub ratiob2: Vec<i16>,
    pub ratiob7: Vec<i16>,
    pub intratiob1: Vec<i16>,
    pub intratiob2: Vec<i16>,
    pub intratiob7: Vec<i16>,
    pub slpratiob1: Vec<i16>,
    pub slpratiob2: Vec<i16>,
    pub slpratiob7: Vec<i16>,
    pub andwi: Vec<i16>,
    pub sndwi: Vec<i16>,
    /// Water vapor and ozone scale factors [MODIS, VIIRS].
    pub wv_scale: f64,
    pub oz_scale: f64,
    pub wv_default: f64,
    pub oz_default: f64,
}

/// Complete surface reflectance result for a scene.
pub struct SurfaceReflectanceResult {
    /// Surface reflectance per band. Each is nlines x nsamps, scaled int16.
    pub sr_bands: Vec<Array2<i16>>,
    /// Brightness temperature per band (Landsat only). Scaled uint16.
    pub bt_bands: Vec<Array2<u16>>,
    /// Aerosol optical thickness. nlines x nsamps, scaled int16.
    pub aerosol: Array2<i16>,
    /// QA flags. nlines x nsamps.
    pub qa: Array2<u8>,
}

/// Main entry point: compute surface reflectance for an entire scene.
pub fn compute_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<f32>],   // TOA reflectance per band (float)
    bt_bands: &[ArrayView2<f32>],    // Brightness temperature (Landsat)
    solar_zenith: &ArrayView2<f32>,  // degrees
    solar_azimuth: &ArrayView2<f32>, // degrees
    view_zenith: &ArrayView2<f32>,   // degrees
    view_azimuth: &ArrayView2<f32>,  // degrees
    qa_band: &ArrayView2<u16>,       // level-1 QA
    lut: &LookupTables,
    aux: &AuxiliaryData,
    scene_center_lat: f64,
    scene_center_lon: f64,
    use_orig_aero: bool,
) -> SurfaceReflectanceResult {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_surface_reflectance_result_structure() {
        // Verify the result type can be constructed
        let nlines = 10;
        let nsamps = 10;
        let result = SurfaceReflectanceResult {
            sr_bands: vec![Array2::zeros((nlines, nsamps)); 8],
            bt_bands: vec![Array2::zeros((nlines, nsamps)); 2],
            aerosol: Array2::zeros((nlines, nsamps)),
            qa: Array2::zeros((nlines, nsamps)),
        };
        assert_eq!(result.sr_bands.len(), 8);
        assert_eq!(result.qa.shape(), &[nlines, nsamps]);
    }
}
```

- [ ] **Step 2: Implement compute_surface_reflectance**

This is a large function. The implementation should follow the C code flow:

1. Map scene center to CMG grid, extract atmospheric params (pressure, ozone, water vapor)
2. Compute gas transmission coefficients at scene center
3. Pre-compute polynomial coefficients for each band (for atmcorlamb2_new)
4. Initial atmospheric correction at scene center with default AOT
5. Aerosol retrieval loop over window centers (parallelized with rayon)
6. Fix invalid aerosols, interpolate
7. Final per-pixel correction with retrieved aerosol
8. Scale to output integers

Due to the size, this should be implemented incrementally across several sub-steps, tested at each stage. The full implementation is the integration of all previously built modules.

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
    scene_center_lat: f64,
    scene_center_lon: f64,
    use_orig_aero: bool,
) -> SurfaceReflectanceResult {
    let nlines = toa_bands[0].nrows();
    let nsamps = toa_bands[0].ncols();
    let npix = nlines * nsamps;
    let nrefl = sensor.num_refl_bands();
    let aero_window = sensor.aerosol_window();
    let half_aero_window = sensor.half_aerosol_window();
    let indices = sensor.band_indices();

    // --- Step 1: Scene-center atmospheric state ---
    let ycmg = (89.975 - scene_center_lat) * 20.0;
    let xcmg = (179.975 + scene_center_lon) * 20.0;
    let lcmg = ycmg.round() as usize;
    let scmg = xcmg.round() as usize;

    let dem_pix = lcmg * DEM_NBLON + scmg;
    let cmg_pix = lcmg * CMG_NBLON + scmg;

    let dem_val = if dem_pix < aux.dem.len() { aux.dem[dem_pix] } else { 0 };
    let pressure = if dem_val != -9999 {
        ATMOS_PRES_0 * (-(dem_val as f64) * ONE_DIV_8500).exp()
    } else {
        ATMOS_PRES_0
    };

    let uwv = if cmg_pix < aux.wv.len() && aux.wv[cmg_pix] != 0 {
        aux.wv[cmg_pix] as f64 / aux.wv_scale
    } else {
        aux.wv_default
    };

    let uoz = if cmg_pix < aux.oz.len() && aux.oz[cmg_pix] != 0 {
        aux.oz[cmg_pix] as f64 / aux.oz_scale
    } else {
        aux.oz_default
    };

    // Scene-center geometry
    let xts_center = solar_zenith[[nlines / 2, nsamps / 2]] as f64;
    let xmus_center = (xts_center * DEG2RAD).cos();

    // --- Step 2: Allocate output arrays ---
    let mut sr_f32: Vec<Vec<f32>> = (0..nrefl).map(|_| vec![0.0f32; npix]).collect();
    let mut taero = vec![DEFAULT_AERO; npix];
    let mut teps = vec![DEFAULT_EPS; npix];
    let mut ipflag = vec![0u8; npix];

    // Flatten QA for indexing
    let qa_flat: Vec<u16> = qa_band.iter().copied().collect();

    // --- Step 3: Aerosol inversion at window centers ---
    // (This would be the full aerosol retrieval loop)
    // For each window center, call subaeroret_new with the appropriate
    // band ratios, geometry, and atmospheric parameters.
    // This is the most complex part and ties together all sub-modules.

    // [Aerosol retrieval loop - uses subaeroret_new for each window center]
    // [Fix invalid aerosols - uses fix_invalid_aerosols]
    // [Interpolate aerosols - uses aerosol_interp]

    fix_invalid_aerosols(
        &mut taero, &mut teps, &mut ipflag,
        nlines, nsamps, aero_window, half_aero_window,
        sensor.fix_aerosol_window(), sensor.half_fix_aerosol_window(),
        sensor.min_clear_pix(),
    );

    aerosol_interp(
        &mut taero, &mut teps, &mut ipflag, &qa_flat,
        nlines, nsamps, aero_window, half_aero_window,
    );

    // --- Step 4: Final per-pixel correction ---
    // Apply atmospheric correction with retrieved per-pixel aerosol
    // using rayon for parallelism.

    // --- Step 5: Scale and pack output ---
    let sr_bands: Vec<Array2<i16>> = sr_f32.iter().map(|band| {
        let scaled: Vec<i16> = band.iter().map(|&v| {
            ((v as f64 + BAND_OFFSET_REFL) * MULT_FACTOR_REFL) as i16
        }).collect();
        Array2::from_shape_vec((nlines, nsamps), scaled).unwrap()
    }).collect();

    let bt_out: Vec<Array2<u16>> = bt_bands.iter().map(|bt| {
        let scaled: Vec<u16> = bt.iter().map(|&v| {
            ((v as f64 + BAND_OFFSET_TH) * MULT_FACTOR_TH) as u16
        }).collect();
        Array2::from_shape_vec((nlines, nsamps), scaled).unwrap()
    }).collect();

    let aerosol_arr: Vec<i16> = taero.iter().zip(qa_flat.iter()).map(|(&a, &qa)| {
        if qa == 0 { AERO_FILL } else { (a * MULT_FACTOR_AERO) as i16 }
    }).collect();

    SurfaceReflectanceResult {
        sr_bands,
        bt_bands: bt_out,
        aerosol: Array2::from_shape_vec((nlines, nsamps), aerosol_arr).unwrap(),
        qa: Array2::from_shape_vec((nlines, nsamps), ipflag).unwrap(),
    }
}
```

Note: The aerosol retrieval loop within `compute_surface_reflectance` is the most complex part and will need careful implementation matching the C code flow. The skeleton above shows the structure; the full loop body will be filled in during implementation following the patterns established in Tasks 5-8.

- [ ] **Step 3: Update lib.rs**

```rust
pub mod aerosol;
pub mod atmospheric;
pub mod constants;
pub mod correction;
pub mod gas_transmission;
pub mod geometry;
pub mod lut;
pub mod rayleigh;
pub mod sensor;
pub mod utils;
```

- [ ] **Step 4: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_core && cargo test`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add lasrc_core/src/correction.rs lasrc_core/src/lib.rs
git commit -m "feat: implement top-level surface reflectance orchestrator"
```

---

## Task 10: PyO3 Bindings

**Files:**
- Modify: `lasrc_py/src/lib.rs`
- Modify: `lasrc_py/Cargo.toml`

Expose the Rust compute functions to Python via PyO3 with NumPy array conversion.

- [ ] **Step 1: Implement PyO3 module with LUT and correction bindings**

Update `lasrc_py/src/lib.rs`:

```rust
use numpy::{PyArray2, PyReadonlyArray2, IntoPyArray};
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use lasrc_core::constants::*;
use lasrc_core::correction::{AuxiliaryData, SurfaceReflectanceResult, compute_surface_reflectance};
use lasrc_core::lut::LookupTables;
use lasrc_core::sensor::*;

#[pyclass]
struct PyLookupTables {
    inner: LookupTables,
}

#[pymethods]
impl PyLookupTables {
    #[new]
    fn new(
        rolutt: Vec<f64>,
        transt: Vec<f64>,
        sphalbt: Vec<f64>,
        normext: Vec<f64>,
        tsmax: Vec<f64>,
        tsmin: Vec<f64>,
        nbfic: Vec<f64>,
        nbfi: Vec<i32>,
        ttv: Vec<f64>,
        tts: Vec<f64>,
        nsr_bands: usize,
    ) -> PyResult<Self> {
        if tts.len() != NSOLAR_ZEN_VALS {
            return Err(PyValueError::new_err("tts must have NSOLAR_ZEN_VALS elements"));
        }
        let mut tts_arr = [0.0f64; NSOLAR_ZEN_VALS];
        tts_arr.copy_from_slice(&tts);
        Ok(Self {
            inner: LookupTables {
                rolutt, transt, sphalbt, normext,
                tsmax, tsmin, nbfic, nbfi, ttv, tts: tts_arr,
                nsr_bands,
            },
        })
    }
}

#[pyclass]
struct PyAuxiliaryData {
    inner: AuxiliaryData,
}

#[pymethods]
impl PyAuxiliaryData {
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        dem: Vec<i16>, wv: Vec<i16>, oz: Vec<i16>,
        ratiob1: Vec<i16>, ratiob2: Vec<i16>, ratiob7: Vec<i16>,
        intratiob1: Vec<i16>, intratiob2: Vec<i16>, intratiob7: Vec<i16>,
        slpratiob1: Vec<i16>, slpratiob2: Vec<i16>, slpratiob7: Vec<i16>,
        andwi: Vec<i16>, sndwi: Vec<i16>,
        wv_scale: f64, oz_scale: f64,
        wv_default: f64, oz_default: f64,
    ) -> Self {
        Self {
            inner: AuxiliaryData {
                dem, wv, oz,
                ratiob1, ratiob2, ratiob7,
                intratiob1, intratiob2, intratiob7,
                slpratiob1, slpratiob2, slpratiob7,
                andwi, sndwi,
                wv_scale, oz_scale, wv_default, oz_default,
            },
        }
    }
}

#[pyfunction]
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
    scene_center_lat: f64,
    scene_center_lon: f64,
    use_orig_aero: bool,
) -> PyResult<PyObject> {
    let sensor: Box<dyn Sensor> = match sensor_name {
        "LANDSAT_8" => Box::new(Landsat8),
        "LANDSAT_9" => Box::new(Landsat9),
        "SENTINEL_2A" => Box::new(Sentinel2A),
        "SENTINEL_2B" => Box::new(Sentinel2B),
        "SENTINEL_2C" => Box::new(Sentinel2C),
        _ => return Err(PyValueError::new_err(format!("Unknown sensor: {sensor_name}"))),
    };

    let toa_views: Vec<_> = toa_bands.iter().map(|a| a.as_array()).collect();
    let bt_views: Vec<_> = bt_bands.iter().map(|a| a.as_array()).collect();

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
        scene_center_lat,
        scene_center_lon,
        use_orig_aero,
    );

    // Return as dict of numpy arrays
    let dict = pyo3::types::PyDict::new(py);
    let sr_list: Vec<_> = result.sr_bands.into_iter()
        .map(|a| a.into_raw_vec_and_offset().0.into_pyarray(py))
        .collect();
    dict.set_item("sr_bands", sr_list)?;
    dict.set_item("aerosol", result.aerosol.into_raw_vec_and_offset().0.into_pyarray(py))?;
    dict.set_item("qa", result.qa.into_raw_vec_and_offset().0.into_pyarray(py))?;
    Ok(dict.into())
}

#[pymodule]
fn lasrc(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;
    m.add_class::<PyLookupTables>()?;
    m.add_class::<PyAuxiliaryData>()?;
    m.add_function(wrap_pyfunction!(process_surface_reflectance, m)?)?;
    Ok(())
}
```

- [ ] **Step 2: Verify build**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_py && pip install -e . 2>&1 | tail -5`
Expected: Builds successfully

- [ ] **Step 3: Quick smoke test from Python**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance && python -c "import lasrc; print(lasrc.__version__)"`
Expected: `0.1.0`

- [ ] **Step 4: Commit**

```bash
git add lasrc_py/src/lib.rs
git commit -m "feat: implement PyO3 bindings for LUT, auxiliary data, and surface reflectance"
```

---

## Task 11: Python Auxiliary Data Loading

**Files:**
- Create: `lasrc_py/python/lasrc/aux.py`
- Create: `lasrc_py/python/lasrc/sensors.py`

- [ ] **Step 1: Create sensors.py with sensor configurations**

```python
"""Sensor configuration dictionaries for I/O and metadata."""

SENSORS = {
    "LANDSAT_8": {
        "name": "LANDSAT_8",
        "satellite": "LANDSAT_8",
        "instrument": "OLI_TIRS",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band9",
        ],
        "thm_band_names": ["bt_band10", "bt_band11"],
        "input_band_names": [
            "band1", "band2", "band3", "band4",
            "band5", "band6", "band7", "band9",
        ],
        "input_thm_names": ["band10", "band11"],
        "resolution_m": 30.0,
        "nsr_bands": 8,
    },
    "LANDSAT_9": {
        "name": "LANDSAT_9",
        "satellite": "LANDSAT_9",
        "instrument": "OLI_TIRS",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band9",
        ],
        "thm_band_names": ["bt_band10", "bt_band11"],
        "input_band_names": [
            "band1", "band2", "band3", "band4",
            "band5", "band6", "band7", "band9",
        ],
        "input_thm_names": ["band10", "band11"],
        "resolution_m": 30.0,
        "nsr_bands": 8,
    },
    "SENTINEL_2A": {
        "name": "SENTINEL_2A",
        "satellite": "SENTINEL-2A",
        "instrument": "MSI",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band8",
            "sr_band8a", "sr_band11", "sr_band12",
        ],
        "thm_band_names": [],
        "input_band_names": [
            "B01", "B02", "B03", "B04", "B05", "B06",
            "B07", "B08", "B8A", "B11", "B12",
        ],
        "resolution_m": 10.0,
        "nsr_bands": 11,
    },
}
```

- [ ] **Step 2: Create aux.py with auxiliary data loading**

```python
"""Load auxiliary data from HDF4/5 files into NumPy arrays."""

import numpy as np


def load_lut_from_hdf(angle_hdf_path: str, intref_path: str,
                      transm_path: str, sphera_path: str,
                      nsr_bands: int) -> dict:
    """Load 6S lookup tables from HDF files.

    Returns dict with keys matching PyLookupTables constructor args.
    """
    import h5py

    with h5py.File(angle_hdf_path, "r") as f:
        tsmax = np.array(f["tsmax"], dtype=np.float64).ravel()
        tsmin = np.array(f["tsmin"], dtype=np.float64).ravel()
        nbfic = np.array(f["nbfic"], dtype=np.float64).ravel()
        nbfi = np.array(f["nbfi"], dtype=np.int32).ravel()
        ttv = np.array(f["ttv"], dtype=np.float64).ravel()
        tts = np.array(f["tts"], dtype=np.float64).ravel()

    with h5py.File(intref_path, "r") as f:
        rolutt = np.array(f["rolutt"], dtype=np.float64).ravel()

    with h5py.File(transm_path, "r") as f:
        transt = np.array(f["transt"], dtype=np.float64).ravel()

    with h5py.File(sphera_path, "r") as f:
        sphalbt = np.array(f["sphalbt"], dtype=np.float64).ravel()
        normext = np.array(f["normext"], dtype=np.float64).ravel()

    return {
        "rolutt": rolutt.tolist(),
        "transt": transt.tolist(),
        "sphalbt": sphalbt.tolist(),
        "normext": normext.tolist(),
        "tsmax": tsmax.tolist(),
        "tsmin": tsmin.tolist(),
        "nbfic": nbfic.tolist(),
        "nbfi": nbfi.tolist(),
        "ttv": ttv.tolist(),
        "tts": tts.tolist(),
        "nsr_bands": nsr_bands,
    }


def load_auxiliary_data(aux_path: str, aux_source: str = "VIIRS") -> dict:
    """Load DEM, water vapor, ozone, and band ratio auxiliary data.

    Args:
        aux_path: Path to auxiliary HDF file.
        aux_source: "VIIRS" or "MODIS".

    Returns dict with keys matching PyAuxiliaryData constructor args.
    """
    if aux_source == "VIIRS":
        return _load_viirs_aux(aux_path)
    else:
        return _load_modis_aux(aux_path)


def _load_viirs_aux(aux_path: str) -> dict:
    import h5py

    wv_scale = 200.0
    oz_scale = 400.0
    wv_default = 2.5
    oz_default = 0.3

    with h5py.File(aux_path, "r") as f:
        base = "/HDFEOS/GRIDS/VIIRS_CMG/Data Fields/"
        wv = np.array(f[base + "Average_Column_WV_CMG"], dtype=np.int16).ravel()
        oz = np.array(f[base + "Average_Column_Ozone_CMG"], dtype=np.int16).ravel()

    return {
        "wv": wv.tolist(),
        "oz": oz.tolist(),
        "wv_scale": wv_scale,
        "oz_scale": oz_scale,
        "wv_default": wv_default,
        "oz_default": oz_default,
    }


def load_dem(dem_path: str) -> list:
    """Load CMG DEM from HDF4 file."""
    from pyhdf.SD import SD, SDC

    hdf = SD(dem_path, SDC.READ)
    dem = hdf.select("CMG_DEM").get().astype(np.int16).ravel()
    hdf.end()
    return dem.tolist()


def load_ratio_file(ratio_path: str) -> dict:
    """Load NDWI and band ratio data from HDF4 file."""
    from pyhdf.SD import SD, SDC

    hdf = SD(ratio_path, SDC.READ)
    result = {}
    for name in ["ratiob1", "ratiob2", "ratiob7",
                  "intratiob1", "intratiob2", "intratiob7",
                  "slpratiob1", "slpratiob2", "slpratiob7",
                  "andwi", "sndwi"]:
        result[name] = np.array(hdf.select(name).get(), dtype=np.int16).ravel().tolist()
    hdf.end()
    return result
```

- [ ] **Step 3: Verify imports**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance && python -c "from lasrc.sensors import SENSORS; print(list(SENSORS.keys()))"`
Expected: `['LANDSAT_8', 'LANDSAT_9', 'SENTINEL_2A']`

- [ ] **Step 4: Commit**

```bash
git add lasrc_py/python/lasrc/aux.py lasrc_py/python/lasrc/sensors.py
git commit -m "feat: implement Python auxiliary data loading and sensor configs"
```

---

## Task 12: Python I/O and Pipeline

**Files:**
- Create: `lasrc_py/python/lasrc/io.py`
- Create: `lasrc_py/python/lasrc/pipeline.py`
- Create: `lasrc_py/python/lasrc/cli.py`
- Modify: `lasrc_py/python/lasrc/__init__.py`

- [ ] **Step 1: Create io.py**

```python
"""Scene I/O: read input scenes and write output products."""

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_bounds


def read_landsat_scene(scene_dir: str | Path) -> dict:
    """Read a Landsat Level-1 scene from GeoTIFF files.

    Returns dict with keys: toa_bands, bt_bands, qa_band, angles, metadata.
    """
    scene_dir = Path(scene_dir)
    # Find the MTL file
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
            # Apply gain/bias and solar angle correction
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

    # Read QA band
    qa_file = _find_band_file(band_files, "QA_PIXEL")
    with rasterio.open(qa_file) as src:
        qa_band = src.read(1).astype(np.uint16)

    # Read angle bands
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
            # Parse reflectance/radiance calibration
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

    nlines = result["sr_bands"][0].shape[0]
    nsamps = result["sr_bands"][0].shape[1]

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

    nlines = result["sr_bands"][0].shape[0]
    nsamps = result["sr_bands"][0].shape[1]
    nbands = len(result["sr_bands"])

    cog_profile = profile.copy()
    cog_profile.update(
        driver="GTiff",
        dtype="int16",
        count=nbands + 2,  # SR bands + aerosol + QA
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
```

- [ ] **Step 2: Create pipeline.py**

```python
"""Top-level scene processing pipeline."""

from pathlib import Path

import numpy as np

from lasrc.aux import load_auxiliary_data, load_dem, load_lut_from_hdf, load_ratio_file
from lasrc.io import read_landsat_scene, write_cog_output, write_espa_output
from lasrc.sensors import SENSORS

import lasrc as _lasrc


def process_scene(
    input_path: str | Path,
    aux_dir: str | Path,
    output_path: str | Path,
    sensor_name: str = "LANDSAT_8",
    output_format: str = "cog",
    use_orig_aero: bool = False,
    aux_source: str = "VIIRS",
) -> None:
    """Process a scene from TOA to surface reflectance.

    Args:
        input_path: Path to input scene directory.
        aux_dir: Path to auxiliary data directory.
        output_path: Path for output file/directory.
        sensor_name: One of LANDSAT_8, LANDSAT_9, SENTINEL_2A, etc.
        output_format: "cog" or "espa".
        use_orig_aero: Use original aerosol algorithm (slower, uses full LUT).
        aux_source: "VIIRS" or "MODIS".
    """
    sensor_config = SENSORS[sensor_name]
    aux_dir = Path(aux_dir)

    # Read input scene
    scene = read_landsat_scene(input_path)

    # Load LUT
    lut_data = load_lut_from_hdf(
        str(aux_dir / "ANGLE_NEW.hdf"),
        str(aux_dir / "RES_LUT_V3.0-LANDSAT.hdf"),
        str(aux_dir / "TRANS_LUT_V3.0-LANDSAT.hdf"),
        str(aux_dir / "AERO_LUT_V3.0-LANDSAT.hdf"),
        nsr_bands=sensor_config["nsr_bands"],
    )
    lut = _lasrc.PyLookupTables(**lut_data)

    # Load auxiliary data
    aux_data = load_auxiliary_data(
        str(aux_dir / "VIIRS_CMG_DAILY.hdf"),
        aux_source=aux_source,
    )
    dem = load_dem(str(aux_dir / "CMGDEM.hdf"))
    ratios = load_ratio_file(str(aux_dir / "ratiomapndwiexp.hdf"))
    aux = _lasrc.PyAuxiliaryData(dem=dem, **aux_data, **ratios)

    # Get scene center lat/lon from profile
    profile = scene["profile"]
    transform = profile["transform"]
    nlines = profile["height"]
    nsamps = profile["width"]
    center_x = transform.c + nsamps / 2 * transform.a
    center_y = transform.f + nlines / 2 * transform.e

    # Prepare angle arrays
    angles = scene["angles"]
    if "sza" in angles:
        sza = angles["sza"]
        saa = angles["saa"]
        vza = angles["vza"]
        vaa = angles["vaa"]
    else:
        sza = np.full((nlines, nsamps), 90.0 - angles["sun_elevation"], dtype=np.float32)
        saa = np.full((nlines, nsamps), 0.0, dtype=np.float32)
        vza = np.full((nlines, nsamps), 0.0, dtype=np.float32)
        vaa = np.full((nlines, nsamps), 0.0, dtype=np.float32)

    # Run Rust surface reflectance computation
    result = _lasrc.process_surface_reflectance(
        sensor_name=sensor_name,
        toa_bands=scene["toa_bands"],
        bt_bands=scene["bt_bands"],
        solar_zenith=sza,
        solar_azimuth=saa,
        view_zenith=vza,
        view_azimuth=vaa,
        qa_band=scene["qa_band"],
        lut=lut,
        aux=aux,
        scene_center_lat=center_y,
        scene_center_lon=center_x,
        use_orig_aero=use_orig_aero,
    )

    # Reshape flat arrays back to 2D
    for i in range(len(result["sr_bands"])):
        result["sr_bands"][i] = result["sr_bands"][i].reshape(nlines, nsamps)
    result["aerosol"] = result["aerosol"].reshape(nlines, nsamps)
    result["qa"] = result["qa"].reshape(nlines, nsamps)

    # Write output
    if output_format == "espa":
        write_espa_output(output_path, result, scene["metadata"], profile, sensor_config)
    else:
        write_cog_output(output_path, result, scene["metadata"], profile, sensor_config)
```

- [ ] **Step 3: Create cli.py**

```python
"""Command-line interface for LaSRC surface reflectance processing."""

import click


@click.command()
@click.option("--input", "input_path", required=True, type=click.Path(exists=True),
              help="Path to input scene directory")
@click.option("--aux-dir", required=True, type=click.Path(exists=True),
              help="Path to auxiliary data directory")
@click.option("--output", "output_path", required=True, type=click.Path(),
              help="Path for output file or directory")
@click.option("--sensor", default=None,
              type=click.Choice(["LANDSAT_8", "LANDSAT_9", "SENTINEL_2A", "SENTINEL_2B", "SENTINEL_2C"]),
              help="Sensor name (auto-detected if omitted)")
@click.option("--output-format", default="cog", type=click.Choice(["cog", "espa"]),
              help="Output format (default: cog)")
@click.option("--use-orig-aero", is_flag=True, default=False,
              help="Use original aerosol algorithm (slower)")
@click.option("--aux-source", default="VIIRS", type=click.Choice(["VIIRS", "MODIS"]),
              help="Auxiliary data source (default: VIIRS)")
def main(input_path, aux_dir, output_path, sensor, output_format, use_orig_aero, aux_source):
    """LaSRC: Compute surface reflectance from satellite imagery."""
    from lasrc.pipeline import process_scene

    if sensor is None:
        sensor = _detect_sensor(input_path)

    click.echo(f"Processing {input_path} with sensor {sensor}")
    process_scene(
        input_path=input_path,
        aux_dir=aux_dir,
        output_path=output_path,
        sensor_name=sensor,
        output_format=output_format,
        use_orig_aero=use_orig_aero,
        aux_source=aux_source,
    )
    click.echo(f"Output written to {output_path}")


def _detect_sensor(input_path: str) -> str:
    """Auto-detect sensor from input scene metadata."""
    from pathlib import Path
    p = Path(input_path)
    name = p.name.upper()
    if "LC08" in name or "LO08" in name:
        return "LANDSAT_8"
    elif "LC09" in name or "LO09" in name:
        return "LANDSAT_9"
    elif "S2A" in name:
        return "SENTINEL_2A"
    elif "S2B" in name:
        return "SENTINEL_2B"
    elif "S2C" in name:
        return "SENTINEL_2C"
    else:
        raise click.ClickException(f"Cannot detect sensor from {name}. Use --sensor flag.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Update __init__.py**

```python
from lasrc.lasrc import __version__, PyLookupTables, PyAuxiliaryData, process_surface_reflectance

__all__ = [
    "__version__",
    "PyLookupTables",
    "PyAuxiliaryData",
    "process_surface_reflectance",
]
```

- [ ] **Step 5: Verify CLI installs**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_py && pip install -e . && lasrc --help`
Expected: Shows usage with --input, --aux-dir, --output options

- [ ] **Step 6: Commit**

```bash
git add lasrc_py/python/lasrc/
git commit -m "feat: implement Python I/O, pipeline orchestration, and CLI"
```

---

## Task 13: Integration Testing Setup

**Files:**
- Create: `tests/test_geometry.py`
- Create: `tests/test_pipeline.py`
- Create: `tests/conftest.py`

- [ ] **Step 1: Create conftest.py with shared fixtures**

```python
"""Shared test fixtures for integration tests."""

import numpy as np
import pytest


@pytest.fixture
def small_scene():
    """Generate a small synthetic test scene (10x10 pixels)."""
    nlines, nsamps = 10, 10
    rng = np.random.default_rng(42)

    toa_bands = [rng.uniform(0.05, 0.3, (nlines, nsamps)).astype(np.float32) for _ in range(8)]
    bt_bands = [rng.uniform(270.0, 310.0, (nlines, nsamps)).astype(np.float32) for _ in range(2)]
    qa_band = np.ones((nlines, nsamps), dtype=np.uint16)
    sza = np.full((nlines, nsamps), 30.0, dtype=np.float32)
    saa = np.full((nlines, nsamps), 150.0, dtype=np.float32)
    vza = np.full((nlines, nsamps), 0.0, dtype=np.float32)
    vaa = np.full((nlines, nsamps), 0.0, dtype=np.float32)

    return {
        "toa_bands": toa_bands,
        "bt_bands": bt_bands,
        "qa_band": qa_band,
        "sza": sza,
        "saa": saa,
        "vza": vza,
        "vaa": vaa,
        "nlines": nlines,
        "nsamps": nsamps,
    }
```

- [ ] **Step 2: Create test_geometry.py**

```python
"""Validate geometry functions via Python bindings."""

import numpy as np
import pytest


def test_import():
    """Verify the lasrc package imports correctly."""
    import lasrc
    assert hasattr(lasrc, "__version__")


def test_version():
    import lasrc
    assert lasrc.__version__ == "0.1.0"
```

- [ ] **Step 3: Create test_pipeline.py (skeleton for future end-to-end tests)**

```python
"""End-to-end pipeline validation tests.

These tests require actual test data and auxiliary files.
Skip if test data is not available.
"""

import pytest
from pathlib import Path

TEST_DATA_DIR = Path(__file__).parent / "data"
HAS_TEST_DATA = TEST_DATA_DIR.exists() and any(TEST_DATA_DIR.iterdir()) if TEST_DATA_DIR.exists() else False


@pytest.mark.skipif(not HAS_TEST_DATA, reason="Test data not available")
def test_landsat_end_to_end():
    """Process a Landsat test scene and compare output to C reference."""
    pass  # Will be implemented when test data is available


@pytest.mark.skipif(not HAS_TEST_DATA, reason="Test data not available")
def test_output_matches_c_reference():
    """Compare Rust+Python output to C output for numerical validation."""
    pass  # Will compare ESPA-format outputs
```

- [ ] **Step 4: Run tests**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance && python -m pytest tests/ -v`
Expected: test_import and test_version PASS, end-to-end tests SKIPPED

- [ ] **Step 5: Commit**

```bash
git add tests/
git commit -m "feat: add integration test scaffolding"
```

---

## Task 14: Workspace Cargo.toml and Final Build Verification

**Files:**
- Create: `Cargo.toml` (workspace root)
- Modify: `.gitignore`

- [ ] **Step 1: Create workspace Cargo.toml**

```toml
[workspace]
members = ["lasrc_core", "lasrc_py"]
resolver = "2"
```

- [ ] **Step 2: Update .gitignore for Rust and Python artifacts**

Append to `.gitignore`:

```
# Rust
/target/
lasrc_core/target/
lasrc_py/target/

# Python
__pycache__/
*.pyc
*.egg-info/
dist/
build/
.eggs/

# Virtual environments
.venv/
```

- [ ] **Step 3: Full build and test run**

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance && cargo test --workspace`
Expected: All Rust tests PASS

Run: `cd /Users/seanharkins/projects/espa-surface-reflectance/lasrc_py && pip install -e . && python -m pytest ../tests/ -v`
Expected: All Python tests PASS

- [ ] **Step 4: Commit**

```bash
git add Cargo.toml .gitignore
git commit -m "feat: add workspace Cargo.toml and update gitignore for Rust/Python"
```
