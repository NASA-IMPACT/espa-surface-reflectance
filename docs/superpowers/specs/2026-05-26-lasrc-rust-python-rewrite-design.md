# LaSRC Rust+Python Rewrite Design

## Overview

Rewrite the LaSRC atmospheric correction module (~14K lines of C) as a Rust core library with Python orchestration. LEDAPS is out of scope and remains unchanged.

The goal is a cleaner, more maintainable codebase that replaces the C+ESPA library stack with:
- **Rust** for performance-critical computation (LUT interpolation, aerosol retrieval, per-pixel atmospheric correction)
- **Python** for I/O, auxiliary data loading, orchestration, and CLI

## Scope

- **In scope**: LaSRC C code conversion (Landsat 8/9 and Sentinel-2 support)
- **Out of scope**: LEDAPS, auxiliary data download scripts (updatelads.py etc.), the 6S radiative transfer model itself

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| ESPA format dependency | Replace with COG (rasterio) | Decouples from ESPA C library stack |
| Auxiliary data format | Keep HDF4/5, read in Python | Avoids reformatting existing data |
| Rust/Python interface | PyO3 via maturin | Tight integration, pip-installable |
| Rust structure | Single crate with internal lib boundary | `lasrc_core` (pure Rust) wrapped by `lasrc_py` (PyO3) |
| Sensor abstraction | Unified trait with dynamic dispatch | Reduces duplication between Landsat/Sentinel correction code |
| Build approach | Bottom-up | Validate numerics at each layer before composing |
| Validation | Match C output first, then iterate | Phase correctness then improvement |

## Project Structure

```
espa-surface-reflectance/
├── lasrc_core/                    # Pure Rust library crate
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs                 # Public API
│       ├── sensor.rs              # Sensor trait + Landsat/Sentinel impls
│       ├── lut.rs                 # LUT data structures & interpolation
│       ├── aerosol.rs             # Aerosol retrieval & spatial interpolation
│       ├── atmospheric.rs         # Per-pixel atmospheric correction
│       ├── geometry.rs            # UTM/degree conversion, angle math
│       └── utils.rs               # Quick select, polynomial coefficients
│
├── lasrc_py/                      # PyO3 wrapper crate + Python package
│   ├── Cargo.toml                 # Depends on lasrc_core
│   ├── pyproject.toml             # maturin build config
│   ├── src/
│   │   └── lib.rs                 # PyO3 bindings, NumPy array conversion
│   └── python/
│       └── lasrc/
│           ├── __init__.py
│           ├── io.py              # Rasterio-based I/O (COG + ESPA format)
│           ├── aux.py             # Auxiliary data loading (HDF4/5 -> arrays)
│           ├── cli.py             # Command-line interface (click)
│           ├── pipeline.py        # Top-level orchestration
│           └── sensors.py         # Sensor configs (band names, resolutions)
│
├── tests/                         # Integration tests & validation data
│   ├── data/                      # Small test scenes
│   └── test_validation.py         # Numerical comparison with C output
│
├── ledaps/                        # Unchanged
├── scripts/                       # Unchanged
└── docs/
```

## Rust Core Library (`lasrc_core`)

### Sensor Trait

```rust
pub trait Sensor {
    fn band_count(&self) -> usize;
    fn reflectance_bands(&self) -> &[BandConfig];
    fn thermal_bands(&self) -> &[BandConfig];
    fn native_resolution(&self, band: usize) -> f64;
    fn output_resolution(&self) -> f64;
    fn aerosol_window_size(&self) -> usize;
    fn band_indices(&self) -> BandIndices;
}
```

Concrete types: `Landsat8`, `Landsat9`, `Sentinel2A`, `Sentinel2B`, `Sentinel2C`. Used via dynamic dispatch (`&dyn Sensor`) since sensor config lookups are not on the hot path.

### LUT Module (`lut.rs`)

Stores pre-computed 6S radiative transfer lookup tables. Current C structure: `[7 pressure levels][22 AOT values][8000+ solar zenith angles]`.

- Flat `Vec<f64>` storage with strides for efficient multi-dimensional indexing
- Trilinear interpolation across pressure, AOT, and angle dimensions
- Data loaded in Python from HDF, passed to Rust as arrays via `from_arrays()` constructor

```rust
pub struct LookupTable {
    data: Vec<f64>,
    pressure_levels: Vec<f64>,
    aot_values: Vec<f64>,
    angle_values: Vec<f64>,
}

impl LookupTable {
    pub fn from_arrays(data: &[f64], pressures: &[f64], aots: &[f64], angles: &[f64]) -> Self;
    pub fn interpolate(&self, pressure: f64, aot: f64, angle: f64) -> f64;
}
```

### Aerosol Retrieval (`aerosol.rs`)

- Takes per-pixel band ratios (NDVI, NDII from visible/SWIR), viewing geometry, and the LUT
- Returns aerosol optical thickness (AOT) per retrieval window (3x3 for Landsat, 6x6 for Sentinel)
- Spatial interpolation fills gaps between retrieval points
- Two algorithm paths: original full inversion and semi-empirical (default)

### Atmospheric Correction (`atmospheric.rs`)

The main compute function, operates on full-scene arrays:

```rust
pub fn compute_surface_reflectance(
    sensor: &dyn Sensor,
    toa_bands: &[ArrayView2<i16>],
    solar_zenith: &ArrayView2<f32>,
    solar_azimuth: &ArrayView2<f32>,
    view_zenith: &ArrayView2<f32>,
    view_azimuth: &ArrayView2<f32>,
    lut: &LookupTable,
    aux: &AuxiliaryData,
) -> SurfaceReflectanceResult;
```

Per-pixel correction formula:
`surface_refl = (TOA_refl - path_radiance) / (transmittance * spherical_albedo_term)`

Also generates QA flags (cloud, shadow, water, aerosol confidence).

Uses `rayon` for parallel iteration over pixels (replacing current OpenMP).

### Utilities (`geometry.rs`, `utils.rs`)

- UTM to lat/lon conversion (for locating auxiliary data pixels)
- Quick select (median finding for aerosol window statistics)
- Polynomial evaluation for spectral adjustments

## Python Layer

### I/O (`lasrc/io.py`)

- **Input**: Reads Landsat/Sentinel Level-1 GeoTIFF or ESPA XML+ENVI products via rasterio. Loads spectral bands, angle bands, and QA into NumPy arrays. Extracts metadata (gain/offset, sun elevation, acquisition date) from MTL (Landsat) or SAFE manifest (Sentinel).
- **Output (COG)**: Writes Cloud-Optimized GeoTIFF with proper CRS, transform, band descriptions, nodata values.
- **Output (ESPA)**: Writes ESPA internal format (flat binary bands + ENVI headers + XML metadata) using `numpy.ndarray.tofile()` and `lxml`. Matches C output structure exactly for validation. Enabled via `--output-format espa` CLI flag.

### Auxiliary Data (`lasrc/aux.py`)

Reads auxiliary datasets, returns NumPy arrays passed to Rust:

- **LUT**: 6S pre-computed tables from HDF, reshaped to Rust-expected layout
- **VIIRS/MODIS AOT**: Daily CMG aerosol grids, subsetted to scene extent
- **Water vapor / Ozone**: HDF grids, subsetted to scene extent
- **DEM**: Coarse CMGDEM (0.05 degree), subsetted and resampled to scene grid

Uses `h5py` for HDF5 and `pyhdf` for HDF4.

### Pipeline (`lasrc/pipeline.py`)

```python
def process_scene(input_path, aux_dir, output_path, sensor_name, output_format="cog"):
    sensor = get_sensor(sensor_name)
    metadata = read_metadata(input_path, sensor)
    toa_bands = read_spectral_bands(input_path, sensor, metadata)
    angles = read_angle_bands(input_path, sensor)
    lut = load_lut(aux_dir, sensor)
    aux_data = load_auxiliary(aux_dir, metadata)

    result = lasrc.compute_surface_reflectance(
        sensor, toa_bands, angles, lut, aux_data
    )

    if output_format == "cog":
        write_cog(output_path, result, metadata)
    else:
        write_espa(output_path, result, metadata)
```

### CLI (`lasrc/cli.py`)

Single entry point, auto-detects sensor from input metadata:

```
lasrc process --input /path/to/scene --aux-dir /path/to/aux --output /path/to/output [--output-format cog|espa]
```

Replaces the separate `do_lasrc_landsat.py` and `do_lasrc_sentinel.py` scripts.

## Dependencies

### Rust (`lasrc_core`)
- `ndarray` - Multi-dimensional arrays
- `rayon` - Parallel iteration

### Rust (`lasrc_py`)
- `pyo3` - Python bindings
- `numpy` (pyo3 crate) - Zero-copy NumPy array conversion

### Python
- `rasterio` - GeoTIFF/COG I/O
- `numpy` - Array handling
- `h5py` - HDF5 reading
- `pyhdf` - HDF4 reading
- `lxml` - XML metadata (ESPA format output)
- `click` - CLI
- `maturin` - Build system (dev)

### Build
- `maturin` with PyO3 backend
- Single `pip install .` from `lasrc_py/` builds everything

## Testing & Validation

### Unit Tests (Rust)
- **LUT interpolation**: Known values, boundary conditions, compared against C output
- **Aerosol retrieval**: Known band ratios + geometry -> AOT, compared against C
- **Atmospheric correction**: Single-pixel tests with known inputs/expected outputs from C
- **Utilities**: UTM/degree conversion, quick_select, polynomial evaluation

### Integration Tests (Python)
- **Round-trip I/O**: Read scene, write back, verify no data loss
- **End-to-end validation**: Process small test scene (~100x100 pixels) through both C and Rust+Python pipelines, compare output arrays
- **Tolerance**: Exact match for integer-scaled outputs (int16 reflectance). The C code's integer scaling absorbs minor floating point path differences.
- **Sensor parity**: Same scene through Landsat and Sentinel paths, verify consistency

### Validation Data
- Small cropped test scene included in repo (or script to generate one)
- Expected C output stored for regression testing
- CI runs comparison on every PR

## Build Order (Bottom-Up)

1. `lasrc_core` utilities (geometry, quick_select, polynomial)
2. `lasrc_core` LUT interpolation
3. `lasrc_core` aerosol retrieval
4. `lasrc_core` atmospheric correction (Landsat first, then Sentinel via trait)
5. `lasrc_py` PyO3 bindings
6. Python auxiliary data loading
7. Python I/O (input reading, ESPA output for validation)
8. Python pipeline + CLI
9. End-to-end validation against C output
10. COG output support
