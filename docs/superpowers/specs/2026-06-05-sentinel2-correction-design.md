# Sentinel-2 Surface Reflectance Correction — Design Spec

## Goal

Add Sentinel-2 support to the Rust+Python LaSRC port, reusing the existing core modules (LUT interpolation, atmospheric correction, aerosol retrieval, gas transmission, Rayleigh, geometry) with Sentinel-specific configuration and orchestration. Process all 13 bands (including B09 and B10). Use the semi-empirical aerosol approach only.

## Architecture

The Sentinel-2 path shares all Rust core algorithms with Landsat. The differences are configuration (band count, coefficients, window sizes) and orchestration (fill detection, aerosol post-processing, NDWI formula, band mapping).

### Files to modify

| File | Change |
|------|--------|
| `lasrc_core/src/constants.rs` | Add `TAURAY_SENTINEL_ALL` (13 bands), `LAMBDA_SENTINEL_ALL`, and 6 gas coefficient arrays for all 13 Sentinel bands |
| `lasrc_core/src/sensor.rs` | Fill in Sentinel `todo!()` stubs for `tauray()` and `gas_coefficients()`, update band config to 13 bands |
| `lasrc_core/src/correction.rs` | Add Sentinel orchestration: 13-band loop, ANY-band fill detection, B09 special handling, B8A/B12 NDWI, B04 reference band, 6x6 aerosol window, Sentinel aerosol post-processing, aerosol QA on B01 |
| `lasrc_core/src/aerosol.rs` | Add three Sentinel-specific functions: `aerosol_interp_sentinel`, `ipflag_expand_failed_sentinel`, `aero_avg_failed_sentinel` |

### Files to create

| File | Purpose |
|------|---------|
| `lasrc_py/python/lasrc/pipeline_sentinel.py` | Sentinel-2 pipeline: reads SAFE archive, resamples bands, calls Rust core, writes output |
| `lasrc_py/python/lasrc/io_sentinel.py` | SAFE archive I/O: read JP2 bands, extract metadata from MTD_MSIL1C.xml and MTD_TL.xml, resample 20m/60m to 10m |
| `lasrc_py/run_test_sentinel.py` | Test runner for the S2B test granule |

### Files to modify (Python)

| File | Change |
|------|---------|
| `lasrc_py/src/lib.rs` | Add PyO3 bindings for Sentinel sensor type and 13-band correction |
| `lasrc_py/compare_outputs.py` | Support Sentinel-2 band names and dimensions |

## Data Flow

### 1. Input reading (Python — `io_sentinel.py`)

Read directly from `.SAFE` archive:
- JP2 bands at native resolutions: B02/B03/B04/B08 at 10m (10980x10980), B05/B06/B07/B8A/B11/B12 at 20m (5490x5490), B01/B09/B10 at 60m (1830x1830)
- Resample 20m and 60m bands to 10m using nearest-neighbor (matching C's `convert_to_10m`)
- Extract from `MTD_TL.xml`: sun zenith/azimuth (scene-center mean), view zenith/azimuth (scene-center mean)
- Extract from `MTD_MSIL1C.xml` or tile metadata: EPSG code, UL corner coordinates, quantification value, radiometric offset
- Unscale TOA reflectance: `toa = (DN + offset) * (1/quantification_value)`
- Output: 13 flat float32 arrays, each 10980x10980

### 2. Fill detection (Rust)

A pixel is fill if ANY of its 13 bands equals the fill value (0 in the raw DN before unscaling). This differs from Landsat which uses QA_PIXEL bits. No cloud detection — Sentinel path does not set cloud QA bits.

### 3. Initial atmospheric correction (Rust)

For each of the 13 bands:
- Compute `atmcorlamb2` once at scene center with default AOT=0.05, eps=-1.0
- Exception: Band B09 (water vapor) uses `tgo=1.0, roatm=0.0, ttatmg=1.0, satm=0.0`
- Apply fast path to every non-fill pixel: `roslamb = (toa - tgo*roatm) / (tgo*ttatmg + satm*roslamb)`

### 4. Semi-empirical coefficient precomputation (Rust)

Same algorithm as Landsat, applied to all 13 bands:
- Loop 22 AOT values, call `atmcorlamb2` for each, store roatm/ttatmg/satm arrays
- Find monotonically increasing roatm range (roatm_iaMax)
- Fit 3rd-order polynomial coefficients
- Store per-band tgo, normext_p0a3

### 5. Aerosol retrieval (Rust)

Step through image in 6x6 windows (SAERO_WINDOW=6):
- Compute per-pixel lat/lon via `utm_to_deg`
- Look up CMG auxiliary data (bilinear interpolation): band ratios, slope/intercept
- NDWI calculation: `(B8A - B12*0.5) / (B8A + B12*0.5)` — uses the initial-corrected sband values
- Band ratios: `erelc[B01] = ndwi*slprb1 + intrb1`, `erelc[B02] = ndwi*slprb2 + intrb2`, `erelc[B04] = 1.0`, `erelc[B12] = ndwi*slprb7 + intrb7`
- TOA averaging within 6x6 window for bands B01, B02, B04, B12
- Reference band: B04 (red), not B05 like Landsat
- Run `subaeroret_new` at eps=1.0, 1.75, 2.5, find optimal eps via quadratic minimization
- Residual test: `residual < 0.015 + 0.005*corf + 0.10*troatm[B12]`
- NDVI validation: compute SR for B8A (`ros5`) and B04 (`ros4`), check `ros5 > 0.1 && (ros5-ros4)/(ros5+ros4) > 0`
- Water pixel retest: use B01/B04/B8A/B12 with `erelc` all set to 1.0, eps=1.5, residual test `> 0.010 + 0.005*corf || ros1 < 0`
- Copy aerosol/eps to all pixels within the 6x6 window

### 6. Aerosol post-processing (Rust — three new functions)

Replace Landsat's `fix_invalid_aerosols` with three Sentinel-specific functions:

**`aerosol_interp_sentinel`**: For each 6x6 window UL pixel, if it's not valid, interpolate from neighboring valid window-UL pixels. Uses the surrounding 8 window-UL neighbors.

**`ipflag_expand_failed_sentinel`**: Expand the "failed" flag to neighboring pixels — if any of the 8 neighbors of a valid pixel is failed, mark the valid pixel as failed too.

**`aero_avg_failed_sentinel`**: For each failed pixel, average the aerosol values from surrounding valid (non-failed, non-fill) pixels within a search window (SFIX_AERO_WINDOW). If not enough valid pixels found, use default aerosol.

### 7. Final per-pixel correction (Rust)

For each of the 13 bands, for each non-fill pixel:
- Call `atmcorlamb2_new` with per-pixel AOT and eps
- Band B10: just copy TOA values (no atmospheric correction)
- On B01: compute aerosol QA bits by comparing initial SR with final SR
- Clamp result to `[MIN_VALID_REFL, MAX_VALID_REFL]`

### 8. Output scaling and writing (Rust + Python)

- SR bands: `uint16 = clamp(roundf((sr + offset) * mult), 0, 65535)` — same as Landsat
- Aerosol: `int16 = clamp(roundf(aot * 1000), 0, 5000)`, fill = -9999
- Aerosol QA: `uint8` ipflag values directly
- All arrays 10980x10980
- Python writes to flat `.img` files in ESPA format

## Sentinel-2 Band Mapping

All 13 bands processed (PROC_ALL_BANDS equivalent):

| Index | Band | Wavelength (um) | Native Res | Notes |
|-------|------|-----------------|------------|-------|
| 0 | B01 | 0.443 | 60m | Coastal aerosol, aerosol QA reference |
| 1 | B02 | 0.490 | 10m | Blue |
| 2 | B03 | 0.560 | 10m | Green |
| 3 | B04 | 0.665 | 10m | Red, aerosol retrieval reference band |
| 4 | B05 | 0.705 | 20m | Red edge 1 |
| 5 | B06 | 0.740 | 20m | Red edge 2 |
| 6 | B07 | 0.783 | 20m | Red edge 3 |
| 7 | B08 | 0.842 | 10m | NIR |
| 8 | B8A | 0.865 | 20m | NIR narrow, used in NDVI/NDWI tests |
| 9 | B09 | 0.945 | 60m | Water vapor, special handling (tgo=1) |
| 10 | B10 | 1.375 | 60m | Cirrus, copy TOA to SR (no correction) |
| 11 | B11 | 1.610 | 20m | SWIR1 |
| 12 | B12 | 2.190 | 20m | SWIR2, used in NDWI and aerosol retrieval |

## Gas Coefficients (13 bands, from gascoef-msi.ASC)

```
oztransa:  [-0.00264691, -0.0272572, -0.0986512, -0.0500348, -0.0204295, -0.0108641, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]
wvtransa:  [2.29849e-27, 2.29849e-27, 0.000777307, 0.00361051, 0.0141249, 0.0137067, 0.00410217, 0.0285871, 0.000390755, 0.00001, 0.01, 0.000640155, 0.018006]
wvtransb:  [0.999742, 0.999742, 0.891099, 0.754895, 0.75596, 0.763497, 0.74117, 0.578722, 0.900899, 0.45818, 1.0, 0.943712, 0.647517]
ogtransa1: [4.91586e-20, 4.91586e-20, 4.91586e-20, 4.91586e-20, 5.3367e-06, 4.91586e-20, 9.03583e-05, 1.64109e-09, 1.90458e-05, 4.91586e-20, 7.62429e-06, 0.0212751, 0.0243065]
ogtransb0: [0.000197019, 0.000197019, 0.000197019, 0.000197019, -0.980313, 0.000197019, 0.0265393, 1.0e-10, 0.0322844, 0.000197019, 0.000197019, 0.000197019, 0.000197019]
ogtransb1: [9.57011e-16, 9.57011e-16, 9.57011e-16, 9.57011e-16, 1.33639, 9.57011e-16, 0.0532256, 1.0e-10, -0.0219907, 9.57011e-16, -0.216849, 0.0116062, 0.0604312]
```

## Tauray (13 bands, from tauray-msi.ASC)

```
[0.23432, 0.15106, 0.09102, 0.04535, 0.03584, 0.02924, 0.02338, 0.01847, 0.01560, 0.01092, 0.00243, 0.00128, 0.00037]
```

## Key Differences from Landsat Path

| Aspect | Landsat | Sentinel-2 |
|--------|---------|------------|
| Bands | 8 reflectance + 2 thermal | 13 reflectance, no thermal |
| Output resolution | 30m | 10m (resampled from mixed) |
| Fill detection | QA_PIXEL bit | ANY band == fill value |
| Cloud detection | Yes (cloud QA bits) | No |
| Aerosol window | LAERO_WINDOW (variable) | 6x6 fixed |
| Reference band | B05 (NIR) | B04 (red) |
| NDWI formula | (B5 - B7) / (B5 + B7) | (B8A - B12*0.5) / (B8A + B12*0.5) |
| Band ratios | erelc on B1, B2, B4=1.0, B7 | erelc on B01, B02, B04=1.0, B12 |
| Water test bands | B1/B4/B5/B7 | B01/B04/B8A/B12 |
| NDVI test | (ros5-ros4)/(ros5+ros4) > 0 | (ros5-ros4)/(ros5+ros4) > 0, ros5=B8A, ros4=B04 |
| Residual threshold | 0.015 + 0.005*corf | 0.015 + 0.005*corf + 0.10*troatm[B12] |
| Water residual | 0.010 + 0.005*corf | 0.010 + 0.005*corf |
| Aerosol post-processing | `fix_invalid_aerosols` (3-pass fill) | `aerosol_interp_sentinel` + `ipflag_expand_failed_sentinel` + `aero_avg_failed_sentinel` |
| Special bands | None | B09: tgo=1 bypass; B10: copy TOA |
| LUT directory | LDCMLUT/ | MSILUT/ |
| Input format | ESPA .img from Landsat | SAFE archive (JP2) |

## Testing and Verification

### C reference output

The C code has already been run on the same S2B granule. Reference output lives in:
```
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band1.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band2.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band3.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band4.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band5.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band6.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band7.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band8.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band8a.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band9.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band10.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band11.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_band12.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_aerosol.img
test_data/c_output/S2B_MSI_L1C_T38PNC_20260124_20260124_sr_aerosol_qa.img
```

All files are 10980x10980. SR bands are uint16, aerosol is int16 (fill=-9999), aerosol QA is uint8.

### Verification process

After each implementation step that produces output, run the comparison against the C reference:

1. **Run the Sentinel-2 pipeline** on the test granule:
   ```bash
   source .venv/bin/activate
   python lasrc_py/run_test_sentinel.py --format espa
   ```

2. **Compare Rust output against C reference** using the extended comparison script:
   ```bash
   python lasrc_py/compare_outputs.py --sensor sentinel \
       --rust-dir test_data/output_espa_s2 \
       --c-dir test_data/c_output \
       --c-prefix S2B_MSI_L1C_T38PNC_20260124_20260124_
   ```

3. **Inspect the comparison table** for all 15 output files. The table should report per-band: % of pixels that differ, median/mean/P95/P99/max absolute difference, and center-pixel difference.

4. **Iterate on discrepancies**: When a band shows unexpected diffs, investigate using the same approach as the Landsat port — check for f32/f64 precision mismatches, formula differences, band index mapping errors, and fill handling. Fix issues one at a time and re-compare after each fix.

### Acceptance criteria

The Rust output should match the C reference with accuracy comparable to the Landsat port:
- **Median diff**: 1-2 scaled integer units (0.0001-0.0002 reflectance)
- **P95 diff**: under 5 units
- **P99 diff**: under 30 units
- **Aerosol**: under 7% of pixels differing, median diff of 1-2
- **Aerosol QA**: under 1% of pixels differing
- Larger diffs in high-AOT/cloud regions are expected (same as Landsat) due to `subaeroret_new` convergence differences
- No systematic biases (center pixel should be exact or within 1-2 units)
