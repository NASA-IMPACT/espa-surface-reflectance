//! Top-level correction orchestrator that ties together all compute modules
//! into the main `compute_surface_reflectance` function.
//!
//! Ported from the C LaSRC main processing loop.

use ndarray::{Array2, ArrayView2};

use crate::aerosol::{aerosol_interp, fix_invalid_aerosols, subaeroret_new};
use crate::atmospheric::{atmcorlamb2, atmcorlamb2_new, AtmCorrCoefficients};
use crate::constants::*;
use crate::gas_transmission::GasCoefficients;
use crate::geometry::{utm_to_deg, SpaceDef};
use crate::lut::LookupTables;
use crate::sensor::Sensor;
use crate::utils::get_3rd_order_poly_coeff;

/// Auxiliary atmospheric data for the scene, loaded by Python.
pub struct AuxiliaryData {
    /// DEM elevation (meters). [DEM_NBLAT x DEM_NBLON]
    pub dem: Vec<i16>,
    /// Water vapor (scaled). [CMG_NBLAT x CMG_NBLON]
    pub wv: Vec<i16>,
    /// Ozone (scaled). [CMG_NBLAT x CMG_NBLON]
    pub oz: Vec<i16>,
    /// Band ratio arrays for aerosol retrieval
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
    pub wv_scale: f64,
    pub oz_scale: f64,
    pub wv_default: f64,
    pub oz_default: f64,
}

/// Complete surface reflectance result for a scene.
pub struct SurfaceReflectanceResult {
    /// Surface reflectance per band, scaled int16
    pub sr_bands: Vec<Array2<i16>>,
    /// Brightness temperature (Landsat only), scaled uint16
    pub bt_bands: Vec<Array2<u16>>,
    /// AOT, scaled int16
    pub aerosol: Array2<i16>,
    /// QA flags
    pub qa: Array2<u8>,
}

/// CMG grid position with four surrounding cell indices and bilinear weights.
struct CmgPosition {
    /// Flat indices of the four surrounding cells: [row][col], [row][col+1],
    /// [row+1][col], [row+1][col+1].
    idx: [usize; 4],
    /// Bilinear weights for the four cells (sum to 1.0).
    w: [f64; 4],
}

/// Map a latitude/longitude to a CMG grid position with bilinear weights.
///
/// The CMG grid covers -90..90 lat and -180..180 lon at 0.05-degree resolution.
/// Uses cell-center coordinates (89.975, 179.975) to match the C code.
fn latlon_to_cmg(lat: f64, lon: f64) -> CmgPosition {
    // Fractional grid coordinates (cell-center aligned, matching C code)
    let ycmg = (89.975 - lat) * 20.0;
    let xcmg = (179.975 + lon) * 20.0;

    // Integer grid indices (truncation matches C code's (int) cast)
    let row = (ycmg as isize).clamp(0, (CMG_NBLAT - 1) as isize) as usize;
    let col = (xcmg as isize).clamp(0, (CMG_NBLON - 1) as isize) as usize;

    // Next row/col with wrapping at edges (matching C code)
    let row1 = if row >= CMG_NBLAT - 1 { 0 } else { row + 1 };
    let col1 = if col >= CMG_NBLON - 1 { 0 } else { col + 1 };

    // Fractional offsets for bilinear interpolation
    let u = ycmg - row as f64;
    let v = xcmg - col as f64;
    let u = u.clamp(0.0, 1.0);
    let v = v.clamp(0.0, 1.0);

    CmgPosition {
        idx: [
            row * CMG_NBLON + col,
            row * CMG_NBLON + col1,
            row1 * CMG_NBLON + col,
            row1 * CMG_NBLON + col1,
        ],
        w: [
            (1.0 - u) * (1.0 - v),
            (1.0 - u) * v,
            u * (1.0 - v),
            u * v,
        ],
    }
}

/// Bilinear interpolation over 4 values with the given weights.
fn bilerp(vals: [f64; 4], w: &[f64; 4]) -> f64 {
    vals[0] * w[0] + vals[1] * w[1] + vals[2] * w[2] + vals[3] * w[3]
}

/// Compute surface pressure from DEM elevation using the barometric formula.
fn pressure_from_elevation(elevation_m: f64) -> f64 {
    ATMOS_PRES_0 * (-elevation_m * ONE_DIV_8500).exp()
}

/// Extract atmospheric parameters (pressure, ozone, water vapor) at a given lat/lon
/// using bilinear interpolation over the 4 surrounding CMG grid cells.
fn extract_atm_params(aux: &AuxiliaryData, lat: f64, lon: f64) -> (f64, f64, f64) {
    let cmg = latlon_to_cmg(lat, lon);

    // DEM -> pressure: convert each corner's elevation to pressure, then interpolate.
    // Fill value (-9999) defaults to sea-level pressure.
    let pres_vals: [f64; 4] = std::array::from_fn(|i| {
        let idx = cmg.idx[i];
        if idx < aux.dem.len() && aux.dem[idx] != -9999 {
            pressure_from_elevation(aux.dem[idx] as f64)
        } else {
            ATMOS_PRES_0
        }
    });
    let pressure = bilerp(pres_vals, &cmg.w);

    // Water vapor: fill/zero values replaced with default DN before interpolation,
    // then unscaled after.
    let wv_vals: [f64; 4] = std::array::from_fn(|i| {
        let idx = cmg.idx[i];
        if idx < aux.wv.len() && aux.wv[idx] > 0 {
            aux.wv[idx] as f64
        } else {
            aux.wv_default * aux.wv_scale // default in DN space
        }
    });
    let uwv = bilerp(wv_vals, &cmg.w) / aux.wv_scale;

    // Ozone: same fill handling as water vapor.
    let oz_vals: [f64; 4] = std::array::from_fn(|i| {
        let idx = cmg.idx[i];
        if idx < aux.oz.len() && aux.oz[idx] > 0 {
            aux.oz[idx] as f64
        } else {
            aux.oz_default * aux.oz_scale
        }
    });
    let uoz = bilerp(oz_vals, &cmg.w) / aux.oz_scale;

    (pressure, uoz, uwv)
}

/// Look up a single band ratio at one CMG cell, falling back to slope/intercept.
fn ratio_at_cell(
    ratio: &[i16],
    intratio: &[i16],
    slpratio: &[i16],
    idx: usize,
    lat: f64,
) -> f64 {
    if idx < ratio.len() && ratio[idx] > 0 {
        ratio[idx] as f64 / 1000.0
    } else if idx < intratio.len() {
        let intr = intratio[idx] as f64 / 1000.0;
        let slp = if idx < slpratio.len() {
            slpratio[idx] as f64 / 1000.0
        } else {
            0.0
        };
        intr + slp * lat
    } else {
        0.5
    }
}

/// Look up band ratios from auxiliary data using bilinear interpolation.
///
/// Returns (rb1, rb2, rb7) as floats.
fn lookup_band_ratios(
    aux: &AuxiliaryData,
    lat: f64,
    lon: f64,
) -> (f64, f64, f64) {
    let cmg = latlon_to_cmg(lat, lon);

    let rb1_vals: [f64; 4] = std::array::from_fn(|i| {
        ratio_at_cell(&aux.ratiob1, &aux.intratiob1, &aux.slpratiob1, cmg.idx[i], lat)
    });
    let rb2_vals: [f64; 4] = std::array::from_fn(|i| {
        ratio_at_cell(&aux.ratiob2, &aux.intratiob2, &aux.slpratiob2, cmg.idx[i], lat)
    });
    let rb7_vals: [f64; 4] = std::array::from_fn(|i| {
        ratio_at_cell(&aux.ratiob7, &aux.intratiob7, &aux.slpratiob7, cmg.idx[i], lat)
    });

    (
        bilerp(rb1_vals, &cmg.w),
        bilerp(rb2_vals, &cmg.w),
        bilerp(rb7_vals, &cmg.w),
    )
}

/// Check if a pixel is water based on auxiliary andwi/sndwi grids.
///
/// Uses the primary (top-left) CMG cell for the discrete water/land decision,
/// matching the C code's use of ratio_pix11 for the NDWI threshold check.
fn is_water_pixel(aux: &AuxiliaryData, lat: f64, lon: f64) -> bool {
    let cmg = latlon_to_cmg(lat, lon);
    let idx = cmg.idx[0]; // primary cell

    let andwi_val = if idx < aux.andwi.len() { aux.andwi[idx] } else { 0 };
    let sndwi_val = if idx < aux.sndwi.len() { aux.sndwi[idx] } else { 0 };

    andwi_val > 0 && sndwi_val > 0
}

/// Build the expected band-ratio array (erelc) for aerosol retrieval.
///
/// For land pixels, erelc encodes the expected ratio of each band's surface
/// reflectance to the reference band. Bands not used are set to -1.0.
///
/// For Landsat, the reference band is the red band (index 3).
fn build_erelc(
    sensor: &dyn Sensor,
    rb1: f64,
    rb2: f64,
    rb7: f64,
    is_water: bool,
) -> Vec<f64> {
    let nbands = sensor.num_refl_bands();
    let bi = sensor.band_indices();
    let mut erelc = vec![-1.0; nbands];

    if is_water {
        // For water: use all active bands with zero ratios (minimize absolute SR)
        erelc[bi.coastal] = 0.0;
        erelc[bi.blue] = 0.0;
        erelc[bi.green] = 0.0;
        erelc[bi.red] = 0.0;
        erelc[bi.nir] = 0.0;
        erelc[bi.swir1] = 0.0;
        erelc[bi.swir2] = 0.0;
    } else {
        // For land: set band ratios relative to the red band
        // The red band is the reference (iband1), erelc for it is not used
        // directly but erelc for other bands encodes ratio to red.
        erelc[bi.coastal] = rb1;
        erelc[bi.blue] = rb2;
        erelc[bi.red] = 1.0; // reference band is always 1.0
        erelc[bi.swir1] = rb7;
    }

    erelc
}

/// Pre-compute polynomial coefficients for the fast atmospheric correction path.
///
/// For each band, evaluate `atmcorlamb2()` at multiple AOT values using the
/// scene-center geometry, extract roatm/ttatmg/satm values, then fit cubic
/// polynomials.
#[allow(clippy::too_many_arguments)]
fn precompute_coefficients(
    sensor: &dyn Sensor,
    lut: &LookupTables,
    gas_coeff: &[GasCoefficients],
    xts: f64,
    xtv: f64,
    xmus: f64,
    xmuv: f64,
    xfi: f64,
    cosxfi: f64,
    pressure: f64,
    uoz: f64,
    uwv: f64,
) -> (Vec<AtmCorrCoefficients>, Vec<f64>, Vec<f64>) {
    let nbands = sensor.num_refl_bands();
    let tauray = sensor.tauray();

    let mut coefficients = Vec::with_capacity(nbands);
    let mut tgo_arr = Vec::with_capacity(nbands);
    let mut normext_p0a3 = Vec::with_capacity(nbands);

    for iband in 0..nbands {
        // Arrays to collect roatm, ttatmg, satm at each AOT value
        let mut roatm_vals = vec![0.0f64; NAOT_VALS];
        let mut ttatmg_vals = vec![0.0f64; NAOT_VALS];
        let mut satm_vals = vec![0.0f64; NAOT_VALS];
        let mut tgo_band = 0.0f64;

        for iaot in 0..NAOT_VALS {
            let raot = AOT550NM[iaot];
            let result = atmcorlamb2(
                lut,
                &gas_coeff[iband],
                tauray[iband],
                iband,
                xts,
                xtv,
                xmus,
                xmuv,
                xfi,
                cosxfi,
                raot,
                pressure,
                uoz,
                uwv,
                0.0, // rotoa placeholder
                DEFAULT_EPS,
            );
            roatm_vals[iaot] = result.roatm;
            ttatmg_vals[iaot] = result.ttatmg;
            satm_vals[iaot] = result.satm;
            if iaot == 0 {
                tgo_band = result.tgo;
            }
        }

        // Fit cubic polynomials
        let roatm_coef = get_3rd_order_poly_coeff(&AOT550NM, &roatm_vals);
        let ttatmg_coef = get_3rd_order_poly_coeff(&AOT550NM, &ttatmg_vals);
        let satm_coef = get_3rd_order_poly_coeff(&AOT550NM, &satm_vals);

        // Find the last AOT index where roatm is still monotonically increasing.
        // The polynomial fit is unreliable beyond this point, so we clamp AOT
        // to this value during aerosol retrieval. Matches C code's roatm_iaMax logic.
        let mut ia_max = NAOT_VALS - 1;
        for ia in 1..NAOT_VALS {
            if roatm_vals[ia] - roatm_vals[ia - 1] <= 1.0e-5 {
                ia_max = ia - 1;
                break;
            }
        }
        let roatm_upper = AOT550NM[ia_max];

        coefficients.push(AtmCorrCoefficients {
            roatm_upper,
            roatm_coef,
            ttatmg_coef,
            satm_coef,
        });

        tgo_arr.push(tgo_band);

        // Extract normext at reference pressure index (ip=0) and reference AOT (iaot=3)
        // normext layout: [iband * NPRES_VALS * NAOT_VALS + ip * NAOT_VALS + iaot]
        let ip_ref = 0;
        let iaot_ref = 3;
        let normext_idx = iband * NPRES_VALS * NAOT_VALS + ip_ref * NAOT_VALS + iaot_ref;
        let ne = if normext_idx < lut.normext.len() {
            lut.normext[normext_idx]
        } else {
            1.0
        };
        normext_p0a3.push(ne);
    }

    (coefficients, tgo_arr, normext_p0a3)
}

/// Compute surface reflectance for an entire scene.
///
/// This is the top-level orchestrator that:
/// 1. Maps scene center to CMG grid, extracts atmospheric params
/// 2. Pre-computes polynomial coefficients for each band
/// 3. Runs aerosol retrieval loop over window centers
/// 4. Fixes invalid aerosols, interpolates
/// 5. Applies final per-pixel correction with retrieved aerosol
/// 6. Scales to output integers
#[allow(clippy::too_many_arguments)]
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
    let (nlines, nsamps) = toa_bands[0].dim();
    let nbands = sensor.num_refl_bands();
    let lambda = sensor.lambda();
    let aero_window = sensor.aerosol_window();
    let half_aero_window = sensor.half_aerosol_window();

    // ── Step 1: Scene-center atmospheric state ──
    let center_line_geo = (nlines / 2) as i32;
    let center_samp_geo = (nsamps / 2) as i32;
    let (scene_center_lat, scene_center_lon) =
        utm_to_deg(space_def, center_line_geo, center_samp_geo);
    let (pressure, uoz, uwv) = extract_atm_params(aux, scene_center_lat, scene_center_lon);

    // Scene-center geometry (use center pixel angles)
    let center_line = nlines / 2;
    let center_samp = nsamps / 2;
    let xts_center = solar_zenith[(center_line, center_samp)] as f64;
    let xtv_center = view_zenith[(center_line, center_samp)] as f64;
    let xmus_center = (xts_center * DEG2RAD).cos();
    let xmuv_center = (xtv_center * DEG2RAD).cos();
    let xfi_center = (view_azimuth[(center_line, center_samp)]
        - solar_azimuth[(center_line, center_samp)])
        .abs() as f64;
    let xfi_center = if xfi_center > 180.0 {
        360.0 - xfi_center
    } else {
        xfi_center
    };
    let cosxfi_center = (xfi_center * DEG2RAD).cos();

    // Build gas coefficients per band from sensor-specific constants.
    let gas_coeff: Vec<GasCoefficients> = sensor.gas_coefficients();

    // ── Step 2: Pre-compute polynomial coefficients at scene center ──
    let (atm_coeff, tgo_arr, normext_p0a3) = precompute_coefficients(
        sensor,
        lut,
        &gas_coeff,
        xts_center,
        xtv_center,
        xmus_center,
        xmuv_center,
        xfi_center,
        cosxfi_center,
        pressure,
        uoz,
        uwv,
    );

    // Extract per-band polynomial coefficient arrays for subaeroret_new
    let roatm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.roatm_coef).collect();
    let ttatmg_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.ttatmg_coef).collect();
    let satm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.satm_coef).collect();
    let roatm_ia_max: Vec<f64> = atm_coeff.iter().map(|c| c.roatm_upper).collect();

    // Threshold values for negative SR check (per band)
    let tth = vec![1.0e-3; nbands];

    // ── Step 3: Allocate working arrays ──
    let npix = nlines * nsamps;
    let mut taero = vec![DEFAULT_AERO; npix];
    let mut teps = vec![DEFAULT_EPS; npix];
    let mut ipflag = vec![0u8; npix];

    // Flatten QA band for 1D access
    let qa_flat: Vec<u16> = qa_band.iter().copied().collect();

    // ── Step 4: Aerosol retrieval at window centers ──
    let bi = sensor.band_indices();
    let iband1 = bi.red; // reference band for aerosol retrieval

    // Iterate over window centers
    let mut iline = half_aero_window;
    while iline < nlines {
        let mut isamp = half_aero_window;
        while isamp < nsamps {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels
            if qa_flat[pix] == INPUT_FILL {
                ipflag[pix] = 1u8 << IPFLAG_FILL;
                isamp += aero_window;
                continue;
            }

            // Get TOA reflectance at this pixel for all bands
            let troatm: Vec<f64> = (0..nbands)
                .map(|ib| toa_bands[ib][(iline, isamp)] as f64)
                .collect();

            // Check for valid TOA data (skip if any band is fill)
            let has_valid_data = troatm.iter().all(|&v| v > -0.5 && v < 2.0);
            if !has_valid_data {
                ipflag[pix] = 1u8 << IPFLAG_FAILED;
                isamp += aero_window;
                continue;
            }

            // Compute per-pixel lat/lon from image coordinates
            let (pixel_lat, pixel_lon) = utm_to_deg(space_def, iline as i32, isamp as i32);

            // Look up band ratios and water flag
            let (rb1, rb2, rb7) = lookup_band_ratios(aux, pixel_lat, pixel_lon);
            let is_water = is_water_pixel(aux, pixel_lat, pixel_lon);

            // Build erelc array
            let erelc = build_erelc(sensor, rb1, rb2, rb7, is_water);

            // Determine eps values to try
            let eps_values: &[f64] = if is_water {
                &[WATER_EPS]
            } else {
                &[LOW_EPS, WATER_EPS, MOD_EPS, HIGH_EPS]
            };

            // Try each eps, pick the one with lowest residual
            let mut best_result = None;
            let mut best_residual = f64::MAX;

            for &eps in eps_values {
                let result = subaeroret_new(
                    is_water,
                    iband1,
                    &erelc,
                    &troatm,
                    &tgo_arr,
                    &roatm_ia_max,
                    &roatm_coef,
                    &ttatmg_coef,
                    &satm_coef,
                    &normext_p0a3,
                    lambda,
                    eps,
                    0, // iaots start
                    &tth,
                );

                if result.residual < best_residual {
                    best_residual = result.residual;
                    best_result = Some(result);
                }
            }

            if let Some(result) = best_result {
                taero[pix] = result.raot;
                teps[pix] = result.eps;
                ipflag[pix] = 1u8 << IPFLAG_CLEAR;

                // Set aerosol QA bits based on residual
                if is_water {
                    ipflag[pix] |= 1u8 << IPFLAG_WATER;
                }
                if result.residual < LOW_AERO_THRESH {
                    ipflag[pix] |= (1u8 << AERO1_QA) | (1u8 << AERO2_QA);
                } else if result.residual < AVG_AERO_THRESH {
                    ipflag[pix] |= 1u8 << AERO1_QA;
                }
            } else {
                ipflag[pix] = 1u8 << IPFLAG_FAILED;
            }

            isamp += aero_window;
        }
        iline += aero_window;
    }

    // ── Step 5: Fix invalid aerosols and interpolate ──
    fix_invalid_aerosols(
        &mut taero,
        &mut teps,
        &mut ipflag,
        nlines,
        nsamps,
        aero_window,
        half_aero_window,
        sensor.fix_aerosol_window(),
        sensor.half_fix_aerosol_window(),
        sensor.min_clear_pix(),
    );

    aerosol_interp(
        &mut taero,
        &mut teps,
        &mut ipflag,
        &qa_flat,
        nlines,
        nsamps,
        aero_window,
        half_aero_window,
    );

    // ── Step 6: Final per-pixel atmospheric correction ──
    let mut sr_f32: Vec<Array2<f64>> = (0..nbands)
        .map(|_| Array2::zeros((nlines, nsamps)))
        .collect();

    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels
            if qa_flat[pix] == INPUT_FILL {
                continue;
            }

            let raot = taero[pix];
            let eps = teps[pix];

            for iband in 0..nbands {
                let rotoa = toa_bands[iband][(iline, isamp)] as f64;

                let roslamb = atmcorlamb2_new(
                    &atm_coeff[iband],
                    tgo_arr[iband],
                    iband,
                    raot,
                    normext_p0a3[iband],
                    rotoa,
                    lambda,
                    eps,
                );

                // Clamp to valid range
                let roslamb = roslamb.clamp(MIN_VALID_REFL, MAX_VALID_REFL);
                sr_f32[iband][(iline, isamp)] = roslamb;
            }
        }
    }

    // ── Step 7: Scale to output integers ──
    let sr_bands: Vec<Array2<i16>> = sr_f32
        .iter()
        .map(|band| {
            band.mapv(|v| {
                ((v + BAND_OFFSET_REFL) * MULT_FACTOR_REFL)
                    .round()
                    .clamp(i16::MIN as f64, i16::MAX as f64) as i16
            })
        })
        .collect();

    // Scale brightness temperature bands
    let bt_out: Vec<Array2<u16>> = bt_bands
        .iter()
        .map(|bt| {
            let arr: Array2<u16> = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
                let val = bt[(i, j)] as f64;
                if val < MIN_VALID_TH || val > MAX_VALID_TH {
                    0u16
                } else {
                    ((val - BAND_OFFSET_TH) * MULT_FACTOR_TH)
                        .round()
                        .clamp(0.0, u16::MAX as f64) as u16
                }
            });
            arr
        })
        .collect();

    // Scale aerosol to int16
    let aerosol = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        let pix = i * nsamps + j;
        if qa_flat[pix] == INPUT_FILL {
            AERO_FILL
        } else {
            (taero[pix] * MULT_FACTOR_AERO)
                .round()
                .clamp(i16::MIN as f64, i16::MAX as f64) as i16
        }
    });

    // Build QA output
    let qa = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        let pix = i * nsamps + j;
        ipflag[pix]
    });

    SurfaceReflectanceResult {
        sr_bands,
        bt_bands: bt_out,
        aerosol,
        qa,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_surface_reflectance_result_structure() {
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

    #[test]
    fn test_latlon_to_cmg_equator_prime_meridian() {
        let cmg = latlon_to_cmg(0.0, 0.0);
        // Primary cell (idx[0]) should be near row=1799, col=3599
        // (89.975 - 0) * 20 = 1799.5 -> row 1799
        // (179.975 + 0) * 20 = 3599.5 -> col 3599
        let row = cmg.idx[0] / CMG_NBLON;
        let col = cmg.idx[0] % CMG_NBLON;
        assert_eq!(row, 1799);
        assert_eq!(col, 3599);
        // Weights should sum to 1.0
        let wsum: f64 = cmg.w.iter().sum();
        assert!((wsum - 1.0).abs() < 1e-10);
        // Fractional position ~0.5 in both directions
        assert!(cmg.w[0] > 0.2 && cmg.w[0] < 0.3); // ~0.25
    }

    #[test]
    fn test_latlon_to_cmg_corners() {
        // Top-left: 90N, 180W -> row 0, col 0, no fractional offset
        let cmg = latlon_to_cmg(90.0, -180.0);
        let row = cmg.idx[0] / CMG_NBLON;
        let col = cmg.idx[0] % CMG_NBLON;
        assert_eq!(row, 0);
        assert_eq!(col, 0);
        // At the corner, almost all weight should be on idx[0]
        assert!(cmg.w[0] > 0.9);

        // Bottom-right: 90S, 180E (clamped to grid)
        let cmg = latlon_to_cmg(-90.0, 180.0);
        let row = cmg.idx[0] / CMG_NBLON;
        assert!(row >= CMG_NBLAT - 2);
    }

    #[test]
    fn test_bilerp_uniform() {
        // All same value -> interpolation returns that value
        let w = [0.25, 0.25, 0.25, 0.25];
        assert!((bilerp([5.0, 5.0, 5.0, 5.0], &w) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_bilerp_weights() {
        // All weight on first cell
        let w = [1.0, 0.0, 0.0, 0.0];
        assert!((bilerp([10.0, 20.0, 30.0, 40.0], &w) - 10.0).abs() < 1e-10);
        // All weight on last cell
        let w = [0.0, 0.0, 0.0, 1.0];
        assert!((bilerp([10.0, 20.0, 30.0, 40.0], &w) - 40.0).abs() < 1e-10);
    }

    #[test]
    fn test_pressure_from_elevation() {
        // Sea level
        let p0 = pressure_from_elevation(0.0);
        assert!((p0 - ATMOS_PRES_0).abs() < 0.01);

        // Higher elevation = lower pressure
        let p1 = pressure_from_elevation(1000.0);
        assert!(p1 < p0);
        assert!(p1 > 0.0);
    }

    #[test]
    fn test_auxiliary_data_struct() {
        let aux = AuxiliaryData {
            dem: vec![0; 10],
            wv: vec![0; 10],
            oz: vec![0; 10],
            ratiob1: vec![0; 10],
            ratiob2: vec![0; 10],
            ratiob7: vec![0; 10],
            intratiob1: vec![0; 10],
            intratiob2: vec![0; 10],
            intratiob7: vec![0; 10],
            slpratiob1: vec![0; 10],
            slpratiob2: vec![0; 10],
            slpratiob7: vec![0; 10],
            andwi: vec![0; 10],
            sndwi: vec![0; 10],
            wv_scale: 0.001,
            oz_scale: 0.001,
            wv_default: 0.5,
            oz_default: 0.3,
        };
        assert_eq!(aux.dem.len(), 10);
        assert_eq!(aux.wv_scale, 0.001);
    }
}
