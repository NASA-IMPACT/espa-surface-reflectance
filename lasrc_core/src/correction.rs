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
    /// Surface reflectance per band, scaled uint16 (fill = 0)
    pub sr_bands: Vec<Array2<u16>>,
    /// Brightness temperature (Landsat only), scaled uint16
    pub bt_bands: Vec<Array2<u16>>,
    /// AOT, scaled int16 (fill = -9999, valid range [0, 5000])
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

/// Compute modified slope and intercept for a single CMG corner pixel.
///
/// Matches C code logic at lines 1109-1128: if band ratios are out of range
/// (water-like), use fixed intercepts; if sndwi < 200 (land), use ratio as
/// intercept; otherwise keep existing slope/intercept.
///
/// Returns (slope, intercept) already unscaled by 0.001.
fn modified_slope_intercept(
    ratiob1_val: i16,
    ratiob2_val: i16,
    sndwi_val: i16,
    slpratio_val: i16,
    intratio_val: i16,
    ratiob_val: i16,
    default_water_intercept: i16,
) -> (f64, f64) {
    let rb1 = ratiob1_val as f64 * 0.001;
    let rb2 = ratiob2_val as f64 * 0.001;

    if rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1 {
        // Water-like: fixed intercepts, zero slopes
        (0.0, default_water_intercept as f64 * 0.001)
    } else if sndwi_val < 200 {
        // Land: slope=0, intercept=ratio value
        (0.0, ratiob_val as f64 * 0.001)
    } else {
        // Keep existing slope/intercept values from the file
        (slpratio_val as f64 * 0.001, intratio_val as f64 * 0.001)
    }
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
) -> (Vec<AtmCorrCoefficients>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let nbands = sensor.num_refl_bands();
    let tauray = sensor.tauray();
    let lambda = sensor.lambda();
    // max_band_idx: last band with valid wavelength for Angstrom scaling.
    // For Landsat: DNL_BAND7 = 6 (bands 0-6 have wavelengths), band 7 (SRL_BAND9) does not.
    // For Sentinel: DNS_BAND12 = last index.
    let max_band_idx = lambda.len() - 1;

    let mut coefficients = Vec::with_capacity(nbands);
    let mut tgo_arr = Vec::with_capacity(nbands);
    let mut normext_p0a3 = Vec::with_capacity(nbands);
    // Climatological atmospheric params at fixed AOT=0.05 (aot550nm[1]),
    // used to do a simplified first-pass SR and then reconstruct TOA.
    let mut btgo = Vec::with_capacity(nbands);
    let mut broatm = Vec::with_capacity(nbands);
    let mut bttatmg = Vec::with_capacity(nbands);
    let mut bsatm = Vec::with_capacity(nbands);

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
                lambda,
                max_band_idx,
                HIGH_EPS, // C uses eps=2.5 for precomputation
            );
            roatm_vals[iaot] = result.roatm;
            ttatmg_vals[iaot] = result.ttatmg;
            satm_vals[iaot] = result.satm;
            if iaot == 0 {
                tgo_band = result.tgo;
            }
            // Save climatological params at iaot=1 (AOT=0.05)
            if iaot == 1 {
                btgo.push(result.tgo);
                broatm.push(result.roatm);
                bttatmg.push(result.ttatmg);
                bsatm.push(result.satm);
            }
        }

        // Find the last AOT index where roatm is still monotonically increasing.
        // C: iaMaxTemp starts at 1, loops ia=1..NAOT_VALS-1, sets NAOT_VALS-1
        // at the last iteration, breaks when diff <= ESPA_EPSILON.
        let mut ia_max = 1usize;
        for ia in 1..NAOT_VALS {
            if ia == NAOT_VALS - 1 {
                ia_max = NAOT_VALS - 1;
            }
            if roatm_vals[ia] - roatm_vals[ia - 1] > 1.0e-5 {
                continue;
            } else {
                ia_max = ia - 1;
                break;
            }
        }
        let roatm_upper = AOT550NM[ia_max];

        // Fit cubic polynomials.
        // C fits roatm to only the first ia_max points (where monotonically
        // increasing), but ttatmg and satm use all NAOT_VALS points.
        let roatm_coef = get_3rd_order_poly_coeff(&AOT550NM[..ia_max], &roatm_vals[..ia_max]);
        let ttatmg_coef = get_3rd_order_poly_coeff(&AOT550NM, &ttatmg_vals);
        let satm_coef = get_3rd_order_poly_coeff(&AOT550NM, &satm_vals);

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

    (coefficients, tgo_arr, normext_p0a3, btgo, broatm, bttatmg, bsatm)
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
    let (atm_coeff, tgo_arr, normext_p0a3, btgo, broatm, bttatmg, bsatm) = precompute_coefficients(
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

    // ── Step 3: Allocate working arrays ──
    let npix = nlines * nsamps;
    let mut taero = vec![DEFAULT_AERO; npix];
    let mut teps = vec![DEFAULT_EPS; npix];
    let mut ipflag = vec![0u8; npix];

    // Flatten QA band for 1D access
    let qa_flat: Vec<u16> = qa_band.iter().copied().collect();

    let bi = sensor.band_indices();
    let iband1 = bi.red; // reference band for aerosol retrieval

    // Threshold values for negative SR check (per band)
    let tth = &LANDSAT_TTH[..nbands];
    let tth_water = &LANDSAT_TTH_WATER[..nbands];

    // ── Step 4: Climatological per-pixel atmospheric correction ──
    // Simplified first-pass SR using scene-center atmospheric params at
    // fixed AOT=0.05. This matches C code lines 648-691: for each pixel,
    // roslamb = (sband[ib][i] - tgo*roatm) / (tgo*ttatmg + satm*(sband[ib][i] - tgo*roatm))
    // where sband[ib][i] is TOA/cos(sza).
    //
    // This MUST be computed before the aerosol retrieval because the
    // retrieval uses climatological SR bands 5 (NIR) and 7 (SWIR2) for NDWI.
    let mut sband: Vec<Vec<f64>> = (0..nbands)
        .map(|_| vec![0.0f64; npix])
        .collect();

    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;

            if is_fill_pixel(qa_flat[pix]) {
                if iline == 0 || isamp == 0 {
                    // Initialize fill flag (C does this for band 1 only)
                }
                ipflag[pix] = 1u8 << IPFLAG_FILL;
                continue;
            }

            let xmus = (solar_zenith[(iline, isamp)] as f64 * DEG2RAD).cos();

            for iband in 0..nbands {
                // TOA / cos(SZA), clamped to valid range
                let raw_toa = toa_bands[iband][(iline, isamp)] as f64;
                let rotoa = (raw_toa / xmus).clamp(MIN_VALID_REFL, MAX_VALID_REFL);

                // Simplified atmospheric correction using scene-center params
                let tgo_x_roatm = btgo[iband] * broatm[iband];
                let tgo_x_ttatmg = btgo[iband] * bttatmg[iband];
                let roslamb = {
                    let num = rotoa - tgo_x_roatm;
                    num / (tgo_x_ttatmg + bsatm[iband] * num)
                };
                sband[iband][pix] = roslamb.clamp(MIN_VALID_REFL, MAX_VALID_REFL);
            }
        }
    }

    // ── Step 5: Aerosol retrieval at window centers ──
    // Precompute the quadratic fit constants for eps optimization.
    // C code lines 963-969.
    let eps1 = LOW_EPS;
    let eps2 = MOD_EPS;
    let eps3 = HIGH_EPS;
    let xa = eps1 * eps1 - eps3 * eps3; // -5.25
    let xd = ((eps2 * eps2 - eps3 * eps3) as i32) as f64; // integer truncation: -3.0
    let xb = eps1 - eps3; // -1.5
    let xe = eps2 - eps3; // -0.75

    // Iterate over window centers
    let mut iline = half_aero_window;
    while iline < nlines {
        let mut isamp = half_aero_window;
        while isamp < nsamps {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels (keep current fill-skip behavior)
            if is_fill_pixel(qa_flat[pix]) {
                ipflag[pix] = 1u8 << IPFLAG_FILL;
                isamp += aero_window;
                continue;
            }

            // Get TOA reflectance at this pixel for the needed bands,
            // divided by cos(solar zenith) to match C code's TOA normalization.
            // These are the original TOA values (aerob1, aerob2, aerob4, aerob5, aerob7
            // in the C code).
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

            // For each of the 4 CMG corners, compute modified slope/intercept
            // for bands 1, 2, 7. The C code mutates the aux arrays, but we
            // compute the modified values locally.
            let corner_indices = [ratio_pix11, ratio_pix12, ratio_pix21, ratio_pix22];

            // Helper to safely read aux array value
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

                // Band 1
                let (s, i) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob1, cidx),
                    safe_read(&aux.intratiob1, cidx),
                    rb1_val,
                    550, // default water intercept for band 1
                );
                slp_b1[ci] = s;
                int_b1[ci] = i;

                // Band 2
                let (s, i) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob2, cidx),
                    safe_read(&aux.intratiob2, cidx),
                    rb2_val,
                    600, // default water intercept for band 2
                );
                slp_b2[ci] = s;
                int_b2[ci] = i;

                // Band 7
                let (s, i) = modified_slope_intercept(
                    rb1_val, rb2_val, sndwi_val,
                    safe_read(&aux.slpratiob7, cidx),
                    safe_read(&aux.intratiob7, cidx),
                    safe_read(&aux.ratiob7, cidx),
                    2000, // default water intercept for band 7
                );
                slp_b7[ci] = s;
                int_b7[ci] = i;
            }

            // Bilinearly interpolate slopes and intercepts
            let slprb1 = bilerp(slp_b1, &cmg.w);
            let intrb1 = bilerp(int_b1, &cmg.w);
            let slprb2 = bilerp(slp_b2, &cmg.w);
            let intrb2 = bilerp(int_b2, &cmg.w);
            let slprb7 = bilerp(slp_b7, &cmg.w);
            let intrb7 = bilerp(int_b7, &cmg.w);

            // Compute NDWI from climatological SR bands 5 (NIR) and 7 (SWIR2)
            let sr_nir = sband[bi.nir][pix];
            let sr_swir2 = sband[bi.swir2][pix];
            let sr_swir2_half = sr_swir2 * 0.5;
            let denom = sr_nir + sr_swir2_half;
            let mut xndwi = if denom.abs() > 1.0e-10 {
                (sr_nir - sr_swir2_half) / denom
            } else {
                0.0
            };

            // Clamp NDWI using andwi/sndwi thresholds from CMG (uses ratio_pix11)
            let andwi_val = safe_read(&aux.andwi, ratio_pix11);
            let sndwi_val = safe_read(&aux.sndwi, ratio_pix11);
            let ndwi_th1 = (andwi_val as f64 + 2.0 * sndwi_val as f64) * 0.001;
            let ndwi_th2 = (andwi_val as f64 - 2.0 * sndwi_val as f64) * 0.001;
            if xndwi > ndwi_th1 {
                xndwi = ndwi_th1;
            }
            if xndwi < ndwi_th2 {
                xndwi = ndwi_th2;
            }

            // Initialize erelc and troatm arrays
            let mut erelc = vec![-1.0f64; nbands];
            let mut troatm = vec![0.0f64; nbands];

            // Compute band ratios from NDWI, slopes, and intercepts
            erelc[bi.coastal] = xndwi * slprb1 + intrb1;
            erelc[bi.blue] = xndwi * slprb2 + intrb2;
            erelc[bi.red] = 1.0;
            erelc[bi.swir2] = xndwi * slprb7 + intrb7;

            // Set TOA reflectance values for the needed bands
            troatm[bi.coastal] = toa_over_cos(bi.coastal);
            troatm[bi.blue] = toa_over_cos(bi.blue);
            troatm[bi.red] = toa_over_cos(bi.red);
            troatm[bi.swir2] = toa_over_cos(bi.swir2);

            // === Eps optimization: 3 retrievals at eps1=1.0, eps2=1.75, eps3=2.5 ===
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

            // Quadratic fit for optimal eps
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
                // epsmin >= HIGH_EPS
                (eps3, sraot3, residual3)
            };

            teps[pix] = eps;
            taero[pix] = raot;

            // corf = raot / xmus_center for !use_orig_aero
            let corf = raot / xmus_center;

            // === Post-retrieval validation ===
            if residual < (0.015 + 0.005 * corf + 0.10 * troatm[bi.swir2]) {
                // Check NIR (band 5) and red (band 4) to compute NDVI
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
                    ipflag[pix] |= 1u8 << IPFLAG_CLEAR;
                } else {
                    ipflag[pix] |= 1u8 << IPFLAG_WATER;
                }
            } else {
                ipflag[pix] |= 1u8 << IPFLAG_WATER;
            }

            // === Water retest ===
            if ipflag[pix] & (1u8 << IPFLAG_WATER) != 0 {
                // Water band ratios: all active bands set to 1.0
                let mut water_erelc = vec![-1.0f64; nbands];
                water_erelc[bi.coastal] = 1.0;
                water_erelc[bi.red] = 1.0;
                water_erelc[bi.nir] = 1.0;
                water_erelc[bi.swir2] = 1.0;

                // Water TOA values
                let mut water_troatm = vec![0.0f64; nbands];
                water_troatm[bi.coastal] = toa_over_cos(bi.coastal);
                water_troatm[bi.red] = toa_over_cos(bi.red);
                water_troatm[bi.nir] = toa_over_cos(bi.nir);
                water_troatm[bi.swir2] = toa_over_cos(bi.swir2);

                let water_result = subaeroret_new(
                    true, iband1, &water_erelc, &water_troatm,
                    &tgo_arr, &roatm_ia_max, &roatm_coef, &ttatmg_coef,
                    &satm_coef, &normext_p0a3, lambda, WATER_EPS, 0, tth_water,
                );

                teps[pix] = WATER_EPS;
                taero[pix] = water_result.raot;
                let water_corf = water_result.raot / xmus_center;

                // Validate: check band 1 reflectance
                let ros1 = atmcorlamb2_new(
                    &atm_coeff[bi.coastal], tgo_arr[bi.coastal], bi.coastal,
                    water_result.raot, normext_p0a3[bi.coastal],
                    toa_over_cos(bi.coastal), lambda, WATER_EPS,
                );

                if water_result.residual > (0.010 + 0.005 * water_corf) || ros1 < 0.0 {
                    // Not valid water, clear all QA bits
                    ipflag[pix] = 0;
                } else {
                    // Valid water pixel
                    ipflag[pix] = (1u8 << IPFLAG_CLEAR) | (1u8 << IPFLAG_WATER);
                }
            }

            isamp += aero_window;
        }
        iline += aero_window;
    }

    // ── Step 6: Fix invalid aerosols and interpolate ──
    fix_invalid_aerosols(
        &mut taero,
        &mut teps,
        &ipflag,
        nlines,
        nsamps,
        aero_window,
        half_aero_window,
        sensor.fix_aerosol_window(),
        sensor.half_fix_aerosol_window(),
        sensor.min_clear_pix(),
    );

    // Interpolate taero, then teps (matching C: two separate calls)
    aerosol_interp(
        &mut taero,
        &mut ipflag,
        &qa_flat,
        nlines,
        nsamps,
        aero_window,
        half_aero_window,
    );
    aerosol_interp(
        &mut teps,
        &mut ipflag,
        &qa_flat,
        nlines,
        nsamps,
        aero_window,
        half_aero_window,
    );

    // ── Step 7: Final per-pixel atmospheric correction ──
    // For each pixel, reconstruct TOA from climatological SR, then apply
    // atmcorlamb2_new with per-pixel retrieved aerosol.
    // C code: rotoa = (rsurf * bttatmg[ib] / (1 - bsatm[ib] * rsurf) + broatm[ib]) * btgo[ib]
    let mut sr_f32: Vec<Array2<f64>> = (0..nbands)
        .map(|_| Array2::zeros((nlines, nsamps)))
        .collect();

    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels
            if is_fill_pixel(qa_flat[pix]) {
                continue;
            }

            let raot = taero[pix];
            let eps = teps[pix];

            for iband in 0..nbands {
                // Reconstruct TOA from climatological SR
                let rsurf = sband[iband][pix];
                let rotoa = (rsurf * bttatmg[iband]
                    / (1.0 - bsatm[iband] * rsurf)
                    + broatm[iband])
                    * btgo[iband];

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

                // Set aerosol QA bits on the coastal aerosol band (band 0)
                // using |rsurf - roslamb| as the aerosol level indicator.
                // Matches C code compute_landsat_refl.c lines 1758-1780.
                if iband == bi.coastal {
                    let tmpf = (rsurf - roslamb).abs();
                    if tmpf <= LOW_AERO_THRESH {
                        // Low aerosol: set AERO1 only
                        ipflag[pix] |= 1u8 << AERO1_QA;
                    } else if tmpf < AVG_AERO_THRESH {
                        // Average aerosol: set AERO2 only
                        ipflag[pix] |= 1u8 << AERO2_QA;
                    } else {
                        // High aerosol: set both AERO1 and AERO2
                        ipflag[pix] |= (1u8 << AERO1_QA) | (1u8 << AERO2_QA);
                    }
                }
            }
        }
    }

    // ── Step 8: Scale to output integers ──
    let sr_bands: Vec<Array2<u16>> = sr_f32
        .iter()
        .map(|band| {
            Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
                let pix = i * nsamps + j;
                if is_fill_pixel(qa_flat[pix]) {
                    0u16
                } else {
                    ((band[(i, j)] + BAND_OFFSET_REFL) * MULT_FACTOR_REFL)
                        .round()
                        .clamp(0.0, u16::MAX as f64) as u16
                }
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

    // Scale aerosol to int16 (valid range [0, 5000], fill = -9999)
    let aerosol = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        let pix = i * nsamps + j;
        if is_fill_pixel(qa_flat[pix]) {
            AERO_FILL
        } else {
            (taero[pix] * MULT_FACTOR_AERO)
                .round()
                .clamp(0.0, 5000.0) as i16
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
