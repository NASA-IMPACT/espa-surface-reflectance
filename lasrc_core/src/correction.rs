//! Top-level correction orchestrator that ties together all compute modules
//! into the main `compute_surface_reflectance` function.
//!
//! Ported from the C LaSRC main processing loop.

use ndarray::{Array2, ArrayView2};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::aerosol::{
    aerosol_interp, fix_invalid_aerosols, subaeroret_new,
    aerosol_interp_sentinel, ipflag_expand_failed_sentinel, aero_avg_failed_sentinel,
};
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

/// Result from a single aerosol window-center retrieval, used for parallel collect-and-scatter.
struct AerosolWindowResult {
    pix: usize,
    taero: f32,
    teps: f32,
    ipflag: u8,
}

struct SentinelAerosolWindowResult {
    win_i: usize,
    win_j: usize,
    curr_pix: usize,
    taero: f32,
    teps: f32,
    ipflag: u8,
    is_fill: bool,
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
    // C declares `float xcmg, ycmg` — truncate to f32 before computing
    // integer indices and fractional offsets so that cell selection and
    // bilinear weights match C exactly.
    let ycmg = ((89.975 - lat) * 20.0) as f32;
    let xcmg = ((179.975 + lon) * 20.0) as f32;

    // Integer grid indices (truncation matches C code's (int) cast on float)
    let row = (ycmg as isize).clamp(0, (CMG_NBLAT - 1) as isize) as usize;
    let col = (xcmg as isize).clamp(0, (CMG_NBLON - 1) as isize) as usize;

    // Next row/col with wrapping at edges (matching C code)
    let row1 = if row >= CMG_NBLAT - 1 { 0 } else { row + 1 };
    let col1 = if col >= CMG_NBLON - 1 { 0 } else { col + 1 };

    // C: float u = (ycmg - lcmg), float v = (xcmg - scmg)
    let u = (ycmg - row as f32).clamp(0.0, 1.0);
    let v = (xcmg - col as f32).clamp(0.0, 1.0);

    // C: float weights for bilinear interpolation
    let one_minus_u = 1.0f32 - u;
    let one_minus_v = 1.0f32 - v;

    CmgPosition {
        idx: [
            row * CMG_NBLON + col,
            row * CMG_NBLON + col1,
            row1 * CMG_NBLON + col,
            row1 * CMG_NBLON + col1,
        ],
        w: [
            (one_minus_u * one_minus_v) as f64,
            (one_minus_u * v) as f64,
            (u * one_minus_v) as f64,
            (u * v) as f64,
        ],
    }
}

/// Bilinear interpolation over 4 values with the given weights.
fn bilerp(vals: [f64; 4], w: &[f64; 4]) -> f64 {
    vals[0] * w[0] + vals[1] * w[1] + vals[2] * w[2] + vals[3] * w[3]
}

/// Bilinear interpolation in f32 precision, matching C's `float` arithmetic.
/// The CMG weights are already f32-precision (stored in f64); values are truncated
/// to f32 before interpolation.
fn bilerp_f32(vals: [f64; 4], w: &[f64; 4]) -> f32 {
    let v0 = vals[0] as f32;
    let v1 = vals[1] as f32;
    let v2 = vals[2] as f32;
    let v3 = vals[3] as f32;
    let w0 = w[0] as f32;
    let w1 = w[1] as f32;
    let w2 = w[2] as f32;
    let w3 = w[3] as f32;
    v0 * w0 + v1 * w1 + v2 * w2 + v3 * w3
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

/// Scene-center atmospheric parameter extraction matching C's init_sr_refl.
///
/// Unlike `extract_atm_params` (which uses truncation + bilinear interpolation),
/// this uses `roundf` to nearest CMG pixel and single-pixel lookup, matching the
/// C code's behavior in init_sr_refl (lines 224-258).
fn extract_atm_params_scene_center(aux: &AuxiliaryData, lat: f64, lon: f64) -> (f64, f64, f64) {
    let ycmg = (89.975 - lat) * 20.0;
    let xcmg = (179.975 + lon) * 20.0;

    // C: lcmg = (int) roundf(ycmg)
    let lcmg = (ycmg as f32).round() as isize;
    let scmg = (xcmg as f32).round() as isize;

    let lcmg = lcmg.clamp(0, CMG_NBLAT as isize) as usize;
    let scmg = scmg.clamp(0, CMG_NBLON as isize) as usize;

    let cmg_pix = lcmg * CMG_NBLON + scmg;

    // Water vapor: single pixel lookup
    let uwv = if cmg_pix < aux.wv.len() && aux.wv[cmg_pix] > 0 {
        aux.wv[cmg_pix] as f64 / aux.wv_scale
    } else {
        aux.wv_default
    };

    // Ozone: single pixel lookup
    let uoz = if cmg_pix < aux.oz.len() && aux.oz[cmg_pix] > 0 {
        aux.oz[cmg_pix] as f64 / aux.oz_scale
    } else {
        aux.oz_default
    };

    // Pressure from DEM: single pixel lookup
    // C: dem_pix = lcmg * DEM_NBLON + scmg (same grid as CMG)
    let pressure = if cmg_pix < aux.dem.len() && aux.dem[cmg_pix] != -9999 {
        pressure_from_elevation(aux.dem[cmg_pix] as f64)
    } else {
        ATMOS_PRES_0
    };

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
    precomp_eps: f64,
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
                precomp_eps, // Landsat=2.5, Sentinel=-1.0
            );
            roatm_vals[iaot] = result.roatm;
            ttatmg_vals[iaot] = result.ttatmg;
            satm_vals[iaot] = result.satm;
            if iaot == 0 {
                tgo_band = result.tgo;
            }
            // Save climatological params at iaot=1 (AOT=0.05)
            // C: btgo/broatm/bttatmg/bsatm are float[] — truncate to f32
            if iaot == 1 {
                btgo.push(result.tgo as f32 as f64);
                broatm.push(result.roatm as f32 as f64);
                bttatmg.push(result.ttatmg as f32 as f64);
                bsatm.push(result.satm as f32 as f64);
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
            // C: float subtraction compared against ESPA_EPSILON (double).
            // roatm_arr values are float, so the subtraction is float precision.
            let diff_f32 = roatm_vals[ia] as f32 - roatm_vals[ia - 1] as f32;
            if (diff_f32 as f64) > 1.0e-5 {
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

        // C: float tgo_arr[] — truncate to f32 precision
        tgo_arr.push(tgo_band as f32 as f64);

        // Extract normext at reference pressure index (ip=0) and reference AOT (iaot=3)
        // normext layout: [iband * NPRES_VALS * NAOT_VALS + ip * NAOT_VALS + iaot]
        // C: float normext_p0a3_arr[] — truncate to f32 precision
        let ip_ref = 0;
        let iaot_ref = 3;
        let normext_idx = iband * NPRES_VALS * NAOT_VALS + ip_ref * NAOT_VALS + iaot_ref;
        let ne = if normext_idx < lut.normext.len() {
            lut.normext[normext_idx] as f32 as f64
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
    num_threads: Option<usize>,
) -> SurfaceReflectanceResult {
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0))
        .build()
        .expect("Failed to build Rayon thread pool");
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
        HIGH_EPS, // Landsat uses eps=2.5 for precomputation
    );

    // Extract per-band polynomial coefficient arrays for subaeroret_new
    let roatm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.roatm_coef).collect();
    let ttatmg_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.ttatmg_coef).collect();
    let satm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.satm_coef).collect();
    let roatm_ia_max: Vec<f64> = atm_coeff.iter().map(|c| c.roatm_upper).collect();

    // ── Step 3: Allocate working arrays ──
    // C: float *taero, float *teps — use f32 to match
    let npix = nlines * nsamps;
    let mut taero = vec![DEFAULT_AERO as f32; npix];
    let mut teps = vec![DEFAULT_EPS as f32; npix];
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
    // C: float **sband — use f32 to match
    let mut sband: Vec<Vec<f32>> = (0..nbands)
        .map(|_| vec![0.0f32; npix])
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
                // C: sband[ib][i] = toa / cos(sza) — float
                let raw_toa = toa_bands[iband][(iline, isamp)] as f64;
                let rotoa = (raw_toa / xmus).clamp(MIN_VALID_REFL, MAX_VALID_REFL) as f32;

                // C: float arithmetic throughout climatological correction
                // tgo_x_roatm = tgo * roatm (float * float = float)
                let tgo_x_roatm = btgo[iband] as f32 * broatm[iband] as f32;
                let tgo_x_ttatmg = btgo[iband] as f32 * bttatmg[iband] as f32;
                let roslamb: f32 = {
                    let num = rotoa - tgo_x_roatm;
                    num / (tgo_x_ttatmg + bsatm[iband] as f32 * num)
                };
                sband[iband][pix] = roslamb.clamp(MIN_VALID_REFL as f32, MAX_VALID_REFL as f32);
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

    // Pre-collect all window center coordinates
    let window_centers: Vec<(usize, usize)> = {
        let mut coords = Vec::new();
        let mut iline = half_aero_window;
        while iline < nlines {
            let mut isamp = half_aero_window;
            while isamp < nsamps {
                coords.push((iline, isamp));
                isamp += aero_window;
            }
            iline += aero_window;
        }
        coords
    };

    // Parallel aerosol retrieval at window centers
    let aero_results: Vec<AerosolWindowResult> = pool.install(|| {
        window_centers.par_iter().map(|&(iline, isamp)| {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels
            if is_fill_pixel(qa_flat[pix]) {
                return AerosolWindowResult {
                    pix,
                    taero: 0.0,
                    teps: 0.0,
                    ipflag: 1u8 << IPFLAG_FILL,
                };
            }

            // Get TOA reflectance at this pixel for the needed bands,
            // divided by cos(solar zenith) to match C code's TOA normalization.
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
            // for bands 1, 2, 7.
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
            let sr_nir = sband[bi.nir][pix] as f64;
            let sr_swir2 = sband[bi.swir2][pix] as f64;
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
            // C: float erelc[NSR_BANDS], float troatm[NSR_BANDS]
            // Truncate to f32 to match C's float precision
            let mut erelc = vec![-1.0f64; nbands];
            let mut troatm = vec![0.0f64; nbands];

            // Compute band ratios from NDWI, slopes, and intercepts
            // C stores these in float arrays, so truncate to f32 precision
            erelc[bi.coastal] = (xndwi * slprb1 + intrb1) as f32 as f64;
            erelc[bi.blue] = (xndwi * slprb2 + intrb2) as f32 as f64;
            erelc[bi.red] = 1.0;
            erelc[bi.swir2] = (xndwi * slprb7 + intrb7) as f32 as f64;

            // Set TOA reflectance values for the needed bands
            // C: troatm[] is float, so truncate
            troatm[bi.coastal] = toa_over_cos(bi.coastal) as f32 as f64;
            troatm[bi.blue] = toa_over_cos(bi.blue) as f32 as f64;
            troatm[bi.red] = toa_over_cos(bi.red) as f32 as f64;
            troatm[bi.swir2] = toa_over_cos(bi.swir2) as f32 as f64;

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

            let mut result_taero = raot as f32;
            let mut result_teps = eps as f32;

            // corf = raot / xmus_center for !use_orig_aero
            let corf = raot / xmus_center;

            // === Post-retrieval validation ===
            let mut result_ipflag: u8 = 0;
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
                    result_ipflag |= 1u8 << IPFLAG_CLEAR;
                } else {
                    result_ipflag |= 1u8 << IPFLAG_WATER;
                }
            } else {
                result_ipflag |= 1u8 << IPFLAG_WATER;
            }

            // === Water retest ===
            if result_ipflag & (1u8 << IPFLAG_WATER) != 0 {
                // Water band ratios: all active bands set to 1.0
                let mut water_erelc = vec![-1.0f64; nbands];
                water_erelc[bi.coastal] = 1.0;
                water_erelc[bi.red] = 1.0;
                water_erelc[bi.nir] = 1.0;
                water_erelc[bi.swir2] = 1.0;

                // Water TOA values (C: float troatm[])
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

                result_teps = WATER_EPS as f32;
                result_taero = water_result.raot as f32;
                let water_corf = water_result.raot / xmus_center;

                // Validate: check band 1 reflectance
                let ros1 = atmcorlamb2_new(
                    &atm_coeff[bi.coastal], tgo_arr[bi.coastal], bi.coastal,
                    water_result.raot, normext_p0a3[bi.coastal],
                    toa_over_cos(bi.coastal), lambda, WATER_EPS,
                );

                if water_result.residual > (0.010 + 0.005 * water_corf) || ros1 < 0.0 {
                    // Not valid water, clear all QA bits
                    result_ipflag = 0;
                } else {
                    // Valid water pixel
                    result_ipflag = (1u8 << IPFLAG_CLEAR) | (1u8 << IPFLAG_WATER);
                }
            }

            AerosolWindowResult {
                pix,
                taero: result_taero,
                teps: result_teps,
                ipflag: result_ipflag,
            }
        }).collect()
    });

    // Scatter results back to the output arrays
    for r in &aero_results {
        taero[r.pix] = r.taero;
        teps[r.pix] = r.teps;
        ipflag[r.pix] = r.ipflag;
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
    let sr_flat: Vec<Vec<f64>> = (0..nbands)
        .map(|_| vec![0.0f64; npix])
        .collect();

    pool.install(|| {
        (0..nlines).into_par_iter().for_each(|iline| {
            for isamp in 0..nsamps {
                let pix = iline * nsamps + isamp;

                // Skip fill pixels
                if is_fill_pixel(qa_flat[pix]) {
                    continue;
                }

                let raot = taero[pix] as f64;
                let eps = teps[pix] as f64;

                for iband in 0..nbands {
                    // Reconstruct TOA from climatological SR
                    let rsurf = sband[iband][pix]; // f32
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

                    // SAFETY: Each iline is processed by exactly one thread
                    // (Rayon's into_par_iter guarantees this). Within a thread,
                    // pix = iline * nsamps + isamp produces unique indices for
                    // that row, so no two threads write to the same index.
                    unsafe {
                        let ptr = sr_flat[iband].as_ptr() as *mut f64;
                        *ptr.add(pix) = roslamb;
                    }
                }
            }
        });
    });

    // Serial post-pass: set aerosol QA bits on the coastal aerosol band
    // using |rsurf - roslamb| as the aerosol level indicator.
    // Matches C code compute_landsat_refl.c lines 1758-1780.
    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;
            if is_fill_pixel(qa_flat[pix]) {
                continue;
            }

            let rsurf = sband[bi.coastal][pix] as f64;
            let roslamb = sr_flat[bi.coastal][pix];
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

    // ── Step 8: Scale to output integers ──
    // C uses float (f32) arithmetic for output scaling:
    //   float tmpf = (sband[band][pix] + offset_value) * mult_value;
    //   out_band[pix] = roundf(tmpf);
    let offset_f32 = BAND_OFFSET_REFL as f32;
    let mult_f32 = MULT_FACTOR_REFL as f32;
    let sr_bands: Vec<Array2<u16>> = sr_flat
        .iter()
        .map(|band| {
            Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
                let pix = i * nsamps + j;
                if is_fill_pixel(qa_flat[pix]) {
                    0u16
                } else {
                    let sband_f32 = band[pix] as f32;
                    let tmpf = (sband_f32 + offset_f32) * mult_f32;
                    tmpf.round().clamp(0.0, u16::MAX as f32) as u16
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
    // C uses float arithmetic: tmpf = (aero[pix] + 0.0f) * 1000.0f; roundf(tmpf)
    let aero_mult_f32 = MULT_FACTOR_AERO as f32;
    let aerosol = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        let pix = i * nsamps + j;
        if is_fill_pixel(qa_flat[pix]) {
            AERO_FILL
        } else {
            let tmpf = taero[pix] * aero_mult_f32;
            tmpf.round().clamp(0.0, 5000.0) as i16
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

/// Compute surface reflectance for a Sentinel-2 scene.
///
/// Similar to `compute_surface_reflectance` (Landsat) but with key differences:
/// - Scene-center scalar angles (not per-pixel grids)
/// - Fill detection: ANY band with TOA == 0.0 marks pixel as fill
/// - TOA is NOT divided by cos(SZA)
/// - Aerosol window starts at (0,0) stepping by SAERO_WINDOW=6
/// - TOA averaging within 6x6 windows for aerosol retrieval
/// - Sentinel-specific eps optimization with resepsmin validation
/// - B09: atmospheric bypass; B10: copy TOA directly
/// - Sentinel-specific post-processing functions
#[allow(clippy::too_many_arguments)]
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
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0))
        .build()
        .expect("Failed to build Rayon thread pool");
    let (nlines, nsamps) = toa_bands[0].dim();
    let nbands = sensor.num_refl_bands(); // 13
    let lambda = sensor.lambda();
    let aero_window = sensor.aerosol_window(); // 6
    let npix = nlines * nsamps;

    // ── Step 1: Fill detection from input QA band ──
    // Bit 0 = fill (set by reader based on DN==0 for any band)
    let qaband: Vec<u16> = (0..npix)
        .map(|pix| {
            let (i, j) = (pix / nsamps, pix % nsamps);
            qa_band[(i, j)]
        })
        .collect();

    // ── Step 2: Scene-center geometry (scalar angles) ──
    let xts = solar_zenith;
    let xtv = view_zenith;
    let xmus = (xts * DEG2RAD).cos();
    let xmuv = (xtv * DEG2RAD).cos();
    let xfi = (view_azimuth - solar_azimuth).abs();
    let xfi = if xfi > 180.0 { 360.0 - xfi } else { xfi };
    let cosxfi = (xfi * DEG2RAD).cos();

    // ── Step 3: Scene-center atmospheric state ──
    // Use rounded nearest-neighbor CMG lookup matching C's init_sr_refl
    let center_line_geo = (nlines / 2) as i32;
    let center_samp_geo = (nsamps / 2) as i32;
    let (scene_center_lat, scene_center_lon) =
        utm_to_deg(space_def, center_line_geo, center_samp_geo);
    let (pressure, uoz, uwv) = extract_atm_params_scene_center(aux, scene_center_lat, scene_center_lon);

    // Build gas coefficients per band from sensor-specific constants.
    let gas_coeff: Vec<GasCoefficients> = sensor.gas_coefficients();

    // ── Step 4: Pre-compute polynomial coefficients ──
    // C uses eps=-1.0 for Sentinel precomputation (no wavelength scaling)
    let (atm_coeff, tgo_arr, normext_p0a3, _btgo, _broatm, _bttatmg, _bsatm) =
        precompute_coefficients(
            sensor, lut, &gas_coeff,
            xts, xtv, xmus, xmuv, xfi, cosxfi,
            pressure, uoz, uwv,
            -1.0, // Sentinel uses eps=-1.0 for precomputation
        );

    // Extract per-band polynomial coefficient arrays for subaeroret_new
    let roatm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.roatm_coef).collect();
    let ttatmg_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.ttatmg_coef).collect();
    let satm_coef: Vec<[f64; NCOEF]> = atm_coeff.iter().map(|c| c.satm_coef).collect();
    let roatm_ia_max: Vec<f64> = atm_coeff.iter().map(|c| c.roatm_upper).collect();

    // ── Step 5: Allocate working arrays ──
    let mut taero = vec![DEFAULT_AERO as f32; npix];
    let mut teps = vec![DEFAULT_EPS as f32; npix];
    let mut ipflag = vec![0u8; npix];

    let bi = sensor.band_indices();
    let iband1 = bi.red; // reference band for aerosol retrieval (DNS_BAND4)

    // Threshold values for aerosol retrieval
    let tth = &SENTINEL_TTH[..nbands];

    // ── Step 6: Climatological per-pixel atmospheric correction ──
    // Simplified first-pass SR using scene-center atmospheric params at
    // fixed AOT=0.05. Sentinel does NOT divide TOA by cos(SZA).
    let mut sband: Vec<Vec<f32>> = (0..nbands)
        .map(|_| vec![0.0f32; npix])
        .collect();

    let tauray = sensor.tauray();
    let max_band_idx = lambda.len() - 1;

    for iband in 0..nbands {
        // C calls atmcorlamb2 with raot=0.05, eps=-1.0 for each band.
        // eps=-1.0 means NO wavelength scaling (mraot550nm = raot directly).
        let (tgo_x_roatm, tgo_x_ttatmg, satm_val) = if iband == DNS_BAND9 {
            (0.0f32, 1.0f32, 0.0f32)
        } else {
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
                0.05,     // raot550nm
                pressure,
                uoz,
                uwv,
                0.0,      // rotoa (unused for coefficient extraction)
                lambda,
                max_band_idx,
                -1.0,     // eps=-1.0: no wavelength scaling
            );
            (
                result.tgo as f32 * result.roatm as f32,
                result.tgo as f32 * result.ttatmg as f32,
                result.satm as f32,
            )
        };

        for pix in 0..npix {
            let (i, j) = (pix / nsamps, pix % nsamps);
            if is_fill_pixel(qaband[pix]) {
                continue;
            }

            // Sentinel TOA is NOT divided by cos(SZA), used directly
            let rotoa = toa_bands[iband][(i, j)];

            let roslamb: f32 = {
                let num = rotoa - tgo_x_roatm;
                num / (tgo_x_ttatmg + satm_val * num)
            };
            sband[iband][pix] = roslamb;
        }
    }

    // ── Step 7: Aerosol retrieval loop ──
    // Precompute the quadratic fit constants for eps optimization.
    let eps1 = LOW_EPS;
    let eps2 = MOD_EPS;
    let eps3 = HIGH_EPS;
    let xa = eps1 * eps1 - eps3 * eps3;
    let xd = ((eps2 * eps2 - eps3 * eps3) as i32) as f64; // integer truncation
    let xb = eps1 - eps3;
    let xe = eps2 - eps3;

    // Pre-collect window center coordinates: Sentinel starts at (0,0), stepping by aero_window
    let window_centers: Vec<(usize, usize)> = {
        let mut coords = Vec::new();
        let mut wi = 0;
        while wi < nlines {
            let mut wj = 0;
            while wj < nsamps {
                coords.push((wi, wj));
                wj += aero_window;
            }
            wi += aero_window;
        }
        coords
    };

    // Parallel aerosol retrieval at window centers
    let aero_results: Vec<SentinelAerosolWindowResult> = pool.install(|| {
        window_centers.par_iter().map(|&(win_i, win_j)| {
            let curr_pix = win_i * nsamps + win_j;

            // Skip fill pixels
            if is_fill_pixel(qaband[curr_pix]) {
                return SentinelAerosolWindowResult {
                    win_i,
                    win_j,
                    curr_pix,
                    taero: 0.0,
                    teps: 0.0,
                    ipflag: 1u8 << IPFLAG_FILL,
                    is_fill: true,
                };
            }

            // Compute per-pixel lat/lon from image coordinates
            let (pixel_lat, pixel_lon) = utm_to_deg(space_def, win_i as i32, win_j as i32);

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

            // Bilinearly interpolate slopes and intercepts
            // C: float slprb1/intrb1 etc. — use f32 bilinear interpolation
            let slprb1 = bilerp_f32(slp_b1, &cmg.w);
            let intrb1 = bilerp_f32(int_b1, &cmg.w);
            let slprb2 = bilerp_f32(slp_b2, &cmg.w);
            let intrb2 = bilerp_f32(int_b2, &cmg.w);
            let slprb7 = bilerp_f32(slp_b7, &cmg.w);
            let intrb7 = bilerp_f32(int_b7, &cmg.w);

            // Compute NDWI from climatological SR: B8A and B12
            // C: xndwi = ((double)sband[8A] - (double)(sband[12]*0.5)) /
            //            ((double)sband[8A] + (double)(sband[12]*0.5))
            // computed in double then stored in float
            let sr_nir = sband[DNS_BAND8A][curr_pix] as f64;
            let sr_swir2 = sband[DNS_BAND12][curr_pix] as f64;
            let sr_swir2_half = sr_swir2 * 0.5;
            let denom = sr_nir + sr_swir2_half;
            let mut xndwi: f32 = if denom.abs() > 1.0e-10 {
                ((sr_nir - sr_swir2_half) / denom) as f32
            } else {
                0.0f32
            };

            // Clamp NDWI using andwi/sndwi thresholds from CMG
            // C: float ndwi_th1/th2
            let andwi_val = safe_read(&aux.andwi, ratio_pix11);
            let sndwi_val = safe_read(&aux.sndwi, ratio_pix11);
            let ndwi_th1 = ((andwi_val as f64 + 2.0 * sndwi_val as f64) * 0.001) as f32;
            let ndwi_th2 = ((andwi_val as f64 - 2.0 * sndwi_val as f64) * 0.001) as f32;
            if xndwi > ndwi_th1 {
                xndwi = ndwi_th1;
            }
            if xndwi < ndwi_th2 {
                xndwi = ndwi_th2;
            }

            // Initialize erelc and troatm arrays — C uses float for both
            let mut erelc = vec![-1.0f64; nbands];
            let mut troatm = vec![0.0f64; nbands];

            // Compute band ratios from NDWI
            erelc[DNS_BAND1] = (xndwi * slprb1 + intrb1) as f32 as f64;
            erelc[DNS_BAND2] = (xndwi * slprb2 + intrb2) as f32 as f64;
            erelc[DNS_BAND4] = 1.0;
            erelc[DNS_BAND12] = (xndwi * slprb7 + intrb7) as f32 as f64;

            // Average TOA in 6x6 window for aerosol retrieval bands
            // Sentinel does NOT divide by cos(SZA)
            // C accumulates in float, so use f32 accumulators for precision match
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
                    if is_fill_pixel(qaband[win_pix]) {
                        continue;
                    }
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

            // === Eps optimization: 3 retrievals at eps1=1.0, eps2=1.75, eps3=2.5 ===
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

            // Sentinel eps quadratic optimization (C lines 1430-1454)
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

            // Run final subaeroret_new at chosen eps
            let result_final = subaeroret_new(
                false, iband1, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
                &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
                lambda, eps, iaots, tth,
            );
            let raot = result_final.raot;
            let residual = result_final.residual;

            let mut result_taero = raot as f32;
            let mut result_teps = eps as f32;
            let corf = raot / xmus;

            // === Post-retrieval validation ===
            let mut result_ipflag: u8 = 0;
            if residual < (0.015 + 0.005 * corf + 0.10 * troatm[DNS_BAND12]) {
                // Average TOA in NxN window for B8A
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

                // Average TOA in NxN window for B04
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
                    // Clear pixel with valid aerosol retrieval
                    result_ipflag |= 1u8 << IPFLAG_CLEAR;
                } else {
                    result_ipflag = 1u8 << IPFLAG_WATER;
                }
            } else {
                result_ipflag = 1u8 << IPFLAG_WATER;
            }

            // === Water retest ===
            if result_ipflag & (1u8 << IPFLAG_WATER) != 0 {
                // Reset erelc and troatm for water retrieval
                let mut water_erelc = vec![-1.0f64; nbands];
                let mut water_troatm = vec![0.0f64; nbands];

                // Average TOA for water bands B01/B04/B8A/B12 (f32 precision)
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

                // Water band ratios: all 1.0
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
                result_teps = WATER_EPS as f32;
                result_taero = water_result.raot as f32;
                let water_corf = water_result.raot / xmus;

                // Validate: check band 1 reflectance
                let ros1 = atmcorlamb2_new(
                    &atm_coeff[DNS_BAND1], tgo_arr[DNS_BAND1], DNS_BAND1,
                    water_result.raot, normext_p0a3[DNS_BAND1],
                    water_troatm[DNS_BAND1], lambda, WATER_EPS,
                );

                if water_result.residual > (0.010 + 0.005 * water_corf) || ros1 < 0.0 {
                    // Not valid water — mark as failed
                    result_ipflag = 1u8 << IPFLAG_FAILED;
                } else {
                    // Valid water pixel
                    result_ipflag = (1u8 << IPFLAG_WATER) | (1u8 << IPFLAG_CLEAR);
                }
            }

            SentinelAerosolWindowResult {
                win_i,
                win_j,
                curr_pix,
                taero: result_taero,
                teps: result_teps,
                ipflag: result_ipflag,
                is_fill: false,
            }
        }).collect()
    });

    // Scatter results back to the output arrays, broadcasting to 6x6 windows
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

    // ── Step 8: Post-processing ──
    aerosol_interp_sentinel(aero_window, &qaband, &mut ipflag, &mut taero, nlines, nsamps);
    ipflag_expand_failed_sentinel(&mut ipflag, nlines, nsamps);
    aero_avg_failed_sentinel(&qaband, &mut ipflag, &mut taero, &mut teps, nlines, nsamps);

    // ── Step 9: Final per-pixel atmospheric correction ──
    let sr_flat: Vec<Vec<f64>> = (0..nbands)
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

                    // SAFETY: Each iline is processed by exactly one thread
                    // (Rayon's into_par_iter guarantees this). Within a thread,
                    // pix = iline * nsamps + isamp produces unique indices for
                    // that row, so no two threads write to the same index.
                    unsafe {
                        let ptr = sr_flat[iband].as_ptr() as *mut f64;
                        *ptr.add(pix) = val;
                    }
                }
            }
        });
    });

    // Serial post-pass: set aerosol QA bits on B01 (not B00 like Landsat)
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

    // ── Step 10: Scale to output integers ──
    let offset_f32 = BAND_OFFSET_REFL as f32;
    let mult_f32 = MULT_FACTOR_REFL as f32;
    let sr_bands: Vec<Array2<u16>> = sr_flat
        .iter()
        .map(|band| {
            Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
                let pix = i * nsamps + j;
                if is_fill_pixel(qaband[pix]) {
                    0u16
                } else {
                    let sband_f32 = band[pix] as f32;
                    let tmpf = (sband_f32 + offset_f32) * mult_f32;
                    tmpf.round().clamp(0.0, u16::MAX as f32) as u16
                }
            })
        })
        .collect();

    // No brightness temperature bands for Sentinel
    let bt_out: Vec<Array2<u16>> = Vec::new();

    // Scale aerosol to int16
    let aero_mult_f32 = MULT_FACTOR_AERO as f32;
    let aerosol = Array2::from_shape_fn((nlines, nsamps), |(i, j)| {
        let pix = i * nsamps + j;
        if is_fill_pixel(qaband[pix]) {
            AERO_FILL
        } else {
            let tmpf = taero[pix] * aero_mult_f32;
            tmpf.round().clamp(0.0, 5000.0) as i16
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
