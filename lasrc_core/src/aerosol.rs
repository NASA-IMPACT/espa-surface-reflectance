//! Aerosol retrieval routines ported from C LaSRC code.
//!
//! Implements `subaeroret_new`, `aerosol_interp`, and `fix_invalid_aerosols`.

use crate::atmospheric::{AtmCorrCoefficients, atmcorlamb2_new};
use crate::constants::*;

/// Result of aerosol retrieval for a single pixel.
#[derive(Debug, Clone)]
pub struct AerosolResult {
    /// Retrieved aerosol optical thickness at 550 nm
    pub raot: f64,
    /// Residual (fit quality)
    pub residual: f64,
    /// Angstrom exponent
    pub eps: f64,
    /// AOT LUT index corresponding to the retrieved AOT
    pub iaots: usize,
}

/// Semi-empirical aerosol retrieval via spectral band ratios.
///
/// Ported from `subaeroret_new()` in C `subaeroret.c`.
/// Matches the C algorithm exactly: convergence loop that stops when
/// residual starts increasing, then 3-point parabolic refinement.
///
/// # Arguments
/// * `is_water`        - Whether the pixel is classified as water
/// * `iband1`          - Index of the reference band
/// * `erelc`           - Expected band ratio coefficients
/// * `troatm`          - TOA reflectance per band
/// * `tgo_arr`         - Gas transmission per band
/// * `roatm_ia_max`    - Atmospheric reflectance at maximum AOT per band
/// * `roatm_coef`      - Polynomial coefficients for atmospheric reflectance, per band
/// * `ttatmg_coef`     - Polynomial coefficients for total transmittance, per band
/// * `satm_coef`       - Polynomial coefficients for spherical albedo, per band
/// * `normext_p0a3`    - Normalised aerosol extinction at reference AOT, per band
/// * `lambda`          - Wavelength table (µm), one entry per band
/// * `eps`             - Angstrom exponent
/// * `iaots`           - Starting AOT LUT index
/// * `tth`             - Threshold surface reflectance per band (negative check)
///
/// # Returns
/// [`AerosolResult`] with the best-fit AOT, residual, eps, and LUT index.
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
    let nband = erelc.len();

    // Helper to compute atmospheric correction for one band.
    // Returns f32 roslamb matching C's `float roslamb`.
    let atm_corr = |ib: usize, raot: f32| -> f32 {
        let coeff = AtmCorrCoefficients {
            roatm_upper: roatm_ia_max[ib],
            roatm_coef: roatm_coef[ib],
            ttatmg_coef: ttatmg_coef[ib],
            satm_coef: satm_coef[ib],
        };
        // atmcorlamb2_new already uses f32 internally and returns f64;
        // truncate back to f32 to match C's `float roslamb`.
        atmcorlamb2_new(&coeff, tgo_arr[ib], ib, raot as f64, normext_p0a3[ib], troatm[ib], lambda, eps) as f32
    };

    // Helper to compute residual at a given raot550nm value.
    // Returns (residual_f32, ros1_f64, testth).
    //
    // C precision: roslamb is float, ros1 is double, *residual is float,
    // point_error is double, erelc[] is float.
    let compute_residual = |raot: f32| -> (f32, f64, bool) {
        // Reference band correction — roslamb is float, ros1 is double
        let roslamb = atm_corr(iband1, raot);
        let mut testth = (roslamb as f64) - tth[iband1] < 0.0;
        let ros1: f64 = roslamb as f64; // C: double ros1 = (float)roslamb

        // C: *residual is float — accumulate into f32
        let mut residual: f32 = 0.0;
        let mut nbval = 0usize;

        // For water: include iband1 in loop. For land: skip iband1.
        if is_water {
            for ib in 0..nband {
                if erelc[ib] > 0.0 {
                    let roslamb = atm_corr(ib, raot);
                    if (roslamb as f64) - tth[ib] < 0.0 {
                        testth = true;
                    }
                    // C: *residual += roslamb*roslamb (float += float*float)
                    residual += roslamb * roslamb;
                    nbval += 1;
                }
            }
        } else {
            for ib in 0..nband {
                // C: if (ib != iband1 && erelc[ib] > 0.0)
                if ib != iband1 && erelc[ib] > 0.0 {
                    let roslamb = atm_corr(ib, raot);
                    if (roslamb as f64) - tth[ib] < 0.0 {
                        testth = true;
                    }
                    // C: point_error = roslamb - erelc[ib] * ros1 (double)
                    // roslamb is float promoted to double, erelc[ib] is float promoted to double
                    let point_error: f64 = roslamb as f64 - (erelc[ib] as f32 as f64) * ros1;
                    // C: *residual += point_error * point_error (float += double)
                    residual += (point_error * point_error) as f32;
                    nbval += 1;
                }
            }
        }

        if nbval > 0 {
            // C: *residual = sqrt(*residual) / nbval
            // sqrt promotes float to double, divides, truncates back to float
            residual = ((residual as f64).sqrt() / nbval as f64) as f32;
        }

        (residual, ros1, testth)
    };

    // C variables: residual1/2 are double, raot1/2 are double
    let mut residual1 = 2000.0f64;
    let mut residual2 = 1000.0f64;
    let mut iaot1 = 0usize;
    let mut iaot2 = 0usize;
    let mut raot1 = 0.0001f64;
    let mut raot2 = 1.0e-6f64;

    // First iteration at iaots
    // C: float raot550nm = aot550nm[iaot]
    let mut iaot = iaots;
    let mut raot550nm: f32 = AOT550NM[iaot] as f32;
    let (mut residual, _ros1, mut testth) = compute_residual(raot550nm);

    // Convergence loop: increment iaot, stop when residual starts increasing
    // or testth is triggered.
    // C: while ((iaot < NAOT_VALS) && (*residual < residual1) && (!testth))
    // Note: *residual is float, residual1 is double — C promotes float to double for comparison
    iaot += 1;
    while iaot < NAOT_VALS && (residual as f64) < residual1 && !testth {
        // Shift history: residual1/2 are double, store float-precision residual promoted to double
        residual2 = residual1;
        iaot2 = iaot1;
        raot2 = raot1;
        residual1 = residual as f64; // C: residual1 = *residual (float→double)
        raot1 = raot550nm as f64;    // C: raot1 = raot550nm (float→double)
        iaot1 = iaot;

        raot550nm = AOT550NM[iaot] as f32;
        let (new_res, _ros1, new_testth) = compute_residual(raot550nm);
        residual = new_res;
        testth = new_testth;

        iaot += 1;
    }

    // Parabolic refinement
    // C: *raot and raotsaved are float; xa, xb, raotmin are double; residualm is double
    let mut final_raot: f32;
    let final_residual: f64;
    let final_iaots;

    if iaot <= 1 {
        // No convergence achieved — use the current AOT value
        // C: if (iaot == 1) { *raot = raot550nm; *iaots unchanged }
        final_raot = raot550nm;
        final_residual = residual as f64;
        final_iaots = iaots; // C does NOT modify *iaots in this case
    } else {
        // C: *raot = raot550nm; raotsaved = *raot;
        final_raot = raot550nm;
        let raotsaved: f32 = final_raot;

        // 3-point quadratic fit using double precision (matching C's double xa, xb, raotmin)
        // C: xa = (residual1 - *residual)*(raot2 - *raot)
        // *residual is float (promoted to double), *raot is float (promoted to double)
        let xa_fit = (residual1 - residual as f64) * (raot2 - final_raot as f64);
        let xb_fit = (residual2 - residual as f64) * (raot1 - final_raot as f64);
        let denom = xa_fit - xb_fit;
        let raotmin = if denom.abs() > 1e-20 {
            0.5 * (xa_fit * (raot2 + final_raot as f64) - xb_fit * (raot1 + final_raot as f64)) / denom
        } else {
            final_raot as f64
        };

        let raotmin = if raotmin < 0.01 || raotmin > 4.0 {
            final_raot as f64
        } else {
            raotmin
        };

        // Evaluate residual at the parabolic minimum
        // C: raot550nm = raotmin (double→float truncation)
        let (residualm_f32, _ros1, _testth) = compute_residual(raotmin as f32);
        let mut residualm: f64 = residualm_f32 as f64;
        // C: *raot = raot550nm (which is raotmin truncated to float)
        final_raot = raotmin as f32;

        // Compare against the three stored values and pick the best
        // C: if (residualm > *residual) { residualm = *residual; *raot = raotsaved; }
        if residualm > residual as f64 {
            residualm = residual as f64;
            final_raot = raotsaved;
        }
        if residualm > residual1 {
            residualm = residual1;
            final_raot = raot1 as f32; // C: *raot = raot1 (double→float)
        }
        if residualm > residual2 {
            residualm = residual2;
            final_raot = raot2 as f32; // C: *raot = raot2 (double→float)
        }
        final_residual = residualm;

        // C: *iaots = MAX((iaot2 - 3), 0)
        // Special case for water: if iaot == 1, iaots = 0
        if is_water && iaot == 1 {
            final_iaots = 0;
        } else {
            final_iaots = if iaot2 >= 3 { iaot2 - 3 } else { 0 };
        }
    }

    AerosolResult {
        raot: final_raot as f64,
        residual: final_residual,
        eps,
        iaots: final_iaots,
    }
}

/// Bilinear interpolation of a per-pixel array from window-center values.
///
/// Ported from `aerosol_interp_landsat()` in C LaSRC code (`aero_interp.c`).
/// Called once for taero, once for teps, matching the C calling convention.
///
/// For each non-center, non-fill pixel, find the four surrounding window centers
/// and bilinearly interpolate.  Center pixels are skipped (they keep their
/// existing values and ipflag).  Non-center pixels get ipflag **assigned** to
/// `IPFLAG_INTERP_WINDOW` (clearing all other bits), with `IPFLAG_WATER` OR'd
/// in if any surrounding center was water.  Fill pixels are cleaned up at the
/// end to have only `IPFLAG_FILL`.
pub fn aerosol_interp(
    data: &mut [f32],
    ipflag: &mut [u8],
    qa_band: &[u16],
    nlines: usize,
    nsamps: usize,
    aero_window: usize,
    half_aero_window: usize,
) {
    let aero_step: f64 = 1.0 / aero_window as f64;
    let water_bit = 1u8 << IPFLAG_WATER;

    for line in 0..nlines {
        // C: center_line = (int)(line * aero_step) * aero_window + half_aero_window
        let center_line = ((line as f64 * aero_step) as isize) as usize * aero_window
            + half_aero_window;

        // Fractional distance and neighbor line
        // C: float yaero, xaero — use f32 to match
        let yaero = ((line as f64 - center_line as f64) * aero_step) as f32;
        let u_signed = yaero - (yaero as i32) as f32; // C: yaero - (int)yaero
        let center_line1 = if u_signed < 0.0 {
            if center_line >= aero_window {
                center_line - aero_window
            } else {
                center_line
            }
        } else {
            let cl1 = center_line + aero_window;
            if cl1 >= nlines.saturating_sub(1) {
                center_line
            } else {
                cl1
            }
        };
        let u: f32 = u_signed.abs();

        for samp in 0..nsamps {
            let pix = line * nsamps + samp;

            // Skip fill pixels
            if is_fill_pixel(qa_band[pix]) {
                continue;
            }

            // C: center_samp = (int)(samp * aero_step) * aero_window + half_aero_window
            let center_samp = ((samp as f64 * aero_step) as isize) as usize * aero_window
                + half_aero_window;

            // Skip center pixels — they already have their values
            if samp == center_samp && line == center_line {
                continue;
            }

            // Fractional distance and neighbor sample (C: float xaero)
            let xaero = ((samp as f64 - center_samp as f64) * aero_step) as f32;
            let v_signed = xaero - (xaero as i32) as f32;
            let center_samp1 = if v_signed < 0.0 {
                if center_samp >= aero_window {
                    center_samp - aero_window
                } else {
                    center_samp
                }
            } else {
                let cs1 = center_samp + aero_window;
                if cs1 >= nsamps.saturating_sub(1) {
                    center_samp
                } else {
                    cs1
                }
            };
            let v: f32 = v_signed.abs();

            // Four corner center pixels
            let pix11 = center_line * nsamps + center_samp;
            let pix12 = center_line * nsamps + center_samp1;
            let pix21 = center_line1 * nsamps + center_samp;
            let pix22 = center_line1 * nsamps + center_samp1;

            // Bilinear interpolation in f32 (C: float arithmetic on float arrays)
            let a11 = data[pix11];
            let a12 = data[pix12];
            let a21 = data[pix21];
            let a22 = data[pix22];
            data[pix] = a11
                + u * (a21 - a11)
                + v * (a12 - a11)
                + u * v * (a11 - a12 - a21 + a22);

            // Set ipflag: clear everything, set INTERP_WINDOW
            ipflag[pix] = 1u8 << IPFLAG_INTERP_WINDOW;

            // If any corner center was water, mark this pixel as water too
            if (ipflag[pix11] & water_bit != 0)
                || (ipflag[pix12] & water_bit != 0)
                || (ipflag[pix21] & water_bit != 0)
                || (ipflag[pix22] & water_bit != 0)
            {
                ipflag[pix] |= water_bit;
            }
        }
    }

    // Clean up fill pixels: ensure they only have IPFLAG_FILL
    let npix = nlines * nsamps;
    for pix in 0..npix {
        if is_fill_pixel(qa_band[pix]) {
            ipflag[pix] = 1u8 << IPFLAG_FILL;
        }
    }
}

/// Three-pass local averaging to fix pixels with invalid aerosol retrievals.
///
/// Ported from `fix_invalid_aerosols_landsat()` and `fill_with_local_average_landsat()`
/// in C LaSRC code (`aero_interp.c`).
///
/// Only iterates over **center pixels** (stepping by `aero_window`, starting at
/// `half_aero_window`).  Uses a separate `smflag` array to track filled pixels
/// (C Landsat never sets IPFLAG_FIXED).
///
/// - Pass 1 (forward):  require `min_clear_pix` valid neighbors, don't use filled.
/// - Pass 2 (forward):  require 1 valid pixel, also allow previously-filled pixels.
/// - Pass 3 (reverse):  same as pass 2 but iterate in reverse order.
///
/// # Arguments
/// * `taero`                - Aerosol optical thickness (modified in place)
/// * `teps`                 - Angstrom exponent (modified in place)
/// * `ipflag`               - Processing flags (NOT modified — Landsat never sets IPFLAG_FIXED)
/// * `nlines`               - Number of image lines
/// * `nsamps`               - Number of image samples
/// * `aero_window`          - Aerosol retrieval window size (NxN)
/// * `half_aero_window`     - Half of aerosol window
/// * `fix_aero_window`      - Search window for fixing invalid retrievals (WxW)
/// * `half_fix_aero_window` - Half of fix window
/// * `min_clear_pix`        - Minimum number of clear neighbors required in pass 1
pub fn fix_invalid_aerosols(
    taero: &mut [f32],
    teps: &mut [f32],
    ipflag: &[u8],
    nlines: usize,
    nsamps: usize,
    aero_window: usize,
    half_aero_window: usize,
    _fix_aero_window: usize,
    half_fix_aero_window: usize,
    min_clear_pix: usize,
) {
    let valid_bit = 1u8 << IPFLAG_CLEAR; // lasrc_qa_is_valid_aerosol_retrieval
    let npix = nlines * nsamps;
    let mut smflag = vec![false; npix];

    // C: window_offset = LHALF_FIX_AERO_WINDOW - half_aero_window
    let window_offset = half_fix_aero_window as isize - half_aero_window as isize;
    let step = aero_window as isize;

    // Inner function matching fill_with_local_average_landsat
    let fill_pass = |forward: bool,
                     required_clear: usize,
                     use_filled: bool,
                     taero: &mut [f32],
                     teps: &mut [f32],
                     smflag: &mut [bool]| {
        let (start_line, start_samp, step_val): (isize, isize, isize) = if forward {
            (half_aero_window as isize, half_aero_window as isize, step)
        } else {
            // C: start = ((int)round((double)nlines / aero_window) - 1) * aero_window + half_aero_window
            let sl = ((nlines as f64 / aero_window as f64).round() as isize - 1)
                * aero_window as isize
                + half_aero_window as isize;
            let ss = ((nsamps as f64 / aero_window as f64).round() as isize - 1)
                * aero_window as isize
                + half_aero_window as isize;
            (sl, ss, -step)
        };

        let mut line = start_line;
        while line > 0 && line < nlines as isize {
            let mut samp = start_samp;
            while samp > 0 && samp < nsamps as isize {
                let curr_pix = line as usize * nsamps + samp as usize;

                // Skip fill pixels
                if ipflag[curr_pix] & (1u8 << IPFLAG_FILL) != 0 {
                    samp += step_val;
                    continue;
                }

                // Skip already-valid or already-filled
                if (ipflag[curr_pix] & valid_bit != 0) || smflag[curr_pix] {
                    samp += step_val;
                    continue;
                }

                // Search WxW window around current center pixel, stepping by aero_window
                // C uses double for sum accumulators
                let mut sum_aero: f64 = 0.0;
                let mut sum_eps: f64 = 0.0;
                let mut nbclrpix = 0usize;

                let mut iline = line - window_offset;
                while iline <= line + window_offset {
                    if iline < 0 || iline >= nlines as isize {
                        iline += step;
                        continue;
                    }
                    let ilpix = iline as usize * nsamps;

                    let mut isamp = samp - window_offset;
                    while isamp <= samp + window_offset {
                        if isamp < 0 || isamp >= nsamps as isize {
                            isamp += step;
                            continue;
                        }
                        let ipix = ilpix + isamp as usize;

                        if (ipflag[ipix] & valid_bit != 0)
                            || (use_filled && smflag[ipix])
                        {
                            nbclrpix += 1;
                            sum_aero += taero[ipix] as f64;
                            sum_eps += teps[ipix] as f64;
                        }
                        isamp += step;
                    }
                    iline += step;
                }

                if nbclrpix >= required_clear {
                    // C: taero[curr_pix] = sum_aero / nbclrpix;
                    // double / int -> double, assigned to float
                    taero[curr_pix] = (sum_aero / nbclrpix as f64) as f32;
                    teps[curr_pix] = (sum_eps / nbclrpix as f64) as f32;
                    smflag[curr_pix] = true;
                } else {
                    taero[curr_pix] = DEFAULT_AERO as f32;
                    teps[curr_pix] = DEFAULT_EPS as f32;
                }

                samp += step_val;
            }
            line += step_val;
        }
    };

    // Pass 1: forward, require min_clear_pix, don't use filled
    fill_pass(true, min_clear_pix, false, taero, teps, &mut smflag);

    // Pass 2: forward, require 1, use filled
    fill_pass(true, 1, true, taero, teps, &mut smflag);

    // Pass 3: reverse, require 1, use filled
    fill_pass(false, 1, true, taero, teps, &mut smflag);
}

/// Bilinear interpolation of aerosol from window UL corners for Sentinel-2.
///
/// Ported from `aerosol_interp_sentinel()` in C `aero_interp.c:527-637`.
pub fn aerosol_interp_sentinel(
    aero_window: usize,
    qaband: &[u16],
    _ipflag: &mut [u8],
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

/// Expand failed aerosol pixels to surrounding area for Sentinel-2.
///
/// Ported from `ipflag_expand_failed_sentinel()` in C `aero_interp.c:658-722`.
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

/// Average aerosol values for failed Sentinel-2 pixels.
///
/// Ported from `aero_avg_failed_sentinel()` in C `aero_interp.c:741-927`.
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

    // Copy averaged values back ONLY for FAILED pixels (not fill).
    // C: only replaces pixels with IPFLAG_FAILED set, and marks them IPFLAG_FIXED.
    for pix in 0..npixels {
        if ipflag[pix] & (1u8 << IPFLAG_FAILED) != 0
            && ipflag[pix] & (1u8 << IPFLAG_FILL) == 0
        {
            taero[pix] = taeros[pix];
            teps[pix] = tepss[pix];
            ipflag[pix] |= 1u8 << IPFLAG_FIXED;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_subaeroret_new_returns_valid_aot() {
        let nband = 7;
        let erelc = vec![1.0, 0.8, 0.0, 0.5, 0.0, 0.0, 0.2];
        let troatm = vec![0.15, 0.14, 0.12, 0.11, 0.08, 0.04, 0.03];
        let tgo_arr = vec![0.95; nband];
        let roatm_ia_max = vec![0.4; nband];
        let roatm_coef: Vec<[f64; NCOEF]> = (0..nband).map(|_| [0.0, 0.0, 0.15, 0.01]).collect();
        let ttatmg_coef: Vec<[f64; NCOEF]> = (0..nband).map(|_| [0.0, 0.0, -0.1, 0.85]).collect();
        let satm_coef: Vec<[f64; NCOEF]> = (0..nband).map(|_| [0.0, 0.0, 0.05, 0.005]).collect();
        let normext_p0a3 = vec![1.0; nband];
        let tth = vec![1.0e-3, 1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-4];
        let result = subaeroret_new(
            false, 0, &erelc, &troatm, &tgo_arr, &roatm_ia_max,
            &roatm_coef, &ttatmg_coef, &satm_coef, &normext_p0a3,
            &LAMBDA_LANDSAT, 1.5, 0, &tth,
        );
        assert!(result.raot >= 0.01, "AOT should be >= 0.01: {}", result.raot);
        assert!(result.raot <= 5.0, "AOT should be <= 5.0: {}", result.raot);
        assert!(result.residual >= 0.0);
    }

    #[test]
    fn test_aerosol_interp_fills_gaps() {
        let nlines = 9;
        let nsamps = 9;
        let npix = nlines * nsamps;
        let mut taero = vec![0.0f32; npix];
        let mut teps = vec![0.0f32; npix];
        let mut ipflag = vec![0u8; npix];
        let qa_band = vec![2u16; npix]; // bit 0 = 0 means not fill
        for i in (1..nlines).step_by(3) {
            for j in (1..nsamps).step_by(3) {
                let pix = i * nsamps + j;
                taero[pix] = 0.1;
                teps[pix] = 1.5;
                ipflag[pix] = 1 << IPFLAG_CLEAR;
            }
        }
        aerosol_interp(&mut taero, &mut ipflag, &qa_band, nlines, nsamps, 3, 1);
        aerosol_interp(&mut teps, &mut ipflag, &qa_band, nlines, nsamps, 3, 1);
        // Check that non-center pixels got interpolated values
        let non_center_pix = 0 * nsamps + 0; // pixel (0,0)
        assert!(taero[non_center_pix] > 0.0, "Pixel (0,0) should have interpolated aerosol");
    }

    #[test]
    fn test_fix_invalid_aerosols_fills() {
        let nlines = 9;
        let nsamps = 9;
        let npix = nlines * nsamps;
        let mut taero = vec![0.0f32; npix];
        let mut teps = vec![0.0f32; npix];
        let mut ipflag = vec![0u8; npix];
        // Set some valid centers (aero_window=3, half=1, so centers at 1,4,7)
        for i in (1..nlines).step_by(3) {
            for j in (1..nsamps).step_by(3) {
                let pix = i * nsamps + j;
                taero[pix] = 0.2;
                teps[pix] = 1.5;
                ipflag[pix] = 1 << IPFLAG_CLEAR;
            }
        }
        // Clear one center to simulate invalid retrieval
        let bad_pix = 4 * nsamps + 4;
        ipflag[bad_pix] = 0;
        taero[bad_pix] = 0.0;

        fix_invalid_aerosols(&mut taero, &mut teps, &ipflag, nlines, nsamps, 3, 1, 15, 7, 4);
        // Landsat never sets IPFLAG_FIXED — but taero should be filled
        assert!(taero[bad_pix] > 0.0, "Should have filled value");
    }

    #[test]
    fn test_aerosol_interp_sentinel_basic() {
        let nlines = 12;
        let nsamps = 12;
        let npix = nlines * nsamps;
        let mut taero = vec![0.0f32; npix];
        let mut ipflag = vec![0u8; npix];
        let qaband = vec![0u16; npix]; // no fill (bit 0 = 0)

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
}
