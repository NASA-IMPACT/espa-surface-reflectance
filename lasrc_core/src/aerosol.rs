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
/// Ported from `subaeroret_new()` in C `lut_subr.c`.
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

    // Track the best result across all AOT steps
    let mut best_raot = AOT550NM[iaots];
    let mut best_residual = f64::MAX;
    let mut best_iaots = iaots;

    let mut prev_residual = f64::MAX;

    for iaot in iaots..NAOT_VALS {
        let raot = AOT550NM[iaot];

        // Build AtmCorrCoefficients for the reference band
        let coeff1 = AtmCorrCoefficients {
            roatm_upper: roatm_ia_max[iband1],
            roatm_coef: roatm_coef[iband1],
            ttatmg_coef: ttatmg_coef[iband1],
            satm_coef: satm_coef[iband1],
        };

        // Correct the reference band
        let ros1 = atmcorlamb2_new(
            &coeff1,
            tgo_arr[iband1],
            iband1,
            raot,
            normext_p0a3[iband1],
            troatm[iband1],
            lambda,
            eps,
        );

        // Test if surface reflectance is negative (below threshold)
        let testth = ros1 - tth[iband1] < 0.0;

        // Accumulate residual across active bands
        let mut residual = 0.0;
        let mut nbval: usize = 0;

        for ib in 0..nband {
            if erelc[ib] < 0.0 {
                continue;
            }

            let coeff_ib = AtmCorrCoefficients {
                roatm_upper: roatm_ia_max[ib],
                roatm_coef: roatm_coef[ib],
                ttatmg_coef: ttatmg_coef[ib],
                satm_coef: satm_coef[ib],
            };

            let roslamb = atmcorlamb2_new(
                &coeff_ib,
                tgo_arr[ib],
                ib,
                raot,
                normext_p0a3[ib],
                troatm[ib],
                lambda,
                eps,
            );

            if is_water {
                residual += roslamb * roslamb;
            } else {
                let diff = roslamb - erelc[ib] * ros1;
                residual += diff * diff;
            }
            nbval += 1;
        }

        if nbval > 0 {
            residual = residual.sqrt() / nbval as f64;
        }

        // Track the global minimum
        if residual < best_residual {
            best_residual = residual;
            best_raot = raot;
            best_iaots = iaot;
        }

        // Break if we passed the minimum (residual is increasing or threshold triggered)
        if testth || residual >= prev_residual {
            break;
        }

        prev_residual = residual;
    }

    // --- Parabolic refinement using the three LUT points surrounding the minimum ---
    // Re-evaluate residuals at the three surrounding AOT indices to fit a quadratic.

    // Use the three surrounding AOT steps around best_iaots for refinement
    let ia_lo = if best_iaots > 0 { best_iaots - 1 } else { best_iaots };
    let ia_hi = if best_iaots + 1 < NAOT_VALS { best_iaots + 1 } else { best_iaots };

    let pts: [(f64, f64); 3] = {
        let mut arr = [(0.0f64, 0.0f64); 3];
        for (k, &ia) in [ia_lo, best_iaots, ia_hi].iter().enumerate() {
            let raot_pt = AOT550NM[ia];
            let coeff1 = AtmCorrCoefficients {
                roatm_upper: roatm_ia_max[iband1],
                roatm_coef: roatm_coef[iband1],
                ttatmg_coef: ttatmg_coef[iband1],
                satm_coef: satm_coef[iband1],
            };
            let ros1 = atmcorlamb2_new(
                &coeff1,
                tgo_arr[iband1],
                iband1,
                raot_pt,
                normext_p0a3[iband1],
                troatm[iband1],
                lambda,
                eps,
            );

            let mut res = 0.0;
            let mut nbval: usize = 0;
            for ib in 0..nband {
                if erelc[ib] < 0.0 {
                    continue;
                }
                let coeff_ib = AtmCorrCoefficients {
                    roatm_upper: roatm_ia_max[ib],
                    roatm_coef: roatm_coef[ib],
                    ttatmg_coef: ttatmg_coef[ib],
                    satm_coef: satm_coef[ib],
                };
                let roslamb = atmcorlamb2_new(
                    &coeff_ib,
                    tgo_arr[ib],
                    ib,
                    raot_pt,
                    normext_p0a3[ib],
                    troatm[ib],
                    lambda,
                    eps,
                );
                if is_water {
                    res += roslamb * roslamb;
                } else {
                    let diff = roslamb - erelc[ib] * ros1;
                    res += diff * diff;
                }
                nbval += 1;
            }
            if nbval > 0 {
                res = res.sqrt() / nbval as f64;
            }
            arr[k] = (raot_pt, res);
        }
        arr
    };

    let (x1, y1) = pts[0];
    let (x2, y2) = pts[1];
    let (x3, y3) = pts[2];

    // Fit quadratic y = a*x^2 + b*x + c through three points.
    // Only attempt if the three x-values are distinct.
    if (x3 - x1).abs() > 1e-12 && (x2 - x1).abs() > 1e-12 && (x3 - x2).abs() > 1e-12 {
        // Using the standard divided-difference / Lagrange approach
        let denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
        if denom.abs() > 1e-20 {
            let a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
            let b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;

            // Vertex of parabola at x = -b / (2*a)
            if a > 0.0 {
                let raotmin = -b / (2.0 * a);
                if raotmin >= 0.01 && raotmin <= 4.0 {
                    // Evaluate residual at parabolic minimum
                    let coeff1 = AtmCorrCoefficients {
                        roatm_upper: roatm_ia_max[iband1],
                        roatm_coef: roatm_coef[iband1],
                        ttatmg_coef: ttatmg_coef[iband1],
                        satm_coef: satm_coef[iband1],
                    };
                    let ros1_min = atmcorlamb2_new(
                        &coeff1,
                        tgo_arr[iband1],
                        iband1,
                        raotmin,
                        normext_p0a3[iband1],
                        troatm[iband1],
                        lambda,
                        eps,
                    );
                    let mut res_min = 0.0;
                    let mut nbval: usize = 0;
                    for ib in 0..nband {
                        if erelc[ib] < 0.0 {
                            continue;
                        }
                        let coeff_ib = AtmCorrCoefficients {
                            roatm_upper: roatm_ia_max[ib],
                            roatm_coef: roatm_coef[ib],
                            ttatmg_coef: ttatmg_coef[ib],
                            satm_coef: satm_coef[ib],
                        };
                        let roslamb = atmcorlamb2_new(
                            &coeff_ib,
                            tgo_arr[ib],
                            ib,
                            raotmin,
                            normext_p0a3[ib],
                            troatm[ib],
                            lambda,
                            eps,
                        );
                        if is_water {
                            res_min += roslamb * roslamb;
                        } else {
                            let diff = roslamb - erelc[ib] * ros1_min;
                            res_min += diff * diff;
                        }
                        nbval += 1;
                    }
                    if nbval > 0 {
                        res_min = res_min.sqrt() / nbval as f64;
                    }
                    if res_min < best_residual {
                        best_residual = res_min;
                        best_raot = raotmin;
                        // Find the nearest LUT index for iaots
                        let mut nearest = best_iaots;
                        let mut min_diff = (AOT550NM[best_iaots] - raotmin).abs();
                        for ia in 0..NAOT_VALS {
                            let d = (AOT550NM[ia] - raotmin).abs();
                            if d < min_diff {
                                min_diff = d;
                                nearest = ia;
                            }
                        }
                        best_iaots = nearest;
                    }
                }
            }
        }
    }

    // Ensure raot is within valid physical range
    let final_raot = best_raot.clamp(0.01, 4.0);

    AerosolResult {
        raot: final_raot,
        residual: best_residual,
        eps,
        iaots: best_iaots,
    }
}

/// Bilinear interpolation of aerosol properties from window-center retrievals.
///
/// Ported from `aerosol_interp()` in C LaSRC code.
///
/// For each pixel not at a window center, find the four surrounding window centers
/// and bilinearly interpolate `taero` and `teps`. Skip fill pixels.
/// Interpolated pixels are marked with [`IPFLAG_INTERP_WINDOW`].
///
/// # Arguments
/// * `taero`            - Aerosol optical thickness array (nlines * nsamps), modified in place
/// * `teps`             - Angstrom exponent array (nlines * nsamps), modified in place
/// * `ipflag`           - QA/processing-flag array (nlines * nsamps), modified in place
/// * `qa_band`          - QA input band (fill == 0 pixels are skipped)
/// * `nlines`           - Number of image lines
/// * `nsamps`           - Number of image samples
/// * `aero_window`      - Full aerosol window size (pixels)
/// * `half_aero_window` - Half aerosol window size (pixels)
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
    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;

            // Skip fill pixels
            if qa_band[pix] == INPUT_FILL as u16 {
                continue;
            }

            // Find the surrounding window center row indices.
            // Window centers are at half_aero_window, half_aero_window + aero_window, ...
            // i.e. row = half_aero_window + k * aero_window  for integer k >= 0.
            let center_row0 = if iline < half_aero_window {
                half_aero_window
            } else {
                // Largest center row <= iline
                let steps = (iline - half_aero_window) / aero_window;
                half_aero_window + steps * aero_window
            };
            let center_row1 = center_row0 + aero_window;

            let center_col0 = if isamp < half_aero_window {
                half_aero_window
            } else {
                let steps = (isamp - half_aero_window) / aero_window;
                half_aero_window + steps * aero_window
            };
            let center_col1 = center_col0 + aero_window;

            // Clamp to image bounds
            let r0 = center_row0.min(nlines - 1);
            let r1 = center_row1.min(nlines - 1);
            let c0 = center_col0.min(nsamps - 1);
            let c1 = center_col1.min(nsamps - 1);

            // Collect corner values (use only valid IPFLAG_CLEAR pixels)
            let flag_bit = 1u8 << IPFLAG_CLEAR;
            let corners = [
                (r0, c0),
                (r0, c1),
                (r1, c0),
                (r1, c1),
            ];

            // Check if any corner has valid data
            let valid_corners: Vec<(usize, usize)> = corners
                .iter()
                .filter(|&&(r, c)| ipflag[r * nsamps + c] & flag_bit != 0)
                .cloned()
                .collect();

            if valid_corners.is_empty() {
                continue;
            }

            // Bilinear interpolation weights.
            // If r0 == r1 or c0 == c1 (edge of image), degenerate to nearest valid.
            // Use signed arithmetic to avoid usize underflow when the pixel is before
            // the first window center (iline < r0 or isamp < c0).
            let dr = if r1 > r0 {
                (iline as isize - r0 as isize).max(0) as f64 / (r1 - r0) as f64
            } else {
                0.0
            };
            let dc = if c1 > c0 {
                (isamp as isize - c0 as isize).max(0) as f64 / (c1 - c0) as f64
            } else {
                0.0
            };

            // Weighted sum over the four corners, using only valid ones.
            let weights = [
                (1.0 - dr) * (1.0 - dc), // (r0, c0)
                (1.0 - dr) * dc,           // (r0, c1)
                dr * (1.0 - dc),           // (r1, c0)
                dr * dc,                   // (r1, c1)
            ];

            let mut sum_aero = 0.0;
            let mut sum_eps = 0.0;
            let mut sum_w = 0.0;

            for (k, &(r, c)) in corners.iter().enumerate() {
                if ipflag[r * nsamps + c] & flag_bit != 0 {
                    let w = weights[k];
                    sum_aero += w * taero[r * nsamps + c];
                    sum_eps += w * teps[r * nsamps + c];
                    sum_w += w;
                }
            }

            if sum_w > 0.0 {
                taero[pix] = sum_aero / sum_w;
                teps[pix] = sum_eps / sum_w;
                ipflag[pix] |= 1u8 << IPFLAG_INTERP_WINDOW;
            }
        }
    }
}

/// Three-pass local averaging to fix pixels with invalid aerosol retrievals.
///
/// Ported from `fix_invalid_aerosols()` in C LaSRC code.
///
/// - Pass 1 (forward):  require `min_clear_pix` valid (IPFLAG_CLEAR) neighbors.
/// - Pass 2 (forward):  require 1 valid pixel; also allow IPFLAG_FIXED pixels.
/// - Pass 3 (reverse):  same as pass 2 but iterate in reverse order.
///
/// Fixed pixels are marked with [`IPFLAG_FIXED`].
///
/// # Arguments
/// * `taero`                - Aerosol optical thickness (modified in place)
/// * `teps`                 - Angstrom exponent (modified in place)
/// * `ipflag`               - Processing flags (modified in place)
/// * `nlines`               - Number of image lines
/// * `nsamps`               - Number of image samples
/// * `aero_window`          - Aerosol retrieval window size
/// * `half_aero_window`     - Half of aerosol window
/// * `fix_aero_window`      - Search window for fixing invalid retrievals
/// * `half_fix_aero_window` - Half of fix window
/// * `min_clear_pix`        - Minimum number of clear neighbors required in pass 1
pub fn fix_invalid_aerosols(
    taero: &mut [f64],
    teps: &mut [f64],
    ipflag: &mut [u8],
    nlines: usize,
    nsamps: usize,
    _aero_window: usize,
    _half_aero_window: usize,
    fix_aero_window: usize,
    half_fix_aero_window: usize,
    min_clear_pix: usize,
) {
    let clear_bit = 1u8 << IPFLAG_CLEAR;
    let fixed_bit = 1u8 << IPFLAG_FIXED;
    let npix = nlines * nsamps;

    // ---- Pass 1: forward, require min_clear_pix IPFLAG_CLEAR neighbors ----
    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;
            if ipflag[pix] & clear_bit != 0 || ipflag[pix] & fixed_bit != 0 {
                continue; // already valid
            }

            let row_lo = iline.saturating_sub(half_fix_aero_window);
            let row_hi = (iline + half_fix_aero_window + 1).min(nlines);
            let col_lo = isamp.saturating_sub(half_fix_aero_window);
            let col_hi = (isamp + half_fix_aero_window + 1).min(nsamps);

            let mut sum_aero = 0.0;
            let mut sum_eps_val = 0.0;
            let mut count = 0usize;

            for r in row_lo..row_hi {
                for c in col_lo..col_hi {
                    let n = r * nsamps + c;
                    if ipflag[n] & clear_bit != 0 {
                        sum_aero += taero[n];
                        sum_eps_val += teps[n];
                        count += 1;
                    }
                }
            }

            if count >= min_clear_pix {
                taero[pix] = sum_aero / count as f64;
                teps[pix] = sum_eps_val / count as f64;
                ipflag[pix] |= fixed_bit;
            }
        }
    }

    // ---- Pass 2: forward, require 1 pixel (CLEAR or FIXED) ----
    for iline in 0..nlines {
        for isamp in 0..nsamps {
            let pix = iline * nsamps + isamp;
            if ipflag[pix] & clear_bit != 0 || ipflag[pix] & fixed_bit != 0 {
                continue;
            }

            let row_lo = iline.saturating_sub(half_fix_aero_window);
            let row_hi = (iline + half_fix_aero_window + 1).min(nlines);
            let col_lo = isamp.saturating_sub(half_fix_aero_window);
            let col_hi = (isamp + half_fix_aero_window + 1).min(nsamps);

            let mut sum_aero = 0.0;
            let mut sum_eps_val = 0.0;
            let mut count = 0usize;

            for r in row_lo..row_hi {
                for c in col_lo..col_hi {
                    let n = r * nsamps + c;
                    if ipflag[n] & clear_bit != 0 || ipflag[n] & fixed_bit != 0 {
                        sum_aero += taero[n];
                        sum_eps_val += teps[n];
                        count += 1;
                    }
                }
            }

            if count >= 1 {
                taero[pix] = sum_aero / count as f64;
                teps[pix] = sum_eps_val / count as f64;
                ipflag[pix] |= fixed_bit;
            }
        }
    }

    // ---- Pass 3: reverse, same as pass 2 ----
    for idx in (0..npix).rev() {
        let pix = idx;
        if ipflag[pix] & clear_bit != 0 || ipflag[pix] & fixed_bit != 0 {
            continue;
        }

        let iline = pix / nsamps;
        let isamp = pix % nsamps;

        let row_lo = iline.saturating_sub(half_fix_aero_window);
        let row_hi = (iline + half_fix_aero_window + 1).min(nlines);
        let col_lo = isamp.saturating_sub(half_fix_aero_window);
        let col_hi = (isamp + half_fix_aero_window + 1).min(nsamps);

        let mut sum_aero = 0.0;
        let mut sum_eps_val = 0.0;
        let mut count = 0usize;

        for r in row_lo..row_hi {
            for c in col_lo..col_hi {
                let n = r * nsamps + c;
                if ipflag[n] & clear_bit != 0 || ipflag[n] & fixed_bit != 0 {
                    sum_aero += taero[n];
                    sum_eps_val += teps[n];
                    count += 1;
                }
            }
        }

        if count >= 1 {
            taero[pix] = sum_aero / count as f64;
            teps[pix] = sum_eps_val / count as f64;
            ipflag[pix] |= fixed_bit;
        }
    }

    // Suppress unused-parameter warnings for window params not needed in the
    // inner loops (the C code uses fix_aero_window only for bounds clamping
    // which we do via saturating_sub / min).
    let _ = fix_aero_window;
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
        let mut taero = vec![0.0f64; npix];
        let mut teps = vec![0.0f64; npix];
        let mut ipflag = vec![0u8; npix];
        let qa_band = vec![1u16; npix];
        for i in (1..nlines).step_by(3) {
            for j in (1..nsamps).step_by(3) {
                let pix = i * nsamps + j;
                taero[pix] = 0.1;
                teps[pix] = 1.5;
                ipflag[pix] = 1 << IPFLAG_CLEAR;
            }
        }
        aerosol_interp(&mut taero, &mut teps, &mut ipflag, &qa_band, nlines, nsamps, 3, 1);
        // Check that non-center pixels got interpolated values
        let non_center_pix = 0 * nsamps + 0; // pixel (0,0)
        assert!(taero[non_center_pix] > 0.0, "Pixel (0,0) should have interpolated aerosol");
    }

    #[test]
    fn test_fix_invalid_aerosols_fills() {
        let nlines = 9;
        let nsamps = 9;
        let npix = nlines * nsamps;
        let mut taero = vec![0.0f64; npix];
        let mut teps = vec![0.0f64; npix];
        let mut ipflag = vec![0u8; npix];
        // Set some valid centers
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

        fix_invalid_aerosols(&mut taero, &mut teps, &mut ipflag, nlines, nsamps, 3, 1, 15, 7, 4);
        assert!(ipflag[bad_pix] & (1 << IPFLAG_FIXED) != 0, "Should be fixed");
        assert!(taero[bad_pix] > 0.0, "Should have filled value");
    }
}
