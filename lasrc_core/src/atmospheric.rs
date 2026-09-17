//! Atmospheric correction routines ported from `atmcorlamb2` and `atmcorlamb2_new`
//! in C `lut_subr.c`.

use crate::constants::*;
use crate::gas_transmission::{GasCoefficients, compute_gas_transmission};
use crate::lut::{LookupTables, LutIndices};
use crate::rayleigh::rayleigh_reflectance;

/// Result of a full atmospheric correction.
#[derive(Debug, Clone)]
pub struct AtmCorrResult {
    /// Surface reflectance
    pub roslamb: f64,
    /// Other gas transmittance
    pub tgo: f64,
    /// Atmospheric intrinsic reflectance (corrected)
    pub roatm: f64,
    /// Total atmospheric transmittance
    pub ttatmg: f64,
    /// Spherical albedo
    pub satm: f64,
    /// Rayleigh reflectance
    pub xrorayp: f64,
}

/// Pre-fitted polynomial coefficients for the fast atmospheric correction path.
#[derive(Debug, Clone)]
pub struct AtmCorrCoefficients {
    /// Upper bound for the scaled AOT
    pub roatm_upper: f64,
    /// Cubic polynomial coefficients for atmospheric reflectance [a3, a2, a1, a0]
    pub roatm_coef: [f64; NCOEF],
    /// Cubic polynomial coefficients for total transmittance [a3, a2, a1, a0]
    pub ttatmg_coef: [f64; NCOEF],
    /// Cubic polynomial coefficients for spherical albedo [a3, a2, a1, a0]
    pub satm_coef: [f64; NCOEF],
}

/// Evaluate a cubic polynomial using f32 arithmetic to match C's float precision.
///
/// C computes: `coef[3] + coef[2]*x + coef[1]*x^2 + coef[0]*x^3` using `float`.
/// We cast to f32 for the evaluation, then return as f64.
#[inline]
fn eval_cubic(c: &[f64; NCOEF], x: f64) -> f64 {
    let xf = x as f32;
    let xf_sq = xf * xf;
    let xf_cube = xf_sq * xf;
    let result = c[3] as f32
        + c[2] as f32 * xf
        + c[1] as f32 * xf_sq
        + c[0] as f32 * xf_cube;
    result as f64
}

/// Fast polynomial atmospheric correction (pre-fitted coefficients path).
///
/// Ported from `atmcorlamb2_new()` in C `lut_subr.c`.
///
/// # Arguments
/// * `coeff`            - Pre-fitted polynomial coefficients for this band/geometry
/// * `tgo`              - Other gas (ozone * other) transmittance
/// * `iband`            - Band index into `lambda`
/// * `raot550nm`        - Aerosol optical thickness at 550 nm
/// * `normext_ib_0_3`   - Normalised aerosol extinction for this band at the reference AOT
/// * `rotoa`            - TOA reflectance
/// * `lambda`           - Wavelength table (µm), one entry per band
/// * `eps`              - Angstrom exponent
///
/// # Returns
/// Surface reflectance `roslamb`.
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
    // C uses float throughout atmcorlamb2_new; match that precision.
    // C: static const float lambda_sf = 1/0.55; -- the division is done in
    // double and only the result is stored as float (1 ULP above the f32
    // division), which pow() then amplifies.
    let lambda_sf: f32 = (1.0f64 / 0.55f64) as f32;

    let mut mraot550nm: f32 = if eps < 0.0 || iband >= lambda.len() {
        raot550nm as f32
    } else {
        // C: pow() uses double precision even though operands start as float.
        // The float product (lambda * lambda_sf) is promoted to double for pow(),
        // then the result (double) is multiplied with the float quotient (promoted
        // to double), and the final result truncated back to float.
        let base = (lambda[iband] as f32 * lambda_sf) as f64;
        let power = base.powf(-(eps as f64));
        ((raot550nm as f32 / normext_ib_0_3 as f32) as f64 * power) as f32
    };

    if mraot550nm >= coeff.roatm_upper as f32 {
        mraot550nm = coeff.roatm_upper as f32;
    }

    // Evaluate cubic polynomials (eval_cubic already uses f32 internally)
    let roatm = eval_cubic(&coeff.roatm_coef, mraot550nm as f64);
    let ttatmg = eval_cubic(&coeff.ttatmg_coef, mraot550nm as f64);
    let satm = eval_cubic(&coeff.satm_coef, mraot550nm as f64);

    // Atmospheric correction — must match C's formula exactly:
    //   *roslamb = rotoa - tgo*roatm;
    //   *roslamb /= tgo*ttatmg + satm*(*roslamb);
    // Note: C's atmcorlamb2_new uses a DIFFERENT formula from atmcorlamb2.
    // The fast path multiplies by tgo instead of dividing.
    let roslamb_f = rotoa as f32 - (tgo as f32) * (roatm as f32);
    let result = roslamb_f / ((tgo as f32) * (ttatmg as f32) + (satm as f32) * roslamb_f);
    result as f64
}

/// Full LUT-based atmospheric correction.
///
/// Ported from `atmcorlamb2()` in C `lut_subr.c`.
///
/// # Arguments
/// * `lut`          - Pre-loaded lookup tables
/// * `gas_coeff`    - Band-specific gas transmission coefficients
/// * `tauray_band`  - Rayleigh optical depth for this band
/// * `iband`        - Band index
/// * `xts`          - Solar zenith angle (degrees)
/// * `xtv`          - View zenith angle (degrees)
/// * `xmus`         - Cosine of solar zenith
/// * `xmuv`         - Cosine of view zenith
/// * `xfi`          - Relative azimuth angle (degrees)
/// * `cosxfi`       - Cosine of relative azimuth
/// * `raot550nm`    - Aerosol optical thickness at 550 nm
/// * `pressure`     - Surface pressure (mbar)
/// * `uoz`          - Ozone amount (cm-atm)
/// * `uwv`          - Water vapour amount (g/cm²)
/// * `rotoa`        - TOA reflectance
/// * `lambda`       - Wavelength table (µm), one entry per band
/// * `max_band_idx` - Maximum valid band index for wavelength scaling
/// * `eps`          - Angstrom exponent
///
/// # Returns
/// [`AtmCorrResult`] with surface reflectance and intermediate quantities.
#[allow(clippy::too_many_arguments)]
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
    lambda: &[f64],
    max_band_idx: usize,
    eps: f64,
) -> AtmCorrResult {
    // C uses float (f32) throughout atmcorlamb2. Truncate all intermediate
    // values to f32 after calling f64 helper functions to match C precision.

    // C: all angle, pressure and gas amount inputs are float
    let (xts, xtv, xmus, xmuv, xfi) = (xts as f32, xtv as f32, xmus as f32, xmuv as f32, xfi as f32);
    let (pres, uoz, uwv) = (pressure as f32, uoz as f32, uwv as f32);

    // Modify AOT based on Angstrom coefficient and wavelength.
    // C (atmcorlamb2, unlike atmcorlamb2_new): static const double lambda_sf = 1/0.55;
    let lambda_sf: f64 = 1.0 / 0.55;
    let mraot550nm: f32 = if eps < 0.0 || iband > max_band_idx {
        raot550nm as f32
    } else {
        let normext_idx = iband * NPRES_VALS * NAOT_VALS + 3;
        let normext_val: f32 = if normext_idx < lut.normext.len() {
            lut.normext[normext_idx] as f32
        } else {
            1.0f32
        };
        // C: mraot550nm = (raot550nm / normext[indx]) * (pow((lambda[iband] * lambda_sf), -eps));
        // float / float, then * double pow(float * double, -float)
        let base = lambda[iband] as f32 as f64 * lambda_sf;
        let power = base.powf(-eps);
        ((raot550nm as f32 / normext_val) as f64 * power) as f32
    };

    // C: atm_pres = pres * ONE_DIV_ATMOS_PRES_0;  (float * double literal)
    let atm_pres: f32 = (pres as f64 * ONE_DIV_ATMOS_PRES_0) as f32;

    // Rayleigh optical depth scaled to surface pressure
    let xtaur: f32 = (tauray_band as f32) * atm_pres;

    let xrorayp: f32 =
        rayleigh_reflectance(xfi as f64, xmuv as f64, xmus as f64, xtaur as f64) as f32;

    // Find pressure and AOT indices into the LUT using modified AOT
    let indices: LutIndices = lut.find_indices(pres as f64, mraot550nm as f64);

    let satm: f32 = lut.interp_spherical_albedo(&indices, iband, pres as f64, mraot550nm as f64) as f32;
    let roatm_raw: f32 = lut.interp_atmospheric_reflectance(
        &indices, iband, pres as f64, mraot550nm as f64, xts as f64, xtv as f64,
        xmus as f64, xmuv as f64, cosxfi,
    ) as f32;

    // Downward and upward transmittances. C calls comptrans with the solar
    // (xtsmin, xtsstep) grid for xts but the view (xtvmin, xtvstep) grid for xtv,
    // while indexing the same tts table.
    let xtts: f32 = lut.interp_transmission(
        &indices, iband, pres as f64, mraot550nm as f64, xts as f64, XTS_MIN as f32, XTS_STEP as f32,
    ) as f32;
    let xttv: f32 = lut.interp_transmission(
        &indices, iband, pres as f64, mraot550nm as f64, xtv as f64, XTV_MIN as f32, XTV_STEP_C,
    ) as f32;
    let ttatm: f32 = xtts * xttv;

    let gt = compute_gas_transmission(gas_coeff, xmus, xmuv, uoz, uwv, atm_pres);
    let tgo: f32 = gt.tgo;
    let tgwv: f32 = gt.tgwv;
    let tgwvhalf: f32 = gt.tgwvhalf;

    // Apply water-vapour half-path correction to atmospheric reflectance
    let roatm_corrected: f32 = (roatm_raw - xrorayp) * tgwvhalf + xrorayp;

    // Total transmittance including gas absorption
    let ttatmg: f32 = ttatm * tgwv;

    // Solve for surface reflectance
    let xroslamb: f32 = rotoa as f32 / tgo - roatm_corrected;
    let roslamb: f32 = xroslamb / (ttatmg + satm * xroslamb);

    AtmCorrResult {
        roslamb: roslamb as f64,
        tgo: tgo as f64,
        roatm: roatm_corrected as f64,
        ttatmg: ttatmg as f64,
        satm: satm as f64,
        xrorayp: xrorayp as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atmcorlamb2_new_low_aot() {
        let coeff = AtmCorrCoefficients {
            roatm_upper: 0.4,
            roatm_coef: [0.0, 0.0, 0.1, 0.01],
            ttatmg_coef: [0.0, 0.0, -0.1, 0.9],
            satm_coef: [0.0, 0.0, 0.05, 0.01],
        };
        let roslamb = atmcorlamb2_new(
            &coeff, 0.95, 0, 0.05, 1.0, 0.15, &LAMBDA_LANDSAT, 1.5,
        );
        // With tgo=0.95 the gas-absorption correction inflates the result slightly
        // above rotoa=0.15; allow up to 0.20 while still bounding it to a physically
        // reasonable small positive range.
        assert!(roslamb < 0.20, "roslamb={roslamb}");
        assert!(roslamb > -0.2, "roslamb={roslamb}");
    }

    #[test]
    fn test_atmcorlamb2_new_higher_aot_differs() {
        let coeff = AtmCorrCoefficients {
            roatm_upper: 0.5,
            roatm_coef: [0.0, 0.0, 0.2, 0.01],
            ttatmg_coef: [0.0, 0.0, -0.15, 0.9],
            satm_coef: [0.0, 0.0, 0.08, 0.01],
        };
        let sr1 = atmcorlamb2_new(&coeff, 0.95, 0, 0.05, 1.0, 0.15, &LAMBDA_LANDSAT, 1.5);
        let sr2 = atmcorlamb2_new(&coeff, 0.95, 0, 0.50, 1.0, 0.15, &LAMBDA_LANDSAT, 1.5);
        assert!((sr1 - sr2).abs() > 0.001, "Different AOT should give different SR");
    }
}
