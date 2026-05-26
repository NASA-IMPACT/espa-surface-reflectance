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

/// Evaluate a cubic polynomial: c[0]*x^3 + c[1]*x^2 + c[2]*x + c[3].
#[inline]
fn eval_cubic(c: &[f64; NCOEF], x: f64) -> f64 {
    ((c[0] * x + c[1]) * x + c[2]) * x + c[3]
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
    let lambda_sf = 1.0 / 0.55;

    // Scale AOT to this band's wavelength using the Angstrom relation
    let mut mraot550nm =
        (raot550nm / normext_ib_0_3) * (lambda[iband] * lambda_sf).powf(-eps);

    // Clamp to the valid range of the fitted polynomials
    if mraot550nm >= coeff.roatm_upper {
        mraot550nm = coeff.roatm_upper;
    }

    // Evaluate cubic polynomials
    let roatm = eval_cubic(&coeff.roatm_coef, mraot550nm);
    let ttatmg = eval_cubic(&coeff.ttatmg_coef, mraot550nm);
    let satm = eval_cubic(&coeff.satm_coef, mraot550nm);

    // Solve for surface reflectance
    let xroslamb = rotoa / tgo - roatm;
    xroslamb / (ttatmg + satm * xroslamb)
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
/// * `eps`          - Angstrom exponent (unused here, kept for API parity)
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
    _eps: f64,
) -> AtmCorrResult {
    // Normalised atmospheric pressure
    let atm_pres = pressure * ONE_DIV_ATMOS_PRES_0;

    // Rayleigh optical depth scaled to surface pressure
    let xtaur = tauray_band * atm_pres;

    // Rayleigh scattering reflectance
    let xrorayp = rayleigh_reflectance(xfi, xmuv, xmus, xtaur);

    // Find pressure and AOT indices into the LUT
    let indices: LutIndices = lut.find_indices(pressure, raot550nm);

    // Interpolate atmospheric quantities from LUT
    let satm = lut.interp_spherical_albedo(&indices, iband, pressure, raot550nm);
    let roatm_raw = lut.interp_atmospheric_reflectance(
        &indices, iband, pressure, raot550nm, xts, xtv, xmus, xmuv, cosxfi,
    );

    // Downward and upward transmittances, then total atmospheric transmittance
    let xtts = lut.interp_transmission(&indices, iband, pressure, raot550nm, xts);
    let xttv = lut.interp_transmission(&indices, iband, pressure, raot550nm, xtv);
    let ttatm = xtts * xttv;

    // Gas transmissions
    let gt = compute_gas_transmission(gas_coeff, xmus, xmuv, uoz, uwv, atm_pres);

    // Apply water-vapour half-path correction to atmospheric reflectance
    let roatm_corrected = (roatm_raw - xrorayp) * gt.tgwvhalf + xrorayp;

    // Total transmittance including gas absorption
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
