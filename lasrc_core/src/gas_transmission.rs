/// Coefficients for computing gas transmissions for a given band.
pub struct GasCoefficients {
    pub ogtransa1: f64,
    pub ogtransb0: f64,
    pub ogtransb1: f64,
    pub wvtransa: f64,
    pub wvtransb: f64,
    pub oztransa: f64,
}

/// Gas transmission values for a given atmospheric state.
pub struct GasTransmission {
    /// Ozone transmission
    pub tgoz: f32,
    /// Water vapor transmission (full path)
    pub tgwv: f32,
    /// Water vapor transmission (solar half-path)
    pub tgwvhalf: f32,
    /// Other gas transmission
    pub tgog: f32,
    /// Combined gas transmission (tgog * tgoz)
    pub tgo: f32,
}

/// Compute gas transmissions for a given band and atmospheric state.
///
/// Ported from `comptg()` in C `lut_subr.c`. Inputs and outputs are float like
/// C; the gas coefficients are double in C, so expressions using them (and the
/// math.h calls) are evaluated in double before being stored as float.
///
/// # Arguments
/// * `coeff` - Band-specific gas coefficients
/// * `xmus` - Cosine of solar zenith angle
/// * `xmuv` - Cosine of view zenith angle
/// * `uoz` - Ozone amount (cm-atm)
/// * `uwv` - Water vapor amount (g/cm^2)
/// * `atm_pres` - Normalized atmospheric pressure (pres/1013)
pub fn compute_gas_transmission(
    coeff: &GasCoefficients,
    xmus: f32,
    xmuv: f32,
    uoz: f32,
    uwv: f32,
    atm_pres: f32,
) -> GasTransmission {
    // C: float m = 1.0 / xmus + 1.0 / xmuv;  (1.0 is double)
    let m: f32 = (1.0 / xmus as f64 + 1.0 / xmuv as f64) as f32;

    // C: *tgoz = exp(oztransa[iband] * m * uoz);  (oztransa is double)
    let tgoz = (coeff.oztransa * m as f64 * uoz as f64).exp() as f32;

    // C: float x = m * uwv; float a = wvtransa; float b = wvtransb;
    let x: f32 = m * uwv;
    let a: f32 = coeff.wvtransa as f32;
    let b: f32 = coeff.wvtransb as f32;

    // C: *tgwv = exp(-a * pow(x, b));
    let tgwv = if x as f64 > 1e-06 {
        (-(a as f64) * (x as f64).powf(b as f64)).exp() as f32
    } else {
        1.0
    };

    // C: x *= 0.5; *tgwvhalf = exp(-a * pow(x, b));
    let xhalf: f32 = (x as f64 * 0.5) as f32;
    let tgwvhalf = if xhalf as f64 > 1e-06 {
        (-(a as f64) * (xhalf as f64).powf(b as f64)).exp() as f32
    } else {
        1.0
    };

    // C: *tgog = -(ogtransa1[iband] * atm_pres) *
    //        pow(m, exp(-(ogtransb0[iband] + ogtransb1[iband] * atm_pres)));
    //    *tgog = exp(*tgog);
    let exponent = (-(coeff.ogtransb0 + coeff.ogtransb1 * atm_pres as f64)).exp();
    let tgog_arg = (-(coeff.ogtransa1 * atm_pres as f64) * (m as f64).powf(exponent)) as f32;
    let tgog = (tgog_arg as f64).exp() as f32;

    // C: *tgo = tgog * tgoz;
    let tgo = tgog * tgoz;

    GasTransmission {
        tgoz,
        tgwv,
        tgwvhalf,
        tgog,
        tgo,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_coeff() -> GasCoefficients {
        GasCoefficients {
            ogtransa1: 0.00014,
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
        assert!(
            (gt.tgwv - 1.0).abs() < 1e-6,
            "No water vapor: tgwv should be ~1.0"
        );
        assert!((gt.tgwvhalf - 1.0).abs() < 1e-6);
    }
}
