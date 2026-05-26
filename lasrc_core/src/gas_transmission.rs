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
    pub tgoz: f64,
    /// Water vapor transmission (full path)
    pub tgwv: f64,
    /// Water vapor transmission (solar half-path)
    pub tgwvhalf: f64,
    /// Other gas transmission
    pub tgog: f64,
    /// Combined gas transmission (tgog * tgoz)
    pub tgo: f64,
}

/// Compute gas transmissions for a given band and atmospheric state.
///
/// Ported from `comptg()` in C `lut_subr.c`.
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
    xmus: f64,
    xmuv: f64,
    uoz: f64,
    uwv: f64,
    atm_pres: f64,
) -> GasTransmission {
    let m = 1.0 / xmus + 1.0 / xmuv;

    // Ozone transmission
    let tgoz = (coeff.oztransa * m * uoz).exp();

    // Water vapor transmission (full path)
    let x = m * uwv;
    let tgwv = if x > 1e-06 {
        (-coeff.wvtransa * x.powf(coeff.wvtransb)).exp()
    } else {
        1.0
    };

    // Water vapor transmission (solar half-path only)
    let xhalf = uwv / xmus;
    let tgwvhalf = if xhalf > 1e-06 {
        (-coeff.wvtransa * xhalf.powf(coeff.wvtransb)).exp()
    } else {
        1.0
    };

    // Other gas transmission
    let exponent = (-(coeff.ogtransb0 + coeff.ogtransb1 * atm_pres)).exp();
    let tgog = (-(coeff.ogtransa1 * atm_pres) * m.powf(exponent)).exp();

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
