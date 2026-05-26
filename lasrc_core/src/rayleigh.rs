use std::f64::consts::PI;

/// Compute Rayleigh scattering reflectance using the depolarization factor
/// and a 10-coefficient Legendre polynomial approximation.
///
/// Ported from `local_chand()` in C `lut_subr.c`.
///
/// # Arguments
/// * `xphi` - Relative azimuth angle (degrees)
/// * `xmuv` - Cosine of view zenith angle
/// * `xmus` - Cosine of solar zenith angle
/// * `xtau` - Rayleigh optical depth
pub fn rayleigh_reflectance(xphi: f64, xmuv: f64, xmus: f64, xtau: f64) -> f64 {
    // Depolarization factor derived from xdep = 0.0279
    // xfd = xdep / (2.0 - xdep) = 0.014147355
    // xfd = (1.0 - xfd) / (1.0 + 2.0 * xfd) = 0.958725777
    let xfd = 0.958_725_777_f64;

    let xmus2 = xmus * xmus;
    let xmuv2 = xmuv * xmuv;

    // Phase function components
    let xph1 = 1.0 + (3.0 * xmus2 - 1.0) * (3.0 * xmuv2 - 1.0) * xfd * 0.125;
    let xph3_base = (1.0 - xmus2) * (1.0 - xmuv2);
    let xph2 = -xmus * xmuv * xph3_base.sqrt() * xfd * 0.75;
    let xph3 = xph3_base * xfd * 0.1875;

    // Azimuth terms — note xcosf2 uses negated cosine (single phi) per C source
    let phios = xphi * PI / 180.0;
    let xcosf2 = -phios.cos();
    let xcosf3 = (2.0 * phios).cos();

    // Single-scattering term: removes xmus factor to avoid extra division
    let xitm_ss = (1.0 - (-xtau * (1.0 / xmus + 1.0 / xmuv)).exp()) / (4.0 * (xmus + xmuv));
    let xp1 = xph1 * xitm_ss;
    let xp2 = xph2 * xitm_ss;
    let xp3 = xph3 * xitm_ss;

    // Multiple-scattering correction term
    let xitm_ms = (1.0 - (-xtau / xmus).exp()) * (1.0 - (-xtau / xmuv).exp());
    let cfonc1 = xph1 * xitm_ms;
    let cfonc2 = xph2 * xitm_ms;
    let cfonc3 = xph3 * xitm_ms;

    // Legendre polynomial coefficients
    let as0: [f64; 10] = [
        0.33243832,
        -6.777104e-02,
        0.16285370,
        1.577425e-03,
        -0.30924818,
        -1.240906e-02,
        -0.10324388,
        3.241678e-02,
        0.11493334,
        -3.503695e-02,
    ];
    let as1: [f64; 2] = [0.19666292, -5.439061e-02];
    let as2: [f64; 2] = [0.14545937, -2.910845e-02];

    // Polynomial basis terms
    let xlntau = xtau.ln();
    let pl = [
        1.0,
        xlntau,
        xmus + xmuv,
        xlntau * (xmus + xmuv),
        xmus * xmuv,
        xlntau * xmus * xmuv,
        xmus2 + xmuv2,
        xlntau * (xmus2 + xmuv2),
        xmus2 * xmuv2,
        xlntau * xmus2 * xmuv2,
    ];

    let fs0: f64 = as0.iter().zip(pl.iter()).map(|(a, p)| a * p).sum();
    let fs1 = pl[0] * as1[0] + pl[1] * as1[1];
    let fs2 = pl[0] * as2[0] + pl[1] * as2[1];

    let xitot1 = xp1 + cfonc1 * fs0;
    let xitot2 = xp2 + cfonc2 * fs1;
    let xitot3 = xp3 + cfonc3 * fs2;

    xitot1 + 2.0 * (xitot2 * xcosf2 + xitot3 * xcosf3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rayleigh_positive() {
        let rray = rayleigh_reflectance(90.0, 0.95, 0.5, 0.17);
        assert!(rray > 0.0, "Rayleigh reflectance should be positive: {rray}");
    }

    #[test]
    fn test_rayleigh_increases_with_tau() {
        let r1 = rayleigh_reflectance(90.0, 0.95, 0.5, 0.05);
        let r2 = rayleigh_reflectance(90.0, 0.95, 0.5, 0.25);
        assert!(r2 > r1, "More optical depth = more scattering: r1={r1}, r2={r2}");
    }

    #[test]
    fn test_rayleigh_range() {
        let rray = rayleigh_reflectance(45.0, 0.7, 0.6, 0.24);
        assert!(rray > 0.0 && rray < 0.5, "rray={rray}");
    }
}
