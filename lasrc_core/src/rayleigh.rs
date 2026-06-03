/// Compute Rayleigh scattering reflectance using the depolarization factor
/// and a 10-coefficient Legendre polynomial approximation.
///
/// Ported from `local_chand()` in C `lut_subr.c`.
/// C uses float variables throughout, but cos/exp/sqrt/log are double (math.h).
/// We use f32 for local variables and f64 for transcendental functions,
/// matching C's implicit float->double promotion and double->float truncation.
///
/// # Arguments
/// * `xphi` - Relative azimuth angle (degrees)
/// * `xmuv` - Cosine of view zenith angle
/// * `xmus` - Cosine of solar zenith angle
/// * `xtau` - Rayleigh optical depth
pub fn rayleigh_reflectance(xphi: f64, xmuv: f64, xmus: f64, xtau: f64) -> f64 {
    let xphi = xphi as f32;
    let xmuv = xmuv as f32;
    let xmus = xmus as f32;
    let xtau = xtau as f32;

    let xfd: f32 = 0.958_725_777;

    let xmus2: f32 = xmus * xmus;
    let xmuv2: f32 = xmuv * xmuv;

    // Phase function components
    let xph1: f32 = 1.0 + (3.0 * xmus2 - 1.0) * (3.0 * xmuv2 - 1.0) * xfd * 0.125;
    let xph3_base: f32 = (1.0 - xmus2) * (1.0 - xmuv2);
    // C: xph2 = -xmus * xmuv * sqrt(xph3); — sqrt is double
    let xph2: f32 = -xmus * xmuv * (xph3_base as f64).sqrt() as f32;
    let xph2: f32 = xph2 * xfd * 0.75;
    let xph3: f32 = xph3_base * xfd * 0.1875;

    // C: phios = xphi * DEG2RAD; (float)
    // C: xcosf2 = -cos(phios); xcosf3 = cos(2.0 * phios); — cos is double
    let deg2rad: f32 = std::f32::consts::PI / 180.0;
    let phios: f32 = xphi * deg2rad;
    let xcosf2: f32 = -(phios as f64).cos() as f32;
    let xcosf3: f32 = (2.0f64 * phios as f64).cos() as f32;

    // C: xitm = (1.0 - exp(-xtau * (1.0/xmus + 1.0/xmuv))) / (4*(xmus+xmuv));
    // Note: 1.0 is double, exp is double, result stored as float
    let xitm_ss: f32 = ((1.0 - (-(xtau as f64) * (1.0 / xmus as f64 + 1.0 / xmuv as f64)).exp())
        / (4.0 * (xmus as f64 + xmuv as f64))) as f32;
    let xp1: f32 = xph1 * xitm_ss;
    let xp2: f32 = xph2 * xitm_ss;
    let xp3: f32 = xph3 * xitm_ss;

    // C: xitm = (1.0 - exp(-xtau/xmus)) * (1.0 - exp(-xtau/xmuv));
    let xitm_ms: f32 = ((1.0 - (-(xtau as f64) / xmus as f64).exp())
        * (1.0 - (-(xtau as f64) / xmuv as f64).exp())) as f32;
    let cfonc1: f32 = xph1 * xitm_ms;
    let cfonc2: f32 = xph2 * xitm_ms;
    let cfonc3: f32 = xph3 * xitm_ms;

    // C: xlntau = log(xtau); — log is double, xlntau is float
    let xlntau: f32 = (xtau as f64).ln() as f32;

    let as0: [f32; 10] = [
        0.33243832, -6.777104e-02, 0.16285370, 1.577425e-03,
        -0.30924818, -1.240906e-02, -0.10324388, 3.241678e-02,
        0.11493334, -3.503695e-02,
    ];
    let as1: [f32; 2] = [0.19666292, -5.439061e-02];
    let as2: [f32; 2] = [0.14545937, -2.910845e-02];

    let pl: [f32; 10] = [
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

    let mut fs0: f32 = 0.0;
    for i in 0..10 {
        fs0 += pl[i] * as0[i];
    }
    let fs1: f32 = pl[0] * as1[0] + pl[1] * as1[1];
    let fs2: f32 = pl[0] * as2[0] + pl[1] * as2[1];

    let xitot1: f32 = xp1 + cfonc1 * fs0;
    let xitot2: f32 = xp2 + cfonc2 * fs1;
    let xitot3: f32 = xp3 + cfonc3 * fs2;

    let result: f32 = xitot1 + 2.0 * (xitot2 * xcosf2 + xitot3 * xcosf3);
    result as f64
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
