//! Differential test: C `local_chand` vs `lasrc_core::rayleigh::rayleigh_reflectance`.
//!
//! Molecular (Rayleigh) reflectance. Uses cos/exp internally (libm both sides),
//! and the Rust port runs f64 vs C float, so assert a tight tolerance rather
//! than bit-exact. Feeds the atmospheric correction's Rayleigh term.
#![cfg(feature = "lut")]

use lasrc_core::rayleigh::rayleigh_reflectance as rust_chand;
use lasrc_cref::c_local_chand;
use proptest::prelude::*;

proptest! {
    #[test]
    fn matches_c(
        xphi in 0.0f32..180.0,   // relative azimuth (deg)
        xmuv in 0.05f32..1.0,    // cos(view zenith), away from grazing
        xmus in 0.05f32..1.0,    // cos(solar zenith)
        xtau in 0.0f32..0.5,     // molecular optical depth
    ) {
        let c = c_local_chand(xphi, xmuv, xmus, xtau);
        let r = rust_chand(xphi as f64, xmuv as f64, xmus as f64, xtau as f64) as f32;

        prop_assume!(c.is_finite() && r.is_finite());
        let tol = 1e-5 * c.abs().max(1e-4);
        prop_assert!(
            (r - c).abs() <= tol,
            "rayleigh differ: rust={} c={} diff={} (xphi={} xmuv={} xmus={} xtau={})",
            r, c, (r - c).abs(), xphi, xmuv, xmus, xtau
        );
    }
}
