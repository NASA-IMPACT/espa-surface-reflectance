//! Differential test: C `get_3rd_order_poly_coeff` vs
//! `lasrc_core::utils::get_3rd_order_poly_coeff`.
//!
//! Cubic fit (LU-solved normal equations, f32 throughout in C). Feeds the
//! roatm/ttatmg/satm poly coefficients consumed by atmcorlamb2_new /
//! subaeroret_new, so a divergence here would cascade even though those
//! kernels are proven faithful. aot = the 22-point AOT550NM grid (distinct ->
//! well-conditioned design matrix); atm is fuzzed.
//!
//! Tolerance (not bit-exact): the LU accumulation is a long f32 op chain where
//! a faithful port can still differ by a last ULP; a gross mismatch (wrong
//! order/formula/stride) blows past this.
#![cfg(feature = "lut")]

use lasrc_core::utils::get_3rd_order_poly_coeff as rust_poly;
use lasrc_cref::{c_get_3rd_order_poly_coeff, AOT550NM};
use proptest::prelude::*;

proptest! {
    #[test]
    fn matches_c(atm in proptest::collection::vec(-1.0f32..1.0, AOT550NM.len())) {
        let aot = AOT550NM.to_vec();

        let c = c_get_3rd_order_poly_coeff(&aot, &atm);

        let aot64: Vec<f64> = aot.iter().map(|&x| x as f64).collect();
        let atm64: Vec<f64> = atm.iter().map(|&x| x as f64).collect();
        let r = rust_poly(&aot64, &atm64);

        for k in 0..4 {
            let (rk, ck) = (r[k] as f32, c[k]);
            prop_assume!(rk.is_finite() && ck.is_finite());
            let tol = 1e-4 * ck.abs().max(rk.abs()).max(1e-3);
            prop_assert!(
                (rk - ck).abs() <= tol,
                "coeff[{}] differ: rust={} c={} diff={}", k, rk, ck, (rk - ck).abs()
            );
        }
    }
}
