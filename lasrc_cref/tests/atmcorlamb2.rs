//! Differential test: C `atmcorlamb2_new` vs `lasrc_core::atmospheric::atmcorlamb2_new`.
//!
//! The polynomial fast-path correction. C picks its `lambda` table + max band
//! index from `Sat_t`; the Rust port takes `lambda` explicitly. We pin
//! sat = Landsat 8 and pass LAMBDA_LANDSAT so both use the identical table
//! (verified equal, len 7 = NREFLL_BANDS).
//!
//! Two regimes:
//!  - eps < 0: modification branch skipped -> pure f32 arithmetic -> bit-exact.
//!  - eps >= 0: includes `pow()` (double, libm both sides) -> assert tight
//!    tolerance, since a last-ULP double difference can survive truncation.
//!
//! Only built with `--features lut` (needs the C leaf compiled against conda
//! headers).
#![cfg(feature = "lut")]

use lasrc_core::atmospheric::{atmcorlamb2_new as rust_atmcorlamb2_new, AtmCorrCoefficients};
use lasrc_core::constants::LAMBDA_LANDSAT;
use lasrc_cref::{c_atmcorlamb2_new, SAT_LANDSAT_8};
use proptest::prelude::*;

/// Physically-plausible input ranges keep the correction well-conditioned
/// (the formula divides by `tgo*ttatmg + satm*roslamb`, which random inputs
/// can drive to zero and blow up f32-vs-f64 rounding at the singularity).
#[derive(Debug, Clone)]
struct Inputs {
    tgo: f32,
    roatm_upper: f32,
    roatm_coef: [f32; 4],
    ttatmg_coef: [f32; 4],
    satm_coef: [f32; 4],
    raot550nm: f32,
    iband: i32,
    normext: f32,
    rotoa: f32,
    eps: f32,
}

fn coef(lo: f32, hi: f32) -> impl Strategy<Value = [f32; 4]> {
    // cubic..constant terms; leading (cubic) terms kept small like real fits.
    [(lo..hi), (lo..hi), (lo..hi), (0.0f32..1.0)].prop_map(|[a, b, c, d]| [a, b, c, d])
}

fn inputs(eps: impl Strategy<Value = f32>) -> impl Strategy<Value = Inputs> {
    (
        0.5f32..1.0,          // tgo
        0.05f32..0.5,         // roatm_upper
        coef(-0.05, 0.05),    // roatm_coef
        coef(0.0, 0.3),       // ttatmg_coef
        coef(-0.05, 0.05),    // satm_coef
        0.0f32..0.5,          // raot550nm
        0i32..=6,             // iband (0..NREFLL_BANDS-1)
        0.5f32..2.0,          // normext
        0.0f32..1.0,          // rotoa
        eps,
    )
        .prop_map(|(tgo, ru, rc, tc, sc, raot, iband, nx, rotoa, eps)| Inputs {
            tgo, roatm_upper: ru, roatm_coef: rc, ttatmg_coef: tc, satm_coef: sc,
            raot550nm: raot, iband, normext: nx, rotoa, eps,
        })
}

fn run_rust(i: &Inputs) -> f32 {
    let coeff = AtmCorrCoefficients {
        roatm_upper: i.roatm_upper as f64,
        roatm_coef: i.roatm_coef.map(|v| v as f64),
        ttatmg_coef: i.ttatmg_coef.map(|v| v as f64),
        satm_coef: i.satm_coef.map(|v| v as f64),
    };
    rust_atmcorlamb2_new(
        &coeff,
        i.tgo as f64,
        i.iband as usize,
        i.raot550nm as f64,
        i.normext as f64,
        i.rotoa as f64,
        &LAMBDA_LANDSAT,
        i.eps as f64,
    ) as f32
}

fn run_c(i: &Inputs) -> f32 {
    c_atmcorlamb2_new(
        SAT_LANDSAT_8, i.tgo, i.roatm_upper,
        &i.roatm_coef, &i.ttatmg_coef, &i.satm_coef,
        i.raot550nm, i.iband, i.normext, i.rotoa, i.eps,
    )
}

// Denominator well-conditioned: skip near-singular draws where any tiny
// rounding difference is meaninglessly amplified.
fn well_conditioned(c: f32, r: f32) -> bool {
    c.is_finite() && r.is_finite() && c.abs() < 10.0 && r.abs() < 10.0
}

proptest! {
    // eps < 0: pow() branch skipped -> pure f32 -> expect bit-exact.
    #[test]
    fn matches_c_no_pow(i in inputs(-2.0f32..0.0)) {
        let c = run_c(&i);
        let r = run_rust(&i);
        prop_assume!(well_conditioned(c, r));
        prop_assert_eq!(
            r.to_bits(), c.to_bits(),
            "mismatch (eps<0): rust={} c={} inputs={:?}", r, c, i
        );
    }

    // eps >= 0: includes double pow() -> allow a tight tolerance.
    #[test]
    fn matches_c_with_pow(i in inputs(0.0f32..2.0)) {
        let c = run_c(&i);
        let r = run_rust(&i);
        prop_assume!(well_conditioned(c, r));
        let tol = 1e-5 * c.abs().max(1e-3);
        prop_assert!(
            (r - c).abs() <= tol,
            "mismatch (eps>=0): rust={} c={} diff={} inputs={:?}",
            r, c, (r - c).abs(), i
        );
    }
}
