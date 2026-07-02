//! Differential test: C `subaeroret_new` vs `lasrc_core::aerosol::subaeroret_new`.
//!
//! THE documented precision-gap source (docs/superpowers .../lasrc-next-steps).
//! The convergence loop picks a discrete AOT index from a residual comparison,
//! so any f32-vs-f64 rounding difference can flip the selected index and jump
//! the retrieved `raot` to a different discrete value. This test is expected to
//! LOCALIZE that gap if the Rust port's float emulation is not yet exact -- a
//! red result here is signal, not a harness bug.
//!
//! Mapping pinned (see src/lib.rs): sat=Landsat8, N=7 bands (DNL_BAND1..7),
//! C `roatm_iaMax` int indices vs Rust AOT550NM-resolved values, C-internal
//! `tth` replicated as LANDSAT_TTH, lambda=LAMBDA_LANDSAT. First cut fixes
//! eps<0 to skip the pow() path and isolate the convergence arithmetic.
#![cfg(feature = "lut")]

use lasrc_core::aerosol::subaeroret_new as rust_subaeroret_new;
use lasrc_core::constants::LAMBDA_LANDSAT;
use lasrc_cref::{c_subaeroret_new, AOT550NM, LANDSAT_TTH, LANDSAT_TTH_WATER, SAT_LANDSAT_8};
use proptest::prelude::*;

const N: usize = 7; // DNL_BAND1..=DNL_BAND7

fn coef_row(lo: f32, hi: f32) -> impl Strategy<Value = [f32; 4]> {
    [(lo..hi), (lo..hi), (lo..hi), (0.0f32..1.0)]
}

fn coef_rows(lo: f32, hi: f32) -> impl Strategy<Value = Vec<[f32; 4]>> {
    proptest::collection::vec(coef_row(lo, hi), N)
}

/// One property case: run C + Rust with identical inputs and assert the
/// decisive outputs (converged AOT index, retrieved AOT) agree bit-for-bit.
#[allow(clippy::too_many_arguments)]
fn check_case(
    water: bool,
    iband1: usize,
    erelc: &[f32],
    troatm: &[f32],
    tgo: &[f32],
    ia_max: &[i32],
    roatm_coef: &[[f32; 4]],
    ttatmg_coef: &[[f32; 4]],
    satm_coef: &[[f32; 4]],
    normext: &[f32],
    iaots_in: i32,
    eps: f32,
) -> Result<(), TestCaseError> {
    // ---- C side (int indices, C-internal tth) ----
    let c = c_subaeroret_new(
        SAT_LANDSAT_8, water, iband1 as i32,
        erelc, troatm, tgo, ia_max,
        roatm_coef, ttatmg_coef, satm_coef, normext,
        iaots_in, eps,
    );

    // ---- Rust side (values widened to f64; ia_max resolved via AOT550NM) ----
    let f64v = |v: &[f32]| v.iter().map(|&x| x as f64).collect::<Vec<_>>();
    let rows = |r: &[[f32; 4]]| r.iter().map(|c| c.map(|x| x as f64)).collect::<Vec<_>>();
    let roatm_ia_max: Vec<f64> = ia_max.iter().map(|&i| AOT550NM[i as usize] as f64).collect();
    // C picks tth by water flag internally; feed the matching table.
    let tth_src = if water { &LANDSAT_TTH_WATER } else { &LANDSAT_TTH };
    let tth: Vec<f64> = tth_src[..N].iter().map(|&x| x as f64).collect();

    let r = rust_subaeroret_new(
        water, iband1,
        &f64v(erelc), &f64v(troatm), &f64v(tgo),
        &roatm_ia_max, &rows(roatm_coef), &rows(ttatmg_coef), &rows(satm_coef),
        &f64v(normext), &LAMBDA_LANDSAT, eps as f64, iaots_in as usize, &tth,
    );

    prop_assume!(c.raot.is_finite() && r.raot.is_finite());
    prop_assert_eq!(r.iaots as i32, c.iaots, "iaots differ: rust={} c={}", r.iaots, c.iaots);
    prop_assert_eq!(
        (r.raot as f32).to_bits(), c.raot.to_bits(),
        "raot differ: rust={} c={} eps={}", r.raot as f32, c.raot, eps
    );
    Ok(())
}

proptest! {
    // eps < 0 -> pow() branch skipped; pure f32 convergence arithmetic.
    #[test]
    fn matches_c_no_pow(
        water in any::<bool>(),
        iband1 in 0usize..N,
        erelc in proptest::collection::vec(0.05f32..1.0, N),
        troatm in proptest::collection::vec(0.0f32..0.5, N),
        tgo in proptest::collection::vec(0.5f32..1.0, N),
        ia_max in proptest::collection::vec(0i32..22, N),
        roatm_coef in coef_rows(-0.05, 0.05),
        ttatmg_coef in coef_rows(0.0, 0.3),
        satm_coef in coef_rows(-0.05, 0.05),
        normext in proptest::collection::vec(0.5f32..2.0, N),
        iaots_in in 0i32..4,
        eps in -2.0f32..0.0,
    ) {
        check_case(water, iband1, &erelc, &troatm, &tgo, &ia_max,
            &roatm_coef, &ttatmg_coef, &satm_coef, &normext, iaots_in, eps)?;
    }
}

proptest! {
    // eps >= 0 -> includes the pow() modification; where the docs expect the
    // f32/f64 gap. More cases to hunt rare residual-driven index flips.
    #![proptest_config(ProptestConfig { cases: 4096, ..ProptestConfig::default() })]
    #[test]
    fn matches_c_with_pow(
        water in any::<bool>(),
        iband1 in 0usize..N,
        erelc in proptest::collection::vec(0.05f32..1.0, N),
        troatm in proptest::collection::vec(0.0f32..0.5, N),
        tgo in proptest::collection::vec(0.5f32..1.0, N),
        ia_max in proptest::collection::vec(0i32..22, N),
        roatm_coef in coef_rows(-0.05, 0.05),
        ttatmg_coef in coef_rows(0.0, 0.3),
        satm_coef in coef_rows(-0.05, 0.05),
        normext in proptest::collection::vec(0.5f32..2.0, N),
        iaots_in in 0i32..4,
        eps in 0.0f32..2.0,
    ) {
        check_case(water, iband1, &erelc, &troatm, &tgo, &ia_max,
            &roatm_coef, &ttatmg_coef, &satm_coef, &normext, iaots_in, eps)?;
    }
}
