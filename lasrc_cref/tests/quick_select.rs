//! Differential test: C `quick_select` vs `lasrc_core::utils::quick_select`.
//!
//! Both are line-for-line median-of-three quickselect ports. They agree
//! bit-for-bit for odd-length arrays. They intentionally DIVERGE for even
//! lengths: C selects index (n-1)/2 (lower-middle) while the Rust port selects
//! n/2 (upper-middle). Real LaSRC aerosol windows are always odd (3x3=9,
//! 15x15=225), so the property test fuzzes odd lengths; the even-length
//! divergence is pinned by an explicit witness below so it is documented rather
//! than hidden. This is the kind of off-by-one the harness exists to surface.

use lasrc_core::utils::quick_select as rust_quick_select;
use lasrc_cref::c_quick_select;
use proptest::prelude::*;

/// Finite f32 values only. NaN has no total order, which makes the result of
/// any selection algorithm implementation-defined and the comparison
/// meaningless; infinities are excluded for the same robustness reason.
fn finite_f32() -> impl Strategy<Value = f32> {
    use proptest::num::f32::{NORMAL, SUBNORMAL, ZERO};
    NORMAL | SUBNORMAL | ZERO
}

proptest! {
    #[test]
    fn matches_c_for_odd_lengths(
        // Odd length n = 2*half + 1 in [1, 401], covering the real window sizes.
        half in 0usize..=200,
        values in proptest::collection::vec(finite_f32(), 401),
    ) {
        let n = 2 * half + 1;
        let mut a_rust: Vec<f32> = values[..n].to_vec();
        let mut a_c = a_rust.clone();

        let r = rust_quick_select(&mut a_rust);
        let c = c_quick_select(&mut a_c);

        // Pure comparison/swap algorithm: expect bit-exact agreement.
        prop_assert_eq!(
            r.to_bits(), c.to_bits(),
            "median mismatch at n={}: rust={} c={}", n, r, c
        );
    }
}

/// Explicit witness for the known, intentional even-length divergence. For a
/// sorted distinct array of length 4, C returns the 2nd-smallest (index
/// (n-1)/2 = 1) and the Rust port returns the 3rd-smallest (index n/2 = 2).
/// Recorded here so the behavior is captured; LaSRC never calls this with even n.
#[test]
fn even_length_divergence_is_documented() {
    let mut a_c = [10.0f32, 20.0, 30.0, 40.0];
    let mut a_rust = a_c;

    let c = c_quick_select(&mut a_c);
    let r = rust_quick_select(&mut a_rust);

    assert_eq!(c, 20.0, "C selects lower-middle, index (n-1)/2");
    assert_eq!(r, 30.0, "Rust selects upper-middle, index n/2");
}
