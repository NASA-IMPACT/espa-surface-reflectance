//! Differential test: C `aerosol_interp_landsat` (via cref shim) vs
//! `lasrc_core::aerosol::aerosol_interp`.
//!
//! Bilinear fill of taero/ipflag from window centers -- the routine the docs
//! flag for fill-boundary artifacts. Both run on identical grids; compare the
//! full taero (bit-exact, pure f32 bilinear) and ipflag arrays. Fuzzes the QA
//! fill mask to exercise the boundary/fill handling.
#![cfg(feature = "lut")]

use lasrc_core::aerosol::aerosol_interp as rust_aerosol_interp;
use lasrc_cref::c_aerosol_interp_landsat;
use proptest::prelude::*;

const NL: usize = 12; // multiple of the window so centers tile cleanly
const NS: usize = 12;
const AW: usize = 3; // LAERO_WINDOW
const HALF: usize = 1; // LHALF_AERO_WINDOW
const N: usize = NL * NS;

proptest! {
    #[test]
    fn matches_c(
        taero0 in proptest::collection::vec(0.0f32..2.0, N),
        ipf0 in proptest::collection::vec(0u8..4, N),
        // 0 = valid, 1 = fill (level-1 QA bit 0)
        qa in proptest::collection::vec(0u16..2, N),
    ) {
        let mut t_c = taero0.clone();
        let mut i_c = ipf0.clone();
        c_aerosol_interp_landsat(&qa, &mut i_c, &mut t_c, NL, NS, AW, HALF);

        let mut t_r = taero0.clone();
        let mut i_r = ipf0.clone();
        rust_aerosol_interp(&mut t_r, &mut i_r, &qa, NL, NS, AW, HALF);

        for p in 0..N {
            prop_assert_eq!(
                t_r[p].to_bits(), t_c[p].to_bits(),
                "taero[{}] differ: rust={} c={}", p, t_r[p], t_c[p]
            );
            prop_assert_eq!(
                i_r[p], i_c[p],
                "ipflag[{}] differ: rust={} c={}", p, i_r[p], i_c[p]
            );
        }
    }
}
