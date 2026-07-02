//! Differential test harness linking the original C LaSRC leaf functions
//! against the `lasrc_core` Rust port.
//!
//! Each `extern "C"` block mirrors a function in `lasrc/c_version/src`, compiled
//! by `build.rs`. The safe wrappers here give the property tests in `tests/` a
//! Rust-shaped call so they can feed identical inputs to both implementations
//! and assert agreement.

use std::os::raw::c_int;
#[cfg(feature = "lut")]
use std::ffi::CString;
#[cfg(feature = "lut")]
use std::os::raw::c_char;

// ---- Stage 0: zero-dependency leaves ----

extern "C" {
    /// C `quick_select(float *arr, int n)` from quick_select.c. Returns the
    /// median (the (n-1)/2-th order statistic) and reorders `arr` in place.
    fn quick_select(arr: *mut f32, n: c_int) -> f32;
}

/// Median via the original C `quick_select`. Reorders `arr` in place, as the C
/// does. Panics on an empty slice (the C contract requires n >= 1).
pub fn c_quick_select(arr: &mut [f32]) -> f32 {
    assert!(!arr.is_empty(), "quick_select requires a non-empty array");
    // SAFETY: `arr` is a valid, non-empty, contiguous f32 buffer of length
    // `arr.len()`, which matches the (ptr, n) contract of the C function.
    unsafe { quick_select(arr.as_mut_ptr(), arr.len() as c_int) }
}

// ---- Stage 1: ESPA/HDF-dependent leaves (feature "lut") ----

/// Sat_t values from common.h: Landsat 8 = 0, Landsat 9 = 1, Sentinel-2 = 2.
/// C `atmcorlamb2_new` uses this to pick the internal `lambda` table and the
/// max band index; the Rust port takes `lambda` explicitly instead.
#[cfg(feature = "lut")]
pub const SAT_LANDSAT_8: c_int = 0;

#[cfg(feature = "lut")]
extern "C" {
    /// C `atmcorlamb2_new` from lut_subr.c (polynomial fast path). `roatm_coef`,
    /// `ttatmg_coef`, `satm_coef` are `float[NCOEF]` (NCOEF = 4). `roslamb` is the
    /// single output.
    fn atmcorlamb2_new(
        sat: c_int,
        tgo: f32,
        roatm_upper: f32,
        roatm_coef: *const f32,
        ttatmg_coef: *const f32,
        satm_coef: *const f32,
        raot550nm: f32,
        iband: c_int,
        normext_ib_0_3: f32,
        rotoa: f32,
        roslamb: *mut f32,
        eps: f32,
    );
}

/// Surface reflectance via the original C `atmcorlamb2_new` (poly fast path).
#[cfg(feature = "lut")]
#[allow(clippy::too_many_arguments)]
pub fn c_atmcorlamb2_new(
    sat: c_int,
    tgo: f32,
    roatm_upper: f32,
    roatm_coef: &[f32; 4],
    ttatmg_coef: &[f32; 4],
    satm_coef: &[f32; 4],
    raot550nm: f32,
    iband: c_int,
    normext_ib_0_3: f32,
    rotoa: f32,
    eps: f32,
) -> f32 {
    let mut roslamb = 0.0f32;
    // SAFETY: the three coef slices are exactly NCOEF=4 long as C expects, and
    // `roslamb` is a valid out-pointer.
    unsafe {
        atmcorlamb2_new(
            sat, tgo, roatm_upper,
            roatm_coef.as_ptr(), ttatmg_coef.as_ptr(), satm_coef.as_ptr(),
            raot550nm, iband, normext_ib_0_3, rotoa, &mut roslamb, eps,
        );
    }
    roslamb
}

/// AOT value table `aot550nm[NAOT_VALS]` hardcoded in C `subaeroret_new`. The C
/// `roatm_iaMax[ib]` are indices into this; the Rust port takes the resolved
/// values, so the test maps `AOT550NM[idx]` when feeding the Rust side.
#[cfg(feature = "lut")]
pub const AOT550NM: [f32; 22] = [
    0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.3, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0,
];

/// Landsat `tth` thresholds hardcoded in C `subaeroret_new` (land / water).
/// C selects internally from sat+water; the Rust port takes `tth` explicitly.
#[cfg(feature = "lut")]
pub const LANDSAT_TTH: [f32; 8] =
    [1.0e-3, 1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-4, 0.0];
#[cfg(feature = "lut")]
pub const LANDSAT_TTH_WATER: [f32; 8] =
    [1.0e-3, 1.0e-3, 0.0, 1.0e-3, 1.0e-3, 0.0, 1.0e-4, 0.0];

#[cfg(feature = "lut")]
extern "C" {
    /// C `subaeroret_new` from subaeroret.c. 2D coef args are `float[N][NCOEF]`,
    /// passed here as row-major flat pointers. `raot`/`residual` are outputs;
    /// `iaots` is in/out (starting AOT index in, converged index out).
    fn subaeroret_new(
        sat: c_int,
        water: bool,
        iband1: c_int,
        erelc: *const f32,
        troatm: *const f32,
        tgo_arr: *const f32,
        roatm_ia_max: *const c_int,
        roatm_coef: *const f32,
        ttatmg_coef: *const f32,
        satm_coef: *const f32,
        normext_p0a3: *const f32,
        raot: *mut f32,
        residual: *mut f32,
        iaots: *mut c_int,
        eps: f32,
    );
}

#[cfg(feature = "lut")]
extern "C" {
    /// C `get_3rd_order_poly_coeff` from poly_coeff.c. Fits a cubic to
    /// (aot[i], atm[i]) via LU-solved normal equations; writes coeff[NCOEF=4].
    fn get_3rd_order_poly_coeff(
        aot: *const f32,
        atm: *const f32,
        n_atm: c_int,
        coeff: *mut f32,
    );
}

/// Cubic fit coefficients via the original C `get_3rd_order_poly_coeff`.
#[cfg(feature = "lut")]
pub fn c_get_3rd_order_poly_coeff(aot: &[f32], atm: &[f32]) -> [f32; 4] {
    assert_eq!(aot.len(), atm.len());
    let mut coeff = [0.0f32; 4];
    // SAFETY: aot/atm are n_atm long; coeff is NCOEF=4 as C writes.
    unsafe {
        get_3rd_order_poly_coeff(aot.as_ptr(), atm.as_ptr(), aot.len() as c_int, coeff.as_mut_ptr());
    }
    coeff
}

#[cfg(feature = "lut")]
extern "C" {
    /// C `local_chand` from lut_subr.c: molecular (Rayleigh) reflectance.
    fn local_chand(xphi: f32, xmuv: f32, xmus: f32, xtau: f32, xrray: *mut f32);
}

/// Rayleigh reflectance via the original C `local_chand`.
#[cfg(feature = "lut")]
pub fn c_local_chand(xphi: f32, xmuv: f32, xmus: f32, xtau: f32) -> f32 {
    let mut xrray = 0.0f32;
    // SAFETY: xrray is a valid out-pointer; all other args are by-value.
    unsafe { local_chand(xphi, xmuv, xmus, xtau, &mut xrray); }
    xrray
}

/// Output of the C `subaeroret_new` aerosol retrieval.
#[cfg(feature = "lut")]
pub struct CAeroResult {
    pub raot: f32,
    pub residual: f32,
    pub iaots: i32,
}

/// Run the original C `subaeroret_new`. `roatm_ia_max` are the C integer
/// indices into `AOT550NM`. Coef slices are `[[f32; 4]; N]` row-major.
#[cfg(feature = "lut")]
#[allow(clippy::too_many_arguments)]
pub fn c_subaeroret_new(
    sat: c_int,
    water: bool,
    iband1: c_int,
    erelc: &[f32],
    troatm: &[f32],
    tgo_arr: &[f32],
    roatm_ia_max: &[c_int],
    roatm_coef: &[[f32; 4]],
    ttatmg_coef: &[[f32; 4]],
    satm_coef: &[[f32; 4]],
    normext_p0a3: &[f32],
    iaots_in: i32,
    eps: f32,
) -> CAeroResult {
    let mut raot = 0.0f32;
    let mut residual = 0.0f32;
    let mut iaots = iaots_in;
    // SAFETY: all slices are >= the band range C indexes (start..=DNL_BAND7);
    // coef slices are row-major [N][NCOEF]; out-pointers are valid.
    unsafe {
        subaeroret_new(
            sat, water, iband1,
            erelc.as_ptr(), troatm.as_ptr(), tgo_arr.as_ptr(),
            roatm_ia_max.as_ptr(),
            roatm_coef.as_ptr() as *const f32,
            ttatmg_coef.as_ptr() as *const f32,
            satm_coef.as_ptr() as *const f32,
            normext_p0a3.as_ptr(),
            &mut raot, &mut residual, &mut iaots, eps,
        );
    }
    CAeroResult { raot, residual, iaots }
}

// ---- LUT loading cross-check: C readluts vs Python lasrc.aux.load_lut ----
// Landsat-sized (nsr_bands = 8). C readluts writes only bands 0..7 for Landsat,
// so these buffers match the Python load_lut output layout exactly.

/// LUT grid sizes for a Landsat load (see common.h).
#[cfg(feature = "lut")]
pub mod lutdim {
    pub const NSR: usize = 8; // NSRL_BANDS
    pub const NPRES: usize = 7;
    pub const NAOT: usize = 22;
    pub const NSOLAR: usize = 8000;
    pub const NSUNANGLE: usize = 22;
    pub const NVIEW_ZEN: usize = 20;
    pub const NSOLAR_ZEN: usize = 22;
    pub const ANGLE: usize = NVIEW_ZEN * NSOLAR_ZEN; // tsmax/tsmin/ttv/nbfic/nbfi
    pub const ROLUTT: usize = NSR * NPRES * NAOT * NSOLAR;
    pub const TRANST: usize = NSR * NPRES * NAOT * NSUNANGLE;
    pub const SPHNORM: usize = NSR * NPRES * NAOT; // sphalbt / normext
}

#[cfg(feature = "lut")]
extern "C" {
    fn readluts(
        sat: c_int,
        tsmax: *mut f32,
        tsmin: *mut f32,
        ttv: *mut f32,
        tts: *mut f32,
        nbfic: *mut f32,
        nbfi: *mut f32,
        indts: *mut i32,
        rolutt: *mut f32,
        transt: *mut f32,
        sphalbt: *mut f32,
        normext: *mut f32,
        xtsstep: f32,
        xtsmin: f32,
        anglehdf: *const c_char,
        intrefnm: *const c_char,
        transmnm: *const c_char,
        spheranm: *const c_char,
    ) -> c_int;
}

/// All arrays produced by C `readluts` (Landsat sizes).
#[cfg(feature = "lut")]
pub struct CLuts {
    pub tsmax: Vec<f32>,
    pub tsmin: Vec<f32>,
    pub ttv: Vec<f32>,
    pub tts: Vec<f32>,
    pub nbfic: Vec<f32>,
    pub nbfi: Vec<f32>,
    pub indts: Vec<i32>,
    pub rolutt: Vec<f32>,
    pub transt: Vec<f32>,
    pub sphalbt: Vec<f32>,
    pub normext: Vec<f32>,
}

/// Load the LUTs via the original C `readluts` (Landsat). Paths are the same
/// four files Python `load_lut` consumes (angle, intref/RES, transm, sphera).
#[cfg(feature = "lut")]
pub fn c_readluts(
    angle: &str,
    intref: &str,
    transm: &str,
    sphera: &str,
    xtsstep: f32,
    xtsmin: f32,
) -> CLuts {
    use lutdim::*;
    let mut o = CLuts {
        tsmax: vec![0.0; ANGLE],
        tsmin: vec![0.0; ANGLE],
        ttv: vec![0.0; ANGLE],
        tts: vec![0.0; NSOLAR_ZEN],
        nbfic: vec![0.0; ANGLE],
        nbfi: vec![0.0; ANGLE],
        indts: vec![0; NSUNANGLE],
        rolutt: vec![0.0; ROLUTT],
        transt: vec![0.0; TRANST],
        sphalbt: vec![0.0; SPHNORM],
        normext: vec![0.0; SPHNORM],
    };
    let ca = CString::new(angle).unwrap();
    let ci = CString::new(intref).unwrap();
    let ct = CString::new(transm).unwrap();
    let cs = CString::new(sphera).unwrap();
    // SAFETY: buffers are sized per common.h for a Landsat (nsr_bands=8) load;
    // path CStrings are nul-terminated as C's char[] expects.
    let status = unsafe {
        readluts(
            SAT_LANDSAT_8,
            o.tsmax.as_mut_ptr(), o.tsmin.as_mut_ptr(), o.ttv.as_mut_ptr(),
            o.tts.as_mut_ptr(), o.nbfic.as_mut_ptr(), o.nbfi.as_mut_ptr(),
            o.indts.as_mut_ptr(), o.rolutt.as_mut_ptr(), o.transt.as_mut_ptr(),
            o.sphalbt.as_mut_ptr(), o.normext.as_mut_ptr(),
            xtsstep, xtsmin,
            ca.as_ptr(), ci.as_ptr(), ct.as_ptr(), cs.as_ptr(),
        )
    };
    assert_eq!(status, 0, "C readluts failed (status {status})");
    o
}
