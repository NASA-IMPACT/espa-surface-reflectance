//! Compile the original C LaSRC leaf functions for differential testing.
//!
//! The C sources live in this repo at ../lasrc/c_version/src. They are split
//! into staged groups by dependency weight:
//!
//!   Stage 0 (always): zero-dependency leaves that include only their own
//!     header -- compile with any C toolchain, no conda env needed.
//!   Stage 1+ (feature "lut"): leaves that pull in common.h -> hdf.h and
//!     lut_subr.h -> espa_*.h. These headers come from the espa-product-formatter
//!     conda package; point $PREFIX (or $CREF_PREFIX) at the activated pixi env.
//!
//! Floating-point note: the C is compiled with -ffp-contract=off so the
//! reference does not fuse multiply-adds. Run inside the activated conda env so
//! `cc` inherits the same $CC/$CFLAGS used to build the shipped binary, keeping
//! the reference codegen faithful to production.

use std::env;
use std::path::PathBuf;

fn main() {
    let c_src = PathBuf::from("../lasrc/c_version/src");
    println!("cargo:rerun-if-changed=build.rs");

    let mut build = cc::Build::new();
    build.include(&c_src);
    build.flag_if_supported("-ffp-contract=off");

    // ---- Stage 0: zero-dependency leaves (compile anywhere) ----
    let stage0 = ["quick_select.c"];
    for f in stage0 {
        build.file(c_src.join(f));
        println!("cargo:rerun-if-changed={}", c_src.join(f).display());
    }

    // ---- Stage 1+: ESPA/HDF-dependent leaves (need a conda/pixi env) ----
    if env::var_os("CARGO_FEATURE_LUT").is_some() {
        let prefix = env::var("CREF_PREFIX")
            .or_else(|_| env::var("PREFIX"))
            .expect(
                "feature `lut` requires $PREFIX or $CREF_PREFIX to point at a conda \
                 env containing espa-product-formatter + hdf headers/libs",
            );
        build.include(format!("{prefix}/include"));
        // libxml2 headers live in their own subdir (espa build: XML2INC=.../include/libxml2).
        build.include(format!("{prefix}/include/libxml2"));

        let stage1 = ["poly_coeff.c", "lut_subr.c", "subaeroret.c"];
        for f in stage1 {
            build.file(c_src.join(f));
            println!("cargo:rerun-if-changed={}", c_src.join(f).display());
        }

        // lut_subr.c's (test-unused) table readers reference HDF4/HDF5; link
        // them so those symbols resolve. Adjust the list if the env's library
        // names differ.
        println!("cargo:rustc-link-search=native={prefix}/lib");
        // readluts uses ESPA error_handler (lib_espa_common.a, pulls libxml2).
        println!("cargo:rustc-link-arg=-l_espa_common");
        // HDF4/5 for the table reads. Wrap in --no-as-needed: libmfhdf's Hoffset
        // lives in libdf, which our code never references directly, so the
        // default --as-needed would drop libdf and break at runtime. Raw link
        // args (the `-as-needed` link-lib modifier is nightly-only).
        println!("cargo:rustc-link-arg=-Wl,--no-as-needed");
        for lib in ["mfhdf", "df", "hdf5", "hdf5_hl"] {
            println!("cargo:rustc-link-arg=-l{lib}");
        }
        println!("cargo:rustc-link-arg=-Wl,--as-needed");
        for lib in ["xml2", "jpeg", "z"] {
            println!("cargo:rustc-link-arg=-l{lib}");
        }
    }

    build.compile("lasrc_cref_c");
}
