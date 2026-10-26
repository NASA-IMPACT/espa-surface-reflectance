//! LUT IO cross-check, part 1: dump the arrays C `readluts` produces so a Python
//! script can diff them against `lasrc.aux.load_lut`. This is the Python-vs-C IO
//! bridge -- readluts + load_lut consume the SAME four files, so byte-diffing
//! their output arrays catches reshape/stride/band-remap bugs in the Python
//! loader (the class of divergence the leaf-kernel tests cannot see).
//!
//! Skips unless the four aux paths + a dump dir are given via env:
//!   CREF_ANGLE, CREF_INTREF, CREF_TRANSM, CREF_SPHERA, CREF_DUMP_DIR
//!   (optional CREF_XTSSTEP, CREF_XTSMIN; default 4.0 / 0.0)
//! Needs --features lut and the real Landsat aux files mounted.
#![cfg(feature = "lut")]

use lasrc_cref::c_readluts;
use std::io::Write;
use std::path::Path;

fn write_f32(dir: &Path, name: &str, v: &[f32]) {
    let mut f = std::fs::File::create(dir.join(name)).unwrap();
    // Raw native-endian f32 (both build + compare hosts are little-endian).
    let bytes = unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) };
    f.write_all(bytes).unwrap();
}

fn write_i32(dir: &Path, name: &str, v: &[i32]) {
    let mut f = std::fs::File::create(dir.join(name)).unwrap();
    let bytes = unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) };
    f.write_all(bytes).unwrap();
}

#[test]
fn dump_c_readluts() {
    let (angle, intref, transm, sphera, dump) = match (
        std::env::var("CREF_ANGLE"),
        std::env::var("CREF_INTREF"),
        std::env::var("CREF_TRANSM"),
        std::env::var("CREF_SPHERA"),
        std::env::var("CREF_DUMP_DIR"),
    ) {
        (Ok(a), Ok(i), Ok(t), Ok(s), Ok(d)) => (a, i, t, s, d),
        _ => {
            eprintln!("SKIP dump_c_readluts: set CREF_ANGLE/INTREF/TRANSM/SPHERA/DUMP_DIR");
            return;
        }
    };
    let xtsstep = std::env::var("CREF_XTSSTEP").ok().and_then(|s| s.parse().ok()).unwrap_or(4.0f32);
    let xtsmin = std::env::var("CREF_XTSMIN").ok().and_then(|s| s.parse().ok()).unwrap_or(0.0f32);

    let dir = Path::new(&dump);
    std::fs::create_dir_all(dir).unwrap();

    let l = c_readluts(&angle, &intref, &transm, &sphera, xtsstep, xtsmin);

    write_f32(dir, "rolutt.bin", &l.rolutt);
    write_f32(dir, "transt.bin", &l.transt);
    write_f32(dir, "sphalbt.bin", &l.sphalbt);
    write_f32(dir, "normext.bin", &l.normext);
    write_f32(dir, "tsmax.bin", &l.tsmax);
    write_f32(dir, "tsmin.bin", &l.tsmin);
    write_f32(dir, "ttv.bin", &l.ttv);
    write_f32(dir, "tts.bin", &l.tts);
    write_f32(dir, "nbfic.bin", &l.nbfic);
    write_f32(dir, "nbfi.bin", &l.nbfi);
    write_i32(dir, "indts.bin", &l.indts);

    eprintln!("wrote C readluts dumps to {}", dir.display());
}
