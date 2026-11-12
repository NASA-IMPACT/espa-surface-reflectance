# GDAS ancillary integration for LaSRC

This branch (`feature/gdas-integration`) adds NOAA GDAS (GRIB) water vapor and ozone as an
alternative to VIIRS/MODIS ancillary data in LaSRC. It is the code the HLS team validated
in 2026.

| | |
|---|---|
| Code base | Eric Vermote's changeset, PR #19, last commit `b271310` (26 Jun 2026), on top of release v3.5.1.0 (`6051e02`) |
| Version string | `3.5.2.0 (Collection 2)` |
| Branch changes on top of `b271310` | one housekeeping commit (file permissions restored, commented-out debug prints removed; no functional change) and this `docs/gdas_integration/` folder |
| Known issues | [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) |

## 1. What changed in the C code

All changes are in `lasrc/c_version/src/`.

| File | Change |
|---|---|
| `grib.c`, `grib.h` | New. GRIB1 reader (legacy `gdas1.PGrbF00.*` files). |
| `PGE11_grib2.c`, `PGE11_grib2.h`, `PGE11_grb2.h` | New. GRIB2 reader (current GDAS `pgrb2.1p00` files). About 5,000 lines of vendored decoder code. |
| `ad_hoc_XML_reader.c` | New. Reads the Sentinel-2 scene-centre time, so the code knows the acquisition time to interpolate GDAS to. |
| `lut_subr.c` | New `open_gb_oz_wv_datasets()`. It reads `PWAT` and `TOZNE` from the two GDAS files, interpolates them linearly to the acquisition time, and resamples the 1° grid to the 0.05° CMG grid that LaSRC uses (bilinear: `get_cart_linesamp()`, `get_cart_interp()`). It also adds `Grb2_file()`, which detects GRIB1 or GRIB2 from the file header. |
| `lasrc.c`, `get_args.c`, `lasrc.h` | New `--aux2` option. GRIB is detected from the file contents, not the file name. GRIB files are read from `$LASRC_AUX_DIR/LADS/` (no year sub-folder). |
| `compute_refl_subr.c`, `compute_landsat_refl.c`, `compute_sentinel_refl.c` | Pass the second aux file name and the acquisition date through to `read_auxiliary_files()`. |
| `common.h` | New `GRIB` entry in `aux_src_t` with its scale factors and defaults. The water-vapor scale factor for GRIB is the same as for VIIRS/MODIS, so GDAS inherits their factor-of-2 behaviour (see KNOWN_ISSUES #1). Version set to 3.5.2.0. |
| `compute_sentinel_refl.c` | Bug fix: the fill check now uses `qaband[curr_pix]`, not `qaband[i]`. |
| `input.c` | Sentinel-2 acquisition date/time is taken from the scene-centre time, not `T00:00:00`. |
| `output.h` | Reflectance output encoding changed from `0.0000275 / -0.2` to `0.0001 / 0.0`. **This changes every SR output, not just GDAS runs** (see KNOWN_ISSUES #5). |

The Python wrappers in `lasrc/c_version/scripts/` (`do_lasrc_sentinel.py`, `do_lasrc_landsat.py`) are
**not** updated in this changeset: they still call LaSRC with a single `--aux`. The wrapper the
validation runs used is kept for reference in [`reference/`](reference/) (section 4).

## 2. How to run LaSRC with GDAS

```
lasrc --xml=<scene>.xml --aux=<gdas file before> --aux2=<gdas file after> --verbose
```

- **Two GDAS analysis files:** the 6-hourly cycles before and after the acquisition time (cycles are at 00, 06, 12 and 18 UTC). LaSRC checks that:
  - the two files are exactly 6 h apart
  - they fall on the same or adjacent days
  - the acquisition time lies between them.

  Otherwise it stops with a "GRIB file staging problem" error.
- **Use `--aux2=FILE` with `=`:** the option is declared `optional_argument`, so `--aux2 FILE` is silently ignored (KNOWN_ISSUES #8).
- **File location:** files must be in `$LASRC_AUX_DIR/LADS/` directly (not `LADS/<year>/`). Any file name works, because the type is detected from the GRIB header. The validation runs named them `gdas1.PGrbF00.YYMMDD.HHz`.
- **Variables used:**
  - `PWAT`, precipitable water in kg/m²
  - `TOZNE`, total ozone in Dobson units
  - both at level "entire atmosphere (considered as a single layer)".
- **Static aux files** are still required as before: `MSILUT/`, `CMGDEM.hdf`, `ratiomapndwiexp.hdf`, `monthly_avgs`. A missing one causes "Could not find cmgdemnm/rationm".

### Where the GDAS files come from

- **Source:** NOAA GFS/GDAS on AWS Open Data, `s3://noaa-gfs-bdp-pds/`. No credentials are needed.
- **Path from 22 Mar 2021 (GFS v16) onward:** `gdas.YYYYMMDD/HH/atmos/gdas.tHHz.pgrb2.1p00.f000`.
- **Path before 22 Mar 2021:** `gdas.YYYYMMDD/HH/gdas.tHHz.pgrb2.1p00.f000` (no `atmos/` folder).
- **Product:** 1° global grid, analysis (`f000`). A cycle is usually posted about 7 h after its nominal time. The second file for a scene therefore arrives roughly 7–13 h after the overpass, and that sets the earliest time a scene can be processed.
- **Upcoming change:** NOAA has proposed GFS/GDAS v17. Folder and file names are expected to change, and continuity of the 1° product should be confirmed before the ingest is frozen.

## 3. Validation summary (HLS team, 2026)

All runs used this code base (3.5.2.0) inside the HLS Sentinel/Landsat containers on UAH Meteor.
VIIRS and GDAS runs used identical inputs; only the ancillary source differed.

| Test | Result |
|---|---|
| S30, 522 granules (Jan and Jul 2025), final HLS products | 4,697 of 4,698 band checks within the Ju et al. (2025) thresholds; mean GDAS − VIIRS differences 0.3–3 HLS DN |
| L30, 614 granules (Jan and Jul 2025) | Included in the same GDAS-vs-VIIRS comparison (LaSRC output stage, Landsat test container) |
| Runtime | 2,027 s (GDAS) vs 2,029 s (VIIRS) per granule; no extra cost |
| Time series, 4 RadCalNet sites (Railroad Valley, La Crau, Baotou sand, Gobabeb), 7 bands, Jan 2024–Dec 2025, 185–258 paired scenes per site | No significant trend in the GDAS − VIIRS difference at any site or band (Mann-Kendall p ≥ 0.25). Largest trend +0.39 DN/yr (Gobabeb B12) = 0.6% of the Ju threshold per year. Largest mean offset 7.4 DN (Gobabeb B12) = 12% of the threshold. |
| Reproduction of operational HLS | At Railroad Valley, our VIIRS run matched operational HLS S30 within 0.5 DN on 104 of 129 clear dates (2024–2025) |

The validation compared GDAS with VIIRS **with the known issues below present in both runs** where they
apply. Fixing issues #2, #3, #4 or #6 changes GDAS outputs only, so after any fix the GDAS-vs-VIIRS
comparison should be repeated on a small reference set.

Unrelated to GDAS: an aerosol-repair artifact over snow/ice (smoke-like plumes in B01–B04) appears with
both ancillary sources. A candidate fix (J. Ju, `docs/jju-review` branch) is under test.

## 4. Reference scripts (documentation only, not wired in)

[`reference/`](reference/) holds the scripts used for the Meteor validation runs. They are kept to
show how the pieces fit together. They are not part of the build and contain site-specific paths.

| File | What it shows |
|---|---|
| `reference/do_lasrc_sentinel_gdas.py` | Modified `do_lasrc_sentinel.py`. It picks the GDAS file pair bracketing the scene time and calls `lasrc --aux=... --aux2=...`. It adds `LASRC_AUX_MODE=modis\|viirs\|grib` (forces a source) and `GRIB_STARTING_DATE` (date from which GDAS is used automatically), and streams LaSRC output to the log. |
| `reference/prepare_grib_aux.sh` | Stages the two GDAS `pgrb2.1p00.f000` files for a granule from a local mirror into `$LASRC_AUX_DIR/LADS/` as `gdas1.PGrbF00.YYMMDD.HHz`. It only knows the post-March-2021 `atmos/` folder layout. |

## 5. Notes for the Rust port

- **Water vapor factor of 2 (KNOWN_ISSUES #1):** reproduce it on purpose, behind a named constant, so Rust output matches C and operations. Whether to fix it is a separate science decision.
- **GRIB decoding:** prefer an established decoder (ecCodes bindings or a maintained GRIB2 crate) over porting the 5,000 lines of `PGE11_grib2.c`. Only `PWAT` and `TOZNE` from the 1° analysis are needed.
- **Make the GDAS steps small, testable functions:** time interpolation, grid reversal and resampling to the CMG grid. Unit-test them on the cases where the C code goes wrong: the antimeridian, the poles, half-cell alignment, and a 1° longitude check against a known field (KNOWN_ISSUES #2–#4).
- **Regression test:** keep a small set of reference granules (S30 and L30, VIIRS and GDAS) with C outputs, and compare Rust outputs band by band against the Ju et al. (2025) thresholds, as in the validation.
- **Output encoding:** decide the output encoding (`0.0001/0` vs `0.0000275/-0.2`) explicitly before porting (KNOWN_ISSUES #5).
