# Known issues in the GDAS changeset

Found while reviewing `b271310`. Line numbers refer to this branch (`feature/gdas-integration`).
Each entry can be copied into a GitHub issue.

**Status: drafts.** Nothing here has been fixed on this branch. Items marked "changes output" will
shift GDAS results slightly once fixed, so each fix should come with a before/after comparison on a
small reference set (see README section 3).

| # | Issue | Severity | Fixing changes output? |
|---|---|---|---|
| 1 | Water vapor is half its true value (same as the existing VIIRS/MODIS path) | Known behaviour | Yes, for all sources; science decision |
| 2 | GDAS grid is shifted 1° in longitude | High | Yes (GDAS) |
| 3 | Wrong neighbour at the 179°–180° longitude wrap | Medium | Yes (GDAS, ~1° strip) |
| 4 | Half-cell (0.025°) offset when resampling to the CMG grid | Low | Yes (GDAS, negligible) |
| 5 | SR output encoding changed to 0.0001/0: negative reflectance becomes fill | High | Yes (all sources) |
| 6 | Uninitialised `swapped` in the time interpolation | High (undefined behaviour) | Possibly (GDAS) |
| 7 | GRIB run without `--aux2` uses an unset file name | Medium | No (crash/robustness) |
| 8 | `--aux2` only works as `--aux2=FILE` | Low | No |
| 9 | Missing checks when a GRIB variable is not found | Medium | No (robustness) |
| 10 | Minor: truncation, format string, memory | Low | Truncation only |

---

## 1. Water vapor is half its true value (same as the existing VIIRS/MODIS path)

**Where:** `common.h:301` (`WV_SCALE_FACTOR {200, 200, 200}`), `lut_subr.c:3122`
(`wv[l] = (uint16)(value*10.0)`), `compute_refl_subr.c:248` (`uwv = wv / WV_SCALE_FACTOR`).

**What happens:**
- GDAS `PWAT` is in kg/m², and 1 g/cm² = 10 kg/m², so the correct value is kg/m² ÷ 10.
- The code stores kg/m² × 10 and later divides by 200, which gives kg/m² ÷ 20: half the true value.
- VIIRS has the same factor of 2: the file stores g/cm² × 100 and the code divides by 200.
- For VIIRS and MODIS the factor of 2 is documented as deliberate. Commit `ebea983` (G. Schmidt, 2023) says it "reintroduced the WV scale factor bug for both MODIS and VIIRS", and the comment on `WV_SCALE_FACTOR` in `common.h` reads "USGS OPS with MODIS bug, and added bug for VIIRS".
- For GDAS there is **no record of intent**. The GDAS changeset (`a64c9cf`, `b271310`) only adds a third entry, `200`, to `WV_SCALE_FACTOR` and stores GDAS as kg/m² × 10, which is the same encoding as the VIIRS file (g/cm² × 100). Neither the commit messages nor the code comments mention the factor of 2. The result is that GDAS goes through exactly the same halving as VIIRS, so GDAS and VIIRS stay consistent, but whether this was a conscious choice should be confirmed with Eric.

**Recommendation:** confirm with Eric, then keep it for now and document it. Reproduce it deliberately in the Rust port,
behind a named constant. Fixing it changes HLS reflectance for every ancillary source, so it needs
a separate science decision (and possibly reprocessing).

## 2. GDAS grid is shifted 1° in longitude

**Where:** `PGE11_grib2.c:289` (`reverse_data`) and `grib.c:352`: `jj = j + (nx / 2) - 1;`
together with `lut_subr.c` `get_cart_linesamp()`.

**What happens:**
- GDAS 1° rows start at 0°E. The reader rotates each row so it "starts at −180", but with the `- 1` the rotated row starts at **−179°** (column 0 = −179°, column 359 = +180°).
- `get_cart_linesamp()` assumes column 0 = −180° (`isamp = (lon + 180) / lonres`).
- Net effect: the value used at longitude L is the GDAS value at **L + 1°**, about 111 km at the equator.

**Evidence:** a direct simulation of the two routines (row values set to their own longitude)
returns 16.12 for a query at 15.12, and 1.0 for a query at 0.0.

**How to confirm on real data:** `lut_subr.c` already has a commented-out block that writes the
resampled `wv`/`oz` arrays to `whatever.hdf`. Enable it for one scene, then compare with the same
GDAS file read by an independent tool (e.g. `xarray` + `cfgrib`, or `wgrib2`).

**Suggested fix:** use `jj = j + nx / 2` (assuming the row starts at 0°E), or derive the first
longitude from the GRIB grid definition rather than assuming it.

## 3. Wrong neighbour at the 179°–180° longitude wrap

**Where:** `lut_subr.c:3199-3200` (`get_cart_interp`):
```c
i3 = (iline + 1) * samps + isamp + 1;
i4 = iline * samps + isamp + 1;
```

**What happens:**
- When `isamp` is the last column, `isamp + 1` runs past the end of the row into column 0 of the *next* row, 1° further south, instead of wrapping to column 0 of the same row.
- At the southern edge, `i3` can fall outside the array and is then set to 0 (no error), which pulls the interpolated value toward zero.

**Suggested fix:** wrap the column index explicitly, `isamp1 = (isamp + 1) % samps`, and clamp the
row index at the poles rather than substituting 0.

## 4. Half-cell offset when resampling to the CMG grid

**Where:** `lut_subr.c:3106` (`CMGlat = 90.0 + CMGgeoincr`) and `lut_subr.c:3112`
(`CMGlon = -180.0 - CMGgeoincr`).

**What happens:**
- The CMG grid is "pixel is area": cell centres are at 89.975 / −179.975.
- The loop samples at 90.0 / −180.0, the cell corners, an offset of 0.025° (half a CMG cell).
- GDAS is "pixel is point" (grid points at whole degrees, so its upper-left cell corner is −180.5 / 90.5).

**Impact:** 2.5% of a GDAS cell; negligible for these smooth fields, but worth fixing for correctness.

**Suggested fix:** start at `90.0 - CMGgeoincr / 2` and `-180.0 + CMGgeoincr / 2`.

## 5. SR output encoding changed: negative reflectance becomes fill

**Where:** `output.h:13-16`; used in `output.c:357-362` (XML metadata) and `output.c:729-750`
(conversion to stored values).

**What happens:**
- `SCALE_FACTOR_REFL` / `OFFSET_REFL` changed from `0.0000275 / -0.2` to `0.0001 / 0.0`.
- The output is `UINT16` and the minimum is clipped to 0 (`if (min_value < 0) min_value = 0`).
- With offset 0, any reflectance ≤ 0 (for example water in SWIR, deep shadow) is stored as 0, which is also `FILL_VALUE` (`output.h:8`).
- Before, values down to −0.2 were stored.

This is unrelated to GDAS and affects every SR output (Landsat and Sentinel-2).

**Action:** confirm with Eric whether this was intended for HLS. If not, restore
`0.0000275 / -0.2`. (The 2026 Landsat aerosol-repair review test on Meteor restored the old
encoding with a local patch.)

## 6. Uninitialised `swapped` in the time interpolation

**Where:** `lut_subr.c:2906` (declaration), `lut_subr.c:3027` (the initialisation, commented out),
`lut_subr.c:3064` (`if (swapped == 0)`).

**What happens:**
- `swapped` is read without ever being set. In C this is undefined behaviour.
- If it holds a non-zero value, the two GDAS times receive each other's weights.
- The swap-detection code it belonged to is commented out.

**Suggested fix:** initialise `swapped = 0`, or remove the variable and its branch entirely. LaSRC
already rejects file pairs in the wrong order.

## 7. GRIB run without `--aux2` uses an unset file name

**Where:** `lasrc.c:349` and `lasrc.c:410`; `aux2nm` is passed to `open_gb_oz_wv_datasets()`.

**What happens:** if the first aux file is GRIB and `--aux2` is not given, `aux2nm` is never filled
in, but the GRIB reader still opens it.

**Suggested fix:** make `--aux2` required when the first file is GRIB, and fail with a clear message.

## 8. `--aux2` only works as `--aux2=FILE`

**Where:** `get_args.c:58`: `{"aux2", optional_argument, 0, 'b'}`.

**What happens:** with `optional_argument`, `getopt_long` only accepts `--aux2=FILE`. The form
`--aux2 FILE` leaves the option empty and treats `FILE` as a stray argument.

**Suggested fix:** use `required_argument` (`--aux2` is only given when it has a value).

## 9. Missing checks when a GRIB variable is not found

**Where:** `lut_subr.c` `open_gb_oz_wv_datasets()`: `grib_ret` (e.g. line 2940) is never checked.

**What happens:** if `PWAT` or `TOZNE` is missing, or the level name differs, the data array stays
NULL and is used anyway. The `"entire atmosphere (considered as a single layer)"` level name is also
hard-coded and may change with GFS v17.

**Suggested fix:** check every read and fail with the variable and level name. Log which GRIB
message was used.

## 10. Minor

- **Truncation:** `lut_subr.c:3119, 3122` use `(uint8)(value*0.40)` and `(uint16)(value*10.0)`, which truncate rather than round (bias up to 1 stored unit: 2.5 DU ozone, 0.1 kg/m² water vapor). Use `roundf()`.
- **Format string:** `input.c:1521` `sprintf (temp, gmeta->scene_center_time);` passes data as a format string. Use `snprintf(temp, sizeof temp, "%s", ...)`.
- **Memory:** the static `tmpfltarray` buffers in `open_gb_oz_wv_datasets()` are allocated by the readers and never freed. This is harmless for one scene per process, but matters if the code becomes a library.
- **Debug output:** several unconditional `printf` lines in `open_gb_oz_wv_datasets()` (acquisition time, weights) should go through the normal verbose logging.
- **Commented-out code in `common.h:38-40`:** test aerosol-window settings (`LAERO_WINDOW 1`) are left in a comment. Production values are active; the comment could be removed.
