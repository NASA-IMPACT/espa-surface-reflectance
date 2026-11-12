# Known issues in LaSRC and the GDAS changeset

Line numbers refer to this branch (`feature/gdas-integration`). Each entry can be copied into a
GitHub issue.

- **Issues 1–10** were found while reviewing the GDAS changeset (`b271310`).
- **Issues 11–18** were reported by J. Ju (LaSRC code review, Nov 2026). They are **not** caused by GDAS: they are already present in release v3.5.1.0, the version HLS runs today, and affect VIIRS and GDAS runs equally. They are listed here so the Rust port can decide, for each one, whether to reproduce or fix it.

**Status: drafts.** Nothing here has been fixed on this branch. Items marked "changes output" will
shift GDAS results slightly once fixed, so each fix should come with a before/after comparison on a
small reference set (see README section 3).

| # | Issue | Severity | Fixing changes output? |
|---|---|---|---|
| 1 | Water vapor is half its true value (same as the existing VIIRS/MODIS path) | Known behaviour (inherited) | Yes, for all sources; science decision |
| 2 | GDAS grid is shifted 1° in longitude | High | Yes (GDAS) |
| 3 | Wrong neighbour at the 179°–180° longitude wrap | Medium | Yes (GDAS, ~1° strip) |
| 4 | Half-cell (0.025°) offset when resampling to the CMG grid | Low | Yes (GDAS, negligible) |
| 5 | SR output encoding changed to 0.0001/0: negative reflectance becomes fill | High | Yes (all sources) |
| 6 | Uninitialised `swapped` in the time interpolation | High (undefined behaviour) | Possibly (GDAS) |
| 7 | GRIB run without `--aux2` uses an unset file name | Medium | No (crash/robustness) |
| 8 | `--aux2` only works as `--aux2=FILE` | Low | No |
| 9 | Missing checks when a GRIB variable is not found | Medium | No (robustness) |
| 10 | Minor: truncation, format string, memory | Low | Truncation only |
| 11 | Landsat green band wavelength is 0.585 µm, OLI band 3 is centred near 0.561 µm | Open science question | Yes (L30 green) |
| 12 | View-zenith transmittance uses the wrong LUT angle grid | Medium | Yes (mainly S30) |
| 13 | Landsat aerosol interpolation reads a window centre outside the image | Medium (out-of-bounds read, rare) | Yes (rare L30 edge rows/columns) |
| 14 | Quadratic coefficient for the Angstrom fit is truncated to an integer | Medium | Yes (S30 and L30) |
| 15 | TOA and intermediate SR clamped to [−0.2, 1.6] without a flag | Design defect | Yes (bright / high-sun-zenith pixels) |
| 16 | CMG row/column clamped one past the end of the array | Low (unreachable in practice) | No |
| 17 | Pixel-centre line offset uses `i - 0.5` instead of `i + 0.5` | Low (negligible) | Negligible |
| 18 | Landsat SR assumes a nadir view for the atmospheric terms | Design limitation, under investigation | Yes (L30, off-nadir pixels) |

---

## 1. Water vapor is half its true value (same as the existing VIIRS/MODIS path)

**Where:** `common.h:301` (`WV_SCALE_FACTOR {200, 200, 200}`), `lut_subr.c:3122`
(`wv[l] = (uint16)(value*10.0)`), `compute_refl_subr.c:248` (`uwv = wv / WV_SCALE_FACTOR`).

**What happens:**

GDAS `PWAT` is in kg/m², and 1 g/cm² = 10 kg/m², so the correct conversion is kg/m² ÷ 10. The GDAS
code stores kg/m² × 10 and later divides by `WV_SCALE_FACTOR` = 200, which gives kg/m² ÷ 20: half the
true value.

VIIRS does exactly the same. In the VIIRS ancillary file (checked on `VJ104ANC.A2025011`) the
`Coarse Resolution Water Vapor` variable is `uint16` with `scale_factor = 0.01` and `units = g/cm^2`,
so the true value is the stored number ÷ 100 (a stored 145 means 1.45 g/cm²). LaSRC reads the stored
numbers unscaled (`H5Dread(..., H5T_NATIVE_USHORT, ..., wv)` in `lut_subr.c`) and then divides by 200,
so it uses 0.725 g/cm² instead of 1.45.

| | Value in the internal `wv` array | Divided by | Result | True value | Ratio |
|---|---|---|---|---|---|
| VIIRS | g/cm² × 100 (as stored in the file) | 200 | g/cm² ÷ 2 | g/cm² | ½ |
| GDAS | kg/m² × 10 = g/cm² × 100 | 200 | g/cm² ÷ 2 | g/cm² | ½ |

The GDAS `× 10` puts GDAS water vapor into the same units as the VIIRS file, so both go through the
same ÷ 200 and both end up at half the true value. This is why switching from VIIRS to GDAS does not
change HLS outputs through this route.

**History:**

1. **MODIS:** USGS operations already divided MODIS water vapor by 200 instead of 100. The commented-out lines in `common.h` show this was known (`{200, 100} /** USGS OPS with MODIS bug **/`, and a fixed `{100, 100}` that was never switched on).
2. **Jan 2023 (`18875d4`):** VIIRS support was added with the correct divisor, 100.
3. **Aug 2023 (`3fc9b0a`, `ebea983`, G. Schmidt):** VIIRS was changed to 200 as well, "reintroduced the WV scale factor bug for both MODIS and VIIRS", so the MODIS-to-VIIRS switch would not change operational outputs. The `common.h` comment reads "USGS OPS with MODIS bug, and added bug for VIIRS". This is the behaviour of release v3.5.1.0, the version HLS runs today.
4. **2026 (`a64c9cf`, `b271310`, E. Vermote):** the GDAS changeset adds a third entry, `200`, for GRIB, so GDAS inherits the same halving. Neither the commit messages nor the code comments mention the factor of 2, so whether this was a conscious choice for GDAS should be confirmed with Eric.

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

---

# Issues reported by J. Ju (pre-existing in LaSRC v3.5.1.0)

Source: J. Ju, "LaSRC Review" slides (Nov 2026). Each item was checked against the code on this
branch; notes added in that check are marked **Check:**.

## 11. Landsat green band wavelength is 0.585 µm, OLI band 3 is centred near 0.561 µm

**Where:** `lut_subr.c:24-25`
```c
float landsat_lambda[NREFLL_BANDS] = {0.443, 0.480, 0.585, 0.655, 0.865, 1.61, 2.2};
```
used in `lut_subr.c:108` to scale the aerosol optical thickness to each band:
`mraot550nm = raot550nm / normext * pow(lambda[iband] / 0.55, -eps)`.

**Evidence:**
- **Rayleigh check (points to 0.561 µm):** the green Rayleigh optical thickness in LaSRC's `tauray-ldcm.ASC` is 0.09070. The standard approximation `tauray = 0.008569 λ⁻⁴ (1 + 0.0113 λ⁻² + 0.00013 λ⁻⁴)` gives 0.0897 at 0.561 µm and 0.0757 at 0.585 µm. **Check:** recomputed, same values.
- **Aerosol LUT check (points to 0.585 µm):** the Angstrom exponent derived from `TRANSLUTV3.0-ldcmb3-URBANCLEAN-V2.0.ASCII` looks reasonable with 0.585 µm.

The two lookup tables seem to have been built with different green wavelengths. For comparison, the
Sentinel-2 green band uses 0.560 µm (`sentinel_lambda`, `lut_subr.c:29-35`).

**Impact:** L30 green band only, through the band-scaled AOT. For an Angstrom exponent of 1.5,
using 0.585 instead of 0.561 lowers the green-band AOT by about 6%.

**Action:** ask Eric which wavelength each LUT was generated for before changing anything.

## 12. View-zenith transmittance uses the wrong LUT angle grid

**Where:** `lut_subr.c:321-322` (in `atmcorlamb2`)
```c
comptrans (ip1, ip2, iaot1, iaot2, xtv, mraot550nm, iband, pres, tpres,
    aot550nm, transt, xtvstep, xtvmin, tts, &xttv);
```
with `xtvmin = 2.84090`, `xtvstep = 3.68` set in `compute_refl_subr.c:160-161`.

**What happens:**
- `xtvmin`/`xtvstep` describe the view-angle grid of the atmospheric-reflectance LUT.
- The transmittance LUT (`TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII`) uses one grid for both sun and view angles: 0° to 80° in 4° steps (`xtsmin = 0`, `xtsstep = 4`).
- `comptrans()` picks the row with `(xtv - 2.84) / 3.68` but interpolates with the 4° grid (`xmts = (xtv - tts[its]) * 0.25`).
- **Check:** for view zeniths ≥ 2.84° the code therefore extrapolates from the wrong pair of rows. For example, at 6° it uses the 0° and 4° rows with weight 1.5 instead of interpolating between 4° and 8°. Below 2.84° the result is correct.

**Impact:**
- **Sentinel-2:** affected, because it uses per-pixel view zenith (up to ~11°). The error is small, since transmittance changes slowly with angle near nadir.
- **Landsat:** no impact by default, because Landsat SR uses a view zenith of 0. **Check:** with `--use_orig_aero_alg`, the Landsat aerosol retrieval uses per-pixel view zenith and would be affected.

**Suggested fix:** pass `xtsstep, xtsmin` (4, 0) for the view-angle call as well.

## 13. Landsat aerosol interpolation reads a window centre outside the image

**Where:** `aero_interp.c:97` and `aero_interp.c:147-148` (`aerosol_interp_landsat`)
```c
center_line = (int) (line * aero_step) * aero_window + half_aero_window;
center_samp = (int) (samp * aero_step) * aero_window + half_aero_window;
```

**What happens:**
- Landsat retrieves aerosol on a 3 × 3 window grid, at centres 1, 4, 7, … (`compute_landsat_refl.c:981-999`).
- When `nlines % 3 == 1` (or `nsamps % 3 == 1`), the last line (or column) maps to a window centre at index `nlines` (or `nsamps`), outside the image, where no aerosol was retrieved.
- The interpolation uses that value anyway.
- **Check:** for the last line this reads past the end of the `taero` array (an out-of-bounds memory read, not just a stale value). For the last column it reads the first pixel of the next row.

**Impact:** L30 only, and rare: it needs the last line or column of the scene to contain valid
(non-fill) data. Sentinel-2 uses a different routine and is not affected.

**Suggested fix:** clamp the window centre to the last centre inside the image.

## 14. Quadratic coefficient for the Angstrom fit is truncated to an integer

**Where:** `compute_landsat_refl.c:968` and **Check:** also `compute_sentinel_refl.c:1036`
```c
xd = (int) ((eps2 * eps2) - (eps3 * eps3));
```

**What happens:**
- With `LOW_EPS = 1.0`, `MOD_EPS = 1.75`, `HIGH_EPS = 2.5`, the true value is −3.1875, but the cast stores −3.
- `xd` is one coefficient of the quadratic used to find the Angstrom exponent with the smallest residual (`coefa`, `coefb`, `compute_landsat_refl.c:1374-1375`).
- The neighbouring coefficients (`xa`, `xb`, `xe`) are not truncated.

**Impact:** shifts the retrieved Angstrom exponent, and hence the aerosol model, for both L30 and S30.

**Suggested fix:** remove the `(int)` cast in both files.

## 15. TOA and intermediate SR clamped to [−0.2, 1.6] without a flag

**Where:**
- TOA: `compute_landsat_refl.c:150-158` (`compute_landsat_toa_refl`)
- SR: `compute_landsat_refl.c:686-689` and `:1784-1787`; **Check:** also `compute_sentinel_refl.c:1910-1913`

**What happens:**
- TOA reflectance is divided by cos(SZA) and then clamped to `[MIN_VALID_REFL, MAX_VALID_REFL]` = [−0.2, 1.6].
- From the Landsat input range [−0.09998, 1.2107], the value can exceed 1.6 for SZA > ~41° and fall below −0.2 for SZA > ~60°.
- The clamped values are not flagged, and they feed the aerosol retrieval and the SR calculation.
- Intermediate SR is clamped the same way.
- Only the final output needs clamping. Also, the code comment says "scaled" where it should say "unscaled".

**Impact:** bright targets (snow, cloud edges) and high-sun-zenith scenes.

**Suggested fix:** keep full-range values internally, clamp only when writing output, and flag clamped pixels.

## 16. CMG row/column clamped one past the end of the array

**Where:** `compute_landsat_refl.c:841-847`, `compute_refl_subr.c:236-242`, `compute_sentinel_refl.c:909-915`
```c
else if (lcmg >= CMG_NBLAT) lcmg = CMG_NBLAT;   /* should be CMG_NBLAT - 1 */
else if (scmg >= CMG_NBLON) scmg = CMG_NBLON;   /* should be CMG_NBLON - 1 */
```

**Impact:** none in practice. The index only reaches the limit for latitudes south of −89.975° or
an invalid longitude. It would be an out-of-bounds read if it ever happened.

**Check:** other places in the same files already use the correct `CMG_NBLAT - 1` / `CMG_NBLON - 1`
(e.g. `compute_landsat_refl.c:1072-1077`).

**Suggested fix:** use `- 1` everywhere.

## 17. Pixel-centre line offset uses `i - 0.5` instead of `i + 0.5`

**Where:** `compute_landsat_refl.c:803` and `:1033`; **Check:** also `compute_sentinel_refl.c:868` and `:1082`
```c
img.l = i - 0.5;   /* strictly should be i + 0.5, matching img.s = j + 0.5 */
```

**Impact:** negligible. The position used for the ancillary lookup is shifted by at most one
pixel in the line direction (≤ 30 m for Landsat), far smaller than the 0.05° (~5 km) CMG cell.

**Suggested fix:** `img.l = i + 0.5`.

## 18. Landsat SR assumes a nadir view for the atmospheric terms (under investigation)

**Where:** `compute_landsat_refl.c:589-600` (view zenith initialised to 0, used as `xtv_center` in
`atmcorlamb2`). Sentinel-2 uses per-pixel view angles.

**What happens:** the Rayleigh (molecular) scattering term depends on the view angle through the
scattering phase function. J. Ju's comparison of the Rayleigh phase function at the observed angle
versus nadir (solar zenith 50°) gives ratios of about:
- 1.03–1.04 for backward views of −4° to −8°
- about 0.94 at +8° for L8/L9
- about 0.92 at +10° for S2.

Treating every Landsat pixel as nadir therefore misstates the Rayleigh contribution toward the
swath edges.

**Impact:** L30 only, growing toward the swath edges and strongest in the blue bands, where Rayleigh
scattering dominates. J. Ju listed this among the points to raise with Eric.

**Action:** discuss with Eric. A fix would mean using per-pixel view angles for Landsat, as is done
for Sentinel-2.
