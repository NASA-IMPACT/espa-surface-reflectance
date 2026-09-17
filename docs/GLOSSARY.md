# LaSRC glossary

The C variable names carried into the Rust port are terse and mostly undocumented.
This is what they mean, with the C declaration to check against. File references are
to `lasrc/c_version/src/`.

## Naming conventions

The names come from the 6S radiative transfer code (Fortran heritage, French
authors), so they are built from a few prefixes rather than words.

| Prefix | Meaning | Examples |
|---|---|---|
| `ro` | rho, reflectance | `rotoa` (rho TOA), `roatm` (rho atmosphere), `roslamb` (rho surface, Lambertian), `ros1`/`ros4`/`ros5` (rho surface per band), `rolutt` (rho LUT), `xrorayp` (rho Rayleigh path) |
| `x` | plain variable prefix | `xts`/`xtv` (theta sun/view), `xmus`/`xmuv` (mu = cosine of theta), `xfi` (phi, azimuth), `xlntau` (ln tau) |
| `t` | transmission, or table | `tgo`/`tgoz`/`tgwv` (transmission of other gases/ozone/water vapor), `ttatmg` (total transmission atmosphere + gas); but `tts`/`ttv`/`tpres` are tables |
| `s` | spherical albedo | `satm`, `sphalbt` |
| `tau` | optical depth | `tauray` (Rayleigh), `raot550nm` (aerosol at 550 nm) |

Leftover French: `nbfi` is "nombre de phi" (number of azimuth angles), `nbfic`
its cumulative form, and `cfonc1/2/3` in `local_chand` are "fonction"
coefficients.

Two names mislead: `troatm` holds TOA reflectance (a table of rho TOA) even
though its C comment says atmospheric reflectance, and `xndwi` is not the
textbook NDWI (see below).

## Angles and geometry

| Name | Meaning |
|---|---|
| `xts` | Solar zenith angle, degrees. For Landsat the scene-center coefficients use the metadata value (90 - MTL `SUN_ELEVATION`); the per-pixel TOA step uses the SZA band. |
| `xtv` | View (observation) zenith angle, degrees. Hard-coded to 0.0 for Landsat in `init_sr_refl` (`compute_refl_subr.c`); read from metadata for Sentinel-2. |
| `xmus`, `xmuv` | `cos(xts)`, `cos(xtv)`. |
| `xfi`, `cosxfi` | Relative azimuth between sun and view, degrees, and its cosine. 0.0 for Landsat. |
| `scaa` / `cscaa` | Scattering angle (degrees) and its cosine, from `xmus`, `xmuv`, `cosxfi`. |
| `its`, `itv` | Indices into the solar / view angle LUT grids. |
| `img.l`, `img.s` | Fractional image line/sample fed to ESPA `from_space`. The aerosol window uses `(i - 0.5, j + 0.5)`. |

## Reflectance chain

| Name | Meaning |
|---|---|
| `rotoa` | Top-of-atmosphere reflectance, already divided by `cos(SZA)`. |
| `roslamb` | Lambertian surface reflectance, the output of an atmospheric correction call. |
| `roatm` | Intrinsic atmospheric (path) reflectance. |
| `ttatmg` | Total atmospheric transmission, including gas absorption. |
| `satm` | Atmospheric spherical albedo. |
| `tgo` | Other-gas transmittance, `tgog * tgoz`. |
| `tgog`, `tgoz`, `tgwv`, `tgwvhalf` | Transmittance for other gases, ozone, water vapor, and water vapor over the solar half-path. |
| `xrorayp` | Rayleigh (molecular) path reflectance, from `local_chand`. |
| `xtaur` | Rayleigh optical depth scaled to surface pressure. |
| `sband` | Working per-band reflectance array: TOA first, then climatological SR, then final SR (C overwrites it in place). |
| `btgo`, `broatm`, `bttatmg`, `bsatm` | Per-band scene-center values at the climatological AOT (0.05) used for the first-pass correction. |
| `troatm` | Per-band TOA reflectance handed to the aerosol retrieval. |
| `aerob1/2/4/5/7` | Saved TOA reflectance for the bands the retrieval needs. |

## Aerosol retrieval

| Name | Meaning |
|---|---|
| `raot`, `raot550nm` | Aerosol optical thickness at 550 nm (the quantity being retrieved). |
| `mraot550nm` | AOT scaled to the band wavelength by the Angstrom exponent, clamped to `roatm_upper`. |
| `eps` | Angstrom exponent (spectral dependence of AOT). `LOW_EPS` 1.0, `MOD_EPS` 1.75, `HIGH_EPS` 2.5, `WATER_EPS` 1.5. |
| `epsmin` | Exponent minimising the residual, from a parabola through the three trial eps values. |
| `taero`, `teps` | Per-pixel retrieved AOT and eps grids. |
| `erelc` | Expected per-band reflectance ratio relative to the reference band (red). `-1.0` means "band unused"; water sets all active bands to 1.0. |
| `xndwi` | Water index `(NIR - SWIR2/2) / (NIR + SWIR2/2)` from climatological SR, clamped to `andwi +/- 2*sndwi`. Drives `erelc`. |
| `ros1`, `ros4`, `ros5` | Surface reflectance of the reference band, red, and NIR during validation. |
| `residual` | Fit quality being minimised: RMS of `roslamb - erelc[ib]*ros1` over the active bands (for water, of `roslamb` itself). |
| `corf` | `raot / xmus_center`, the aerosol impact term that scales the validation thresholds. |
| `ipflag` | Per-pixel QA bits: fill, valid retrieval (`IPFLAG_CLEAR`), water, window-interpolated, plus the two aerosol-level bits. |
| `iaots` | Starting AOT index carried between retrieval calls. |

## Lookup tables (6S)

| Name | Meaning |
|---|---|
| `rolutt` | Intrinsic reflectance table, `[band][pressure][AOT][solar geometry]`. |
| `transt` | Transmission table, `[band][pressure][AOT][sun angle]`. |
| `sphalbt` | Spherical albedo table, `[band][pressure][AOT]`. |
| `normext` | Aerosol extinction normalised at 550 nm, same shape as `sphalbt`. |
| `tsmax`, `tsmin` | Maximum / minimum scattering angle per `[view][solar]` grid cell. |
| `nbfi`, `nbfic` | Number of azimuth angles per grid cell, and its cumulative sum. |
| `tts`, `ttv` | Sun angle and view angle tables. |
| `indts` | Offsets into `rolutt` per sun angle. |
| `tpres`, `aot550nm` | Pressure levels and the 22 AOT grid values the LUTs are sampled on. |
| `roatm_coef`, `ttatmg_coef`, `satm_coef` | Cubic fits of `roatm`/`ttatmg`/`satm` versus AOT at scene-center geometry, used by the fast path. |
| `roatm_iaMax` / `roatm_upper` | Last AOT index where `roatm` still increases; the fit and the AOT clamp use it. |
| `normext_p0a3` | `normext[band][pressure 0][AOT index 3]`, the reference extinction for Angstrom scaling. |

## Ancillary inputs

| Name | Meaning |
|---|---|
| `uoz`, `uwv` | Total column ozone (cm-atm) and water vapor (g/cm^2). |
| `pres`, `tp` | Surface pressure from the CMG DEM, and its per-pixel grid. |
| `atm_pres` | `pres / 1013`, the pressure ratio used for Rayleigh and gas terms. |
| `twvi`, `tozi` | Per-pixel interpolated water vapor and ozone (original aerosol path only). |
| `lcmg`, `scmg`, `xcmg`, `ycmg` | Line/sample and fractional position in the 0.05 degree CMG climatology grid. |
| `andwi`, `sndwi` | Mean and standard deviation of the NDWI climatology, used to clamp `xndwi`. |
| `ratiob1/2/7`, `slpratiob*`, `intratiob*` | Band-ratio climatology and its per-cell slope/intercept versus NDWI. |

## Functions

| Name | Meaning |
|---|---|
| `atmcorlamb2` | Full atmospheric correction: interpolates the LUTs for one band and geometry. |
| `atmcorlamb2_new` | Fast path: evaluates the pre-fitted cubics instead of the LUTs. |
| `subaeroret` / `subaeroret_new` | Aerosol inversion: sweeps the AOT grid, then refines with a parabola. `_new` uses the fitted coefficients. |
| `comproatm`, `comptrans`, `compsalb`, `comptg` | LUT interpolation for `roatm`, transmission, spherical albedo, and gas transmittance. |
| `local_chand` | Rayleigh reflectance (Chandrasekhar). |
| `aerosol_interp_landsat` | Bilinear spread of window-center AOT to every pixel. |
| `fix_invalid_aerosols_landsat` | Replaces failed retrievals with a local average of valid neighbours, in three passes. |
| `find_closest_non_fill` | Finds the nearest non-fill pixel when a window center is fill. |
| `readluts` | Loads all the 6S tables. |

## Abbreviations

| Short | Long |
|---|---|
| AOT | Aerosol optical thickness |
| TOA / SR | Top of atmosphere / surface reflectance |
| CMG | Climate Modeling Grid (0.05 degree global grid) |
| LUT | Lookup table |
| NDWI | Normalized difference water index (here the `xndwi` variant above) |
| eps | Angstrom exponent |
| ESPA | EROS Science Processing Architecture (the XML + raw binary format) |
| GCTP | General Cartographic Transformation Package (the projection library) |
| 6S | Second Simulation of a Satellite Signal in the Solar Spectrum (source of the LUTs) |
