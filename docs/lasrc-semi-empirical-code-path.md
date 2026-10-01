# LaSRC's semi-empirical code path

A ~10 minute tour of how LaSRC turns top-of-atmosphere (TOA) reflectance into surface
reflectance (SR) on its default, "semi-empirical" path (`use_orig_aero=false`). This is the path
NASA HLS runs in production. It applies to both Landsat and Sentinel-2; differences are noted
where they matter. For variable names (`roatm`, `erelc`, `troatm`, ...) see `docs/GLOSSARY.md`.

### The big picture

The atmosphere changes what the sensor sees in three ways: it adds scattered light (path
reflectance, `roatm`), it attenuates light on the way down and back up (transmission,
`ttatmg`, and absorption by other gases, `tgo`), and it reflects some surface-reflected light
back down (spherical albedo, `satm`). Given those terms, TOA reflectance can be inverted to
surface reflectance assuming a Lambertian surface:

$$
\begin{aligned}
x &= \rho_{\mathrm{TOA}} - t_{go}\,\rho_{\mathrm{atm}} \\
\rho_{\mathrm{sfc}} &= \frac{x}{t_{go}\,T_{\mathrm{atm}} + S_{\mathrm{atm}}\,x}
\end{aligned}
$$

where $\rho_{\mathrm{atm}}$ is `roatm`, $T_{\mathrm{atm}}$ is `ttatmg`, $S_{\mathrm{atm}}$ is
`satm`, and $t_{go}$ is `tgo`.

The terms depend on geometry, gases, surface pressure, and aerosol. LaSRC treats everything but
aerosol as known, from metadata and ancillary data. Aerosol is retrieved from the image itself:
its optical thickness at 550 nm (AOT) and its spectral slope (the Angstrom exponent, `eps`).

The default path makes two simplifications that set it apart from `use_orig_aero=true`:

- **One atmosphere per scene**: water vapor, ozone and pressure come from the scene center and
  are shared by every pixel. Only aerosol varies across the scene.
- **Semi-empirical evaluation**: instead of interpolating the radiative transfer lookup tables
  every time, LaSRC fits each atmospheric term as a cubic polynomial in AOT once per scene,
  then evaluates those polynomials everywhere else.

The steps, in order:

1. Gather inputs: TOA reflectance, fill mask, geometry, and ancillary data
2. Set up the scene-level geometry and atmosphere
3. Fit the atmospheric terms as cubic functions of AOT
4. Compute a quick first-pass SR
5. Retrieve aerosol in small windows
6. Validate each window, retrying failures as water
7. Fill in failed windows and spread aerosol to every pixel
8. Apply the final per-pixel correction

### 1. Inputs

- **TOA reflectance** for the reflective bands. Landsat applies the MTL gains and offsets;
  Sentinel-2 applies `RADIO_ADD_OFFSET` and `QUANTIFICATION_VALUE`, with its 20 m and 60 m
  bands replicated to 10 m. Landsat thermal bands and Sentinel-2 B10 (cirrus) are carried
  through but not corrected.
- **Fill mask**: a pixel is fill if any band is fill.
- **Geometry**: sun and view angles from the product metadata.
- **Daily water vapor and ozone**: VIIRS (or MODIS) on a global 0.05 degree grid (the "CMG").
- **Static ancillary data**:
  - radiative transfer lookup tables for path reflectance, transmission, spherical albedo and
    aerosol extinction, indexed by angles, pressure and AOT
  - a CMG elevation model, converted to surface pressure
  - the **band ratio map**: a climatology, per CMG cell, of the ratios of coastal, blue and
    SWIR surface reflectance to red surface reflectance, as a function of NDWI. This is what
    makes aerosol retrieval possible over land.

### 2. Scene-level geometry and atmosphere

LaSRC picks one geometry and one atmosphere for the whole scene:

- **Geometry**: solar zenith (`xts`, cosine `xmus`), view zenith (`xtv`, `xmuv`), and relative
  azimuth (`xfi`, `cosxfi`). Landsat uses the metadata solar zenith with view zenith and
  relative azimuth fixed at 0; Sentinel-2 uses the scene-mean sun and view angles.
- **Atmosphere**: the scene center is converted to lat/long, and water vapor, ozone and pressure
  are read from the single nearest CMG cell.

### 3. Fitting the atmospheric terms

For each band and each of the 22 AOT levels in the lookup tables (0.01 to 5.0), LaSRC runs the
full lookup-table evaluation, `atmcorlamb2`, with the scene geometry and atmosphere. That gives
22 samples each of `roatm`, `ttatmg` and `satm` as a function of AOT. A least-squares cubic is
fitted to each (`get_3rd_order_poly_coeff`). `roatm` is fitted only up to the AOT where it
stops increasing, and higher AOT values are clamped to that point later. `tgo` does not depend
on AOT and is kept as one value per band.

This is the "semi-empirical" part: the terms come from radiative transfer physics, but from here
on they are applied through fitted curves. Everything after this step that needs "the
atmospheric terms for this band at this AOT" calls `atmcorlamb2_new`. It first converts the
retrieved 550 nm AOT, $\tau_{550}$, to the AOT that band's polynomials are indexed by:

$$
\tau_b = \min\left( \frac{\tau_{550}}{n_b} \left( \frac{\lambda_b}{0.55} \right)^{-\epsilon},\ \tau_{\max,b} \right)
$$

where $\lambda_b$ is the band's center wavelength in micrometers, $\epsilon$ is the Angstrom
exponent (`eps`), $n_b$ is the band's aerosol extinction relative to 550 nm from the lookup
tables (`normext`), and $\tau_{\max,b}$ is the AOT where the `roatm` fit stops (`roatm_upper`).
The $(\lambda_b / 0.55)^{-\epsilon}$ factor is the Angstrom law for the band's aerosol optical
depth; dividing by $n_b$ maps that back onto the lookup tables' AOT axis, which the cubics were
fitted on. It then evaluates the three cubics at $\tau_b$ and applies the Lambertian
correction above.

### 4. First-pass ("climatological") surface reflectance

The aerosol retrieval needs a rough idea of the surface type, so LaSRC first corrects the image
at a fixed, low AOT of 0.05. Landsat divides TOA by the cosine of the per-pixel solar zenith
before this; Sentinel-2 does not.

Only a few bands of this first pass are used later: the NIR and SWIR bands for NDWI at each
aerosol window, and on Sentinel-2 the coastal band for the aerosol QA test.

### 5. Aerosol retrieval in windows

The scene is tiled into small windows: 3x3 pixels around a center for Landsat, 6x6 pixels
anchored at the upper-left corner for Sentinel-2. For each window:

1. **Locate it.** Convert the window position to lat/long.
2. **Expected surface ratios.** From the band ratio map, bilinearly interpolate (between the 4
   surrounding CMG cells) the slope and intercept of the blue, coastal and SWIR ratios. Combine
   them with the first-pass NDWI, clamped to the map's local NDWI range, to get `erelc`: the
   expected surface reflectance of each retrieval band relative to the red band. Where the map
   is unreliable, the slopes are set to 0 and fixed ratios are used.
3. **Observed TOA.** Average the window's TOA in the retrieval bands (`troatm`), skipping fill.
4. **Invert.** `subaeroret_new` searches for the AOT at which the corrected surface reflectances
   best match the expected ratios, and returns that AOT with a fit residual. It does this for
   three values of `eps` (1.0, 1.75, 2.5), fits a parabola through the three residuals to pick
   the best `eps`, and runs a final retrieval at that value.

Why this works: aerosol scattering is much stronger at short wavelengths than in the SWIR.
Correcting with too little AOT leaves the coastal and blue surface reflectances too bright
relative to red; too much makes them too dark. The band ratio map supplies what those
ratios should be for this kind of surface, so the AOT that reproduces them is the estimate. The
residual measures how well any AOT could match them, which is why it is used to judge the
retrieval next.

The window's result is stored at its anchor pixel and copied to the window's other non-fill
pixels.

### 6. Validating windows and the water retest

A retrieval counts as valid land only if:

- its residual is below a threshold that grows with aerosol load (`corf`), and
- the corrected NIR and red reflectances look like vegetation or land: NIR above 0.1 and a
  positive NDVI.

Otherwise the window is retried as water: the same inversion with all band ratios set to 1 and
`eps` fixed at 1.5. A good fit with non-negative coastal reflectance makes it valid water;
anything else marks the window as failed. These outcomes are recorded per window in `ipflag`
(fill, clear, water, failed, ...), which becomes the aerosol QA band.

### 7. Filling in and spreading aerosol

Failed windows need an aerosol value, and every pixel needs a smooth aerosol field.

- **Landsat**: `fix_invalid_aerosols` fills failed windows from nearby valid ones, then
  `aerosol_interp` bilinearly interpolates AOT and `eps` between window centers.
- **Sentinel-2**:
  - `aerosol_interp_sentinel` interpolates AOT from each window toward its right and lower
    neighbors.
  - `ipflag_expand_failed_sentinel` treats pixels within 12 px of a failed window as failed too.
  - `aero_avg_failed_sentinel` replaces each failed pixel with the average of valid pixels in a
    61x61 neighborhood, repeating outward until everything is filled.

The result is an AOT (`taero`) and `eps` (`teps`) value for every non-fill pixel.

### 8. Final per-pixel correction

Every non-fill pixel and band is corrected with `atmcorlamb2_new` at that pixel's AOT and
`eps`, clamped to the valid reflectance range. The outputs are:

- surface reflectance per band
- AOT per pixel
- the aerosol QA band: the window outcome from step 6 plus a low, medium or high aerosol level,
  set from how much the final correction differs from the first-pass SR in a reference band

### How `use_orig_aero=true` differs

For comparison, the alternative path changes steps 2, 3, 5, 6 and 8:

- water vapor, ozone and pressure are interpolated to every pixel from the 4 surrounding CMG
  cells, instead of one scene-center value
- step 3 is skipped; the inversion (`subaeroret`) and the final correction (`atmcorlamb2`)
  evaluate the full lookup tables directly, using each window's or pixel's own atmosphere

Steps 1, 4 and 7 are the same in both paths. Step 4 uses the scene-center atmosphere even
when the flag is on.

### Known issues

Two behaviors in this path look like bugs:

- The Sentinel-2 fill test runs on converted TOA rather than the DN, so a valid DN of exactly
  1000 (TOA 0, possible since processing baseline 04.00 added a -1000 offset) is treated as
  fill ([#28](https://github.com/NASA-IMPACT/espa-surface-reflectance/issues/28)).
- The Sentinel-2 aerosol interpolation in step 7 reads fill pixels as AOT 0, pulling AOT down
  near fill edges ([#29](https://github.com/NASA-IMPACT/espa-surface-reflectance/issues/29)).
