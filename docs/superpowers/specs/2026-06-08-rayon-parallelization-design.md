# Rayon Parallelization Design

## Goal

Add Rayon-based parallelism to the hot loops in both Landsat and Sentinel-2 surface reflectance pipelines, targeting reasonable speedup on multi-core machines. The C code uses OpenMP for these same loops.

## Acceptance Criteria

- Landsat and Sentinel output remain within existing comparison thresholds (no new regressions beyond current C-vs-Rust diffs)
- Optional `num_threads` parameter controls thread count; default uses all available cores
- No `unsafe` code required

## Architecture

Row-chunked `par_iter` (Rayon Approach A). Parallelize by rows using `par_chunks_mut` on output arrays or by collecting independent work items into a Vec and using `par_iter`. Each thread gets a contiguous block of work with no shared mutable state.

Rayon is already a dependency in `lasrc_core/Cargo.toml`.

## Thread Pool Configuration

Add an optional `num_threads: Option<usize>` parameter to:
- `compute_surface_reflectance()` (Rust Landsat entry point in correction.rs)
- PyO3 binding: `process_surface_reflectance()` (lasrc_py/src/lib.rs)
- `compute_sentinel_surface_reflectance()` (Rust Sentinel entry point in correction.rs)
- PyO3 binding: `process_sentinel_surface_reflectance()` (lasrc_py/src/lib.rs)

If provided, build a Rayon `ThreadPool` via `ThreadPoolBuilder::new().num_threads(n).build()` and run parallel sections inside `pool.install(|| ...)`. If `None`, use Rayon's default global pool (all available cores).

The PyO3 Python bindings pass through an optional `num_threads` keyword argument.

## Loops to Parallelize

### 1. Aerosol Window-Center Retrieval (highest impact)

**Landsat:** correction.rs Step 5 (~lines 543-809)
**Sentinel:** correction.rs Step 7 (~lines 1135-1470)

**Current structure:** Nested `while` loops stepping by `aero_window` over lines and samples. Each window center runs 3-4 `subaeroret_new` calls (the heaviest per-pixel computation). Writes to `taero[pix]`, `teps[pix]`, `ipflag[pix]` at non-overlapping window-center indices.

**Strategy:** Collect-and-scatter.
1. Pre-collect all valid window-center coordinates into a `Vec<(usize, usize)>` (line, samp pairs), skipping fill pixels.
2. Use `par_iter` over this Vec. Each iteration computes aerosol retrieval and returns `(pix_index, taero_val, teps_val, ipflag_val)`.
3. Scatter results back into `taero`, `teps`, `ipflag` in a serial pass.

**Read-only inputs:** `sband`, `aux`, `atm_coeff`, `toa_bands`, `tgo_arr`, `normext_p0a3`, `roatm_coef`, `ttatmg_coef`, `satm_coef`, `roatm_ia_max`, `qa_flat`/`qaband`, `solar_zenith`, `space_def`.

**Why collect-and-scatter:** Avoids `unsafe` code entirely. The window-center count is small (scene_pixels / aero_window^2), so the Vec allocation is negligible.

### 2. Final Per-Pixel Atmospheric Correction (second highest impact)

**Landsat:** correction.rs Step 7 (~lines 853-910)
**Sentinel:** correction.rs Step 9 (~lines 1477-1535)

**Current structure:** Nested `for` loops over all pixels. Each pixel calls `atmcorlamb2_new` per band. Writes to `sr_f32[band][(line, samp)]` and sets aerosol QA bits on `ipflag` for the coastal band.

**Strategy:** Split into parallel SR computation + serial QA post-pass.
1. Restructure the loop to iterate over rows. Use `par_chunks_mut` to give each thread a contiguous block of row indices.
2. Each thread computes SR for all bands for its rows, writing to `sr_f32[band][(line, samp)]`. No `ipflag` mutation in this pass.
3. After the parallel loop completes, a serial post-pass over all pixels compares `sband[coastal][pix]` vs `sr_f32[coastal][(line, samp)]` and sets the aerosol QA bits on `ipflag`.

**Read-only inputs:** `taero`, `teps`, `sband`, `atm_coeff`, `tgo_arr`, `normext_p0a3`, `btgo`, `broatm`, `bttatmg`, `bsatm`, `qa_flat`/`qaband`, `toa_bands`.

**Note on Sentinel:** The Sentinel final correction loop iterates band-first (`for iband in 0..nbands { for iline... }`). This can stay band-first with `par_chunks_mut` on the inner row loop, or be restructured to pixel-first. Band-first with parallel rows is simpler and avoids restructuring.

### 3. Climatological SR (moderate impact)

**Landsat:** correction.rs Step 4 (~lines 500-530)
**Sentinel:** correction.rs Step 6 (~lines 1065-1123)

**Current structure:** Nested `for` loops over lines/samples. Each pixel divides TOA by cos(SZA) and applies climatological atmospheric correction. Writes to `sband[band][pix]`. Landsat also marks fill pixels in `ipflag`.

**Strategy:** Pre-pass fill marking + parallel SR.
1. **Pre-pass (serial):** Iterate over all pixels, set `ipflag[pix] = 1 << IPFLAG_FILL` for fill pixels. This is already partially done in the existing code and is cheap.
2. **Parallel SR (Landsat):** Use `par_chunks_mut` over row ranges, same pattern as the final correction loop. Each thread processes a range of rows, reading `qa_flat`, `solar_zenith`, `toa_bands` and writing to `sband[band][pix]`. Skip fill pixels by reading `qa_flat` directly (no `ipflag` dependency).
3. **Parallel SR (Sentinel):** Same approach but no SZA division. The band-first outer loop stays; the inner pixel loop parallelizes over rows.

**Read-only inputs:** `qa_flat`/`qaband`, `solar_zenith`, `toa_bands`, `btgo`, `broatm`, `bttatmg`, `bsatm` (or per-band atmcorlamb2 results for Sentinel).

### 4. Aerosol Interpolation (deferred)

**Functions:** `aerosol_interp`, `aerosol_interp_sentinel` in aerosol.rs

**Not parallelized in this iteration.** The bilinear interpolation reads from 4 neighboring window-center values and writes to the current pixel. Since a row can read from a different row that another thread might be writing, this is not race-free without restructuring (e.g., double-buffering or read-only snapshot). The computation per pixel is cheap (one bilinear interpolation), so the payoff is low.

Can revisit if profiling shows this is a bottleneck.

## Loops That Stay Serial

- **Aerosol post-processing** (`fix_invalid_aerosols`, `ipflag_expand_failed`, `aero_avg_failed`): Iterative neighbor-dependent passes. Complex to parallelize, low cost.
- **Output scaling** (Step 8/10): Already uses `Array2::from_shape_fn` which could potentially use `par_bridge`, but the per-pixel cost is trivial (one multiply + round).
- **LUT precomputation** (Steps 1-3): Runs once per scene, negligible cost.

## Implementation Order

1. **Thread pool plumbing** — Add `num_threads` parameter, build pool, wire through to both pipelines
2. **Aerosol window-center retrieval** — Highest impact, collect-and-scatter pattern
3. **Final per-pixel correction** — Split QA post-pass, `par_chunks_mut` over rows
4. **Climatological SR** — Pre-pass fill marking, parallel row iteration
5. **Verify** — Rebuild, regenerate both Landsat and Sentinel outputs, confirm comparison thresholds hold

## Testing Strategy

- Run existing `compare_outputs.py` for both Landsat and Sentinel after each parallelized loop to catch regressions
- Test with `num_threads=1` to verify serial-equivalent output (should be bit-identical to current)
- Test with default threads to verify parallel output is within thresholds
- Run `cargo test -p lasrc_core` to verify unit tests still pass
