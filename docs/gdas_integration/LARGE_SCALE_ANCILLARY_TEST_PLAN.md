# Plan: large-scale test of low-latency ancillary sources in LaSRC

Follow-up to the 3-granule smoke test in
[ANCILLARY_SOURCE_COMPARISON.md](ANCILLARY_SOURCE_COMPARISON.md). This test reruns the HLS
validation granule set (S30 and L30) with each low-latency ancillary scenario and compares the
results with the validated GDAS-analysis run.

**Status:** plan; phase 0 (inventory) mostly done, see [section 7.1](#71-phase-0-results). Decisions
are listed in [section 8](#8-decisions).

## 1. Goal

Answer, on the full validation set, whether replacing the GDAS analyses with a lower-latency
NOAA source changes LaSRC surface reflectance by a meaningful amount:

- per granule and band, the difference from the GDAS-analysis baseline (V0)
- pass rates against the Ju et al. (2025) Table 7 thresholds, immediate and final stages
- where and when differences occur (season, latitude, water vapor amount)

## 2. Scenarios

Same definitions as the smoke test. C1 is the 6-hourly cycle before the acquisition,
C2 = C1 + 6 h and CP = C1 − 6 h. All files are NOAA NODD 1° `pgrb2.1p00`.

| Run | File valid at C1 | File valid at C2 | Ancillary ready after acquisition | New runs needed |
|---|---|---|---|---|
| **V0** baseline | GDAS C1 f000 | GDAS C2 f000 | 7–13 h | **No**: reuse the validated GDAS outputs |
| **S1** | GDAS C1 f000 | GDAS C1 f006 | 1–7 h | Yes |
| **S2** | GFS C1 f000 | GFS C1 f006 | 0–3.6 h | Yes |
| **S3** | GFS CP f006 | GFS CP f012 | 0 (before the overpass) | Yes |

Forecast files are relabelled so that reference time = valid time and forecast hour = 0, because
LaSRC reads only the reference time. The code is unchanged (`b271310`) and runs in the same HLS
containers used for the GDAS-vs-VIIRS validation.

## 3. Granule set

| Product | Granules | Dates | Inputs | V0 (GDAS) outputs |
|---|---|---|---|---|
| S30 | 605 SAFE IDs (522 granules after twins) | Jan and Jul 2025 | Sentinel-2 L1C SAFE archive used for the validation | Validated S30 GDAS workspace: 559 product folders, 2.8 TB (to be matched against the granule list) |
| L30 | 614 | Jan and Jul 2025 | Landsat L1 inputs used for the validation | Landsat `b271310` workspace, `outputs/grib/`: **614 granules, 2.8 TB** (VIIRS runs in `outputs/viirs/`, also 614) |

The Landsat wrapper used for the validation (`do_lasrc_landsat.py`, workspace copy) already selects
a `gdas1.PGrbF00.*` pair from `LADS/` and accepts `LASRC_AUX_MODE=grib`, so the same per-run
`LADS/` staging works for L30.

## 4. Ancillary data

| Data | Needed for | Status | Size |
|---|---|---|---|
| GDAS f000, Jan and Jul 2025 | V0, S1 (C1 file) | **Already mirrored** for the validation | — |
| GDAS f006 | S1 | Download | ~264 files |
| GFS f000, f006 | S2 | Download | ~528 files |
| GFS f006, f012 (from CP) | S3 | Download | ~528 files (shared with S2 where cycles coincide) |

- Cycles: every 00/06/12/18Z cycle in Jan and Jul 2025, plus 31 Dec 2024, 1 Feb, 30 Jun and 1 Aug
  for cycles that cross a day boundary: about 66 days × 4 cycles.
- Only PWAT and TOZNE are downloaded (`.idx` byte ranges, about 150 KB per file): **about 1,050
  files, ~160 MB**, or ~330 MB with relabelled copies.
- Stored once in a shared **ancillary pool** by cycle. Each run's `LADS/` links to the pool rather
  than copying files.

## 5. Storage and compute estimates

Measured: one S30 run in debug mode is 6.1 GB (smoke test: final COGs 0.25 GB, immediate LaSRC
`.img` bands ~3 GB, intermediate HDFs ~2.8 GB). One L30 run in debug mode is ~4.6 GB (existing GDAS
runs: 2.8 TB / 614). Trimmed and final-only sizes for L30 are estimates.

| Output kept per run | S30 (~1,815 new runs) | L30 (~1,842 new runs) | Total |
|---|---|---|---|
| Full debug (current) | ~11 TB | ~8.5 TB | ~20 TB |
| **Final + immediate (chosen)**: final products, immediate LaSRC SR bands, logs, ancillary record | ~3.6–5 TB | ~2.8–3.5 TB | **~6.5–8.5 TB** |
| Final products only | ~0.45 TB | ~0.7 TB | ~1.2 TB |

Available: **139 TB free** on the HPC scratch file system (2% used), no per-user quota reported.
Storage does not constrain the retention choice.

Plus ~0.3 TB scratch at 30 concurrent jobs. "Trimmed" keeps both validation stages.

**Compute:** S30 ~35 min per run (4 CPUs, 32 GB requested, ~7 GB peak), so ~1,060 job-hours,
~36 h wall at the 30-job limit. L30 is expected to be similar or shorter. About 3 days in total if
every granule is run.

## 6. Screening step (recommended before running LaSRC at scale)

LaSRC applies one water vapor and one ozone value per scene: a single 0.05° CMG cell, stored as
integers in 0.1 kg/m² and 2.5 DU steps. The smoke test showed that an independent script
reproduces these values from the GRIB files (11 of 12 runs).

For every granule and scenario, the screening step computes the stored values LaSRC would use. It
needs only:

- the scene-centre lat/lon and the time weights, which are already printed in each V0 LaSRC log
- the ancillary pool

Output: a table of granules where S1, S2 or S3 would give LaSRC a different stored water vapor or
ozone than V0, and by how many steps. Granules where all scenarios match V0 will give identical SR,
so they can be reported without running LaSRC. This takes minutes instead of days and shows
up front where differences can occur.

## 7. Phases

| Phase | Work | Deliverable | Check before moving on |
|---|---|---|---|
| **0. Inventory** | Free space and quota on the HPC scratch; granule lists; V0 outputs present for every S30 and L30 granule; Landsat workspace, container and binary; L30 output size per run | Go / no-go note ([7.1](#71-phase-0-results)) | All V0 outputs found; storage covers the chosen retention level |
| **1. Ancillary pool** | One job: download ~1,050 PWAT/TOZNE extracts, relabel forecasts, write a manifest (file, cycle, forecast hour, MD5) | `ancillary_pool/` + manifest | Every needed cycle present; spot-check decoded times |
| **2. Screening** | Compute LaSRC's stored `uwv`/`uoz` for V0, S1, S2, S3 for every granule | Screening table (CSV) | Screened V0 values equal the printed V0 values |
| **3. Pilot** | 10 S30 + 10 L30 granules × S1–S3, mostly from granules flagged in phase 2 | Pilot outputs and logs | Ancillary record correct in every log; trimming works; runtimes as expected |
| **4. Full run** | Chunked Slurm arrays (≤ 30 jobs), final + immediate outputs kept, ancillary record per run | ~3,650 runs (all granules × S1–S3) + ancillary table per product | All jobs exit 0; output counts match; ancillary table complete |
| **5. Validation** | Existing chunked validation CLI (`prepare_lasrc_validation_run.py` → chunks → merge): each scenario vs V0, optionally vs VIIRS; immediate and final | Merged metrics and plots per scenario | Merge complete; V0-vs-VIIRS reproduces the earlier validation |
| **6. Report** | Extend [ANCILLARY_SOURCE_COMPARISON.md](ANCILLARY_SOURCE_COMPARISON.md) with the large-scale results | Doc update on this branch | Team review |

### 7.1 Phase 0 results

| Check | Result |
|---|---|
| Free space | 139 TB free of 141 TB on the HPC scratch file system; no quota reported |
| S30 V0 outputs | 559 product folders, 2.8 TB in the validated S30 GDAS workspace. To do: match against the 605 SAFE IDs / 522 granules (twins, reruns) |
| L30 V0 outputs | `outputs/grib/` in the Landsat `b271310` workspace: 614 date/path-row folders, 2.8 TB. Same workspace holds `outputs/viirs/` (614) |
| L30 run setup | Chunked Slurm arrays (`run_chunk_landsat_apptainer.slurm`, `submit_all_grib_chunks.sh`) with `LASRC_AUX_MODE` per run and `OUTPUT_DIR` per mode; reusable for S1–S3 by pointing each run at its own `LADS/` |
| GRIB mode | The Landsat workspace's `lasrc_aux/LADS/` already holds staged `gdas1.PGrbF00.*` files, so the L30 GDAS runs used the same GRIB path |

### Run details

- **Per run:** the run-specific `LADS/` contains only its two files, named `gdas1.PGrbF00.YYMMDD.HHz` by
  valid time. Static aux files are symlinks. The GDAS staging step is replaced by a stub, because the
  pair is already staged.
- **Logging:** each log starts with an ancillary record (source path, size, MD5, and the GRIB
  reference, forecast and valid times of both files). The live log is written to node-local disk and
  mirrored to shared storage, because editing a log on the share from another machine broke running
  jobs in the smoke test.
- **Validation jobs** are submitted with `--export=NONE`. An inherited conda environment made jobs die
  at start-up in the smoke test.

## 8. Decisions

| # | Topic | Decision |
|---|---|---|
| 1 | Scope | **Run every granule** for S1, S2 and S3 (S30: 605 SAFE IDs, L30: 614). The screening step (section 6) is still produced as a cross-check of the runs, not as a filter. |
| 2 | Retention | **Final products and immediate LaSRC outputs only.** Intermediate HDFs and other debug files are deleted after each run. The scene-centre values LaSRC used (`uwv`, `uoz`, `pres`), the scene-centre lat/lon, the two GRIB files (source, MD5, valid times) and the time weights are written to a per-run record and collected into one table per product. |
| 3 | Validation reference | Each scenario is compared with **V0** (the existing GDAS-analysis run), which answers the latency question directly. The same validation is also run against the **VIIRS**-based product so the results sit next to the earlier GDAS-vs-VIIRS numbers. *(To confirm.)* |
| 4 | GFS resolution | 1°, as in the smoke test. *(Optional 0.25° follow-up, given NWS PNS 26-30.)* |

## 9. Risks

| Risk | Mitigation |
|---|---|
| Storage runs out mid-run | Trimmed retention; per-chunk checks of free space; clean scratch per job |
| Missing V0 output for some granules | Phase 0 inventory; rerun V0 only where missing |
| L30 workspace differs from the validated setup | Pilot includes L30; compare pilot V0 rerun with the stored V0 for 2 granules |
| Ancillary cycle missing on NODD | Manifest check in phase 1; NOMADS or NCAR RDA as fallback |
| Results driven by the 1° shift (KNOWN_ISSUES #2) | Report it alongside; repeat after a fix to #2 |
