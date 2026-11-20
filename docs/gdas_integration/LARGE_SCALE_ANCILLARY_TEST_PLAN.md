# Plan: large-scale test of low-latency ancillary sources in LaSRC

Follow-up to the 3-granule smoke test in
[ANCILLARY_SOURCE_COMPARISON.md](ANCILLARY_SOURCE_COMPARISON.md). This test reruns the HLS
validation granule set (S30 and L30) with each low-latency ancillary scenario and compares the
results with the validated GDAS-analysis run.

**Status:** plan, not started. Open decisions are listed in [section 8](#8-open-decisions).

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
| S30 | 605 SAFE IDs (522 granules after twins) | Jan and Jul 2025 | Sentinel-2 L1C SAFE archive used for the validation | Validated GDAS workspace (to be confirmed in phase 0) |
| L30 | 614 | Jan and Jul 2025 | Landsat L1 inputs used for the validation | Landsat `b271310` GDAS workspace (location to be confirmed in phase 0) |

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

Measured in the smoke test: one S30 run in debug mode is 6.1 GB (final COGs 0.25 GB, immediate
LaSRC `.img` bands ~3 GB, intermediate HDFs ~2.8 GB). L30 sizes to be measured in phase 0.

| Output kept per run | S30 (~1,815 new runs) | L30 (~1,842 new runs) | Total |
|---|---|---|---|
| Full debug (current) | ~11 TB | ~3–5 TB | ~15 TB |
| **Trimmed** (recommended): final COGs + immediate core bands B01–B04, B8A, B11, B12 + logs | ~3.6 TB | ~2–3 TB | **~6 TB** |
| Final COGs only | ~0.45 TB | ~0.3 TB | ~0.8 TB |

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
| **0. Inventory** | Free space and quota on the HPC scratch; granule lists; V0 outputs present for every S30 and L30 granule; Landsat workspace, container and binary; L30 output size per run | Go / no-go note | All V0 outputs found; storage covers the chosen retention level |
| **1. Ancillary pool** | One job: download ~1,050 PWAT/TOZNE extracts, relabel forecasts, write a manifest (file, cycle, forecast hour, MD5) | `ancillary_pool/` + manifest | Every needed cycle present; spot-check decoded times |
| **2. Screening** | Compute LaSRC's stored `uwv`/`uoz` for V0, S1, S2, S3 for every granule | Screening table (CSV) | Screened V0 values equal the printed V0 values |
| **3. Pilot** | 10 S30 + 10 L30 granules × S1–S3, mostly from granules flagged in phase 2 | Pilot outputs and logs | Ancillary record correct in every log; trimming works; runtimes as expected |
| **4. Full run** | Chunked Slurm arrays (≤ 30 jobs), trimmed outputs, ancillary record in every log | ~3,650 runs (or only flagged granules) | All jobs exit 0; output counts match |
| **5. Validation** | Existing chunked validation CLI (`prepare_lasrc_validation_run.py` → chunks → merge): each scenario vs V0, optionally vs VIIRS; immediate and final | Merged metrics and plots per scenario | Merge complete; V0-vs-VIIRS reproduces the earlier validation |
| **6. Report** | Extend [ANCILLARY_SOURCE_COMPARISON.md](ANCILLARY_SOURCE_COMPARISON.md) with the large-scale results | Doc update on this branch | Team review |

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

## 8. Open decisions

1. **Scope:** run every granule, or screen first and run only granules where a scenario changes
   LaSRC's stored values? (Recommended: screen first.)
2. **Retention:** full debug, trimmed (recommended) or final COGs only?
3. **Validation reference:** V0 only, or V0 and VIIRS?
4. **GFS resolution:** 1° for consistency with GDAS (as in the smoke test), or add 0.25° for one
   scenario, given the proposed removal of GFS 1° (NWS PNS 26-30)?

## 9. Risks

| Risk | Mitigation |
|---|---|
| Storage runs out mid-run | Trimmed retention; per-chunk checks of free space; clean scratch per job |
| Missing V0 output for some granules | Phase 0 inventory; rerun V0 only where missing |
| L30 workspace differs from the validated setup | Pilot includes L30; compare pilot V0 rerun with the stored V0 for 2 granules |
| Ancillary cycle missing on NODD | Manifest check in phase 1; NOMADS or NCAR RDA as fallback |
| Results driven by the 1° shift (KNOWN_ISSUES #2) | Report it alongside; repeat after a fix to #2 |
