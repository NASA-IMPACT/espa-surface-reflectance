# Ancillary-input selection and checks (low-latency ancillary test)

Copies of the scripts used on Meteor (`/rtmp/tvo/LaSRC_ancillary_latency_20261120/scripts/`).

| Script | What it does |
|---|---|
| `select_pair.py` | The rule: which two GRIB files a granule gets in each scenario (V0, S1, S2, S3), from its acquisition time |
| `audit_ancillary_inputs.py` | Checks every run record against that rule: names, sources, what LaSRC read, checksums, GRIB time, time weights, distinct files per scenario |
| `plot_scenario_inputs.py` | Independent check from the records alone: decodes the model, cycle and forecast hour each run used, compares with the rule, and plots it |

Try the rule without any data:

```bash
python3 select_pair.py --acq 2025-01-14T14:29:34 --scenario S3 --pool /POOL --gdas-mirror /MIRROR
# gdas1.PGrbF00.250114.12z /POOL/gfs/20250114/06/gfs.t06z.1p00.f006.valid2025011412.grb2
# gdas1.PGrbF00.250114.18z /POOL/gfs/20250114/06/gfs.t06z.1p00.f012.valid2025011418.grb2
```

(It also prints `MISSING ...` and exits with 3 because the example paths do not exist.)
Details and results: [../ANCILLARY_INPUT_VERIFICATION.md](../ANCILLARY_INPUT_VERIFICATION.md), section 6.
