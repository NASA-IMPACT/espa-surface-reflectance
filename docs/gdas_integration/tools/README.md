# Ancillary-input selection and checks (low-latency ancillary test)

Copies of the scripts used on Meteor (`/rtmp/tvo/LaSRC_ancillary_latency_20261120/scripts/`).

| Script | What it does |
|---|---|
| `select_pair.py` | The rule: which two GRIB files a granule gets in each scenario (V0, S1, S2, S3), from its acquisition time |
| `audit_ancillary_inputs.py` | Checks every run record against that rule: names, sources, what LaSRC read, checksums, GRIB time, time weights, distinct files per scenario |
| `plot_scenario_inputs.py` | Independent check from the records alone: decodes the model, cycle and forecast hour each run used, compares with the rule, and plots it |
| `scene_centre_values.py` | GRIB2 reader (numpy) and reproduction of LaSRC's water vapor / ozone steps: time weights, bilinear onto the CMG grid, scene-centre cell, truncation |
| `gas_sensitivity.py` | Band sensitivity curves from LaSRC's gas-transmittance coefficients; predicted vs measured SR change |
| `slide_figures.py` | Result figures: maps vs VIIRS / V0, boxplots, storage effect, SR change vs input change (%) |
| `plot_spatial_interp.py` | How 4 GRIB points become the 0.05° CMG cell at the scene centre; optional nearest-vs-bilinear comparison |
| `plot_fields_to_scene.py` | One granule: V0 fields, S1–S3 % differences, and scene-centre value → stored value (`--grib` for 1° cells) |
| `world_land_110m.json` | Coastlines used by the map figures |
| `validation_boxplots.py` | Boxplots per band and summary / failure tables from the validation CLI metrics (one product and stage) |
| `slide_maps_s30.py` | S30 B03 maps vs VIIRS and change vs V0 (granules at MGRS tile centres) |

Try the rule without any data:

```bash
python3 select_pair.py --acq 2025-01-14T14:29:34 --scenario S3 --pool /POOL --gdas-mirror /MIRROR
# gdas1.PGrbF00.250114.12z /POOL/gfs/20250114/06/gfs.t06z.1p00.f006.valid2025011412.grb2
# gdas1.PGrbF00.250114.18z /POOL/gfs/20250114/06/gfs.t06z.1p00.f012.valid2025011418.grb2
```

(It also prints `MISSING ...` and exits with 3 because the example paths do not exist.)
Details and results: [../ANCILLARY_INPUT_VERIFICATION.md](../ANCILLARY_INPUT_VERIFICATION.md), section 6.

Sensitivity figures: [../LASRC_GAS_SENSITIVITY.md](../LASRC_GAS_SENSITIVITY.md), section 8.
