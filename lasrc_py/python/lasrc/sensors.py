"""Sensor configuration dictionaries for I/O and metadata."""

SENSORS = {
    "LANDSAT_8": {
        "name": "LANDSAT_8",
        "satellite": "LANDSAT_8",
        "instrument": "OLI_TIRS",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band9",
        ],
        "thm_band_names": ["bt_band10", "bt_band11"],
        "input_band_names": [
            "band1", "band2", "band3", "band4",
            "band5", "band6", "band7", "band9",
        ],
        "input_thm_names": ["band10", "band11"],
        "resolution_m": 30.0,
        "nsr_bands": 8,
    },
    "LANDSAT_9": {
        "name": "LANDSAT_9",
        "satellite": "LANDSAT_9",
        "instrument": "OLI_TIRS",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band9",
        ],
        "thm_band_names": ["bt_band10", "bt_band11"],
        "input_band_names": [
            "band1", "band2", "band3", "band4",
            "band5", "band6", "band7", "band9",
        ],
        "input_thm_names": ["band10", "band11"],
        "resolution_m": 30.0,
        "nsr_bands": 8,
    },
    "SENTINEL_2A": {
        "name": "SENTINEL_2A",
        "satellite": "SENTINEL-2A",
        "instrument": "MSI",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band8",
            "sr_band8a", "sr_band9", "sr_band10", "sr_band11", "sr_band12",
        ],
        "thm_band_names": [],
        "input_band_names": [
            "B01", "B02", "B03", "B04", "B05", "B06",
            "B07", "B08", "B8A", "B09", "B10", "B11", "B12",
        ],
        "resolution_m": 10.0,
        "nsr_bands": 13,
    },
    "SENTINEL_2B": {
        "name": "SENTINEL_2B",
        "satellite": "SENTINEL-2B",
        "instrument": "MSI",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band8",
            "sr_band8a", "sr_band9", "sr_band10", "sr_band11", "sr_band12",
        ],
        "thm_band_names": [],
        "input_band_names": [
            "B01", "B02", "B03", "B04", "B05", "B06",
            "B07", "B08", "B8A", "B09", "B10", "B11", "B12",
        ],
        "resolution_m": 10.0,
        "nsr_bands": 13,
    },
    "SENTINEL_2C": {
        "name": "SENTINEL_2C",
        "satellite": "SENTINEL-2C",
        "instrument": "MSI",
        "refl_band_names": [
            "sr_band1", "sr_band2", "sr_band3", "sr_band4",
            "sr_band5", "sr_band6", "sr_band7", "sr_band8",
            "sr_band8a", "sr_band9", "sr_band10", "sr_band11", "sr_band12",
        ],
        "thm_band_names": [],
        "input_band_names": [
            "B01", "B02", "B03", "B04", "B05", "B06",
            "B07", "B08", "B8A", "B09", "B10", "B11", "B12",
        ],
        "resolution_m": 10.0,
        "nsr_bands": 13,
    },
}
