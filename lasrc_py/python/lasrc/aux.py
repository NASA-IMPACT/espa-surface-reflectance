"""Load LaSRC look-up tables and ancillary data into NumPy arrays.

This module is the single source of truth for reading the LaSRC auxiliary
inputs. The on-disk formats are:

  - ANGLE_NEW.hdf            : HDF4, datasets TSMAX/TSMIN/NBFIC/NBFI/TTV/TTS/INDTS
  - RES_LUT_*.hdf            : HDF4, per-band datasets NRLUT_BAND_<n> [NSOLAR, NAOT, NPRES]
  - TRANS_LUT_*.ASCII        : plain ASCII transmission table
  - AERO_LUT_*.ASCII         : plain ASCII spherical-albedo / extinction table
  - V*04ANC.A*.h5            : HDF5 daily VIIRS CMG water vapor + ozone
  - CMGDEM.hdf               : HDF4 averaged elevation
  - ratiomapndwiexp.hdf      : HDF4 band-ratio / NDWI climatology
"""

import numpy as np

# 6S LUT grid dimensions (shared by all sensors).
NSOLAR = 8000     # solar/view geometry samples per (band, pressure, aot)
NAOT = 22         # aerosol optical thickness steps
NPRES = 7         # pressure levels
NSUNANGLE = 22    # sun-angle slots allocated for the transmission table
NSUNANGLE_FILE = 21  # sun-angle lines actually present per pressure block


# -- 6S look-up tables -----------------------------------------------------------

def load_angle_lut(path: str) -> dict:
    """Load the angle LUT (ANGLE_NEW.hdf) from HDF4.

    Returns a dict with keys tsmax/tsmin/nbfic/nbfi/ttv/tts/indts, matching the
    corresponding PyLookupTables constructor arguments.
    """
    from pyhdf.SD import SD, SDC

    hdf = SD(str(path), SDC.READ)
    result = {}
    for name in ["TSMAX", "TSMIN", "NBFIC", "TTV", "TTS"]:
        result[name.lower()] = (
            hdf.select(name).get().astype(np.float64).ravel().tolist()
        )
    result["nbfi"] = hdf.select("NBFI").get().astype(np.int32).ravel().tolist()
    result["indts"] = hdf.select("INDTS").get().astype(np.int32).ravel().tolist()
    hdf.end()
    return result


def load_rolutt(path: str, band_names: list[str]) -> list:
    """Load the intrinsic reflectance LUT (RES_LUT_*.hdf) from HDF4.

    HDF4 stores one dataset per band (named in ``band_names``), each shaped
    [NSOLAR, NAOT, NPRES]. The Rust core expects a flat layout of
    band x NPRES x NAOT x NSOLAR.
    """
    from pyhdf.SD import SD, SDC

    nsr_bands = len(band_names)
    hdf = SD(str(path), SDC.READ)
    rolutt = np.zeros(nsr_bands * NPRES * NAOT * NSOLAR, dtype=np.float64)

    for ib, ds_name in enumerate(band_names):
        data = hdf.select(ds_name).get().astype(np.float64)  # [NSOLAR, NAOT, NPRES]
        for ip in range(NPRES):
            for ia in range(NAOT):
                base = (ib * NPRES * NAOT * NSOLAR
                        + ip * NAOT * NSOLAR
                        + ia * NSOLAR)
                rolutt[base:base + NSOLAR] = data[:, ia, ip]

    hdf.end()
    return rolutt.tolist()


def load_aero_ascii(path: str, nsr_bands: int) -> tuple[list, list]:
    """Load spherical albedo and normalized extinction from the AERO ASCII file.

    Format: per band, NPRES pressure blocks, NAOT lines per block
    (aot sphalbt normext). Returns (sphalbt, normext) flat lists in
    band x NPRES x NAOT layout.
    """
    sphalbt = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)
    normext = np.zeros(nsr_bands * NPRES * NAOT, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        idx += 1  # band header line
        for ip in range(NPRES):
            idx += 1  # pressure header line
            for ia in range(NAOT):
                parts = lines[idx].split()
                base = ib * NPRES * NAOT + ip * NAOT + ia
                sphalbt[base] = float(parts[1])
                normext[base] = float(parts[2])
                idx += 1

    return sphalbt.tolist(), normext.tolist()


def load_trans_ascii(path: str, nsr_bands: int) -> list:
    """Load the transmission LUT from the TRANS ASCII file.

    Format: per band, NPRES pressure blocks, NSUNANGLE_FILE sun-angle lines per
    block, each line a sun angle followed by NAOT transmission values. The Rust
    core allocates NSUNANGLE slots per (band, pressure, aot) but only the first
    NSUNANGLE_FILE are populated from the file. Returns a flat list in
    band x NPRES x NAOT x NSUNANGLE layout.
    """
    transt = np.zeros(nsr_bands * NPRES * NAOT * NSUNANGLE, dtype=np.float64)

    with open(path) as f:
        lines = f.readlines()

    idx = 0
    for ib in range(nsr_bands):
        idx += 1  # band header line
        for ip in range(NPRES):
            idx += 1  # pressure header line
            for isun in range(NSUNANGLE_FILE):
                parts = lines[idx].split()
                # parts[0] is the sun angle; parts[1:] are per-AOT transmission.
                for ia in range(NAOT):
                    base = (ib * NPRES * NAOT * NSUNANGLE
                            + ip * NAOT * NSUNANGLE
                            + ia * NSUNANGLE
                            + isun)
                    transt[base] = float(parts[1 + ia])
                idx += 1

    return transt.tolist()


def load_lut(angle_path: str, intref_path: str, transm_path: str,
             sphera_path: str, band_names: list[str]) -> dict:
    """Load all 6S look-up tables and return PyLookupTables constructor kwargs.

    Parameters
    ----------
    angle_path : path to ANGLE_NEW.hdf (HDF4)
    intref_path : path to RES_LUT_*.hdf (HDF4)
    transm_path : path to TRANS_LUT_*.ASCII
    sphera_path : path to AERO_LUT_*.ASCII
    band_names : RES_LUT dataset names for this sensor (see SENSORS[...]['lut_band_names'])
    """
    nsr_bands = len(band_names)
    angle_data = load_angle_lut(angle_path)
    rolutt = load_rolutt(intref_path, band_names)
    sphalbt, normext = load_aero_ascii(sphera_path, nsr_bands)
    transt = load_trans_ascii(transm_path, nsr_bands)

    return {
        "rolutt": rolutt,
        "transt": transt,
        "sphalbt": sphalbt,
        "normext": normext,
        **angle_data,
        "nsr_bands": nsr_bands,
    }


# -- Atmospheric ancillary data --------------------------------------------------

def load_auxiliary_data(aux_path: str, aux_source: str = "VIIRS") -> dict:
    """Load water vapor and ozone grids.

    Args:
        aux_path: Path to the daily auxiliary file.
        aux_source: "VIIRS" or "MODIS".

    Returns dict with keys matching the PyAuxiliaryData constructor args
    (wv, oz, wv_scale, oz_scale, wv_default, oz_default).
    """
    if aux_source == "VIIRS":
        return _load_viirs_aux(aux_path)
    return _load_modis_aux(aux_path)


def _load_viirs_aux(aux_path: str) -> dict:
    """Load water vapor + ozone from a daily VIIRS CMG HDF5 file.

    Scale factors (VIIRS): WV / 200 -> g/cm^2, OZ / 400 -> cm-atm.
    """
    import h5py

    base = "/HDFEOS/GRIDS/VIIRS_CMG/Data Fields/"
    with h5py.File(str(aux_path), "r") as f:
        wv = np.array(f[base + "Coarse Resolution Water Vapor"], dtype=np.int16).ravel()
        oz = np.array(f[base + "Coarse Resolution Ozone"], dtype=np.int16).ravel()

    return {
        "wv": wv.tolist(),
        "oz": oz.tolist(),
        "wv_scale": 200.0,
        "oz_scale": 400.0,
        "wv_default": 2.5,
        "oz_default": 0.3,
    }


def _load_modis_aux(aux_path: str) -> dict:
    """Placeholder for MODIS auxiliary data loading."""
    raise NotImplementedError("MODIS auxiliary loading not yet implemented")


def load_dem(dem_path: str) -> list:
    """Load the CMG DEM ('averaged elevation') from HDF4."""
    from pyhdf.SD import SD, SDC

    hdf = SD(str(dem_path), SDC.READ)
    dem = hdf.select("averaged elevation").get().astype(np.int16).ravel()
    hdf.end()
    return dem.tolist()


def load_ratio_file(ratio_path: str) -> dict:
    """Load NDWI and band-ratio climatology from HDF4.

    HDF4 dataset names differ from the internal field names. The mapping
    (from the C reference code) maps Landsat band numbers to the internal
    ratio slots used by the aerosol retrieval.
    """
    from pyhdf.SD import SD, SDC

    name_map = {
        "ratiob1": "average ratio b9",
        "ratiob2": "average ratio b10",
        "ratiob7": "average ratio b7",
        "intratiob1": "inter ratiob9",
        "intratiob2": "inter ratiob10",
        "intratiob7": "inter ratiob7",
        "slpratiob1": "slope ratiob9",
        "slpratiob2": "slope ratiob10",
        "slpratiob7": "slope ratiob7",
        "andwi": "average ndvi",
        "sndwi": "standard ndvi",
    }

    hdf = SD(str(ratio_path), SDC.READ)
    result = {}
    for field_name, dataset_name in name_map.items():
        result[field_name] = np.array(
            hdf.select(dataset_name).get(), dtype=np.int16
        ).ravel().tolist()
    hdf.end()
    return result
