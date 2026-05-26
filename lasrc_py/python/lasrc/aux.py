"""Load auxiliary data from HDF4/5 files into NumPy arrays."""

import numpy as np


def load_lut_from_hdf(angle_hdf_path: str, intref_path: str,
                      transm_path: str, sphera_path: str,
                      nsr_bands: int) -> dict:
    """Load 6S lookup tables from HDF files.

    Returns dict with keys matching PyLookupTables constructor args.
    """
    import h5py

    with h5py.File(angle_hdf_path, "r") as f:
        tsmax = np.array(f["tsmax"], dtype=np.float64).ravel()
        tsmin = np.array(f["tsmin"], dtype=np.float64).ravel()
        nbfic = np.array(f["nbfic"], dtype=np.float64).ravel()
        nbfi = np.array(f["nbfi"], dtype=np.int32).ravel()
        ttv = np.array(f["ttv"], dtype=np.float64).ravel()
        tts = np.array(f["tts"], dtype=np.float64).ravel()

    with h5py.File(intref_path, "r") as f:
        rolutt = np.array(f["rolutt"], dtype=np.float64).ravel()

    with h5py.File(transm_path, "r") as f:
        transt = np.array(f["transt"], dtype=np.float64).ravel()

    with h5py.File(sphera_path, "r") as f:
        sphalbt = np.array(f["sphalbt"], dtype=np.float64).ravel()
        normext = np.array(f["normext"], dtype=np.float64).ravel()

    return {
        "rolutt": rolutt.tolist(),
        "transt": transt.tolist(),
        "sphalbt": sphalbt.tolist(),
        "normext": normext.tolist(),
        "tsmax": tsmax.tolist(),
        "tsmin": tsmin.tolist(),
        "nbfic": nbfic.tolist(),
        "nbfi": nbfi.tolist(),
        "ttv": ttv.tolist(),
        "tts": tts.tolist(),
        "nsr_bands": nsr_bands,
    }


def load_auxiliary_data(aux_path: str, aux_source: str = "VIIRS") -> dict:
    """Load DEM, water vapor, ozone, and band ratio auxiliary data.

    Args:
        aux_path: Path to auxiliary HDF file.
        aux_source: "VIIRS" or "MODIS".

    Returns dict with keys matching PyAuxiliaryData constructor args.
    """
    if aux_source == "VIIRS":
        return _load_viirs_aux(aux_path)
    else:
        return _load_modis_aux(aux_path)


def _load_viirs_aux(aux_path: str) -> dict:
    import h5py

    wv_scale = 200.0
    oz_scale = 400.0
    wv_default = 2.5
    oz_default = 0.3

    with h5py.File(aux_path, "r") as f:
        base = "/HDFEOS/GRIDS/VIIRS_CMG/Data Fields/"
        wv = np.array(f[base + "Average_Column_WV_CMG"], dtype=np.int16).ravel()
        oz = np.array(f[base + "Average_Column_Ozone_CMG"], dtype=np.int16).ravel()

    return {
        "wv": wv.tolist(),
        "oz": oz.tolist(),
        "wv_scale": wv_scale,
        "oz_scale": oz_scale,
        "wv_default": wv_default,
        "oz_default": oz_default,
    }


def _load_modis_aux(aux_path: str) -> dict:
    """Placeholder for MODIS auxiliary data loading."""
    raise NotImplementedError("MODIS auxiliary loading not yet implemented")


def load_dem(dem_path: str) -> list:
    """Load CMG DEM from HDF4 file."""
    from pyhdf.SD import SD, SDC

    hdf = SD(dem_path, SDC.READ)
    dem = hdf.select("CMG_DEM").get().astype(np.int16).ravel()
    hdf.end()
    return dem.tolist()


def load_ratio_file(ratio_path: str) -> dict:
    """Load NDWI and band ratio data from HDF4 file."""
    from pyhdf.SD import SD, SDC

    hdf = SD(ratio_path, SDC.READ)
    result = {}
    for name in ["ratiob1", "ratiob2", "ratiob7",
                  "intratiob1", "intratiob2", "intratiob7",
                  "slpratiob1", "slpratiob2", "slpratiob7",
                  "andwi", "sndwi"]:
        result[name] = np.array(hdf.select(name).get(), dtype=np.int16).ravel().tolist()
    hdf.end()
    return result
