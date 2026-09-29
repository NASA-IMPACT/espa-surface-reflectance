"""Tests for input validation at the Python/Rust binding boundary."""

import numpy as np
import pytest

import lasrc

NSOLAR_ZEN_VALS = 22


def _placeholder_lut() -> "lasrc.PyLookupTables":
    f64 = np.zeros(1, dtype=np.float64)
    i32 = np.zeros(1, dtype=np.int32)
    return lasrc.PyLookupTables(
        rolutt=f64, transt=f64, sphalbt=f64, normext=f64, tsmax=f64, tsmin=f64,
        nbfic=f64, nbfi=i32, ttv=f64, tts=np.zeros(NSOLAR_ZEN_VALS), indts=i32,
        nsr_bands=13,
    )


def _placeholder_aux() -> "lasrc.PyAuxiliaryData":
    i16 = np.zeros(1, dtype=np.int16)
    return lasrc.PyAuxiliaryData(
        dem=i16, wv=i16, oz=i16, ratiob1=i16, ratiob2=i16, ratiob7=i16,
        intratiob1=i16, intratiob2=i16, intratiob7=i16, slpratiob1=i16,
        slpratiob2=i16, slpratiob7=i16, andwi=i16, sndwi=i16,
        wv_scale=1.0, oz_scale=1.0, wv_default=0.0, oz_default=0.0,
    )


def test_sentinel_rejects_non_contiguous_qa_band():
    toa_bands = [np.zeros((4, 4), dtype=np.float32) for _ in range(13)]
    qa_band = np.zeros((4, 8), dtype=np.uint16)[:, ::2]
    assert not qa_band.flags.c_contiguous

    with pytest.raises(ValueError, match="qa_band"):
        lasrc.process_sentinel_surface_reflectance(
            sensor_name="SENTINEL_2A", toa_bands=toa_bands,
            solar_zenith=30.0, solar_azimuth=150.0, view_zenith=0.0, view_azimuth=0.0,
            qa_band=qa_band, lut=_placeholder_lut(), aux=_placeholder_aux(),
            ul_corner_x=0.0, ul_corner_y=0.0, pixel_size_x=10.0, pixel_size_y=10.0,
            utm_zone=14,
        )
