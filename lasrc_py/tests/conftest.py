"""Shared test fixtures for integration tests."""

import numpy as np
import pytest


@pytest.fixture
def small_scene():
    """Generate a small synthetic test scene (10x10 pixels)."""
    nlines, nsamps = 10, 10
    rng = np.random.default_rng(42)

    toa_bands = [rng.uniform(0.05, 0.3, (nlines, nsamps)).astype(np.float32) for _ in range(8)]
    bt_bands = [rng.uniform(270.0, 310.0, (nlines, nsamps)).astype(np.float32) for _ in range(2)]
    qa_band = np.ones((nlines, nsamps), dtype=np.uint16)
    sza = np.full((nlines, nsamps), 30.0, dtype=np.float32)
    saa = np.full((nlines, nsamps), 150.0, dtype=np.float32)
    vza = np.full((nlines, nsamps), 0.0, dtype=np.float32)
    vaa = np.full((nlines, nsamps), 0.0, dtype=np.float32)

    return {
        "toa_bands": toa_bands,
        "bt_bands": bt_bands,
        "qa_band": qa_band,
        "sza": sza,
        "saa": saa,
        "vza": vza,
        "vaa": vaa,
        "nlines": nlines,
        "nsamps": nsamps,
    }
