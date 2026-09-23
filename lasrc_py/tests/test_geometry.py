"""Validate geometry functions via Python bindings."""

import numpy as np
import pytest


def test_import():
    """Verify the lasrc package imports correctly."""
    import lasrc
    assert hasattr(lasrc, "__version__")


def test_version():
    import lasrc
    assert lasrc.__version__ == "0.1.0"
