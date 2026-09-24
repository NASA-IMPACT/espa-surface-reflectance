"""Tests for Sentinel-2 SAFE file discovery and metadata parsing."""

from pathlib import Path

import numpy as np
import pytest

from lasrc.io.sentinel import (
    _dn_to_toa,
    _find_band_file,
    _find_tile_metadata,
    _parse_angles,
    _parse_quantification,
)

MTD_TL_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<n1:Level-1C_Tile_ID xmlns:n1="https://psd-15.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd">
  <n1:Geometric_Info>
    <Tile_Angles>
      {sun}
      <Mean_Viewing_Incidence_Angle_List>
        {view}
      </Mean_Viewing_Incidence_Angle_List>
    </Tile_Angles>
  </n1:Geometric_Info>
</n1:Level-1C_Tile_ID>
"""

SUN = """<Mean_Sun_Angle>
        <ZENITH_ANGLE unit="deg">66.3360005107259</ZENITH_ANGLE>
        <AZIMUTH_ANGLE unit="deg">163.18021256875</AZIMUTH_ANGLE>
      </Mean_Sun_Angle>"""

VIEW = """<Mean_Viewing_Incidence_Angle bandId="0">
          <ZENITH_ANGLE unit="deg">3.37889124840365</ZENITH_ANGLE>
          <AZIMUTH_ANGLE unit="deg">213.175195014188</AZIMUTH_ANGLE>
        </Mean_Viewing_Incidence_Angle>
        <Mean_Viewing_Incidence_Angle bandId="4">
          <ZENITH_ANGLE unit="deg">4.0</ZENITH_ANGLE>
          <AZIMUTH_ANGLE unit="deg">220.0</AZIMUTH_ANGLE>
        </Mean_Viewing_Incidence_Angle>"""


def _make_safe(tmp_path: Path, granules: list[str], sun: str = SUN, view: str = VIEW) -> Path:
    safe = tmp_path / "S2C_MSIL1C_TEST.SAFE"
    safe.mkdir()
    (safe / "MTD_MSIL1C.xml").write_text("<root/>")
    for g in granules:
        gdir = safe / "GRANULE" / g
        gdir.mkdir(parents=True)
        (gdir / "MTD_TL.xml").write_text(MTD_TL_TEMPLATE.format(sun=sun, view=view))
    return safe


def test_find_tile_metadata_in_granule_dir(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_T14TLS_A002070_20250127T174637"])
    assert _find_tile_metadata(safe) == (
        safe / "GRANULE" / "L1C_T14TLS_A002070_20250127T174637" / "MTD_TL.xml"
    )


def test_find_tile_metadata_missing_raises(tmp_path):
    safe = _make_safe(tmp_path, [])
    with pytest.raises(FileNotFoundError, match="MTD_TL.xml"):
        _find_tile_metadata(safe)


def test_find_tile_metadata_multiple_raises(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A", "L1C_B"])
    with pytest.raises(ValueError, match="MTD_TL.xml"):
        _find_tile_metadata(safe)


def test_parse_angles_matches_espa_metadata_chain(tmp_path):
    # C LaSRC reads angles from the ESPA XML: each value is parsed into a
    # float, written with "%f", and read back into a float. The view angle is
    # the first Mean_Viewing_Incidence_Angle only (bandId 0 here), not a mean.
    safe = _make_safe(tmp_path, ["L1C_A"])
    angles = _parse_angles(_find_tile_metadata(safe))
    assert angles == {
        "solar_zenith": 66.33599853515625,
        "solar_azimuth": 163.18020629882812,
        "view_zenith": 3.3788909912109375,
        "view_azimuth": 213.17520141601562,
    }


@pytest.mark.parametrize(
    "sun,view,missing",
    [
        ("", VIEW, "Mean_Sun_Angle"),
        (SUN.replace("<ZENITH_ANGLE unit=\"deg\">66.3360005107259</ZENITH_ANGLE>", ""),
         VIEW, "Mean_Sun_Angle"),
        (SUN, "", "Mean_Viewing_Incidence_Angle"),
    ],
)
def test_parse_angles_missing_element_raises(tmp_path, sun, view, missing):
    safe = _make_safe(tmp_path, ["L1C_A"], sun=sun, view=view)
    with pytest.raises(ValueError, match=missing):
        _parse_angles(_find_tile_metadata(safe))


def test_parse_angles_rejects_product_metadata(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A"])
    with pytest.raises(ValueError, match="Mean_Sun_Angle"):
        _parse_angles(safe / "MTD_MSIL1C.xml")


def _add_jp2(safe: Path, relpath: str) -> Path:
    path = safe / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def test_find_band_file_ignores_qi_masks(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A"])
    _add_jp2(safe, "GRANULE/L1C_A/QI_DATA/MSK_QUALIT_B01.jp2")
    _add_jp2(safe, "GRANULE/L1C_A/QI_DATA/MSK_DETFOO_B01.jp2")
    band = _add_jp2(safe, "GRANULE/L1C_A/IMG_DATA/T14TLS_20250127T174641_B01.jp2")
    assert _find_band_file(safe, "B01") == band


def test_find_band_file_does_not_confuse_b8a_with_b08(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A"])
    _add_jp2(safe, "GRANULE/L1C_A/IMG_DATA/T14TLS_20250127T174641_B8A.jp2")
    b08 = _add_jp2(safe, "GRANULE/L1C_A/IMG_DATA/T14TLS_20250127T174641_B08.jp2")
    assert _find_band_file(safe, "B08") == b08


def test_find_band_file_missing_raises(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A"])
    _add_jp2(safe, "GRANULE/L1C_A/QI_DATA/MSK_QUALIT_B01.jp2")
    with pytest.raises(FileNotFoundError, match="B01"):
        _find_band_file(safe, "B01")


def test_find_band_file_multiple_raises(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A", "L1C_B"])
    _add_jp2(safe, "GRANULE/L1C_A/IMG_DATA/T14TLS_20250127T174641_B01.jp2")
    _add_jp2(safe, "GRANULE/L1C_B/IMG_DATA/T14TLS_20250127T174641_B01.jp2")
    with pytest.raises(ValueError, match="B01"):
        _find_band_file(safe, "B01")


MTD_MSIL1C_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<n1:Level-1C_User_Product xmlns:n1="https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd">
  <n1:General_Info>
    <Product_Image_Characteristics>
      {quant}
      <Radiometric_Offset_List>
        {offsets}
      </Radiometric_Offset_List>
    </Product_Image_Characteristics>
  </n1:General_Info>
</n1:Level-1C_User_Product>
"""

QUANT = '<QUANTIFICATION_VALUE unit="none">10000</QUANTIFICATION_VALUE>'


def _offsets(band_ids) -> str:
    return "\n".join(
        f'<RADIO_ADD_OFFSET band_id="{i}">{-1000 - i}</RADIO_ADD_OFFSET>' for i in band_ids
    )


def _write_mtd(tmp_path: Path, quant: str = QUANT, offsets: str = _offsets(range(13))) -> Path:
    path = tmp_path / "MTD_MSIL1C.xml"
    path.write_text(MTD_MSIL1C_TEMPLATE.format(quant=quant, offsets=offsets))
    return path


def test_parse_quantification_per_band_offsets(tmp_path):
    quant, offsets = _parse_quantification(_write_mtd(tmp_path))
    assert quant == 10000.0
    assert offsets == [-1000.0 - i for i in range(13)]


def test_parse_quantification_missing_offsets_raises(tmp_path):
    with pytest.raises(ValueError, match="RADIO_ADD_OFFSET"):
        _parse_quantification(_write_mtd(tmp_path, offsets=""))


def test_parse_quantification_missing_value_raises(tmp_path):
    with pytest.raises(ValueError, match="QUANTIFICATION_VALUE"):
        _parse_quantification(_write_mtd(tmp_path, quant=""))


@pytest.mark.parametrize(
    "band_ids",
    [range(12), [*range(13), 0], [*range(12), 13]],
    ids=["incomplete", "duplicate", "out-of-range"],
)
def test_parse_quantification_bad_offsets_raise(tmp_path, band_ids):
    with pytest.raises(ValueError, match="RADIO_ADD_OFFSET"):
        _parse_quantification(_write_mtd(tmp_path, offsets=_offsets(band_ids)))


def test_dn_to_toa_matches_c_float_arithmetic():
    # C LaSRC: toa = (DN + add_offset) * scale_factor in float, where
    # scale_factor is 1/QUANTIFICATION_VALUE written to the ESPA XML as
    # "%10.8f". Dividing by 10000 instead differs by 1 ULP for these DNs.
    toa = _dn_to_toa(np.array([4.0, 8.0], dtype=np.float32), -1000.0, 10000.0)
    assert toa.dtype == np.float32
    assert toa.tolist() == [-0.09959999471902847, -0.09919999539852142]
