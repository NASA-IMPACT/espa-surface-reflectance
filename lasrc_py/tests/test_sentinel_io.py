"""Tests for Sentinel-2 SAFE metadata parsing."""

from pathlib import Path

import pytest

from lasrc.io.sentinel import _find_tile_metadata, _parse_angles

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
          <ZENITH_ANGLE unit="deg">3.0</ZENITH_ANGLE>
          <AZIMUTH_ANGLE unit="deg">210.0</AZIMUTH_ANGLE>
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


def test_parse_angles(tmp_path):
    safe = _make_safe(tmp_path, ["L1C_A"])
    angles = _parse_angles(_find_tile_metadata(safe))
    assert angles["solar_zenith"] == pytest.approx(66.3360005107259)
    assert angles["solar_azimuth"] == pytest.approx(163.18021256875)
    assert angles["view_zenith"] == pytest.approx(3.5)
    assert angles["view_azimuth"] == pytest.approx(215.0)


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
