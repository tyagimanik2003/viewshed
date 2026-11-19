"""Upstream GDAL behaviour the engine works around. If one of these starts failing after a
GDAL upgrade, the quirk was fixed upstream: update core.py and docs/gdal_notes.md."""
from __future__ import annotations

import math

import numpy as np
from osgeo import gdal

from conftest import CURV, R_EARTH


def _farthest_visible(dem, tmp_path, **kw):
    out = tmp_path / "vs.tif"
    gdal.Run("raster", "viewshed", input=str(dem), output=str(out), max_distance=8000,
             curvature_coefficient=CURV, overwrite=True, **kw).Finalize()
    ds = gdal.Open(str(out))
    a, gt = ds.GetRasterBand(1).ReadAsArray(), ds.GetGeoTransform()
    c, r = np.meshgrid(np.arange(a.shape[1]) + 0.5, np.arange(a.shape[0]) + 0.5)
    return float(np.hypot(gt[0] + c * gt[1], gt[3] + r * gt[5])[a == 255].max())


def horizon(h):
    return math.sqrt(2 * R_EARTH * h / CURV)


def test_height_inside_position_is_honoured(flat_aeqd, tmp_path):
    assert abs(_farthest_visible(flat_aeqd, tmp_path, position=[0.0, 0.0, 0.5]) - horizon(0.5)) < 60


def test_quirk_separate_height_option_is_ignored(flat_aeqd, tmp_path):
    """GDAL 3.13.3: -z/--height is ignored; the observer stays at the 2 m default."""
    d = _farthest_visible(flat_aeqd, tmp_path, position=[0.0, 0.0], height=0.5)
    assert abs(d - horizon(2.0)) < 60, f"farthest visible {d:.0f} m; fixed upstream? horizon(0.5)={horizon(0.5):.0f}"
