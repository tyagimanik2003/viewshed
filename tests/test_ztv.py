"""ZTV (multi-observer visibility counts) and web-map export, checked against known answers.

All inputs are synthetic and local, so these run offline.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from osgeo import gdal, osr
from pyproj import Geod, Transformer

from conftest import CURV, ELEV, OBS_LAT, OBS_LON, R_EARTH
from viewshed import INVISIBLE, VISIBLE, viewshed, ztv
from viewshed.export import viewshed_lut, viewshed_overlay, ztv_classes, ztv_lut, ztv_overlay

RADIUS = 5000.0
GEOD = Geod(ellps="WGS84")


@pytest.fixture
def plateau(tmp_path):
    """Flat EPSG:4326 plateau, 1 arc-second cells, +-11 km around the test location."""
    res = 1 / 3600
    half_lon = 11000 / (111320 * math.cos(math.radians(OBS_LAT)))
    half_lat = 11000 / 110950
    ncol, nrow = int(2 * half_lon / res), int(2 * half_lat / res)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    p = tmp_path / "plateau.tif"
    ds = gdal.GetDriverByName("GTiff").Create(str(p), ncol, nrow, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((OBS_LON - ncol / 2 * res, res, 0, OBS_LAT + nrow / 2 * res, 0, -res))
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(np.full((nrow, ncol), ELEV, np.float32))
    ds = None
    return p


def cell_lonlat(res):
    gt = res.geotransform
    h, w = res.count.shape
    c, r = np.meshgrid(np.arange(w) + 0.5, np.arange(h) + 0.5)
    t = Transformer.from_crs(res.crs_wkt, "EPSG:4326", always_xy=True)
    return t.transform(gt[0] + c * gt[1], gt[3] + r * gt[5])


def east(lon, lat, metres):
    lo, la, _ = GEOD.fwd(lon, lat, 90, metres)
    return lo, la


def test_single_observer_ztv_equals_viewshed(plateau):
    """One observer: the ZTV grid is the viewshed grid, and count == 1 exactly where visible."""
    z = ztv([(OBS_LON, OBS_LAT)], RADIUS, observer_height=1.5, target_height=0, dem=plateau)
    v = viewshed(OBS_LON, OBS_LAT, RADIUS, observer_height=1.5, dem=plateau)
    assert z.geotransform == pytest.approx(v.geotransform)
    assert np.array_equal(z.count == 1, v.codes == VISIBLE)
    assert np.array_equal(z.analysed, np.isin(v.codes, (VISIBLE, INVISIBLE)))


def test_counts_follow_the_curvature_horizon(plateau):
    """Two observers 3 km apart on a flat plateau: a cell sees observer i iff it lies within that
    observer's curvature horizon (and the radius), so the count is analytic."""
    h = 1.5
    horizon = math.sqrt(2 * R_EARTH * h / CURV)  # ~4.72 km, inside the 5 km radius
    a = east(OBS_LON, OBS_LAT, -1500)
    b = east(OBS_LON, OBS_LAT, 1500)
    z = ztv([(*a, h), (*b, h)], RADIUS, target_height=0, dem=plateau)
    lon, lat = cell_lonlat(z)
    d = [GEOD.inv(np.full(lon.shape, o[0]), np.full(lat.shape, o[1]), lon, lat)[2] for o in (a, b)]
    expected = sum((di <= horizon).astype(int) for di in d)
    clear = np.all([np.abs(di - horizon) > 60 for di in d], axis=0) & np.all([np.abs(di - RADIUS) > 60 for di in d], axis=0)
    assert clear.sum() > 10_000
    assert (z.count[clear] == expected[clear]).mean() > 0.999
    assert set(np.unique(z.count[clear])) == {0, 1, 2}


def test_area_by_count_adds_up(plateau):
    obs = [east(OBS_LON, OBS_LAT, dx) for dx in (-2000, 0, 2000)]
    z = ztv(obs, RADIUS, observer_height=20, dem=plateau)
    areas = z.area_by_count()
    assert set(areas) == {1, 2, 3}
    assert sum(areas.values()) == pytest.approx(z.visible_km2)
    assert z.visible_km2 <= z.analysed_km2
    # Each observer's own visible area, summed, counts every cell once per observer it sees.
    assert len(z.observer_visible_km2) == 3
    assert sum(z.observer_visible_km2) == pytest.approx(sum(k * a for k, a in areas.items()))


def test_out_writes_counts_geotiff(plateau, tmp_path):
    out = tmp_path / "ztv.tif"
    z = ztv([(OBS_LON, OBS_LAT, 10)], 3000, dem=plateau, out=out)
    ds = gdal.Open(str(out))
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    assert band.GetNoDataValue() == 65535
    assert np.array_equal(arr[z.analysed], z.count[z.analysed])
    assert (arr[~z.analysed] == 65535).all()
    ds = None


@pytest.mark.parametrize("observers", [[], [(OBS_LON,)], [(200.0, 10.0)]])
def test_bad_observers_are_rejected(observers, plateau):
    with pytest.raises(ValueError):
        ztv(observers, RADIUS, dem=plateau)


def test_observer_without_elevation_is_rejected(plateau):
    far = east(OBS_LON, OBS_LAT, 30_000)  # off the synthetic plateau
    with pytest.raises(ValueError, match="no elevation"):
        ztv([(OBS_LON, OBS_LAT), far], 2000, dem=plateau)


def _png_rgba(png: bytes, tmp_path):
    p = tmp_path / "overlay.png"
    p.write_bytes(png)
    ds = gdal.Open(str(p))
    rgba = np.dstack([ds.GetRasterBand(b + 1).ReadAsArray() for b in range(4)])
    ds = None
    return rgba


def test_viewshed_overlay_lines_up_with_the_observer(plateau, tmp_path):
    """The overlay's corners bracket the observer, its centre is drawn as visible, and its
    corners (outside the circular range) are transparent."""
    v = viewshed(OBS_LON, OBS_LAT, 3000, observer_height=1.5, dem=plateau)
    png, corners = viewshed_overlay(v)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    (w, n), (e, _), (_, s), _ = corners
    assert w < OBS_LON < e and s < OBS_LAT < n
    rgba = _png_rgba(png, tmp_path)
    h, wd, _ = rgba.shape
    col = int((OBS_LON - w) / (e - w) * wd)
    # Web Mercator rows are not linear in latitude; use the projected y for the observer's row.
    to_m = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    yn, yo, ys = to_m.transform([w, w, w], [n, OBS_LAT, s])[1]
    row = int((yn - yo) / (yn - ys) * h)
    assert tuple(rgba[row, col]) == tuple(viewshed_lut()[VISIBLE])
    assert rgba[0, 0, 3] == 0 and rgba[-1, -1, 3] == 0


def test_ztv_classes_and_ramp(plateau, tmp_path):
    obs = [east(OBS_LON, OBS_LAT, dx) for dx in (-1000, 1000)]
    z = ztv(obs, 3000, observer_height=1.5, target_height=0, dem=plateau)
    cls = ztv_classes(z)
    assert set(np.unique(cls)) <= {0, 1, 2, 3, 255}
    assert (cls[~z.analysed & ~z.nodata] == 0).all()
    lut = ztv_lut(2)
    assert lut[0, 3] == 0 and lut[2, 3] > 0 and lut[3, 3] > 0
    assert lut[3, :3].astype(int).sum() < lut[2, :3].astype(int).sum()  # more visible -> darker
    png, _ = ztv_overlay(z)
    assert _png_rgba(png, tmp_path)[..., 3].any()
