"""Engine specifications checked against geometry that has a known answer.

All inputs are synthetic and local, so these run offline.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from osgeo import gdal, osr

from conftest import ELEV, HORIZON_M, OBS_LAT, OBS_LON
from viewshed import INVISIBLE, NODATA, OUT_OF_RANGE, VISIBLE, viewshed

RADIUS = 5000.0


def grid_xy(res):
    gt = res.geotransform
    h, w = res.codes.shape
    c, r = np.meshgrid(np.arange(w) + 0.5, np.arange(h) + 0.5)
    x, y = gt[0] + c * gt[1], gt[3] + r * gt[5]
    return np.hypot(x, y), (np.degrees(np.arctan2(x, y)) + 360) % 360


def test_lat_lon_dem_gives_physical_viewshed(flat_geographic):
    """Handed a raw EPSG:4326 (lat/lon) DEM, the engine must still give the physical answer."""
    res = viewshed(OBS_LON, OBS_LAT, RADIUS, dem=flat_geographic)
    dist, _ = grid_xy(res)
    near = dist < 0.9 * HORIZON_M
    ring = (dist > 1.01 * HORIZON_M) & (dist < 0.99 * RADIUS)
    far = dist > 1.02 * RADIUS
    assert (res.codes[near] == VISIBLE).mean() > 0.99
    assert ring.sum() > 50 and (res.codes[ring] == INVISIBLE).mean() > 0.95
    assert (res.codes[far] == OUT_OF_RANGE).all()


def test_output_codes_are_distinct(flat_geographic):
    res = viewshed(OBS_LON, OBS_LAT, RADIUS, dem=flat_geographic)
    assert len({VISIBLE, INVISIBLE, OUT_OF_RANGE, NODATA}) == 4
    assert set(np.unique(res.codes)) <= {VISIBLE, INVISIBLE, OUT_OF_RANGE, NODATA}


def test_analysed_area_is_a_true_circle(flat_geographic):
    res = viewshed(OBS_LON, OBS_LAT, RADIUS, dem=flat_geographic)
    assert res.analysed_km2 == pytest.approx(math.pi * (RADIUS / 1000) ** 2, rel=0.01)


def test_sector_through_north_is_honoured(flat_geographic):
    res = viewshed(OBS_LON, OBS_LAT, RADIUS, dem=flat_geographic, start_angle=300, end_angle=60)
    dist, az = grid_xy(res)
    inside = ((az > 305) | (az < 55)) & (dist > 200) & (dist < 0.9 * HORIZON_M)
    outside = (az > 70) & (az < 290) & (dist > 200)
    assert (res.codes[inside] == VISIBLE).mean() > 0.99
    assert not (res.codes[outside] == VISIBLE).any()
    assert res.analysed_km2 == pytest.approx(math.pi * 25 * 120 / 360, rel=0.03)


def test_min_distance_is_honoured(flat_geographic):
    res = viewshed(OBS_LON, OBS_LAT, RADIUS, dem=flat_geographic, min_distance=1000)
    dist, _ = grid_xy(res)
    assert not (res.codes[dist < 950] == VISIBLE).any()
    assert (res.codes[(dist > 1100) & (dist < 4000)] == VISIBLE).mean() > 0.99


def test_missing_elevation_is_nodata_not_sea_level(tmp_path):
    """A hole in the DEM must come out as NODATA, never as visible/invisible 0 m terrain."""
    res_deg = 1 / 3600
    n = 600
    arr = np.full((n, n), ELEV, np.float32)
    arr[250:280, 310:340] = -9999  # void ~0.6-1.5 km north-east of the observer, as in raw SRTM
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    p = tmp_path / "void.tif"
    ds = gdal.GetDriverByName("GTiff").Create(str(p), n, n, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((OBS_LON - n / 2 * res_deg, res_deg, 0, OBS_LAT + n / 2 * res_deg, 0, -res_deg))
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).SetNoDataValue(-9999)
    ds.GetRasterBand(1).WriteArray(arr)
    ds = None
    res = viewshed(OBS_LON, OBS_LAT, 2000, dem=p)
    assert (res.codes == NODATA).sum() > 100
    assert not (res.dem[res.codes != NODATA] == 0).any()


def test_observer_outside_dem_is_rejected(flat_geographic):
    with pytest.raises(ValueError, match="no elevation data at the observer"):
        viewshed(OBS_LON + 1, OBS_LAT, 1000, dem=flat_geographic)


@pytest.mark.parametrize("kwargs", [dict(start_angle=10), dict(end_angle=10), dict(min_distance=6000)])
def test_bad_arguments_are_rejected(kwargs, flat_geographic):
    with pytest.raises(ValueError):
        viewshed(OBS_LON, OBS_LAT, RADIUS, dem=flat_geographic, **kwargs)
