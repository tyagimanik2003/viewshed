"""Shared fixtures: synthetic flat DEMs in geographic and projected CRSs.

A flat plateau (4000 m, like Lahaul) makes the expected viewshed analytic:
with observer height h, target height 0 and GDAL's curvature coefficient
c = 1 - 1/7, a cell at ground distance d is visible iff d <= sqrt(2 R h / c).
For h = 1.5 m that horizon is about 4.72 km.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, osr

gdal.UseExceptions()
sys.path.insert(0, str(Path(__file__).parent))

OBS_LON, OBS_LAT = 77.015, 32.329  # a Himalayan test location (Kangra-Lahaul ridge)
ELEV = 4000.0
OBS_H = 1.5
CURV = 1 - 1 / 7  # GDAL default: curvature with standard refraction
R_EARTH = 6371008.8  # GDAL uses the CRS's semi-major axis; difference is << 1 cell here
HORIZON_M = math.sqrt(2 * R_EARTH * OBS_H / CURV)


def _write(path: Path, arr: np.ndarray, gt, wkt: str, nodata=None) -> Path:
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    band.WriteArray(arr)
    ds.FlushCache()
    ds = None
    return path


@pytest.fixture
def flat_geographic(tmp_path) -> Path:
    """SRTM-like EPSG:4326 grid, 1 arc-second cells, +-7 km around the observer."""
    res = 1 / 3600
    half_lon = 7000 / (111320 * math.cos(math.radians(OBS_LAT)))
    half_lat = 7000 / 110950
    ncol, nrow = int(2 * half_lon / res), int(2 * half_lat / res)
    gt = (OBS_LON - ncol / 2 * res, res, 0, OBS_LAT + nrow / 2 * res, 0, -res)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    return _write(tmp_path / "flat_4326.tif", np.full((nrow, ncol), ELEV, np.float32), gt, srs.ExportToWkt())


@pytest.fixture
def flat_aeqd(tmp_path) -> Path:
    """Same plateau in an azimuthal-equidistant projection centred on the observer (30 m cells)."""
    res, half = 30.0, 7000.0
    n = int(2 * half / res) + 1
    gt = (-n / 2 * res, res, 0, n / 2 * res, 0, -res)
    srs = osr.SpatialReference()
    srs.ImportFromProj4(f"+proj=aeqd +lat_0={OBS_LAT} +lon_0={OBS_LON} +datum=WGS84 +units=m +no_defs")
    return _write(tmp_path / "flat_aeqd.tif", np.full((n, n), ELEV, np.float32), gt, srs.ExportToWkt())


@pytest.fixture
def flat_utm(tmp_path) -> Path:
    """Same plateau in UTM 43N (EPSG:32643), 30 m cells."""
    from pyproj import Transformer

    x0, y0 = Transformer.from_crs(4326, 32643, always_xy=True).transform(OBS_LON, OBS_LAT)
    res, n = 30.0, 401
    gt = (x0 - n / 2 * res, res, 0, y0 + n / 2 * res, 0, -res)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32643)
    return _write(tmp_path / "flat_utm.tif", np.full((n, n), ELEV, np.float32), gt, srs.ExportToWkt())


def ground_distance_grid(path: Path) -> np.ndarray:
    """Geodesic distance (m) from the observer to every cell centre of a raster."""
    from pyproj import CRS, Geod, Transformer

    ds = gdal.Open(str(path))
    gt = ds.GetGeoTransform()
    cols, rows = np.meshgrid(np.arange(ds.RasterXSize) + 0.5, np.arange(ds.RasterYSize) + 0.5)
    x = gt[0] + cols * gt[1]
    y = gt[3] + rows * gt[5]
    to_ll = Transformer.from_crs(CRS.from_wkt(ds.GetProjection()), 4326, always_xy=True)
    lon, lat = to_ll.transform(x, y)
    _, _, dist = Geod(ellps="WGS84").inv(np.full(lon.shape, OBS_LON), np.full(lat.shape, OBS_LAT), lon, lat)
    return dist


def observer_xy(path: Path) -> tuple[float, float]:
    """Observer position in the raster's own CRS units."""
    from pyproj import CRS, Transformer

    ds = gdal.Open(str(path))
    t = Transformer.from_crs(4326, CRS.from_wkt(ds.GetProjection()), always_xy=True)
    return t.transform(OBS_LON, OBS_LAT)
