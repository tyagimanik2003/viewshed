"""Zone of theoretical visibility (ZTV): from how many observers is each cell visible?

The usual application is visual impact assessment: the "observers" are wind turbines, masts or
buildings at their tip/roof height, and the target is a person (eye height about 2 m). A sightline
is symmetric, so the cells where a 2 m target is visible from a turbine tip are exactly the cells
from which a person can see that tip. The result counts, per cell, how many observers are visible.

All observers share one azimuthal-equidistant grid centred on their geographic centroid. Over the
spans this function accepts (see MAX_SPAN_M), the projection's distance error stays below about
1 part in 10,000.
"""
from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
from osgeo import gdal
from pyproj import Transformer

from .core import _NODATA_ARG, CURVATURE_COEFF, INVISIBLE, NODATA, OUT_OF_RANGE, VISIBLE
from .dem import aeqd_grid, aeqd_proj, geodesic_bbox

gdal.UseExceptions()

MAX_SPAN_M = 100_000.0  # max distance of any observer from the centroid
COUNT_NODATA = 65535    # GeoTIFF nodata for cells that were not analysed


@dataclass
class ZTVResult:
    """Per-cell visible-observer counts on an AEQD grid centred on the observers' centroid."""
    count: np.ndarray        # uint16: number of observers visible from the cell
    analysed: np.ndarray     # bool: within range of at least one observer, with elevation data
    dem: np.ndarray
    geotransform: tuple
    crs_wkt: str
    center: tuple[float, float]
    observers: list[tuple[float, float, float]]   # (lon, lat, height above ground)
    observer_xy: list[tuple[float, float]]        # observer positions in grid metres
    radius: float
    target_height: float
    path: Path | None = None
    params: dict = field(default_factory=dict)
    observer_visible_km2: list[float] = field(default_factory=list)  # area each observer is visible from

    @property
    def cell_area_km2(self) -> float:
        return abs(self.geotransform[1] * self.geotransform[5]) / 1e6

    @property
    def nodata(self) -> np.ndarray:
        return self.dem <= -32767

    @property
    def analysed_km2(self) -> float:
        return float(self.analysed.sum() * self.cell_area_km2)

    @property
    def visible_km2(self) -> float:
        """Area from which at least one observer is visible."""
        return float(((self.count > 0) & self.analysed).sum() * self.cell_area_km2)

    def area_by_count(self) -> dict[int, float]:
        """{k: km2 from which exactly k observers are visible}, for k = 1 .. number of observers."""
        c = self.count[self.analysed]
        return {k: float((c == k).sum() * self.cell_area_km2) for k in range(1, len(self.observers) + 1)}

    def summary(self) -> str:
        a = self.analysed_km2
        pct = 100 * self.visible_km2 / a if a else float("nan")
        return (f"{len(self.observers)} observers, radius {self.radius / 1000:g} km, target height "
                f"{self.target_height} m: at least one visible from {self.visible_km2:.2f} km2 of "
                f"{a:.2f} km2 analysed ({pct:.1f}%)")


def _normalise(observers, observer_height: float) -> list[tuple[float, float, float]]:
    out = []
    for o in observers:
        if len(o) == 2:
            lon, lat, h = o[0], o[1], observer_height
        elif len(o) == 3:
            lon, lat, h = o
        else:
            raise ValueError("each observer must be (lon, lat) or (lon, lat, height)")
        if not (-180 <= lon <= 180 and -90 <= lat <= 90):
            raise ValueError(f"observer ({lon}, {lat}) is not a valid longitude/latitude")
        out.append((float(lon), float(lat), float(h)))
    if not out:
        raise ValueError("at least one observer is required")
    return out


def _centroid(obs) -> tuple[float, float]:
    """Spherical centroid (mean of unit vectors), so groups across the antimeridian work."""
    lon, lat = np.radians([o[0] for o in obs]), np.radians([o[1] for o in obs])
    x, y, z = (np.cos(lat) * np.cos(lon)).mean(), (np.cos(lat) * np.sin(lon)).mean(), np.sin(lat).mean()
    return math.degrees(math.atan2(y, x)), math.degrees(math.atan2(z, math.hypot(x, y)))


def ztv(observers: Sequence[Sequence[float]], radius: float = 10_000.0, *, observer_height: float = 2.0,
        target_height: float = 2.0, dem="glo30", resolution: float = 30.0,
        curvature_coeff: float = CURVATURE_COEFF, out: str | Path | None = None) -> ZTVResult:
    """Count, for every cell, how many observers are visible from it.

    observers: (lon, lat) or (lon, lat, height) per observer, in degrees (WGS84) and metres above
    ground. A per-observer height overrides `observer_height`. For a wind farm, use the turbines'
    tip (or hub) heights as observer heights and a person's eye height as `target_height`.
    radius: analysis distance around each observer, metres.
    out: optional GeoTIFF of the counts (UInt16; cells not analysed are 65535).
    """
    obs = _normalise(observers, observer_height)
    if radius <= 0:
        raise ValueError("radius must be positive")
    lon0, lat0 = _centroid(obs)
    to_grid = Transformer.from_crs("EPSG:4326", aeqd_proj(lon0, lat0), always_xy=True)
    xy = [to_grid.transform(o[0], o[1]) for o in obs]
    if max(math.hypot(x, y) for x, y in xy) > MAX_SPAN_M:
        raise ValueError(f"observers spread too far apart: keep them within {MAX_SPAN_M / 1000:g} km of their centroid")

    # Cell edges at (k + 1/2) * resolution, as in observer_grid, with a two-cell margin.
    res = resolution
    xs, ys = [p[0] for p in xy], [p[1] for p in xy]
    edge_lo = lambda v: (math.floor((v - radius - 2 * res) / res) - 0.5) * res  # noqa: E731
    edge_hi = lambda v: (math.ceil((v + radius + 2 * res) / res) + 0.5) * res  # noqa: E731
    bounds = (edge_lo(min(xs)), edge_lo(min(ys)), edge_hi(max(xs)), edge_hi(max(ys)))

    boxes = [geodesic_bbox(o[0], o[1], radius * 1.05) for o in obs]
    bbox = (min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes))

    with tempfile.TemporaryDirectory(prefix="viewshed_ztv_") as tmp:
        dem_arr, gt, wkt, count, analysed, per_obs = _compute(obs, xy, lon0, lat0, bounds, bbox, Path(tmp), radius,
                                                     target_height, dem, res, curvature_coeff)

    out_path = Path(out) if out else None
    if out_path:
        _write_counts(out_path, np.where(analysed, count, COUNT_NODATA).astype(np.uint16), gt, wkt)
    return ZTVResult(count=count, analysed=analysed, dem=dem_arr, geotransform=gt, crs_wkt=wkt,
                     center=(lon0, lat0), observers=obs, observer_xy=xy, radius=radius,
                     target_height=target_height, path=out_path,
                     params=dict(dem=str(dem), resolution=res, curvature_coeff=curvature_coeff),
                     observer_visible_km2=[n * res * res / 1e6 for n in per_obs])


def _compute(obs, xy, lon0, lat0, bounds, bbox, work: Path, radius, target_height, dem, res, curvature_coeff):
    dem_path = work / "dem_aeqd.tif"
    ds = aeqd_grid(lon0, lat0, bounds, bbox, dem, res, dem_path)
    dem_arr = ds.GetRasterBand(1).ReadAsArray()
    gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
    ds = None
    nodata = dem_arr <= -32767

    count = np.zeros(dem_arr.shape, np.uint16)
    analysed = np.zeros(dem_arr.shape, bool)
    per_obs: list[int] = []
    for i, ((lon, lat, h), (x, y)) in enumerate(zip(obs, xy)):
        col, row = int((x - gt[0]) // gt[1]), int((y - gt[3]) // gt[5])
        if nodata[row, col]:
            raise ValueError(f"observer {i} at ({lon:.5f}, {lat:.5f}) has no elevation data")
        vs_path = work / f"vs_{i}.tif"
        kwargs = dict(
            # Height goes in position=[x, y, h]: the separate --height option is ignored by GDAL
            # 3.12-3.13 (docs/gdal_notes.md).
            input=str(dem_path), output=str(vs_path), position=[x, y, h], target_height=target_height,
            max_distance=radius, curvature_coefficient=curvature_coeff, visible_value=VISIBLE,
            invisible_value=INVISIBLE, out_of_range_value=OUT_OF_RANGE, overwrite=True,
        )
        kwargs[_NODATA_ARG] = NODATA
        gdal.Run("raster", "viewshed", **kwargs).Finalize()
        vs = gdal.Open(str(vs_path))
        codes = vs.GetRasterBand(1).ReadAsArray()
        vgt = vs.GetGeoTransform()
        vs = None
        # GDAL crops the output to the observer's range; place it back on the shared grid.
        c0, r0 = int(round((vgt[0] - gt[0]) / gt[1])), int(round((vgt[3] - gt[3]) / gt[5]))
        win = (slice(r0, r0 + codes.shape[0]), slice(c0, c0 + codes.shape[1]))
        vis = (codes == VISIBLE) & ~nodata[win]
        count[win] += vis
        per_obs.append(int(vis.sum()))
        analysed[win] |= np.isin(codes, (VISIBLE, INVISIBLE))
        vs_path.unlink(missing_ok=True)
    count[nodata] = 0
    analysed &= ~nodata
    return dem_arr, gt, wkt, count, analysed, per_obs


def _write_counts(path: Path, arr: np.ndarray, gt, wkt) -> None:
    ds = gdal.GetDriverByName("GTiff").Create(str(path), arr.shape[1], arr.shape[0], 1, gdal.GDT_UInt16,
                                              options=["COMPRESS=DEFLATE"])
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(COUNT_NODATA)
    band.WriteArray(arr)
    ds.FlushCache()
    ds = None
