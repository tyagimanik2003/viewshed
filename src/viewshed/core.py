"""The viewshed itself: lat/lon in, physically meaningful raster out."""
from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from osgeo import gdal

from .dem import observer_grid

gdal.UseExceptions()

# Output codes. Deliberately distinct, so visible cells can never be confused with no-data or out-of-range.
OUT_OF_RANGE = 0   # beyond max distance, inside min distance, or outside the sector
INVISIBLE = 1
MAYBE_VISIBLE = 2  # reserved for uncertainty-aware output
VISIBLE = 3
NODATA = 255       # no elevation data
CURVATURE_COEFF = 1 - 1 / 7  # GDAL default: curvature with standard refraction (k = 1/7)
# GDAL 3.13 renamed `gdal raster viewshed --dst-nodata` to `--output-nodata`.
_NODATA_ARG = "output_nodata" if tuple(map(int, gdal.__version__.split(".")[:2])) >= (3, 13) else "dst_nodata"


@dataclass
class ViewshedResult:
    """A viewshed on the observer-centred AEQD grid (metres; observer at x=0, y=0)."""
    codes: np.ndarray
    dem: np.ndarray
    geotransform: tuple
    crs_wkt: str
    observer: tuple[float, float]
    observer_height: float
    radius: float
    path: Path | None = None
    params: dict = field(default_factory=dict)

    @property
    def cell_area_km2(self) -> float:
        return abs(self.geotransform[1] * self.geotransform[5]) / 1e6

    @property
    def visible_km2(self) -> float:
        return float((self.codes == VISIBLE).sum() * self.cell_area_km2)

    @property
    def analysed_km2(self) -> float:
        return float(np.isin(self.codes, (VISIBLE, INVISIBLE, MAYBE_VISIBLE)).sum() * self.cell_area_km2)

    @property
    def visible_fraction(self) -> float:
        a = self.analysed_km2
        return self.visible_km2 / a if a else float("nan")

    def summary(self) -> str:
        lon, lat = self.observer
        return (f"Observer {lat:.4f}N {lon:.4f}E, height {self.observer_height} m, radius {self.radius/1000:g} km: "
                f"{self.visible_km2:.2f} km2 visible of {self.analysed_km2:.2f} km2 analysed "
                f"({100 * self.visible_fraction:.1f}%)")


def viewshed(lon: float, lat: float, radius: float = 5000.0, *, observer_height: float = 1.5,
             target_height: float = 0.0, start_angle: float | None = None, end_angle: float | None = None,
             min_distance: float | None = None, dem="glo30", resolution: float = 30.0,
             curvature_coeff: float = CURVATURE_COEFF, out: str | Path | None = None,
             engine: str = "gdal") -> ViewshedResult:
    """Compute what an observer standing at (lon, lat) can see.

    Parameters are in degrees (WGS84) and metres. Sector angles are geodesic azimuths,
    clockwise from true north; the sector runs clockwise from start_angle to end_angle
    (e.g. 300 -> 60 covers north). `dem` is "glo30" (Copernicus, fetched on demand),
    a raster path, or a list of raster paths such as local SRTM tiles.
    `engine` is "gdal" (fast, Wang et al. reference planes) or "r3" (exact reference, slow).
    """
    if engine not in {"gdal", "r3"}:
        raise ValueError("engine must be 'gdal' or 'r3'")
    if (start_angle is None) != (end_angle is None):
        raise ValueError("give both start_angle and end_angle, or neither")
    if min_distance is not None and not 0 <= min_distance < radius:
        raise ValueError("min_distance must be >= 0 and smaller than radius")

    work = Path(tempfile.mkdtemp(prefix="viewshed_"))
    dem_path = work / "dem_aeqd.tif"
    ds = observer_grid(lon, lat, radius, dem=dem, resolution_m=resolution, out_path=dem_path)
    dem_arr = ds.GetRasterBand(1).ReadAsArray()
    ds = None

    params = dict(target_height=target_height, start_angle=start_angle, end_angle=end_angle,
                  min_distance=min_distance, dem=str(dem), resolution=resolution,
                  curvature_coeff=curvature_coeff, engine=engine)
    if engine == "r3":
        from .reference import r3, sector_mask

        gt = gdal.Open(str(dem_path)).GetGeoTransform()
        codes = r3(dem_arr, resolution, observer_height, target_height, radius, curvature_coeff)
        region = sector_mask(codes.shape, resolution, start_angle, end_angle, min_distance, radius)
        codes[~region] = OUT_OF_RANGE
        return ViewshedResult(codes=codes, dem=dem_arr, geotransform=gt, crs_wkt=gdal.Open(str(dem_path)).GetProjection(),
                              observer=(lon, lat), observer_height=observer_height, radius=radius,
                              path=None, params=params)

    out_path = Path(out) if out else work / "viewshed.tif"
    kwargs = dict(
        # Observer height MUST go in position=[x, y, h]: in GDAL 3.13.3 the separate
        # -z/--height option is ignored and the default 2 m is used (docs/gdal_notes.md).
        input=str(dem_path), output=str(out_path), position=[0.0, 0.0, observer_height],
        target_height=target_height, max_distance=radius, curvature_coefficient=curvature_coeff,
        # maybe_visible_value is not passed: it only exists from GDAL 3.13, is unused without an
        # --sd-filename, and its default (2) already equals MAYBE_VISIBLE.
        visible_value=VISIBLE, invisible_value=INVISIBLE,
        out_of_range_value=OUT_OF_RANGE, overwrite=True,
    )
    kwargs[_NODATA_ARG] = NODATA
    if start_angle is not None:
        kwargs.update(start_angle=start_angle % 360, end_angle=end_angle % 360)
    if min_distance:
        kwargs.update(min_distance=min_distance)
    # gdal.Run is the route that honours sector and min-distance options (docs/gdal_notes.md).
    gdal.Run("raster", "viewshed", **kwargs).Finalize()

    vs = gdal.Open(str(out_path))
    codes = vs.GetRasterBand(1).ReadAsArray()
    gt, wkt = vs.GetGeoTransform(), vs.GetProjection()
    vs = None
    codes, dem_arr, gt = _align(codes, gt, dem_arr, dem_path)
    return ViewshedResult(codes=codes, dem=dem_arr, geotransform=gt, crs_wkt=wkt, observer=(lon, lat),
                          observer_height=observer_height, radius=radius, path=out_path, params=params)


def _align(codes, gt, dem_arr, dem_path):
    """GDAL clips the viewshed raster to the range/sector; place it back on the full DEM grid
    so codes and DEM share one geotransform."""
    dgt = gdal.Open(str(dem_path)).GetGeoTransform()
    full = np.full(dem_arr.shape, OUT_OF_RANGE, dtype=np.uint8)
    c0 = int(round((gt[0] - dgt[0]) / dgt[1]))
    r0 = int(round((gt[3] - dgt[3]) / dgt[5]))
    h, w = codes.shape
    full[r0:r0 + h, c0:c0 + w] = codes
    full[dem_arr <= -32767] = NODATA
    return full, dem_arr, dgt
