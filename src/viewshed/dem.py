"""Where elevation data comes from, and how it gets into an observer-centred metric grid.

Viewsheds computed directly on latitude/longitude grids mix degrees and metres. Here every analysis grid is an
azimuthal-equidistant (AEQD) projection centred on the observer: distances and azimuths
from the observer are true, so circles stay circles and sector bearings stay geodesic.
"""
from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Iterable, Sequence

from osgeo import gdal, osr
from pyproj import Geod

gdal.UseExceptions()
GEOD = Geod(ellps="WGS84")
DEM_NODATA = -32768.0

GLO30_URL = (
    "/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{alat:02d}_00_{ew}{alon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{alat:02d}_00_{ew}{alon:03d}_00_DEM.tif"
)


def aeqd_proj(lon: float, lat: float) -> str:
    """PROJ string of the azimuthal-equidistant projection centred on (lon, lat)."""
    return f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs"


def aeqd_srs(lon: float, lat: float) -> osr.SpatialReference:
    srs = osr.SpatialReference()
    srs.ImportFromProj4(aeqd_proj(lon, lat))
    return srs


def geodesic_bbox(lon: float, lat: float, radius_m: float) -> tuple[float, float, float, float]:
    """(west, south, east, north) in degrees enclosing a geodesic circle (accounts for cos(latitude))."""
    lons, lats = [], []
    for az in range(0, 360, 15):
        lo, la, _ = GEOD.fwd(lon, lat, az, radius_m)
        lons.append(lo)
        lats.append(la)
    return min(lons), min(lats), max(lons), max(lats)


def _tiles_1deg(bbox) -> Iterable[tuple[int, int]]:
    w, s, e, n = bbox
    for la in range(math.floor(s), math.floor(n) + 1):
        for lo in range(math.floor(w), math.floor(e) + 1):
            yield la, lo


def cache_dir() -> Path | None:
    """Tile cache: $VIEWSHED_CACHE if set, else ~/.cache/viewshed. Set VIEWSHED_CACHE=off to disable."""
    v = os.environ.get("VIEWSHED_CACHE")
    if v and v.lower() == "off":
        return None
    d = Path(v) if v else Path.home() / ".cache" / "viewshed"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download(url: str, dest: Path) -> bool:
    """Fetch one tile to dest atomically. False if the tile does not exist (e.g. open ocean)."""
    import urllib.error
    import urllib.request

    part = dest.with_suffix(".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(part, "wb") as f:
            while chunk := r.read(1 << 20):
                f.write(chunk)
    except urllib.error.HTTPError as e:
        part.unlink(missing_ok=True)
        if e.code in (403, 404):
            return False
        raise
    part.replace(dest)
    return True


def glo30_sources(bbox) -> list[str]:
    """Copernicus GLO-30 tiles covering bbox: cached locally when a cache is enabled,
    otherwise read in place over HTTP. Tiles are fetched at run time, never bundled."""
    cache = cache_dir()
    srcs = []
    for la, lo in _tiles_1deg(bbox):
        url = GLO30_URL.format(ns="N" if la >= 0 else "S", alat=abs(la), ew="E" if lo >= 0 else "W", alon=abs(lo))
        if cache is not None:
            local = cache / "glo30" / url.rsplit("/", 1)[1]
            missing = local.with_suffix(".missing")
            if local.exists():
                srcs.append(str(local))
            elif not missing.exists():
                local.parent.mkdir(parents=True, exist_ok=True)
                if _download(url.removeprefix("/vsicurl/"), local):
                    srcs.append(str(local))
                else:
                    missing.touch()
            continue
        try:
            gdal.Open(url)  # some tiles (open ocean) do not exist
            srcs.append(url)
        except RuntimeError:
            pass
    if not srcs:
        raise ValueError(f"No Copernicus GLO-30 tiles found for bbox {bbox}")
    return srcs


def resolve_sources(dem: str | Path | Sequence[str | Path], bbox) -> list[str]:
    """dem = 'glo30', one raster path, or a list of raster paths (e.g. local SRTM tiles)."""
    if isinstance(dem, str) and dem.lower() in {"glo30", "copernicus"}:
        return glo30_sources(bbox)
    if isinstance(dem, (str, Path)):
        return [str(dem)]
    return [str(p) for p in dem]


def observer_grid(lon: float, lat: float, radius_m: float, dem="glo30", resolution_m: float = 30.0,
                  out_path: str | Path | None = None) -> gdal.Dataset:
    """Mosaic the sources and reproject them to an AEQD grid centred on the observer.

    Bilinear resampling (reprojection changes elevations; bilinear is the usual mitigation).
    The grid covers radius_m plus a two-cell margin, has an odd number of cells per side, and
    the observer (0, 0) sits exactly on the centre of the middle cell.
    """
    bbox = geodesic_bbox(lon, lat, radius_m * 1.05)
    half = math.ceil((radius_m + 2 * resolution_m) / resolution_m) * resolution_m + resolution_m / 2
    return aeqd_grid(lon, lat, (-half, -half, half, half), bbox, dem, resolution_m, out_path)


def aeqd_grid(lon0: float, lat0: float, bounds: tuple[float, float, float, float], bbox_deg,
              dem="glo30", resolution_m: float = 30.0, out_path: str | Path | None = None) -> gdal.Dataset:
    """Mosaic the DEM sources covering bbox_deg (west, south, east, north) and warp them to the
    AEQD grid centred on (lon0, lat0) with the given (xmin, ymin, xmax, ymax) bounds in metres."""
    sources = resolve_sources(dem, bbox_deg)
    vrt = gdal.BuildVRT("", sources)
    dst = str(out_path) if out_path else ""
    return gdal.Warp(
        dst, vrt, format="GTiff" if out_path else "MEM",
        dstSRS=aeqd_srs(lon0, lat0).ExportToWkt(),
        outputBounds=bounds, xRes=resolution_m, yRes=resolution_m,
        resampleAlg="bilinear", dstNodata=DEM_NODATA, outputType=gdal.GDT_Float32,
        multithread=True,
    )
