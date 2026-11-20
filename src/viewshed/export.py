"""Export results for web maps: a coloured PNG in Web Mercator plus its corner coordinates.

Web maps (MapLibre, Leaflet, OpenLayers) draw image overlays in Web Mercator (EPSG:3857). The
analysis grids here are azimuthal-equidistant, so the classes are reprojected (nearest neighbour,
so no class is invented) before colouring; the PNG then lines up exactly with the basemap.
"""
from __future__ import annotations

import uuid

import numpy as np
from osgeo import gdal
from pyproj import Transformer

from .core import INVISIBLE, NODATA, OUT_OF_RANGE, VISIBLE, ViewshedResult

gdal.UseExceptions()

_TO_LONLAT = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)


def web_overlay(classes: np.ndarray, geotransform, crs_wkt: str, lut: np.ndarray) -> tuple[bytes, list[list[float]]]:
    """Reproject a uint8 class raster to Web Mercator and colour it with `lut` (256 x 4 RGBA, uint8).

    Class 0 must mean "draw nothing": it is also used for cells outside the source grid.
    Returns (png_bytes, corners), where corners are [lon, lat] for the top-left, top-right,
    bottom-right and bottom-left of the image, the order MapLibre's image source expects.
    """
    lut = np.asarray(lut, np.uint8)
    if lut.shape != (256, 4):
        raise ValueError("lut must have shape (256, 4)")
    if lut[0, 3] != 0:
        raise ValueError("class 0 must be fully transparent")
    h, w = classes.shape
    src = gdal.GetDriverByName("MEM").Create("", w, h, 1, gdal.GDT_Byte)
    src.SetGeoTransform(geotransform)
    src.SetProjection(crs_wkt)
    src.GetRasterBand(1).WriteArray(classes.astype(np.uint8))
    # Class 0 is nodata on both sides; otherwise GDAL remaps source 0 to 1 to keep it distinct.
    src.GetRasterBand(1).SetNoDataValue(0)
    merc = gdal.Warp("", src, format="MEM", dstSRS="EPSG:3857", resampleAlg="near", dstNodata=0)
    src = None
    warped = merc.GetRasterBand(1).ReadAsArray()
    mgt = merc.GetGeoTransform()
    merc = None

    rgba = lut[warped]
    mem = gdal.GetDriverByName("MEM").Create("", rgba.shape[1], rgba.shape[0], 4, gdal.GDT_Byte)
    for b in range(4):
        mem.GetRasterBand(b + 1).WriteArray(rgba[:, :, b])
    name = f"/vsimem/overlay_{uuid.uuid4().hex}.png"
    gdal.GetDriverByName("PNG").CreateCopy(name, mem, options=["ZLEVEL=6"])
    mem = None
    try:
        f = gdal.VSIFOpenL(name, "rb")
        png = gdal.VSIFReadL(1, gdal.VSIStatL(name).size, f)
        gdal.VSIFCloseL(f)
    finally:
        for p in (name, name + ".aux.xml"):
            if gdal.VSIStatL(p) is not None:
                gdal.Unlink(p)

    x0, y0 = mgt[0], mgt[3]
    x1, y1 = x0 + warped.shape[1] * mgt[1], y0 + warped.shape[0] * mgt[5]
    corners = [list(_TO_LONLAT.transform(x, y)) for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))]
    return bytes(png), corners


def viewshed_lut(visible=(250, 158, 28, 175), invisible=(20, 30, 45, 70), nodata=(88, 80, 236, 150)) -> np.ndarray:
    """Colours for ViewshedResult.codes: visible amber, invisible dim, no data red, rest transparent."""
    lut = np.zeros((256, 4), np.uint8)
    lut[VISIBLE], lut[INVISIBLE], lut[NODATA] = visible, invisible, nodata
    lut[OUT_OF_RANGE] = 0
    return lut


def viewshed_overlay(result: ViewshedResult, lut: np.ndarray | None = None) -> tuple[bytes, list[list[float]]]:
    """PNG overlay of a single viewshed for web maps."""
    return web_overlay(result.codes, result.geotransform, result.crs_wkt, viewshed_lut() if lut is None else lut)


def ztv_classes(result) -> np.ndarray:
    """ZTV counts as overlay classes: 0 not analysed, 1 + count (capped at 254), 255 no data."""
    cls = np.where(result.analysed, np.minimum(result.count, 253).astype(np.int32) + 1, 0)
    cls[result.nodata] = 255
    return cls.astype(np.uint8)


def ztv_lut(n_observers: int, nodata=(88, 80, 236, 150)) -> np.ndarray:
    """Sequential ramp for ZTV counts: none visible dim, then light yellow to deep red as the count grows."""
    stops = np.array([[255, 237, 160], [254, 178, 76], [240, 59, 32], [128, 0, 38]], float)
    lut = np.zeros((256, 4), np.uint8)
    lut[1] = (20, 30, 45, 70)  # analysed, nothing visible
    n = max(1, min(n_observers, 253))
    for k in range(1, n + 1):
        t = 0.0 if n == 1 else (k - 1) / (n - 1)
        pos = t * (len(stops) - 1)
        i = min(int(pos), len(stops) - 2)
        rgb = stops[i] + (stops[i + 1] - stops[i]) * (pos - i)
        lut[k + 1] = (*np.round(rgb).astype(int), 185)
    lut[255] = nodata
    return lut


def ztv_overlay(result, lut: np.ndarray | None = None) -> tuple[bytes, list[list[float]]]:
    """PNG overlay of a ZTV for web maps."""
    lut = ztv_lut(len(result.observers)) if lut is None else lut
    return web_overlay(ztv_classes(result), result.geotransform, result.crs_wkt, lut)
