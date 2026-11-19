"""ViewshedBench engine adapters.

Every adapter takes the SAME observer-centred AEQD DEM (odd-sized, observer on the centre
cell) and returns codes on that grid: VISIBLE / INVISIBLE / OUT_OF_RANGE / NODATA.
Settings are matched across engines; what each engine cannot express is recorded in
ENGINE_NOTES rather than silently ignored.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr

from viewshed import INVISIBLE, NODATA, OUT_OF_RANGE, VISIBLE
from viewshed.reference import r3

gdal.UseExceptions()
K_REFRACTION = 1 / 7  # standard refraction; curved scenario uses curvature coeff 1 - k

ENGINE_NOTES = {
    "r3": "Exact R3 (this package). Declared terrain model in viewshed/reference.py.",
    "gdal": "gdal raster viewshed (Wang et al. 2000 reference planes), via gdal.Run, height in position.",
    "grass": "GRASS r.viewshed (Haverkort, Toma & Zhuang 2009 radial sweep); -c curvature, -r refraction 1/7.",
    "whitebox": "WhiteboxTools Viewshed; no curvature option, no max distance (masked afterwards).",
    "xrspatial": "xarray-spatial viewshed (CPU); no curvature option in the documented API, no max distance.",
}
SUPPORTS_CURVATURE = {"r3", "gdal", "grass"}


def _read(path):
    ds = gdal.Open(str(path))
    return ds.GetRasterBand(1).ReadAsArray(), ds.GetGeoTransform(), ds.GetProjection()


def _dist(shape, cell):
    n = shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    return np.hypot(xx - n // 2, yy - n // 2) * cell


def finish(vis_bool, dem, cell, radius, nodata=-32768.0):
    codes = np.where(vis_bool, VISIBLE, INVISIBLE).astype(np.uint8)
    codes[dem <= nodata + 0.5] = NODATA
    codes[_dist(dem.shape, cell) > radius] = OUT_OF_RANGE
    n = dem.shape[0]
    codes[n // 2, n // 2] = VISIBLE
    return codes


def run_r3(dem_path, radius, h, curved):
    dem, gt, _ = _read(dem_path)
    return r3(dem, gt[1], h, 0.0, radius, (1 - K_REFRACTION) if curved else 0.0)


def run_gdal(dem_path, radius, h, curved):
    dem, gt, _ = _read(dem_path)
    out = Path(tempfile.mkdtemp()) / "vs.tif"
    kw = dict(input=str(dem_path), output=str(out), position=[0.0, 0.0, h],
              max_distance=radius, curvature_coefficient=(1 - K_REFRACTION) if curved else 0.0,
              visible_value=VISIBLE, invisible_value=INVISIBLE, out_of_range_value=OUT_OF_RANGE, overwrite=True)
    # GDAL 3.13 renamed --dst-nodata to --output-nodata
    nodata_kw = "output_nodata" if tuple(map(int, gdal.__version__.split(".")[:2])) >= (3, 13) else "dst_nodata"
    gdal.Run("raster", "viewshed", **kw, **{nodata_kw: NODATA}).Finalize()
    a, ogt, _ = _read(out)
    full = np.full(dem.shape, OUT_OF_RANGE, np.uint8)
    c0, r0 = int(round((ogt[0] - gt[0]) / gt[1])), int(round((ogt[3] - gt[3]) / gt[5]))
    full[r0:r0 + a.shape[0], c0:c0 + a.shape[1]] = a
    return finish(full == VISIBLE, dem, gt[1], radius)


def _bin(name):
    exe = shutil.which(name)
    if exe:
        return exe
    for d in (Path(sys.prefix) / "Library" / "bin", Path(sys.prefix) / "Scripts", Path(sys.prefix)):
        for ext in ("", ".exe", ".bat", ".cmd"):
            if (d / f"{name}{ext}").exists():
                return str(d / f"{name}{ext}")
    raise FileNotFoundError(name)


def run_grass(dem_path, radius, h, curved):
    dem, gt, _ = _read(dem_path)
    work = Path(tempfile.mkdtemp())
    out = work / "vs.tif"
    flags = "cbr" if curved else "b"
    job = work / "job.py"
    job.write_text(textwrap.dedent(f"""
        import grass.script as gs
        gs.run_command("r.in.gdal", input=r"{dem_path}", output="dem", flags="o", quiet=True)
        gs.run_command("g.region", raster="dem", quiet=True)
        gs.run_command("r.viewshed", input="dem", output="vs", coordinates=(0, 0),
                       observer_elevation={h}, target_elevation=0, max_distance={radius},
                       refraction_coeff={K_REFRACTION}, flags="{flags}", quiet=True, overwrite=True)
        gs.run_command("r.out.gdal", input="vs", output=r"{out}", format="GTiff", type="Byte",
                       nodata=255, flags="c", quiet=True, overwrite=True)
    """), encoding="utf-8")
    subprocess.run([_bin("grass"), "--tmp-project", str(dem_path), "--exec", sys.executable, str(job)],
                   check=True, capture_output=True, text=True)
    a, _, _ = _read(out)
    return finish(a == 1, dem, gt[1], radius)


def run_whitebox(dem_path, radius, h, curved):
    dem, gt, wkt = _read(dem_path)
    work = Path(tempfile.mkdtemp())
    pts = work / "obs.shp"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    ds = drv.CreateDataSource(str(pts))
    srs = osr.SpatialReference(); srs.ImportFromWkt(wkt)
    lyr = ds.CreateLayer("obs", srs, ogr.wkbPoint)
    f = ogr.Feature(lyr.GetLayerDefn()); g = ogr.Geometry(ogr.wkbPoint); g.AddPoint(0.0, 0.0); f.SetGeometry(g)
    lyr.CreateFeature(f); ds = None
    out = work / "vs.tif"
    subprocess.run([_bin("whitebox_tools"), "-r=Viewshed", f"--wd={work}", f"--dem={dem_path}",
                    f"--stations={pts}", f"--output={out}", f"--height={h}"],
                   check=True, capture_output=True, text=True)
    a, _, _ = _read(out)
    return finish(a > 0, dem, gt[1], radius)


def run_xrspatial(dem_path, radius, h, curved):
    import xarray as xr
    from xrspatial import viewshed as xv

    dem, gt, _ = _read(dem_path)
    n = dem.shape[0]
    xs = gt[0] + (np.arange(n) + 0.5) * gt[1]
    ys = gt[3] + (np.arange(n) + 0.5) * gt[5]
    da = xr.DataArray(dem.astype(np.float64), dims=("y", "x"), coords={"y": ys, "x": xs})
    res = xv(da, x=0.0, y=0.0, observer_elev=h, target_elev=0.0)
    vis = np.asarray(res.values) != -1  # INVISIBLE sentinel is -1; visible cells hold the view angle
    return finish(vis, dem, gt[1], radius)


ENGINES = {"r3": run_r3, "gdal": run_gdal, "grass": run_grass, "whitebox": run_whitebox, "xrspatial": run_xrspatial}
