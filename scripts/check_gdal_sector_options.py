"""Check (GDAL issue #13451): do sector / min-distance viewshed options work from Python?

GDAL 3.12 added --start-angle/--end-angle/--min-distance (and pitch limits) to
`gdal raster viewshed`. Issue #13451 reports that the C API -- which the classic
Python binding gdal.ViewshedGenerate calls -- drops the angle and pitch settings.
If that holds, the rebuilt package cannot rely on ViewshedGenerate for sectors.

Test: flat 4000 m plateau in an observer-centred AEQD projection, observer 1.5 m,
max distance 4 km, sector 0-90 deg (clockwise from north = north-east quadrant),
min distance 1 km. A route "honours" the options when
  * cells in the NE quadrant between 1.2 and 3.8 km are visible,
  * no cell in the SW quadrant is visible (masked or clipped away), and
  * no cell closer than 0.8 km is visible.

Routes tried:
  A  gdal.ViewshedGenerate(...)              classic Python binding (no angle arguments)
  A2 gdal.ViewshedGenerate(..., options=...) same, passing the options as key=value strings
  B  `gdal raster viewshed` CLI via subprocess
  C  gdal.Run("raster", "viewshed", ...)     Python API for the new CLI algorithms

Usage:  python scripts/check_gdal_sector_options.py [--out sector_check.md]
"""
from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

LON, LAT, ELEV, H = 77.015, 32.329, 4000.0, 1.5
MAXD, MIND, START, END = 4000.0, 1000.0, 0.0, 90.0
VIS, OOR = 255, 100


def make_dem(path: Path) -> None:
    res, n = 30.0, 301
    srs = osr.SpatialReference()
    srs.ImportFromProj4(f"+proj=aeqd +lat_0={LAT} +lon_0={LON} +datum=WGS84 +units=m +no_defs")
    ds = gdal.GetDriverByName("GTiff").Create(str(path), n, n, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((-n / 2 * res, res, 0, n / 2 * res, 0, -res))
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(np.full((n, n), ELEV, np.float32))
    ds = None


def evaluate(out: Path) -> dict:
    ds = gdal.Open(str(out))
    arr = ds.GetRasterBand(1).ReadAsArray()
    gt = ds.GetGeoTransform()
    c, r = np.meshgrid(np.arange(ds.RasterXSize) + 0.5, np.arange(ds.RasterYSize) + 0.5)
    x, y = gt[0] + c * gt[1], gt[3] + r * gt[5]
    dist = np.hypot(x, y)
    az = (np.degrees(np.arctan2(x, y)) + 360) % 360
    ne = (az > 5) & (az < 85) & (dist > 1200) & (dist < 3800)
    sw = (az > 185) & (az < 265) & (dist > 1200) & (dist < 3800)
    close = dist < 800
    res = {
        "ne_visible_frac": float((arr[ne] == VIS).mean()) if ne.any() else float("nan"),
        "sw_cells_in_output": int(sw.sum()),
        "sw_visible": int((arr[sw] == VIS).sum()),
        "close_visible": int((arr[close] == VIS).sum()),
        "output_size": f"{ds.RasterXSize}x{ds.RasterYSize}",
    }
    res["honours_sector"] = res["sw_visible"] == 0 and res["ne_visible_frac"] > 0.99
    res["honours_min_distance"] = res["close_visible"] == 0
    return res


def route_a(dem: Path, out: Path, options=None) -> None:
    ds = gdal.Open(str(dem))
    kw = dict(
        srcBand=ds.GetRasterBand(1), driverName="GTiff", targetRasterName=str(out), creationOptions=[],
        observerX=0.0, observerY=0.0, observerHeight=H, targetHeight=0.0, visibleVal=VIS, invisibleVal=0,
        outOfRangeVal=OOR, noDataVal=254, dfCurvCoeff=1 - 1 / 7, mode=2, maxDistance=MAXD,
    )
    if options is not None:
        kw["options"] = options
    r = gdal.ViewshedGenerate(**kw)
    r = None


def route_b(dem: Path, out: Path) -> None:
    exe = shutil.which("gdal") or str(Path(sys.prefix) / "Library" / "bin" / "gdal.exe")
    cmd = [exe, "raster", "viewshed", "--position", f"0,0,{H}", "--max-distance", str(MAXD),
           "--min-distance", str(MIND), "--start-angle", str(START), "--end-angle", str(END),
           "--out-of-range-value", str(OOR), "--curvature-coefficient", str(1 - 1 / 7),
           "--overwrite", "-q", str(dem), str(out)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def route_c(dem: Path, out: Path) -> None:
    gdal.Run("raster", "viewshed", input=str(dem), output=str(out), position=[0.0, 0.0, H],
             max_distance=MAXD, min_distance=MIND, start_angle=START, end_angle=END,
             out_of_range_value=OOR, curvature_coefficient=1 - 1 / 7, overwrite=True).Finalize()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    args = ap.parse_args()
    tmp = Path(tempfile.mkdtemp())
    dem = tmp / "flat_aeqd.tif"
    make_dem(dem)
    routes = {
        "A  gdal.ViewshedGenerate (no angle args exist)": lambda o: route_a(dem, o),
        "A2 gdal.ViewshedGenerate options=[START_ANGLE=...]": lambda o: route_a(
            dem, o, [f"START_ANGLE={START}", f"END_ANGLE={END}", f"MIN_DISTANCE={MIND}"]),
        "B  `gdal raster viewshed` CLI": lambda o: route_b(dem, o),
        "C  gdal.Run('raster','viewshed')": lambda o: route_c(dem, o),
    }
    rows = []
    for i, (name, fn) in enumerate(routes.items()):
        out = tmp / f"vs_{i}.tif"
        try:
            fn(out)
            r = evaluate(out)
            rows.append((name, "yes" if r["honours_sector"] else "**no**",
                         "yes" if r["honours_min_distance"] else "**no**",
                         f"NE visible {r['ne_visible_frac']:.2f}; SW visible cells {r['sw_visible']} "
                         f"of {r['sw_cells_in_output']}; <0.8 km visible {r['close_visible']}; out {r['output_size']}"))
        except Exception as e:  # report, don't crash: a refused route is a result
            rows.append((name, "error", "error", f"{type(e).__name__}: {str(e).splitlines()[0][:120]}"))

    sig = (gdal.ViewshedGenerate.__doc__ or "").strip().splitlines()[0]
    lines = [
        "# GDAL sector-option check (issue #13451)",
        "",
        f"GDAL {gdal.__version__}. Flat AEQD plateau, observer 1.5 m, max 4 km, "
        f"sector {START:.0f}-{END:.0f} deg clockwise from north, min distance {MIND:.0f} m.",
        "",
        "| Route | Sector honoured | Min distance honoured | Evidence |",
        "|---|---|---|---|",
        *[f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows],
        "",
        f"`ViewshedGenerate` signature: `{sig}`",
    ]
    text = "\n".join(lines)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
