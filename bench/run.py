"""Run ViewshedBench: every engine vs the exact R3 reference, on identical grids.

    python -m bench.run                     # all sites, radii 5 and 10 km
    python -m bench.run --sites 3 --radii 5000

Resumable: results append to outputs/bench/results.jsonl and finished runs are skipped.
Terrain grids and engine outputs are kept under data/bench/ (git-ignored).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("VIEWSHED_CACHE", str(ROOT / "data" / "cache"))

from viewshed import INVISIBLE, VISIBLE  # noqa: E402
from viewshed.dem import observer_grid  # noqa: E402

from .engines import ENGINES, SUPPORTS_CURVATURE  # noqa: E402
from .sites import SITES  # noqa: E402

OUT = ROOT / "outputs" / "bench"
DATA = ROOT / "data" / "bench"
H, RES = 1.5, 30.0
BANDS = [(0, 1000), (1000, 2500), (2500, 5000), (5000, 10000)]


def slug(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def metrics(codes, ref, dem, cell):
    n = ref.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    dist = np.hypot(xx - n // 2, yy - n // 2) * cell
    m = np.isin(ref, (VISIBLE, INVISIBLE)) & np.isin(codes, (VISIBLE, INVISIBLE))
    ev, rv = (codes == VISIBLE) & m, (ref == VISIBLE) & m
    fv, fi = ev & ~rv, ~ev & rv & m
    area = cell * cell / 1e6

    def f1(sel):
        tp, fp, fn = (ev & rv & sel).sum(), (fv & sel).sum(), (fi & sel).sum()
        return round(float(2 * tp / (2 * tp + fp + fn)), 4) if tp + fp + fn else None

    rk = float(rv.sum() * area)
    return dict(
        ref_visible_km2=round(rk, 3), visible_km2=round(float(ev.sum() * area), 3),
        area_bias_pct=round(100 * (ev.sum() - rv.sum()) / rv.sum(), 2) if rv.sum() else None,
        agreement_pct=round(100 * float((ev == rv)[m].mean()), 3), f1=f1(m),
        false_visible_km2=round(float(fv.sum() * area), 3), false_invisible_km2=round(float(fi.sum() * area), 3),
        f1_by_band={f"{a/1000:g}-{b/1000:g}": f1(m & (dist >= a) & (dist < b)) for a, b in BANDS if a < dist[m].max()},
        relief_m=round(float(dem[m].max() - dem[m].min()), 1),
        roughness_m=round(float(np.std(np.diff(dem, axis=1)[m[:, 1:]])), 2),
    )


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", type=int, help="only the first N sites")
    ap.add_argument("--radii", type=float, nargs="+", default=[5000, 10000])
    ap.add_argument("--engines", nargs="+", default=list(ENGINES))
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    (DATA / "grids").mkdir(parents=True, exist_ok=True)
    (DATA / "codes").mkdir(parents=True, exist_ok=True)
    res_file = OUT / "results.jsonl"
    done = set()
    if res_file.exists():
        for line in res_file.read_text(encoding="utf-8").splitlines():
            d = json.loads(line)
            done.add((d["site"], d["radius"], d["scenario"], d["engine"]))

    sites = list(SITES.items())[: a.sites] if a.sites else list(SITES.items())
    for name, (lon, lat, region) in sites:
        for radius in a.radii:
            grid = DATA / "grids" / f"{slug(name)}_{int(radius)}.tif"
            if not grid.exists():
                ds = observer_grid(lon, lat, radius, dem="glo30", resolution_m=RES, out_path=grid)
                ds = None  # close = flush to disk before anyone reads it
            from osgeo import gdal
            ds = gdal.Open(str(grid))  # keep the dataset alive while reading its band
            dem = ds.GetRasterBand(1).ReadAsArray()
            ds = None
            for scenario in ("flat", "curved"):
                curved = scenario == "curved"
                ref_file = DATA / "codes" / f"{slug(name)}_{int(radius)}_{scenario}_r3.npz"
                if ref_file.exists():
                    ref = np.load(ref_file)["codes"]
                else:
                    ref = ENGINES["r3"](grid, radius, H, curved)
                    np.savez_compressed(ref_file, codes=ref)
                for eng in a.engines:
                    if eng == "r3" or (curved and eng not in SUPPORTS_CURVATURE):
                        continue
                    key = (name, radius, scenario, eng)
                    if key in done:
                        continue
                    row = dict(site=name, region=region, lon=lon, lat=lat, radius=radius, scenario=scenario, engine=eng)
                    try:
                        t0 = time.perf_counter()
                        codes = ENGINES[eng](grid, radius, H, curved)
                        row["seconds"] = round(time.perf_counter() - t0, 3)
                        np.savez_compressed(DATA / "codes" / f"{slug(name)}_{int(radius)}_{scenario}_{eng}.npz", codes=codes)
                        row.update(metrics(codes, ref, dem, RES))
                    except Exception as e:  # an engine failing is a result, not a crash
                        row["error"] = f"{type(e).__name__}: {str(e)[:300]}"
                        traceback.print_exc()
                    with res_file.open("a", encoding="utf-8") as f:
                        f.write(json.dumps(row) + "\n")
                    print(f"{name:32s} {radius/1000:>4g} km {scenario:6s} {eng:9s} "
                          f"F1={row.get('f1')} bias={row.get('area_bias_pct')}% {row.get('error', '')[:80]}", flush=True)


if __name__ == "__main__":
    main()
