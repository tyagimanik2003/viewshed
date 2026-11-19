"""Do the GDAL-dependent numbers in the synthetic and bump experiments change between GDAL 3.12.3
(used there, from .env-bench) and 3.13.3 (used for the main benchmark, .env)?

Rebuilds every seeded terrain exactly, runs the GDAL engine with the interpreter's GDAL, and
compares its metrics with the stored 3.12.3 results. Run with .env (GDAL 3.13.3):

    python -m bench.gdal_version_check     -> outputs/gdal_version_check.json
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
from osgeo import gdal

from viewshed import INVISIBLE, VISIBLE
from viewshed.reference import r3

from . import bumps as B
from . import explain as X
from .engines import ENGINES
from .run import ROOT


def metrics(c, ref, m):
    ev, rv = (c == VISIBLE) & m, (ref == VISIBLE) & m
    tp, fp, fn = (ev & rv).sum(), (ev & ~rv).sum(), (~ev & rv).sum()
    return dict(f1=float(2 * tp / (2 * tp + fp + fn)) if tp + fp + fn else 1.0,
                disagree_pct=float(100 * (ev != rv)[m].mean()),
                hidden=int(((c == INVISIBLE) & m).sum()), visible=int(ev.sum()),
                false_hidden=int((rv & ~ev & m).sum()), false_visible=int((~rv & ev & m).sum()))


def main():
    tmp = Path(tempfile.mkdtemp())
    out = {"gdal_version": gdal.__version__, "synthetic": [], "bumps": []}
    stored_b = json.loads((X.OUT / "part_b.json").read_text())
    for row in stored_b:
        h, A, s, seed = row["h"], row["A"], row["s"], row["seed"]
        dem = X.synthetic_dem(A, s, seed)
        p = tmp / f"syn_{h}_{A}_{s}_{seed}.tif"
        X.write_aeqd(dem, p)
        ref, clr = r3(dem, X.CELL, h, 0.0, 5000.0, 0.0, return_clearance=True)
        m = np.isin(ref, (VISIBLE, INVISIBLE)) & np.isfinite(clr)  # same mask as explain.part_b
        new = metrics(ENGINES["gdal"](p, 5000.0, h, False), ref, m)
        out["synthetic"].append(dict(h=h, A=A, s=s, seed=seed, f1_312=row["f1_gdal"], f1_313=new["f1"],
                                     disagree_312=row["disagree_pct_gdal"], disagree_313=new["disagree_pct"]))
    stored_bump = json.loads((B.OUT / "results.json").read_text())
    for row in stored_bump:
        if row["kind"] == "bump":
            dr, dc = B.DIRS[row["direction"]]
            dem = np.full((B.N, B.N), B.BASE)
            dem[B.O + dr, B.O + dc] += row["height"]
        elif row["kind"] == "checker":
            yy, xx = np.mgrid[0:B.N, 0:B.N]
            dem = B.BASE + row["height"] * np.where((yy + xx) % 2 == 0, 1.0, -1.0)
        else:
            dem = B.BASE + row["height"] * np.random.default_rng(row["seed"]).normal(size=(B.N, B.N))
        p = tmp / f"{row['case']}.tif"
        X.write_aeqd(dem, p)
        ref = r3(dem, B.CELL, B.H, 0.0, B.R, 0.0)
        m = B.analysed(dem.shape) & np.isin(ref, (VISIBLE, INVISIBLE))
        new = metrics(ENGINES["gdal"](p, B.R, B.H, False), ref, m)
        out["bumps"].append(dict(case=row["case"], hidden_312=row["gdal_hidden"], hidden_313=new["hidden"],
                                 false_hidden_312=row["gdal_false_hidden"], false_hidden_313=new["false_hidden"],
                                 false_visible_312=row["gdal_false_visible"], false_visible_313=new["false_visible"]))
    syn_diff = max(abs(r["f1_312"] - r["f1_313"]) for r in out["synthetic"])
    bump_diff = sum(r["hidden_312"] != r["hidden_313"] or r["false_hidden_312"] != r["false_hidden_313"]
                    or r["false_visible_312"] != r["false_visible_313"] for r in out["bumps"])
    out["summary"] = dict(synthetic_cases=len(out["synthetic"]), max_abs_f1_difference=syn_diff,
                          bump_cases=len(out["bumps"]), bump_cases_with_any_difference=int(bump_diff))
    (ROOT / "outputs" / "gdal_version_check.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out["summary"], indent=2), "GDAL", gdal.__version__)


if __name__ == "__main__":
    main()
