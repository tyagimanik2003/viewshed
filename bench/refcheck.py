"""Reference-model robustness: do the benchmark conclusions depend on R3's terrain model?

For every site x radius (flat-Earth scenario, observer 1.5 m), R3 is computed under three
terrain models (linear = the benchmark reference, nearest, flatcell). Each engine is scored
against each reference, and the references are compared with each other.

    python -m bench.refcheck            (resumable; outputs/refcheck/)
    python -m bench.refcheck --analyse  (tables + figure only)
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from statistics import median

import numpy as np
from osgeo import gdal

from viewshed import INVISIBLE, VISIBLE
from viewshed.reference import r3

from .run import DATA, ROOT, slug
from .sites import SITES

OUT = ROOT / "outputs" / "refcheck"
OUT.mkdir(parents=True, exist_ok=True)
MODELS = ["linear", "nearest", "flatcell"]
ENG = ["gdal", "whitebox", "xrspatial"]
H, CELL = 1.5, 30.0
MIN_KM2 = 1.0


def score(codes, ref):
    m = np.isin(ref, (VISIBLE, INVISIBLE)) & np.isin(codes, (VISIBLE, INVISIBLE))
    ev, rv = (codes == VISIBLE) & m, (ref == VISIBLE) & m
    tp, fp, fn = (ev & rv).sum(), (ev & ~rv).sum(), (~ev & rv).sum()
    return dict(f1=float(2 * tp / (2 * tp + fp + fn)) if tp + fp + fn else 1.0,
                bias_pct=float(100 * (ev.sum() - rv.sum()) / rv.sum()) if rv.sum() else None,
                agree_pct=float(100 * (ev == rv)[m].mean()), ref_km2=float(rv.sum() * CELL * CELL / 1e6))


def run():
    res = OUT / "results.jsonl"
    done = {(d["site"], d["radius"]) for d in map(json.loads, res.read_text(encoding="utf-8").splitlines())} if res.exists() else set()
    for name, (lon, lat, region) in SITES.items():
        for radius in (5000, 10000):
            if (name, radius) in done:
                continue
            ds = gdal.Open(str(DATA / "grids" / f"{slug(name)}_{radius}.tif"))
            dem = ds.GetRasterBand(1).ReadAsArray()
            ds = None
            refs = {"linear": np.load(DATA / "codes" / f"{slug(name)}_{radius}_flat_r3.npz")["codes"]}
            for mdl in ("nearest", "flatcell"):
                f = DATA / "codes" / f"{slug(name)}_{radius}_flat_r3{mdl}.npz"
                if f.exists():
                    refs[mdl] = np.load(f)["codes"]
                else:
                    refs[mdl] = r3(dem, CELL, H, 0.0, radius, 0.0, model=mdl)
                    np.savez_compressed(f, codes=refs[mdl])
            row = dict(site=name, region=region, radius=radius)
            for a in MODELS:
                for b in MODELS:
                    if a < b:
                        row[f"ref_{a}_vs_{b}"] = score(refs[b], refs[a])  # how B scores against A
            for e in ENG:
                codes = np.load(DATA / "codes" / f"{slug(name)}_{radius}_flat_{e}.npz")["codes"]
                row[e] = {mdl: score(codes, refs[mdl]) for mdl in MODELS}
            with res.open("a", encoding="utf-8") as f:
                f.write(json.dumps(row) + "\n")
            print(f"{name[:34]:34s} {radius//1000:>2} km  refs km2: " +
                  " ".join(f"{mdl} {row['gdal'][mdl]['ref_km2']:.2f}" for mdl in MODELS) +
                  "  | GDAL F1: " + " ".join(f"{row['gdal'][mdl]['f1']:.3f}" for mdl in MODELS), flush=True)


def analyse():
    rows = [json.loads(l) for l in (OUT / "results.jsonl").read_text(encoding="utf-8").splitlines()]
    ok = [r for r in rows if r["gdal"]["linear"]["ref_km2"] >= MIN_KM2]
    plains = lambda r: r["region"] == "plains"
    out = {"runs": len(rows), "runs_in_stats": len(ok)}
    lines = ["| Engine | Reference model | All: median F1 | All: median bias | Plains: median F1 | Plains: median bias | Mountains*: median F1 |",
             "|---|---|---|---|---|---|---|"]
    for e in ENG:
        for mdl in MODELS:
            def stat(sel, key):
                v = [r[e][mdl][key] for r in ok if sel(r) and r[e][mdl][key] is not None]
                return median(v) if v else float("nan")
            mount = lambda r: r["region"] in ("high Himalaya", "high hills", "Ladakh plateau")
            s = dict(all_f1=stat(lambda r: True, "f1"), all_bias=stat(lambda r: True, "bias_pct"),
                     plains_f1=stat(plains, "f1"), plains_bias=stat(plains, "bias_pct"), mount_f1=stat(mount, "f1"))
            out[f"{e}|{mdl}"] = s
            lines.append(f"| {e} | {mdl} | {s['all_f1']:.3f} | {s['all_bias']:+.1f}% | {s['plains_f1']:.3f} | "
                         f"{s['plains_bias']:+.1f}% | {s['mount_f1']:.3f} |")
    lines += ["", "| Reference vs reference | All: median F1 | All: median bias | Plains: median F1 | Plains: median bias |",
              "|---|---|---|---|---|"]
    for a in MODELS:
        for b in MODELS:
            if a < b:
                k = f"ref_{a}_vs_{b}"
                v = lambda sel, key: median([r[k][key] for r in ok if sel(r) and r[k][key] is not None])
                s = dict(all_f1=v(lambda r: True, "f1"), all_bias=v(lambda r: True, "bias_pct"),
                         plains_f1=v(plains, "f1"), plains_bias=v(plains, "bias_pct"))
                out[k] = s
                lines.append(f"| {b} scored against {a} | {s['all_f1']:.3f} | {s['all_bias']:+.1f}% | {s['plains_f1']:.3f} | {s['plains_bias']:+.1f}% |")
    # Does each engine's bias keep its sign across references? Share of runs per engine/model with bias < 0.
    for e in ENG:
        for mdl in MODELS:
            v = [r[e][mdl]["bias_pct"] for r in ok if r[e][mdl]["bias_pct"] is not None]
            out[f"{e}|{mdl}|share_negative_bias"] = float(np.mean([x < 0 for x in v]))
    # Is the "plains are hardest" ranking stable? Spearman between reference models of per-run engine F1.
    (OUT / "summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    figure(ok)


def figure(ok):
    import matplotlib.pyplot as plt

    from .explain import COLOR, INK, LABEL, MUTED, style
    bench = [json.loads(l) for l in (ROOT / "outputs" / "bench" / "results.jsonl").read_text(encoding="utf-8").splitlines()]
    relief = {(b["site"], int(b["radius"])): b["relief_m"] for b in bench if b["engine"] == "gdal" and b["scenario"] == "flat"}
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), dpi=130, sharey=True)
    for ax, mdl in zip(axs, MODELS):
        style(ax)
        for e in ENG:
            ax.scatter([relief[(r["site"], r["radius"])] for r in ok], [r[e][mdl]["f1"] for r in ok], s=14, color=COLOR[e],
                       edgecolor="white", lw=0.6, label=LABEL[e], alpha=0.85)
        ax.set_xscale("log")
        ax.set_title(f"reference = {mdl}", fontsize=10, color=INK)
        ax.set_xlabel("relief within radius (m)", fontsize=9, color=MUTED)
    axs[0].set_ylabel("engine F1 vs reference", fontsize=9, color=MUTED)
    axs[0].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Do the conclusions depend on the reference's terrain model?", fontsize=11, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "f1_by_reference_model.png"); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    if not a.analyse:
        run()
    analyse()
