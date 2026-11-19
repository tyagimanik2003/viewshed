"""Why do engines split on flat land? Test the grazing-sightline hypothesis.

Hypothesis: engines disagree where the exact sightline clearance |c| is tiny (a few cm),
and flat terrain has far more such cells. On open flat ground c ~ h * cell / (D - cell).

Part A  real sites: disagreement rate vs |clearance|, and per-site near-tie fraction vs spread.
Part B  synthetic terrain sweep: large-scale relief A x micro-roughness s x observer height h.
Part C  each engine's effective threshold: the clearance offset d (cm) for which
        "visible iff clearance > d" best reproduces that engine.

    python -m bench.explain          (run in .env-bench: needs all engines)
Writes outputs/explain/*.json and *.png.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal, osr

from viewshed import INVISIBLE, VISIBLE
from viewshed.reference import r3

from .engines import ENGINES
from .run import DATA, ROOT, slug
from .sites import SITES

gdal.UseExceptions()
OUT = ROOT / "outputs" / "explain"
OUT.mkdir(parents=True, exist_ok=True)
ENG = ["gdal", "whitebox", "xrspatial"]
LABEL = {"gdal": "GDAL", "whitebox": "WhiteboxTools", "xrspatial": "xarray-spatial"}
COLOR = {"gdal": "#2a78d6", "whitebox": "#eb6834", "xrspatial": "#1baf7a"}  # validated slots 1-3
INK, MUTED, GRID = "#1f2328", "#5b6470", "#e4e7ea"
BINS = [0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 1e9]
CELL, H = 30.0, 1.5


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)


def read(path):
    ds = gdal.Open(str(path))
    a = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    return a


def threshold_offset(clr, eng_visible, mask):
    """Best d (m): engine ~ 'visible iff clearance > d'. Positive = stricter than exact."""
    ds_ = np.linspace(-0.5, 0.5, 201)
    c, v = clr[mask], eng_visible[mask]
    near = np.abs(c) < 1.0  # only cells where the offset can matter
    c, v = c[near], v[near]
    if c.size == 0:
        return None, None
    acc = [((c > d) == v).mean() for d in ds_]
    i = int(np.argmax(acc))
    return float(ds_[i]), float(acc[i])


# ----------------------------------------------------------------------------- Part A
def part_a():
    pooled = {e: np.zeros((2, len(BINS) - 1)) for e in ENG}  # [disagree, total] per bin
    per_site, offsets = [], {e: [] for e in ENG}
    for name in SITES:
        for radius in (5000, 10000):
            grid = DATA / "grids" / f"{slug(name)}_{radius}.tif"
            dem = read(grid)
            ref, clr = r3(dem, CELL, H, 0.0, radius, 0.0, return_clearance=True)
            m = np.isin(ref, (VISIBLE, INVISIBLE)) & np.isfinite(clr)
            row = dict(site=name, radius=radius,
                       near_tie_10cm=float((np.abs(clr[m]) < 0.10).mean()),
                       near_tie_50cm=float((np.abs(clr[m]) < 0.50).mean()),
                       r3_visible_km2=float((ref[m] == VISIBLE).sum() * CELL * CELL / 1e6))
            areas = []
            for e in ENG:
                codes = np.load(DATA / "codes" / f"{slug(name)}_{radius}_flat_{e}.npz")["codes"]
                mm = m & np.isin(codes, (VISIBLE, INVISIBLE))
                dis = (codes == VISIBLE) != (ref == VISIBLE)
                idx = np.digitize(np.abs(clr[mm]), BINS) - 1
                np.add.at(pooled[e][0], idx, dis[mm])
                np.add.at(pooled[e][1], idx, 1)
                row[f"disagree_pct_{e}"] = float(100 * dis[mm].mean())
                areas.append(float((codes[mm] == VISIBLE).sum() * CELL * CELL / 1e6))
                d, acc = threshold_offset(clr, codes == VISIBLE, mm)
                offsets[e].append(d)
                row[f"offset_cm_{e}"] = None if d is None else round(100 * d, 1)
            row["engine_spread_pct"] = float(100 * (max(areas) - min(areas)) / row["r3_visible_km2"]) if row["r3_visible_km2"] else None
            per_site.append(row)
            print(f"A {name[:28]:28s} {radius/1000:>3.0f} km near-tie<10cm {100*row['near_tie_10cm']:5.1f}%  "
                  f"spread {row['engine_spread_pct']:.1f}%  offsets cm " + " ".join(str(row[f'offset_cm_{e}']) for e in ENG), flush=True)
    rates = {e: (pooled[e][0] / np.maximum(pooled[e][1], 1)).tolist() for e in ENG}
    counts = pooled[ENG[0]][1].tolist()
    res = dict(bins=BINS, disagreement_rate=rates, cells_per_bin=counts, per_site=per_site,
               median_offset_cm={e: round(100 * float(np.median([o for o in offsets[e] if o is not None])), 1) for e in ENG})
    (OUT / "part_a.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    return res


# ----------------------------------------------------------------------------- Part B
def synthetic_dem(A, s, seed, n=339):
    rng = np.random.default_rng(seed)
    y, x = (np.mgrid[0:n, 0:n] - n // 2) * CELL
    field = np.zeros((n, n))
    for _ in range(30):
        cx, cy = rng.uniform(-6000, 6000, 2)
        w = rng.uniform(300, 1500)
        field += rng.normal() * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * w * w))
    field = (field - field.mean()) / (field.std() + 1e-12)
    return 300.0 + A * field + s * rng.normal(size=(n, n))


def write_aeqd(arr, path):
    n = arr.shape[0]
    srs = osr.SpatialReference()
    srs.ImportFromProj4("+proj=aeqd +lat_0=30 +lon_0=76 +datum=WGS84 +units=m +no_defs")
    ds = gdal.GetDriverByName("GTiff").Create(str(path), n, n, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((-n / 2 * CELL, CELL, 0, n / 2 * CELL, 0, -CELL))
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(arr.astype(np.float32))
    ds = None


def part_b():
    rows = []
    tmp = Path(tempfile.mkdtemp())
    for h in (1.5, 10.0):
        for A in (0.0, 2.0, 10.0, 50.0, 250.0):
            for s in (0.0, 0.1, 0.5, 2.0):
                for seed in (0, 1):
                    dem = synthetic_dem(A, s, seed)
                    p = tmp / f"syn_{h}_{A}_{s}_{seed}.tif"
                    write_aeqd(dem, p)
                    ref, clr = r3(dem, CELL, h, 0.0, 5000.0, 0.0, return_clearance=True)
                    m = np.isin(ref, (VISIBLE, INVISIBLE)) & np.isfinite(clr)
                    rv = (ref == VISIBLE) & m
                    row = dict(h=h, A=A, s=s, seed=seed, near_tie_10cm=float((np.abs(clr[m]) < 0.10).mean()),
                               r3_visible_frac=float(rv.sum() / m.sum()))
                    areas = []
                    for e in ENG:
                        c = ENGINES[e](p, 5000.0, h, False)
                        ev = (c == VISIBLE) & m
                        tp, fp, fn = (ev & rv).sum(), (ev & ~rv).sum(), (~ev & rv).sum()
                        row[f"f1_{e}"] = float(2 * tp / (2 * tp + fp + fn)) if tp + fp + fn else 1.0
                        row[f"disagree_pct_{e}"] = float(100 * (ev != rv)[m].mean())
                        areas.append(ev.sum())
                    row["engine_spread_pct"] = float(100 * (max(areas) - min(areas)) / max(rv.sum(), 1))
                    rows.append(row)
                    print(f"B h={h:>4} A={A:>5} s={s:>3} seed={seed} near-tie {100*row['near_tie_10cm']:5.1f}% "
                          f"spread {row['engine_spread_pct']:6.1f}% disagree " +
                          " ".join(f"{row[f'disagree_pct_{e}']:.2f}" for e in ENG), flush=True)
    (OUT / "part_b.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return rows


# ----------------------------------------------------------------------------- Figures
def figures(a, b):
    # 1. disagreement rate vs |clearance| (real sites, pooled)
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=130)
    style(ax)
    mids = [np.sqrt(max(lo, 0.01) * min(hi, 20)) for lo, hi in zip(BINS[:-1], BINS[1:])]
    for e in ENG:
        ax.plot(mids, [100 * r for r in a["disagreement_rate"][e]], marker="o", ms=5, lw=2, color=COLOR[e], label=LABEL[e])
    ax.set_xscale("log")
    ax.set_xlabel("|sightline clearance| from the exact reference (m), log scale", fontsize=9, color=MUTED)
    ax.set_ylabel("cells where engine disagrees (%)", fontsize=9, color=MUTED)
    ax.set_title("Disagreement rises steeply as the sightline grazes the ground", fontsize=11, color=INK)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "disagreement_vs_clearance.png"); plt.close(fig)

    # (A near-tie-fraction scatter was dropped: Part B showed near-ties alone do not cause disagreement.)

    # 3. synthetic: spread as a function of relief A and roughness s (h = 1.5)
    As, Ss = sorted({r["A"] for r in b}), sorted({r["s"] for r in b})
    grid = np.zeros((len(Ss), len(As)))
    for i, s in enumerate(Ss):
        for j, A in enumerate(As):
            v = [r["engine_spread_pct"] for r in b if r["h"] == 1.5 and r["A"] == A and r["s"] == s]
            grid[i, j] = np.mean(v)
    fig, ax = plt.subplots(figsize=(6.4, 3.9), dpi=130)
    im = ax.imshow(grid, cmap="Blues", origin="lower", aspect="auto")
    for i in range(len(Ss)):
        for j in range(len(As)):
            ax.text(j, i, f"{grid[i, j]:.0f}%", ha="center", va="center", fontsize=9,
                    color="white" if grid[i, j] > grid.max() * 0.55 else INK)
    ax.set_xticks(range(len(As)), [f"{a:g}" for a in As]); ax.set_yticks(range(len(Ss)), [f"{s:g}" for s in Ss])
    ax.set_xlabel("large-scale relief, std (m)", fontsize=9, color=MUTED)
    ax.set_ylabel("micro-roughness, std (m)", fontsize=9, color=MUTED)
    ax.set_title("Engine spread on synthetic terrain (observer 1.5 m)", fontsize=11, color=INK)
    fig.colorbar(im, ax=ax, label="spread (% of exact visible area)")
    fig.tight_layout(); fig.savefig(OUT / "synthetic_spread_heatmap.png"); plt.close(fig)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-a", action="store_true", help="reuse outputs/explain/part_a.json")
    ap.add_argument("--only-a", action="store_true", help="rerun Part A + site terrain + prediction, reuse Part B")
    args = ap.parse_args()
    if args.only_a:
        a = part_a()
        real_site_terrain()
        rows = predicted_vs_observed()
        figures(a, json.loads((OUT / "part_b.json").read_text()))
        big = [r for r in rows if r["r3_visible_km2"] >= 3]
        lo = lambda v: np.log10(max(v, 0.1))
        print("log-log r (>=3 km2):", round(float(np.corrcoef([lo(r["pred_spread_pct"]) for r in big],
                                                               [lo(r["spread_pct"]) for r in big])[0, 1]), 2), "n =", len(big))
        return
    a = json.loads((OUT / "part_a.json").read_text()) if args.skip_a else part_a()
    b = part_b()
    figures(a, b)
    print("median effective offsets (cm, + = stricter than exact):", a["median_offset_cm"])



def predicted_vs_observed():
    """Place each real site on the synthetic map (nearest relief/roughness level in log space,
    observer 1.5 m) and compare the predicted engine spread with the observed one."""
    b = json.loads((OUT / "part_b.json").read_text())
    real = json.loads((OUT / "real_site_terrain.json").read_text())
    As, Ss = sorted({r["A"] for r in b}), sorted({r["s"] for r in b})
    G = np.array([[np.mean([x["engine_spread_pct"] for x in b if x["h"] == 1.5 and x["A"] == A and x["s"] == s])
                   for A in As] for s in Ss])
    la, ls = np.log10(np.array(As) + 1), np.log10(np.array(Ss) + 1)

    def interp(A, s):  # bilinear in log10(value + 1) between the synthetic levels
        x = np.clip(np.log10(A + 1), la[0], la[-1]); y = np.clip(np.log10(s + 1), ls[0], ls[-1])
        j = min(np.searchsorted(la, x, side="right") - 1, len(la) - 2)
        i = min(np.searchsorted(ls, y, side="right") - 1, len(ls) - 2)
        tx, ty = (x - la[j]) / (la[j + 1] - la[j]), (y - ls[i]) / (ls[i + 1] - ls[i])
        return float((1 - tx) * (1 - ty) * G[i, j] + tx * (1 - ty) * G[i, j + 1]
                     + (1 - tx) * ty * G[i + 1, j] + tx * ty * G[i + 1, j + 1])

    rows = [{**r, "pred_spread_pct": round(interp(r["large_relief_std_m"], r["roughness_m"]), 1)} for r in real]
    (OUT / "predicted_vs_observed.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    fig, ax = plt.subplots(figsize=(6.8, 5.2), dpi=130)
    style(ax)
    lim = [0.05, 100]
    ax.plot(lim, lim, color=MUTED, lw=1, ls="--", zorder=1)
    ax.text(30, 50, "1 : 1", fontsize=8, color=MUTED, rotation=38)
    for r in rows:
        small = r["r3_visible_km2"] < 3
        ax.scatter(max(r["pred_spread_pct"], 0.06), max(r["spread_pct"], 0.06), s=38, zorder=3,
                   color="#9aa3ad" if small else "#2a78d6", edgecolor="white", lw=1.2)
        if r["r3_visible_km2"] >= 3 and (r["site"].startswith(("Leh", "Kurukshetra")) or r["spread_pct"] > 30):
            ax.annotate(f"{r['site'].split(',')[0].split(' (')[0]} {r['radius']//1000} km", (max(r["pred_spread_pct"], 0.06), r["spread_pct"]),
                        xytext=(5, 3), textcoords="offset points", fontsize=8, color=INK)
    ax.scatter([], [], color="#2a78d6", s=38, label="real site")
    ax.scatter([], [], color="#9aa3ad", s=38, label="real site, < 3 km² visible (noisy %)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("spread predicted from synthetic terrain with the same relief & roughness (%)", fontsize=9, color=MUTED)
    big = [r for r in rows if r["r3_visible_km2"] >= 3]
    lo = lambda v: np.log10(max(v, 0.1))
    rr = float(np.corrcoef([lo(r["pred_spread_pct"]) for r in big], [lo(r["spread_pct"]) for r in big])[0, 1])
    ax.text(0.98, 0.04, f"r = {rr:.2f} (log-log, {len(big)} runs with >= 3 km² visible)", transform=ax.transAxes,
            ha="right", fontsize=8, color=MUTED)
    ax.set_ylabel("spread observed at the real site (%)", fontsize=9, color=MUTED)
    ax.set_title("Relief + roughness predict where engines will disagree", fontsize=11, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "predicted_vs_observed.png"); plt.close(fig)
    return rows


def _box(a, k):
    p = k // 2
    c = np.pad(np.pad(a, p, mode="edge").cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / (k * k)


def real_site_terrain():
    """Large-scale relief (std of 15-cell box-smoothed DEM) and cell-scale roughness
    (std of DEM minus 3x3 mean, white-noise equivalent) within each site's radius."""
    a = json.loads((OUT / "part_a.json").read_text())
    sp = {(r["site"], r["radius"]): r for r in a["per_site"]}
    out = []
    for name in SITES:
        for radius in (5000, 10000):
            d = read(DATA / "grids" / f"{slug(name)}_{radius}.tif").astype(float)
            n = d.shape[0]
            yy, xx = np.mgrid[0:n, 0:n]
            m = np.hypot(xx - n // 2, yy - n // 2) * CELL <= radius
            r = sp[(name, radius)]
            out.append(dict(site=name, radius=radius, large_relief_std_m=round(float(_box(d, 15)[m].std()), 1),
                            roughness_m=round(float((d - _box(d, 3))[m].std() / 0.943), 2),
                            spread_pct=round(r["engine_spread_pct"], 1) if r["engine_spread_pct"] is not None else None,
                            r3_visible_km2=round(r["r3_visible_km2"], 2)))
    (OUT / "real_site_terrain.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


if __name__ == "__main__":
    main()
