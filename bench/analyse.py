"""Summarise ViewshedBench results: tables + figures in outputs/bench/.

    python -m bench.analyse
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
import numpy as np

from viewshed import INVISIBLE, VISIBLE

from .engines import ENGINE_NOTES
from .run import DATA, OUT, slug

ENGINE_LABEL = {"gdal": "GDAL", "grass": "GRASS r.viewshed", "whitebox": "WhiteboxTools", "xrspatial": "xarray-spatial"}
ORDER = ["gdal", "grass", "whitebox", "xrspatial"]
RADIUS_COLOR = {5000.0: "#2a78d6", 10000.0: "#eb6834"}  # categorical slots 1-2 (validated palette)
INK, MUTED, GRID = "#1f2328", "#5b6470", "#e4e7ea"
MIN_VISIBLE_KM2 = 1.0


def load():
    rows = [json.loads(l) for l in (OUT / "results.jsonl").read_text(encoding="utf-8").splitlines()]
    return [r for r in rows if "error" not in r], [r for r in rows if "error" in r]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)


def summary(rows):
    by = defaultdict(list)
    for r in rows:
        by[(r["scenario"], r["engine"])].append(r)
    table = []
    for (sc, eng), rs in sorted(by.items(), key=lambda kv: (kv[0][0], ORDER.index(kv[0][1]) if kv[0][1] in ORDER else 9)):
        # Relative metrics are unstable when almost nothing is visible (e.g. Ludhiana: 0.17 km2),
        # so percentage statistics use only runs with >= MIN_VISIBLE_KM2 visible in the reference.
        big = [r for r in rs if r["ref_visible_km2"] >= MIN_VISIBLE_KM2] or rs
        table.append(dict(
            scenario=sc, engine=eng, runs=len(rs), runs_in_pct_stats=len(big),
            median_agreement_pct=round(median(r["agreement_pct"] for r in rs), 2),
            worst_agreement_pct=round(min(r["agreement_pct"] for r in rs), 2),
            median_abs_diff_km2=round(median(abs(r["visible_km2"] - r["ref_visible_km2"]) for r in rs), 3),
            worst_abs_diff_km2=round(max(abs(r["visible_km2"] - r["ref_visible_km2"]) for r in rs), 3),
            median_f1=round(median(r["f1"] for r in big), 4), worst_f1=round(min(r["f1"] for r in big), 4),
            median_bias_pct=round(median(r["area_bias_pct"] for r in big), 2),
            worst_abs_bias_pct=round(max(abs(r["area_bias_pct"]) for r in big), 2),
            median_seconds=round(median(r["seconds"] for r in rs), 3),
        ))
    return table


def fig_f1_vs_relief(rows, scenario, path):
    engs = [e for e in ORDER if any(r["engine"] == e and r["scenario"] == scenario for r in rows)]
    fig, axs = plt.subplots(1, len(engs), figsize=(3.6 * len(engs), 3.8), dpi=130, sharey=True, squeeze=False)
    ymin = min(r["f1"] for r in rows if r["scenario"] == scenario)
    for ax, eng in zip(axs[0], engs):
        style(ax)
        for rad, col in RADIUS_COLOR.items():
            rs = [r for r in rows if r["engine"] == eng and r["scenario"] == scenario and r["radius"] == rad]
            ax.scatter([r["relief_m"] for r in rs], [r["f1"] for r in rs], s=34, color=col,
                       edgecolor="white", linewidth=1.2, label=f"{rad/1000:g} km radius", zorder=3)
        ax.set_title(ENGINE_LABEL.get(eng, eng), fontsize=10, color=INK)
        ax.set_xlabel("relief within radius (m)", fontsize=9, color=MUTED)
        ax.set_ylim(max(0, ymin - 0.02), 1.002)
    axs[0][0].set_ylabel("F1 of visible cells vs exact R3", fontsize=9, color=MUTED)
    axs[0][0].legend(fontsize=8, frameon=False, loc="lower left")
    fig.suptitle(f"Agreement with the exact reference ({scenario}-earth scenario)", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_bias(rows, path):
    combos = [(sc, e) for sc in ("flat", "curved") for e in ORDER if any(r["engine"] == e and r["scenario"] == sc for r in rows)]
    fig, ax = plt.subplots(figsize=(8, 0.55 * len(combos) + 1.4), dpi=130)
    style(ax)
    ax.axvline(0, color=MUTED, lw=1)
    for i, (sc, eng) in enumerate(combos):
        rs = [r for r in rows if r["engine"] == eng and r["scenario"] == sc]
        for r in rs:
            ax.scatter(r["area_bias_pct"], i, s=26, color=RADIUS_COLOR[r["radius"]], edgecolor="white", linewidth=1, zorder=3)
        m = median(r["area_bias_pct"] for r in rs)
        ax.plot([m, m], [i - 0.3, i + 0.3], color=INK, lw=2, zorder=4)
    ax.set_yticks(range(len(combos)), [f"{ENGINE_LABEL[e]} · {sc}" for sc, e in combos], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("visible area vs exact R3 (%)   ← under-counts · over-counts →", fontsize=9, color=MUTED)
    ax.set_title("How much each engine mis-states the visible area (bar = median)", fontsize=11, color=INK)
    for rad, col in RADIUS_COLOR.items():
        ax.scatter([], [], color=col, s=26, label=f"{rad/1000:g} km")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_maps(site, radius, scenario, path):
    from osgeo import gdal

    from viewshed.plot import hillshade

    grid = DATA / "grids" / f"{slug(site)}_{int(radius)}.tif"
    ds = gdal.Open(str(grid))
    dem = ds.GetRasterBand(1).ReadAsArray()
    ds = None
    ref = np.load(DATA / "codes" / f"{slug(site)}_{int(radius)}_{scenario}_r3.npz")["codes"]
    engs = [e for e in ORDER if (DATA / "codes" / f"{slug(site)}_{int(radius)}_{scenario}_{e}.npz").exists()]
    fig, axs = plt.subplots(1, len(engs) + 1, figsize=(3.6 * (len(engs) + 1), 3.9), dpi=120)
    hs = hillshade(dem, 30.0)
    h = radius / 1000
    ext = [-h, h, -h, h]
    for ax, (label, codes) in zip(axs, [("Exact R3", ref)] + [(ENGINE_LABEL[e], np.load(DATA / "codes" / f"{slug(site)}_{int(radius)}_{scenario}_{e}.npz")["codes"]) for e in engs]):
        ax.imshow(hs, cmap="gray", vmin=0, vmax=1, extent=ext)
        rgba = np.zeros((*codes.shape, 4))
        if label == "Exact R3":
            rgba[codes == VISIBLE] = (0.98, 0.62, 0.11, 0.65)
        else:
            m = np.isin(ref, (VISIBLE, INVISIBLE)) & np.isin(codes, (VISIBLE, INVISIBLE))
            rgba[m & (codes == VISIBLE) & (ref == VISIBLE)] = (0.98, 0.62, 0.11, 0.35)
            rgba[m & (codes == VISIBLE) & (ref != VISIBLE)] = (0.16, 0.47, 0.84, 0.95)
            rgba[m & (codes != VISIBLE) & (ref == VISIBLE)] = (0.89, 0.29, 0.28, 0.95)
        ax.imshow(rgba, extent=ext, interpolation="nearest")
        ax.plot(0, 0, marker="*", ms=11, mec="white", mfc="black")
        ax.set_title(label, fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{site}, {radius/1000:g} km, {scenario}-earth. Blue = engine says visible but hidden; red = engine says hidden but visible",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    rows, errors = load()
    table = summary(rows)
    (OUT / "summary.json").write_text(json.dumps(dict(summary=table, errors=errors, engines=ENGINE_NOTES), indent=2), encoding="utf-8")
    lines = ["| Scenario | Engine | Runs | Cells agreeing (median / worst) | abs area diff km2 (median / worst) | "
             "F1 (median / worst)* | Area bias (median / worst abs)* | Median time |",
             "|---|---|---|---|---|---|---|---|"]
    lines += [f"| {t['scenario']} | {ENGINE_LABEL.get(t['engine'], t['engine'])} | {t['runs']} | "
              f"{t['median_agreement_pct']:.2f}% / {t['worst_agreement_pct']:.2f}% | "
              f"{t['median_abs_diff_km2']:.2f} / {t['worst_abs_diff_km2']:.2f} | {t['median_f1']:.3f} / {t['worst_f1']:.3f} | "
              f"{t['median_bias_pct']:+.1f}% / {t['worst_abs_bias_pct']:.1f}% | {t['median_seconds']:.2f} s |" for t in table]
    lines += ["", f"*Runs with at least {MIN_VISIBLE_KM2:g} km2 visible in the reference only."]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    if errors:
        print(f"\n{len(errors)} failed runs, e.g.: {errors[0]['engine']}: {errors[0]['error'][:200]}")
    for sc in ("flat", "curved"):
        if any(r["scenario"] == sc for r in rows):
            fig_f1_vs_relief(rows, sc, OUT / f"f1_vs_relief_{sc}.png")
    fig_bias(rows, OUT / "area_bias.png")
    fig_maps("Leh, Ladakh", 5000.0, "flat", OUT / "maps_leh_5km_flat.png")
    region_table(rows)


def region_table(rows):
    """Median F1 / area bias per landscape (flat-Earth scenario, runs with >= 1 km2 visible)."""
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r["scenario"] != "flat" or r["ref_visible_km2"] < MIN_VISIBLE_KM2:
            continue
        reg = r["region"].replace("Ladakh plateau", "Ladakh plateau")
        by[reg][r["engine"]].append(r)
    order = ["high Himalaya", "high hills", "Ladakh plateau", "mountains", "foothills", "plateau", "plains"]
    lines = ["| Landscape | Runs | " + " | ".join(f"{ENGINE_LABEL[e]} F1 / bias" for e in ORDER if e != "grass") + " |",
             "|---|---|" + "---|" * (len(ORDER) - 1)]
    out = {}
    for reg in sorted(by, key=lambda k: order.index(k) if k in order else 99):
        cells, out[reg] = [], {}
        for e in ORDER:
            if e == "grass":
                continue
            rs = by[reg].get(e, [])
            if rs:
                f1, b = median(x["f1"] for x in rs), median(x["area_bias_pct"] for x in rs)
                out[reg][e] = dict(runs=len(rs), median_f1=round(f1, 4), median_bias_pct=round(b, 2))
                cells.append(f"{f1:.3f} / {b:+.1f}%")
            else:
                cells.append("–")
        n = max(len(v) for v in by[reg].values())
        lines.append(f"| {reg} | {n} | " + " | ".join(cells) + " |")
    (OUT / "by_region.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "by_region.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
