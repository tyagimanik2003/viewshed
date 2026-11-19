"""Summarise observer-position sensitivity (outputs/observer/results.jsonl).

    python -m bench.observer_analyse
"""
from __future__ import annotations

import json
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from .explain import INK, MUTED, style
from .observer import OUT, SHIFTS

BLUE, GREY = "#2a78d6", "#9aa3ad"
MIN_KM2 = 1.0  # below this, percentages are dominated by a handful of cells


def main():
    rows = [json.loads(l) for l in (OUT / "results.jsonl").read_text(encoding="utf-8").splitlines()]
    rows = [r for r in rows if r["area_change_pct"] is not None]
    ok = [r for r in rows if r["base_visible_km2"] >= MIN_KM2]
    sites = sorted({r["site"] for r in rows})
    table = []
    for d in SHIFTS:
        a = np.array([r["area_change_pct"] for r in ok if r["shift_m"] == d])
        ious = [r["iou"] for r in ok if r["shift_m"] == d and r.get("iou") is not None]
        table.append(dict(shift_m=d, runs=int(a.size), median_area_change_pct=round(float(np.median(a)), 2),
                          p90_area_change_pct=round(float(np.percentile(a, 90)), 2), max_area_change_pct=round(float(a.max()), 1),
                          median_iou=round(float(np.median(ious)), 3) if ious else None,
                          p10_iou=round(float(np.percentile(ious, 10)), 3) if ious else None))
    per_site = defaultdict(dict)
    for r in ok:
        per_site[r["site"]].setdefault(r["shift_m"], []).append(r["area_change_pct"])
        per_site[r["site"]]["slope"] = r["observer_slope_deg"]
        per_site[r["site"]]["visible"] = r["base_visible_km2"]
    site_rows = [dict(site=s, slope=v["slope"], visible_km2=v["visible"],
                      median_change_10m=float(np.median(v[10])) if 10 in v else None) for s, v in per_site.items()]
    sl = np.array([r["slope"] for r in site_rows if r["median_change_10m"] is not None])
    ch = np.array([r["median_change_10m"] for r in site_rows if r["median_change_10m"] is not None])
    rho = float(np.corrcoef(np.argsort(np.argsort(sl)), np.argsort(np.argsort(ch)))[0, 1]) if sl.size > 2 else None
    summary = dict(sites=len(sites), sites_in_stats=len(per_site), table=table, spearman_slope_vs_change10=round(rho, 3) if rho else None,
                   per_site=site_rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    for t in table:
        print(t)
    print("Spearman(observer slope, median area change at 10 m) =", summary["spearman_slope_vs_change10"])

    # Figure 1: area change vs shift distance
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=130)
    style(ax)
    rng = np.random.default_rng(0)
    for i, d in enumerate(SHIFTS):
        a = np.array([r["area_change_pct"] for r in ok if r["shift_m"] == d])
        ax.scatter(i + rng.uniform(-0.18, 0.18, a.size), a, s=8, color=BLUE, alpha=0.35, edgecolor="none", zorder=2)
        ax.plot([i - 0.28, i + 0.28], [np.median(a)] * 2, color=INK, lw=2.4, zorder=3)
        ax.text(i + 0.32, np.median(a), f"{np.median(a):.1f}%", fontsize=8, color=INK, va="center")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(range(len(SHIFTS)), [f"{d} m" for d in SHIFTS])
    ax.set_xlabel("observer moved by (8 directions each)", fontsize=9, color=MUTED)
    ax.set_ylabel("change in visible area (%)", fontsize=9, color=MUTED)
    ax.set_title(f"A GPS-sized shift changes the exact viewshed ({len(per_site)} sites, 5 km, bar = median)", fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "area_change_vs_shift.png"); plt.close(fig)

    # Figure 2: sensitivity vs observer slope
    fig, ax = plt.subplots(figsize=(6.6, 4.2), dpi=130)
    style(ax)
    ax.scatter(sl, ch, s=32, color=BLUE, edgecolor="white", lw=1.1, zorder=3)
    for r in site_rows:
        if r["median_change_10m"] is not None and r["median_change_10m"] > np.percentile(ch, 85):
            ax.annotate(r["site"].split(",")[0][:22], (r["slope"], r["median_change_10m"]), xytext=(4, 3),
                        textcoords="offset points", fontsize=7.5, color=INK)
    ax.set_xlabel("ground slope where the observer stands (°)", fontsize=9, color=MUTED)
    ax.set_ylabel("median area change for a 10 m shift (%)", fontsize=9, color=MUTED)
    ax.set_title(f"Is a steep standing spot more sensitive? (Spearman ρ = {summary['spearman_slope_vs_change10']})", fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "change_vs_slope.png"); plt.close(fig)


if __name__ == "__main__":
    main()
