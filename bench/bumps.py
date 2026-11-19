"""Bump test: how big a shadow does each engine draw behind small obstacles on flat ground?

Shadows larger than R3's mean the engine is conservative (under-counts visible area);
smaller means generous (over-counts). Flat Earth, observer 1.5 m, 30 m cells, radius 3 km.

Cases
  bump   one raised cell, height b, ~600 m from the observer, in three directions:
         E (along a grid axis), NE (diagonal), ENE (2:1 offset, between axes)
  checker  checkerboard of +/- s
  noise    white-noise roughness s (6 seeds)

    python -m bench.bumps        (run in .env-bench: needs all engines)
Writes outputs/bumps/*.json and *.png.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from viewshed import INVISIBLE, VISIBLE
from viewshed.reference import r3

from .engines import ENGINES
from .explain import COLOR, GRID, INK, LABEL, MUTED, style, write_aeqd
from .run import ROOT

OUT = ROOT / "outputs" / "bumps"
OUT.mkdir(parents=True, exist_ok=True)
ENG = ["gdal", "whitebox", "xrspatial"]
CELL, H, R, N = 30.0, 1.5, 3000.0, 211
O = N // 2
BASE = 300.0
DIRS = {"E (grid axis)": (0, 20), "NE (diagonal)": (-14, 14), "ENE (between axes)": (-9, 18)}  # (drow, dcol)
HEIGHTS = [0.1, 0.3, 1.0, 3.0]


def analysed(shape):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.hypot(xx - O, yy - O) * CELL <= R


def run_case(dem, label, tmp):
    p = tmp / f"{label}.tif"
    write_aeqd(dem, p)
    ref = r3(dem, CELL, H, 0.0, R, 0.0)
    m = analysed(dem.shape) & np.isin(ref, (VISIBLE, INVISIBLE))
    out = {"r3": ref}
    row = {"case": label, "r3_hidden": int((ref[m] == INVISIBLE).sum()), "r3_visible": int((ref[m] == VISIBLE).sum())}
    for e in ENG:
        c = ENGINES[e](p, R, H, False)
        out[e] = c
        mm = m & np.isin(c, (VISIBLE, INVISIBLE))
        ev, rv = c == VISIBLE, ref == VISIBLE
        row[f"{e}_hidden"] = int((c[mm] == INVISIBLE).sum())
        row[f"{e}_false_hidden"] = int((rv & ~ev & mm).sum())   # engine shadows ground R3 sees
        row[f"{e}_false_visible"] = int((~rv & ev & mm).sum())  # engine sees ground R3 shadows
        row[f"{e}_visible_bias_pct"] = round(100 * ((ev & mm).sum() - (rv & mm).sum()) / max((rv & mm).sum(), 1), 3)
    return row, out


def main():
    tmp = Path(tempfile.mkdtemp())
    rows, maps = [], {}
    for dname, (dr, dc) in DIRS.items():
        for b in HEIGHTS:
            dem = np.full((N, N), BASE)
            dem[O + dr, O + dc] += b
            row, out = run_case(dem, f"bump_{dname.split()[0]}_{b}", tmp)
            row.update(kind="bump", direction=dname, height=b, distance_m=round(float(np.hypot(dr, dc) * CELL)))
            rows.append(row)
            if b == 1.0:
                maps[dname] = (out, (dr, dc))
            print(f"bump {dname:20s} b={b:>4} m: shadow cells R3 {row['r3_hidden']:4d} | "
                  + " | ".join(f"{LABEL[e]} {row[f'{e}_hidden']:4d} (+{row[f'{e}_false_hidden']}/-{row[f'{e}_false_visible']})" for e in ENG), flush=True)
    for s in (0.1, 0.5):
        yy, xx = np.mgrid[0:N, 0:N]
        dem = BASE + s * np.where((yy + xx) % 2 == 0, 1.0, -1.0)
        row, _ = run_case(dem, f"checker_{s}", tmp)
        row.update(kind="checker", height=s)
        rows.append(row)
        print(f"checker s={s}: visible-area bias % " + " ".join(f"{LABEL[e]} {row[f'{e}_visible_bias_pct']:+.2f}" for e in ENG), flush=True)
    for s in (0.1, 0.5, 2.0):
        for seed in range(6):
            dem = BASE + s * np.random.default_rng(seed).normal(size=(N, N))
            row, _ = run_case(dem, f"noise_{s}_{seed}", tmp)
            row.update(kind="noise", height=s, seed=seed)
            rows.append(row)
            print(f"noise s={s} seed={seed}: visible-area bias % " + " ".join(f"{LABEL[e]} {row[f'{e}_visible_bias_pct']:+.2f}" for e in ENG), flush=True)
    (OUT / "results.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    figures(rows, maps)


MIN_VISIBLE_CELLS = 1000  # below this, % bias is dominated by a handful of cells


def figures(rows, maps):
    if maps:
        bump_maps(maps)
    bump_ratio(rows)
    roughness_bias(rows)


def bump_maps(maps):
    # 1. shadow behind a 1 m bump, zoomed, per direction x engine
    cols = ["r3"] + ENG
    fig, axs = plt.subplots(len(maps), len(cols), figsize=(3.2 * len(cols), 3.0 * len(maps)), dpi=120)
    for i, (dname, (out, (dr, dc))) in enumerate(maps.items()):
        ref = out["r3"]
        r0, c0 = O + dr, O + dc
        rs, cs = slice(max(r0 - 45, 0), min(r0 + 45, N)), slice(max(c0 - 45, 0), min(c0 + 45, N))
        for j, k in enumerate(cols):
            ax = axs[i, j]
            c = out[k]
            img = np.ones((*ref.shape, 3))
            if k == "r3":
                img[ref == INVISIBLE] = (0.35, 0.38, 0.42)
            else:
                both = (c == INVISIBLE) & (ref == INVISIBLE)
                img[both] = (0.35, 0.38, 0.42)
                img[(c == INVISIBLE) & (ref == VISIBLE)] = (0.89, 0.29, 0.28)   # extra shadow
                img[(c == VISIBLE) & (ref == INVISIBLE)] = (0.16, 0.47, 0.84)   # missing shadow
            ax.imshow(img[rs, cs], interpolation="nearest")
            ax.plot(c0 - cs.start, r0 - rs.start, marker="^", color="black", ms=6)
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title("Exact R3" if k == "r3" else LABEL[k], fontsize=10, color=INK)
            if j == 0:
                ax.set_ylabel(dname, fontsize=9, color=MUTED)
    fig.suptitle("Shadow behind a single 1 m bump, 600 m away (grey = shadow both agree on; red = extra shadow; blue = missing shadow)",
                 fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "bump_shadows.png"); plt.close(fig)


def bump_ratio(rows):
    # 2. shadow size relative to R3 vs bump height, faceted by direction
    fig, axs = plt.subplots(1, len(DIRS), figsize=(4.0 * len(DIRS), 3.6), dpi=130, sharey=True)
    for ax, dname in zip(axs, DIRS):
        style(ax)
        rs = [r for r in rows if r["kind"] == "bump" and r["direction"] == dname]
        for e in ENG:
            ax.plot([r["height"] for r in rs], [r[f"{e}_hidden"] / max(r["r3_hidden"], 1) for r in rs],
                    marker="o", ms=5, lw=2, color=COLOR[e], label=LABEL[e])
        ax.axhline(1, color=MUTED, lw=1, ls="--")
        ax.set_xscale("log")
        ax.set_title(dname, fontsize=10, color=INK)
        ax.set_xlabel("bump height (m)", fontsize=9, color=MUTED)
    axs[0].set_ylabel("shadow area / exact shadow area", fontsize=9, color=MUTED)
    axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Above 1 = conservative (draws too much shadow); below 1 = generous", fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(OUT / "bump_shadow_ratio.png"); plt.close(fig)


def roughness_bias(rows):
    # 3. signed visible-area bias on roughness (cases with too little visible ground are excluded)
    rows = [r for r in rows if r["kind"] == "bump" or r["r3_visible"] >= MIN_VISIBLE_CELLS]
    cases = [("checker", 0.1), ("checker", 0.5), ("noise", 0.1), ("noise", 0.5), ("noise", 2.0)]
    fig, ax = plt.subplots(figsize=(8, 3.6), dpi=130)
    style(ax)
    w = 0.26
    for k, e in enumerate(ENG):
        vals = [np.mean([r[f"{e}_visible_bias_pct"] for r in rows if r["kind"] == kd and r["height"] == s]) for kd, s in cases]
        ax.bar(np.arange(len(cases)) + (k - 1) * w, vals, width=w - 0.02, color=COLOR[e], label=LABEL[e])
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xticks(range(len(cases)), [f"{kd} ±{s:g} m" for kd, s in cases], fontsize=9)
    ax.set_ylabel("visible area vs exact (%)", fontsize=9, color=MUTED)
    ax.set_title(f"On rough flat ground: which engines over- or under-count (runs with >= {MIN_VISIBLE_CELLS} visible cells)",
                 fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "roughness_bias.png"); plt.close(fig)


if __name__ == "__main__":
    main()
