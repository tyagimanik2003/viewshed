"""Observer-position sensitivity: how much does a small GPS-sized shift change the viewshed?

For each site, the exact R3 viewshed (5 km, observer 1.5 m, curved Earth) is recomputed with
the observer moved d metres in 8 directions. End-to-end on a 30 m grid, as a user would
experience it: the terrain is re-gridded around the new position (bilinear), so the shift
changes both the eye position and the interpolated ground under it.

  area change  |visible(shifted) - visible(original)| / visible(original), all d
  IoU          overlap of the two visible sets, only for d = 30, 60 m (whole-cell shifts,
               so both grids line up exactly and no resampling error enters the comparison)

    python -m bench.observer           (resumable; outputs/observer/results.jsonl)
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from pyproj import Geod

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("VIEWSHED_CACHE", str(ROOT / "data" / "cache"))

from viewshed import INVISIBLE, VISIBLE, viewshed  # noqa: E402

from .sites import NAMED, RANDOM  # noqa: E402

OUT = ROOT / "outputs" / "observer"
OUT.mkdir(parents=True, exist_ok=True)
GEOD = Geod(ellps="WGS84")
R, H, CELL = 5000.0, 1.5, 30.0
SHIFTS = [5, 10, 20, 30, 60]
AZ = list(range(0, 360, 45))


def sites():
    s = dict(NAMED)
    seen = set()
    for name, v in RANDOM.items():          # first random site of each tile
        tile = name.split("-")[0]
        if tile not in seen:
            seen.add(tile)
            s[name] = v
    return s


def local_slope_deg(dem):
    o = dem.shape[0] // 2
    w = dem[o - 1:o + 2, o - 1:o + 2]
    gx = ((w[0, 2] + 2 * w[1, 2] + w[2, 2]) - (w[0, 0] + 2 * w[1, 0] + w[2, 0])) / (8 * CELL)
    gy = ((w[2, 0] + 2 * w[2, 1] + w[2, 2]) - (w[0, 0] + 2 * w[0, 1] + w[0, 2])) / (8 * CELL)
    return float(np.degrees(np.arctan(np.hypot(gx, gy))))


def shift_codes(codes, dx_cells, dy_cells):
    """Express a grid centred on the shifted observer in the original grid (whole-cell shifts)."""
    out = np.zeros_like(codes)
    n = codes.shape[0]
    # shifted grid origin is (dx, dy) metres east/north of the original -> columns +dx, rows -dy
    c, r = int(dx_cells), int(-dy_cells)
    src = codes[max(0, -r):n - max(0, r), max(0, -c):n - max(0, c)]
    out[max(0, r):max(0, r) + src.shape[0], max(0, c):max(0, c) + src.shape[1]] = src
    return out


def main():
    res_file = OUT / "results.jsonl"
    done = set()
    if res_file.exists():
        done = {(d["site"], d["shift_m"], d["azimuth"]) for d in map(json.loads, res_file.read_text(encoding="utf-8").splitlines())}
    for name, (lon, lat, region) in sites().items():
        base = viewshed(lon, lat, R, observer_height=H, engine="r3")
        bv = base.codes == VISIBLE
        ba = float(bv.sum() * CELL * CELL / 1e6)
        slope = local_slope_deg(base.dem)
        for d in SHIFTS:
            for az in AZ:
                if (name, d, az) in done:
                    continue
                lon2, lat2, _ = GEOD.fwd(lon, lat, az, d)
                s = viewshed(lon2, lat2, R, observer_height=H, engine="r3")
                sa = float((s.codes == VISIBLE).sum() * CELL * CELL / 1e6)
                row = dict(site=name, region=region, lon=lon, lat=lat, shift_m=d, azimuth=az, base_visible_km2=round(ba, 4),
                           shifted_visible_km2=round(sa, 4), area_change_pct=round(100 * abs(sa - ba) / ba, 3) if ba else None,
                           observer_slope_deg=round(slope, 2))
                if d % CELL == 0 and az % 90 == 0:     # whole-cell shift along an axis: exact alignment
                    dx = round(d * np.sin(np.radians(az)) / CELL)
                    dy = round(d * np.cos(np.radians(az)) / CELL)
                    sv = shift_codes(s.codes, dx, dy) == VISIBLE
                    both = np.isin(base.codes, (VISIBLE, INVISIBLE)) & np.isin(shift_codes(s.codes, dx, dy), (VISIBLE, INVISIBLE))
                    inter, union = (bv & sv & both).sum(), ((bv | sv) & both).sum()
                    row["iou"] = round(float(inter / union), 4) if union else None
                with res_file.open("a", encoding="utf-8") as f:
                    f.write(json.dumps(row) + "\n")
            print(f"{name[:36]:36s} shift {d:>2} m: median area change "
                  f"{np.median([r for r in [json.loads(l)['area_change_pct'] for l in res_file.read_text(encoding='utf-8').splitlines() if json.loads(l)['site']==name and json.loads(l)['shift_m']==d] if r is not None]):.1f}%",
                  flush=True)


if __name__ == "__main__":
    main()
