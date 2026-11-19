"""Stratified random benchmark sites: K random observers inside each chosen 1-degree tile.

Random placement (fixed seed) avoids cherry-picking. Each candidate's 10 km circle must lie
inside its tile, have < 1% missing elevation, and not be dominated by perfectly flat
patches (a crude water / fill-flat filter). Rejected candidates are redrawn.

    python -m bench.sample_sites        -> bench/sites_random.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("VIEWSHED_CACHE", str(ROOT / "data" / "cache"))

from viewshed.dem import observer_grid  # noqa: E402

# (lat0, lon0): landscape label, for reading only. Analysis uses measured relief/roughness.
TILES = {
    (32, 77): "high Himalaya (Lahaul)", (32, 78): "high Himalaya (Spiti)", (31, 78): "high Himalaya (Kinnaur)",
    (30, 79): "high Himalaya (Garhwal)", (33, 76): "high Himalaya (Zanskar)", (31, 77): "high hills (Kullu-Shimla)",
    (30, 78): "high hills (lower Garhwal)", (34, 77): "Ladakh plateau", (32, 79): "Ladakh plateau (Changthang)",
    (32, 76): "foothills (Kangra)", (31, 76): "foothills (Mandi)", (30, 77): "foothills (Sirmaur-Dehradun)",
    (29, 76): "plains (Haryana)", (30, 75): "plains (Punjab)", (28, 73): "plains (Thar)",
    (26, 80): "plains (Gangetic, Kanpur)", (17, 75): "plateau (Deccan)", (10, 76): "mountains (Western Ghats)",
    (27, 88): "high Himalaya (Sikkim)",
}
PER_TILE, SEED, MARGIN, CHECK_R = 3, 2026, 0.13, 10000.0


def acceptable(dem: np.ndarray) -> bool:
    n = dem.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    m = np.hypot(xx - n // 2, yy - n // 2) * 30 <= CHECK_R
    if (dem[m] <= -32767).mean() > 0.01:
        return False
    flat = (np.abs(np.diff(dem, axis=1))[m[:, 1:]] < 1e-3).mean()  # identical neighbours: water/fill
    return flat < 0.3


def main():
    rng = np.random.default_rng(SEED)
    sites = {}
    for (lat0, lon0), label in TILES.items():
        k = tries = 0
        while k < PER_TILE and tries < 15:
            tries += 1
            lat = round(float(rng.uniform(lat0 + MARGIN, lat0 + 1 - MARGIN)), 4)
            lon = round(float(rng.uniform(lon0 + MARGIN, lon0 + 1 - MARGIN)), 4)
            ds = observer_grid(lon, lat, CHECK_R, dem="glo30")
            dem = ds.GetRasterBand(1).ReadAsArray()
            ds = None
            if not acceptable(dem):
                print(f"  rejected {lat},{lon} ({label})")
                continue
            k += 1
            sites[f"R{lat0:02d}{lon0:03d}-{k} {label}"] = [lon, lat, label.split(" (")[0]]
            print(f"site {len(sites):2d}: {lat:.4f} N {lon:.4f} E  {label}", flush=True)
    out = Path(__file__).with_name("sites_random.json")
    out.write_text(json.dumps(sites, indent=1), encoding="utf-8")
    print(f"{len(sites)} sites -> {out}")


if __name__ == "__main__":
    main()
