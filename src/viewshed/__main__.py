"""Command line:  python -m viewshed LON LAT [--radius 5000] [--sector 300 60] [--png map.png]"""
from __future__ import annotations

import argparse

from .core import viewshed


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="viewshed", description="What can an observer at LON, LAT see?")
    ap.add_argument("lon", type=float)
    ap.add_argument("lat", type=float)
    ap.add_argument("--radius", type=float, default=5000, help="metres (default 5000)")
    ap.add_argument("--height", type=float, default=1.5, help="observer height above ground, m (default 1.5)")
    ap.add_argument("--target-height", type=float, default=0.0)
    ap.add_argument("--sector", type=float, nargs=2, metavar=("START", "END"),
                    help="azimuths in degrees clockwise from north; sector runs clockwise START -> END")
    ap.add_argument("--min-distance", type=float)
    ap.add_argument("--dem", nargs="+", default=["glo30"], help="'glo30' (default) or raster path(s)")
    ap.add_argument("--resolution", type=float, default=30.0, help="analysis cell size, m")
    ap.add_argument("--out", help="output GeoTIFF (AEQD grid centred on the observer)")
    ap.add_argument("--png", help="quick-look map")
    a = ap.parse_args(argv)

    dem = a.dem[0] if len(a.dem) == 1 else a.dem
    res = viewshed(a.lon, a.lat, a.radius, observer_height=a.height, target_height=a.target_height,
                   start_angle=a.sector[0] if a.sector else None, end_angle=a.sector[1] if a.sector else None,
                   min_distance=a.min_distance, dem=dem, resolution=a.resolution, out=a.out)
    print(res.summary())
    if a.out:
        print(f"GeoTIFF: {res.path}")
    if a.png:
        from .plot import plot

        plot(res, a.png)
        print(f"Map: {a.png}")


if __name__ == "__main__":
    main()
