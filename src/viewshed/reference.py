"""Exact reference viewshed (R3, Franklin & Ray 1994), vectorised with numpy.

This is the yardstick for ViewshedBench: slow (O(n^1.5)) but with a fully declared
terrain model, so "accuracy vs R3" has a precise meaning.

Declared terrain model
----------------------
* Elevations are point samples at cell centres. The observer is the centre of its cell.
* The line of sight runs from the observer eye (ground + observer_height) to the target
  (ground + target_height), both at cell centres.
* The LOS is tested at EVERY crossing with a row line or a column line through cell
  centres. At a crossing, terrain height is linearly interpolated between the two
  adjacent cell centres on that line.
* Earth curvature and refraction: every height at ground distance d is lowered by
  curvature_coeff * d^2 / (2 R), with R the ellipsoid semi-major axis (as GDAL does).
* A target is visible iff its elevation angle is >= the maximum elevation angle of all
  crossings strictly between observer and target (ties count as visible).
* Missing elevation (NaN) at a crossing is ignored; a missing target is NODATA.

Alternative terrain models (model=...), used to test whether conclusions depend on the
reference's own assumptions. Everything else above stays the same:
* "linear"   (default) as described above.
* "nearest"  same crossings, but terrain height is the NEAREST cell centre (no interpolation).
* "flatcell" each cell is a flat square at its centre elevation (a "table-top" model). The LOS
             is tested at every crossing of a cell BOUNDARY (half-integer lines), where the
             terrain height is the higher of the two cells sharing that boundary. Boundaries of
             the target's own cell are excluded, so a cell never hides itself.
"""
from __future__ import annotations

import numpy as np

from .core import CURVATURE_COEFF, INVISIBLE, NODATA, OUT_OF_RANGE, VISIBLE

WGS84_A = 6378137.0


def sector_mask(shape, cell, start_angle=None, end_angle=None, min_distance=None, max_distance=None):
    """True where a cell centre lies inside the analysis region (observer at the centre cell)."""
    rows, cols = shape
    oi, oj = rows // 2, cols // 2
    di, dj = np.mgrid[0:rows, 0:cols]
    x, y = (dj - oj) * cell, (oi - di) * cell
    dist = np.hypot(x, y)
    m = np.ones(shape, bool)
    if max_distance is not None:
        m &= dist <= max_distance
    if min_distance:
        m &= dist >= min_distance
    if start_angle is not None:
        az = (np.degrees(np.arctan2(x, y)) + 360) % 360
        s, e = start_angle % 360, end_angle % 360
        m &= ((az - s) % 360) <= ((e - s) % 360)
    return m


def r3(dem: np.ndarray, cell: float, observer_height: float = 1.5, target_height: float = 0.0,
       max_distance: float | None = None, curvature_coeff: float = CURVATURE_COEFF,
       radius_earth: float = WGS84_A, nodata: float | None = -32768.0,
       return_clearance: bool = False, model: str = "linear"):
    """Exact R3 viewshed of the centre cell of `dem` (odd-sized, square cells of size `cell`).

    With return_clearance=True, also returns the sightline clearance in metres: how far the
    target could be lowered (if positive) or would have to be raised (if negative) before its
    visibility flips. |clearance| near 0 means a grazing sightline, where any small modelling
    difference between engines decides the answer. Cells next to the observer get +inf.
    """
    z = dem.astype(np.float64).copy()
    if nodata is not None:
        z[dem <= nodata + 0.5] = np.nan
    rows, cols = z.shape
    oi, oj = rows // 2, cols // 2
    k_curv = curvature_coeff / (2.0 * radius_earth)
    z0 = z[oi, oj] + observer_height

    di, dj = np.mgrid[0:rows, 0:cols]
    di, dj = di - oi, dj - oj
    dist = np.hypot(di, dj) * cell
    with np.errstate(invalid="ignore", divide="ignore"):
        tgt = (z + target_height - k_curv * dist**2 - z0) / dist
    maxslope = np.full(z.shape, -np.inf)
    adj, adi = np.abs(dj), np.abs(di)
    if model not in ("linear", "nearest", "flatcell"):
        raise ValueError("model must be 'linear', 'nearest' or 'flatcell'")

    if model == "flatcell":
        _flatcell_crossings(z, di, dj, adi, adj, oi, oj, dist, z0, k_curv, maxslope)
    else:
        nearest = model == "nearest"
        # Crossings with column lines (x = integer offsets): interpolate along rows.
        for k in range(1, int(adj.max())):
            m = adj > k
            t = k / adj[m]
            col = oj + np.sign(dj[m]) * k
            r = oi + di[m] * t
            maxslope[m] = np.fmax(maxslope[m], _slope(z, r, col, t * dist[m], z0, k_curv, True, nearest))
        # Crossings with row lines (y = integer offsets): interpolate along columns.
        for k in range(1, int(adi.max())):
            m = adi > k
            t = k / adi[m]
            row = oi + np.sign(di[m]) * k
            c = oj + dj[m] * t
            maxslope[m] = np.fmax(maxslope[m], _slope(z, row, c, t * dist[m], z0, k_curv, False, nearest))

    codes = np.where(tgt >= maxslope - 1e-12, VISIBLE, INVISIBLE).astype(np.uint8)
    codes[oi, oj] = VISIBLE
    codes[np.isnan(z)] = NODATA
    if max_distance is not None:
        codes[dist > max_distance] = OUT_OF_RANGE
    if return_clearance:
        with np.errstate(invalid="ignore"):
            clearance = (tgt - maxslope) * dist
        clearance[oi, oj] = np.inf
        clearance[codes == OUT_OF_RANGE] = np.nan
        clearance[codes == NODATA] = np.nan
        return codes, clearance
    return codes


def _slope(z, r, c, d, z0, k_curv, along_rows, nearest=False):
    """Elevation angle (as slope) of terrain at fractional positions on grid lines."""
    if along_rows:  # c is an exact column index, r fractional
        c = c.astype(np.intp)
        if nearest:  # exactly halfway between two centres: take the higher (symmetric tie rule)
            lo, hi = np.ceil(r - 0.5).astype(np.intp), np.floor(r + 0.5).astype(np.intp)
            h = np.fmax(z[lo, c], z[hi, c])
        else:
            r0 = np.floor(r).astype(np.intp)
            f = r - r0
            r1 = np.minimum(r0 + 1, z.shape[0] - 1)
            h = z[r0, c] * (1 - f) + np.where(f > 0, z[r1, c] * f, 0.0)
    else:           # r is an exact row index, c fractional
        r = r.astype(np.intp)
        if nearest:
            lo, hi = np.ceil(c - 0.5).astype(np.intp), np.floor(c + 0.5).astype(np.intp)
            h = np.fmax(z[r, lo], z[r, hi])
        else:
            c0 = np.floor(c).astype(np.intp)
            f = c - c0
            c1 = np.minimum(c0 + 1, z.shape[1] - 1)
            h = z[r, c0] * (1 - f) + np.where(f > 0, z[r, c1] * f, 0.0)
    return (h - k_curv * d**2 - z0) / d


def _flatcell_crossings(z, di, dj, adi, adj, oi, oj, dist, z0, k_curv, maxslope):
    """Table-top model: test the LOS at every cell-boundary crossing; the blocking height is
    the higher of the two cells sharing the boundary. Boundaries of the target cell are skipped."""
    ti, tj = oi + di, oj + dj
    # Column boundaries x = oj + s*(k - 0.5), between columns oj + s*(k-1) and oj + s*k.
    for k in range(1, int(adj.max()) + 1):
        m = adj >= k
        t = (k - 0.5) / adj[m]
        s = np.sign(dj[m])
        ca, cb = oj + s * (k - 1), oj + s * k
        r = oi + di[m] * t
        # A crossing exactly at a cell corner touches four cells: use all of them (symmetric).
        rl, rh = np.ceil(r - 0.5).astype(np.intp), np.floor(r + 0.5).astype(np.intp)
        h = np.fmax(np.fmax(z[rl, ca], z[rl, cb]), np.fmax(z[rh, ca], z[rh, cb]))
        d = t * dist[m]
        sl = (h - k_curv * d**2 - z0) / d
        sl[(cb == tj[m]) & ((rl == ti[m]) | (rh == ti[m]))] = -np.inf
        maxslope[m] = np.fmax(maxslope[m], sl)
    # Row boundaries y = oi + s*(k - 0.5), between rows oi + s*(k-1) and oi + s*k.
    for k in range(1, int(adi.max()) + 1):
        m = adi >= k
        t = (k - 0.5) / adi[m]
        s = np.sign(di[m])
        ra, rb = oi + s * (k - 1), oi + s * k
        c = oj + dj[m] * t
        cl, ch = np.ceil(c - 0.5).astype(np.intp), np.floor(c + 0.5).astype(np.intp)
        h = np.fmax(np.fmax(z[ra, cl], z[rb, cl]), np.fmax(z[ra, ch], z[rb, ch]))
        d = t * dist[m]
        sl = (h - k_curv * d**2 - z0) / d
        sl[(rb == ti[m]) & ((cl == tj[m]) | (ch == tj[m]))] = -np.inf
        maxslope[m] = np.fmax(maxslope[m], sl)
