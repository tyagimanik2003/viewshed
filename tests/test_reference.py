"""The exact R3 reference must match geometry exactly, and obey invariants on any terrain."""
from __future__ import annotations

import math

import numpy as np
import pytest

from conftest import OBS_LAT, OBS_LON
from viewshed import INVISIBLE, OUT_OF_RANGE, VISIBLE, viewshed
from viewshed.reference import WGS84_A, r3

CELL = 30.0
CURV = 1 - 1 / 7


def rough_terrain(n=161, seed=0):
    """Smooth random mountains: sum of Gaussian bumps, up to ~600 m relief."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:n, 0:n]
    z = np.full((n, n), 3500.0)
    for _ in range(25):
        cy, cx = rng.uniform(0, n, 2)
        s = rng.uniform(4, 25)
        z += rng.uniform(-250, 400) * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * s * s))
    return z


def test_flat_plateau_horizon_matches_theory():
    n = 401
    codes = r3(np.full((n, n), 4000.0), CELL, observer_height=1.5, curvature_coeff=CURV)
    row = codes[n // 2, n // 2:]  # due east
    last_visible = np.nonzero(row == VISIBLE)[0].max() * CELL
    theory = math.sqrt(2 * WGS84_A * 1.5 / CURV)
    assert abs(last_visible - theory) <= CELL


def test_wall_shadow_starts_where_geometry_says():
    """Flat earth (no curvature), eye 10 m, wall 5 m high at 300 m.
    Ground behind the wall is hidden until the ray over the wall top lands: 10/D = 5/300 -> D = 600 m."""
    n = 101
    dem = np.zeros((n, n))
    o = n // 2
    dem[:, o + 10] = 5.0  # wall 10 cells = 300 m east
    codes = r3(dem, CELL, observer_height=10.0, curvature_coeff=0.0)
    east = codes[o, o:]
    dist = np.arange(east.size) * CELL
    assert (east[dist <= 300] == VISIBLE).all()
    assert (east[(dist > 300) & (dist < 600)] == INVISIBLE).all()
    assert (east[dist >= 600] == VISIBLE).all()  # exactly 600 m: grazing ray counts as visible


def test_clearance_matches_wall_geometry():
    """Eye 10 m, flat ground, 5 m wall at 300 m. The steepest blocker of a ground target at D is
    either the wall top (slope -5/300) or the ground cell just before the target (slope
    -10/(D - cell)). Clearance = D * (target slope - steepest blocker slope).

    The second case is the flat-land effect: on open flat ground clearance is only
    h * cell / (D - cell), i.e. ~0.4 m at 750 m here and ~4.5 cm at 1 km for a 1.5 m observer,
    so centimetre-level modelling differences decide visibility."""
    n = 101
    dem = np.zeros((n, n))
    o = n // 2
    dem[:, o + 10] = 5.0
    codes, clr = r3(dem, CELL, observer_height=10.0, curvature_coeff=0.0, return_clearance=True)
    for k in (12, 15, 18, 25, 40):            # 360, 450, 540, 750, 1200 m east
        d = k * CELL
        blocker = max(-5 / 300, -10 / (d - CELL))
        expected = d * (-10 / d - blocker)
        assert clr[o, o + k] == pytest.approx(expected, abs=1e-6)
        assert (codes[o, o + k] == VISIBLE) == (expected >= -1e-9)
    assert clr[o, o + 25] == pytest.approx(10 * CELL / (25 * CELL - CELL), abs=1e-6)  # ~0.42 m


@pytest.mark.parametrize("b", [0.3, 1.0])
def test_single_bump_shadow_length_on_axis(b):
    """Flat ground, eye h = 1.5 m, one cell raised by b at d = 600 m due east. The ray to a
    ground target at D passes the bump at height h(1 - d/D), so the target is hidden for
    d < D < d*h/(h-b): 750 m for b = 0.3, 1800 m for b = 1.0."""
    h, n = 1.5, 161
    o = n // 2
    dem = np.full((n, n), 300.0)
    dem[o, o + 20] += b
    codes = r3(dem, CELL, observer_height=h, curvature_coeff=0.0)
    east = codes[o, o:]
    dist = np.arange(east.size) * CELL
    end = 600 * h / (h - b)
    assert (east[(dist > 600) & (dist < end - 1e-6)] == INVISIBLE).all()
    assert (east[dist >= end - 1e-6] == VISIBLE).all()
    assert (east[dist <= 600] == VISIBLE).all()


MODELS = ["linear", "nearest", "flatcell"]


@pytest.mark.parametrize("model", MODELS)
def test_flat_plateau_horizon_all_models(model):
    n = 401
    codes = r3(np.full((n, n), 4000.0), CELL, observer_height=1.5, curvature_coeff=CURV, model=model)
    last_visible = np.nonzero(codes[n // 2, n // 2:] == VISIBLE)[0].max() * CELL
    assert abs(last_visible - math.sqrt(2 * WGS84_A * 1.5 / CURV)) <= CELL


def test_flatcell_wall_shadow():
    """Table-top model: the 5 m wall cell spans 285-315 m; its far top edge at 315 m sets the
    shadow, which ends where 10/D = 5/315, i.e. D = 630 m (linear model: 600 m)."""
    n = 101
    dem = np.zeros((n, n))
    o = n // 2
    dem[:, o + 10] = 5.0
    codes = r3(dem, CELL, observer_height=10.0, curvature_coeff=0.0, model="flatcell")
    east = codes[o, o:]
    dist = np.arange(east.size) * CELL
    assert (east[dist <= 300] == VISIBLE).all()
    assert (east[(dist > 300) & (dist < 630)] == INVISIBLE).all()
    assert (east[dist >= 630] == VISIBLE).all()


@pytest.mark.parametrize("model", ["nearest", "flatcell"])
def test_alternative_models_rotation_and_monotonicity(model):
    dem = rough_terrain(seed=7)
    assert (np.rot90(r3(dem, CELL, model=model)) == r3(np.rot90(dem), CELL, model=model)).mean() > 0.999
    low, high = r3(dem, CELL, observer_height=2.0, model=model), r3(dem, CELL, observer_height=20.0, model=model)
    assert not ((low == VISIBLE) & (high != VISIBLE)).any()


def test_unknown_model_rejected():
    with pytest.raises(ValueError):
        r3(np.zeros((11, 11)), CELL, model="cubic")


def test_clearance_sign_agrees_with_codes():
    dem = rough_terrain(seed=4)
    codes, clr = r3(dem, CELL, return_clearance=True)
    ok = np.isfinite(clr)
    assert ((clr[ok] >= -1e-6) == (codes[ok] == VISIBLE)).all()


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_rotation_invariance(seed):
    dem = rough_terrain(seed=seed)
    a = np.rot90(r3(dem, CELL))
    b = r3(np.rot90(dem), CELL)
    assert (a == b).mean() > 0.999


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_raising_the_observer_never_hides_anything(seed):
    dem = rough_terrain(seed=seed)
    low, high = r3(dem, CELL, observer_height=2.0), r3(dem, CELL, observer_height=20.0)
    assert not ((low == VISIBLE) & (high != VISIBLE)).any()
    assert (high == VISIBLE).sum() > (low == VISIBLE).sum()


def test_engine_r3_passes_the_same_specs_as_gdal(flat_geographic):
    res = viewshed(OBS_LON, OBS_LAT, 5000, dem=flat_geographic, engine="r3")
    assert res.analysed_km2 == pytest.approx(math.pi * 25, rel=0.01)
    gdal_res = viewshed(OBS_LON, OBS_LAT, 5000, dem=flat_geographic, engine="gdal")
    both = (res.codes != OUT_OF_RANGE) & (gdal_res.codes != OUT_OF_RANGE)
    assert (res.codes[both] == gdal_res.codes[both]).mean() > 0.99


def test_engine_r3_sector(flat_geographic):
    res = viewshed(OBS_LON, OBS_LAT, 5000, dem=flat_geographic, engine="r3", start_angle=300, end_angle=60)
    assert res.analysed_km2 == pytest.approx(math.pi * 25 / 3, rel=0.03)


def test_unknown_engine_rejected(flat_geographic):
    with pytest.raises(ValueError):
        viewshed(OBS_LON, OBS_LAT, 5000, dem=flat_geographic, engine="magic")
