# ViewshedBench

A reproducible benchmark of viewshed engines against an exact reference. It lives in `bench/`
and is separate from the `viewshed` package.

## What it measures

Every engine receives the **same** terrain grid: Copernicus GLO-30, reprojected to a 30 m
azimuthal-equidistant grid centred on the observer, with the observer on the centre cell. Each
engine's visible/invisible classification is compared with the exact R3 reference in
`viewshed.reference`, using agreement, F1 of the visible class, signed area difference and
absolute area difference.

| Engine | Adapter | Notes |
|---|---|---|
| R3 (reference) | `viewshed.reference.r3` | Exact under a declared terrain model: `linear` (default), `nearest`, `flatcell` |
| GDAL ≥ 3.12 | `gdal raster viewshed` via `gdal.Run` | Observer height passed in `position` (see `gdal_notes.md`) |
| GRASS r.viewshed 8.5 | GRASS session | Curvature `-c`, refraction `-r` (k = 1/7) |
| WhiteboxTools 2.4 | `Viewshed` | No curvature option; radius masked afterwards |
| xarray-spatial 0.10 | `xrspatial.viewshed` | No curvature option; radius masked afterwards |

Scenarios: **flat** (no Earth curvature; all engines) and **curved** (curvature with refraction
k = 1/7; engines that support it). Observer height 1.5 m, target height 0, radii 5 and 10 km.

## Sites

`bench/sites.py`: 14 named sites spanning relief classes in India, plus stratified random sites
from `bench/sites_random.json` (three random positions in each of 19 one-degree tiles, generated
by `python -m bench.sample_sites` with a fixed seed).

## Environments

- `environment.yml`: the package, GDAL and the test suite.
- `environment-bench.yml`: adds GRASS, WhiteboxTools and xarray-spatial. On Windows GRASS may
  not run; see `troubleshooting.md`.

## Running it

```bash
python -m bench.run                      # all sites, radii 5 and 10 km, all engines in the environment
python -m bench.run --sites 3 --radii 5000 --engines gdal
python -m bench.analyse                  # tables and figures in outputs/bench/
```

Runs are resumable: results append to `outputs/bench/results.jsonl` and finished runs are
skipped. Terrain tiles are cached under `data/cache/`; grids and engine outputs under
`data/bench/`. Both folders are git-ignored.

On Linux, `.github/workflows/grass-bench.yml` runs the R3 reference, GRASS and GDAL for all sites.

## Further experiments

| Module | Question |
|---|---|
| `bench.explain` | Per-cell sightline clearance; synthetic terrains varying relief, roughness and observer height |
| `bench.bumps` | Shadows behind single obstacles and on rough flat ground |
| `bench.observer` | Sensitivity of the exact viewshed to small shifts of the observer |
| `bench.refcheck` | Sensitivity of engine scores to the reference's terrain model |
| `bench.gdal_version_check` | Whether GDAL-dependent results change between GDAL versions |

## Reproducibility notes

- Record GDAL and PROJ versions with any result: grids built with different versions can differ
  slightly, which matters where sightlines graze the terrain.
- Percent-based scores are unstable when very little is visible. `bench.analyse` reports them
  only for runs with at least 1 km² visible in the reference, and reports absolute areas for all
  runs.
