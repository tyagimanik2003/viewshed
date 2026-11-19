# Notes on GDAL's viewshed behaviour

The GDAL backend of this package calls `gdal.Run("raster", "viewshed", ...)`. These notes
record GDAL behaviours that affect correctness, and how the package handles them. Each one is
covered by a test, so a GDAL upgrade that changes the behaviour shows up as a test failure.

## 1. `-z/--height` is ignored (GDAL 3.12.3 and 3.13.3)

When the observer position is given with two coordinates, `gdal raster viewshed` overwrites the
`--height` value with a 2 m default (`apps/gdalalg_raster_viewshed.cpp`: the `m_opts.observer.z = 2`
fallback runs even when `--height` was set explicitly). Only the third value of `--position`
is honoured.

Flat 4,000 m plateau, azimuthal-equidistant grid at 30 m, curvature coefficient 1 − 1/7.
"Farthest visible" is the distance of the last visible cell; theory is √(2Rh / c).

| Observer height h | `--height h` | `--position x,y,h` | Theory |
|---|---|---|---|
| 0.5 m | 5,476 m | 2,748 m | 2,728 m |
| 1.5 m | 5,476 m | 4,745 m | 4,725 m |
| 10 m | 5,476 m | 8,485 m (grid edge) | 12,199 m |

5,476 m is the horizon of a 2 m observer, whatever `--height` says.

**Handling:** the package always passes the observer height as the third value of `position`.
`tests/test_gdal_quirks.py` pins both behaviours.

## 2. Sector and minimum-distance options from Python

GDAL ≥ 3.12 supports `--start-angle`, `--end-angle` and `--min-distance`. Whether they take effect
depends on the route used (GDAL 3.13.3; flat plateau, sector 0–90°, min distance 1 km):

| Route | Sector honoured | Min distance honoured |
|---|---|---|
| `gdal.ViewshedGenerate(...)` (classic binding; no angle arguments) | no | no |
| `gdal.ViewshedGenerate(..., options=["START_ANGLE=..."])` | yes | no |
| `gdal raster viewshed` command line | yes | yes |
| `gdal.Run("raster", "viewshed", ...)` | yes | yes |

**Handling:** the package uses `gdal.Run`. Reproduce with `python scripts/check_gdal_sector_options.py`.
See also GDAL issue #13451.

## 3. Option names differ between versions

- `--dst-nodata` (3.12) was renamed `--output-nodata` (3.13). The package picks the right name
  for the installed version.
- `--maybe-visible-value` exists only from 3.13. The package doesn't pass it; its default (2)
  equals `viewshed.MAYBE_VISIBLE`.

## 4. Geographic (latitude/longitude) DEMs

GDAL's viewshed assumes a projected grid: on a latitude/longitude DEM, horizontal distances are
in degrees while heights are in metres, so curvature and maximum distance are wrong. The package
always reprojects to an azimuthal-equidistant grid centred on the observer, and
`tests/test_engine.py::test_lat_lon_dem_gives_physical_viewshed` checks the result.

## 5. Results depend on the GDAL/PROJ version used to build grids

Reprojecting the same Copernicus GLO-30 tiles with GDAL 3.12.3 / PROJ 9.7.1 and with
GDAL 3.13.3 / PROJ 9.9.0 gives grids that differ by centimetres on average. At sites where many
sightlines graze the terrain, this can change the visible area noticeably. Record the GDAL and
PROJ versions alongside any published viewshed result.
