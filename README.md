# viewshed

[![tests](https://github.com/tyagimanik2003/viewshed/actions/workflows/tests.yml/badge.svg)](https://github.com/tyagimanik2003/viewshed/actions/workflows/tests.yml)

Correct, reproducible viewsheds from latitude/longitude, anywhere on Earth: *what can an
observer standing here see?*

Give it a longitude, latitude and radius. It fetches the terrain, builds a metric grid centred on
the observer, and returns a viewshed you can check: distances and azimuths are true, the circle
is a circle, and visible cells can never be confused with missing data.

```bash
viewshed 77.015 32.329 --radius 5000 --png map.png
```

```python
from viewshed import viewshed

r = viewshed(77.015, 32.329, radius=5000, observer_height=1.5, start_angle=300, end_angle=60)
print(r.summary())        # visible area, analysed area, percentage
r.codes                   # 3 visible, 1 invisible, 0 out of range/sector, 255 no data
```

> **Status:** alpha (0.1.0.dev0). The API may change. See [CHANGELOG.md](CHANGELOG.md).

## Why another viewshed tool?

Mature engines (GDAL, GRASS, WhiteboxTools, xarray-spatial) compute viewsheds quickly. They also
assume a projected grid and leave terrain acquisition, projection and output conventions to the
user, and viewsheds computed on latitude/longitude grids are silently wrong. This package focuses
on getting the workflow right, and on being able to show that it is right:

- **Lat/lon in, metric geometry inside.** Each run uses an azimuthal-equidistant grid centred on
  the observer, so radius and sector azimuths are geodesic.
- **Terrain on demand.** Copernicus GLO-30 tiles are downloaded and cached automatically. Local
  DEMs (e.g. SRTM tiles) also work.
- **An exact reference.** `engine="r3"` runs an exact R3 line-of-sight algorithm with a declared
  terrain model, and can report each cell's sightline clearance.
- **Documented GDAL behaviour.** Known pitfalls (an ignored height option, version-specific option
  names) are handled and pinned by tests. See [docs/gdal_notes.md](docs/gdal_notes.md).
- **Tested against geometry.** The tests check analytic answers: the curvature horizon on a flat
  plateau, shadows behind walls and single obstacles, true circular area, rotation invariance,
  and monotonicity in observer height.

## Installation

GDAL comes from conda-forge, because binary GDAL wheels are not available on PyPI for every
platform.

```bash
git clone https://github.com/tyagimanik2003/viewshed.git
cd viewshed
conda env create -f environment.yml
conda activate viewshed
pip install --no-deps -e .
python -m pytest
```

Requires GDAL ≥ 3.12 (tested with 3.12.3 and 3.13.3).

## Usage

### Command line

```bash
viewshed LON LAT [--radius M] [--height M] [--target-height M]
                 [--sector START END] [--min-distance M]
                 [--dem glo30 | FILE ...] [--resolution M]
                 [--out viewshed.tif] [--png map.png]
```

Sector azimuths are degrees clockwise from true north; the sector runs clockwise from START to
END (e.g. `--sector 300 60` covers north).

### Python

```python
from viewshed import viewshed
from viewshed.plot import plot          # needs matplotlib

r = viewshed(lon, lat, radius=10_000,
             observer_height=10,            # e.g. a tower
             dem=["n32_e077.tif"],          # local DEMs instead of GLO-30
             engine="gdal")                 # or "r3" (exact, slower)
plot(r, "map.png")
```

`ViewshedResult` contains the codes, the DEM grid, the geotransform and CRS (azimuthal
equidistant, observer at x = y = 0), and helpers such as `visible_km2` and `summary()`.

### Exact reference

```python
from viewshed.reference import r3

codes, clearance = r3(dem_array, cell=30.0, observer_height=1.5,
                      max_distance=5000, model="linear", return_clearance=True)
```

`model` is `"linear"` (heights interpolated between cell centres), `"nearest"` or `"flatcell"`
(cells as flat squares). `clearance` is how far each target could be lowered (+) or must be
raised (−) before its visibility changes.

## ViewshedBench

`bench/` contains a reproducible benchmark that compares GDAL, GRASS r.viewshed, WhiteboxTools and
xarray-spatial against the exact reference on identical grids at sites across Indian terrain.
See [docs/benchmark.md](docs/benchmark.md).

## Data and licences

DEM tiles are never stored in this repository; they are fetched at run time. Copernicus GLO-30:
© DLR e.V. 2010–2014 and © Airbus Defence and Space GmbH 2014–2018, provided under COPERNICUS by
the European Union and ESA. Other DEMs have their own licences (e.g. FABDEM is non-commercial).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Troubleshooting notes are in
[docs/troubleshooting.md](docs/troubleshooting.md).

## Licence

Not yet chosen. Until a licence is added, all rights are reserved.
