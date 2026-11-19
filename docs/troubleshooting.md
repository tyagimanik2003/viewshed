# Troubleshooting

## GRASS r.viewshed fails on Windows (conda-forge)

With conda-forge `grass 8.5.0` on Windows, every GRASS module can exit with code 0xC0000135
(DLL not found). `ldd` on a module such as `g.proj.exe` shows:

    libomp140.x86_64.dll => not found

This is Microsoft's LLVM OpenMP runtime, which ships with the Visual Studio C++ build tools and
is not included in conda-forge's `vc14_runtime`. A second, smaller issue: `grass.bat` sets
`PROJ_LIB=` (empty); setting `PROJ_DATA` to `<env>/Library/share/proj` works around it.

Options:
1. Install the Visual Studio C++ build tools (provides `libomp140.x86_64.dll`).
2. Install GRASS through OSGeo4W, which bundles its own runtimes.
3. Run on Linux. The GitHub Actions workflow `.github/workflows/grass-bench.yml` does this.

## The viewshed looks too generous

If you call GDAL directly, check that the observer height is passed as the third value of the
position and not with `--height` (see `docs/gdal_notes.md`, section 1).

## A viewshed on a latitude/longitude DEM looks wrong

Use `viewshed.viewshed()`, which reprojects automatically, or reproject to a metric projection
before calling GDAL yourself (see `docs/gdal_notes.md`, section 4).

## First run is slow

The first run downloads Copernicus GLO-30 tiles (about 40 MB per 1° tile) into the cache
(`$VIEWSHED_CACHE`, default `~/.cache/viewshed`). Later runs at the same location take about a
second. Set `VIEWSHED_CACHE=off` to read tiles over HTTP without caching.
