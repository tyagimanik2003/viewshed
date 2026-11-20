# Changelog

## 0.1.0.dev0 (unreleased)

First development version.

- `viewshed.viewshed()`: viewsheds from longitude/latitude with an automatic, observer-centred
  azimuthal-equidistant grid; circular radius, sectors (geodesic azimuths), minimum distance,
  observer and target heights, Earth curvature and refraction.
- Copernicus GLO-30 fetched on demand with a local tile cache; local DEM files also accepted.
- Distinct output codes for visible, invisible, out-of-range and no-data cells.
- Two engines: GDAL (`gdal raster viewshed` via `gdal.Run`) and an exact R3 reference with three
  declared terrain models (`linear`, `nearest`, `flatcell`) and per-cell sightline clearance.
- `viewshed.ztv()`: zone of theoretical visibility, giving per-cell counts of visible observers
  (e.g. wind turbines) on a shared grid, areas by count, the area each observer is visible from,
  and GeoTIFF output.
- `viewshed.export`: Web Mercator PNG overlays with corner coordinates for web maps.
- `viewshed()` raises a clear error when the observer has no elevation data (e.g. outside a local DEM).
- Command line: `viewshed LON LAT [--radius] [--sector] [--height] [--png]`.
- Works with GDAL 3.12 and 3.13, including workarounds for documented GDAL behaviour
  (`docs/gdal_notes.md`).
- ViewshedBench (`bench/`): benchmark of GDAL, GRASS r.viewshed, WhiteboxTools and
  xarray-spatial against the exact reference.
- Test suite against analytic geometry; CI on Linux.
- Released under the MIT licence.
