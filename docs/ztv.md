# Zone of theoretical visibility (ZTV)

A ZTV maps, for every location around a development, how many of its elements (turbines, masts,
buildings) could be seen from there. It is the standard first step of a landscape and visual
impact assessment.

## Method

`viewshed.ztv(observers, radius, target_height=...)`:

1. **Shared grid.** All observers share one azimuthal-equidistant grid centred on their spherical
   centroid. Cell edges sit at half-cell offsets, as for single viewsheds. The grid covers every
   observer's radius plus a two-cell margin, and DEM tiles are fetched for that whole area.
2. **One viewshed per observer.** Each observer runs through `gdal raster viewshed` on the shared
   grid, with its own height, the common target height, Earth curvature and refraction
   (coefficient 1 − 1/7), and the radius as maximum distance.
3. **Counting.** A cell's count is the number of observers for which it is visible. A cell is
   *analysed* if it lies within range of at least one observer and has elevation data.

## Interpreting the heights

Line of sight is symmetric: the segment from a turbine tip to a person's eyes is the same segment
in both directions. So "the tip sees a 2 m target at this cell" is exactly "a person standing at
this cell sees the tip".

- **Observer height:** the part of the structure whose visibility matters, usually blade-tip
  height. Hub height answers "is the nacelle visible?".
- **Target height:** viewer eye height, conventionally 1.5–2 m.

## Assumptions and limits

- **Bare-earth versus surface model.** Copernicus GLO-30 is a surface model: it includes
  buildings and some vegetation, but at 30 m resolution it smooths them out. Treat the result as
  a *theoretical* ZTV, which is how ZTVs are normally presented. Screening by local features
  needs field work or higher-resolution data.
- **Visibility is binary.** A count says whether any part of the tip is above the horizon, not
  how much of the turbine shows or how prominent it is.
- **Spread of observers.** They must lie within 100 km of their centroid. Over that span the
  shared projection distorts distances by less than about 1 part in 10,000.
- **Engine.** Counts use GDAL's viewshed algorithm. Engines differ in how they treat terrain
  between cell centres, and the differences are largest on flat, rough land (see the benchmark).

## Validation

`tests/test_ztv.py` checks that:
- with one observer, the ZTV equals the single viewshed exactly;
- on a flat plateau, each count matches the analytic curvature horizon of every observer;
- areas by count add up to the visible area;
- the GeoTIFF output round-trips.

On real Himalayan terrain, the counts also matched the sum of five independently computed single
viewsheds, sampled on the ZTV grid, for 99.4% of cells. The remaining cells differed by one, at
the boundaries between the two grids.
