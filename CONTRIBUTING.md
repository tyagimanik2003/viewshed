# Contributing

Bug reports, test cases and pull requests are welcome.

## Development setup

```bash
conda env create -f environment.yml
conda activate viewshed
pip install --no-deps -e .
python -m pytest
```

## Guidelines

- **Tests against geometry.** New behaviour should come with a test whose expected answer can be
  derived independently (an analytic horizon, a wall shadow, an invariance), not only a
  regression snapshot.
- **No silent changes to results.** If a change alters viewshed output, say so in the pull
  request and update `CHANGELOG.md`.
- **Record versions.** Results depend on GDAL/PROJ versions; include them in bug reports.
- **Keep DEM data out of the repository.** Tiles are fetched at run time; licences differ per DEM.
- **Disclose AI assistance** in pull requests, stating what was generated and what was reviewed.

## Reporting a bug

Include the coordinates, radius and options, the DEM source, `python -c "from osgeo import gdal;
print(gdal.__version__)"`, and what you expected versus what you got.
