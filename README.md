# SOCIO-FUNCTIONAL CLUSTERING

Code to detect socio-functional areas by integrating morphological analysis 
(MomePy) and urban activity data (Foursquare OS Places).

## Description

The code processes building and street geometries (from cadastral sources or 
OSM) and computes morphological metrics via MomePy tessellations. It then 
integrates Foursquare LBSN data to derive activity diversity indices at the 
tessellation-cell level. Two clustering procedures are implemented:

- **Part 1 – Morphological clustering**: based solely on urban form metrics
- **Part 2 – Socio-functional clustering**: combines morphological metrics 
  with activity variables (counts, check-ins, diversity indices)

The code is based on [MomePy](https://github.com/pysal/momepy).

## Required libraries

- `geopandas`, `osmnx`, `momepy`, `libpysal`
- `clustergram`, `scipy`, `numpy`, `pandas`
- `matplotlib`, `bokeh`, `tqdm`, `shapely`

## Required input files

- Excel file with the list of places (`Case-studies.xlsx`), including city 
  name, country, and number of clusters for each part
- Building geometries: cadastral GeoJSON or downloaded from OSM
- Street network: GeoJSON or downloaded from OSM
- Foursquare places dataset per city (GeoJSON with `catName1`, 
  `Activity_Type`, and `checkins:2023` fields) — Part 2 only

## Usage

Set `INPUT_FOLDER`, `RESULTS_FOLDER` and `LOCAL_CRS` at the top of the 
script. Run Part 1 independently of Part 2. Results are exported as GeoJSON 
files and cluster score plots as PNG.

## Citation

If you use this code, please cite:

Nolasco-Cirugeda, A., López-Baeza, J., García-Pérez, S., & García-Mayor, C. (2026). Shaping socio-functional urban identity: When morphometrics meet social media data to reveal neighbourhood structure and cores. *Cities*. https://doi.org/10.1016/j.cities.2026.107373

Code repository: https://github.com/mappingame/socio-functional-clustering
Code DOI: https://doi.org/10.5281/zenodo.21258671

## Contributors

[Jesús López-Baeza](https://github.com/JesusLopezBaeza)
[Sergio García-Pérez](https://github.com/s-garciap)
[Almudena Nolasco-Cirugeda](https://github.com/mappingame)

