"""
Socio-Functional Clustering Procedure Using Morphological and Activity Data

This code integrates MomePy with Foursquare OS Places data for the detection
of socio-functional areas at the urban tessellation level.

General description:
    The code searches from the list of places or specific files already
    downloaded with building and street grid geometries, or downloads them
    from OSM. It calculates morphological clusters (Part 1) or
    socio-functional clusters (Part 2) on the basis of MomePy, incorporating
    data from Foursquare OS Places.

    The code is based on https://github.com/pysal/momepy

Usage:
    Set the INPUT_FOLDER, RESULTS_FOLDER and LOCAL_CRS variables below to
    match your local environment before running.

    Part 1 (morphological clustering) and Part 2 (socio-functional clustering)
    can be run independently.

Input files required:
    - Excel file with list of places: Case-studies.xlsx
      Columns: City, Country, Form Clusters (Part 1), FS Clusters (Part 2)
    - Building geometries: cadastral GeoJSON or downloaded from OSM
    - Street network: GeoJSON or downloaded from OSM
    - Foursquare places dataset per city (Part 2 only):
      GeoJSON with fields: catName1, Activity_Type, checkins:2023

Citation:
    If you use this code, please cite:
    Nolasco-Cirugeda, A., López-Baeza, J., García-Pérez, S., & García-Mayor, C. (2026).
    Shaping socio-functional urban identity: When morphometrics meet social media data
    to reveal neighbourhood structure and cores. Cities.
    https://doi.org/10.1016/j.cities.2026.107373
    Code: https://doi.org/10.5281/zenodo.XXXXXXX
"""

# =============================================================================
# CONFIGURATION — edit these paths before running
# =============================================================================

INPUT_FOLDER = "input/"       # folder containing Case-studies.xlsx and GeoJSONs
RESULTS_FOLDER = "results/"   # folder where outputs will be saved
LOCAL_CRS = 25830             # projected CRS for your study area (EPSG code)
                              # e.g. 25830 for UTM zone 30N (Spain)
                              # e.g. 5514 for Czech Republic / Central Europe

# =============================================================================
# PART 1 — MORPHOLOGICAL CLUSTERING
# =============================================================================

import os
import warnings

import geopandas as gpd
import libpysal
import momepy
import osmnx as ox
import pandas as pd
from clustergram import Clustergram
import matplotlib.pyplot as plt
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union
import numpy as np
from math import log
from tqdm import tqdm
import scipy.stats as stats
from concurrent.futures import ProcessPoolExecutor

print('Libraries imported.')

places_list = pd.read_excel(os.path.join(INPUT_FOLDER, "Case-studies.xlsx"))

print('List of cities imported.')


def process_city_morphology(row):
    place = str(row['City'] + ', ' + row['Country'])
    k_number_of_clusters = int(row['Form Clusters'])
    print(f"{place} — starting morphological clustering for {k_number_of_clusters} clusters.")

    # INPUT FILES
    # Buildings: from local GeoJSON (cadastral or similar) or downloaded from OSM
    buildings_path = os.path.join(INPUT_FOLDER, f"{place} Input-buildings.geojson")
    if os.path.isfile(buildings_path):
        buildings = gpd.read_file(buildings_path).to_crs(LOCAL_CRS)
        if "REFCAT" in buildings.columns:
            buildings = buildings.dissolve(by="REFCAT").reset_index()
        buildings = buildings.explode().reset_index(drop=True)
        buildings = buildings[buildings.geom_type == "Polygon"].reset_index(drop=True)
        buildings = buildings[["geometry"]].to_crs(LOCAL_CRS)
        buildings["uID"] = range(len(buildings))
    else:
        buildings = ox.features_from_place(place, tags={'building': True})
        buildings = buildings.explode().reset_index(drop=True)
        buildings = buildings[buildings.geom_type == "Polygon"].reset_index(drop=True)
        buildings = buildings[["geometry"]].to_crs(LOCAL_CRS)
        buildings["uID"] = range(len(buildings))
        buildings.to_crs(4326).to_file(buildings_path, driver="GeoJSON")
        buildings = buildings.to_crs(LOCAL_CRS)
        print(f"{place} — buildings downloaded from OSM.")

    # Streets: from local GeoJSON or downloaded from OSM
    streets_path = os.path.join(INPUT_FOLDER, f"{place} Input-streets.geojson")
    if os.path.isfile(streets_path):
        streets = gpd.read_file(streets_path).to_crs(LOCAL_CRS)
        streets = momepy.remove_false_nodes(streets)
        streets = streets[["geometry"]]
        streets["nID"] = range(len(streets))
        streets = streets.to_crs(LOCAL_CRS)
    else:
        osm_graph = ox.graph_from_place(place, network_type='drive')
        osm_graph = ox.projection.project_graph(osm_graph, to_crs=LOCAL_CRS)
        streets = ox.graph_to_gdfs(osm_graph, nodes=False, edges=True,
                                   node_geometry=False, fill_edge_geometry=True)
        streets = momepy.remove_false_nodes(streets)
        streets = streets[["geometry"]]
        streets["nID"] = range(len(streets))
        streets.to_crs(4326).to_file(streets_path, driver="GeoJSON")
        streets = streets.to_crs(LOCAL_CRS)
        print(f"{place} — streets downloaded from OSM.")

    print('Streets and buildings imported.')

    # CALCULATE METRICS
    tess_path = os.path.join(RESULTS_FOLDER, f"{place} Results-tesselation.geojson")
    bld_path  = os.path.join(RESULTS_FOLDER, f"{place} Results-buildings.geojson")
    str_path  = os.path.join(RESULTS_FOLDER, f"{place} Results-streets.geojson")

    if os.path.isfile(tess_path) and os.path.isfile(bld_path) and os.path.isfile(str_path):
        tessellation = gpd.read_file(tess_path).to_crs(LOCAL_CRS)
        buildings    = gpd.read_file(bld_path).to_crs(LOCAL_CRS)
        streets      = gpd.read_file(str_path).to_crs(LOCAL_CRS)
        graph = momepy.gdf_to_nx(streets)
        graph = momepy.node_degree(graph)
        nodes, streets = momepy.nx_to_gdf(graph)
        queen_1 = libpysal.weights.contiguity.Queen.from_dataframe(
            tessellation, ids="uID", silence_warnings=True)
        queen_3 = momepy.sw_high(k=3, weights=queen_1)
    else:
        limit = momepy.buffered_limit(buildings, 300)
        tessellation = momepy.Tessellation(
            buildings, "uID", limit, verbose=False, segment=1).tessellation
        print('Tessellation created.')

        buildings = buildings.sjoin_nearest(streets, max_distance=1000, how="left")
        buildings = buildings.drop_duplicates("uID").drop(columns="index_right")
        tessellation = tessellation.merge(buildings[['uID', 'nID']], on='uID', how='left')
        print('Streets linked.')

        buildings["building-area"]    = buildings.area
        tessellation["tile-area"]     = tessellation.area
        streets["segment-length"]     = streets.length
        print('Dimensions measured.')

        buildings['building-rectangular-eri'] = momepy.EquivalentRectangularIndex(buildings).series
        buildings['building-elongation']      = momepy.Elongation(buildings).series
        tessellation['tile-convexity']        = momepy.Convexity(tessellation).series
        streets["street-linearity"]           = momepy.Linearity(streets).series
        print('Shape measured.')

        buildings["building-shared_walls"] = momepy.SharedWallsRatio(buildings).series
        queen_1 = libpysal.weights.contiguity.Queen.from_dataframe(
            tessellation, ids="uID", silence_warnings=True)
        tessellation["tile-covered_area"] = momepy.CoveredArea(
            tessellation, queen_1, "uID", verbose=False).series

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

        queen_3      = momepy.sw_high(k=3, weights=queen_1)
        buildings_q1 = libpysal.weights.contiguity.Queen.from_dataframe(
            buildings, silence_warnings=True)
        buildings['interbuilding_distance'] = momepy.MeanInterbuildingDistance(
            buildings, queen_1, 'uID', queen_3, verbose=False).series
        buildings['building-adjacency'] = momepy.BuildingAdjacency(
            buildings, queen_3, 'uID', buildings_q1, verbose=False).series

        profile = momepy.StreetProfile(streets, buildings)
        streets["street-width"]           = profile.w
        streets["street-width_deviation"] = profile.wd
        streets["street-openness"]        = profile.o
        print('Spatial distribution measured.')

        tessellation['tile-area-ratio'] = momepy.AreaRatio(
            tessellation, buildings, 'tile-area', 'building-area', 'uID').series
        print('Intensity measured.')

        graph = momepy.gdf_to_nx(streets)
        graph = momepy.node_degree(graph)
        graph = momepy.closeness_centrality(graph, radius=400, distance="mm_len")
        graph = momepy.meshedness(graph, radius=400, distance="mm_len")
        graph = momepy.betweenness_centrality(
            graph, name='betweenness_metric_n', mode='nodes', weight='mm_len')
        nodes, streets = momepy.nx_to_gdf(graph)
        print('Connectivity measured.')

        tessellation.to_crs(4326).to_file(tess_path, driver="GeoJSON")
        buildings.to_crs(4326).to_file(bld_path, driver="GeoJSON")
        streets.to_crs(4326).to_file(str_path, driver="GeoJSON")
        print(f"{place} — results exported.")

    print('Metrics calculated.')

    buildings["nodeID"] = momepy.get_node_id(buildings, nodes, streets, "nodeID", "nID")
    merged = tessellation.merge(buildings.drop(columns=['nID', 'geometry']), on='uID')
    merged = merged.merge(streets.drop(columns='geometry'), on='nID', how='left')
    merged = merged.merge(nodes.drop(columns='geometry'), on='nodeID', how='left')

    percentiles = []
    for column in merged.columns.drop(
            ["uID", "nodeID", "nID", 'mm_len', 'node_start', 'node_end', "geometry"]):
        perc = momepy.Percentiles(merged, column, queen_3, "uID", verbose=False).frame
        perc.columns = [f"{column}_" + str(x) for x in perc.columns]
        percentiles.append(perc)

    percentiles_joined = pd.concat(percentiles, axis=1)
    standardized = (percentiles_joined - percentiles_joined.mean()) / percentiles_joined.std()

    cgram = Clustergram(range(1, int(k_number_of_clusters + 1)), n_init=10, random_state=42)
    cgram.fit(standardized.fillna(0))

    fig, axs = plt.subplots(2, figsize=(7, 4), sharex=True)
    cgram.silhouette_score().plot(
        xlabel="Number of clusters (k)", ylabel="Silhouette score", ax=axs[0])
    cgram.calinski_harabasz_score().plot(
        xlabel="Number of clusters (k)", ylabel="Calinski-Harabasz score", ax=axs[1])

    plot_path = os.path.join(RESULTS_FOLDER, f"{place} Plot scores (form-only).png")
    plt.savefig(plot_path)
    print(f"{'*'*10}  {place} finished.  {'*'*10}")


print('Starting morphological clustering for all places.')

if __name__ == "__main__":
    with ProcessPoolExecutor() as executor:
        executor.map(process_city_morphology,
                     [row for _, row in places_list.iterrows()])

print('All places done — Part 1.')


# =============================================================================
# PART 2 — SOCIO-FUNCTIONAL CLUSTERING
# =============================================================================

import geopandas
import pandas

print('Libraries imported.')

places_list = pandas.read_excel(os.path.join(INPUT_FOLDER, "Case-studies.xlsx"))

print('List of cities imported.')

for index, row in places_list.iterrows():
    city = str(row['City'])
    place = str(row['City'] + ', ' + row['Country'])
    k_number_of_clusters = int(row['FS Clusters'])
    print(f"{place} — starting socio-functional clustering for {k_number_of_clusters} clusters.")

    # INPUT FILES
    # Buildings
    buildings_path = os.path.join(INPUT_FOLDER, f"{place} Input-buildings.geojson")
    if os.path.isfile(buildings_path):
        buildings = geopandas.read_file(buildings_path).to_crs(LOCAL_CRS)
        if "REFCAT" in buildings.columns:
            buildings = buildings.dissolve(by="REFCAT").reset_index()
        buildings = buildings.explode().reset_index(drop=True)
        buildings = buildings[buildings.geom_type == "Polygon"].reset_index(drop=True)
        buildings = buildings[["geometry"]].to_crs(LOCAL_CRS)
        buildings["uID"] = range(len(buildings))
    else:
        buildings = osmnx.features_from_place(place, tags={'building': True})
        buildings = buildings.explode().reset_index(drop=True)
        buildings = buildings[buildings.geom_type == "Polygon"].reset_index(drop=True)
        buildings = buildings[["geometry"]].to_crs(LOCAL_CRS)
        buildings["uID"] = range(len(buildings))
        buildings.to_crs(4326).to_file(buildings_path, driver="GeoJSON")
        buildings = buildings.to_crs(LOCAL_CRS)
        print(f"{place} — buildings downloaded from OSM.")

    # Streets
    streets_path = os.path.join(INPUT_FOLDER, f"{place} Input-streets.geojson")
    if os.path.isfile(streets_path):
        streets = geopandas.read_file(streets_path).to_crs(LOCAL_CRS)
        streets = momepy.remove_false_nodes(streets)
        streets = streets[["geometry"]]
        streets["nID"] = range(len(streets))
        streets = streets.to_crs(LOCAL_CRS)
    else:
        osm_graph = osmnx.graph_from_place(place, network_type='drive')
        osm_graph = osmnx.projection.project_graph(osm_graph, to_crs=LOCAL_CRS)
        streets = osmnx.graph_to_gdfs(osm_graph, nodes=False, edges=True,
                                      node_geometry=False, fill_edge_geometry=True)
        streets = momepy.remove_false_nodes(streets)
        streets = streets[["geometry"]]
        streets["nID"] = range(len(streets))
        streets.to_crs(4326).to_file(streets_path, driver="GeoJSON")
        streets = streets.to_crs(LOCAL_CRS)
        print(f"{place} — streets downloaded from OSM.")

    print('Streets and buildings imported.')

    # CALCULATE METRICS
    tess_path = os.path.join(RESULTS_FOLDER, f"{place} Results-tesselation.geojson")
    bld_path  = os.path.join(RESULTS_FOLDER, f"{place} Results-buildings.geojson")
    str_path  = os.path.join(RESULTS_FOLDER, f"{place} Results-streets.geojson")

    if os.path.isfile(tess_path) and os.path.isfile(bld_path) and os.path.isfile(str_path):
        tessellation = geopandas.read_file(tess_path).to_crs(LOCAL_CRS)
        buildings    = geopandas.read_file(bld_path).to_crs(LOCAL_CRS)
        streets      = geopandas.read_file(str_path).to_crs(LOCAL_CRS)
        graph = momepy.gdf_to_nx(streets)
        graph = momepy.node_degree(graph)
        nodes, streets = momepy.nx_to_gdf(graph)
        queen_1 = libpysal.weights.contiguity.Queen.from_dataframe(
            tessellation, ids="uID", silence_warnings=True)
        queen_3 = momepy.sw_high(k=3, weights=queen_1)
    else:
        limit = momepy.buffered_limit(buildings, 300)
        tessellation = momepy.Tessellation(
            buildings, "uID", limit, verbose=False, segment=1)
        tessellation = tessellation.tessellation
        print('Tessellation created.')

        buildings = buildings.sjoin_nearest(streets, max_distance=1000, how="left")
        buildings = buildings.drop_duplicates("uID").drop(columns="index_right")
        tessellation = tessellation.merge(buildings[['uID', 'nID']], on='uID', how='left')
        print('Streets linked.')

        buildings["building-area"]    = buildings.area
        tessellation["tile-area"]     = tessellation.area
        streets["segment-length"]     = streets.length
        print('Dimensions measured.')

        buildings['building-rectangular-eri'] = momepy.EquivalentRectangularIndex(buildings).series
        buildings['building-elongation']      = momepy.Elongation(buildings).series
        tessellation['tile-convexity']        = momepy.Convexity(tessellation).series
        streets["street-linearity"]           = momepy.Linearity(streets).series
        print('Shape measured.')

        buildings["building-shared_walls"] = momepy.SharedWallsRatio(buildings).series
        queen_1 = libpysal.weights.contiguity.Queen.from_dataframe(
            tessellation, ids="uID", silence_warnings=True)
        tessellation["tile-covered_area"] = momepy.CoveredArea(
            tessellation, queen_1, "uID", verbose=False).series

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

        queen_3      = momepy.sw_high(k=3, weights=queen_1)
        buildings_q1 = libpysal.weights.contiguity.Queen.from_dataframe(
            buildings, silence_warnings=True)
        buildings['interbuilding_distance'] = momepy.MeanInterbuildingDistance(
            buildings, queen_1, 'uID', queen_3, verbose=False).series
        buildings['building-adjacency'] = momepy.BuildingAdjacency(
            buildings, queen_3, 'uID', buildings_q1, verbose=False).series

        profile = momepy.StreetProfile(streets, buildings)
        streets["street-width"]           = profile.w
        streets["street-width_deviation"] = profile.wd
        streets["street-openness"]        = profile.o
        print('Spatial distribution measured.')

        tessellation['tile-area-ratio'] = momepy.AreaRatio(
            tessellation, buildings, 'tile-area', 'building-area', 'uID').series
        print('Intensity measured.')

        graph = momepy.gdf_to_nx(streets)
        graph = momepy.node_degree(graph)
        graph = momepy.closeness_centrality(graph, radius=400, distance="mm_len")
        graph = momepy.meshedness(graph, radius=400, distance="mm_len")
        graph = momepy.betweenness_centrality(
            graph, name='betweenness_metric_n', mode='nodes', weight='mm_len')
        nodes, streets = momepy.nx_to_gdf(graph)
        print('Connectivity measured.')

        buildings.to_crs(4326).to_file(bld_path, driver="GeoJSON")
        streets.to_crs(4326).to_file(str_path, driver="GeoJSON")
        print(f"{place} — results exported.")

    print('Metrics calculated.')

    # IMPORT FOURSQUARE DATA
    fs_tess_path = os.path.join(RESULTS_FOLDER, f"{place} Results-FS-tesselation.geojson")

    if os.path.isfile(fs_tess_path):
        tessellation = geopandas.read_file(fs_tess_path).to_crs(LOCAL_CRS)
    else:
        points_path = os.path.join(INPUT_FOLDER, f"{city}_Points.geojson")
        points = geopandas.read_file(points_path).to_crs(4326)
        grid   = tessellation.to_crs(4326)
        print('Foursquare files imported.')

        points = geopandas.clip(points, grid)
        join   = geopandas.sjoin(points, grid, how="inner", op="within")
        print(f"Points clipped to {place}.")

        counts = join.groupby(["index_right"]).size()
        grid["CNT_points_total"] = grid.index.map(counts).fillna(0)

        checkin_sums = join.groupby(["index_right"])["checkins:2023"].sum()
        grid["SUM_checkins_2023"] = grid.index.map(checkin_sums).fillna(0)

        plcats = join.groupby(["index_right"])["catName1"].nunique()
        grid["CNTU_catID"] = grid.index.map(plcats).fillna(0)

        for act_type, col_prefix in [
            ('Consumption', 'Consumption'),
            ('Social', 'Social'),
            ('Income', 'Income'),
            ('Leisure', 'Leisure'),
        ]:
            cnt = join[join['Activity_Type'] == act_type].groupby(['index_right']).size()
            grid[f'CNT_{col_prefix}'] = grid.index.map(cnt).fillna(0)
            sm  = join[join['Activity_Type'] == act_type].groupby(
                ['index_right'])["checkins:2023"].sum()
            grid[f'SUM_{col_prefix}'] = grid.index.map(sm).fillna(0)

        print('General Foursquare statistics calculated.')

        diversity_results = {}
        diversity_indices = [
            'simpsons', 'shannon', 'shannon_wiener',
            'invsimpson', 'gini_simpson', 'renyi_entropy'
        ]

        for idx in tqdm(diversity_indices):
            counts_div = join.groupby(
                ['index_right', 'catName1'])["index_right"].size().reset_index(name='count')
            table = counts_div.pivot(
                index='index_right', columns='catName1', values='count').fillna(0)
            proportions = table.div(table.sum(axis=1), axis=0)

            if idx == 'simpsons':
                diversity = 1 - (proportions ** 2).sum(axis=1)
            elif idx == 'shannon':
                diversity = -(proportions * np.log(proportions)).sum(axis=1)
            elif idx == 'shannon_wiener':
                shannon   = -(proportions * np.log2(proportions)).sum(axis=1)
                diversity = 2 ** shannon - 1
            elif idx == 'invsimpson':
                diversity = 1 / (proportions ** 2).sum(axis=1)
            elif idx == 'gini_simpson':
                sq = (proportions ** 2).sum(axis=1)
                diversity = 1 - (sq / (sq + ((1 - proportions) ** 2).sum(axis=1)))
            elif idx == 'renyi_entropy':
                alpha     = 2
                diversity = (1 / (1 - alpha)) * np.log((proportions ** alpha).sum(axis=1))

            diversity_results[idx] = diversity

        results      = pandas.DataFrame(diversity_results)
        tessellation = pandas.concat([grid, results], axis=1).to_crs(LOCAL_CRS)
        tessellation.to_crs(4326).to_file(fs_tess_path, driver="GeoJSON")

    print('Foursquare diversity statistics calculated.')

    # CALCULATE SOCIO-FUNCTIONAL TYPES
    buildings["nodeID"] = momepy.get_node_id(buildings, nodes, streets, "nodeID", "nID")
    merged = tessellation.merge(buildings.drop(columns=['nID', 'geometry']), on='uID')
    merged = merged.merge(streets.drop(columns='geometry'), on='nID', how='left')
    merged = merged.merge(nodes.drop(columns='geometry'), on='nodeID', how='left')

    percentiles = []
    for column in merged.columns.drop(
            ["uID", "nodeID", "nID", 'mm_len', 'node_start', 'node_end', "geometry"]):
        perc = momepy.Percentiles(merged, column, queen_3, "uID", verbose=False).frame
        perc.columns = [f"{column}_" + str(x) for x in perc.columns]
        percentiles.append(perc)

    percentiles_joined = pandas.concat(percentiles, axis=1)
    standardized = (percentiles_joined - percentiles_joined.mean()) / percentiles_joined.std()

    cgram = Clustergram(
        range(1, int(k_number_of_clusters + 1)), n_init=10, random_state=42)
    cgram.fit(standardized.fillna(0))

    fig, axs = plt.subplots(2, figsize=(7, 4), sharex=True)
    cgram.silhouette_score().plot(
        xlabel="Number of clusters (k)", ylabel="Silhouette score", ax=axs[0])
    cgram.calinski_harabasz_score().plot(
        xlabel="Number of clusters (k)", ylabel="Calinski-Harabasz score", ax=axs[1])

    plot_path = os.path.join(RESULTS_FOLDER, f"{place} Plot scores (foursquare-full).png")
    plt.savefig(plot_path)
    print(f"{'*'*10}  {place} finished.  {'*'*10}")

print('All places done — Part 2.')
