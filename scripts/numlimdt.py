# -*- coding: utf-8 -*-
"""
Created on Wed May  5 11:22:17 2021

@author: danie
"""

import netCDF4 as nc
import geopandas as gpd
from shapely.geometry import Point

fn = r"..\modellen\boezemmodel\fm\DFM_OUTPUT_boezemmodel\boezemmodel_map.nc"
ds = nc.Dataset(fn)

# inlezen waterhoogtes uit NetCDF
numlimdt_gdf = gpd.GeoDataFrame(
    data={"numlimdt": ds["mesh1d_Numlimdt"][:].max(axis=0),
          "geometry": [
              Point(coords) for coords in zip(ds["mesh1d_node_x"], ds["mesh1d_node_y"])
              ]}
    )

numlimdt_gdf.crs = "epsg:28992"
numlimdt_gdf.to_file("numlimdt.shp")
