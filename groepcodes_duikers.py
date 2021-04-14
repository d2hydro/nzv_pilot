# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:05:46 2021

@author: danie
"""
import geopandas as gpd
import pandas as pd
from pathlib import Path

import hydrotools

data_path = Path(r".\data").absolute().resolve()
beheerregister = data_path.joinpath("beheerregister")
excelbestanden = data_path.joinpath("xlsx")

invoerbestanden = {
    "modelgebied": "Bemalings_en_lozingsgebieden.shp",
    "branches": "Afaanvoervakken_DHYDRO_v5.shp",
    "profielpunten": "profielpunten_3feb2021_actueel.shp",
    "bruggen": "Bruggen_(doorstroomopeningen).shp",
    "duikers": "duikers.shp",
    "sluizen": "sluizen.shp",
    "sifons": "sifons.shp",
    "stuwen": "stuwen.shp",
    "inlaten": "inlaten.shp",
    "gemalen": "gemalen_v2.shp",
    "peilgebieden": "Peilgebieden.shp"
}

principe_profielen_df = pd.read_csv(
    excelbestanden.joinpath("principe_profielen.csv"),
    index_col=0)

hydamo = hydrotools.load_model(Path(r"./hydamo_model/boezemmodel_v4.pickle"))

all_culverts_gdf = gpd.read_file(beheerregister.joinpath(invoerbestanden["duikers"]))
groups = {"code":[],
          "object_laag":[],
          "group_code": []}
group_id = 1

for code, row in hydamo.culverts.iterrows():
    centroid = row["geometry"].centroid
    bodembreedte = min(principe_profielen_df.loc[row["branch_id"]]["bottomwidth"],
                       25)
    buffer = row["geometry"].centroid.buffer(bodembreedte)
    gdf = all_culverts_gdf.loc[all_culverts_gdf["geometry"].centroid.within(buffer)]
    if len(gdf) > 1:
        groups["code"] += gdf["CODE"].to_list()
        groups["object_laag"] += ["culverts"] * len(gdf)
        groups["group_code"] += [f"kw_groep_{group_id:03d}"] * len(gdf)
        group_id += 1

groups_df = pd.DataFrame(data=groups)

groups_df.to_csv(excelbestanden.joinpath("duiker_groepen.csv"), index=False)
