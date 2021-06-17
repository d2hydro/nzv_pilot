"""
Created on Mon Sep 21 23:18:13 2020

@author: danie
"""

import geopandas as gpd
import pandas as pd
from delft3dfmpy import HyDAMO
from delft3dfmpy.converters.hydamo_to_dflowfm import roughness_gml
from delft3dfmpy.core.geometry import find_nearest_branch
from pathlib import Path
import re
import numpy as np
import csv

from shapely.geometry import LineString, Point
from shapely.ops import snap
import pickle

#%%
hydamo = HyDAMO()

ATTRIBUTES = ["crosssections",
              "bridges",
              "culverts",
              "orifices",
              "weirs",
              "gemalen",
              "pumps"]

def _valid_pprof(prof_def):
    prof_def["slope"] = max(0.1, prof_def["slope"])
    prof_def["bottomwidth"] = max(0.5, prof_def["bottomwidth"])
    prof_def["maximumflowwidth"] = max(prof_def["bottomwidth"],
                                       prof_def["maximumflowwidth"])
    
    return(prof_def)

def _make_list(item):
    if not isinstance(item, list):
        item = [item]
    return item


def _filter(gdf, attribute_filter):
    if isinstance(attribute_filter, dict):
        for key, value in attribute_filter.items():
            value = _make_list(value)
            gdf = gdf[gdf[key].isin(value)]

        return gdf
    else:
       raise IOError('attribute_filter should be dictionary') 

def _cut_line(line, distance):
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]
        
def _get_nodes(gdf):
    coords = gdf["geometry"].apply(lambda x: x.coords[0]).to_list()
    coords += gdf["geometry"].apply(lambda x: x.coords[-1]).to_list()
    nodes = gpd.GeoSeries([Point(i) for i in list(set(coords))])
    return nodes


def read_file(path,
              hydamo_attribute,
              index_col=None,
              keep_indices=None,
              attribute_filter=None,
              snap_to_branches=None,
              keep_columns=None,
              column_mapping=None,
              z_coord=False
              ):
        """
        Read any OGR supported feature-file to match hydamo-property.

        A mask file can be specified to clip the selection.

        Parameters
        ----------
        path : str or Path
            Path to the feature file
        hydamo_attribute : HyDAMO property
            property to map to (HyDAMO.branches, HyDAMO.crosssections)
        index_col: str
            column to be used to identify rows to keep
        keep_indices: List[str]
            values in index_col to keep after filtering & snapping
        attribute_filter: dict
            dict with lists or strings of the format {'column_name': [values to keep]}
        snap_to_branches: dict
            snap to a geodataframe with LineStrings. Keep only geometries that snap
            to a LineString with a certain attribute value.
            Format: {'branches': GeoDataFrame,
                     attribute_filter: {'column_name': [values to keep]}}
        column_mapping: dict
            dict for renaming input colunns to required columns
        Result: GeoDataFrame matching the HyDAMO property
        """
        gdf = gpd.read_file(path)
        gdf.columns = gdf.columns.str.lower()

        index_col = next(
                (k.lower() for k, v in column_mapping.items() if v.lower() == "code"),
                None)

        if index_col is None:
            index_col = "code"
            
            

        # filter by gdf_snap
        if snap_to_branches:
            distance = snap_to_branches["distance"]
            branches = snap_to_branches["branches"]
            find_nearest_branch(branches, gdf, maxdist=distance)
            gdf = gdf.loc[gdf["branch_offset"].notna()]
            gdf.loc[:, "hydromodel"] = gdf.apply(
                (lambda x: branches.loc[x["branch_id"]]["hydromodel"]),
                axis=1)

            if snap_to_branches["attribute_filter"]:
                snap_attribute_filter = {
                    key.lower(): value for key, value in snap_to_branches[
                        "attribute_filter"].items()}
                gdf = _filter(gdf, snap_attribute_filter)

        # filter by attribute
        if attribute_filter:
            attribute_filter = {
                key.lower(): value for key, value in attribute_filter.items()}
            gdf = _filter(gdf,
                          attribute_filter)

        # map to hydamo columns
        if column_mapping:
            column_mapping = {
                key.lower(): value.lower() for key, value in column_mapping.items()
                }
            gdf.rename(columns=column_mapping, inplace=True)

        # drop all columns not needed
        required_columns = getattr(hydamo, hydamo_attribute).required_columns.copy()

        if keep_columns:
            required_columns += [col.lower() for col in keep_columns]

        if hydamo_attribute == 'crosssections':
            if z_coord:
                required_columns += ['z', 'order']

        drop_cols = [
            col for col in gdf.columns if col not in required_columns + ['geometry']
            ]
        if len(drop_cols) > 0:
            gdf = gdf.drop(drop_cols, axis=1)

        return gdf


def to_file(model, hydamo_attribute, length=False, path=Path('.')):
    """Convert hydamo class to shape-file."""
    path = Path(path)
    hydamo_class = getattr(model, hydamo_attribute)
    if not hydamo_class.empty:
        data = {col: hydamo_class[col].values for col in hydamo_class.columns}

        if length:
            data = data = {**data, 'length': hydamo_class['geometry'].length.values}

        gpd.GeoDataFrame(data=data).to_file(path.joinpath(f'{hydamo_attribute}.shp'))


def snap_ends(gdf, tolerance, digits=None):
    """Snap all end-vertices within a specified tolerance."""
    sindex = gdf.sindex
    snapped = []
    for index, row in gdf.iterrows():
        # rough selection on index
        buffer_geom = row['geometry'].buffer(tolerance)
        # precise selection on distance < tolerance
        gdf_selec = gdf.iloc[list(sindex.intersection(buffer_geom.bounds))].copy()
        gdf_selec['distance'] = gdf_selec.distance(row['geometry'])
        gdf_selec = gdf_selec.loc[gdf_selec['distance'] < tolerance]
        # only snap to features that will not be modified
        gdf_selec = gdf_selec.loc[gdf_selec.index.isin(snapped)]
        # snapping to remaining objects
        geom = row['geometry']
        # round digits (optionally)
        if digits:
            geom = LineString([[round(coord, ndigits=digits) for
                               coord in coords] for coords in geom.coords])

        if not gdf_selec.empty:
            geom_coords = list(geom.coords)

            for _, row_selec in gdf_selec.iterrows():
                for dst_vert in [0, -1]:
                    for src_vert in [0, -1]:
                        geom_coords[dst_vert] = snap(
                            Point(geom_coords[dst_vert]),
                            Point(row_selec['geometry'].coords[src_vert]),
                            tolerance=1).coords[0]
            geom = LineString(geom_coords)
        # write feature in original GeoDataFrame
        gdf.loc[index, 'geometry'] = geom
        # mark index as snapped
        snapped += [index]

    return gdf


def filter_model(model, attribute_filter=None, geometry=None):
    """Filter a hydamo model on an attribute filter on branches."""
    drop_branches = []

    if attribute_filter:
        attribute_filter = {
            key.lower(): value for key, value in attribute_filter.items()}
        for key, value in attribute_filter.items():
            drop_branches += list(
                model.branches.loc[
                    model.branches[key] != value].index)

    if geometry:
        drop_branches += list(model.branches.loc[
            ~model.branches.intersects(geometry)].index)

    drop_branches = list(set(drop_branches))

    model.branches = model.branches.loc[~model.branches.index.isin(drop_branches)]
    for attribute in ATTRIBUTES:
        hydamo_class = getattr(model, attribute)
        if 'branch_id' in hydamo_class.columns:
            hydamo_class.set_data(hydamo_class.loc[
                ~hydamo_class['branch_id'].isin(drop_branches)],
                index_col="code",
                check_columns=True,
                check_geotype=True)

    return model


def export_shapes(model, path=Path('.')):
    """Export a hydamo class to shape-files."""
    path = Path(path)
    path.mkdir(exist_ok=True)
    for attribute in ATTRIBUTES:
        to_file(model, attribute, length=False, path=path)

    to_file(model, "branches", length=True, path=path)


def save_model(model, file_name=Path('model.pickle')):
    """Save the model as a pickle."""
    file_name = Path(file_name)
    parent = file_name.parent
    parent.mkdir(exist_ok=True)
    with open(file_name, 'wb') as dst:
        pickle.dump(model, dst, protocol=pickle.HIGHEST_PROTOCOL)


def load_model(file_name):
    """Load the model from a pickle."""
    with open(file_name, 'rb') as src:
        model = pickle.load(src)
    return model


def get_trapeziums(gdf,
                   index,
                   bottom_width,
                   bottom_level,
                   waterlevel_width,
                   slope_left,
                   slope_right,
                   roughnesstype,
                   roughnessvalue):
    """Return trapezium profiles for branches."""
    gdf = gdf.set_index(index)
    definitions = {}
    for idx, row in gdf.iterrows():
        slope = (row[slope_left] + row[slope_right]) / 2
        maximumflowwidth = row[waterlevel_width] + (2 * slope)
        bottomwidth = row[bottom_width]
        bottomlevel = row[bottom_level]
        definitions[idx] = dict(slope=slope,
                                bottomwidth=bottomwidth,
                                bottomlevel=bottomlevel,
                                maximumflowwidth=maximumflowwidth,
                                roughnesstype=row[roughnesstype],
                                roughnessvalue=row[roughnessvalue]
                                )
    return pd.DataFrame.from_dict(definitions, orient="index")


def add_trapeziums(dfmmodel, principe_profielen_df, closed=False):
    """Add trapezium profiles on branches with missing crosssections."""
    xs = dfmmodel.crosssections
    for branch in xs.get_branches_without_crosssection():
        if branch in principe_profielen_df.index:
            prof_def = _valid_pprof(dict(principe_profielen_df.loc[branch]))
            chainage = dfmmodel.network.branches.loc[branch]['geometry'].length / 2
            definition = f"PPRO_{branch}"
            xs.add_crosssection_location(branch,
                                         chainage,
                                         definition,
                                         shift=prof_def["bottomlevel"]
                                         )

            xs.add_trapezium_definition(
                name=definition,
                slope=prof_def["slope"],
                maximumflowwidth=prof_def["maximumflowwidth"],
                bottomwidth=prof_def["bottomwidth"],
                closed=closed,
                roughnesstype=roughness_gml[int(prof_def["roughnesstype"])],
                roughnessvalue=float(prof_def["roughnessvalue"]))

    return dfmmodel


def filter_to_other_object(row, object_gdf, max_distance):
    """Filter HyDAMO-class-objects within distance to another object-class."""
    gdf = object_gdf.loc[
        object_gdf["geometry"].centroid.distance(row["geometry"]) < max_distance
        ]

    if not gdf.empty:
        gdf = gdf.loc[gdf["branch_id"] == row["branch_id"]]

    return gdf.empty

def move_end_nodes(branches_gdf, move_lines_gdf, threshold):
    #%% add start & end node of linestrings
    branches_gdf.loc[:, "start_node"] = branches_gdf["geometry"].apply(
        lambda x: Point(x.coords[0])
        )
    branches_gdf.loc[:, "end_node"] = branches_gdf["geometry"].apply(
        lambda x: Point(x.coords[-1])
        )
    
    #%% add start & end node of linestrings
    modified_rows = []
    for _, row in move_lines_gdf.iterrows():
        from_node = row["geometry"].coords[0]
        to_node = row["geometry"].coords[-1]
        from_poly = Point(from_node).buffer(threshold)
        
        # add to_node at beginning of LineString when start_node intersects from_node
        rows_select = branches_gdf[
            branches_gdf["start_node"].within(from_poly)
            ].index.to_list()
        branches_gdf.loc[rows_select, "geometry"] = branches_gdf.loc[
            rows_select, "geometry"
            ].apply(lambda x: LineString([to_node] + list(x.coords)))
        modified_rows += rows_select
    
        # extend LineString with to_node when start_node intersects to_node
        rows_select = branches_gdf[
            branches_gdf["end_node"].within(from_poly)
            ].index.to_list()
        branches_gdf.loc[rows_select, "geometry"] = branches_gdf.loc[
            rows_select, "geometry"
            ].apply(lambda x: LineString(list(x.coords) + [to_node]))
        modified_rows += rows_select
        
    #%% remove all lines with a length < treshold between new startand end_nodes
    branches_gdf.loc[modified_rows, "start_end_dist"] = branches_gdf.loc[
        modified_rows, "geometry"].apply(lambda x: Point(x.coords[0]).distance(Point(x.coords[-1])))
    
    branches_gdf = branches_gdf.loc[
        (branches_gdf["start_end_dist"] > threshold) | (branches_gdf["start_end_dist"].isna())]
    
    branches_gdf.drop(["start_node", "end_node", "start_end_dist"], axis=1, inplace=True)
    
    return branches_gdf


class Sobek(object):
    """Sobek project-class."""    
    def __init__(self, path):
        path = Path(path)
        self.path = path
        self.cases = self._read_sbk_cases(path.joinpath("caselist.cmt"))
        
    def _read_sbk_cases(self, path):
        cases = {}
        with open(path, "r") as src:
            for line in src.readlines():
                line_split = re.split(r" '", line.replace("'\n", ""))
                cases[line_split[1]] = line_split[0]
        return cases

    def list_cases(self):
        """Return a list of cases."""
        return list(self.cases.keys())  
    
    def read_network(self, case):
        netw_path = self.path.joinpath(self.cases[case],"NETWORK.NTW")
        rows = []
        with open(netw_path) as src:
            for idx, line in enumerate(src):
                if "*" in line:
                    break
                elif idx > 0:
                    rows.append(line)
        df = pd.DataFrame([
            i for i in csv.reader(rows, skipinitialspace=True)
            ])
        return df

    def read_branches(self, case, pattern=".*", code_col="code"):
        """Read network-branches from Sobek case."""

        branches = {}
        nodes = {}

        # read network nodes and links
        netw_tp_path = self.path.joinpath(self.cases[case], "NETWORK.TP")
        
        #read branches to populate dictionary
        branches_raw = re.finditer(r"BRCH([\S\s]*?)brch",
                                   netw_tp_path.read_text(),
                                   flags=re.MULTILINE)

        for branch_raw in branches_raw:
            branch = branch_raw.group(0)
            branch_id = re.search(r"id '(\S*)'", branch).group(1)
            if re.match(pattern, branch_id):
                from_node = re.search(r"bn '(\S*)'", branch).group(1)
                to_node = re.search(r"en '(\S*)'", branch).group(1)
                branches[branch_id] = {"from_node": from_node,
                                       "to_node": to_node}

        # read nodes for start and end-coordinates
        nodes_raw = re.finditer(r"NODE([\S\s]*?)node",
                                netw_tp_path.read_text(),
                                flags=re.MULTILINE)

        for node_raw in nodes_raw:
            node = node_raw.group(0)
            node_id = re.search(r"id '(\S*)'", node).group(1)
            node_x = float(re.search(r"px ([+-]?([0-9]*[.])?[0-9]+)", node).group(1))
            node_y = float(re.search(r"py ([+-]?([0-9]*[.])?[0-9]+)", node).group(1))
            nodes[node_id] = Point(node_x, node_y)
        
        # read cp to add the line vertices
        netw_cp_path = self.path.joinpath(self.cases[case], "NETWORK.CP")
        branches_raw = re.finditer(r"BRCH([\S\s]*?)brch",
                                   netw_cp_path.read_text(),
                                   flags=re.MULTILINE)
        

        for branch_raw in branches_raw:
            branch = branch_raw.group(0)
            branch_id = re.search(r"id '(\S*)'", branch).group(1)
            if branch_id in branches.keys():
                table_str = re.findall(r"TBLE\n([\S\s]*?)\ntble", branch)[0]
                table_list = table_str.split("<\n")
                
                sumDistance = 0.0
                coord_list = [nodes[branches[branch_id]["from_node"]].coords[0]]
                for cp in table_list:
                    distance, angle = cp.split()[0:2]
                    distance = (float(distance) - sumDistance) * 2
                    angle = np.deg2rad(90 - float(angle))
                    x = coord_list[-1][0] + float(distance) * np.cos(angle)
                    y = coord_list[-1][1] + float(distance) * np.sin(angle)
                    coord_list += [(x, y)]
                    sumDistance += distance
                coord_list[-1] = nodes[branches[branch_id]["to_node"]].coords[0]
                branches[branch_id]["geometry"] = LineString(coord_list)

        gdf = gpd.GeoDataFrame.from_dict(branches, orient='index')
        gdf[code_col] = gdf.index

        return gdf

    def read_profiles(self, case, roughnesstype, roughnessvalue, pattern=".*"):
        """Read crosssections from Sobek case."""
        crosssections = {}
        def_ids = []

        # read dat-file as iterator
        dat_path = self.path.joinpath(self.cases[case], "Profile.dat")
        dats_raw = re.finditer(r"CRSN([\S\s]*?)crsn",
                               dat_path.read_text(),
                               flags=re.MULTILINE)

        # populate a list of profile locations
        for prof_dat_raw in dats_raw:
            prof_dat = prof_dat_raw.group(0)
            dat_id = re.search(r"id '(\S*)'", prof_dat).group(1)
            def_id = re.search(r"di '(\S*)'", prof_dat).group(1)
            if re.match(pattern, dat_id):
                crosssections[dat_id] = {"definition": def_id}
                def_ids += [def_id]

        # read network.cr file as iterator
        netw_path = self.path.joinpath(self.cases[case], "NETWORK.CR")
        netw_crs_raw = re.finditer(r"CRSN([\S\s]*?)crsn",
                                   netw_path.read_text(),
                                   flags=re.MULTILINE)

        for crs_raw in netw_crs_raw:
            crs_def = crs_raw.group(0)
            dat_id = re.search(r"id '(\S*)'", crs_def).group(1)
            if dat_id in crosssections.keys():
                branch_id = re.search(r"ci '(\S*)'", crs_def).group(1)
                chainage = float(re.search(r"lc ([+-]?([0-9]*[.])?[0-9]+)", crs_def).group(1))
                crosssections[dat_id].update({"branch": branch_id,
                                              "chainage": chainage})


        # read def-file as iterator
        def_path = self.path.joinpath(self.cases[case], "Profile.def")
        defs_raw = re.finditer(r"CRDS([\S\s]*?)crds",
                               def_path.read_text(),
                               flags=re.MULTILINE)

        # add profile definition
        for prof_def_raw in defs_raw:
            prof_def = prof_def_raw.group(0)

            def_id = re.search(r"id '(\S*)'", prof_def).group(1)
            if def_id in def_ids:
                dat_id = next(
                    k for k, v in crosssections.items() if v["definition"] == def_id
                    )
                def_type = int(re.search("ty ([0-9])", prof_def).group(1))
                if def_type in [0, 10]:
                    table_str = re.findall(r"TBLE\n([\S\s]*?)\ntble", prof_def)[0]
                    table_list = table_str.split("<\n")
                    if def_type == 0:
                        z = np.array([float(i.split()[0]) for i in table_list])
                        z = np.concatenate([np.flip(z), z])
                        w = np.array([float(i.split()[1]) for i in table_list])
                        y = np.concatenate([np.flip(-w / 2), w / 2])
                        yz = yz_fixer(np.array([list(i) for i in zip(y, z)]))
                        thalweg = 0
                    elif def_type == 10:
                        yz = yz_fixer(np.array(
                            [[float(i.split(" ")[0]), float(i.split(" ")[1])] for i in table_list]
                            ))
                        thalweg = (yz[-1][0] - yz[0][0])/2
                crosssections[dat_id].update({"yz": yz,
                                              "thalweg": thalweg,
                                              "roughnesstype": roughnesstype,
                                              "roughnessvalue": roughnessvalue})
        return crosssections
  
    
    def read_rr_laterals(self, case):
        """Read 3B boundaries from Sobek case."""
        bounds = {}

        # get levels
        bound_path = self.path.joinpath(self.cases[case], "BOUND3B.3B")
        bounds_raw = re.finditer(r"BOUN([\S\s]*?)boun",
                               bound_path.read_text(),
                               flags=re.MULTILINE)
        
        for bound_raw in bounds_raw:
            bound = bound_raw.group(0)
            bound_id = re.search(r"id '(\S*)'", bound).group(1)
            bound_level = re.search(r"bl 0 ([+-]?([0-9]*[.])?[0-9]+)", bound).group(1)
            bounds[bound_id] = {"level":bound_level}
            
        # get chainage
        netw_path = self.path.joinpath(self.cases[case], "NETWORK.CN")
        netw_cns_raw = re.finditer(r"FLBX([\S\s]*?)flbx",
                                  netw_path.read_text(),
                                  flags=re.MULTILINE)
        
        for netw_cn_raw in netw_cns_raw:
            netw_cn = netw_cn_raw.group(0)
            bound_id = re.search(r"id '(\S*)'", netw_cn).group(1)
            branch_id = re.search(r"ci '(\S*)'", netw_cn).group(1)
            chainage = float(re.search(r"lc ([+-]?([0-9]*[.])?[0-9]+)", netw_cn).group(1))
            bounds[bound_id].update({"branch": branch_id,
                                     "chainage": chainage})
        
        # get xy
        df = self.read_network(case)
        start_nodes_df = df[df[19] == "SBK_SBK-3B-REACH"][[14,21,22]].copy()
        start_nodes_df.columns = ["node_id", "x", "y"]
        end_nodes_df = df[df[32] == "SBK_SBK-3B-REACH"][[27,34,35]].copy()
        end_nodes_df.columns = ["node_id", "x", "y"]
        
        df = end_nodes_df.append(start_nodes_df).drop_duplicates(subset='node_id', keep='first')
        df.loc[:, "x"] = df["x"].astype(float)
        df.loc[:, "y"] = df["y"].astype(float)
        bounds_df = df.set_index("node_id").to_dict(orient="index")
        for k,v in bounds_df.items():
            bounds[k].update(v)
            
        return bounds
    
    def copy_rr(self, case, target_dir, fnm_path):
        # get file_list
        target_dir = Path(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)
        files = []
        with open(fnm_path) as src:
            for line in src:
                if not line.replace(" ","")[0] == "*":
                    if "*" in line:
                        files.append(line[0:line.find("*")].replace(" ","").replace("'",""))
                    else:
                        files.append(line).replace(" ","").replace("'","")
                        
        # now get that files copied
        for file in files:
            file_path = self.path.joinpath(self.cases[case], file)
            to_path = Path(target_dir).joinpath(file)
            if file_path.exists():
                file_path.copy(to_path)
                
        Path(fnm_path).copy(Path(target_dir).joinpath(fnm_path.name))
        
        
def test_sbk_profiles():
    sobek = Sobek(r"c:\SK215003\TKI3_NZV.lit")
    case = "Boezemmodel 0D1D"
    pattern = "pStor.*"
    roughnesstype="Manning"
    roughnessvalue=0.04
    return sobek.read_profiles(case, roughnesstype, roughnessvalue, pattern=pattern)

def test_sbk_branches():
    sobek = Sobek(r"c:\SK215003\TKI3_NZV.lit")
    case = "Boezemmodel 0D1D"
    pattern = "rStor.*"
    return sobek.read_branches(case, pattern=pattern)



#%%
def merge_dummy_branches(gdf, dummy_gdf, append_dummies=False):
    gdf.reset_index(inplace=True, drop=True)
    nodes = _get_nodes(gdf)
    for index, row in dummy_gdf.iterrows():
        #print(index)
        point = Point(row["geometry"].coords[-1])
        if nodes.distance(point).min() >= 1:
            # cut branch on point
            buffer = point.buffer(1)
            intersects_gdf = gdf[gdf["geometry"].intersects(buffer)]
            if not intersects_gdf.empty:
                branch = gdf[gdf["geometry"].intersects(buffer)].iloc[0]
                gdf_index = branch.name.copy()
                line = branch["geometry"]
                distance = line.project(point)
                line1, line2 = _cut_line(line, distance)
    
                # modify original gdf and append new line
                gdf.at[gdf_index, "geometry"] = line1
                gdf.at[gdf_index, "code"] = f"{branch['code']}_1"
    
                # modify original gdf and append new line
                branch.name =  len(gdf)
                branch["geometry"] = line2
                branch["code"] = f"{branch['code']}_2"
    
                gdf = gdf.append(branch)
                # print(gdf["code"].to_list())
                nodes = _get_nodes(gdf)
            else:
                print(f"warning {row.name} does not intersect a line")
    
    #merge the dummy_gdf and reset index
    drop_cols = [i for i in dummy_gdf if i not in gdf.columns]
    dummy_gdf.drop(columns=drop_cols, axis=1)
    if append_dummies:
        gdf = gdf.append(dummy_gdf)
    gdf.reset_index(inplace=True, drop=True)

    return gdf

def test_merge_dummies(gdf,dummy_gdf):
    dummy_gdf = dummy_gdf.loc[dummy_gdf["code"].isin(["rStorGFE04919", 
                                                      "rStorGFE04927",
                                                      "rStorGFE04929",
                                                      "rStorGFE04931"])]
    gdf = gdf.loc[gdf["code"].isin(["OAF004749",
                                    "OAF004750"])]
    
    gdf_out = merge_dummy_branches(gdf, dummy_gdf)
    
    gdf_out.to_file("branches_and_dummies_test.shp")
    
def yz_fixer(yz):
    y = np.array([i[0] for i in yz])
    removals = np.all([y[1:-1] == y[2:],y[1:-1] == y[0:-2]], axis=0)
    return np.array([list(i) for idx, i in enumerate(yz[1:-1]) if not removals[idx]])


def write_rr_boundaries(rr_writer):
    output_dir = Path(rr_writer.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    filepath = Path(rr_writer.output_dir).joinpath('BoundaryConditions.bc')
    header = {'majorVersion':'1', 'minorVersion':'0', 'fileType':'boundConds'}
    with open(filepath, 'w') as f:                        
        rr_writer._write_dict(f, header, 'General','\n')            
        for _, dct in rr_writer.rrmodel.external_forcings.boundary_nodes.items():                
            temp = {"name":''+dct['id'], 'function':'constant','quantity':'water_level','unit':'m'} 
            rr_writer._write_dict(f,temp,'Boundary','    0\n\n')

def generate_meteo_series(mm_day,
                          start_datetime,
                          end_datetime,
                          timedelta=pd.Timedelta(hours=1)):
    timestamps = int(((end_datetime - start_datetime)/ timedelta) + 1.5)
    value = mm_day/(pd.Timedelta(days=1)/timedelta)
    index = [start_datetime + timedelta * i for i in range(0,int(timestamps))]
    return pd.Series(data=[value]* len(index), index=index)