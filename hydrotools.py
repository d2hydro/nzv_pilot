# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 23:18:13 2020

@author: danie
"""

import geopandas as gpd
from delft3dfmpy import HyDAMO
from pathlib import Path

from shapely.geometry import Polygon, LineString, Point
from shapely.ops import snap


hydamo = HyDAMO()

def _filter(gdf, attribute_filter):
    if isinstance(attribute_filter, dict):
        for key, value in attribute_filter.items():
            if not isinstance(value, list):
                value = [value]
            gdf = gdf[gdf[key].isin(value)]
            
        return gdf
    else:
       raise IOError('attribute_filter should be dictionary') 


def read_file(path,
              hydamo_attribute,
              index_col=None,
              attribute_filter=None,
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
            column to be used as index (after optional mapping)
        attribute_filter: dict
            dict with lists or strings of the format {'column_name': [values to keep]}
        column_mapping: dict
            dict for renaming input colunns to required columns
            
        Result: GeoDataFrame matching the HyDAMO property
        """
        gdf = gpd.read_file(path)
        gdf.columns = gdf.columns.str.lower()
        
        #filter by attribute
        if attribute_filter:
            attribute_filter = {key.lower():value for key,value in attribute_filter.items()}
            gdf = _filter(gdf, attribute_filter)
            
        #map to hydamo columns
        if column_mapping:
            column_mapping = {key.lower():value.lower() for key,value in column_mapping.items()}
            gdf.rename(columns=column_mapping, inplace=True)
            
        # drop all columns not needed 
        required_columns = getattr(hydamo,hydamo_attribute).required_columns
        
        if hydamo_attribute =='crosssections':
            required_columns += ['order', 'category']
            if z_coord:
                required_columns += ['z']
        
        drop_cols = [col for col in gdf.columns if not col in required_columns + ['geometry']]
        
        if len(drop_cols) > 0:
            gdf = gdf.drop(drop_cols, axis=1)
        
        return gdf
    
def to_file(model, hydamo_attribute, length=False, path=Path('.')):
    '''converts hydamo class to shape-file'''
    
    path = Path(path)
    hydamo_class = getattr(model,hydamo_attribute)
    data = {col:hydamo_class[col].values for col in hydamo_class.columns}
    
    if length:
        data = data = {**data, 'length': hydamo_class['geometry'].length.values}
        
    gpd.GeoDataFrame(data=data).to_file(path.joinpath(f'{hydamo_attribute}.shp'))
    
def snap_ends(gdf, tolerance, digits=None):
    ''' function to snap all end-vertices of a GeoDataFrame with LineStrings within a specified tolerance'''
    
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
                #unfortunately this line doesn't work (always):
                #geom = snap(geom, row_selec['geometry'], snap_tolerance)
                
                #but this does
                for dst_vert in [0,-1]:
                    for src_vert in [0,-1]:
                        geom_coords[dst_vert] = snap(Point(geom_coords[dst_vert]),
                                                     Point(row_selec['geometry'].coords[src_vert]), 
                                                     tolerance=1).coords[0]
#                        print(f'{_}: {geom_coords[-1]} {row_selec["geometry"].coords[src_vert]}')
            geom = LineString(geom_coords)
        # write feature in original GeoDataFrame
        gdf.loc[index, 'geometry'] = geom
        # mark index as snapped
        snapped += [index]
        
    return gdf