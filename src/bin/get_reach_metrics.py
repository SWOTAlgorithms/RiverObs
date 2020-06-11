#!/usr/bin/env python
'''
Copyright (c) 2020-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

# A collection of functions for the RiverObs metanalysis
Author(s): Cassie Stuurman
'''
import numpy as np
import pdb
from scipy.io import loadmat

def get_quadratic_fit(nodes, wse):
    fit = np.polyfit(nodes, wse, 2)
    quad_coefficient = round(fit[0], 3)
    lin_coefficient = round(fit[1], 3)
    return quad_coefficient, lin_coefficient

def get_lake_proximity(reach, reach_lat, reach_lon):
    lake_data = loadmat('/u/swot-fn-r0/fore/rivers-lakes-closeness/lakes_near_rivers/matches.mat')
    reach_ids = lake_data['reach_id']
    reach_lats = lake_data['reach_lat']
    reach_lons = lake_data['reach_lon']
    lake_dist = lake_data['min_dist']
    proximity = lake_dist[np.where(reach_ids==reach)]
    if not proximity: # theres no matching reach id
        proximity = lake_dist[np.where(reach_lats==reach_lat)]
    if not proximity: # theres no matching latitude
        # find minimum distance reach
        try:
            lat_diff = abs(reach_lats - reach_lat[0])
            lon_diff = abs(reach_lons - reach_lon[0])
            diff = lat_diff + lon_diff
            min = lake_dist[np.where(diff == np.nanmin(diff))]
            proximity = round(np.min(min), 3)
        except IndexError: # there's something wrong with the input lat/lon. Temp fix since the new database will solve this
            proximity=np.nan

    return proximity

