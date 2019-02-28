#!/usr/bin/env python
"""
Pre-processes a GDEM for use with RiverObs validation studies

Author (s): Alex Fore
"""
import shutil
import netCDF4
import argparse
import logging
import scipy.stats
import scipy.ndimage
import numpy as np

import RiverObs.ReachDatabase

LOGGER = logging.getLogger(__name__)

def select_river_label(type, type_label, gdem_x, gdem_y, reaches):
    """
    Picks the label that is nearest to the reach lons/lats.
    """
    LOGGER.info('select_river_label')
    reach_x, reach_y = np.array([]), np.array([])
    for item in reaches:
        reach_x = np.append(reach_x, item.x)
        reach_y = np.append(reach_y, item.y)

    uniq_labels = np.unique(type_label[type == 1])
    type_dist = 9999999*np.ones(uniq_labels.shape)

    for ilabel, uniq_label in enumerate(uniq_labels[:20]):
        these_x = gdem_x[type_label == uniq_label]
        these_y = gdem_y[type_label == uniq_label]
        delta2 = (
            (these_x[:, np.newaxis] - reach_x)**2 +
            (these_y[:, np.newaxis] - reach_y)**2)
        min_d2 = delta2.min(axis=0)
        type_dist[ilabel] = np.mean(np.sqrt(min_d2))

        LOGGER.debug(
            "Feature label, dist, num: {} {} {}".format(
                uniq_label, type_dist[ilabel], len(these_x)))

    return uniq_labels[type_dist.argmin()]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_gdem_file', help='Input GDEM file')
    parser.add_argument('out_gdem_file', help='Output GDEM file')
    parser.add_argument('reachdb_path', help='reach DB path/file')
    parser.add_argument(
        '-l', '--log-level', type=str, default="warning",
        help="logging level, one of: debug info warning error")
    args = parser.parse_args()

    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

    with netCDF4.Dataset(args.in_gdem_file, 'r') as ifp:
        in_type = ifp.variables['landtype'][:]
        lat = ifp.variables['latitude'][:]
        lon = ifp.variables['longitude'][:]

    type = np.zeros(in_type.shape)
    type[in_type == 1] = 1

    # Extract prior node locations, first make a LatLonRegion
    llbox = RiverObs.ReachDatabase.LatLonRegion(
        [lon.min(), lat.min(), lon.max(), lat.max()])

    # project GDEM coordinates
    gdem_x, gdem_y = llbox.proj(lon, lat)

    # Extract Reaches
    reaches = RiverObs.ReachDatabase.ReachExtractor(args.reachdb_path, llbox)

    # Labeled in descending order of counts
    type_label, num_labels = scipy.ndimage.label(type)

    # Get land and river labels.
    land_label = scipy.stats.mode(type_label[type == 0]).mode[0]
    river_label = select_river_label(type, type_label, gdem_x, gdem_y, reaches)

    # Water that is not in the largest water features
    water_not_main_label = np.logical_and(
        type_label != river_label, type_label != land_label)

    out_type = type.copy()
    out_type[water_not_main_label] = 0

    shutil.copy(args.in_gdem_file, args.out_gdem_file)
    with netCDF4.Dataset(args.out_gdem_file, 'a') as ofp:
        ifp.variables['landtype'][:] = out_type

if __name__ == "__main__":
    main()
