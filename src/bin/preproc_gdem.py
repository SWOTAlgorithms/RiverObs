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

def select_river_labels(type, type_label, gdem_x, gdem_y, reaches):
    """
    Picks the labels that are nearest to the reach lons/lats.  For each reach
    finds the closest label.
    """
    LOGGER.info('select_river_label')

    min_compare = 3200 # 10m * 3.125m * 3200 == 2*200m*250m (2 nodes)
    max_compare = 100000

    uniq_labels = np.unique(type_label[type == 1])
    cnts, _ = np.histogram(type_label[type == 1], uniq_labels)

    idxsort = np.argsort(cnts)[::-1]
    uniq_labels = uniq_labels[idxsort]
    cnts = cnts[idxsort]

    labels = []
    for ireach, reach in enumerate(reaches):
        type_dist = 9999999*np.ones(uniq_labels.shape)
        for ilabel, uniq_label in enumerate(uniq_labels):

            # probably not it if only 1000 pixels in feature
            if cnts[ilabel] < min_compare:
                continue

            these_x = gdem_x[type_label == uniq_label]
            these_y = gdem_y[type_label == uniq_label]

            if cnts[ilabel] > max_compare:
                these_x, these_y = np.random.permutation(
                    np.array([these_x, these_y]))[:, :max_compare]

            delta2 = (
                (these_x[:, np.newaxis] - reach.x)**2 +
                (these_y[:, np.newaxis] - reach.y)**2)
            min_d2 = delta2.min(axis=0)
            type_dist[ilabel] = np.mean(np.sqrt(min_d2))

            LOGGER.debug(
                "reach, label, dist, num: {} {} {} {}".format(
                    ireach, uniq_label, type_dist[ilabel], len(these_x)))

        this_label = uniq_labels[type_dist.argmin()]
        if this_label not in labels:
            labels.append(this_label)

    return labels

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_gdem_file', help='Input GDEM file')
    parser.add_argument('out_gdem_file', help='Output GDEM file')
    parser.add_argument('reachdb_path', help='reach DB path/file')
    parser.add_argument(
        '-l', '--log-level', type=str, default="warning",
        help="logging level, one of: debug info warning error")
    parser.add_argument('--plot', default=False, action='store_true')
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
    land_label = 0
    river_labels = select_river_labels(
        type, type_label, gdem_x, gdem_y, reaches)

    # Water that is not in the largest water features
    water_not_main_label = np.logical_and(
        np.isin(type_label, river_labels, invert=True),
        type_label != land_label)

    out_type = type.copy()
    out_type[water_not_main_label] = 0

    if args.plot:
        import matplotlib.pyplot as plt

        reach_lon, reach_lat = np.array([]), np.array([])
        for item in reaches:
            reach_lon = np.append(reach_lon, item.lon)
            reach_lat = np.append(reach_lat, item.lat)

        figure, axis = plt.subplots()
        axis.plot(reach_lon, reach_lat, 'rs')
        axis.plot(lon[type==1], lat[type==1], 'g.')
        axis.plot(lon[out_type==1], lat[out_type==1], 'k.')
        axis.set_xlabel('longitude')
        axis.set_ylabel('latitude')

        figure, axis = plt.subplots()
        axis.imshow(type)

        plt.show()

    shutil.copy(args.in_gdem_file, args.out_gdem_file)
    with netCDF4.Dataset(args.out_gdem_file, 'a') as ofp:
        ifp.variables['landtype'][:] = out_type

if __name__ == "__main__":
    main()
