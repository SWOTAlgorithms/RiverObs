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
        type = ifp.variables['landtype'][:]
        lat = ifp.variables['latitude'][:]
        lon = ifp.variables['longitude'][:]

    # labeled in descending order of counts
    type_label, num_labels = scipy.ndimage.label(type)

    # Get land and river labels (assume largest connected water feature is the
    # river of interest).
    land_label = scipy.stats.mode(type_label[type == 0]).mode[0]
    river_label = scipy.stats.mode(type_label[type == 1]).mode[0]

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
