#!/usr/bin/env python
"""
Stand-in for RiverObs SDS-like processing

Useage:
    swot_rivertiles2riversp.py ...
"""
import os
import logging
import argparse
import numpy as np

from SWOTRiver.products.riversp import L2HRRiverSP
from SWOTRiver.products.rivertile import L2HRRiverTile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'rivertile_files', type=str, nargs='*',
        help='space serperated list of rivertiles')
    parser.add_argument('river_sp_ncfile', type=str)
    parser.add_argument('--from-shapes', help='From shapefiles',
        action='store_true', default=False)
    parser.add_argument('--shpbasedir', type=str, default=None)
    parser.add_argument(
        '-l', '--log-level', type=str, default="info",
        help="logging level, one of: debug info warning error")
    args = parser.parse_args()

    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

    if args.from_shapes:
        rivertile_files = list(zip(
            args.rivertile_files[0::2], args.rivertile_files[1::2]))
    else:
        rivertile_files = args.rivertile_files

    river_tiles = []
    for ii, rivertile_file in enumerate(rivertile_files):
        if args.from_shapes:
            river_tile = L2HRRiverTile.from_shapes(
                rivertile_file[0], rivertile_file[1])
        else:
            river_tile = L2HRRiverTile.from_ncfile(rivertile_file)

        river_tiles.append(river_tile)

    # write river sp
    river_sp = L2HRRiverSP.from_rivertiles(river_tiles)

    # write shapefile version of river sp
    if args.shpbasedir is not None:
        if not os.path.isdir(args.shpbasedir):
            os.mkdir(args.shpbasedir)
        river_sp.nodes.write_shapes(
            os.path.join(args.shpbasedir, 'nodes.shp'))
        river_sp.reaches.write_shapes(
            os.path.join(args.shpbasedir, 'reaches.shp'))

    if args.from_shapes:
        fill_value = river_sp.nodes.VARIABLES['time_str']['fill_value']
        river_sp.nodes.time_str = (
            np.ones(river_sp.nodes.time.shape)*fill_value)
        river_sp.reaches.time_str = (
            np.ones(river_sp.reaches.time.shape)*fill_value)
    river_sp.to_ncfile(args.river_sp_ncfile)

if __name__ == "__main__":
    main()
