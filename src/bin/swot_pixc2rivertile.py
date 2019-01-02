#!/usr/bin/env python
"""
Stand-in for RiverObs SDS-like processing

Useage:
swot_pixc2rivertile.py l2pixc rivertile.nc pixcvector.nc config.rdf

Optional args:
--shpbasedir dirname    -- writes shapefiles in dirname/reaches dirname/nodes
--sensor-file sensor.nc -- gets sensor info from sensor.nc


template config file:

width_db_file             (-) = None
use_width_db              (-) = False
shape_file_root           (-) = /u/onde-r0/fore/data/River_Prior_Database/ADT_priordatabase_reaches_polylines/NA_reaches_data_discharge_depth_chn_grdc_revised_GCS
class_list                (-) = [2, 3, 4, 5, 6]
use_fractional_inundation (-) = [True, True, False, False, False]
use_segmentation          (-) = [False, True, True, True, True]
use_heights               (-) = [False, False, True, False, False]
min_points                (-) = 100
clip_buffer               (-) = 20.0
ds                        (-) = 300.0
refine_centerline         (-) = False
smooth                    (-) = 0.01
alpha                     (-) = 1
max_iter                  (-) = 1
scalar_max_width          (-) = 600.0
minobs                    (-) = 10
trim_ends                 (-) = False
fit_types                 (-) = ['OLS', 'WLS', 'RLM']
min_fit_points            (-) = 3
do_improved_geolocation   (-) = True
geolocation_method        (-) = taylor
height_agg_method         (-) = weight
area_agg_method           (-) = composite

Config file just has processing parameters, no filenames (shape_file_root
will be overwritten in SDS env with "prior_rivers" in current
working directory by SDS pre-processor).

For using with GDEMs change to these key/value pairs:
class_list                (-) = [1,]
use_fractional_inundation (-) = [False,]
use_segmentation          (-) = [True,]
use_heights               (-) = [True,]

Author (s): Alex Fore
"""
import sys
import os
import ast
import argparse
import netCDF4
import logging

import RDF
import SWOTRiver.Estimate

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_file', help='pixel cloud file')
    parser.add_argument('out_riverobs_file', help='Output NetCDF file')
    parser.add_argument('out_pixc_vector_file', help='Output PIXC vector file')
    parser.add_argument('rdf_file', help='Static config params')
    parser.add_argument('--shpbasedir', type=str, default=None)
    parser.add_argument('--sensor-file', type=str, default=None)
    parser.add_argument(
        '-l', '--log-level', type=str, default="warning",
        help="logging level, one of: debug info warning error")
    args = parser.parse_args()

    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

    config = RDF.RDF()
    config.rdfParse(args.rdf_file)
    config = dict(config)

    # typecast most config values with eval since RDF won't do it for me
    # (excluding strings)
    for key in config.keys():
        if any([key == item for item in [
            'geolocation_method', 'shape_file_root']]):
            continue
        config[key] = ast.literal_eval(config[key])

    l2pixc_to_rivertile = SWOTRiver.Estimate.L2PixcToRiverTile(
            args.pixc_file, args.out_pixc_vector_file, args.sensor_file)

    l2pixc_to_rivertile.load_config(config)
    l2pixc_to_rivertile.do_river_processing()
    l2pixc_to_rivertile.match_pixc_idx()
    l2pixc_to_rivertile.do_improved_geolocation()
    l2pixc_to_rivertile.flag_lakes_pixc()
    l2pixc_to_rivertile.build_products()

    l2pixc_to_rivertile.rivertile_product.to_ncfile(args.out_riverobs_file)
    if args.shpbasedir is not None:
        if not os.path.isdir(args.shpbasedir):
            os.path.mkdir(args.shpbasedir)
        l2pixc_to_rivertile.rivertile_product.nodes.write_shapes(
            os.path.join(args.shpbasedir, 'nodes.shp'))
        l2pixc_to_rivertile.rivertile_product.reaches.write_shapes(
            os.path.join(args.shpbasedir, 'reaches.shp'))

if __name__ == "__main__":
    main()
