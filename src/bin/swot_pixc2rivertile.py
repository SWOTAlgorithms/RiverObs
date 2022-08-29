#!/usr/bin/env python
"""
Stand-in for RiverObs SDS-like processing

Useage:
swot_pixc2rivertile.py l2pixc rivertile.nc pixcvector.nc config.rdf

Optional args:
--shpbas nodebase reachbase -- writes shapefiles in nodebase.shp reachbase.shp
--sensor-file sensor.nc -- gets sensor info from sensor.nc
--gdem-file gdem.nc     -- will make a fake pixel cloud from a gdem and run
                           riverobs on that instead.

template config file:

reach_db_path               (-) = REPLACE_ME
class_list                  (-) = [2, 3, 4, 5, 6, 7]
use_fractional_inundation   (-) = [True, True, False, False, False, False]
use_segmentation            (-) = [False, True, True, True, True, True]
use_heights                 (-) = [False, True, True, False, True, True]
min_points                  (-) = 0
minobs                      (-) = 1
trim_ends                   (-) = False
min_fit_points              (-) = 2
do_improved_geolocation     (-) = True
geolocation_method          (-) = taylor
height_agg_method           (-) = weight
area_agg_method             (-) = composite
slope_method                (-) = bayes
prior_wse_method            (-) = fit
prior_unc_alpha             (-) = 3.0
char_length_tau             (-) = 10000
use_multiple_reaches        (-) = True
use_ext_dist_coef           (-) = True
outlier_method              (-) = piecewise_linear
outlier_abs_thresh          (-) = 1.5
outlier_rel_thresh          (-) = 68
outlier_upr_thresh          (-) = 80
outlier_iter_num            (-) = 30
outlier_breakpoint_min_dist (-) = 0.1
outlier_edge_min_dist       (-) = 0.1
outlier_n_boot              (-) = 10
pixc_quality_handling       (-) = nominal
num_good_sus_pix_thresh_wse  (-) = 1
num_good_sus_pix_thresh_area (-) = 1
use_bright_land              (-) = True
geo_qual_wse_suspect         (-) = 0x0000ffff
geo_qual_wse_degraded        (-) = 0x01ff0000
geo_qual_wse_bad             (-) = 0xfe000000
class_qual_area_suspect      (-) = 0x0000ffff
class_qual_area_degraded     (-) = 0x01ff0000
class_qual_area_bad          (-) = 0xfe000000
sig0_qual_suspect            (-) = 0x01ffffff
sig0_qual_bad                (-) = 0xfe000000

Config file just has processing parameters, no filenames (shape_file_root
will be overwritten in SDS env with "prior_rivers" in current
working directory by SDS pre-processor).

For using with GDEMs change to these key/value pairs:
class_list                (-) = [4, 5, 24]
use_fractional_inundation (-) = [False, False, False]
use_segmentation          (-) = [True, True, True]
use_heights               (-) = [True, True, True]
do_improved_geolocation   (-) = False
slope_method              (-) = first_to_last
use_ext_dist_coef         (-) = False
outlier_method            (-) = None

Author (s): Alex Fore
"""
import sys
import os
import ast
import argparse
import netCDF4
import logging
import subprocess

import RDF
import SWOTRiver.Estimate
from SWOTRiver.products.pixcvec import L2PIXCVector
from SWOTRiver.errors import RiverObsException

LOGGER = logging.getLogger('swot_pixc2rivertile')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_file', help='pixel cloud file')
    parser.add_argument('out_riverobs_file', help='Output NetCDF file')
    parser.add_argument('out_pixc_vector_file', help='Output PIXC vector file')
    parser.add_argument('rdf_file', help='Static config params')
    parser.add_argument('--shpbase', type=str, nargs=2, default=[None, None])
    parser.add_argument(
        '-l', '--log-level', type=str, default="info",
        help="logging level, one of: debug info warning error")
    parser.add_argument(
        '--gdem-file', '-g', type=str, default=None,
        help="GDEM file; if commanded makes a fake pixc from GDEM and runs"+
             "RiverObs on that instead of on pixc_file")
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
        if key in ['geolocation_method', 'reach_db_path', 'height_agg_method',
                   'area_agg_method', 'slope_method', 'outlier_method',
                   'prior_wse_method', 'pixc_quality_handling']:
            if config[key].lower() != 'none':
                continue
        config[key] = ast.literal_eval(config[key])

    pixc_file = args.pixc_file
    if args.gdem_file is not None:
        import fake_pixc_from_gdem
        import tempfile
        pixc_file = tempfile.mktemp()
        fake_pixc_from_gdem.fake_pixc_from_gdem(
            args.gdem_file, args.pixc_file, pixc_file)

    l2pixc_to_rivertile = SWOTRiver.Estimate.L2PixcToRiverTile(
            pixc_file, args.out_pixc_vector_file)

    l2pixc_to_rivertile.load_config(config)

    # generate empty output file on errors
    try:
        l2pixc_to_rivertile.validate_inputs()
        l2pixc_to_rivertile.do_river_processing()
        l2pixc_to_rivertile.match_pixc_idx()
        l2pixc_to_rivertile.do_improved_geolocation()

    except RiverObsException as exception:
        LOGGER.error(
            'Unable to continue river processing: {}'.format(exception))

    l2pixc_to_rivertile.build_products()

    l2pixc_to_rivertile.rivertile_product.to_ncfile(args.out_riverobs_file)
    if args.shpbase[0] is not None:
        l2pixc_to_rivertile.rivertile_product.nodes.write_shapes(
            args.shpbase[0]+'.shp')
        l2pixc_to_rivertile.rivertile_product.reaches.write_shapes(
            args.shpbase[1]+'.shp')

    if args.gdem_file is not None:
        os.remove(pixc_file)

if __name__ == "__main__":
    main()
