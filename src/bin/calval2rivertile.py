#!/usr/bin/env python
'''
Copyright (c) 2022-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import os
import ast
import argparse
import logging

import RDF
import SWOTRiver.Estimate
import SWOTRiver.products.calval
from SWOTRiver.errors import RiverObsException

LOGGER = logging.getLogger('calval2river')

def main():
    """Sample script for running calval data through RiverObs"""
    parser = argparse.ArgumentParser()
    parser.add_argument('gps_profile', help='GPS profile file')
    parser.add_argument('out_pixc_file', help='Output pixc file')
    parser.add_argument('out_riverobs_file', help='Output river NETCDF file')
    parser.add_argument('out_pixcvec_file', help='Output PIXCVec file')
    parser.add_argument('rdf_file', help='Static config params')
    parser.add_argument(
        '-l', '--log-level', type=str, default="info",
        help="logging level, one of: debug info warning error")
    args = parser.parse_args()

    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR
            }[args.log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

    config = RDF.RDF()
    config.rdfParse(args.rdf_file)
    config = dict(config)

    pixc_simple = SWOTRiver.products.calval.SimplePixelCloud.from_any(
        args.gps_profile)
    # TODO: no need to write the file if it already is a SimplePixelCloud
    pixc_simple.to_ncfile(args.out_pixc_file)

    # typecast most config values with eval since RDF won't do it for me
    # (excluding strings)
    for key in config.keys():
        if key in ['geolocation_method', 'reach_db_path', 'height_agg_method',
                   'area_agg_method', 'slope_method', 'outlier_method',
                   'pixc_quality_handling']:
            if config[key].lower() != 'none':
                continue
        config[key] = ast.literal_eval(config[key])


    estimator = SWOTRiver.Estimate.CalValToRiverTile(
        args.out_pixc_file, args.out_pixcvec_file)
    estimator.load_config(config)

    # generate empty output file on errors
    try:
        estimator.do_river_processing()
    except RiverObsException as exception:
        LOGGER.error(
            'Unable to continue river processing: {}'.format(exception))

    # Build and write a rivertile-style output file
    estimator.build_products()
    estimator.rivertile_product.to_ncfile(args.out_riverobs_file)

if __name__ == "__main__":
    main()
