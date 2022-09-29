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
import warnings
import tempfile

import RDF
import SWOTRiver.Estimate
from SWOTRiver.products.calval import \
    SimplePixelCloud, Drifter, PressureTransducers
from SWOTRiver.errors import RiverObsException

LOGGER = logging.getLogger('calval2rivertile')

FORMATS = ['simple_pixc', 'drifter', 'geotiff', 'pt']

def main():
    """Sample script for running calval data through RiverObs"""
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Input calval file to process')
    parser.add_argument(
        'format', help='Input calval file format', choices=FORMATS)
    parser.add_argument(
        'out_riverobs_file', help='Output RiverTile NETCDF file')
    parser.add_argument('out_pixcvec_file', help='Output PIXCVecRiver file')
    parser.add_argument('rdf_file', help='Static config params')
    parser.add_argument(
        'out_pixc_file',
        help='If specified, write reformatted SimplePixelCloud to this file',
        default=None, nargs='?')
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

    if args.format == 'drifter':
        # Use file extension to determine Drifter classmethod to use
        if args.input_file.lower().endswith('.nc'):
            drifter = Drifter.from_ncfile(args.input_file)

        elif args.input_file.lower().endswith('.shp'):
            drifter = Drifter.from_shp(args.input_file)

        elif args.input_file.lower().endswith('.txt'):
            drifter = Drifter.from_native(args.input_file)

        else:
            raise Exception("Unknown drifter file format!")

        pixc_simple = SimplePixelCloud.from_drifter(drifter)

    elif args.format == 'geotiff':
        pixc_simple = SimplePixelCloud.from_geotif(args.input_file)

    elif args.format == 'pt':
        pixc_simple = SimplePixelCloud.from_pressure_transducer(
            PressureTransducers.from_native(args.input_file))

    elif args.format == 'simple_pixc':
        pixc_simple = SimplePixelCloud.from_ncfile(args.input_file)

    output_pixc_file = args.out_pixc_file

    # If args.out_pixc_file not specified use a tempfile and clean it later
    if args.out_pixc_file is None:
        fid, output_pixc_file = tempfile.mkstemp()
    pixc_simple.to_ncfile(output_pixc_file)

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
        output_pixc_file, args.out_pixcvec_file)
    estimator.load_config(config)

    # generate empty output file on errors
    try:
        estimator.do_river_processing()
    except RiverObsException as exception:
        LOGGER.error(
            'Unable to continue river processing: {}'.format(exception))

    # Build and write a rivertile-style output file
    estimator.build_products()

    # Suppress litany of warnings from Product class writer
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        estimator.rivertile_product.to_ncfile(args.out_riverobs_file)

    # Delete tempfile if args.out_pixc_file is not specified
    if args.out_pixc_file is None:
        os.remove(output_pixc_file)

if __name__ == "__main__":
    main()
