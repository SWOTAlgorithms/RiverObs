#!/usr/bin/env python
"""
Stand-in for RiverObs SDS-like processing

Author (s): Alex Fore
"""
import sys
import os
import ast
import argparse
import netCDF4

import RDF
import SWOTRiver.Estimate
import RiverObs.ShapeWriter
import RiverObs.NetCDFReachWriter

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_file', help='pixel cloud file')
    parser.add_argument('out_riverobs_file', help='Output NetCDF file')
    parser.add_argument('out_pixc_vector_file', help='Output PIXC vector file')
    parser.add_argument('rdf_file', help='Static config params')
    parser.add_argument('--shpbasedir', type=str, default=None)
    parser.add_argument('--sensor-file', type=str, default=None)
    args = parser.parse_args()

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

    is_new_pixc = True if args.sensor_file is None else False

    l2pixc_to_rivertile = SWOTRiver.Estimate.L2PixcToRiverTile(
            args.pixc_file, args.out_pixc_vector_file, args.sensor_file)

    l2pixc_to_rivertile.load_config(config)
    l2pixc_to_rivertile.do_river_processing()
    l2pixc_to_rivertile.match_pixc_idx()
    l2pixc_to_rivertile.do_improved_geolocation()

    RiverObs.NetCDFReachWriter.write(
        args.out_riverobs_file, l2pixc_to_rivertile.node_outputs,
        l2pixc_to_rivertile.reach_outputs)

    RiverObs.NetCDFReachWriter.fixup_metadata(args.out_riverobs_file)

    # optional shapefile outputs
    if args.shpbasedir is not None:
        try:
            RiverObs.ShapeWriter.write(
                l2pixc_to_rivertile.reach_collection,
                os.path.join(args.shpbasedir, 'nodes'),
                os.path.join(args.shpbasedir, 'reaches'))

        # No reaches found, skip writing of shapefiles
        except IndexError:
            pass


if __name__ == "__main__":
    main()
