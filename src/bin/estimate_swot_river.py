#!/usr/bin/env python
"""
Runs RiverObs on a swot river.

Most of code lives in SWOTRiver.EstimateSWOTRiver.estimate
"""
import os
import ast
import argparse
import SWOTRiver.EstimateSWOTRiver
import RiverObs.RiverReachWriter
import RDF

def main():
    """
    Does the running of RiverObs on a swot river.
    """
    parser = argparse.ArgumentParser(
        description='Estimate heigt and slope for a SWOT simulation.',
        usage=SWOTRiver.EstimateSWOTRiver.rdf_file_template)

    parser.add_argument('rdf_file', help='Input RDF file')
    parser.add_argument(
        '-f','--format', help="OGR file format  (default 'ESRI Shapefile')",
        default='ESRI Shapefile')
    args = parser.parse_args()

    params = RDF.RDF_to_class(
        SWOTRiver.EstimateSWOTRiver.input_vars, file=args.rdf_file)
    print((params.l2_file))
    print((params.fout_reach))
    print((params.fout_node))
    print((params.fout_index))
    # Reformat some inputs
    print(params.lonmin)

    # Heavy lifting done in SWOTRiver.EstimateSWOTRiver.estimate
    node_variables, reach_variables, reach_collection = \
        SWOTRiver.EstimateSWOTRiver.estimate(params)

    # Write shapefiles
    writer = RiverObs.RiverReachWriter(
        reach_collection, node_variables, reach_variables)
    driver = args.format
    writer.write_nodes_ogr(params.fout_node, driver=driver)
    writer.write_reaches_ogr(params.fout_reach, driver=driver)
    print('Successfuly estimated river heights and slopes')

if __name__ == "__main__":
    main()
