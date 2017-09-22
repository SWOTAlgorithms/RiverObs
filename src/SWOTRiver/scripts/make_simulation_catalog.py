#! /usr/bin/env python

"""
Make a GIS file containing the lat/lon polygons covered by
a SWOT simulation run.
"""

from __future__ import absolute_import, division, print_function

# make sure the libraries are importable

def search_for_libraries():
    """Search for the libraries."""

    import os, os.path
    import sys

    # Try importing the root library
    try:
        from GDALOGRUtilities import OGRWriter
        from shapely.geometry import box
    except:
        sys.stderr.write("Libraries not found. Make sure you are running in the SWOTRiver environment.\n")
        sys.exit(1)

search_for_libraries()

# Imports

from os.path import join, split
import argparse
from glob import glob
from GDALOGRUtilities import OGRWriter
from shapely.geometry import box

def parse_inputs():
    """Create the argument parser."""

    parser = argparse.ArgumentParser(description='Generate GIS file with coverage polygons for a SWOT simulation.')

    parser.add_argument('output_file',help='Output file name')
    parser.add_argument('data_directories',metavar='data_dir',
                        nargs='+',help='Director containing simulation data')

    parser.add_argument('-f','--format',help="OGR file format  (default 'ESRI Shapefile')",
                        default='ESRI Shapefile')
    parser.add_argument('--dlat',type=float,help='latitude tile size (default +1)',
                        default=1.)
    parser.add_argument('--dlon',type=float,help='longitude tile size (default +1)',
                        default=1.)

    args = parser.parse_args()

    return args

def get_bbox_from_files(sim_files,dlat,dlon):
    """Given a list of files, return a list of bounding box
    coordinates suitable for passing to OGRWriter."""

    bounding_boxes = []
    for name in sim_files:
        location = split(name)[-1].split('_')[0]
        EW = location[0]
        NS = location[4]

        if EW.lower() == 'w':
            lonmin = -float(location[1:4])
        else:
            lonmin = float(location[1:4])

        if NS.lower() == 's':
            latmin = -float(location[5:7])
        else:
            latmin = float(location[5:7])

        bounding_boxes.append(box(lonmin,latmin,lonmin+dlon,latmin+dlat))

    return bounding_boxes

def write_catalog(output_file,format,bounding_boxes,sim_files):
    """Write coverage polygons with metadata."""

    fields = {'file':('str',42)}

    writer = OGRWriter(output_file,fields=fields,driver=format,
                       geometry='Polygon')

    for i,bbox in enumerate(bounding_boxes):
        field_record = {'file':split(sim_files[i])[-1]}
        print(split(sim_files[i])[-1][0:7],bbox.wkt)
        writer.add_wkt_feature(bbox.wkt,field_record)
    writer.close()


def main():

    args = parse_inputs()

    # Get a list of the files

    sim_files = []
    for directory in args.data_directories:
        files = glob(join(directory,'*_cycle_*_pass_*.*.nc'))
        if len(files) > 0:
            sim_files += files

    bounding_boxes = get_bbox_from_files(sim_files,args.dlat,args.dlon)
    print(bounding_boxes)

    write_catalog(args.output_file,args.format,bounding_boxes,sim_files)

    print('Successfuly wrote catalog file')


if __name__ == '__main__': main()
