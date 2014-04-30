#! /usr/bin/env python

"""
Given a SWOTL2 file, estimate all of the reaches within it and
output the results to a file
"""

# make sure the libraries are importable

def search_for_libraries():
    """Search for the libraries."""

    import os, os.path
    import sys
        
    # Try importing the root library
    try:
        from SWOTRiver import SWOTRiverEstimator
    except:
        sys.stderr.write("Libraries not found. Make sure you are running in the SWOTRiver environment.\n")
        sys.exit(1)
        
search_for_libraries()

# Imports

from os.path import join, split
import argparse
from glob import glob
from SWOTRiver import SWOTRiverEstimator
from GDALOGRUtilities import OGRWriter
from shapely.geometry import box
from RDF import RDF

RDF_template = """
SWOT L2 DATA INPUTS

good data labels                       = 1 ! List of numbers separated by space
latitude keyword                       = no_layover_latitude
longitude keyword                      = no_layover_longitude
classification keyword                 = no_layover_classification
height measurement keyword             = height
true height keyword                    = water_height
projection                             = laea ! Proj4 projection name
false easting                      (m) = 0
false northing                     (m) = 0
lat_0                            (deg) = None !if no number, then use mean data lat
lon_0                            (deg) = None
ellipsoid                              = WGS84

CENTERLINE DATA BASE INPUTS

data base directory                    = ../../data/nAmerica_GRWDL_river_topo/
shape file root name                   = nAmerica_GRWDL_river_topo
clip to data bounding box              = 1 ! 0/1
clip buffer                      (deg) = 0.1

"""

def parse_inputs():
    """Create the argument parser."""

    parser = argparse.ArgumentParser(description='Generate GIS file with coverage polygons for a SWOT simulation.')

    parser.add_argument('rdf_file',help='Input RDF file')
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
        print split(sim_files[i])[-1][0:7],bbox.wkt
        writer.add_wkt_feature(bbox.wkt,field_record)
    writer.close()
                        

def main():

    args = parse_inputs()

    # Parse the input RDF

    rdf = RDF().rdfParse(args.rdf_file)

    class_list = rdf.int('good data labels')
    if type(class_list) == int:
        class_list = [class_list]
    
    lat_kwd = rdf['latitude keyword']
    lon_kwd = rdf['longitude keyword']
    class_kwd = rdf['classification keyword']
    height_kwd = rdf['height measurement keyword']
    true_height_kwd = rdf['true height keyword']
    proj = rdf['projection']
    x_0 = rdf.float('false easting')
    y_0 = rdf.float('false northing')
    try:
        lat_0 = rdf.float('lat_0')
        lon_0 = rdf.float('lon_0')
    else:
        lat_0 = None
        lon_0 = None
    ellps = rdf['ellipsoid']

    db_dir = rdf['data base directory']
    shape_file_root = join(db_dir,rdf['shape file root name'])
    clip = bool(rdf.int('clip to data bounding box'))
    clip_buffer = rdf.float('clip buffer')
    
    # Get a list of the files

    sim_files = []
    for directory in args.data_directories:
        files = glob(join(directory,'*_cycle_*_pass_*.*.nc'))
        if len(files) > 0:
            sim_files += files

    # Write a catalog of the tiles processed as a GIS file
    
    bounding_boxes = get_bbox_from_files(sim_files,args.dlat,args.dlon)
    write_catalog(args.output_file,args.format,bounding_boxes,sim_files)

    # Process each input file

    for sim_file in sim_files:

        print('Reading file: %s'%sim_file)

        # Initialize the data
        
        estimator = SWOTRiverEstimator(sim_file,class_list=class_list,
                    lat_kwd=lat_kwd, lon_kwd=lon_kwd,
                    class_kwd=class_kwd,
                    height_kwd=height_kwd,true_height_kwd=true_height_kwd,
                    proj=proj,x_0=x_0,y_0=y_0,lat_0=lat_0,lon_0=lon_0,
                    ellps=ellps)
        print('data read')

        # Get the reaches corresponding to these data

        reaches = estimator.get_reaches(shape_file_root, clip=clip,clip_buffer=clip_buffer)
        print('reaches read')

    print('Successfuly estimated river heights and slopes')


if __name__ == '__main__': main()
