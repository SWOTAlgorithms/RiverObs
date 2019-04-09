#! /usr/bin/env python
"""
Given a SWOTL2 file, estimate all of the reaches within it and
output the results to a file
"""

from __future__ import absolute_import, division, print_function

# make sure the libraries are importable


def search_for_libraries():
    """Search for the libraries."""

    import os, os.path
    import sys

    # Try importing the root library
    try:
        from SWOTRiver import SWOTRiverEstimator
    except:
        sys.stderr.write(
            "Libraries not found. Make sure you are running in the SWOTRiver environment.\n"
        )
        sys.exit(1)


#search_for_libraries()

# Imports

from os.path import join, split
import argparse
from glob import glob
from SWOTRiver import SWOTRiverEstimator
from RiverObs import WidthDataBase
#from SWOTRiver import WidthDataBase
from GDALOGRUtilities import OGRWriter
from shapely.geometry import box
from RDF import RDF_to_class, RDF

rdf_file_template = """%(prog)s [-h] rdf_file data_dir [data_dir ...]

Estimate river slopes and height measurements from simulated data.

positional arguments:
  rdf_file    rdf file containing the inputs.
  data_dir    directories containing simulated data.

optional arguments:
  -h, --help   show this help message and exit
  -f, --format OGR file format  (default 'ESRI Shapefile')
  --dlat       latitude tile size (default +1)
  --dlon       longitude tile size (default +1)

A template RDF file follows:

SWOT L2 DATA INPUTS

good data labels                       = 1 ! List of numbers separated by space
latitude keyword                       = no_layover_latitude
longitude keyword                      = no_layover_longitude
classification keyword                 = no_layover_classification
height measurement keyword             = height
true height keyword                    = water_height
no noise height keyword                = no_noise_height
xtrack keyword                         = no_layover_cross_track
projection                             = laea ! Proj4 projection name
false easting                      (m) = 0
false northing                     (m) = 0
lat_0                            (deg) = None !if no number, then use mean data lat
lon_0                            (deg) = None
ellipsoid                              = WGS84

CENTERLINE DATA BASE INPUTS

reach data base directory              = ../../data/nAmerica_GRWDL_river_topo/
shape file root name                   = nAmerica_GRWDL_river_topo
width data base directory              = ../../data/
width db file name                     = nAmerica_GRWDL.h5
clip to data bounding box              = 1 ! 0/1
clip buffer                      (deg) = 0.1
use width db                           = 1 ! 0/1

CENTERLINE REFINEMENT PARAMETERS

centerline spacing                 (m) = 50.
refine centerline                      = 1 ! 0/1
smooth parameter                       = 1.e-2 ! between 0 and 1
alpha                                  = 1 ! <= 1
maximum iterations                     = 1 ! 1 or 2 recommended
scalar maximum width               (m) = 600. ! look for points in this width

ESTIMATION PARAMETERS

subreach size                      (m) = 10.e3
minimum number of points to fit        = 150
step                               (m) = 2.5e3
fit types                              = OLS WLS RLM ! space separate list of OLS WLS RLM
starting along track distance      (m) = 0.
minimum observations per node          = 10

OUTPUT OPTIONS

output file                            = sacramento_reach_estimates.h5
"""

# Input directory for RDF_to_class

input_vars = {
    # SWOT L2 DATA INPUTS
    'good data labels': ('class_list', 'd'),
    'latitude keyword': ('lat_kwd', 's'),
    'longitude keyword': ('lon_kwd', 's'),
    'classification keyword': ('class_kwd', 's'),
    'height measurement keyword': ('height_kwd', 's'),
    'true height keyword': ('true_height_kwd', 's'),
    'no noise height keyword': ('no_noise_height_kwd', 's'),
    'xtrack keyword': ('xtrack_kwd', 's'),
    'projection': ('proj', 's'),
    'false easting': ('x_0', 'f'),
    'false northing': ('y_0', 'f'),
    'lat_0': ('lat_0', 's'),
    'lon_0': ('lon_0', 's'),
    'ellipsoid': ('ellps', 's'),

    # CENTERLINE DATA BASE INPUTS
    'reach data base directory': ('reach_db_dir', 's'),
    'shape file root name': ('shape_file_root', 's'),
    'width data base directory': ('width_db_dir', 's'),
    'width db file name': ('width_db_file', 's'),
    'clip to data bounding box': ('clip', 'd'),
    'clip buffer': ('clip_buffer', 'f'),
    'use width db': ('use_width_db', 'd'),

    # CENTERLINE AND REFINEMENT PARAMETERS
    'centerline spacing': ('ds', 'float'),
    'refine centerline': ('refine_centerline', 'd'),
    'smooth parameter': ('smooth', 'f'),
    'alpha': ('alpha', 'f'),
    'maximum iterations': ('max_iter', 'd'),
    'scalar maximum width': ('scalar_max_width', 'f'),

    # ESTIMATION PARAMETERS
    'subreach size': ('subreach_size', 'f'),
    'minimum number of points to fit': ('min_fit_points', 'd'),
    'step': ('step', 'f'),
    'fit types': ('fit_types', 's'),
    'starting along track distance': ('smin', 'f'),
    'minimum observations per node': ('minobs', 'd'),

    # OUTPUT OPTIONS
    'output file': ('output_file', 's'),
}


def parse_inputs():
    """Create the argument parser."""

    parser = argparse.ArgumentParser(
        description='Estimate heigt and slope for a SWOT simulation.',
        usage=rdf_file_template)

    parser.add_argument('rdf_file', help='Input RDF file')
    parser.add_argument(
        'data_dirs',
        metavar='data_dir',
        nargs='+',
        help='Directories containing simulation data')

    parser.add_argument(
        '-f',
        '--format',
        help="OGR file format  (default 'ESRI Shapefile')",
        default='ESRI Shapefile')
    parser.add_argument(
        '--dlat',
        type=float,
        help='latitude tile size (default +1)',
        default=1.)
    parser.add_argument(
        '--dlon',
        type=float,
        help='longitude tile size (default +1)',
        default=1.)

    args = parser.parse_args()

    return args


def get_bbox_from_files(sim_files, dlat, dlon):
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

        bounding_boxes.append(
            box(lonmin, latmin, lonmin + dlon, latmin + dlat))

    return bounding_boxes


def write_catalog(output_file, format, bounding_boxes, sim_files):
    """Write coverage polygons with metadata."""

    fields = {'file': ('str', 42)}

    writer = OGRWriter(
        output_file, fields=fields, driver=format, geometry='Polygon')

    for i, bbox in enumerate(bounding_boxes):
        field_record = {'file': split(sim_files[i])[-1]}
        print(split(sim_files[i])[-1][0:7], bbox.wkt)
        writer.add_wkt_feature(bbox.wkt, field_record)
    writer.close()


def main():

    args = parse_inputs()

    # Parse the input RDF

    print(args.rdf_file)
    pars = RDF_to_class(input_vars, file=args.rdf_file)

    if type(pars.class_list) == int:
        pars.class_list = [pars.class_list]
    try:
        pars.lat_0 = float(pars.lat_0)
    except:
        pars.lat_0 = None
    try:
        pars.lon_0 = float(pars.lon_0)
    except:
        pars.lon_0 = None

    pars.fit_types = pars.fit_types.split()

    shape_file_root = join(pars.reach_db_dir, pars.shape_file_root)

    # Get a list of the files

    sim_files = []
    for directory in args.data_dirs:
        files = glob(join(directory, '*_cycle_*_pass_*.*.nc'))
        if len(files) > 0:
            sim_files += files

    # Write a catalog of the tiles processed as a GIS file

    ## bounding_boxes = get_bbox_from_files(sim_files,args.dlat,args.dlon)
    ## write_catalog(args.output_file,args.format,bounding_boxes,sim_files)

    # Open the width database, if desired

    if bool(pars.use_width_db):
        width_db_file = join(pars.width_db_dir, pars.width_db_file)
        width_db = WidthDataBase(width_db_file)
        print(('Width data base file opened: %s' % width_db_file))
    else:
        width_db = None

    # Process each input file

    for sim_file in sim_files:

        print(('Reading file: %s' % sim_file))

        # Initialize the data

        for k in pars.__dict__:
            if k not in ['d']:
                print('%s : %s' % (k, pars.__dict__[k]))

        try:
            estimator = SWOTRiverEstimator(
                sim_file,
                class_list=pars.class_list,
                lat_kwd=pars.lat_kwd,
                lon_kwd=pars.lon_kwd,
                class_kwd=pars.class_kwd,
                height_kwd=pars.height_kwd,
                true_height_kwd=pars.true_height_kwd,
                proj=pars.proj,
                x_0=pars.x_0,
                y_0=pars.y_0,
                lat_0=pars.lat_0,
                lon_0=pars.lon_0,
                ellps=pars.ellps)
            print('data read')
        except:
            print('Insufficient number of points. Continuing to next file')
            continue

        # Initialize the output file for append

        estimator.init_output_file(pars.output_file, mode='a')

        # Get the reaches corresponding to these data

        reaches = estimator.get_reaches(
            shape_file_root,
            clip=bool(pars.clip),
            clip_buffer=pars.clip_buffer)
        print(('Number of reaches read: %d' % reaches.nreaches))

        # Set the estimator width data base

        estimator.set_width_db(width_db)
        print('width data base set')

        # Process the reaches

        try:
            """
            estimator.process_reaches(use_width_db=bool(pars.use_width_db),
                                    refine_centerline=bool(pars.refine_centerline),
                                    smooth=pars.smooth,alpha=pars.alpha,
                                    max_iter=pars.max_iter,scalar_max_width=pars.scalar_max_width,
                                    subreach_size=pars.subreach_size,
                                    min_fit_points=pars.min_fit_points,
                                    step=pars.step,fit_types=pars.fit_types,
                                    ds=pars.ds,smin=pars.smin,minobs=pars.minobs,max_width=None)
                                    """
            estimator.process_reaches(
                use_width_db=bool(pars.use_width_db),
                refine_centerline=bool(pars.refine_centerline),
                smooth=pars.smooth,
                alpha=pars.alpha,
                max_iter=pars.max_iter,
                scalar_max_width=pars.scalar_max_width,
                min_fit_points=pars.min_fit_points,
                fit_types=pars.fit_types,
                ds=pars.ds,
                minobs=pars.minobs)
            print('reaches processed')
            #print estimator.river_reach_collection[0]
            #print estimator.river_reach_collection[0].h_n_ave
            #estimator.store.append('reachHeights',estimator.river_reach_collection[0].h_n_ave)
            #estimator.store.append('fits',estimator.fit_collection)
            #estimator.store.append('obs',estimator.river_obs_collection)
        except:
            estimator.store.close()
            print('HDFStore closed')
            print('Something went wrong in process reaches')
            continue

        # Close the output file HDFStore

        estimator.store.close()
        print('HDFStore closed')

    # Close the width data base before exiting

    if width_db != None:
        width_db.h5.close()

    print('Successfuly estimated river heights and slopes')


if __name__ == '__main__': main()
