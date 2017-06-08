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
        sys.stderr.write("Libraries not found. Make sure you are running in the SWOTRiver environment.\n")
        sys.exit(1)

# Imports

import ast
from os.path import join, split
import argparse
from glob import glob
from os.path import join
from SWOTRiver import SWOTRiverEstimator
from RiverObs import RiverReachWriter

from shapely.geometry import box
from RDF import RDF_to_class, RDF

rdf_file_template = """%(prog)s [-h] rdf_file

Estimate river widths, slopes and height measurements from simulated data.

positional arguments:
  rdf_file    rdf file containing the inputs.

optional arguments:
  -h, --help   show this help message and exit
  -f, --format OGR file format  (default 'ESRI Shapefile')

  A template RDF file follows:

! These are the locations of the reach data and the width data base (if desired)

! This is the width database (if desired), no longer needed since prior database has al this info
width_db_file = ../../data/Databases/OhioRightSwath_GRWDL.h5

! This is the path to the reach files (not including the .shp suffix), this is the prior reach database
shape_file_root = ../../data/Databases/OhioRightSwath_GRWDL_river_topo/OhioRightSwath_GRWDL_river_topo

! This is the location of the simulated water file

l2_file = ../../data/OhioRightSwathData/heights/swot_intf_ohio_cycle_0001_pass_0413.RightSwath.EstClass.nc

! This is the location of the output data

fout_reach = ../../data/Results/ohio_cycle_0001_pass_0413.RightSwath.EstClass.NLoc.CClass_reach
fout_node = ../../data/Results/ohio_cycle_0001_pass_0413.RightSwath.EstClass.NLoc.CClass_reach
fout_index = ../../data/Results/ohio_cycle_0001_pass_0413.RightSwath.EstClass.NLoc.CClass_reach

! The bounding box tells what region is of interest.
! The center of this bounding box defines the x,y projection that the nodes are defined on.
! Generally, should make this to just cover the "tile" (or gdem or gdem_dem extent).
! The bounding box is no longer updated based on the data to enable processing gdem as well as pixc
! and get the same node locations

lonmin = -125.
latmin = 30.
lonmax = -10.
latmax = 50.
bounding_box = lonmin,latmin,lonmax,latmax

! The second set of inputs have to do with what data are extracted from the
! water file.

! either 'no_layover_latitude' (a priori lat) or 'latitude' (estimated lat)
lat_kwd = latitude
! either 'no_layover_longitude' (a priori lon) or 'longitude' (estimated lat)
lon_kwd = longitude
! either 'classification' (estimated classification)
! or 'no_layover_classification' (truth classification)
class_kwd = classification
! either 'height' (estimated height) or 'water_height' (truth height)
height_kwd = height

! The third set of inputs have to do with how to use the classification
! to estimate river width

! The list of classes to consider for potential inundation.
! The truth classes are [1], if no_layover_classification' is used.
! If estimated classification is used, the choice depends on whether
! use_fractional_inundation is set.
! If it is not set, either [3,4] or [4] should be used.
! If it is set, [2,3,4] or [3,4] should be used.
class_list = [2,3,4,5]

! If the L2 water file has been updated to contain the fractional
! inundation, this is the name of the variable. If it has not been
! updated or you do not wish to use it, set this to None
fractional_inundation_kwd = 'continuous_classification'

! This corresponds to the clases set above.
! If True, use fractional inundation estimate to get the inundated area for this class.
! If False, assume that this class is fully flooded.
use_fractional_inundation=[True, True, False, False]

! This corresponds to the clases set above.
! if true, assume this class if water for feature segmentation purposes
use_segmentation=[False, True, True, True]

! This corresponds to the clases set above.
! if true, assume this class is good for estimating heights
use_heights=[False, False, True, False]

! This is the minimum number of measurements that the data set must have.
min_points=100

! The fourth set of inputs have to do with the reaches and width data base.

! The clip_buffer is a buffer (in degrees) that is drawn around the data
! bounding box so that the full reach is included and is not broken if
! the river. 0.01 ~ 1km
clip_buffer=0.02


! The fifth set of options has to do with how the data are sampled and
! quantities are estimated

! This option is only possible if you have an a priori estimate
! of width for each width point. It will load that width into
! the centerline for comparison with the estimated data.
use_width_db = True

! This option determines the separation between centerline nodes.
! If set to None, the the spacing in the input reach is used.
! The units are meters. The default is to use the input reach spacing.
ds = 300.

! The next set of options are required if one desires to refine
! the centerline if it does not align well with the data.
! If you do not know what these are, don't change them.
refine_centerline=False ! Set to True if you want to refine the centerline.
smooth=1.e-2
alpha=1.
max_iter=1
! This is how far from the centerline points are accepted
scalar_max_width=600.

! This variable states how many valid data points are required before
! a node can be consired to have a sufficient number of observations.
minobs = 10

! Set this if there seem to be issues with the first/last nodes.
! This can happen sometimes in the near range.
trim_ends = True

! These are the fitting algorithms desired for mean height and slope estimation.
! More than one type of fit can be requested.
! 'OLS': ordinary least square
! 'WLS': weighted least squares
! 'RLM': Robust Linear Model
fit_types=['OLS','WLS','RLM']

! These are the minimum number of points required for a slope fit
min_fit_points = 3

!TODO: should simplify this rdf by using standard defaluts
"""

# Input directory for RDF_to_class

input_vars = {

    'l2_file' : ('l2_file','s'),

    'fout_node' : ('fout_node','s'),
    'fout_reach' : ('fout_reach','s'),
    'fout_index' : ('fout_index','s'),

    'lonmin' : ('lonmin','f'),
    'latmin' : ('latmin','f'),
    'lonmax' : ('lonmax','f'),
    'latmax' : ('latmax','f'),#'proj_lon_0' : ('proj_lon_0','f'),#'proj_lat_0' : ('proj_lat_0','f'),

    'lat_kwd' : ('lat_kwd','s'),
    'lon_kwd' : ('lon_kwd','s'),
    'class_kwd' : ('class_kwd','s'),
    'height_kwd' : ('height_kwd','s'),

    'class_list' : ('class_list','s'),
    'fractional_inundation_kwd' : ('fractional_inundation_kwd','s'),
    'use_fractional_inundation' : ('use_fractional_inundation','s'),
    'use_segmentation' : ('use_segmentation','s'),
    'use_heights' : ('use_heights','s'),
    'min_points' : ('min_points','d'),

    'shape_file_root' : ('shape_file_root','s'),
    'clip_buffer' : ('clip_buffer','f'),

    'use_width_db' : ('use_width_db','s'),
    'width_db_file' : ('width_db_file', 's'),
    'ds' : ('ds','f'),
    'refine_centerline' : ('refine_centerline','s'),
    'smooth' : ('smooth','f'),
    'alpha' : ('alpha','f'),
    'scalar_max_width' : ('scalar_max_width','f'),
    'max_iter' : ('max_iter','d'),
    'minobs' : ('minobs','d'),
    'trim_ends' : ('trim_ends','s'),
    'fit_types' : ('fit_types','s'),
    'min_fit_points' : ('min_fit_points','d'),
    }

def parse_inputs():
    """Create the argument parser."""

    parser = argparse.ArgumentParser(
        description='Estimate heigt and slope for a SWOT simulation.',
        usage=rdf_file_template)

    parser.add_argument('rdf_file',help='Input RDF file')
    parser.add_argument('-f','--format',help="OGR file format  (default 'ESRI Shapefile')",
                        default='ESRI Shapefile')

    args = parser.parse_args()

    return args


def main():

    # make sure you are running in the right environment

    search_for_libraries()

    # Parse the inputs

    args = parse_inputs()

    # Parse the input RDF

    pars = RDF_to_class(input_vars,file=args.rdf_file)

    print((pars.l2_file))
    print((pars.fout_reach))
    print((pars.fout_node))
    print((pars.fout_index))
    # Reformat some inputs

    bounding_box = pars.lonmin,pars.latmin,pars.lonmax,pars.latmax
    lat_0,lon_0=None,None#pars.proj_lat_0,pars.proj_lon_0
    class_list = ast.literal_eval(pars.class_list)
    use_fractional_inundation = ast.literal_eval(pars.use_fractional_inundation)
    use_segmentation = ast.literal_eval(pars.use_segmentation)
    use_heights = ast.literal_eval(pars.use_heights)
    use_width_db = ast.literal_eval(pars.use_width_db)
    refine_centerline = ast.literal_eval(pars.refine_centerline)
    fit_types = ast.literal_eval(pars.fit_types)

    # Read the data and estimate the flooded area.

    river_estimator = SWOTRiverEstimator(pars.l2_file,
                                        bounding_box=bounding_box,
                                        lat_kwd=pars.lat_kwd,
                                        lon_kwd=pars.lon_kwd,
                                        class_kwd=pars.class_kwd,
                                        height_kwd=pars.height_kwd,
                                        class_list=class_list,
                                        fractional_inundation_kwd=pars.fractional_inundation_kwd,
                                        use_fractional_inundation=use_fractional_inundation,
                                        use_segmentation=use_segmentation,
                                        use_heights=use_heights,
                                        min_points=pars.min_points,
                                        verbose=True,store_obs=False,
                                        store_reaches=False,
                                        store_fits=False,
                                        output_file=pars.fout_index,
                                        proj='laea',x_0=0,y_0=0,lat_0=lat_0,lon_0=lon_0)

    # Load the reaches and width data base

    river_estimator.get_reaches(pars.shape_file_root,clip_buffer=pars.clip_buffer)
    if use_width_db:
        river_estimator.get_width_db(pars.width_db_file)

    # Process all of the reaches
    river_reach_collection = river_estimator.process_reaches(scalar_max_width=pars.scalar_max_width,
                    minobs=pars.minobs,min_fit_points=pars.min_fit_points,
                    fit_types=fit_types,
                    use_width_db = use_width_db,
                    ds=pars.ds,
                    refine_centerline=refine_centerline,
                    smooth=pars.smooth,alpha=pars.alpha,max_iter=pars.max_iter)

    # Initialize the output writer
    #print river_reach_collection
    reach_output_variables = list(river_reach_collection[0].metadata.keys())
    # Brent williams May 2017: added extra fields
    node_output_variables = ['lat','lon','x','y','nobs','s',
                                 'w_ptp','w_std','w_area','w_db','area',
                                 'h_n_ave','h_n_std','h_a_ave','h_a_std',
                                 'nobs_h','x_prior','y_prior','node_indx','reach_indx']

    writer = RiverReachWriter(river_reach_collection,
                          node_output_variables,
                          reach_output_variables)

    # Write shapefiles

    driver = args.format
    writer.write_nodes_ogr(pars.fout_node,driver=driver)
    writer.write_reaches_ogr(pars.fout_reach,driver=driver)

    print('Successfuly estimated river heights and slopes')
