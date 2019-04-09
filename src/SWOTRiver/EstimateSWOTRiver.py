"""
Given a SWOTL2 file, estimate all of the reaches within it and
output the results to a file
"""
from __future__ import absolute_import, division, print_function
import ast
import os
import SWOTRiver.SWOTRiverEstimator

rdf_file_template = """%(prog)s [-h] rdf_file

Estimate river widths, slopes and height measurements from simulated data.

positional arguments:
  rdf_file    rdf file containing the inputs.

optional arguments:
  -h, --help   show this help message and exit
  -f, --format OGR file format  (default 'ESRI Shapefile')

  A template RDF file follows (Input after ! is ignored in parsing):

! These are the locations of the reach data and the width data base (if desired)
! This is the width database (if desired)
width_db_file = None !/path/to/RiverObsTestData/GRWDL/nAmerica_GRWDL.h5

! This is the path to the reach files (not including the .shp suffix)
shape_file_root = /path/to/NA_reaches_data_discharge_depth_chn_grdc_revised_GCS

! This is the location of the simulated water file
l2_file = /path/to/swot_heights_ohio_example_v1.Multilook_L2PIXC.nc

! Where to write output files
fout_reach = /path/to/reaches (will make reaches/reaches.shp)
fout_node = /path/to/nodes (will make nodes/nodes.shp)
fout_index = /path/to/index.nc


! Bounding box may be required because the simulator sometimes has anomalous
! latitudes and longitudes. If this is not the case, it can be set to None.
! It only has to be good enough to exclude location anomalies and will be
! updated to the true bounding box by the program.
lonmin =  -83
latmin =  38
lonmax =  -82
latmax =  39
bounding_box = lonmin, latmin, lonmax, latmax

! Which dataset to use for latitudes
lat_kwd = latitude_medium

! Which dataset to use for longitudes
lon_kwd = longitude_medium

! Which dataset to use for classification
class_kwd = classification

! Which dataset to use for heights
height_kwd = height_medium

! The third set of inputs have to do with how to use the classification
! to estimate river width
! The list of classes to consider for potential inundation.
! The truth classes are [1], if no_layover_classification' is used.
! If estimated classification is used, the choice depends on whether
! use_fractional_inundation is set.
! If it is not set, either [3,4] or [4] should be used.
! If it is set, [2,3,4,5] or [3,4,5] should be used.
class_list = [2, 3, 4, 5]

! If the L2 water file has been updated to contain the fractional
! inundation, this is the name of the variable. If it has not been
! updated or you do not wish to use it, set this to None
fractional_inundation_kwd = continuous_classification

! This corresponds to the clases set above.
! If True, use fractional inundation to get the inundated area per class.
! If False, assume that class is fully flooded.
use_fractional_inundation = [True, True, False, False]

! This corresponds to the clases set above.
! if true, assume this class is water for feature segmentation purposes
use_segmentation = [False, True, True, True]

! This corresponds to the clases set above.
! if true, use this class for estimating heights
use_heights = [False, False, True, False]

! This is the minimum number of measurements that the data set must have.
min_points = 100

! The fourth set of inputs have to do with the reaches and width data base.

! The clip_buffer is a buffer (in degrees) that is drawn around the data
! bounding box so that the full reach is included and is not broken if
! the river. 0.01 ~ 1km
clip_buffer = 0.02

! The fifth set of options has to do with how the data are sampled and
! quantities are estimated

! This option is only possible if you have an a priori estimate
! of width for each width point. It will load that width into
! the centerline for comparison with the estimated data.
use_width_db =  False !True

! This option determines the separation between centerline nodes.
! If set to None, the the spacing in the input reach is used.
! The units are meters. The default is to use the input reach spacing.
ds = 300.

! The next set of options are required if one desires to refine
! the centerline if it does not align well with the data.
! If you do not know what these are, don't change them.
refine_centerline = False ! Set to True if you want to refine the centerline.
smooth = 1.e-2
alpha = 1.
max_iter = 1

! This is how far from the centerline points are accepted
scalar_max_width = 600.

! This variable states how many valid data points are required before
! a node can be consired to have a sufficient number of observations.
minobs = 10

! Set this if there seem to be issues with the first/last nodes.
! This can happen sometimes in the near range.
trim_ends = True

! These are the minimum number of points required for a slope fit
min_fit_points = 3
"""

# dict mapping input params for estimate function to RDF key/value pairs
# each key in dict is the RDF value, while the values in dict are a tuple
# with the variable name in struct as 1st item and datatype as second.
# Datatypes: 's' - string; 'd' - int; 'f' - float
input_vars = {
    'l2_file': ('l2_file', 's'),
    'fout_node': ('fout_node', 's'),
    'fout_reach': ('fout_reach', 's'),
    'fout_index': ('fout_index', 's'),
    'lonmin': ('lonmin', 'f'),
    'lonmax': ('lonmax', 'f'),
    'latmin': ('latmin', 'f'),
    'latmax': ('latmax', 'f'),
    'lat_kwd': ('lat_kwd', 's'),
    'lon_kwd': ('lon_kwd', 's'),
    'class_kwd': ('class_kwd', 's'),
    'height_kwd': ('height_kwd', 's'),
    'xtrack_kwd': ('xtrack_kwd', 's'),
    'class_list': ('class_list', 's'),
    'fractional_inundation_kwd': ('fractional_inundation_kwd', 's'),
    'use_fractional_inundation': ('use_fractional_inundation', 's'),
    'use_segmentation': ('use_segmentation', 's'),
    'use_heights': ('use_heights', 's'),
    'min_points': ('min_points', 'd'),
    'shape_file_root': ('shape_file_root', 's'),
    'clip_buffer': ('clip_buffer', 'f'),
    'use_width_db': ('use_width_db', 's'),
    'width_db_file': ('width_db_file', 's'),
    'ds': ('ds', 'f'),
    'refine_centerline': ('refine_centerline', 's'),
    'smooth': ('smooth', 'f'),
    'alpha': ('alpha', 'f'),
    'scalar_max_width': ('scalar_max_width', 'f'),
    'max_iter': ('max_iter', 'd'),
    'minobs': ('minobs', 'd'),
    'trim_ends': ('trim_ends', 's'),
    'min_fit_points': ('min_fit_points', 'd'),
    'subsample_factor': ('subsample_factor', 'd')
}


def estimate(params):
    """
    Runs RiverObs using commanded param structure of parameters either
    from RDF file or generated by another process.

    Returns iterable of reaches "reach_collection" to caller
    """
    bounding_box = params.lonmin, params.latmin, params.lonmax, params.latmax
    lat_0, lon_0 = None, None

    class_list = ast.literal_eval(params.class_list)
    use_fractional_inundation = ast.literal_eval(
        params.use_fractional_inundation)
    use_segmentation = ast.literal_eval(params.use_segmentation)
    use_heights = ast.literal_eval(params.use_heights)
    use_width_db = ast.literal_eval(params.use_width_db)
    refine_centerline = ast.literal_eval(params.refine_centerline)

    # Read the data and estimate the flooded area.
    river_estimator = SWOTRiver.SWOTRiverEstimator(
        params.l2_file,
        bounding_box=bounding_box,
        lat_kwd=params.lat_kwd,
        lon_kwd=params.lon_kwd,
        class_kwd=params.class_kwd,
        height_kwd=params.height_kwd,
        class_list=class_list,
        xtrack_kwd=params.xtrack_kwd,
        fractional_inundation_kwd=params.fractional_inundation_kwd,
        use_fractional_inundation=use_fractional_inundation,
        use_segmentation=use_segmentation,
        use_heights=use_heights,
        min_points=params.min_points,
        store_obs=False,
        store_reaches=False,
        output_file=params.fout_index,
        proj='laea',
        x_0=0,
        y_0=0,
        lat_0=lat_0,
        lon_0=lon_0,
        subsample_factor=params.subsample_factor)

    # Load the reaches and width data base

    river_estimator.get_reaches(
        params.shape_file_root, clip_buffer=params.clip_buffer)

    if use_width_db:
        river_estimator.get_width_db(params.width_db_file)

    # Process all of the reaches
    reach_collection = river_estimator.process_reaches(
        scalar_max_width=params.scalar_max_width,
        minobs=params.minobs,
        min_fit_points=params.min_fit_points,
        use_width_db=use_width_db,
        ds=params.ds,
        refine_centerline=refine_centerline,
        smooth=params.smooth,
        alpha=params.alpha,
        max_iter=params.max_iter)

    return reach_collection
