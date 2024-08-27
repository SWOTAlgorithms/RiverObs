"""
Given a SWOTL2 file, fit all of the reaches observed and output results.
"""

from __future__ import absolute_import, division, print_function
import os

import scipy.ndimage
import numpy as np
import netCDF4 as nc
import collections
import scipy.stats
import scipy.ndimage
import statsmodels.api
import logging
import contextlib
import warnings
import bottleneck
import piecewise_regression
from scipy.ndimage.morphology import binary_dilation

import RiverObs.ReachDatabase
import SWOTWater.aggregate
import SWOTRiver.discharge
from SWOTWater.products.constants import FILL_VALUES
from .SWOTL2 import SWOTL2
from RiverObs import WidthDataBase
from RiverObs import IteratedRiverObs
from RiverObs import RiverNode
from RiverObs import RiverReach
from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9
from SWOTRiver.errors import RiverObsException
from Centerline.Centerline import Centerline, CenterLineException

LOGGER = logging.getLogger(__name__)

REACH_WSE_SYS_UNCERT = 0.09  # m
SLOPE_SYS_UNCERT = 0.000003  # m/m

NEAR_RANGE_XTRACK_THRESHOLD = 10000
FAR_RANGE_XTRACK_TRESHOLD = 60000
WATER_FRAC_SUSPECT_THRESHOLD = 3

MIN_VALID_WSE = -500
MAX_VALID_WSE = 8000

TOO_FEW_PIXELS_THRESHOLD = 10

BLOCKING_WIDTH_THRESHOLD_FRACTION = 0.90

class SWOTRiverEstimator(SWOTL2):
    """
    Given a SWOTL2 file, fit all of the reaches observed and output results.

    This class is derived from the SWOTL2 class.

    This class contains four basic components:
    1. It is a subclass of the SWOTL2, and therefore is a LatLonRegion and
    also contains SWOT data.
    2. For each reach of interest, it keeps a dictionary RiverObs instances
    for each reach.
    3. It keeps a dictionary of RiverReach instances, containing the results
    of height/slope fitting and width estimation.

    Keeping copies for each observation, although useful for visualization,
    can be a tax on memory for a large number of reaches and observation,
    so that storing can be turned off, if so desired.

    Parameters
    ----------

    height_kwd : str
        name of netcdf variable to use as a measurement.
    true_height_kwd: str
        name of variable to use as truth for comparison.
    no_noise_height_kwd : str
        name for no noise measurement
    xtrack_kwd : str
        name of cross-track distance variable.
    bounding_box : 4-tuple or array_like
        If the bounding_box is provided, select only the data in the
        bounding box, otherwise, get the bounding box from the data
        itself. bounding_box should be of the form
        (lonmin,latmin,lonmax,latmax)
    lon_kwd, lat_kwd : str
        netcdf names for the longitudes and latitudes to be used
        for georeferencing the data set. The a priori geolocation can be
        retrieved using 'no_layover_longitude' and 'no_layover_latitude'.
        The estimated latitude and longitude can be retrieved using
        'longitude' and 'latitude', respectively.
    class_kwd : str, default 'classification'
        the netcdf name of the classification layer to use.
        'classification' will return the result of applying a classifier
        to the simulated data. 'no_layover_classification' can be used as the
        truth
    class_list : list, default [2,3,4,5,6,7]
        a list of the class labels for what is considered good water data.
        This should be defined as any classes that might contain water, even
        if only partially. If truth data are desired (class_kwd =
        'no_layover_classification'), the this should be set to [1].
        The default has interior water pixels(4), border water pixels (3),
        border land pixels (2), dark water (5), low-coherence edge water (6),
        and low-coherence interior water (7).
    fractional_inundation_kwd : str, default 'continuous_classification'
        Netcdf keyword containing the inundation fraction. If None, the no
        inundation fraction is used. The inundation fraction is, in theory,
        a number between 0 and 1 estimating the fraction of the pixel covered
        with water. In practice, because of noise, it can be outside the
        bounds and even be negative! It will produced an ensemble
        unbiased estimate if the class mean cross sections are known.
    use_fractional_inundation : bool list,
        default [True, True, False, False, False, False]
        For which classes should the inundation fraction be used instead of a
        binary classification. The default is to assume that interior pixels
        are 100% water, but to use both land and water edge pixels partially
        by using the fractional inundation kwd. Low-coherence classes (6,7) may
        have spurious water fractions and thus the fractional inundation is not
        used.
    use_segmentation :  bool list,
        default [False, True, True, True, True, True]
        Selects which classes to assume as water for segmentation purposes.
    use_heights : bool list,
    default [False, True, True, True, True, True]
        Selects which classes to use for estimating heights
        TO-DO: validate low-coherence class heights for this parameter
    min_points : int
        Minimum number of good PIXC pixels required in order for the RiverTile
        processor to output a populated RiverTile product. Default is 0.
    height_kwd : str, default 'height'
        These are the heights to extract from the water file.
    true_height_kwd : str, default 'water_height'
        These are the reference heights to use for performance assessment
    no_noise_height_kwd : str, default 'no_noise_height'
        These are the a priori heights used for phase unwrapping
    xtrack_kwd : str, 'no_layover_cross_track'
        This is the distance to the nadir track.
    xtrack_res_kwd : str, default 'flat_pixel_resolution'
        This is the cross-track dimension of the radar pixel. The netcdf file
        has an estimate 'no_layover_ground_resolution' that takes into account
        a priori surface slopes to derive an incidence angle. However, this
        runs into trouble if water
    trim_ends : bool, default False
        If True, the first and last nodes are not used, to avoid potential
        edge effects.
    store_obs : bool, default True
        If True, store each RiverObs instance in a dictionary.
    height_agg_method : str, default 'weight'
        Method for aggregating pixel heights to node. By default, pixels are
        weighted using a weighted mean using the inverse height variance from
        the pixel cloud as the weights.
    area_agg_method : str, default 'composite'
        Method for aggregating pixel-level areas to node. The default
        'composite' method uses the full pixel area for interior water and the
        water fraction for edge water pixels.
    slope_method : str, default 'weighted'
        Algorithm for computing reach-level slope. Default is 'weighted'. May
        also compute using 'first_to_last' or 'unweighted' fit to the reach.
    bayes_slope_use_all_nodes : bool, default True
        If True, add unpopulated nodes to current reach for bayes
        reconstruction.
    use_ext_dist_coef : bool, default True
        Flag for using the extreme distance coefficient to define
        search thresholds for the pixel assignment. This value is defined in
        the prior river database (SWORD) for river reaches near water bodies.
    outlier_method : str, default 'None'.
        This describes which method to use for reach-level wse outlier
        rejection. Options are 'None', 'iterative_linear' and 'piecewise_linear'
    outlier_abs_thresh : float, default 2
        The threshold for absolute error in outlier rejection. Nodes with error
        greater than this value will be detected as outliers. Default is 2 m.
    outlier_rel_thresh : float (0-100), default 68
        The threshold for relative error in outlier rejection. Nodes with error
        greater than this value will be detected as outliers. Default is 68%.
    outlier_upr_thresh : float (0-100), default 80
        The percentage of nodes that must be kept for a given reach following
        outlier detection and masking. By default, 80% of nodes must not be
        masked (i.e. the maximum proportion of nodes that may be masked is
        20%, so 80% of the reach node information is preserved).
    outlier_iter_num   : int, default 2
        The number of iterations of fitting, residual checking, and outlier
        masking that will be performed by the outlier detection algorithm.
    outlier_min_dist_btw_bps : float (0-1), default 0.1
        Minimum distance between breakpoints, as a proportion of the input data
        range for piecewise linear outlier detection algorithm.
    outlier_min_dist_to_edge : float (0-1), default 0.1
        Minimum distance from edge to a breakpoint, as a proportion of the
        input data range for piecewise linear outlier detection algorithm.
    pixc_qual_handling : str, default None
        Describes approach for pixel cloud quality handling. Options are
        'None', which ignores the pixel cloud quality data, 'nominal', which
        uses the standard quality handling described in ATBD, and
        'experimental', which is a testing state for new features.
    geo_qual_wse_suspect: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC geolocation qual map
        to a suspect wse measurement.
    geo_qual_wse_degraded: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC geolocation qual map
        to a degraded wse measurement.
    geo_qual_wse_bad: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC geolocation qual map
        to a bad wse measurement.
    class_qual_area_suspect: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC classification qual
        map to a suspect area measurement.
    class_qual_area_degraded: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC classification qual
        map to a degraded area measurement.
    class_qual_area_bad: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC classification qual
        map to a bad area measurement.
    sig0_qual_suspect: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC sig0 qual map to a
        suspect sig0 measurement.
    sig0_qual_bad: 32 bit hex, default None
        Bitwise word describing which bits from the PIXC sig0 qual map to a
        bad sig0 measurement.
    num_good_sus_pix_thresh_wse: Minimum number of good and suspect pixels
        required to use only good/suspect data in node aggregation for wse.
    num_good_sus_pix_thresh_area: Minimum number of good and suspect pixels
        required to use only good/suspect data in node aggregation for area.
    use_bright_land: Boolean value, if True we use bright land in pixel
        assignment and aggregations (mark as suspect), if False do not use
        in pixel assignment.

    The final set of keywords are projection options for pyproj. See Notes.
    Notes
    -----

    A full list of projection options to set plus explanations of their
    meaning can be found here: https://trac.osgeo.org/proj/wiki/GenParms

    The default projection is Lambert equiarea, which has a proj4 string with
    the following parameters:

    +proj=laea
    +lat_0=Latitude at projection center, set to bounding box center lat
    +lon_0=Longitude at projection center, set to bounding box center lon
    +x_0=False Easting, set to 0
    +y_0=False Northing, set to 0
    """

    xtrack_res_max = 200.  # Maximum allowed cross-track range resolution
    platform_height = 970.e3
    earth_radius = 6378e3

    def __init__(self,
                 swotL2_file,
                 bounding_box=None,
                 lat_kwd='no_layover_latitude',
                 lon_kwd='no_layover_longitude',
                 class_list=[2, 3, 4, 5, 6, 7],
                 class_kwd='classification',
                 rngidx_kwd='range_index',
                 aziidx_kwd='azimuth_index',
                 fractional_inundation_kwd='water_frac',
                 fractional_inundation_uncert_kwd='water_frac_uncert',
                 use_fractional_inundation=[True, True, False,
                                            False, False, False],
                 use_segmentation=[False, True, True, True, True, True],
                 use_heights=[False, False, True, False, True, True],
                 min_points=0,
                 height_kwd='height',
                 trim_ends=False,
                 store_obs=True,
                 xtrack_kwd='no_layover_cross_track',
                 sig0_kwd='sig0',
                 sig0_uncert_kwd='sig0_uncert',
                 ifgram_kwd='interferogram',
                 power1_kwd='power_plus_y',
                 power2_kwd='power_minus_y',
                 phase_noise_std_kwd='phase_noise_std',
                 dh_dphi_kwd='dheight_dphase',
                 dlat_dphi_kwd='dlatitude_dphase',
                 dlon_dphi_kwd='dlongitude_dphase',
                 num_rare_looks_kwd='eff_num_rare_looks',
                 num_med_looks_kwd='eff_num_medium_looks',
                 looks_to_efflooks_kwd='looks_to_efflooks',
                 false_detection_rate_kwd='false_detection_rate',
                 missed_detection_rate_kwd='missed_detection_rate',
                 darea_dheight_kwd = 'darea_dheight',
                 geoid_kwd='geoid',
                 solid_earth_tide_kwd='solid_earth_tide',
                 load_tide_fes_kwd='load_tide_fes',
                 load_tide_got_kwd='load_tide_got',
                 pole_tide_kwd='pole_tide',
                 bright_land_flag_kwd='bright_land_flag',
                 classification_qual_kwd='classification_qual',
                 geolocation_qual_kwd='geolocation_qual',
                 sig0_qual_kwd='sig0_qual',
                 layover_impact_kwd='layover_impact',
                 proj='laea',
                 x_0=0,
                 y_0=0,
                 lat_0=None,
                 lon_0=None,
                 ellps='WGS84',
                 output_file=None,
                 subsample_factor=1,
                 height_agg_method='weight',  # [weight, median, uniform, orig]
                 area_agg_method='composite',
                 preseg_dilation_iter=0,
                 slope_method='bayes',
                 bayes_slope_use_all_nodes=True,
                 prior_unc_alpha=1.5,
                 char_length_tau=10000,
                 prior_wse_method= 'fit',
                 use_multiple_reaches=False,
                 use_ext_dist_coef=True,
                 outlier_method=None,
                 outlier_abs_thresh=None,
                 outlier_rel_thresh=None,
                 outlier_upr_thresh=None,
                 outlier_iter_num=None,
                 outlier_breakpoint_min_dist=None,
                 outlier_edge_min_dist=None,
                 outlier_n_boot=None,
                 pixc_qual_handling=None,
                 geo_qual_wse_suspect=None,
                 geo_qual_wse_degraded=None,
                 geo_qual_wse_bad=None,
                 class_qual_area_suspect=None,
                 class_qual_area_degraded=None,
                 class_qual_area_bad=None,
                 sig0_qual_suspect=None,
                 sig0_qual_bad=None,
                 num_good_sus_pix_thresh_wse=None,
                 num_good_sus_pix_thresh_area=None,
                 use_bright_land=None,
                 reach_pct_good_sus_thresh=None,
                 **proj_kwds):

        self.trim_ends = trim_ends
        self.input_file = os.path.split(swotL2_file)[-1]
        self.output_file = output_file  # index file
        self.subsample_factor = subsample_factor
        self.slope_method = slope_method
        self.bayes_slope_use_all_nodes = bayes_slope_use_all_nodes
        self.prior_unc_alpha = prior_unc_alpha
        self.char_length_tau = char_length_tau
        self.prior_wse_method = prior_wse_method
        self.use_multiple_reaches = use_multiple_reaches
        self.use_ext_dist_coef = use_ext_dist_coef
        self.outlier_method = outlier_method
        self.outlier_abs_thresh = outlier_abs_thresh
        self.outlier_rel_thresh = outlier_rel_thresh
        self.outlier_upr_thresh = outlier_upr_thresh
        self.outlier_iter_num = outlier_iter_num
        self.outlier_breakpoint_min_dist = outlier_breakpoint_min_dist
        self.outlier_edge_min_dist = outlier_edge_min_dist
        self.outlier_n_boot = outlier_n_boot
        self.num_good_sus_pix_thresh_wse = num_good_sus_pix_thresh_wse
        self.num_good_sus_pix_thresh_area = num_good_sus_pix_thresh_area
        self.use_bright_land = use_bright_land
        self.reach_pct_good_sus_thresh = reach_pct_good_sus_thresh

        # Classification inputs
        self.class_kwd = class_kwd
        self.class_list = class_list
        self.fractional_inundation_kwd = fractional_inundation_kwd
        self.use_fractional_inundation = use_fractional_inundation

        self.use_segmentation = use_segmentation
        self.use_heights = use_heights
        self.height_agg_method = height_agg_method
        self.area_agg_method = area_agg_method
        self.preseg_dilation_iter = preseg_dilation_iter

        # Initialize the base class
        SWOTL2.__init__(
            self,
            swotL2_file,
            bounding_box=bounding_box,
            class_list=class_list,
            lat_kwd=lat_kwd,
            lon_kwd=lon_kwd,
            class_kwd=class_kwd,
            rngidx_kwd=rngidx_kwd,
            aziidx_kwd=aziidx_kwd,
            min_points=min_points,
            proj=proj,
            x_0=x_0,
            y_0=y_0,
            lat_0=lat_0,
            lon_0=lon_0,
            ellps=ellps,
            subsample_factor=subsample_factor,
            **proj_kwds)

        if self.output_file is not None:
            self.create_index_file()

        if np.ma.is_masked(self.lat):
            mask = self.lat.mask
        else:
            mask = np.zeros(len(self.lat), dtype=bool)
        self.h_noise = self.get(height_kwd)
        if np.ma.is_masked(self.h_noise):
            mask = mask | self.h_noise.mask

        if pixc_qual_handling not in [None, 'nominal', 'experimental']:
            raise Exception('AUX Param key: pixc_qual_handling has unexpected '+
                            'value: %s'%pixc_qual_handling)

        # set quality masks from input PIXC quality flags
        try:
            classification_qual = self.get(classification_qual_kwd)
            geolocation_qual = self.get(geolocation_qual_kwd)
            bright_land_flag = self.get(bright_land_flag_kwd)
        except KeyError:
            LOGGER.warning("PIXC Quality flags not found in input file")
            classification_qual = np.zeros(mask.shape, dtype='u4')
            geolocation_qual = np.zeros(mask.shape, dtype='u4')
            bright_land_flag = np.zeros(mask.shape, dtype='u4')

        # Update mask of pixels to reject from pixel assignment with
        # commanded bitmasks from aux param file applied to classification_qual
        # and geolocation_qual.
        mask = np.logical_or.reduce([
            mask, class_qual_area_bad & classification_qual > 0,
            geo_qual_wse_bad & geolocation_qual > 0])

        # Remove degraded class_qual where there is no prior water
        PIXC_CLASS_QUAL_NO_PRIOR = 2**2
        degraded_class_mask = class_qual_area_degraded \
                              & classification_qual > 0
        no_prior_water_mask = classification_qual \
                              & PIXC_CLASS_QUAL_NO_PRIOR > 0
        degraded_and_no_prior_water = degraded_class_mask & no_prior_water_mask
        degraded_and_prior_water = degraded_class_mask & ~no_prior_water_mask
        # uncomment to make make degraded class_qual good anywhere there is
        # prior water
        # degraded_class_mask[degraded_and_prior_water] = False
        mask = np.logical_or(mask, degraded_and_no_prior_water)

        # HACK in avoidance of bad water_frac values
        water_frac = self.get(fractional_inundation_kwd)
        mask = np.logical_or(mask, water_frac.mask)

        # Exclude bright land as well if self.use_bright_land is False
        if not self.use_bright_land:
            mask = np.logical_or(mask, bright_land_flag > 0)

        # skip NaNs in dheight_dphase
        good = ~mask
        if good.sum() == 0:
            LOGGER.warning("No usable pixels found in input PIXC file")
            raise RiverObsException(
                "No usable pixels found in input PIXC file")
        for key in ['lat', 'lon', 'x', 'y', 'klass', 'h_noise',
                    'img_x', 'img_y']:
            try:
                setattr(self, key, getattr(self, key)[good])
            except TypeError:
                pass

        # modify self.index for subsetting of good pixels
        indicies = np.where(self.index)[0]
        indicies = indicies[good]
        self.index[:] = False
        self.index[indicies] = True

        datasets2load = [
            ['xtrack', xtrack_kwd], ['sig0', sig0_kwd],
            ['sig0_uncert', sig0_uncert_kwd],
            ['water_frac', fractional_inundation_kwd],
            ['water_frac_uncert', fractional_inundation_uncert_kwd],
            ['ifgram', ifgram_kwd],
            ['power1', power1_kwd],
            ['power2', power2_kwd],
            ['phase_noise_std', phase_noise_std_kwd],
            ['dh_dphi', dh_dphi_kwd],
            ['dlat_dphi', dlat_dphi_kwd],
            ['dlon_dphi', dlon_dphi_kwd],
            ['num_rare_looks', num_rare_looks_kwd],
            ['num_med_looks', num_med_looks_kwd],
            ['false_detection_rate', false_detection_rate_kwd],
            ['missed_detection_rate', missed_detection_rate_kwd],
            ['darea_dheight', darea_dheight_kwd],
            ['geoid', geoid_kwd], ['solid_earth_tide', solid_earth_tide_kwd],
            ['load_tide_fes', load_tide_fes_kwd],
            ['load_tide_got', load_tide_got_kwd],
            ['pole_tide', pole_tide_kwd],
            ['layover_impact', layover_impact_kwd],
            ['geolocation_qual', geolocation_qual_kwd],
            ['classification_qual', classification_qual_kwd],
            ['sig0_qual', sig0_qual_kwd],
            ['bright_land_flag', bright_land_flag_kwd],
        ]

        for dset_name, keyword in datasets2load:
            try:
                value = self.get(keyword)
                # hack ifgram re/im parts into complex dtype
                if dset_name == 'ifgram':
                    value = value[:,0] + 1j* value[:,1]
            except KeyError:
                value = None
            setattr(self, dset_name, value)

        # Set PIXC quality masks for use in aggregation using commanded bitmasks
        # from aux param file. Set as attributes of self so they may be used
        # later in node aggregation.
        PIXC_GEO_QUAL_XOVR_SUSPECT = 2**6
        PIXC_GEO_QUAL_XOVR_BAD = 2**23
        PIXC_CLASS_QUAL_RINGING = 2**19

        if self.classification_qual is not None:
            self.is_area_degraded = degraded_class_mask[good]
            self.is_area_suspect = (
                    class_qual_area_suspect & self.classification_qual > 0)
            self.sring_flg = (
                self.classification_qual & PIXC_CLASS_QUAL_RINGING) > 0
        else:
            self.is_area_degraded = np.zeros(self.lat.shape, dtype='bool')
            self.is_area_suspect = np.zeros(self.lat.shape, dtype='bool')
            self.sring_flg = np.zeros(self.lat.shape, dtype='bool')

        if self.geolocation_qual is not None:
            self.is_wse_degraded = (
                geo_qual_wse_degraded & self.geolocation_qual > 0)
            self.is_wse_suspect = (
                geo_qual_wse_suspect & self.geolocation_qual > 0)
            self.is_xovercal_suspect = (
                self.geolocation_qual & PIXC_GEO_QUAL_XOVR_SUSPECT) > 0
            self.is_xovercal_degraded = (
                self.geolocation_qual & PIXC_GEO_QUAL_XOVR_BAD) > 0
        else:
            self.is_wse_degraded = np.zeros(self.lat.shape, dtype='bool')
            self.is_wse_suspect = np.zeros(self.lat.shape, dtype='bool')
            self.is_xovercal_suspect = np.zeros(self.lat.shape, dtype='bool')
            self.is_xovercal_degraded = np.zeros(self.lat.shape, dtype='bool')

        if self.sig0_qual is not None:
            self.is_sig0_bad = sig0_qual_bad & self.sig0_qual > 0
            self.is_sig0_suspect = sig0_qual_suspect & self.sig0_qual > 0
        else:
            self.is_sig0_bad = np.zeros(self.lat.shape, dtype='bool')
            self.is_sig0_suspect = np.zeros(self.lat.shape, dtype='bool')

        try:
            self.looks_to_efflooks = self.getatt(looks_to_efflooks_kwd)
            if self.looks_to_efflooks == 'None':
                self.looks_to_efflooks = 1.75  # set to default value
        except (AttributeError, KeyError):
            self.looks_to_efflooks = None

        # Try to read the pixel area from the L2 file, or compute it
        # from look angle and azimuth spacing, or from azimuth spacing
        # and ground spacing
        try:
            # hopefully already there
            self.pixel_area = self.get('pixel_area')

        except KeyError:
            try:
                # try compute with look angle
                look_angle = self.get('no_layover_look_angle')
                incidence_angle = (look_angle) * (
                    1. + self.platform_height / self.earth_radius)
                range_resolution = float(self.getatt('range_resolution'))
                azimuth_spacing = float(self.getatt('azimuth_spacing'))
                self.xtrack_res = range_resolution / np.sin(
                    np.radians(incidence_angle))
                self.pixel_area = azimuth_spacing * self.xtrack_res

            except KeyError:
                try:
                    # try compute azi / ground spacing (assume GDEM)
                    azimuth_spacing = float(self.getatt('azimuth_spacing'))
                    ground_spacing = float(self.getatt('ground_spacing'))
                    self.pixel_area = (
                        azimuth_spacing * ground_spacing *
                        np.ones(np.shape(self.h_noise)))

                except AttributeError:
                    self.pixel_area = 10.0 * np.zeros(len(self.h_noise))
                    LOGGER.warning(
                        "could not find correct pixel area parameters")

        # need to scale pixel area by the subsampling factor if subsampling
        if self.subsample_factor > 1:
            self.pixel_area = self.pixel_area * self.subsample_factor

        if fractional_inundation_kwd is None:  # all water pixels are inundated
            self.fractional_inundation = None
            self.inundated_area = self.pixel_area

        else:
            self.fractional_inundation = self.get(fractional_inundation_kwd)
            self.inundated_area = self.pixel_area.copy()
            for i, k in enumerate(class_list):
                if use_fractional_inundation[i]:
                    # TO-DO: validate fractional inundation result here
                    index = self.klass == k
                    self.inundated_area[index] = (self.pixel_area[
                        index] * self.fractional_inundation[index])

        if (self.use_segmentation is not None and
            len(self.use_segmentation) == len(self.class_list)):
            # Segment the image into disjoint features
            self.isWater = np.zeros(np.shape(self.h_noise))
            for i, k in enumerate(class_list):
                if self.use_segmentation[i]:
                    index = self.klass == k
                    self.isWater[index] = 1
            self.segment_water_class(self.preseg_dilation_iter)

        else:
            self.seg_label = None

        if len(self.use_heights) == len(self.class_list):
            self.wse_class_flg = np.zeros(np.shape(self.h_noise), dtype='bool')
            for i, k in enumerate(class_list):
                if self.use_heights[i]:
                    index = self.klass == k
                    self.wse_class_flg[index] = True

        else:
            self.wse_class_flg = None

        # Use wse_class_flg and is_wse_degraded to make self.h_flg
        if self.wse_class_flg is None:
            self.h_flg = np.logical_not(self.is_wse_degraded)
        else:
            self.h_flg = np.logical_and(
                self.wse_class_flg, np.logical_not(self.is_wse_degraded))

        # Set area / sig0 flg based on degraded area / bad sig0
        self.area_flg = np.logical_not(self.is_area_degraded)
        self.sig0_flg = np.logical_not(self.is_sig0_bad)

        LOGGER.debug('Data loaded')

        # Initialize the list of reaches to output to product
        self.river_reach_collection = collections.OrderedDict()

    def flatten_interferogram(self):
        """Flattens the pixel cloud interferogram"""
        # range index is self.img_x, azi is self.img_y
        LOGGER.info('flatten_interferogram')
        try:
            import cnes.modules.geoloc.lib.geoloc as geoloc
            from cnes.common.lib.my_variables import \
                GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE

        except ImportError:
            LOGGER.warning(
                "Can't flatten interferogram as swotCNES not found.")
            return

        try:
            tvp_plus_y_antenna_xyz = (
                self.nc.groups['tvp']['plus_y_antenna_x'][:],
                self.nc.groups['tvp']['plus_y_antenna_y'][:],
                self.nc.groups['tvp']['plus_y_antenna_z'][:])
            tvp_minus_y_antenna_xyz = (
                self.nc.groups['tvp']['minus_y_antenna_x'][:],
                self.nc.groups['tvp']['minus_y_antenna_y'][:],
                self.nc.groups['tvp']['minus_y_antenna_z'][:])
            pixc = {'tvp': {'time': self.nc.groups['tvp']['time'][:]},
                    'pixel_cloud': {'illumination_time': 
                     self.nc.groups['pixel_cloud']['illumination_time'][:]}}
        except KeyError:
            LOGGER.warning(
                "Can't flatten interferogram as TVP not in pixel cloud.")
            return

        hgt_2d = np.nan*np.ones((np.max(self.img_y)+1, np.max(self.img_x)+1))
        hgt_2d[self.img_y, self.img_x] = self.h_noise

        # do a 2d median filter on the heights
        # following line raises a Warning (all-NaN slice encountered)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hgt_filt = scipy.ndimage.generic_filter(
                hgt_2d, bottleneck.nanmedian, size=11)

        target_xyz = geoloc.convert_llh2ecef(
            self.lat, self.lon, hgt_filt[self.img_y, self.img_x],
            GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE)

        pixc_tvp_index = SWOTWater.aggregate.get_sensor_index(pixc)
        pixc_wavelength = self.nc.wavelength
        flat_ifgram = SWOTWater.aggregate.flatten_interferogram(
            self.ifgram, tvp_plus_y_antenna_xyz, tvp_minus_y_antenna_xyz,
            target_xyz, pixc_tvp_index[self.index], pixc_wavelength)

        self.ifgram = flat_ifgram

    def segment_water_class(self, preseg_dilation_iter=0):
        """
        do image segmentation algorithm on the water class to label
        unconnected features
        """
        maxX = np.max(self.img_x)
        maxY = np.max(self.img_y)
        cls_img = np.zeros((maxY + 1, maxX + 1))
        cls_img[self.img_y, self.img_x] = self.isWater

        # Do some regularization with morphological operations
        # so that water features very close to each other
        # (with 2*preseg_dilation_iter or fewer land pixels separating them)
        # are given same label
        if preseg_dilation_iter > 0:
            cls_tmp = np.zeros((maxY + 1, maxX + 1))
            cls_tmp[self.img_y, self.img_x] = 1
            cls_tmp = binary_dilation(cls_tmp, iterations=preseg_dilation_iter)
            cls_img[cls_tmp == 1] = 1

        # segment the water class image
        lbl, nlbl = scipy.ndimage.label(cls_img)

        # assign land edge segments (label 0) to nearest water feature
        # (label >0) using grey dilation for this, probably need to re-think
        # how to handle land edge pixels that touch two different
        # segments (this will arbitrarily assign it to the one with
        # the largest label index)
        lbl2 = scipy.ndimage.grey_dilation(lbl, 3)
        lbl_out = lbl.copy()
        lbl_out[lbl == 0] = lbl2[lbl == 0]

        # write out segmentation images to a nc file for testing
        dumpLbl = False
        if dumpLbl:
            M0, M1 = np.shape(lbl)
            with nc.Dataset("segDump.nc", 'w') as f:
                f.createDimension('azimuth', M0)
                f.createDimension('range', M1)
                ol = f.createVariable(
                    'orig_label', 'i4', ('azimuth', 'range'), fill_value=-128)
                fl = f.createVariable(
                    'final_label', 'i4', ('azimuth', 'range'), fill_value=-128)
                ol[:] = lbl[:]
                fl[:] = lbl_out[:]

        # create the segmentation label variable
        self.seg_label = lbl_out[self.img_y, self.img_x]

    def get_reaches(self, reach_db_path, day_of_year=None):
        """Get all of the reaches using a ReachExtractor."""
        LOGGER.debug('get_reaches')
        self.reaches = RiverObs.ReachDatabase.ReachExtractor(
            reach_db_path, self, day_of_year=day_of_year)

        return self.reaches

    def get_width_db(self, width_db_file):
        """Open the width data base for later use."""
        self.width_db = WidthDataBase(width_db_file)

    def set_width_db(self, width_db):
        """Set width data base from an already opened version."""
        self.width_db = width_db

    def get_max_width_from_db(self, reach_idx):
        """
        Get the width associated with a given reach from the width data base.
        """
        return self.width_db.get_river(
            reach_idx,
            columns=['width'],
            asarray=True,
            transpose=False,
            bounding_box=self.bounding_box,
            clip_buffer=self.clip_buffer)

    def process_reaches(self,
                        scalar_max_width=600.,
                        minobs=10,
                        min_fit_points=3,
                        ds=None,
                        refine_centerline=False,
                        smooth=1.e-2,
                        alpha=1.,
                        max_iter=1,
                        max_window_size=10000,
                        min_sigma=1000,
                        window_size_sigma_ratio=5,
                        enhanced=False):
        """
        Process all of the reaches in the data bounding box.

        Parameters
        ----------

        scalar_max_width : float, default 600
            How far away to look for points
        minobs : int, default 10
            Minimum number of observations for valid node.
        min_fit_points : int, default 3
            Minimum number of populated nodes required for height/slope fit
        ds : float, optional
            Separation between centerline nodes (in m). If None, uses reach
            points.
        refine_centerline: bool, default True
            Refine the centerline?
        smooth : float, default 1.e-2
            Centerline smoothing constant (see Centerline)
        alpha : float, default 1
            Centerline refinement update weight
        max_iter : int, default 1
            Maximum number of centerline iterations
        max_window_size : max window for gaussian averaging, default is 10km
        min_sigma : min sigma for gaussian averaging, default is 1km
        window_size_sigma_ratio : default is 5

        Returns
        -------

        Returns a list containing a RiverReach instance for each reach in the
        bounding box.
        """
        LOGGER.info('process_reaches: assigning reaches')
        # assign the reaches
        if self.use_ext_dist_coef:
            river_obs_list, reach_idx_list, ireach_list = \
                self.assign_reaches_ext_dist_coef(
                    scalar_max_width, minobs, ds)
        else:
            river_obs_list, reach_idx_list, ireach_list = \
                self.assign_reaches(scalar_max_width, minobs, ds)


        LOGGER.info('process_reaches: processing nodes')
        river_reach_collection_no_outlier_flag = []
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:

            try:
                # probably should scale this to look some fraction
                # farther than the database width
                max_width = (
                    self.reaches[ireach].metadata['Wmean'] * np.ones(
                        np.shape(self.reaches[i_reach].x)))  # *2.0

            except KeyError:
                max_width = scalar_max_width

            # Ugly way process_reach/process_node uses the data
            self.river_obs = river_obs

            river_reach = self.process_node(
                self.reaches[ireach],
                ireach,
                reach_idx,
                scalar_max_width=scalar_max_width,
                max_width=max_width,
                ds=ds,
                refine_centerline=refine_centerline,
                smooth=smooth,
                alpha=alpha,
                max_iter=max_iter)

            river_reach_collection_no_outlier_flag.append(river_reach)

        # Now iterate over reaches again and do node outlier computations and
        # add node quality flags to river_reach
        reach_zips = zip(river_reach_collection_no_outlier_flag, ireach_list)
        river_reach_collection = []
        for river_reach, ireach in reach_zips:
            # compute node mask to be used for reach aggregation
            if self.use_multiple_reaches:
                (wse, wse_r_u, node_q, this_len, first_node, ss,
                 mask, p_wse) = self.get_multi_reach_height(
                    river_reach_collection_no_outlier_flag, river_reach,
                    ireach)
                ww = 1 / (wse_r_u ** 2)
                reach_mask, valid_wse_mask, wse_outlier_mask = (
                    self.get_reach_mask(ss, wse, ww, node_q, min_fit_points,
                                        first_node, this_len))
                river_reach.mask = reach_mask[first_node:first_node+this_len][
                    river_reach.populated_nodes]
                river_reach.valid_wse_mask = valid_wse_mask[
                                             first_node:first_node+this_len][
                    river_reach.populated_nodes]
                wse_outlier_mask = \
                wse_outlier_mask[first_node:first_node + this_len][
                    river_reach.populated_nodes]

            else:
                ww = 1 / (river_reach.wse_r_u ** 2)
                river_reach.mask, river_reach.valid_wse_mask, wse_outlier_mask \
                    = self.get_reach_mask(river_reach.node_ss, river_reach.wse,
                                        ww, river_reach.node_q,
                                        min_fit_points)

            # Use node mask to set wse_outlier bit in node_q_b
            river_reach.node_q_b[~wse_outlier_mask] |= (
                SWOTRiver.products.rivertile.QUAL_IND_WSE_OUTLIER)
            # update summary qual flag to include outlier nodes
            river_reach.node_q[~wse_outlier_mask] = 3

            river_reach_collection.append(river_reach)

        LOGGER.info('process_reaches: processing reaches')
        out_river_reach_collection = []
        # Now iterate over reaches again and do reach average computations
        reach_zips = zip(
            river_reach_collection, river_obs_list, reach_idx_list,
            ireach_list)
        for river_reach, river_obs, reach_idx, ireach in reach_zips:

            # Ugly way process_reach/process_node uses the data
            self.river_obs = river_obs
            out_river_reach = self.process_reach(
                river_reach_collection, river_reach, self.reaches[ireach],
                ireach, reach_idx, min_fit_points=min_fit_points
            )

            if out_river_reach is not None:
                if enhanced:
                    slp2 = self.compute_enhanced_slope(
                        river_reach_collection, river_reach, ireach,
                        max_window_size=max_window_size,
                        min_sigma=min_sigma,
                        window_size_sigma_ratio=window_size_sigma_ratio,
                        min_fit_points=min_fit_points
                    )
                    enhanced_slope, enhanced_slope_r_u, enhanced_slope_u = slp2
                else:
                    enhanced_slope = MISSING_VALUE_FLT
                    enhanced_slope_r_u = MISSING_VALUE_FLT
                    enhanced_slope_u = MISSING_VALUE_FLT

                # add enhanced slope to river reach outputs
                out_river_reach.metadata['slope2'] = enhanced_slope
                out_river_reach.metadata['slope2_r_u'] = enhanced_slope_r_u
                out_river_reach.metadata['slope2_u'] = enhanced_slope_u
                out_river_reach_collection.append(out_river_reach)

        return out_river_reach_collection

    def assign_reaches(self,
                       scalar_max_width,
                       minobs=10,
                       ds=None,
                       use_ext_dist_coef=False):
        """
        Assigns pixels to nodes for every reach.
        """
        LOGGER.info("assign_reaches")
        # First extract the segmentation labels to keep
        all_dominant_labels = []
        all_ids = []
        all_up_ids = []
        all_dn_ids = []
        for i_reach, reach_idx in enumerate(self.reaches.reach_idx):
            if len(self.reaches[i_reach].x) <= 3:
                continue
            try:
                river_obs = IteratedRiverObs(
                    self.reaches[i_reach], self.x, self.y, ds=ds,
                    seg_label=self.seg_label, max_width=scalar_max_width,
                    minobs=minobs, use_ext_dist_coef=use_ext_dist_coef)
            except CenterLineException as e:
                print("CenterLineException: ", e)
                continue

            if river_obs.dominant_label is not None:
                all_dominant_labels.append(river_obs.dominant_label)
                all_ids.append(reach_idx)
                all_up_ids.append(
                    self.reaches[i_reach].metadata['rch_id_up'][:,0])
                all_dn_ids.append(
                    self.reaches[i_reach].metadata['rch_id_dn'][:,0])

        # Iterate over reaches, assign pixels to nodes
        river_obs_list = []
        reach_idx_list = []
        ireach_list = []
        node_x_list = []
        node_y_list = []
        for i_reach, reach_idx in enumerate(self.reaches.reach_idx):

            LOGGER.debug('assign_reaches: reach %d/%d reach index: %d' %(
                i_reach + 1, self.reaches.nreaches, reach_idx))

            try:
                seg_label = self.seg_label.copy()
            except AttributeError:
                # if self.seg_label is None
                seg_label = None

            # save node x and y in list for min-dist calculation later
            node_x_list.extend(self.reaches[i_reach].x.data.tolist())
            node_y_list.extend(self.reaches[i_reach].y.data.tolist())

            try:
                this_idx = np.where(reach_idx == all_ids)[0][0]
                adjacent_ids = np.concatenate([
                    all_up_ids[this_idx], all_dn_ids[this_idx]])
                adjacent_ids = adjacent_ids[adjacent_ids != 0]

                this_label = all_dominant_labels[this_idx]
                for that_label, that_id in zip(all_dominant_labels, all_ids):
                    if that_id in adjacent_ids:
                        seg_label[self.seg_label == that_label] = this_label

            except IndexError:
                pass

            try:
                river_obs = IteratedRiverObs(
                    self.reaches[i_reach],
                    self.x,
                    self.y,
                    ds=ds,
                    seg_label=seg_label,
                    max_width=scalar_max_width,
                    minobs=minobs,
                    use_ext_dist_coef=use_ext_dist_coef)

            except CenterLineException as e:
                print("CenterLineException: ", e)
                continue

            # Add width per node to centerline and re-init IteratedRiverObs.
            # Search width is the width it uses to always include (1/2 on
            # each side of centerline).
            
            # set the best initial search width for each node/reach
            primary_width = np.maximum(self.reaches[i_reach].max_width,
                                       self.reaches[i_reach].width)
            # fix any potential NaN or zero values from PRD
            primary_width[np.isnan(primary_width)] = np.nanmedian(primary_width)
            primary_width[primary_width <= 0] = np.nanmedian(primary_width)
            search_width = 2 * (
                primary_width * self.reaches[i_reach].wth_coef)
            river_obs.add_centerline_obs(
                self.reaches[i_reach].x, self.reaches[i_reach].y,
                search_width, 'max_width')
            river_obs.reinitialize()

            if len(river_obs.x) == 0:
                LOGGER.debug(('assign_reaches: no observations mapped to nodes '
                              'in this reach'))
                continue

            LOGGER.debug(
                'assign_reaches: {} observations mapped to nodes in this '
                'reach'.format(len(river_obs.x)))
            river_obs_list.append(river_obs)
            reach_idx_list.append(reach_idx)
            ireach_list.append(i_reach)

        # Ensure unique and optimal assignments of pixels to reach.
        tile_centerline = Centerline(node_x_list, node_y_list, k=3)
        tile_i, tile_d, tile_x, tile_y, tile_s, tile_n = \
            tile_centerline.to_centerline(self.x, self.y)
        # set min distance target to smallest distance-to-node for whole tile
        min_dist = tile_d
        reach_ind = -1 * np.ones(self.x.shape, dtype=int)
        cnts_assigned = np.zeros(self.x.shape, dtype=int)
        for ii, river_obs in enumerate(river_obs_list):
            # Get current reach assignment and min distance to node for all
            # pixels assigned to this reach.
            these_reach_inds = reach_ind[river_obs.in_channel]
            these_min_dists = min_dist[river_obs.in_channel]

            # Figure out which ones are better than current assignment
            mask = river_obs.d <= these_min_dists

            # Re-assign the pixels to reaches with a better assignment
            these_reach_inds[mask] = ii
            these_min_dists[mask] = river_obs.d[mask]
            reach_ind[river_obs.in_channel] = these_reach_inds
            min_dist[river_obs.in_channel] = these_min_dists
            cnts_assigned[river_obs.in_channel] += 1

        # iterate over river_obs again to set it so optimized
        for ii, river_obs in enumerate(river_obs_list):
            mask_keep = reach_ind[river_obs.in_channel] == ii
            river_obs.mask_pixels(mask_keep)

        # Iterate through and only return reaches with pixels in them.
        # (don't iterate and modify!)
        river_obs_list_out = []
        reach_idx_list_out = []
        ireach_list_out = []
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:

            # Filter out ghost reaches and reaches with only one node
            if reach_idx % 10 == 6 or len(self.reaches[ireach].x) == 1:
                continue

            if len(river_obs.populated_nodes) > 0:
                river_obs_list_out.append(river_obs)
                reach_idx_list_out.append(reach_idx)
                ireach_list_out.append(ireach)

        return river_obs_list_out, reach_idx_list_out, ireach_list_out

    def assign_reaches_ext_dist_coef(
            self, scalar_max_width, minobs=10, ds=None):
        """
        Does the reach assignments using ext_dist_coef and the outputs from
        assign_reaches.
        """
        LOGGER.info("assign_reaches_ext_dist_coef")
        # Note we do the ext_dist_coef processing out of assign_reaches to
        # get the desired behavior
        river_obs_list, reach_idx_list, ireach_list = self.assign_reaches(
            scalar_max_width, minobs, ds)

        # Iterate through and only return reaches with pixels in them.
        river_obs_list_out = []
        reach_idx_list_out = []
        ireach_list_out = []
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:
            # Only keep pixels which were assigned with the larger extreme
            # distance, and which are also in the smaller extreme distance
            # when scaled by ext_dist_coef. This code replicates the logic
            # in RiverObs.flag_out_channel_and_label method.

            # set the best initial search width for each node/reach
            primary_width = np.maximum(self.reaches[ireach].max_width,
                                       self.reaches[ireach].width)
            # fix any potential NaN or zero values from PRD
            primary_width[np.isnan(primary_width)] = np.nanmedian(primary_width)
            primary_width[primary_width <= 0] = np.nanmedian(primary_width)
            # compute extreme distance threshold
            max_distance = (
                primary_width * self.reaches[ireach].wth_coef)
            extreme_dist = self.reaches[ireach].ext_dist_coef * np.array(
                [abs(river_obs.ds), max_distance]).max(axis=0)

            ds = abs(river_obs.s - river_obs.centerline.s[river_obs.index])
            dn = abs(river_obs.n)

            mask_keep = np.logical_and(ds <= extreme_dist[river_obs.index],
                                       dn <= extreme_dist[river_obs.index])

            river_obs.mask_pixels(mask_keep)

            if len(river_obs.populated_nodes) > 0:
                river_obs_list_out.append(river_obs)
                reach_idx_list_out.append(reach_idx)
                ireach_list_out.append(ireach)

        return river_obs_list_out, reach_idx_list_out, ireach_list_out

    def process_node(self,
                     reach,
                     reach_id,
                     reach_idx=None,
                     scalar_max_width=600.,
                     max_width=None,
                     ds=None,
                     refine_centerline=False,
                     smooth=1.e-2,
                     alpha=1.,
                     max_iter=1):
        """
        Estimate the node quantities for a single reach
        Parameters
        ----------
        reach : Reach instance
            One of the reaches from ReachExtractor.
        reach_id : int
            Index in the list of reaches extracted for this scene.
        reach_idx, int
            Reach index used as pointer to a reach collection that may be
            different than the one used as input (e.g., a global reach
            collection). If using a width database, the retrieved width will
            be associated with this reach index. If None, it will be set equal
            to the reach_id.
        scalar_max_width : float, default 600
            How far away to look for points
        max_width: float or array_like, optional
            Maximum width to use for accepting points. From width database or
            apriori.
        ds : float, optional
            Separation between centerline nodes (in m). If None, uses reach
            points.
        refine_centerline: bool, default True
            Refine the centerline?
        smooth : float, default 1.e-2
            Centerline smoothing constant (see Centerline)
        alpha : float, default 1
            Centerline refinement update weight
        max_iter : int, default 1
            Maximum number of centerline iterations
        Returns
        -------
        Return a RiverReach instance with the node quantities populated but
        not the reach quantities.
        """
        # Refine the centerline, if desired
        # get the number of node in the reach and only refine if there are
        # enough to do spline
        num_nodes = len(np.unique(self.river_obs.index))
        enough_nodes = True if num_nodes - 1 > self.river_obs.k else False
        LOGGER.debug(
            "process_node for reach_id: {}; {} nodes".format(
                reach_idx, num_nodes))

        if refine_centerline and enough_nodes:
            self.river_obs.iterate(
                max_iter=max_iter, alpha=alpha, weights=True, smooth=smooth)

            # Associate the width to the new centerline
            if np.iterable(max_width):
                xw = reach.x
                yw = reach.y
                self.river_obs.add_centerline_obs(xw, yw, max_width,
                                                  'max_width')

            # Reinitialize to the new centerline and max_width
            self.river_obs.reinitialize()
            LOGGER.debug('centerline refined')

        else:
            # Associate the width to the new centerline
            if np.iterable(max_width):
                xw = reach.x
                yw = reach.y
                self.river_obs.add_centerline_obs(xw, yw, max_width,
                                                  'max_width')

        # Exclude beginning and end nodes, if desired
        if self.trim_ends:
            first_node = self.river_obs.populated_nodes[0]
            last_node = self.river_obs.populated_nodes[-1]
            self.river_obs.remove_nodes([first_node, last_node])

        # write out the image coordinates for each node in a netcdf file
        try:
            segOut = self.seg_label[self.river_obs.in_channel]

        except TypeError:
            segOut = None

        # get the prior locations and indices of the nodes
        xw = reach.x
        yw = reach.y
        self.river_obs.add_centerline_obs(xw, yw, xw, 'x_prior')
        self.river_obs.add_centerline_obs(xw, yw, yw, 'y_prior')

        # should probably compute prior lat/lon from prior x,y too
        windex = self.river_obs.centerline_obs['x_prior'].populated_nodes
        x_prior = xw
        x_prior[windex] = self.river_obs.centerline_obs['x_prior'].v
        x_prior = x_prior[self.river_obs.populated_nodes]

        windex = self.river_obs.centerline_obs['y_prior'].populated_nodes
        y_prior = yw
        y_prior[windex] = self.river_obs.centerline_obs['y_prior'].v

        if reach.node_indx is None:
            node_indx = np.arange(len(y_prior))
        else:
            node_indx = reach.node_indx

        node_indx = node_indx[self.river_obs.populated_nodes]
        y_prior = y_prior[self.river_obs.populated_nodes]
        reach_index = np.ones(len(node_indx)) * (reach_idx)

        if self.output_file is not None:
            self.write_index_file(
                self.img_x[self.river_obs.in_channel],
                self.img_y[self.river_obs.in_channel],
                reach.node_indx[self.river_obs.index],
                self.river_obs.d, self.river_obs.s,
                self.river_obs.n, reach_idx,
                segOut, self.lat[self.river_obs.in_channel],
                self.lon[self.river_obs.in_channel],
                self.h_noise[self.river_obs.in_channel],
                self.h_flg[self.river_obs.in_channel],
                self.area_flg[self.river_obs.in_channel],
                self.sig0_flg[self.river_obs.in_channel],
                self.wse_class_flg[self.river_obs.in_channel],
                self.geolocation_qual[self.river_obs.in_channel],
                self.classification_qual[self.river_obs.in_channel],
                self.sig0_qual[self.river_obs.in_channel])

        # Load these datasets from input file and add observations to
        # self.river_obs.
        dsets = [
            'h_noise', 'h_flg', 'wse_class_flg', 'area_flg', 'sig0_flg', 'lon',
            'lat', 'inundated_area', 'klass', 'pixel_area', 'sring_flg',
            'xtrack', 'sig0', 'sig0_uncert', 'water_frac', 'water_frac_uncert',
            'ifgram', 'power1', 'power2', 'phase_noise_std', 'dh_dphi',
            'dlat_dphi', 'dlon_dphi', 'num_rare_looks', 'num_med_looks',
            'false_detection_rate', 'missed_detection_rate', 'darea_dheight',
            'looks_to_efflooks', 'geoid', 'solid_earth_tide', 'load_tide_fes',
            'load_tide_got', 'pole_tide', 'is_area_degraded',
            'is_area_suspect', 'is_wse_degraded', 'is_wse_suspect',
            'is_sig0_bad', 'is_sig0_suspect', 'bright_land_flag',
            'is_xovercal_suspect', 'is_xovercal_degraded', 'layover_impact']

        dsets_to_load = []
        for name in dsets:
            value = getattr(self, name)
            if value is not None:
                dsets_to_load.append(name)
                if name == 'looks_to_efflooks':
                    value = value + np.zeros(np.shape(self.lat))
                self.river_obs.add_obs(name, value)

        # need to get array of land/water edge classes
        # to decode/encode the classification routine 
        # for external call to area agg in RiverNode
        edge_water = np.zeros(np.shape(self.klass))
        for i, k in enumerate(self.class_list):
            if self.use_fractional_inundation[i]:
                # this is actually both land and water edges,
                # but setting to water edge
                edge_water[self.klass == k] = 1

        self.river_obs.add_obs('edge_water', edge_water)
        dsets_to_load.append('edge_water')

        # Adjust heights to geoid and do tide corrections
        # (need to do before load_nodes or it is not updated in nodes)
        try:
            mask = np.logical_and(
                self.river_obs.geoid > -200, self.river_obs.geoid < 200)
            self.river_obs.h_noise[mask] -= (
                self.river_obs.geoid[mask] +
                self.river_obs.solid_earth_tide[mask] +
                self.river_obs.load_tide_fes[mask] +
                self.river_obs.pole_tide[mask])
        except AttributeError:
            LOGGER.warning(
                'Unable to adjust heights to geoid or do tide corrections')

        # add these separately since input/output names don't match
        self.river_obs.add_obs('xobs', self.x)
        self.river_obs.add_obs('yobs', self.y)
        dsets_to_load.append('xobs')
        dsets_to_load.append('yobs')

        self.river_obs.load_nodes(dsets_to_load)

        # Get various node statistics
        nobs = np.asarray(self.river_obs.get_node_stat('count', ''))
        s_median = np.asarray(self.river_obs.get_node_stat('median', 's'))
        x_median = np.asarray(self.river_obs.get_node_stat('median', 'xobs'))
        y_median = np.asarray(self.river_obs.get_node_stat('median', 'yobs'))

        if self.xtrack is not None:
            xtrack_median = np.asarray(
                self.river_obs.get_node_stat('median', 'xtrack'))
        else:
            xtrack_median = np.nan*np.ones(nobs.shape)

        # Compute number of good+sus pixels in wse / area masks
        n_pix_wse_good_sus = np.array(
            self.river_obs.get_node_stat('count_good', 'h_flg'))
        n_pix_area_good_sus = np.array(
            self.river_obs.get_node_stat('count_good', 'area_flg'))

        # Use good/suspect aggregations if number of pixels used is greater or
        # equal to num_good_sus_pix_thresh_(wse/area)
        mask_good_sus_wse = (
            n_pix_wse_good_sus >= self.num_good_sus_pix_thresh_wse)
        mask_good_sus_area = (
            n_pix_area_good_sus >= self.num_good_sus_pix_thresh_area)

        # prints warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lon_median = np.asarray(self.river_obs.get_node_stat(
                'sincos_median', 'lon', goodvar='h_flg'))
            lon_median[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat(
                    'sincos_median', 'lon', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]

            lat_median = np.asarray(self.river_obs.get_node_stat(
                'median', 'lat', goodvar='h_flg'))
            lat_median[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat(
                    'median', 'lat', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]

        ds = np.asarray(self.river_obs.get_node_stat('value', 'ds'))

        # number of good heights
        nobs_h = np.asarray(self.river_obs.get_node_stat('count_good', 'h_flg'))
        nobs_h[~mask_good_sus_wse] = np.asarray(
            self.river_obs.get_node_stat('count_good', 'wse_class_flg')
            )[~mask_good_sus_wse]

        # heights using only "good" heights
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wse = np.asarray(self.river_obs.get_node_stat(
                    'median', 'h_noise', goodvar='h_flg'))
            wse[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat(
                    'median', 'h_noise', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]

            wse_r_u = np.asarray(self.river_obs.get_node_stat(
                'std', 'h_noise', goodvar='h_flg'))
            wse_r_u[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat(
                    'std', 'h_noise', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]

        wse_std = wse_r_u

        # These are area based estimates, the nominal SWOT approach
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            area = np.asarray(
                self.river_obs.get_node_stat('sum', 'inundated_area'))
            area[~mask_good_sus_area] = np.asarray(
                self.river_obs.get_node_stat(
                    'sum', 'inundated_area', goodvar='good')
                )[~mask_good_sus_area]
            # compute the area of all pixels flagged specular ringing
            area_sring = np.asarray(
                self.river_obs.get_node_stat('sum', 'inundated_area',
                goodvar='sring_flg'))

        # area of pixels used to compute heights
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            area_of_ht = np.asarray(self.river_obs.get_node_stat(
                'sum', 'inundated_area', goodvar='h_flg'))
            area_of_ht[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat(
                    'sum', 'inundated_area', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]

        width_area = np.asarray(self.river_obs.get_node_stat(
            'width_area', 'inundated_area'))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                rdr_sig0 = np.asarray(self.river_obs.get_node_stat(
                    'median', 'sig0', goodvar='h_flg'))
            except AttributeError:
                rdr_sig0 = np.nan*np.ones(area.shape)

        # get the aggregated heights and widths with their corresponding
        # uncertainty estimates all in one shot
        if self.height_agg_method == 'orig' and self.area_agg_method == 'orig':
            # copy some of these over for the output products
            area_det = area
            area_u = np.nan*np.ones(area.shape)
            width_u = np.nan*np.ones(area.shape)
            area_det_u = np.nan*np.ones(area.shape)
            rdr_sig0_u = np.nan*np.ones(area.shape)
            latitude_u = np.nan*np.ones(area.shape)
            longitud_u = np.nan*np.ones(area.shape)
        else:
            node_aggs = self.river_obs.get_node_agg(
                height_method=self.height_agg_method,
                area_method=self.area_agg_method, goodvar_wse='h_flg',
                goodvar_area='area_flg', goodvar_sig0='sig0_flg')

            node_aggs_w_degraded = self.river_obs.get_node_agg(
                height_method=self.height_agg_method,
                area_method=self.area_agg_method, goodvar_wse='wse_class_flg',
                goodvar_area='good', goodvar_sig0='sig0_flg')

            latitude_u = node_aggs['lat_u']
            latitude_u[~mask_good_sus_wse] = node_aggs_w_degraded[
                'lat_u'][~mask_good_sus_wse]

            longitud_u = node_aggs['lon_u']
            latitude_u[~mask_good_sus_wse] = node_aggs_w_degraded[
                'lon_u'][~mask_good_sus_wse]

            rdr_sig0 = node_aggs['sig0']
            rdr_sig0_u = node_aggs['sig0_u']

            if self.area_agg_method != 'orig':

                width_area = node_aggs['width_area']
                width_u = node_aggs['width_area_u']
                area = node_aggs['area']
                area_u = node_aggs['area_u']
                area_det = node_aggs['area_det']
                area_det_u = node_aggs['area_det_u']
                area_of_ht = node_aggs['area_of_ht']

                # Use degraded pix as well when not enough good/suspect pix
                width_area[~mask_good_sus_area] = node_aggs_w_degraded[
                    'width_area'][~mask_good_sus_area]
                width_u[~mask_good_sus_area] = node_aggs_w_degraded[
                    'width_area_u'][~mask_good_sus_area]
                area[~mask_good_sus_area] = node_aggs_w_degraded[
                    'area'][~mask_good_sus_area]
                area_u[~mask_good_sus_area] = node_aggs_w_degraded[
                    'area_u'][~mask_good_sus_area]
                area_det[~mask_good_sus_area] = node_aggs_w_degraded[
                    'area_det'][~mask_good_sus_area]
                area_det_u[~mask_good_sus_area] = node_aggs_w_degraded[
                    'area_det_u'][~mask_good_sus_area]
                area_of_ht[~mask_good_sus_area] = node_aggs_w_degraded[
                    'area_of_ht'][~mask_good_sus_area]

            if self.height_agg_method != 'orig':
                wse = node_aggs['h']
                wse[~mask_good_sus_wse] = node_aggs_w_degraded[
                    'h'][~mask_good_sus_wse]

                # height_uncert_std
                wse_std = node_aggs['h_std']
                wse_std[~mask_good_sus_wse] = node_aggs_w_degraded[
                    'h_std'][~mask_good_sus_wse]

                # height_uncert_multilook
                wse_r_u = node_aggs['h_u']
                wse_r_u[~mask_good_sus_wse] = node_aggs_w_degraded[
                    'h_u'][~mask_good_sus_wse]

        # geoid heights and tide corrections weighted by height uncertainty
        try:
            geoid_hght = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean', 'geoid',
                                             goodvar='h_flg',
                                             method=self.height_agg_method))
            geoid_hght[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean', 'geoid',
                                             goodvar='wse_class_flg',
                                             method=self.height_agg_method)
                                             )[~mask_good_sus_wse]
        except AttributeError:
            geoid_hght = np.asarray(self.river_obs.get_node_stat(
                'mean', 'geoid', goodvar = 'h_flg'))

        try:
            solid_tide = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'solid_earth_tide',
                                             goodvar='h_flg',
                                             method=self.height_agg_method))
            solid_tide[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'solid_earth_tide',
                                             goodvar='wse_class_flg',
                                             method=self.height_agg_method)
                                             )[~mask_good_sus_wse]
        except AttributeError:
            solid_tide = np.asarray(self.river_obs.get_node_stat(
                'mean', 'solid_earth_tide', goodvar='h_flg'))

        try:
            load_tidef = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'load_tide_fes',
                                             goodvar='h_flg',
                                             method=self.height_agg_method))
            load_tidef[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'load_tide_fes',
                                             goodvar='wse_class_flg',
                                             method=self.height_agg_method)
                                             )[~mask_good_sus_wse]
        except AttributeError:
            load_tidef = np.asarray(self.river_obs.get_node_stat(
                'mean', 'load_tide_fes', goodvar='h_flg'))

        try:
            load_tideg = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'load_tide_got',
                                             goodvar='h_flg',
                                             method=self.height_agg_method))
            load_tideg[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'load_tide_got',
                                             goodvar='wse_class_flg',
                                             method=self.height_agg_method)
                                             )[~mask_good_sus_wse]
        except AttributeError:
            load_tideg = np.asarray(self.river_obs.get_node_stat(
                'mean', 'load_tide_got', goodvar='h_flg'))

        try:
            pole_tide = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'pole_tide',
                                             goodvar='h_flg',
                                             method=self.height_agg_method))
            pole_tide[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'pole_tide',
                                             goodvar='wse_class_flg',
                                             method=self.height_agg_method)
                                             )[~mask_good_sus_wse]
        except AttributeError:
            pole_tide = np.asarray(self.river_obs.get_node_stat(
                'mean', 'pole_tide', goodvar='h_flg'))

        try:
            layovr_val = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'layover_impact',
                                             goodvar='h_flg',
                                             method=self.height_agg_method))
            layovr_val[~mask_good_sus_wse] = np.asarray(
                self.river_obs.get_node_stat('height_weighted_mean',
                                             'layover_impact',
                                             goodvar='wse_class_flg',
                                             method=self.height_agg_method)
                                             )[~mask_good_sus_wse]
        except AttributeError:
            layovr_val = np.nan*np.ones(area.shape)

        # These are the values from the width database
        width_db = np.ones(self.river_obs.n_nodes, dtype=np.float64) * \
            self.river_obs.missing_value

        try:
            windex = self.river_obs.centerline_obs['max_width'].populated_nodes
            width_db[windex] = self.river_obs.centerline_obs['max_width'].v
            width_db = width_db[self.river_obs.populated_nodes]

        except KeyError:
            width_db = np.ones(
                self.river_obs.n_nodes,
                dtype=np.float64) * self.river_obs.missing_value
            width_db = width_db[self.river_obs.populated_nodes]

        # flag for bad widths
        min_n = np.array(self.river_obs.get_node_stat('min', 'n'))
        max_n = np.array(self.river_obs.get_node_stat('max', 'n'))
        blocking_width = reach.blocking_widths[self.river_obs.populated_nodes]

        # test at 10% inside of blocking width
        test_width = blocking_width * BLOCKING_WIDTH_THRESHOLD_FRACTION
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            is_blocked = np.logical_or(
                np.logical_and(test_width < 0, min_n < test_width),
                np.logical_and(test_width > 0, max_n > test_width))

        lon_prior = reach.lon[self.river_obs.populated_nodes]
        lat_prior = reach.lat[self.river_obs.populated_nodes]

        dark_frac = MISSING_VALUE_FLT * np.ones(area.shape)
        dark_frac[area > 0] = 1 - area_det[area > 0] / area[area > 0]
        dark_frac[np.logical_and(dark_frac > 1, area > 0)] = 1  # clip to <= 1
        sring_frac = MISSING_VALUE_FLT * np.ones(area.shape)
        sring_frac[area > 0] = area_sring[area > 0] / area[area > 0]
        sring_frac[np.logical_and(sring_frac > 1, area > 0)] = 1

        # Compute flow direction relative to along-track
        tangent = self.river_obs.centerline.tangent[
            self.river_obs.populated_nodes]

        fit_xx = np.array([np.ones(x_median.shape), x_median, y_median]).T

        if xtrack_median is not None:
            fit = statsmodels.api.OLS(xtrack_median, fit_xx).fit()
            dxt_dx = fit.params[1]
            dxt_dy = fit.params[2]
            xt_angle = np.arctan2(dxt_dy, dxt_dx)
            at_angle = xt_angle-np.pi/2
            tangent_angle = np.arctan2(tangent[:, 1], tangent[:, 0])
            flow_dir = np.rad2deg(tangent_angle - at_angle) % 360
        else:
            flow_dir = np.nan*np.ones(area.shape)

        prior_s = np.cumsum(reach.node_length)

        # total number of pixels assigned to each node
        n_pix = np.array(self.river_obs.get_node_stat('count', 'x'))

        # Number of pixels that are in wse class types for each node
        n_pix_wse_class = np.array(
            self.river_obs.get_node_stat('count_good', 'wse_class_flg'))

        # Number of pixels used for wse/area/sig0 for each node
        n_pix_wse = n_pix_wse_good_sus
        n_pix_wse[~mask_good_sus_wse] = n_pix_wse_class[~mask_good_sus_wse]

        n_pix_area = n_pix_area_good_sus
        n_pix_area[~mask_good_sus_area] = n_pix[~mask_good_sus_area]

        n_pix_sig0 = np.array(
            self.river_obs.get_node_stat('count_good', 'sig0_flg'))

        # Number of pixels with suspect wse/area/sig0 used in aggregations for
        # each node
        n_pix_wse_suspect = np.array(self.river_obs.get_node_stat(
            'sum', 'is_wse_suspect', goodvar='h_flg'))
        n_pix_area_suspect = np.array(self.river_obs.get_node_stat(
            'sum', 'is_area_suspect', goodvar='area_flg'))
        n_pix_sig0_suspect = np.array(self.river_obs.get_node_stat(
            'sum', 'is_sig0_suspect', goodvar='sig0_flg'))

        # Number of pixels with degraded wse/area that were rejected for each
        # node
        n_pix_wse_degraded = np.array(self.river_obs.get_node_stat(
            'sum', 'is_wse_degraded', goodvar='wse_class_flg'))
        n_pix_area_degraded = np.array(self.river_obs.get_node_stat(
            'sum', 'is_area_degraded'))

        # create empty arrays for the reconst WSE & uncertainty (for later)
        w_opt = np.ones(wse.shape,
                        dtype=np.float64) * self.river_obs.missing_value
        w_opt_r_u = np.ones(wse.shape,
                        dtype=np.float64) * self.river_obs.missing_value

        # Create node_q_b quality bitwise flag
        node_q_b = np.zeros(lat_median.shape, dtype='i4')

        # bit 0 / sig0_qual_suspect
        node_q_b[n_pix_sig0_suspect > 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_SIG0_QUAL_SUSPECT)

        # bit 1 / classification_qual_suspect
        node_q_b[n_pix_area_suspect > 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_CLASS_QUAL_SUSPECT)

        # bit 2 / geolocation_qual_suspect
        node_q_b[n_pix_wse_suspect > 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_GEOLOCATION_QUAL_SUSPECT)

        # bit 3 / water_fraction_suspect
        # If water_frac not an attribute of nodes, skip (GPS profiles)
        with contextlib.suppress(AttributeError):
            water_frac_max = np.array(self.river_obs.get_node_stat(
                'max', 'water_frac', goodvar='area_flg'))
            water_frac_max[~mask_good_sus_area] = np.array(
                self.river_obs.get_node_stat('max', 'water_frac',
                goodvar='good'))[~mask_good_sus_area]
            node_q_b[water_frac_max > WATER_FRAC_SUSPECT_THRESHOLD] |= (
                SWOTRiver.products.rivertile.QUAL_IND_WATER_FRAC_SUSPECT);

        # bit 4 / blocking_width_suspect
        node_q_b[is_blocked] |= (
            SWOTRiver.products.rivertile.QUAL_IND_BLOCK_WIDTH_SUSPECT)

        # bit 7 / bright_land
        # If bright_land_flag not an attribute of nodes, skip (GPS profiles)
        with contextlib.suppress(AttributeError):
            n_pix_bright_land = np.array(
                self.river_obs.get_node_stat('sum', 'bright_land_flag'))
            node_q_b[n_pix_bright_land > 0] |= (
                SWOTRiver.products.rivertile.QUAL_IND_BRIGHT_LAND_SUSPECT)

        # bit (9/few_sig0_observations) / (10/few_area_observations) /
        # (11/few_wse_observations)
        node_q_b[n_pix_sig0 < TOO_FEW_PIXELS_THRESHOLD] |= (
            SWOTRiver.products.rivertile.QUAL_IND_FEW_SIG0_PIX)
        node_q_b[n_pix_area < TOO_FEW_PIXELS_THRESHOLD] |= (
            SWOTRiver.products.rivertile.QUAL_IND_FEW_AREA_PIX)
        node_q_b[n_pix_wse < TOO_FEW_PIXELS_THRESHOLD] |= (
            SWOTRiver.products.rivertile.QUAL_IND_FEW_WSE_PIX)

        # bit 13 / far_range_suspect and bit 14 / near_range_suspect
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            node_q_b[np.abs(xtrack_median) > FAR_RANGE_XTRACK_TRESHOLD] |= (
                SWOTRiver.products.rivertile.QUAL_IND_FAR_RANGE_SUSPECT)
            # bit 14 / near_range_suspect
            node_q_b[np.abs(xtrack_median) < NEAR_RANGE_XTRACK_THRESHOLD] |= (
                SWOTRiver.products.rivertile.QUAL_IND_NEAR_RANGE_SUSPECT)

        # bit 18 / classification_qual_degraded
        # Only set if not enough pixels for good/sus only aggregations and
        # degraded data is actually used (n_pix_area_degraded>0).
        node_q_b[np.logical_and(
            ~mask_good_sus_area, n_pix_area_degraded > 0)] |= (
                SWOTRiver.products.rivertile.QUAL_IND_CLASS_QUAL_DEGRADED)

        # bit 19 / geolocation_qual_degraded
        # Only set if not enough pixels for good/sus only aggregations, and
        # n_pix_wse_degraded is not zero
        node_q_b[np.logical_and(
            ~mask_good_sus_wse, n_pix_wse_degraded > 0)] |= (
                SWOTRiver.products.rivertile.QUAL_IND_GEOLOCATION_QUAL_DEGRADED)

        # bit 23 / wse_outlier will be set in later stage of processing when
        # reach masks are generated)

        # bit 24 / wse_bad
        node_q_b[np.logical_or(wse < MIN_VALID_WSE, wse > MAX_VALID_WSE)] |= (
            SWOTRiver.products.rivertile.QUAL_IND_WSE_BAD)

        # bit 25 / no_sig0_observations
        node_q_b[n_pix_sig0 == 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_NO_SIG0_PIX)

        # bit 26 / no_area_observations
        node_q_b[n_pix_area == 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_NO_SIG0_PIX)

        # bit 27 / no_wse_observations
        node_q_b[n_pix_wse == 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_NO_SIG0_PIX)

        # bit 28 / no_pixels
        node_q_b[n_pix == 0] |= (
            SWOTRiver.products.rivertile.QUAL_IND_NO_SIG0_PIX)

        # Create node_q from node_q_b
        thresh_sus = 1
        thresh_deg = (
            SWOTRiver.products.rivertile.QUAL_IND_NODE_DEGRADED_THRESHOLD)
        thresh_bad = SWOTRiver.products.rivertile.QUAL_IND_NODE_BAD_THRESHOLD

        node_q = np.zeros(lat_median.shape, dtype='i2')
        node_q[node_q_b >= thresh_sus] = 1
        node_q[node_q_b >= thresh_deg] = 2
        node_q[node_q_b >= thresh_bad] = 3

        # create xovr_cal_q
        xovr_cal_q = np.zeros(lat_median.shape, dtype='i2')
        n_pix_xovercal_suspect = np.array(
            self.river_obs.get_node_stat(
                'sum', 'is_xovercal_suspect', goodvar='h_flg'))
        n_pix_xovercal_suspect[~mask_good_sus_wse] = np.array(
            self.river_obs.get_node_stat(
                'sum', 'is_xovercal_suspect', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]
        xovr_cal_q[n_pix_xovercal_suspect > 0] = 1

        n_pix_xovercal_degraded = np.array(
            self.river_obs.get_node_stat(
                'sum', 'is_xovercal_degraded', goodvar='h_flg'))
        n_pix_xovercal_degraded[~mask_good_sus_wse] = np.array(
            self.river_obs.get_node_stat(
                'sum', 'is_xovercal_degraded', goodvar='wse_class_flg')
                )[~mask_good_sus_wse]
        xovr_cal_q[n_pix_xovercal_degraded > 0] = 2

        # if no valid pixels for height, set xovr_cal_q to BAD
        xovr_cal_q[n_pix_wse == 0] = 2

        # initiate river reach mask, fill with Trues
        river_node_mask = np.full(len(wse), True, dtype=bool)

        # type cast node outputs and pack it up for RiverReach constructor
        river_reach_kw_args = {
            'lat': lat_median.astype('float64'),
            'lon': lon_median.astype('float64'),
            'x': x_median.astype('float64'),
            'y': y_median.astype('float64'),
            's': s_median.astype('float64'),
            'ds': ds.astype('float64'),
            'w_area': width_area.astype('float64'),
            'w_db': width_db.astype('float64'),
            'area': area.astype('float64'),
            'area_u': area_u.astype('float64'),
            'area_det': area_det.astype('float64'),
            'area_det_u': area_det_u.astype('float64'),
            'area_of_ht': area_of_ht.astype('float64'),
            'area_sring': area_sring.astype('float64'),
            'wse': wse.astype('float64'),
            'wse_std': wse_std.astype('float64'),
            'wse_r_u': wse_r_u.astype('float64'),
            'w_opt': w_opt.astype('float64'),
            'w_opt_r_u': w_opt_r_u.astype('float64'),
            'nobs': nobs.astype('int32'),
            'nobs_h': nobs_h.astype('int32'),
            'n_good_pix': n_pix_wse.astype('int32'),
            'node_indx': node_indx.astype('int64'),
            'reach_indx': reach_index.astype('int64'),
            'rdr_sig0': rdr_sig0.astype('float64'),
            'rdr_sig0_u': rdr_sig0_u.astype('float64'),
            'latitude_u': latitude_u.astype('float64'),
            'longitud_u': longitud_u.astype('float64'),
            'width_u': width_u.astype('float64'),
            'geoid_hght': geoid_hght.astype('float64'),
            'solid_tide': solid_tide.astype('float64'),
            'load_tidef': load_tidef.astype('float64'),
            'load_tideg': load_tideg.astype('float64'),
            'pole_tide': pole_tide.astype('float64'),
            'node_blocked': is_blocked.astype('uint8'),
            'dark_frac': dark_frac,
            'sring_frac': sring_frac,
            'x_prior': x_prior.astype('float64'),
            'y_prior': y_prior.astype('float64'),
            'lon_prior': lon_prior.astype('float64'),
            'lat_prior': lat_prior.astype('float64'),
            'p_wse': reach.wse[self.river_obs.populated_nodes],
            'p_wse_var': reach.wse_var[self.river_obs.populated_nodes],
            'p_wse_all': reach.wse,
            'p_width': reach.width[self.river_obs.populated_nodes],
            'p_wid_var': reach.width_var[self.river_obs.populated_nodes],
            'p_dist_out': reach.dist_out[self.river_obs.populated_nodes],
            'p_length': reach.node_length[self.river_obs.populated_nodes],
            'grand_id': reach.grod_id[self.river_obs.populated_nodes],
            'n_chan_max': reach.n_chan_max[self.river_obs.populated_nodes],
            'n_chan_mod': reach.n_chan_mod[self.river_obs.populated_nodes],
            'flow_dir': flow_dir.astype('float64'),
            'prior_node_ss': prior_s,
            'node_ss': prior_s[self.river_obs.populated_nodes],
            'populated_nodes': self.river_obs.populated_nodes,
            'ice_clim_f': reach.metadata['iceflag']*np.ones(lat_median.shape),
            'river_name': reach.river_name[self.river_obs.populated_nodes],
            'node_q': node_q, 'node_q_b': node_q_b, 'xovr_cal_q': xovr_cal_q,
            'layovr_val': layovr_val.astype('float64'),
            'mask': river_node_mask
        }

        # Get wse_u from RSS of random wse_r_u and REACH_WSE_SYS_UNCERT
        river_reach_kw_args['wse_u'] = MISSING_VALUE_FLT * np.ones(
            river_reach_kw_args['wse_r_u'].shape)
        mask = river_reach_kw_args['wse_r_u'] > 0
        river_reach_kw_args['wse_u'][mask] = np.sqrt(
            river_reach_kw_args['wse_r_u'][mask]**2 + REACH_WSE_SYS_UNCERT**2
            ).astype('float64')

        if xtrack_median is not None:
            river_reach_kw_args['xtrack'] = xtrack_median.astype('float64')
        else:
            river_reach_kw_args['xtrack'] = None

        # For swotCNES
        river_reach_kw_args['h_n_ave'] = river_reach_kw_args['wse']

        river_reach = RiverReach(**river_reach_kw_args)

        return river_reach

    def process_reach(self, river_reach_collection, river_reach, reach,
                      reach_id, reach_idx=None, min_fit_points=2):
        """
        Estimate the width, height, and slope for one reach.

        Parameters
        ----------
        river_reach : partially populated RiverReach instance with node
            quantities already computed
        reach : Reach instance
            One of the reaches from ReachExtractor.
        reach_id : int
            Index in the list of reaches extracted for this scene.
        reach_idx, int
            Reach index used as pointer to a reach collection that may be
            different than the one used as input (e.g., a global reach
            collection). If using a width database, the retrieved width will
            be associated with this reach index. If None, it will be set equal
            to the reach_id.
        min_fit_points : int, default 2
            Minimum number of populated nodes required for height/slope fit

        Returns
        -------
        Nothing

        Modifies river_reach with the reach quantities (river_reach.metadata)
        """
        # Check to see if there are sufficient number of points for fit
        ngood = len(river_reach.s)
        LOGGER.debug(("process_reach for reach_id: {}; number of nodes: {}"
                      ).format(reach_idx, ngood))

        reach_stats = collections.OrderedDict()
        reach_stats['length'] = np.sum(river_reach.p_length)
        reach_stats['reach_id'] = reach_id
        reach_stats['reach_idx'] = reach_idx

        reach_stats['node_dist'] = np.mean(np.sqrt(
                (river_reach.x-river_reach.x_prior)**2 +
                (river_reach.y-river_reach.y_prior)**2))

        reach_stats['area'] = np.sum(river_reach.area)
        reach_stats['area_u'] = np.sqrt(np.sum(
            river_reach.area_u**2))

        reach_stats['area_det'] = np.sum(river_reach.area_det)
        reach_stats['area_det_u'] = np.sqrt(np.sum(
            river_reach.area_det_u**2))

        reach_stats['area_of_ht'] = np.sum(river_reach.area_of_ht)

        pct_unmasked = 100*np.sum(river_reach.mask) / len(river_reach.mask)
        if pct_unmasked >= self.reach_pct_good_sus_thresh:
            reach_stats['width'] = np.sum(river_reach.area[river_reach.mask])\
                                   /np.sum(river_reach.p_length[river_reach.mask])
            reach_stats['width_u'] = np.sqrt(
                np.sum(river_reach.area_u[river_reach.mask]**2)) / np.sum(
                river_reach.p_length[river_reach.mask])
        else:
            reach_stats['width'] = reach_stats['area'] / reach_stats['length']
            reach_stats['width_u'] = np.sqrt(
                np.sum(river_reach.area_u**2))/ reach_stats['length']
        reach_stats['layovr_val'] = np.sqrt(np.sum(
            river_reach.layovr_val[river_reach.mask]**2))

        reach_stats['loc_offset'] = (
            river_reach.s.mean() - self.river_obs.centerline.s.mean())

        if river_reach.xtrack is not None:
            reach_stats['xtrk_dist'] = np.ma.median(river_reach.xtrack)

        # Along-flow distance for all PRD nodes.
        all_ss = np.cumsum(reach.node_length)

        # Make center of PRD reach the ss == 0 intercept (reference point) for
        # reach WSE.
        all_ss = all_ss - np.mean(all_ss)

        # Along-flow distance for just the observed nodes.
        ss = all_ss[self.river_obs.populated_nodes]

        hh = river_reach.wse
        wse_r_u = river_reach.wse_r_u
        ww = 1/(wse_r_u**2)
        SS = np.c_[ss, np.ones(len(ss), dtype=ss.dtype)]
        mask = river_reach.mask

        if mask.sum() >= min_fit_points:
            # compute slope according to input config method
            if self.slope_method == 'first_to_last':
                reach_stats['slope'] = (
                    hh[mask][0]-hh[mask][-1])/(ss[mask][0]-ss[mask][-1])
                reach_stats['height'] = (
                    np.mean(hh[mask]) + reach_stats['slope'] *
                    (np.mean(all_ss)-np.mean(ss[mask])))

                # TBD on unc quantities for first_to_last method
                reach_stats['slope_r_u'] = MISSING_VALUE_FLT
                reach_stats['height_r_u'] = MISSING_VALUE_FLT
                reach_stats['slope_u'] = SLOPE_SYS_UNCERT
                reach_stats['height_u'] = REACH_WSE_SYS_UNCERT
                reach_stats['slope2_u'] = MISSING_VALUE_FLT
                reach_stats['slope2_r_u'] = MISSING_VALUE_FLT

            elif self.slope_method in ['unweighted', 'weighted']:
                # use weighted fit if commanded and all weights are good
                if self.slope_method == 'weighted' \
                        and all(np.isfinite(ww[mask])):
                    fit = statsmodels.api.WLS(
                        hh[mask], SS[mask], weights=ww[mask]).fit()

                # use unweighted fit
                else:
                    fit = statsmodels.api.OLS(hh[mask], SS[mask]).fit()

                # fit slope is meters per meter
                reach_stats['slope'] = fit.params[0]
                reach_stats['height'] = fit.params[1]

                # use White’s (1980) heteroskedasticity robust standard errors.
                # https://www.statsmodels.org/dev/generated/
                #    statsmodels.regression.linear_model.RegressionResults.html
                reach_stats['slope_r_u'] = fit.HC0_se[0]
                reach_stats['height_r_u'] = fit.HC0_se[1]
                reach_stats['slope_u'] = np.sqrt(
                    SLOPE_SYS_UNCERT**2 + reach_stats['slope_r_u']**2)
                reach_stats['height_u'] = np.sqrt(
                    REACH_WSE_SYS_UNCERT**2 + reach_stats['height_r_u']**2)

            elif self.slope_method == 'bayes':
                if self.bayes_slope_use_all_nodes:
                    # add unpopulated nodes.
                    hh_opt, wse_r_u_opt, _, mask_opt = \
                        self.add_unpopulated_nodes(river_reach)
                    ss_opt = all_ss
                else:
                    # populated nodes only
                    hh_opt, wse_r_u_opt, mask_opt, ss_opt = \
                        hh, wse_r_u, mask, ss
                # get the optimal reconstruction (Bayes estimate)
                wse_opt, wse_opt_r_u, height_u, slope_u,  = \
                    self.optimal_reconstruct(
                        river_reach_collection,
                        river_reach, reach_id,
                        ss_opt, hh_opt,
                        wse_r_u_opt, mask_opt,
                        min_fit_points,
                        method='Bayes',
                    )
                # Store the reconstructed node heights in river reach object.
                # Currently, only store populated nodes (node_indx clipped to
                # nodes with minobs observations as specified by L2_HR_Param
                # file).
                river_reach.w_opt = wse_opt[self.river_obs.populated_nodes]
                river_reach.w_opt_r_u = wse_opt_r_u[
                    self.river_obs.populated_nodes]
                # Use reconstruction height and slope for reach outputs
                dx = ss_opt[0] - ss_opt[-1]  # along-reach dist
                reach_stats['slope'] = (wse_opt[0] - wse_opt[-1]) / dx
                reach_stats['height'] = np.mean(wse_opt)
                reach_stats['slope_r_u'] = slope_u
                reach_stats['height_r_u'] = height_u
                reach_stats['slope_u'] = np.sqrt(
                    SLOPE_SYS_UNCERT**2 + reach_stats['slope_r_u']**2)
                reach_stats['height_u'] = np.sqrt(
                    REACH_WSE_SYS_UNCERT**2 + reach_stats['height_r_u']**2)
                reach_stats['slope2_u'] = MISSING_VALUE_FLT
                reach_stats['slope2_r_u'] = MISSING_VALUE_FLT

        else:
            # insufficient node heights for fit to reach
            reach_stats['slope'] = MISSING_VALUE_FLT
            reach_stats['slope_r_u'] = MISSING_VALUE_FLT
            reach_stats['slope_u'] = MISSING_VALUE_FLT
            reach_stats['height'] = MISSING_VALUE_FLT
            reach_stats['height_r_u'] = MISSING_VALUE_FLT
            reach_stats['height_u'] = MISSING_VALUE_FLT
            reach_stats['slope2_u'] = MISSING_VALUE_FLT
            reach_stats['slope2_r_u'] = MISSING_VALUE_FLT

        reach_stats['n_good_nod'] = mask.sum()
        reach_stats['frac_obs'] = (
            mask.sum() / len(self.river_obs.centerline.s))

        # do fit on geoid heights for reach-level outputs
        gg = river_reach.geoid_hght
        if mask.sum() >= min_fit_points:
            geoid_fit = statsmodels.api.WLS(
                gg[mask], SS[mask], weights=ww[mask]).fit()

            # fit slope is meters per meter
            reach_stats['geoid_slop'] = geoid_fit.params[0]
            reach_stats['geoid_hght'] = geoid_fit.params[1]
        else:
            reach_stats['geoid_slop'] = MISSING_VALUE_FLT
            reach_stats['geoid_hght'] = MISSING_VALUE_FLT

        # trap out of range / missing data
        if reach.metadata['lakeFlag'] < 0 or reach.metadata['lakeFlag'] > 255:
            uint8_flg = 255
        else:
            uint8_flg = reach.metadata['lakeFlag']
        reach_stats['lake_flag'] = uint8_flg
        reach_stats['centerline_lon'] = reach.metadata['centerline_lon']
        reach_stats['centerline_lat'] = reach.metadata['centerline_lat']

        # add fit_height for improved geolocation
        if reach_stats['slope'] != MISSING_VALUE_FLT:
            river_reach.fit_height = (
                reach_stats['height'] + reach_stats['slope'] * ss)
        else:
            river_reach.fit_height = MISSING_VALUE_FLT * np.ones(ss.shape)

        # copy things from the prior DB into reach outputs
        reach_stats['rch_id_up'] = np.array([
            item[0] for item in reach.metadata['rch_id_up']], dtype='i8')
        reach_stats['rch_id_up'][reach_stats['rch_id_up']==0] = \
            MISSING_VALUE_INT9

        reach_stats['rch_id_dn'] = np.array([
            item[0] for item in reach.metadata['rch_id_dn']], dtype='i8')
        reach_stats['rch_id_dn'][reach_stats['rch_id_dn']==0] = \
            MISSING_VALUE_INT9

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dark_frac = (
                1-np.sum(river_reach.area_det)/np.sum(river_reach.area))
            reach_stats['dark_frac'] = min(dark_frac, 1)  # clip to <= 1
            sring_frac = (
                np.sum(river_reach.area_sring)/np.sum(river_reach.area))
            reach_stats['sring_frac'] = min(sring_frac, 1)  # clip to <= 1

        reach_stats['n_reach_up'] = (reach_stats['rch_id_up'] > 0).sum()
        reach_stats['n_reach_dn'] = (reach_stats['rch_id_dn'] > 0).sum()

        reach_stats['p_lon'] = reach.metadata['lon']
        reach_stats['p_lat'] = reach.metadata['lat']
        reach_stats['p_wse'] = reach.metadata['wse']
        reach_stats['p_wse_var'] = reach.metadata['wse_var']
        reach_stats['p_width'] = reach.metadata['width']
        reach_stats['p_wid_var'] = reach.metadata['width_var']
        reach_stats['p_n_nodes'] = reach.metadata['n_nodes']
        reach_stats['p_dist_out'] = reach.metadata['dist_out']
        reach_stats['p_length'] = reach.metadata['reach_length']
        reach_stats['grand_id'] = reach.metadata['grod_id']
        reach_stats['n_chan_max'] = reach.metadata['n_chan_max']
        reach_stats['n_chan_mod'] = reach.metadata['n_chan_mod']
        reach_stats['ice_clim_f'] = reach.metadata['iceflag']
        reach_stats['river_name'] = reach.metadata['river_name']


        # Avoid letting fill value of -9999 from PRD propagate into outputs
        # (these variables are just passed through from PRD to RiverTile).
        def fill_if_was_fill(value, other_fill, fill):
            return value if value != other_fill else fill

        reach_stats['p_low_slp'] = fill_if_was_fill(
            reach.metadata['p_low_slp'], -9999, MISSING_VALUE_INT4)

        dsch_m_uc = reach.metadata['discharge_models']['unconstrained']
        dsch_m_c = reach.metadata['discharge_models']['constrained']

        # Set reach_q_b
        # Start by bitwise OR-ing all the node_q_b of the good (non-outlier)
        # nodes node_q_b qual flag together.
        reach_q_b = np.bitwise_or.reduce(river_reach.node_q_b[mask])

        # bitwise AND of node_q_b used to set bits 26-28
        reach_q_b_and = np.bitwise_and.reduce(river_reach.node_q_b[mask])

        # The following bits are not given by logically OR-ing all the
        # node_q_b flags together.

        # unset bit 1
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_SIG0_QUAL_SUSPECT

        # unset bit 4
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_BLOCK_WIDTH_SUSPECT

        # unset bit 9
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_FEW_SIG0_PIX

        # overwrite bits 10 / 11, and set if more than half of observed nodes
        # have QUAL_IND_FEW_AREA_PIX / QUAL_IND_FEW_WSE_PIX set
        nodes_have_few_area_pix = (
            river_reach.node_q_b[river_reach.valid_wse_mask] &
            SWOTRiver.products.rivertile.QUAL_IND_FEW_AREA_PIX > 0)
        nodes_have_few_wse_pix = (
            river_reach.node_q_b[river_reach.valid_wse_mask] &
            SWOTRiver.products.rivertile.QUAL_IND_FEW_WSE_PIX > 0)

        num_valid_wse_nodes = river_reach.valid_wse_mask.sum()
        if num_valid_wse_nodes > 0:
            frac_few_area_pix = sum(nodes_have_few_area_pix)/num_valid_wse_nodes
            frac_few_wse_pix = sum(nodes_have_few_wse_pix)/num_valid_wse_nodes
        else:
            frac_few_area_pix = 0
            frac_few_wse_pix = 0

        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_FEW_AREA_PIX
        if frac_few_area_pix >= 0.5:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_FEW_AREA_PIX

        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_FEW_WSE_PIX
        if frac_few_wse_pix >= 0.5:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_FEW_WSE_PIX

        # unset bit 15 / re-set it with partially_observed
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_PARTIAL_OBS
        if reach_stats['frac_obs'] < 0.5:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_PARTIAL_OBS

        # unset bits 23 and 24
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_WSE_OUTLIER
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_WSE_BAD

        # overwrite 25 / below_min_fit_points
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_NO_SIG0_PIX
        if mask.sum() < min_fit_points:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_MIN_FIT_POINTS

        # overwrite bit 26 / no_area_observations if ALL node_q_b bit 26 is set
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_NO_AREA_PIX
        if reach_q_b_and & SWOTRiver.products.rivertile.QUAL_IND_NO_AREA_PIX:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_NO_AREA_PIX

        # overwrite bit 27 / no_wse_observations if ALL node_q_b bit 27 is set
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_NO_WSE_PIX
        if reach_q_b_and & SWOTRiver.products.rivertile.QUAL_IND_NO_WSE_PIX:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_NO_WSE_PIX

        # overwrite bit 28 / no_observations if ALL node_q_b bit 28 is set
        reach_q_b &= ~SWOTRiver.products.rivertile.QUAL_IND_NO_OBS
        if reach_q_b_and & SWOTRiver.products.rivertile.QUAL_IND_NO_OBS:
            reach_q_b |= SWOTRiver.products.rivertile.QUAL_IND_NO_OBS

        # Create reach_q from reach_q_b
        reach_q = 0
        thresh_sus = 1
        thresh_deg = (
            SWOTRiver.products.rivertile.QUAL_IND_REACH_DEGRADED_THRESHOLD)
        thresh_bad = SWOTRiver.products.rivertile.QUAL_IND_REACH_BAD_THRESHOLD

        if thresh_sus <= reach_q_b < thresh_deg:
            reach_q = 1
        elif thresh_deg <= reach_q_b < thresh_bad:
            reach_q = 2
        elif reach_q_b >= thresh_bad:
            reach_q = 3
        reach_stats['reach_q_b'] = reach_q_b
        reach_stats['reach_q'] = reach_q
        reach_stats['xovr_cal_q'] = max(
            river_reach.xovr_cal_q[mask], default=2)

        # Compute discharge
        # TODO: slope2 or slope
        discharge_model_values = SWOTRiver.discharge.compute(
            reach, reach_stats['height'], reach_stats['height_u'],
            reach_stats['width'], reach_stats['width_u'],
            reach_stats['slope'], reach_stats['slope_u'], reach_q)

        # Put in discharge_model_values into reach_stats for output
        reach_stats.update(discharge_model_values)

        river_reach.metadata = reach_stats
        return river_reach

    def create_index_file(self):
        """Initializes the pixel cloud vector file"""
        with nc.Dataset(self.output_file, 'w') as ofp:
            ofp.createDimension('points', None)
            ofp.createVariable(
                'range_index', 'i4', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'azimuth_index', 'i4', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'node_id', 'i8', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'reach_id', 'i8', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'segmentation_label', 'i4', 'points',
                fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'distance_to_node', 'f4', 'points',
                fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'along_reach', 'f4', 'points', fill_value=FILL_VALUES['f4'])
            ofp.createVariable(
                'cross_reach', 'f4', 'points', fill_value=FILL_VALUES['f4'])
            ofp.createVariable(
                'latitude_vectorproc', 'f8', 'points',
                fill_value=FILL_VALUES['f8'])
            ofp.createVariable(
                'longitude_vectorproc', 'f8', 'points',
                fill_value=FILL_VALUES['f8'])
            ofp.createVariable(
                'height_vectorproc', 'f8', 'points',
                fill_value=FILL_VALUES['f8'])
            ofp.createVariable(
                'h_flg', 'i4', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'area_flg', 'i4', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'sig0_flg', 'i4', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'wse_class_flg', 'i4', 'points', fill_value=FILL_VALUES['i4'])
            ofp.createVariable(
                'geolocation_qual', 'u4', 'points',
                fill_value=FILL_VALUES['u4'])
            ofp.createVariable(
                'classification_qual', 'u4', 'points',
                fill_value=FILL_VALUES['u4'])
            ofp.createVariable(
                'sig0_qual', 'u4', 'points', fill_value=FILL_VALUES['u4'])


    def write_index_file(self, img_x, img_y, node_index, dst, along_reach,
                         cross_reach, reach_index, seg_lbl, lat, lon,
                         height, h_flg, area_flg, sig0_flg, wse_class_flg,
                         geolocation_qual, classification_qual, sig0_qual):
        """
        Write out the river obs indices for each pixel that get mapped to a
        node as well as the pixel cloud coordinates (range and azimuth, or
        original image coordinate [e.g., gdem along- and cross-track index])
        """
        # append the new data
        with nc.Dataset(self.output_file, 'a') as ofp:
            curr_len = len(ofp.variables['range_index'])
            new_len = curr_len + len(img_x)
            ofp.variables['range_index'][curr_len:new_len] = img_x
            ofp.variables['azimuth_index'][curr_len:new_len] = img_y
            ofp.variables['node_id'][curr_len:new_len] = node_index
            ofp.variables['reach_id'][curr_len:new_len] = reach_index
            ofp.variables['segmentation_label'][curr_len:new_len] = seg_lbl
            ofp.variables['distance_to_node'][curr_len:new_len] = dst
            ofp.variables['along_reach'][curr_len:new_len] = along_reach
            ofp.variables['cross_reach'][curr_len:new_len] = cross_reach
            # for improved geolocation
            ofp.variables['latitude_vectorproc'][curr_len:new_len] = lat
            ofp.variables['longitude_vectorproc'][curr_len:new_len] = lon
            ofp.variables['height_vectorproc'][curr_len:new_len] = height
            ofp.variables['h_flg'][curr_len:new_len] = h_flg.astype('i')
            ofp.variables['area_flg'][curr_len:new_len] = area_flg.astype('i')
            ofp.variables['sig0_flg'][curr_len:new_len] = sig0_flg.astype('i')
            ofp.variables['wse_class_flg'][curr_len:new_len] = (
                wse_class_flg.astype('i'))
            ofp.variables['geolocation_qual'][curr_len:new_len] = (
                geolocation_qual.astype('u4'))
            ofp.variables['classification_qual'][curr_len:new_len] = (
                classification_qual.astype('u4'))
            ofp.variables['sig0_qual'][curr_len:new_len] = (
                sig0_qual.astype('u4'))

        return

    def compute_enhanced_slope(
        self, river_reach_collection, river_reach, ireach,
        max_window_size, min_sigma, window_size_sigma_ratio, min_fit_points=3):
        """
        This function calculate enhanced reach slope from smoothed
        node height using Gaussian moving average.
        For more information, please see Dr. Renato Frasson's paper:
        https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017WR020887

        Inputs:
        river_reach: collection of node / reach quantities for this reach
        river_reach_collection: same for all reaches
        ireach: index into prior DB extracted reaches (self.reaches)
        max_window_size: the max of Gaussian window, default is 10km
        min_sigma : min sigma for gaussian averaging, default is 1km
        window_size_sigma_ratio : default is 5

        Output:
        enhanced_slope: enhanced reach slope
        """
        wse, wse_r_u, _, this_len, first_node, ss, mask, p_wse = \
            self.get_multi_reach_height(
            river_reach_collection, river_reach, ireach)

        if np.sum(river_reach.mask) < min_fit_points:
            enhanced_slope = MISSING_VALUE_FLT
            enhanced_slope_r_u = MISSING_VALUE_FLT
            enhanced_slope_u = MISSING_VALUE_FLT
        else:
            # handle indexing for masked reaches
            this_reach_edges = np.ma.flatnotmasked_edges(
                np.ma.masked_array(river_reach.wse, mask=~river_reach.mask))
            this_reach_len = (river_reach.node_ss[this_reach_edges[1]]
                              - river_reach.node_ss[this_reach_edges[0]])

            first_node_masked = first_node - np.sum(~mask[:first_node])
            last_node_masked = first_node_masked + np.sum(river_reach.mask) - 1

            # window size and sigma for Gaussian averaging
            window_size = np.min([max_window_size, this_reach_len])
            sigma = np.max([min_sigma, window_size / window_size_sigma_ratio])

            # smooth h_n_ave, and get slope
            slope = np.polyfit(ss[mask], wse[mask], 1)[0]
            heights_detrend = wse[mask] - slope * ss[mask]
            heights_smooth, wse_smooth_r_u, weight_matrix = \
                self.gaussian_averaging(
                    ss[mask], heights_detrend, wse_r_u[mask], window_size,
                    sigma)

            heights_smooth = heights_smooth + slope * ss[mask]
            enhanced_slope = (heights_smooth[last_node_masked] -
                              heights_smooth[first_node_masked])/this_reach_len

            # Compute covariance of first and last node smoothed heights
            cov_first_last = (
                weight_matrix[first_node_masked] * wse_r_u[mask] *
                weight_matrix[last_node_masked] * wse_r_u[mask]).sum() / (
                weight_matrix[first_node_masked].sum() *
                weight_matrix[last_node_masked].sum())

            enhanced_slope_r_u = (1/this_reach_len) * np.sqrt(
                    wse_smooth_r_u[first_node_masked]**2 +
                    wse_smooth_r_u[last_node_masked]**2 -
                    2 * cov_first_last)

            enhanced_slope_u = np.sqrt(
                SLOPE_SYS_UNCERT**2 + enhanced_slope_r_u**2)

        if np.isnan(enhanced_slope):
            enhanced_slope = MISSING_VALUE_FLT
            enhanced_slope_r_u = MISSING_VALUE_FLT
            enhanced_slope_u = MISSING_VALUE_FLT
        return enhanced_slope, enhanced_slope_r_u, enhanced_slope_u

    @staticmethod
    def gaussian_averaging(ss, wse, wse_r_u, window_size, sigma):
        """
        Gaussian smoothing of heights using distances
        ss:             along-flow distance
        wse:            water heights
        wse_r_u:        random component of uncertinity of water heights
        window_size:    size of data window to use for averaging
        sigma:          STD of Gaussian used for averaging

        outputs:
        wse_smooth:     smoothed water heights
        wse_smooth_r_u: random component of uncertinity of smoothed water
                        heights
        weight_matrix:  the full matrix of gaussian weights used
        """
        wse_smooth = np.zeros(wse.shape)
        wse_smooth_r_u = np.zeros(wse.shape)
        weight_matrix = np.zeros((wse.shape[0], wse.shape[0]))
        for ii, this_distance in enumerate(ss):

            # get data window
            mask = np.logical_and(
                np.abs(this_distance-ss) <= window_size / 2,
                ~np.isnan(wse))

            weights = scipy.stats.norm.pdf(
                this_distance, ss[mask], sigma)

            wse_smooth[ii] = (
                np.multiply(weights, wse[mask]).sum() /
                weights.sum())

            wse_smooth_r_u[ii] = np.sqrt(
                (np.multiply(weights, wse_r_u[mask])**2).sum() /
                (weights.sum())**2)

            weight_matrix[ii][mask] = weights

        return wse_smooth, wse_smooth_r_u, weight_matrix

    def get_reach_mask(self, ss, hh, ww, node_q, min_fit_points,
                       first_node=None, this_len=None):
        """
        Mask each node in an input reach based on whether or not it has a valid
        wse. Then mask for node-level wse outliers, if self.outlier_method is
        set. Returns a numpy array where good nodes are 1 and bad nodes are 0.

        Inputs:
            :param ss: along-reach distance for each node
            :param hh: height of each node
            :param ww: weights for each node height, to use in outlier
                flagging
            :param node_q: node quality flags for each node
            :param min_fit_points: minimum number of points needed to flag
                outliers
            :param first_node : 1st node index of the current reach, if
                neighbor reaches were used
            :param this_len : the length of current reach if neighbor reaches
                were used
        Outputs:
            :mask: reach mask where good nodes are 1 and bad nodes are 0
            :valid_wse_mask: reach mask where nodes with valid wse are 1,
                is the same as mask if outlier_method is None
        """
        mask = np.logical_and(hh > MIN_VALID_WSE, hh < MAX_VALID_WSE)
        if (ww[mask] == 0).any():
            LOGGER.warning(
                "get_reach_mask: Removing invalid wse_r_u (Inf) values!")
            mask = np.logical_and(mask, ww > 0)

        # Keep a copy of the valid WSE mask before outlier removal
        valid_wse_mask = mask.copy()

        SS = np.c_[ss, np.ones(len(ss), dtype=ss.dtype)]
        # throw away bad nodes upfront with > coarse threshold distance from
        # WLS fit line
        if mask.sum() > min_fit_points:
            fit = statsmodels.api.WLS(hh[mask], SS[mask],
                                      weights=ww[mask]).fit()
            wse_fit = fit.params[1] + fit.params[0] * SS[mask][:, 0]
            # use slope to define the coarse threshold (10 * slope[m/km]),
            # use 20 m if coarse threshold < 20 m
            coarse_thresh = np.maximum(abs(fit.params[0] * 1e3 * 10), 20)
            coarse_mask = abs(wse_fit - hh[mask]) < coarse_thresh
            mask[mask] = coarse_mask

        if first_node is not None and this_len is not None:
            # get reach node number and apply neighbor reach filters
            this_reach_node_num = mask[first_node:first_node+this_len].sum()
            if this_reach_node_num > min_fit_points:
                mask = self.neighbor_reach_filter(SS[:, 0], mask, node_q,
                                               first_node, this_len)
        else:
            this_reach_node_num = mask.sum()

        if (self.outlier_method is not None and
            this_reach_node_num > min_fit_points):
            # use WLS fit slope to define the reach absolute threshold:
            # self.outlier_abs_thresh * slope[m/km]
            # self.outlier_abs_thresh is the absolute outlier threshold for
            # reaches with 1m/km slope, use observed reach slope to scale
            # the threshold up or down, default to 2, tested on the calval
            # sites and a few more tiles
            fit = statsmodels.api.WLS(hh[mask], SS[mask],
                                      weights=ww[mask]).fit()
            reach_outlier_abs_thresh = abs(
                fit.params[0] * 1e3 * self.outlier_abs_thresh)

            if first_node is not None and this_len is not None:
                # use neighbor reaches
                masked_first_node = mask[0:first_node].sum()
                masked_this_len = mask[first_node:first_node+this_len].sum()
                mask = self.flag_outliers(hh[mask], SS[mask], ww[mask],
                                          node_q[mask], mask,
                                          reach_outlier_abs_thresh,
                                          masked_first_node, masked_this_len)

            else:
                # current reach only
                mask = self.flag_outliers(hh[mask], SS[mask], ww[mask],
                                          node_q[mask], mask,
                                          reach_outlier_abs_thresh)

        elif (self.outlier_method is not None and
              this_reach_node_num <= min_fit_points):
            mask = np.zeros(len(hh), dtype=bool)

        wse_outlier_mask = mask.copy()

        if mask.sum() > 0:
            pct_good_sus = 100 * (node_q[mask] < 2).sum() / mask.sum()
        else:
            pct_good_sus = 0

        if pct_good_sus > self.reach_pct_good_sus_thresh:
            mask[node_q >= 2] = False

        return mask, valid_wse_mask, wse_outlier_mask

    def flag_outliers(self, wse, flow_dist, weights, node_q, input_mask,
                      reach_outlier_abs_thresh, first_node=None,
                      this_len=None):
        """
        Flag wse outliers within a reach using iterative WLS slope algorithm
        wse: node wse
        flow_dist: node flow distance
        weights: weight for each node
        node_q: quality flag for each node
        input_mask: overall mask that wse and flow_dist are subset using
        reach_outlier_abs_thresh: reach outlier absolute threshold
        first_node: 1st node index of the current reach, if neighbor reaches
            were used
        this_len: the length of the current reach, if neighbor reaches were
            used

        valid methods: ['iterative_linear', 'piecewise_linear']

        iterative_linear outlier flagging method uses these class attributes:
        - self.outlier_rel_thresh: relative threshold, default 68% percentile
        - self.outlier_upr_thresh: upper limit, i.e. what percent of nodes are
                                   kept at the last iteration, default 80%
        - self.outlier_iter_num:   number of iterations

        piecewise_linear outlier flagging method uses these class attributes:
        - self.outlier_upr_thresh:       upper limit, i.e. what percent of
                                         nodes are kept at the last iteration,
                                         default 80%
        - self.outlier_breakpoint_min_dist: Minimum number of nodes between
                                            breakpoints, default 5
        - self.outlier_edge_min_dist:    Minimum number of nodes from edge to a
                                         breakpoint, default 10
        - self.outlier_n_boot:           The number of times to run the
                                         bootstrap restarting, positive int,
                                         default 10
        - self.outlier_iter_num:         Maximum iterations of Muggeo
                                         algorithms if not converged,
                                         positive int, default 30
        For more information of the piecewise linear fit:
        https://github.com/chasmani/piecewise-regression

        Three of the parameters below are used in both of the outlier flagging
        algorithms:
        - reach_outlier_abs_thresh
        - outlier_upr_thresh
        - outlier_iter_num

        outputs: returns a copy of the input mask array, further subset using
                 the outlier flagging method.
        """
        # Copy input_mask to mask (will be returned as output)
        mask = input_mask.copy()

        pct_good_sus = 100 * (node_q < 2).sum() / mask.sum()
        if pct_good_sus > self.reach_pct_good_sus_thresh:
            weights[node_q >= 2] = 0

        if self.outlier_method == 'iterative_linear':
            current_ind = self.iterative_linear(wse, flow_dist, weights,
                                                self.outlier_iter_num,
                                                self.outlier_rel_thresh,
                                                reach_outlier_abs_thresh,
                                                first_node, this_len)
            mask[mask] = current_ind
            
        elif self.outlier_method == 'piecewise_linear':
            current_ind = None
            # at least to have 10 node (2km) to run piecewise linear
            if len(wse) > 10:
                break_num = int(np.floor(len(wse)/10) - 1)
                while current_ind is None and break_num >=1:
                    current_ind = self.piecewise_linear(
                        flow_dist[:, 0], wse, node_q, break_num,
                        reach_outlier_abs_thresh, first_node, this_len)
                    break_num -= 1

            # if piecewise linear didn't converge at the end, use iterative
            # linear approach, the reaches at this stage are linear, so lower
            # the absolute threshold and iteration number
            if current_ind is None:
                current_ind = self.iterative_linear(wse, flow_dist, weights,
                                                    3, 68,
                                                    reach_outlier_abs_thresh,
                                                    first_node, this_len)
            mask[mask] = current_ind

        else:
            raise NotImplementedError(
                'Outlier flagging method %s is not supported!'%
                self.outlier_method)

        return mask

    @staticmethod
    def add_unpopulated_nodes(river_reach):
        """
        Adds unpopulated nodes to current reach for Bayes height reconstruction

        Parameters
        ----------
        river_reach : partially populated RiverReach instance with node
            quantities already computed

        Outputs
        ----------
        all_hh : wse with unpopulated nodes filled with NaN
        all_wse_r_u : wse_r_u with unpopulated nodes filled with NaN
        all_node_q : node_q with unpopulated nodes filled with NaN
        all_mask : mask with unpopulated nodes filled with NaN
        """
        all_hh = np.empty(river_reach.prior_node_ss.shape) * np.nan
        all_wse_r_u = np.empty(river_reach.prior_node_ss.shape) * np.nan
        all_node_q = np.empty(river_reach.prior_node_ss.shape) * np.nan
        all_mask = np.zeros(river_reach.prior_node_ss.shape, dtype=bool)

        all_hh[river_reach.populated_nodes] = river_reach.wse
        all_wse_r_u[river_reach.populated_nodes] = river_reach.wse_r_u
        all_node_q[river_reach.populated_nodes] = river_reach.node_q
        all_mask[river_reach.populated_nodes] = river_reach.mask
        return all_hh, all_wse_r_u, all_node_q, all_mask


    def get_multi_reach_height(
            self, river_reach_collection, river_reach, ireach):
        """
        Handles the upstream and downstream reaches from the PRD and returns
        heights and lengths over multiple reaches for enhanced/Bayes slope
        calculations. Includes checks for dam reaches and edge node proximity.

        Parameters
        ----------
        river_reach_collection : List of RiverReach instances for tile
        river_reach : partially populated RiverReach instance with node
            quantities already computed
        reach : Reach instance
            One of the reaches from ReachExtractor.
        ireach : int
            Index in the list of reaches extracted for this scene.

        Outputs
        ----------
        wse : node-level wse over (valid) upstream, current, and (valid)
              downstream reaches.
        wse_r_u : node-level wse_r_u over (valid) upstream, current, and
                  (valid) downstream reaches.
        this_len : number of nodes in current reach
        first_node : index of first node in current reach following
                     concatenation of multiple reaches
        ss : node-level distance over (valid) upstream, current, and (valid)
             downstream reaches.
        """
        this_id = river_reach.reach_indx[0]
        other_ids = [item.reach_indx[0] for item in river_reach_collection]

        # get up/dn id from prior db
        prior_s = river_reach.prior_node_ss

        # skip (disconnected lake, dam, unrealizable topology) reaches
        skip_types = [3, 4, 5]

        # use the first good adjacent reach that is not skipped for height
        # smoothing
        # TO-DO: add handling for multiple upstream/downstream reaches
        for up_id_try in self.reaches[ireach].metadata['rch_id_up'][:, 0]:
            if not up_id_try % 10 in skip_types:
                break

        for dn_id_try in self.reaches[ireach].metadata['rch_id_dn'][:, 0]:
            if not dn_id_try % 10 in skip_types:
                break

        # In all these dicts: -1 -- downstream, +1 -- upstream
        prd_rch = {}
        prd_rch[0] = self.reaches.reach[
            np.where(self.reaches.reach_idx == this_id)[0][0]]
        prd_is_good = {-1: False, 1: False}
        prd_delta = {}
        adj_rch = {}

        for side in [-1, 1]:
            for id_try in [dn_id_try, up_id_try]:
                try:
                    # index in observed reaches
                    other_idx = np.where(other_ids == id_try)[0][0]

                    # get PRD reach
                    try_prd_rch = self.reaches.reach[
                        np.where(self.reaches.reach_idx == id_try)[0][0]]

                except IndexError:
                    # cannot find adjacent reach with id_try in either
                    # PRD or observed reaches
                    continue

                # only add when neighbor has more than 5 non-bad nodes
                adj_rch_obs = river_reach_collection[other_idx]
                adj_rch_mask = adj_rch_obs.node_q < 3
                if side == -1 and adj_rch_mask.sum() > 5:
                    # side is downstream of current reach
                    dx = try_prd_rch.x[-1] - prd_rch[0].x[0]
                    dy = try_prd_rch.y[-1] - prd_rch[0].y[0]
                elif side == 1 and adj_rch_mask.sum() > 5:
                    # side is upstream of current reach
                    dx = try_prd_rch.x[0] - prd_rch[0].x[-1]
                    dy = try_prd_rch.y[0] - prd_rch[0].y[-1]
                else:
                    dx = np.nan
                    dy = np.nan

                delta = np.sqrt(dx**2+dy**2)
                if delta < 300:
                    prd_is_good[side] = True
                    prd_delta[side] = delta
                    adj_rch[side] = river_reach_collection[other_idx]
                    prd_rch[side] = try_prd_rch

        # Build up array of data to be smoothed from upstream to
        # downstream.  Adjust along-reach to be cumulative across
        # reach boundaries.
        first_node = 0
        ss = np.array([])
        p_wse = np.array([])
        wse = np.array([])
        wse_r_u = np.array([])
        node_q = np.array([])
        mask = np.array([], dtype=bool)
        # if upstream PRD reach is usable
        if prd_is_good[1]:
            ss = np.concatenate([adj_rch[1].node_ss, ss])
            p_wse = np.concatenate([adj_rch[1].p_wse, p_wse])
            wse = np.concatenate([adj_rch[1].wse, wse])
            wse_r_u = np.concatenate([adj_rch[1].wse_r_u, wse_r_u])
            node_q = np.concatenate([adj_rch[1].node_q, node_q])
            mask = np.concatenate([adj_rch[1].mask, mask])

        if self.bayes_slope_use_all_nodes:
            # add unpopulated nodes for current reach
            this_hh, this_wse_r_u, this_node_q, this_mask = \
                self.add_unpopulated_nodes(river_reach)
            this_ss, this_p_wse = \
                river_reach.prior_node_ss, river_reach.p_wse_all
        else:
            # use populated nodes only
            (this_hh, this_wse_r_u, this_node_q, this_mask, this_ss,
             this_p_wse) = (river_reach.wse, river_reach.wse_r_u,
                            river_reach.node_q, river_reach.mask,
                            river_reach.node_ss, river_reach.p_wse)

        ss = np.concatenate([this_ss, ss+prior_s[-1]])
        p_wse = np.concatenate([this_p_wse, p_wse])
        wse = np.concatenate([this_hh, wse])
        wse_r_u = np.concatenate([this_wse_r_u, wse_r_u])
        node_q = np.concatenate([this_node_q, node_q])
        mask = np.concatenate([this_mask, mask])
        this_len = len(this_ss)

        # if downstream PRD reach is usable
        if prd_is_good[-1]:
            downstream_prior_s = adj_rch[-1].prior_node_ss
            first_node = first_node + len(adj_rch[-1].wse)
            ss = np.concatenate([adj_rch[-1].node_ss,
                                 ss+downstream_prior_s[-1]])
            p_wse = np.concatenate([adj_rch[-1].p_wse, p_wse])
            wse = np.concatenate([adj_rch[-1].wse, wse])
            wse_r_u = np.concatenate([adj_rch[-1].wse_r_u, wse_r_u])
            node_q = np.concatenate([adj_rch[-1].node_q, node_q])
            mask = np.concatenate([adj_rch[-1].mask, mask])

        return (wse, wse_r_u, node_q, this_len, first_node, ss, mask,
                p_wse)

    @staticmethod
    def neighbor_reach_filter(flow_dist, mask, node_q, first_node, this_len):
        """
        filter out nodes in neighbor reaches beyond big gaps (5 nodes)
        the gap includes unobserved nodes, masked nodes, and bad nodes

        Parameters
        ----------
        flow_dist: node flow distance
        mask: node mask
        node_q: node quality flag
        first_node: 1st node index of the current reach, if neighbor reaches
            were used
        this_len: length of current reachh, if neighbor reaches were used

        Outputs
        ----------
        mask: updated node mask after the neighbor reach filter
        """
        # mask out neighbor bad nodes
        node_q_mask = np.ones(len(node_q), dtype=bool)
        node_q_mask[0:first_node] = node_q[0:first_node] < 3
        node_q_mask[first_node+this_len:] = node_q[first_node+this_len:] < 3
        mask = np.logical_and(node_q_mask, mask)
        ss_masked = flow_dist[mask]
        node_dist = np.zeros(len(flow_dist))
        # find and remove the last True element,
        # since distance array is one element less
        dist_mask = mask.copy()
        ind = np.max([i for i, val in enumerate(dist_mask) if val])
        dist_mask[ind] = False
        # compute the distance between populated unmasked nodes
        node_dist[dist_mask] = [abs(t - s) for s, t in
                                zip(ss_masked, ss_masked[1:])]
        obs_gaps = np.where(np.array(node_dist) > 1000)[0]
        # if neighbor has big gaps (> 5 nodes),
        # remove neighbor nodes beyond the gap
        if len(obs_gaps) > 0:
            with contextlib.suppress(ValueError):
                gap_dn = np.max(obs_gaps[obs_gaps < first_node])
                mask[0:gap_dn + 1] = False
            with contextlib.suppress(ValueError):
                gap_up = np.min(obs_gaps[obs_gaps >= first_node+this_len])
                mask[gap_up:] = False
        return mask

    def optimal_reconstruct(
            self,
            river_reach_collection,
            river_reach,
            ireach,
            ss,
            wse,
            wse_r_u,
            mask,
            min_fit_points=2,
            prior_cov_method='exponential',
            full_noise_cov=False,
            method='Bayes'):
        """
        This function estimates the optimal reconstruction estimator under
        certain assumptions given by the options.
        Inputs:
        river_reach_collection : List of RiverReach instances for tile
        river_reach : partially populated RiverReach instance with node
                      quantities already computed
        ireach : int
            Index in the list of reaches extracted for this scene.
        ss         : Node-level distance from prior database
        wse        : Measured Node wse (with masked values for missing nodes)
        wse_r_u    : Node-wise random wse uncertainty
        mask       : True if node-level height is good, False where it is bad
        Options:
        prior_cov_method : independent, exponential
                           This option controls how spatial structure is
                           imposed by the prior.
        full_noise_cov : True, or False
                         This option sets the noise covariance and the sampling
                         operator to assume that all nodes are observed (but
                         weighted appropriately for the unobserved nodes).
        method = Bayes
           Bayes        : Bayes estimate of wse given prior and
                          prior_cov method.
        """

        if self.use_multiple_reaches:
            wse, wse_r_u, _, this_len, first_node, ss, mask, prior_wse = \
                self.get_multi_reach_height(
                    river_reach_collection, river_reach, ireach)
            # get the multi-reach mask
            ss = ss - np.mean(ss)
        else:
            first_node = 0
            this_len = len(ss)

        # define vectors b and c for uncertainty estimates later
        this_reach_mask_b = np.zeros_like(ss)
        this_reach_mask_b[first_node:first_node+this_len] = 1
        this_reach_mask_b = this_reach_mask_b / np.sum(this_reach_mask_b)
        first_and_last_node_c = np.zeros_like(ss)
        reach_length = ss[first_node + this_len - 1] - ss[first_node]
        first_and_last_node_c[first_node] = -1 / reach_length
        first_and_last_node_c[first_node+this_len-1] = 1 / reach_length

        # create a wse prior if flagged
        if self.prior_wse_method == 'fit':
            prior_wse = None
        elif self.prior_wse_method == 'prd':
            if self.use_multiple_reaches:
                # prior_wse already set above
                pass
            else:
                # get prior wse from PRD
                if self.bayes_slope_use_all_nodes:
                    prior_wse = river_reach.p_wse_all
                else:
                    prior_wse = river_reach.p_wse
        else:
            raise Exception('Prior wse method %s is not an implemented option '
                            'for the reconstruction' % self.prior_wse_method)

        if prior_wse is None:
            # If no prior wse is given, use obs wse to make linear fit prior
            ww = 1/wse_r_u**2
            SS = np.c_[ss, np.ones(len(ss), dtype=ss.dtype)]
            wse_fit = statsmodels.api.WLS(wse[mask], SS[mask],
                                          weights=ww[mask]).fit()
            prior_wse = wse_fit.predict(SS)

        # get the sampling operator
        # find where the data is not masked out or NaN
        ind = np.where(np.logical_not(np.logical_and(mask, np.isfinite(wse))))
        # get vector with 1 for valid data elements
        h = np.ones(np.shape(wse))
        h[ind] = 0
        # make full sampling matrix
        H = np.diag(h)
        # now remove the zero rows
        num_removed = 0
        for k, val in enumerate(h):
            if val == 0:
                row = k - num_removed
                H = np.delete(H, row, 0)
                num_removed = num_removed + 1
        # get the covariance matrices
        msk = np.logical_and(mask, np.isfinite(wse))
        Rv = self.get_noise_autocov(wse, wse_r_u, mask, full_noise_cov)
        if prior_cov_method == 'independent':
            # assume no spatial structure
            Ry0 = np.identity(len(ss))
        elif prior_cov_method == 'exponential':
            # get the signal covariance assuming an exponential random process
            Ry0 = np.zeros((len(ss), len(ss)))
            for k, d0 in enumerate(ss):
                t = ss - d0
                Ry0[k, :] = np.exp(-np.abs(t) / self.char_length_tau)
        else:
            raise Exception('Covariance model %s is not an implemented option '
                            'for the reconstruction' % prior_cov_method)
        # scale the covariance to trade-off noise.vs "spectral resolution"
        Ry = Ry0 / np.max(Ry0) * self.prior_unc_alpha ** 2
        # compute the optimal wse reconstruction filter
        if method == 'Bayes':
            # get the bayes estimate
            K, K_bar, A_inv = self.compute_bayes_estimator(Ry, Rv, H)
        else:
            raise Exception('Reconstruction method %s is not an implemented '
                            'option for the reconstruction' % method)
        # handle the missing node wse measurements
        if full_noise_cov:
            # wse_reg = prior_wse.copy()
            wse_reg = np.zeros(np.shape(wse))
            wse_reg[msk] = wse[msk]
        else:
            wse_reg = wse[msk]
        # apply the wse estimator(s) to the measurement term
        wse_out0 = np.matmul(K, wse_reg)
        # apply the prior term
        wse_out = wse_out0 + np.matmul(K_bar, prior_wse)
        wse_out_std = np.sqrt(np.diag(A_inv))  # node-level height uncertainty
        height_u = np.sqrt(this_reach_mask_b @ A_inv @ this_reach_mask_b.T)
        slope_u = np.sqrt(first_and_last_node_c @ A_inv @ np.atleast_2d(
            first_and_last_node_c).T)
        if self.use_multiple_reaches:
            wse_out = wse_out[first_node:first_node+this_len]
            wse_out_std = wse_out_std[first_node:first_node+this_len]
        # Return node WSEs & stdevs, and reach height_u and slope_u as scalars
        return wse_out, wse_out_std, height_u.item(), slope_u.item()

    @staticmethod
    def compute_bayes_estimator(Ry, Rv, H):
        """
        Implements the bayes estimator of the signal/parameters
        Ry is the signal (or parameter) covariance
        Rv is the measurement noise covariance
        H is the sampling operator (or sampling then basis projection operator)
        A_inv is the posterior covariance (post_cov)

        These are the equations implemented
        K = (Ry^-1 + H.T Rv^-1 H)^-1 H.T Rv^-1
        K_bar = (Ry^-1 + H.T Rv^-1 H)^-1 Ry^-1 (if non-zero mean of prior wse)
        A = (Ry^-1 + H.T Rv^-1 H))
        """
        Ry_inv = np.linalg.pinv(Ry)
        Rv_inv = np.linalg.pinv(Rv)
        A = Ry_inv + H.T @ Rv_inv @ H
        post_cov = np.linalg.pinv(A)  # A_inv
        K = post_cov @ H.T @ Rv_inv
        K_bar = post_cov @ Ry_inv

        return K, K_bar, post_cov

    @staticmethod
    def get_noise_autocov(wse, wse_r_u, mask, full=False):
        msk = np.logical_and(mask, np.isfinite(wse))
        if full:
            reg_inf = 1e10
            wse_r_u_reg = np.zeros(np.shape(wse_r_u)) + reg_inf
            wse_r_u_reg[msk] = wse_r_u[msk]
            wse_r_u_reg[wse_r_u_reg < 0] = reg_inf
        else:
            wse_r_u_reg = wse_r_u[msk]
        return np.diag(wse_r_u_reg)

    def iterative_linear(self, wse, flow_dist, weights, outlier_iter_num,
                         outlier_rel_thresh, reach_outlier_abs_thresh,
                         first_node=None, this_len=None):
        """
        Flag wse outliers within a reach using iterative WLS slope algorithm
        wse: node wse
        flow_dist: node flow distance
        weights: weight for each node
        outlier_iter_num: number of iterations
        outlier_rel_thresh: relative threshold, default 68% percentile
        reach_outlier_abs_thresh: reach outlier absolute threshold
        first_node: 1st node index of the current reach, if neighbor reaches
            were used
        this_len: the length of the current reach, if neighbor reaches were
            used

        iterative_linear outlier flagging method uses these class attributes:

        - self.outlier_upr_thresh: upper limit, i.e. what percent of nodes are
                                   kept at the last iteration, default 80%

        outputs: indices of node not marked as outliers
        """
        if np.ma.isMaskedArray(flow_dist):
            flow_dist= flow_dist.data
        if np.ma.isMaskedArray(wse):
            wse = wse.data
        if np.ma.isMaskedArray(weights):
            weights = weights.data
        # initial fit OLS
        fit = statsmodels.api.OLS(wse, flow_dist).fit()
        r_slp = fit.params[0]
        r_wse = fit.params[1]
        wse_fit = r_wse + r_slp * flow_dist[:, 0]

        for i in range(outlier_iter_num):
            metric = abs(wse_fit - wse)
            if np.ma.isMaskedArray(metric):
                metric = metric.data
            metric_one_sigma = np.percentile(
                metric, outlier_rel_thresh)
            metric_upr_thresh = np.percentile(metric, self.outlier_upr_thresh)
            icond = np.logical_or(metric < reach_outlier_abs_thresh,
                                  metric < metric_one_sigma)
            if sum(icond) / len(metric) * 100 <= self.outlier_upr_thresh:
                current_ind = metric < metric_upr_thresh
            else:
                current_ind = icond
            if all(np.isfinite(weights)):
                fit = statsmodels.api.WLS(wse[current_ind],
                                          flow_dist[current_ind],
                                          weights=weights[current_ind],
                                          missing='drop').fit()
                r_slp = fit.params[0]
                r_wse = fit.params[1]
            else:
                fit = statsmodels.api.OLS(wse[current_ind],
                                          flow_dist[current_ind],
                                          missing='drop').fit()
                r_slp = fit.params[0]
                r_wse = fit.params[1]
            wse_fit = r_wse + r_slp * flow_dist[:, 0]

        metric = abs(wse_fit - wse)
        # if neighbor reaches used in outlier flagging and current reach first
        # and last node indices are provided, apply thresholds only on current
        # reach
        if first_node is not None and this_len is not None:
            this_metric_upr_thresh = (
                np.percentile(metric[first_node:first_node+this_len],
                              self.outlier_upr_thresh))
            # if metric at upper threshold (default 80%) is bigger than
            # absolute threshold (default 1.5 m), use upper threshold to avoid
            # flagging too many nodes, else remove > absolute threshold nodes
            if reach_outlier_abs_thresh <= this_metric_upr_thresh:
                current_ind[first_node:first_node+this_len] = (
                        metric[first_node:first_node+this_len] <
                        this_metric_upr_thresh)
            else:
                current_ind[first_node:first_node+this_len] = (
                    metric[first_node:first_node+this_len] <
                    reach_outlier_abs_thresh)
        else:
            metric_upr_thresh = np.percentile(metric, self.outlier_upr_thresh)
            if reach_outlier_abs_thresh <= metric_upr_thresh:
                current_ind = metric < metric_upr_thresh
            else:
                current_ind = metric < reach_outlier_abs_thresh
        return current_ind

    def piecewise_linear(self, x, y, quality, n_breakpoints,
                         reach_outlier_abs_thresh, first_node=None,
                         this_len=None):
        """
        Piecewise linear outlier detector

        Parameters
        ----------
        x : Node-level flow distance
        y : Measured Node wse
        quality : Node quality flag for adding weight to each node
        n_breakpoints : Number of breakpoints
        reach_outlier_abs_thresh: reach outlier absolute threshold
        first_node: 1st node index of the current reach, if neighbor reaches
            were used
        this_len: the length of the current reach, if neighbor reaches were
            used
        self.outlier_breakpoint_min_dist: Minimum number of nodes between
            breakpoints, default 5 (1km)
            In older version this parameter was defined as proportion of a
            reach, e.g. 0.1, 10% percent of the reach
        self.outlier_edge_min_dist: Minimum number of nodes from edge to a
            breakpoint, default 10 (2km)
            In older version this parameter was defined as proportion of a
            reach, e.g. 0.1, 10% percent of the reach
        For more information of the piecewise linear fit:
        https://github.com/chasmani/piecewise-regression

        Outputs
        ----------
        returns the input mask array
        """
        # Remove masks from arrays if present before using quantile,
        # percentile, ...etc
        if np.ma.isMaskedArray(x):
            x = x.data
        if np.ma.isMaskedArray(y):
            y = y.data
        if np.ma.isMaskedArray(quality):
            quality = quality.data

        # duplicate data to make good observations to have more weight,
        # if fraction of good and suspect nodes > reach_pct_good_sus_thresh:
        # good(*3), suspect(*2)
        # else: good (*4), suspect (*3), degraded (*2) and bad nodes (*1)
        good_ind = quality == 0
        suspect_ind = quality == 1
        degraded_ind = quality == 2
        x_good = np.concatenate((x[good_ind], x[good_ind], x[good_ind]),
                                axis=None)
        y_good = np.concatenate((y[good_ind], y[good_ind], y[good_ind]),
                                axis=None)
        x_suspect = np.concatenate((x[suspect_ind], x[suspect_ind]),
                                   axis=None)
        y_suspect = np.concatenate((y[suspect_ind], y[suspect_ind]),
                                   axis=None)
        x_degraded = x[degraded_ind]
        y_degraded = y[degraded_ind]

        pct_good_sus = 100 * (quality < 2).sum() / len(x)
        if pct_good_sus > self.reach_pct_good_sus_thresh:
            x_dup = np.concatenate((x_good, x_suspect), axis=None)
            y_dup = np.concatenate((y_good, y_suspect), axis=None)
        else:
            x_dup = np.concatenate((x, x_good, x_suspect, x_degraded),
                                   axis=None)
            y_dup = np.concatenate((y, y_good, y_suspect, y_degraded),
                                   axis=None)

        # in case the older version of param file is used
        if 0 < self.outlier_breakpoint_min_dist < 1:
            min_dist_between_bps = self.outlier_breakpoint_min_dist
        else:
            # compute min dist between breakpoints (%), default 5 nodes,
            # use the original x array, not the multiplied array
            min_dist_between_bps = self.outlier_breakpoint_min_dist / len(x)
        # in case the older version of param file is used
        if 0 < self.outlier_edge_min_dist < 1:
            min_dist_edge = self.outlier_edge_min_dist
        else:
            # compute min dist (%) between edge to breakpoints, default 10
            # nodes, use the original x array, not the multiplied array
            min_dist_edge = self.outlier_edge_min_dist / len(x)
        min_allowed_bp = np.quantile(x, min_dist_edge)
        max_allowed_bp = np.quantile(x, 1 - min_dist_edge)
        # generate pseudo-random seeds for piecewise_regression package
        seed = int((np.nanmean(y_dup) - MIN_VALID_WSE)*1e6) % 2**32
        np.random.seed(seed)
        start_values = np.linspace(
            min_allowed_bp, max_allowed_bp, num=n_breakpoints)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pw_fit = piecewise_regression.Fit(
                x_dup, y_dup, start_values=start_values,
                n_breakpoints=n_breakpoints,
                min_distance_between_breakpoints=min_dist_between_bps,
                min_distance_to_edge=min_dist_edge,
                n_boot=self.outlier_n_boot,
                max_iterations=self.outlier_iter_num)

        if pw_fit.best_muggeo is not None:
            final_params = pw_fit.best_muggeo.best_fit.raw_params
            breakpoints = pw_fit.best_muggeo.best_fit.next_breakpoints
            # Extract what we need from params etc
            intercept_hat = final_params[0]
            alpha_hat = final_params[1]
            beta_hats = final_params[2:2 + len(breakpoints)]

            # Build the fit plot segment by segment. Betas are defined as
            # difference in gradient from previous section
            y_hat = intercept_hat + alpha_hat * x
            for bp_count in range(len(breakpoints)):
                y_hat += (
                    beta_hats[bp_count] *
                    np.maximum(x - breakpoints[bp_count], 0))
            metric = abs(y - y_hat)
            # if neighbor reaches used in outlier flagging and current reach
            # first and last node indices are provided, apply thresholds only
            # on current reach
            if first_node is not None and this_len is not None:
                this_metric_upr_thresh = np.percentile(
                    metric[first_node:first_node+this_len],
                    self.outlier_upr_thresh)
                # if metric at upper threshold (default 80%) is greater than
                # absolute threshold, use upper threshold to avoid flagging
                # too many nodes, else remove > absolute threshold nodes
                current_ind = np.zeros(len(metric), dtype=bool)
                if reach_outlier_abs_thresh <= this_metric_upr_thresh:
                    current_ind[first_node:first_node+this_len] = (
                        metric[first_node:first_node+this_len] <
                        this_metric_upr_thresh)
                else:
                    current_ind[first_node:first_node+this_len] = (
                            metric[first_node:first_node+this_len] <
                            reach_outlier_abs_thresh)
            else:
                metric_upr_thresh = np.percentile(
                    metric, self.outlier_upr_thresh)
                if reach_outlier_abs_thresh <= metric_upr_thresh:
                    current_ind = metric < metric_upr_thresh
                else:
                    current_ind = metric < reach_outlier_abs_thresh
        else:
            current_ind = None
        return current_ind

