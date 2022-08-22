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
import piecewise_regression

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

from Centerline.Centerline import Centerline, CenterLineException
from scipy.ndimage.morphology import binary_dilation

LOGGER = logging.getLogger(__name__)

REACH_WSE_SYS_UNCERT = 0.09  # m
SLOPE_SYS_UNCERT = 0.000003  # m/m

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
        Minimum distance from edge to a breakpoint, as a proportion of the input data
        range for piecewise linear outlier detection algorithm.
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
                 **proj_kwds):

        self.trim_ends = trim_ends
        self.input_file = os.path.split(swotL2_file)[-1]
        self.output_file = output_file  # index file
        self.subsample_factor = subsample_factor
        self.slope_method = slope_method
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

        self.create_index_file()
        if np.ma.is_masked(self.lat):
            mask = self.lat.mask
        else:
            mask = np.zeros(len(self.lat), dtype=np.bool)
        self.h_noise = self.get(height_kwd)
        if np.ma.is_masked(self.h_noise):
            mask = mask | self.h_noise.mask

        # skip NaNs in dheight_dphase
        good = ~mask
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
            ['pole_tide', pole_tide_kwd]]

        for dset_name, keyword in datasets2load:
            try:
                value = self.get(keyword)
                # hack ifgram re/im parts into complex dtype
                if dset_name is 'ifgram':
                    value = value[:,0] + 1j* value[:,1]
            except KeyError:
                value = None
            setattr(self, dset_name, value)

        try:
            self.looks_to_efflooks = self.getatt(looks_to_efflooks_kwd)
            if self.looks_to_efflooks == 'None':
                self.looks_to_efflooks = 1.75  # set to default value
        except KeyError:
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
                    print("could not find correct pixel area parameters")

        # need to scale pixel area by the subsampling factor if subsampling
        if self.subsample_factor > 1:
            self.pixel_area = self.pixel_area * self.subsample_factor

        if fractional_inundation_kwd is None:  # all water pixels are inundated
            self.fractional_inundation = None
            self.inundated_area = self.pixel_area

        else:
            self.fractional_inundation = self.get(fractional_inundation_kwd)
            self.inundated_area = self.pixel_area
            for i, k in enumerate(class_list):
                if use_fractional_inundation[i]:
                    # TO-DO: validate fractional inundation result here
                    index = self.klass == k
                    self.inundated_area[index] = (self.pixel_area[
                        index] * self.fractional_inundation[index])

        # Segment the image into disjoint features
        if len(self.use_segmentation) == len(self.class_list):
            self.isWater = np.zeros(np.shape(self.h_noise))
            for i, k in enumerate(class_list):
                if self.use_segmentation[i]:
                    index = self.klass == k
                    self.isWater[index] = 1
            self.segment_water_class(self.preseg_dilation_iter)
        else:
            self.seg_label = None

        # set which class pixels to use for heights (should be modified
        # by layover flag later)
        if len(self.use_heights) == len(self.class_list):
            self.h_flg = np.zeros(np.shape(self.h_noise))
            for i, k in enumerate(class_list):
                if self.use_heights[i]:
                    index = self.klass == k
                    self.h_flg[index] = 1
        else:
            self.h_flg = None

        LOGGER.debug('Data loaded')

        # Initialize the list of reaches to output to product
        self.river_reach_collection = collections.OrderedDict()

        self.flatten_interferogram()
        self.nc.close()

    def flatten_interferogram(self):
        """Flattens the pixel cloud interferogram"""
        # range index is self.img_x, azi is self.img_y
        try:
            import cnes.modules.geoloc.lib.geoloc as geoloc
            from cnes.common.lib.my_variables import \
                GEN_RAD_EARTH_EQ, GEN_RAD_EARTH_POLE

        except ImportError:
            print("Can't flatten interferogram as swotCNES not found.")
            return

        try:
            tvp_plus_y_antenna_xyz = (
                self.nc['tvp']['plus_y_antenna_x'][:],
                self.nc['tvp']['plus_y_antenna_y'][:],
                self.nc['tvp']['plus_y_antenna_z'][:])
            tvp_minus_y_antenna_xyz = (
                self.nc['tvp']['minus_y_antenna_x'][:],
                self.nc['tvp']['minus_y_antenna_y'][:],
                self.nc['tvp']['minus_y_antenna_z'][:])
            pixc = {'tvp': {'time': self.nc['tvp']['time'][:]},
                    'pixel_cloud': {'illumination_time': 
                     self.nc['pixel_cloud']['illumination_time'][:]}}
        except KeyError:
            print("Can't flatten interferogram as TVP not in pixel cloud.")
            return

        hgt_2d = np.nan*np.ones((np.max(self.img_y)+1, np.max(self.img_x)+1))
        hgt_2d[self.img_y, self.img_x] = self.h_noise

        # do a 2d median filter on the heights
        hgt_filt = scipy.ndimage.generic_filter(hgt_2d, np.nanmedian, size=11)

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

    def get_reaches(self, reach_db_path, clip=False, clip_buffer=0.1,
                    day_of_year=None):
        """Get all of the reaches using a ReachExtractor."""
        self.clip = clip
        self.clip_buffer = clip_buffer

        self.reaches = RiverObs.ReachDatabase.ReachExtractor(
            reach_db_path, self, clip=clip, clip_buffer=clip_buffer,
            day_of_year=day_of_year)

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
                        use_width_db=False,
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
        use_width_db : bool, default False
            Use the width data base for setting widths?
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
        # assign the reaches
        if self.use_ext_dist_coef:
            river_obs_list, reach_idx_list, ireach_list = \
                self.assign_reaches_two_pass(
                    scalar_max_width, minobs, use_width_db, ds)
        else:
            river_obs_list, reach_idx_list, ireach_list = \
                self.assign_reaches(scalar_max_width, minobs, use_width_db, ds)

        river_reach_collection = []
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:

            if use_width_db:
                max_width = self.get_max_width_from_db(reach_idx)
                LOGGER.debug('max_width read')

            else:
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
                use_width_db=use_width_db,
                max_width=max_width,
                ds=ds,
                refine_centerline=refine_centerline,
                smooth=smooth,
                alpha=alpha,
                max_iter=max_iter)
            # compute node mask to be used for reach aggregation
            ww = 1 / (river_reach.wse_r_u**2)
            river_reach.mask = self.get_reach_mask(river_reach.node_ss,
                                                   river_reach.wse,
                                                   ww,
                                                   min_fit_points)
            river_reach_collection.append(river_reach)

            LOGGER.debug('reach processed')

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
                    enhanced_slope = self.compute_enhanced_slope(
                        river_reach_collection, river_reach, ireach,
                        max_window_size=max_window_size,
                        min_sigma=min_sigma,
                        window_size_sigma_ratio=window_size_sigma_ratio,
                        min_fit_points=min_fit_points
                    )
                else:
                    enhanced_slope = MISSING_VALUE_FLT

                # add enhanced slope to river reach outputs
                out_river_reach.metadata['slope2'] = enhanced_slope
                out_river_reach_collection.append(out_river_reach)

        return out_river_reach_collection

    def assign_reaches(self,
                       scalar_max_width,
                       minobs=10,
                       use_width_db=False,
                       ds=None,
                       second_pass=False):
        """
        Assigns pixels to nodes for every reach.
        """
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
                    minobs=minobs, second_pass=second_pass)
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

            LOGGER.debug('Reach %d/%d Reach index: %d' %(
                i_reach + 1, self.reaches.nreaches, reach_idx))

            seg_label = self.seg_label.copy()

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
                    second_pass=second_pass)

            except CenterLineException as e:
                print("CenterLineException: ", e)
                continue

            # Add width per node to centerline and re-init IteratedRiverObs.
            # Search width is the width it uses to always include (1/2 on
            # each side of centerline).
            search_width = 2 * (
                self.reaches[i_reach].width * self.reaches[i_reach].wth_coef)
            river_obs.add_centerline_obs(
                self.reaches[i_reach].x, self.reaches[i_reach].y,
                search_width, 'max_width')
            river_obs.reinitialize()

            if len(river_obs.x) == 0:
                LOGGER.debug(
                    'No observations mapped to nodes in this reach')
                continue

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

            # set in_channel mask to exlude the nodes to drop
            river_obs.in_channel[reach_ind != ii] = False

            # Drop pixels that were double-assigned to reaches and
            # recompute things set in RiverObs constructor
            river_obs.index = river_obs.index[mask_keep]
            river_obs.d = river_obs.d[mask_keep]
            river_obs.x = river_obs.x[mask_keep]
            river_obs.y = river_obs.y[mask_keep]
            river_obs.s = river_obs.s[mask_keep]
            river_obs.n = river_obs.n[mask_keep]
            river_obs.nedited_data = len(river_obs.d)
            river_obs.populated_nodes, river_obs.obs_to_node_map = \
                river_obs.get_obs_to_node_map(river_obs.index,
                                              river_obs.minobs)

            # Recompute things set in IteratedRiverObs constructor
            river_obs.add_obs('xo', river_obs.xobs)
            river_obs.add_obs('yo', river_obs.yobs)
            river_obs.load_nodes(['xo', 'yo'])

        # Iterate through and only return reaches with no pixels in them.
        # (don't iterate and modify!)
        river_obs_list_out = []
        reach_idx_list_out = []
        ireach_list_out = []
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:

            # Filter out ghost reaches
            if reach_idx % 10 == 6:
                continue

            if len(river_obs.populated_nodes) > 0:
                river_obs_list_out.append(river_obs)
                reach_idx_list_out.append(reach_idx)
                ireach_list_out.append(ireach)

        return river_obs_list_out, reach_idx_list_out, ireach_list_out

    def assign_reaches_two_pass(
            self, scalar_max_width, minobs=10, use_width_db=False, ds=None):
        """
        Does the second pass of reach assignments using the results from the
        first pass.
        """
        river_obs_list1, reach_idx_list1, ireach_list1 = self.assign_reaches(
            scalar_max_width, minobs, use_width_db, ds)

        river_obs_list, reach_idx_list, ireach_list = self.assign_reaches(
            scalar_max_width, minobs, use_width_db, ds, second_pass=True)

        # iterate over reaches in second pass, only keep pixels that were
        # also assigned to the same reach in first pass assignments
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:

            irch1 = np.argwhere(reach_idx_list1==reach_idx)[0][0]
            river_obs1 = river_obs_list1[irch1]

            if (river_obs1.in_channel != river_obs.in_channel).any():
                # Make new in channel mask
                new_in_channel = np.logical_and(
                    river_obs.in_channel, river_obs1.in_channel)

                # subset it to the previous subset of in_channel for mask
                # of pixels to keep that were already kept.
                mask_keep = new_in_channel[river_obs.in_channel]

                # update river_obs in channel mask
                river_obs.in_channel = new_in_channel

                # Drop pixels that were double-assigned to reaches and
                # recompute things set in RiverObs constructor
                river_obs.index = river_obs.index[mask_keep]
                river_obs.d = river_obs.d[mask_keep]
                river_obs.x = river_obs.x[mask_keep]
                river_obs.y = river_obs.y[mask_keep]
                river_obs.s = river_obs.s[mask_keep]
                river_obs.n = river_obs.n[mask_keep]
                river_obs.nedited_data = len(river_obs.d)
                river_obs.populated_nodes, river_obs.obs_to_node_map = \
                    river_obs.get_obs_to_node_map(river_obs.index,
                    river_obs.minobs)

                # Recompute things set in IteratedRiverObs constructor
                river_obs.add_obs('xo', river_obs.xobs)
                river_obs.add_obs('yo', river_obs.yobs)
                river_obs.load_nodes(['xo', 'yo'])

        return river_obs_list, reach_idx_list, ireach_list

    def process_node(self,
                     reach,
                     reach_id,
                     reach_idx=None,
                     scalar_max_width=600.,
                     use_width_db=False,
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
        use_width_db : bool, default False
            Use the width data base for setting widths?
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
        numNodes = len(np.unique(self.river_obs.index))
        enough_nodes = True if numNodes - 1 > self.river_obs.k else False
        LOGGER.debug("numNodes: %d, k: %d" % (numNodes, self.river_obs.k))

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

        self.write_index_file(self.img_x[self.river_obs.in_channel],
                              self.img_y[self.river_obs.in_channel],
                              reach.node_indx[self.river_obs.index],
                              self.river_obs.d, self.river_obs.s,
                              self.river_obs.n, reach_idx,
                              segOut, self.lat[self.river_obs.in_channel],
                              self.lon[self.river_obs.in_channel],
                              self.h_noise[self.river_obs.in_channel])

        # Add the observations
        self.river_obs.add_obs('h_noise', self.h_noise)
        self.river_obs.add_obs('h_flg', (self.h_flg > 0))
        self.river_obs.add_obs('lon', self.lon)
        self.river_obs.add_obs('lat', self.lat)
        self.river_obs.add_obs('xobs', self.x)
        self.river_obs.add_obs('yobs', self.y)
        self.river_obs.add_obs('inundated_area', self.inundated_area)

        dsets_to_load = [
            'h_noise', 'h_flg', 'lon', 'lat', 'xobs', 'yobs', 'inundated_area'
        ]

        other_obs_keys = [
            'xtrack', 'sig0', 'sig0_uncert', 'water_frac', 'water_frac_uncert',
            'ifgram', 'power1', 'power2', 'phase_noise_std', 'dh_dphi',
            'dlat_dphi', 'dlon_dphi', 'num_rare_looks', 'num_med_looks',
            'false_detection_rate', 'missed_detection_rate', 'darea_dheight',
            'looks_to_efflooks', 'geoid', 'solid_earth_tide', 'load_tide_fes',
            'load_tide_got', 'pole_tide']

        for name in other_obs_keys:
            value = getattr(self, name)
            if value is not None:
                if name is 'looks_to_efflooks':
                    value = value + np.zeros(np.shape(self.lat))
                self.river_obs.add_obs(name, value)
                dsets_to_load.append(name)

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

        self.river_obs.add_obs('klass', self.klass)
        dsets_to_load.append('klass')

        self.river_obs.add_obs('pixel_area', self.pixel_area)
        dsets_to_load.append('pixel_area')

        # Adjust heights to geoid and do tide corrections 
        # (need to do before load_nodes or it is not updated in nodes)
        mask = np.logical_and(
            self.river_obs.geoid > -200, self.river_obs.geoid < 200)
        self.river_obs.h_noise[mask] -= (
            self.river_obs.geoid[mask] +
            self.river_obs.solid_earth_tide[mask] +
            self.river_obs.load_tide_fes[mask] +
            self.river_obs.pole_tide[mask])

        self.river_obs.load_nodes(dsets_to_load)
        LOGGER.debug('Observations added to nodes')

        # Get various node statistics
        nobs = np.asarray(self.river_obs.get_node_stat('count', ''))
        s_median = np.asarray(self.river_obs.get_node_stat('median', 's'))
        x_median = np.asarray(self.river_obs.get_node_stat('median', 'xobs'))
        y_median = np.asarray(self.river_obs.get_node_stat('median', 'yobs'))

        if self.xtrack is not None:
            xtrack_median = np.asarray(
                self.river_obs.get_node_stat('median', 'xtrack'))
        else:
            xtrack_median = None

        lon_median = np.asarray(
            self.river_obs.get_node_stat('sincos_median', 'lon',
                                         goodvar='h_flg'))
        lat_median = np.asarray(self.river_obs.get_node_stat('median', 'lat',
                                                             goodvar='h_flg'))
        ds = np.asarray(self.river_obs.get_node_stat('value', 'ds'))

        # number of good heights
        nobs_h = np.asarray(self.river_obs.get_node_stat('countGood', 'h_flg'))

        # heights using only "good" heights
        wse = np.asarray(self.river_obs.get_node_stat(
                'median', 'h_noise', goodvar='h_flg'))
        wse_r_u = np.asarray(
            self.river_obs.get_node_stat('std', 'h_noise', goodvar='h_flg'))
        wse_std = wse_r_u

        # The following are estimates of river width
        width_ptp = np.asarray(self.river_obs.get_node_stat('width_ptp', ''))
        width_std = np.asarray(self.river_obs.get_node_stat('width_std', ''))

        # These are area based estimates, the nominal SWOT approach
        area = np.asarray(
            self.river_obs.get_node_stat('sum', 'inundated_area'))

        rdr_sig0 = np.asarray(self.river_obs.get_node_stat(
            'median', 'sig0', goodvar='h_flg'))

        # area of pixels used to compute heights
        area_of_ht = np.asarray(self.river_obs.get_node_stat(
            'sum', 'inundated_area', goodvar='h_flg'))
        width_area = np.asarray(
            self.river_obs.get_node_stat('width_area', 'inundated_area'))

        # get the aggregated heights and widths with their corresponding
        # uncertainty estimates all in one shot
        if self.height_agg_method is 'orig' and self.area_agg_method is 'orig':
            pass
        else:
            node_aggs = self.river_obs.get_node_agg(
                height_method=self.height_agg_method,
                area_method=self.area_agg_method, goodvar='h_flg')

            latitude_u = node_aggs['lat_u']
            longitud_u = node_aggs['lon_u']
            rdr_sig0 = node_aggs['sig0']
            rdr_sig0_u = node_aggs['sig0_u']

            if self.area_agg_method is not 'orig':
                width_area = node_aggs['width_area']
                width_u = node_aggs['width_area_u']
                area = node_aggs['area']
                area_u = node_aggs['area_u']
                area_det = node_aggs['area_det']
                area_det_u = node_aggs['area_det_u']
                area_of_ht = node_aggs['area']

            if self.height_agg_method is not 'orig':
                wse = node_aggs['h']
                wse_std = node_aggs['h_std']  # height_uncert_std
                wse_r_u = node_aggs['h_u']    # height_uncert_multilook
                area_of_ht = area

        # geoid heights and tide corrections weighted by height uncertainty
        geoid_hght = np.asarray(
            self.river_obs.get_node_stat('height_weighted_mean', 'geoid',
                                         goodvar='h_flg',
                                         method=self.height_agg_method))
        solid_tide = np.asarray(
            self.river_obs.get_node_stat('height_weighted_mean',
                                         'solid_earth_tide',
                                         goodvar='h_flg',
                                         method=self.height_agg_method))
        load_tidef = np.asarray(
            self.river_obs.get_node_stat('height_weighted_mean',
                                         'load_tide_fes',
                                         goodvar='h_flg',
                                         method=self.height_agg_method))
        load_tideg = np.asarray(
            self.river_obs.get_node_stat('height_weighted_mean',
                                         'load_tide_got',
                                         goodvar='h_flg',
                                         method=self.height_agg_method))
        pole_tide = np.asarray(
            self.river_obs.get_node_stat('height_weighted_mean',
                                         'pole_tide',
                                         goodvar='h_flg',
                                         method=self.height_agg_method))

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

        # test at 5% inside of blocking width
        test_width = blocking_width * 0.90
        nanmask = np.isfinite(test_width)  # mask NaN's to avoid warning
        is_blocked = np.empty_like(nanmask)
        is_blocked[nanmask] = np.logical_or(
            np.logical_and(test_width[nanmask] < 0,
                           min_n[nanmask] < test_width[nanmask]),
            np.logical_and(test_width[nanmask] > 0,
                           max_n[nanmask] > test_width[nanmask]))
        is_blocked = np.logical_or(is_blocked, nanmask)

        lon_prior = reach.lon[self.river_obs.populated_nodes]
        lat_prior = reach.lat[self.river_obs.populated_nodes]

        dark_frac = MISSING_VALUE_FLT * np.ones(area.shape)
        dark_frac[area > 0] = 1 - area_det[area > 0] / area[area > 0]

        # Compute flow direction relative to along-track
        tangent = self.river_obs.centerline.tangent[
            self.river_obs.populated_nodes]

        fit_xx = np.array([np.ones(x_median.shape), x_median, y_median]).T
        fit = statsmodels.api.OLS(xtrack_median, fit_xx).fit()
        dxt_dx = fit.params[1]
        dxt_dy = fit.params[2]
        xt_angle = np.arctan2(dxt_dy, dxt_dx)
        at_angle = xt_angle-np.pi/2
        tangent_angle = np.arctan2(tangent[:, 1], tangent[:, 0])
        flow_dir = np.rad2deg(tangent_angle - at_angle) % 360

        prior_s = np.cumsum(reach.node_length)

        # type cast node outputs and pack it up for RiverReach constructor
        river_reach_kw_args = {
            'lat': lat_median.astype('float64'),
            'lon': lon_median.astype('float64'),
            'x': x_median.astype('float64'),
            'y': y_median.astype('float64'),
            'nobs': nobs.astype('int32'),
            's': s_median.astype('float64'),
            'ds': ds.astype('float64'),
            'w_ptp': width_ptp.astype('float64'),
            'w_std': width_std.astype('float64'),
            'w_area': width_area.astype('float64'),
            'w_db': width_db.astype('float64'),
            'area': area.astype('float64'),
            'area_u': area_u.astype('float64'),
            'area_det': area_det.astype('float64'),
            'area_det_u': area_det_u.astype('float64'),
            'area_of_ht': area_of_ht.astype('float64'),
            'wse': wse.astype('float64'),
            'wse_std': wse_std.astype('float64'),
            'wse_r_u': wse_r_u.astype('float64'),
            'nobs_h': nobs_h.astype('int32'),
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
            'x_prior': x_prior.astype('float64'),
            'y_prior': y_prior.astype('float64'),
            'lon_prior': lon_prior.astype('float64'),
            'lat_prior': lat_prior.astype('float64'),
            'p_wse': reach.wse[self.river_obs.populated_nodes],
            'p_wse_var': reach.wse_var[self.river_obs.populated_nodes],
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
            'river_name': reach.river_name[self.river_obs.populated_nodes]
        }

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
        LOGGER.debug(('number of fit points: %d' % ngood))

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

        reach_stats['width'] = np.sum(river_reach.area)/reach_stats['length']
        reach_stats['width_u'] = np.sqrt(np.sum(
            river_reach.area_u**2)) / reach_stats['length']

        reach_stats['loc_offset'] = (
            river_reach.s.mean() - self.river_obs.centerline.s.mean())

        if river_reach.xtrack is not None:
            reach_stats['xtrk_dist'] = np.median(river_reach.xtrack)

        # Along-flow distance for all PRD nodes.
        all_ss = np.cumsum(reach.node_length)

        # Make center of PRD reach the ss == 0 intercept (reference point) for
        # reach WSE.
        all_ss = all_ss - np.mean(all_ss)

        # Along-flow distance for just the observed nodes.
        ss = all_ss[self.river_obs.populated_nodes]

        hh = river_reach.wse
        ww = 1/(river_reach.wse_r_u**2)  # TO DO: validate wse_r_u here
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
                # get the optimal reconstruction (Bayes estimate)
                wse_opt, height_u, slope_u = self.optimal_reconstruct(
                    river_reach_collection,
                    river_reach, reach_id,
                    ss, hh,
                    np.sqrt(1.0 / ww),
                    min_fit_points,
                    method='Bayes',
                )
                # Use reconstruction height and slope for reach outputs
                dx = ss[0] - ss[-1]  # along-reach dist
                reach_stats['slope'] = (wse_opt[0] - wse_opt[-1]) / dx
                reach_stats['height'] = np.mean(wse_opt)
                reach_stats['slope_r_u'] = slope_u * 0.01    # m/m
                reach_stats['height_r_u'] = height_u * 0.01  # m
                reach_stats['slope_u'] = np.sqrt(
                    SLOPE_SYS_UNCERT**2 + reach_stats['slope_r_u']**2)
                reach_stats['height_u'] = np.sqrt(
                    REACH_WSE_SYS_UNCERT**2 + reach_stats['height_r_u']**2)

        else:
            # insufficient node heights for fit to reach
            reach_stats['slope'] = MISSING_VALUE_FLT
            reach_stats['slope_r_u'] = MISSING_VALUE_FLT
            reach_stats['slope_u'] = MISSING_VALUE_FLT
            reach_stats['height'] = MISSING_VALUE_FLT
            reach_stats['height_r_u'] = MISSING_VALUE_FLT
            reach_stats['height_u'] = MISSING_VALUE_FLT

        reach_stats['n_good_nod'] = mask.sum()
        reach_stats['frac_obs'] = (
            mask.sum() / len(self.river_obs.centerline.s))

        # do fit on geoid heights for reach-level outputs
        gg = river_reach.geoid_hght
        # make changes here
        #mask = self.get_reach_mask(ss, gg, ww, min_fit_points)
        if mask.sum() >= min_fit_points:
            geoid_fit = statsmodels.api.WLS(
                gg[mask], SS[mask], weights=ww[mask]).fit()

            # fit slope is meters per meter
            reach_stats['geoid_slop'] = geoid_fit.params[0]
            reach_stats['geoid_hght'] = geoid_fit.params[1]
        else:
            reach_stats['geoid_slop'] = MISSING_VALUE_FLT
            reach_stats['geoid_hght'] = MISSING_VALUE_FLT

        LOGGER.debug('Reach height/slope processing finished')

        # trap out of range / missing data
        if reach.metadata['lakeFlag'] < 0 or reach.metadata['lakeFlag'] > 255:
            uint8_flg = 255
        else:
            uint8_flg = reach.metadata['lakeFlag']
        reach_stats['lake_flag'] = uint8_flg
        reach_stats['centerline_lon'] = reach.metadata['centerline_lon']
        reach_stats['centerline_lat'] = reach.metadata['centerline_lat']

        # Compute discharge
        discharge_model_values = SWOTRiver.discharge.compute(
            reach, reach_stats['height'], reach_stats['width'],
            reach_stats['slope'])

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

        reach_stats['dark_frac'] = (
            1-np.sum(river_reach.area_det)/np.sum(river_reach.area))

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
        reach_stats['p_low_slp'] = reach.metadata['p_low_slp']
        reach_stats['river_name'] = reach.metadata['river_name']

        dsch_m_uc = reach.metadata['discharge_models']['unconstrained']
        dsch_m_c = reach.metadata['discharge_models']['constrained']

        # Avoid letting fill value of -9999 from PRD propagate into outputs
        # (these variables are just passed through from PRD to RiverTile).
        def fill_if_was_fill(value, other_fill, fill):
            return value if value != other_fill else fill
        reach_stats['dschg_msf'] = fill_if_was_fill(
            dsch_m_uc['MetroMan']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_bsf'] = fill_if_was_fill(
            dsch_m_uc['BAM']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_hsf'] = fill_if_was_fill(
            dsch_m_uc['HiVDI']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_osf'] = fill_if_was_fill(
            dsch_m_uc['MOMMA']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_ssf'] = fill_if_was_fill(
            dsch_m_uc['SADS']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_isf'] = fill_if_was_fill(
            dsch_m_uc['SIC4DVar']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)

        reach_stats['dschg_gmsf'] = fill_if_was_fill(
            dsch_m_c['MetroMan']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_gbsf'] = fill_if_was_fill(
            dsch_m_c['BAM']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_ghsf'] = fill_if_was_fill(
            dsch_m_c['HiVDI']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_gosf'] = fill_if_was_fill(
            dsch_m_c['MOMMA']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_gssf'] = fill_if_was_fill(
            dsch_m_c['SADS']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)
        reach_stats['dschg_gisf'] = fill_if_was_fill(
            dsch_m_c['SIC4DVar']['sbQ_rel'].item(), -9999, MISSING_VALUE_FLT)

        # Put in discharge_model_values into reach_stats for output
        reach_stats.update(discharge_model_values)

        river_reach.metadata = reach_stats
        return river_reach

    def flag_outliers(self, wse, flow_dist, weights, mask):
        """
        Flag wse outliers within a reach using iterative WLS slope algorithm
        wse: node wse
        flow_dist: node flow distance
        mask: overall mask that wse and flow_dist are subset using

        valid methods: ['iterative_linear', 'piecewise_linear']

        iterative_linear outlier flagging method uses these class attributes:
        - self.outlier_abs_thresh: absolute threshold, default 2[m]
        - self.outlier_rel_thresh: relative threshold, default 68% percentile
        - self.outlier_upr_thresh: upper limit, i.e. what percent of nodes are
                                   kept at the last iteration, default 80%
        - self.outlier_iter_num:   number of iterations

        piecewise_linear outlier flagging method uses these class attributes:
        - self.outlier_abs_thresh:       absolute threshold, default 1.5[m]
        - self.outlier_upr_thresh:       upper limit, i.e. what percent of
                                         nodes are kept at the last iteration,
                                         default 80%
        - self.outlier_breakpoint_min_dist: Minimum distance between
                                            breakpoints, as a proportion of the
                                            input data range, default 0.1
        - self.outlier_edge_min_dist:    Minimum distance from edge to a
                                         breakpoint, as a proportion of the
                                         input data range, default 0.1
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
        - outlier_abs_thresh
        - outlier_upr_thresh
        - outlier_iter_num

        outputs: returns the input mask array, further subset using the
                 outlier flagging method.
        """
        if self.outlier_method == 'iterative_linear':
            # initial fit OLS
            fit = statsmodels.api.OLS(wse, flow_dist).fit()
            r_slp = fit.params[0]
            r_wse = fit.params[1]
            wse_fit = r_wse + r_slp * flow_dist[:, 0]

            for i in range(self.outlier_iter_num):
                metric = abs(wse_fit - wse)
                metric_one_sigma = np.percentile(
                    metric, self.outlier_rel_thresh)
                upr_threshs = np.zeros(self.outlier_iter_num)
                upr_threshs[-1] = self.outlier_upr_thresh
                frac_keep = np.percentile(metric, upr_threshs[i])
                icond = np.logical_or(metric < self.outlier_abs_thresh,
                                      metric < metric_one_sigma)

                if sum(icond) / len(metric) * 100 <= upr_threshs[i]:
                    current_ind = metric < frac_keep
                else:
                    current_ind = icond
                if i == self.outlier_iter_num - 1 and all(np.isfinite(weights)):
                    # last iteration use WLS
                    fit = statsmodels.api.WLS(wse[current_ind],
                                              flow_dist[current_ind],
                                              weights=weights[current_ind]).fit()
                    r_slp = fit.params[0]
                    r_wse = fit.params[1]
                else:
                    fit = statsmodels.api.OLS(wse[current_ind],
                                              flow_dist[current_ind]).fit()
                    r_slp = fit.params[0]
                    r_wse = fit.params[1]

                wse_fit = r_wse + r_slp * flow_dist[:, 0]
            mask[mask] = current_ind
            
        elif self.outlier_method == 'piecewise_linear':
            current_ind = None
            # at least to have 10 node (2km) to run piecewise linear
            if len(wse) > 10:
                break_num = int(np.floor(len(wse)/10)-1)
                while current_ind is None and break_num >=1:
                    current_ind = self.piecewise_linear(flow_dist[:, 0],
                                                        wse,
                                                        break_num)
                    break_num -= 1

            if current_ind is None:
                fit = statsmodels.api.OLS(wse, flow_dist).fit()
                wse_fit = fit.params[1] + fit.params[0] * flow_dist[:, 0]
                metric = abs(wse_fit - wse)
                metric_upr_thresh = np.percentile(
                    metric, self.outlier_upr_thresh)
                if self.outlier_abs_thresh < metric_upr_thresh:
                    current_ind = metric < metric_upr_thresh
                else: current_ind = metric < self.outlier_abs_thresh
            mask[mask] = current_ind

        else:
            raise NotImplementedError(
                'Outlier flagging method %s is not supported!'%
                self.outlier_method)

        return mask

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

    def write_index_file(self, img_x, img_y, node_index, dst, along_reach,
                         cross_reach, reach_index, seg_lbl, lat, lon,
                         height):
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
        wse, wse_r_u, this_len, first_node, ss, mask, p_wse = \
            self.get_multi_reach_height(
            river_reach_collection, river_reach, ireach)

        if np.sum(river_reach.mask) < min_fit_points:
            enhanced_slope = MISSING_VALUE_FLT
        else:
            # handle indexing for masked reaches
            this_reach_edges = np.ma.flatnotmasked_edges(
                np.ma.masked_array(river_reach.wse, mask=~river_reach.mask))
            this_reach_len = river_reach.node_ss[this_reach_edges[1]] \
                             - river_reach.node_ss[this_reach_edges[0]]
            if first_node == 0:
                # No upstream reach. Do not adjust first node location
                first_node_masked = 0
            else:
                # Adjust first node location based on upstream masked nodes
                first_node_masked = first_node - np.sum(
                    ~mask[:first_node + this_reach_edges[0]])
            last_node_masked = first_node_masked + np.sum(river_reach.mask) - 1

            # window size and sigma for Gaussian averaging
            window_size = np.min([max_window_size, this_reach_len])
            sigma = np.max([min_sigma, window_size / window_size_sigma_ratio])
            # smooth h_n_ave, and get slope
            slope = np.polyfit(ss[mask], wse[mask], 1)[0]
            heights_detrend = wse[mask] - slope * ss[mask]
            heights_smooth = self.gaussian_averaging(
                ss[mask], heights_detrend, window_size, sigma)
            heights_smooth = heights_smooth + slope*(ss[mask] - ss[mask][0])
            enhanced_slope = (heights_smooth[last_node_masked] -
                              heights_smooth[first_node_masked])/this_reach_len

        if np.isnan(enhanced_slope):
            enhanced_slope = MISSING_VALUE_FLT
        return enhanced_slope

    @staticmethod
    def gaussian_averaging(ss, wse, window_size, sigma):
        """
        Gaussian smoothing of heights using distances
        ss:   along-flow distance
        wse:     water heights
        window_size: size of data window to use for averaging
        sigma:       STD of Gaussian used for averaging

        outputs:
        smooth_heights : smoothed elevations
        """
        smooth_heights = np.zeros(wse.shape)
        for ii, this_distance in enumerate(ss):

            # get data window
            mask = np.logical_and(
                np.abs(this_distance-ss) <= window_size / 2,
                ~np.isnan(wse))

            weights = scipy.stats.norm.pdf(
                this_distance, ss[mask], sigma)

            smooth_heights[ii] = (
                np.multiply(weights, wse[mask]).sum() /
                weights.sum())

        return smooth_heights

    def get_reach_mask(self, ss, hh, ww, min_fit_points):
        """
        Mask each node in an input reach based on whether or not it has a valid
        wse. Then mask for node-level wse outliers, if self.outlier_method is
        set. Returns a numpy array where good nodes are 1 and bad nodes are 0.
        :param ss: along-reach distance for each node
        :param hh: height of each node
        :param ww: weights for each node height, to use in outlier flagging
        :param min_fit_points: minimum number of points needed to flag outliers
        :return: reach mask where good nodes are 1 and bad nodes are 0
        """
        mask = np.logical_and(hh > -500, hh < 8000)
        if self.outlier_method is not None and mask.sum() > min_fit_points:
            SS = np.c_[ss, np.ones(len(ss), dtype=ss.dtype)]
            mask = self.flag_outliers(hh[mask], SS[mask], ww[mask], mask)

        elif self.outlier_method is not None and mask.sum() <= min_fit_points:
            mask = np.zeros(len(hh), dtype=bool)
        return mask

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
        other_ids = [
            item.reach_indx[0] for item in river_reach_collection
        ]

        # get up/dn id from prior db
        prior_s = river_reach.prior_node_ss

        # skip (disconnected lake, dam) reaches
        skip_types = [2, 4]

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

                if side == -1:
                    # side is downstream of current reach
                    dx = try_prd_rch.x[-1] - prd_rch[0].x[0]
                    dy = try_prd_rch.y[-1] - prd_rch[0].y[0]
                else:
                    # side is upstream of current reach
                    dx = try_prd_rch.x[0] - prd_rch[0].x[-1]
                    dy = try_prd_rch.y[0] - prd_rch[0].y[-1]

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
        mask = np.array([], dtype=bool)
        # if upstream PRD reach is usable
        if prd_is_good[1]:
            ss = np.concatenate([adj_rch[1].node_ss, ss])
            p_wse = np.concatenate([adj_rch[1].p_wse, p_wse])
            wse = np.concatenate([adj_rch[1].wse, wse])
            wse_r_u = np.concatenate([adj_rch[1].wse_r_u, wse_r_u])
            mask = np.concatenate([adj_rch[1].mask, mask])

        ss = np.concatenate([river_reach.node_ss, ss+prior_s[-1]])
        p_wse = np.concatenate([river_reach.p_wse, p_wse])
        wse = np.concatenate([river_reach.wse, wse])
        wse_r_u = np.concatenate([river_reach.wse_r_u, wse_r_u])
        mask = np.concatenate([river_reach.mask, mask])
        this_len = len(river_reach.wse)

        # if downstream PRD reach is usable
        if prd_is_good[-1]:
            downstream_prior_s = adj_rch[-1].prior_node_ss
            first_node = first_node + len(adj_rch[-1].wse)
            ss = np.concatenate([adj_rch[-1].node_ss,
                                 ss+downstream_prior_s[-1]])
            p_wse = np.concatenate([adj_rch[-1].p_wse, p_wse])
            wse = np.concatenate([adj_rch[-1].wse, wse])
            wse_r_u = np.concatenate([adj_rch[-1].wse_r_u, wse_r_u])
            mask = np.concatenate([adj_rch[-1].mask, mask])

        return wse, wse_r_u, this_len, first_node, ss, mask, p_wse


    def optimal_reconstruct(
            self,
            river_reach_collection,
            river_reach,
            ireach,
            ss,
            wse,
            wse_r_u,
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
        wse        : Measured Node wse (with masked values for missing nodes)
        wse_r_u    : Node-wise random wse uncertainty
        ss         : Node-level distance from prior database
        mask       : True if node-level height is good, False where it is bad
        Options:
        method = Bayes
           Bayes        : Bayes estimate of wse given prior and
                          prior_cov method.
        full_noise_cov : True, or False
                         This option sets the noise covariance and the sampling
                         operator to assume that all nodes are observed (but
                         weighted appropriately for the unobserved nodes).
        prior_wse : The mean wse of the prior we want to impose (defaults to
                    the weighted linear fit to reach if there is no input)
        prior_cov_method : independent, exponential
                           This option controls how spatial structure is
                           imposed by the prior.
        """

        if self.use_multiple_reaches:
            wse, wse_r_u, this_len, first_node, ss, mask, prior_wse = \
                self.get_multi_reach_height(
                    river_reach_collection, river_reach, ireach)
            # get the multi-reach mask
            ss = ss - np.mean(ss)
            end_slice = first_node + this_len
        else:
            first_node = 0
            end_slice = len(ss)

        # define vectors b and c for uncertainty estimates later
        this_reach_mask_b = np.zeros_like(ss)
        this_reach_mask_b[first_node:end_slice] = 1
        first_and_last_node_c = np.zeros_like(ss)
        first_and_last_node_c[first_node] = -1
        first_and_last_node_c[end_slice-1] = 1

        # create a wse prior if flagged
        if self.prior_wse_method == 'fit':
            prior_wse = None
        elif self.prior_wse_method == 'prd':
            if self.use_multiple_reaches:
                # prior_wse already set above
                pass
            else:
                # get prior wse from PRD
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
        height_u = this_reach_mask_b @ A_inv @ np.atleast_2d(
            this_reach_mask_b).T
        slope_u = first_and_last_node_c @ A_inv @ np.atleast_2d(
            first_and_last_node_c).T
        if self.use_multiple_reaches:
            wse_out = wse_out[first_node:end_slice]
        return wse_out, height_u[0], slope_u[0]

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

    def piecewise_linear(self, x, y, n_breakpoints):
        """
        Piecewise linear outlier detector

        Parameters
        ----------
        x : Node-level flow distance
        y : Measured Node wse
        n_breakpoints : Number of breakpoints
        For more information of the piecewise linear fit:
        https://github.com/chasmani/piecewise-regression

        Outputs
        ----------
        returns the input mask array
        """
        pw_fit = piecewise_regression.Fit(
            x, y, n_breakpoints=n_breakpoints,
            min_distance_between_breakpoints=self.outlier_breakpoint_min_dist,
            min_distance_to_edge=self.outlier_edge_min_dist,
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
            metric_upr_thresh = np.percentile(metric, self.outlier_upr_thresh)
            if self.outlier_abs_thresh < metric_upr_thresh:
                current_ind = metric < metric_upr_thresh
            else:
                current_ind = metric < self.outlier_abs_thresh
        else:
            current_ind = None
        return current_ind

