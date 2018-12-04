"""
Given a SWOTL2 file, fit all of the reaches observed and output results.
"""

from __future__ import absolute_import, division, print_function

import os
import scipy.ndimage
import numpy as np
import pandas as pd
import netCDF4 as nc
import collections
import scipy.stats

from .SWOTL2 import SWOTL2
from RiverObs import ReachExtractor
from RiverObs import WidthDataBase
from RiverObs import IteratedRiverObs
from RiverObs import FitRiver
from RiverObs import RiverNode
from RiverObs import RiverReach
from Centerline.Centerline import CenterLineException

class SWOTRiverEstimator(SWOTL2):
    """
    Given a SWOTL2 file, fit all of the reaches observed and output results.

    This class is derived from the SWOTL2 class.

    This class contains four basic components:
    1. It is a subclass of the SWOTL2, and therefor is a LatLonRegion and
    also contains SWOT data.
    2. For each reach of interest, it keeps a dictionary RiverObs instances
    for each reach.
    3. It keeps a dictionary of RiverReach instances, containing the results
    of height/slope fitting and width estimation.
    4. A dictionary of FitRiver instances, containing the full fit diagnostics.

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
        thruth
    class_list : list, default [2,3,4,5]
        a list of the class labels for what is considered good water data.
        This should be piched as any classes that might contain water, even
        if partially. If truth data are desired (class_kwd =
        'no_layover_classification'), the this should be set to [1].
        The default has interior water pixels(4), border water pixels (3),
        and border land pixels (2). This should be used with
        inundation_fraction turned on.
    fractional_inundation_kwd : str, default 'continuous_classification'
        Netcdf keyword containing the inundation fraction. If None, the no
        inundation fraction is used. The inundation fraction is, in theory,
        a number between 0 and 1 estimating the fraction of the pixel covered
        with water. In practice, because of noise, it can be outside the
        bounds and even be negative! It will produced an ensemble
        unbiased estimate if the class mean cross sections are known.
    use_fractional_inundation : bool list, default [True, True, False, False]
        For which classes should the inundation fraction be used instead of a
        binary classification? the default is to assume that interior pixels
        are 100% water, but to use both land and water edge pixels partially
        by using the fractional inundation kwd.
    use_segmentation :  bool list, default [False, True, True, True]
        Selects which classes to assume as water for segmentation purposes
    use_heights : bool list, default [False, False, True, False]
        Selects which classes to use for estimating heights
    min_points : int
        If the number of good points is less than this, raise an exception.
    height_kwd : str, default 'height'
        These are the heights to exract from the water file.
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
    store_reaches : bool, default True
        If True, store each RiverRiver instance in a dictionary.
    store_fits : bool, default True
        If True, store each fit result in a dictionary.
    verbose : bool, default False
        If True, print to sdtout while processing.
    use_segmentation : bool list, default [False, True, True, True]
        Defines which classes should the assumed as water for segmatation
        algorithm to label disjoint features
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
                 class_list=[2, 3, 4, 5],
                 class_kwd='classification',
                 rngidx_kwd='range_index',
                 aziidx_kwd='azimuth_index',
                 fractional_inundation_kwd='continuous_classification',
                 use_fractional_inundation=[True, True, False, False],
                 use_segmentation=[False, True, True, True],
                 use_heights=[False, False, True, False],
                 min_points=100,
                 height_kwd='height',
                 trim_ends=False,
                 store_obs=True,
                 store_reaches=True,
                 store_fits=True,
                 verbose=True,
                 xtrack_kwd='no_layover_cross_track',
                 sig0_kwd='sig0',
                 proj='laea',
                 x_0=0,
                 y_0=0,
                 lat_0=None,
                 lon_0=None,
                 ellps='WGS84',
                 output_file=None,
                 subsample_factor=1,
                 **proj_kwds):

        self.trim_ends = trim_ends
        self.store_obs = store_obs
        self.store_reaches = store_reaches
        self.store_fits = store_fits
        self.verbose = verbose
        self.input_file = os.path.split(swotL2_file)[-1]
        self.output_file = output_file  # index file
        self.subsample_factor = subsample_factor

        # Classification inputs
        self.class_kwd = class_kwd
        self.class_list = class_list
        self.fractional_inundation_kwd = fractional_inundation_kwd
        self.use_fractional_inundation = use_fractional_inundation

        self.use_segmentation = use_segmentation
        self.use_heights = use_heights

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
            verbose=verbose,
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

        try:
            self.xtrack = self.get(xtrack_kwd)
        except KeyError:
            self.xtrack = None

        try:
            self.sig0 = self.get(sig0_kwd)
        except KeyError:
            self.sig0 = None

        good = ~mask
        self.lat = self.lat[good]
        self.lon = self.lon[good]
        self.x = self.x[good]
        self.y = self.y[good]
        self.klass = self.klass[good]
        self.h_noise = self.h_noise[good]
        if self.xtrack is not None:
            self.xtrack = self.xtrack[good]
        if self.sig0 is not None:
            self.sig0 = self.sig0[good]
        self.img_x = self.img_x[good]  # range or x index
        self.img_y = self.img_y[good]  # azimuth or y index

        # set the pixcvec geolocations to the pixel cloud values
        # TODO: update these with height constrained geolocation
        self.lat_vec = self.lat
        self.lon_vec = self.lon
        self.height_vec = self.h_noise

        # Try to read the pixel area from the L2 file, or compute it
        # from look angle and azimuth spacing, or from azimuth spacing
        # and ground spacing

        try:
            # hopefully already there
            self.pixel_area = self.get('pixel_area')

        except KeyError:
            try:
                # try compute with look angle
                look_angle = self.get('no_layover_look_angle')[good]
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
        if (self.subsample_factor > 1):
            self.pixel_area = self.pixel_area * self.subsample_factor
        #
        if fractional_inundation_kwd is None:  # all water pixels are inundated
            self.fractional_inundation = None
            self.inundated_area = self.pixel_area

        else:
            self.fractional_inundation = self.get(fractional_inundation_kwd)
            self.inundated_area = self.pixel_area
            for i, k in enumerate(class_list):
                if use_fractional_inundation[i]:
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
            self.segment_water_class()
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

        if self.verbose:
            print('Data loaded')

        # Initialize the list of observations and reaches
        self.river_obs_collection = collections.OrderedDict()
        self.river_reach_collection = collections.OrderedDict()
        self.fit_collection = collections.OrderedDict()

    def segment_water_class(self):
        """
        do image segmentation algorithm on the water class to label
        unconnected features
        """
        maxX = np.max(self.img_x)
        maxY = np.max(self.img_y)
        cls_img = np.zeros((maxY + 1, maxX + 1))
        cls_img[self.img_y, self.img_x] = self.isWater

        # Do some regularization with morphological operations? (TODO)
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

    def get_reaches(self, shape_file_root, clip=True, clip_buffer=0.1):
        """Get all of the reaches using a ReachExtractor."""
        self.clip = clip
        self.clip_buffer = clip_buffer

        self.reaches = ReachExtractor(
            shape_file_root, self, clip=clip, clip_buffer=clip_buffer)

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
                        fit_types=['OLS', 'WLS', 'RLM'],
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
        fit_types : list
            A list of fit types to perform. Can contain 'OLS','WLS','RLM'.
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

        Returns a list containg a RiverReach instance for each reach in the
        bounding box.
        """
        # assign the reaches
        river_obs_list, reach_idx_list, ireach_list = self.assign_reaches(
            scalar_max_width, minobs, use_width_db, ds)

        river_reach_collection = []
        reach_zips = zip(river_obs_list, reach_idx_list, ireach_list)
        for river_obs, reach_idx, ireach in reach_zips:

            if use_width_db:
                max_width = self.get_max_width_from_db(reach_idx)
                if self.verbose:
                    print('max_width read')

            else:
                max_width = None

            if max_width is None:
                try:
                    # probably should scale this to look some fraction
                    # farther than the database width
                    max_width = (
                        self.reaches[i_reach].metadata['Wmean'] * np.ones(
                            np.shape(self.reaches[i_reach].x)))  # *2.0

                except:
                    max_width = scalar_max_width

            # Ugly way process_reach/process_node uses the data
            self.river_obs = river_obs

            river_reach = self.process_node(
                self.reaches[ireach],
                ireach,
                reach_idx,
                scalar_max_width=scalar_max_width,
                minobs=minobs,
                use_width_db=use_width_db,
                max_width=max_width,
                ds=ds,
                refine_centerline=refine_centerline,
                smooth=smooth,
                alpha=alpha,
                max_iter=max_iter)

            river_reach_collection.append(river_reach)
            if self.store_reaches:
                self.river_reach_collection[ireach] = river_reach

            if self.verbose:
                print('reach pocessed')

        # calculate reach enhanced slope, and add to river_reach_collection
        enhanced_slopes = self.compute_enhanced_slopes(
            river_reach_collection, max_window_size=max_window_size,
            min_sigma=min_sigma,
            window_size_sigma_ratio=window_size_sigma_ratio,
            enhanced=enhanced)

        out_river_reach_collection = []
        # Now iterate over reaches again and do reach average computations
        reach_zips = zip(
            river_reach_collection, river_obs_list, reach_idx_list,
            ireach_list, enhanced_slopes)
        for river_reach, river_obs, reach_idx, ireach, enhanced_slope in\
            reach_zips:

            # Ugly way process_reach/process_node uses the data
            self.river_obs = river_obs

            out_river_reach = self.process_reach(river_reach,
                                                 self.reaches[ireach],
                                                 ireach,
                                                 reach_idx,
                                                 min_fit_points=min_fit_points,
                                                 fit_types=fit_types)

            if out_river_reach is not None:
                # add enhanced slope to river reach outputs
                out_river_reach.metadata['slp_enhncd'] = np.float32(
                    enhanced_slope)

                out_river_reach_collection.append(out_river_reach)

        return out_river_reach_collection

    def assign_reaches(self,
                       scalar_max_width,
                       minobs=10,
                       use_width_db=False,
                       ds=None):
        """
        Assigns pixels to nodes for every reach.
        """
        # Iterate over reaches, assign pixels to nodes
        river_obs_list = []
        reach_idx_list = []
        ireach_list = []
        for i_reach, reach_idx in enumerate(self.reaches.reach_idx):

            if len(self.reaches[i_reach].x) <= 3:
                print("reach does not have enough points",
                      len(self.reaches[i_reach].x))
                continue

            if self.verbose:
                print('Reach %d/%d Reach index: %d' %
                      (i_reach + 1, self.reaches.nreaches, reach_idx))

            try:
                river_obs = IteratedRiverObs(
                    self.reaches[i_reach],
                    self.x,
                    self.y,
                    ds=ds,
                    seg_label=self.seg_label,
                    max_width=scalar_max_width,
                    minobs=minobs,
                    verbose=self.verbose)

            except CenterLineException as e:
                print("CenterLineException: ", e)
                continue

            if len(river_obs.x) == 0:
                print('No observations mapped to nodes in this reach')
                continue

            river_obs_list.append(river_obs)
            reach_idx_list.append(reach_idx)
            ireach_list.append(i_reach)

        # Ensure unique and optimal assignments of pixels to reach.
        min_dist = 9999999 * np.ones(self.x.shape)
        reach_ind = -1 * np.ones(self.x.shape, dtype=int)
        cnts_assigned = np.zeros(self.x.shape, dtype=int)
        for ii, river_obs in enumerate(river_obs_list):

            # Get current reach assingment and min distance to node for all
            # pixels assigned to this reach.
            these_reach_inds = reach_ind[river_obs.in_channel]
            these_min_dists = min_dist[river_obs.in_channel]

            # Figure out which ones are better than current assignment
            mask = river_obs.d < these_min_dists

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
                river_obs.get_obs_to_node_map(river_obs.index, river_obs.minobs)

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
            if len(river_obs.s) > 0:
                river_obs_list_out.append(river_obs)
                reach_idx_list_out.append(reach_idx)
                ireach_list_out.append(ireach)

        return river_obs_list_out, reach_idx_list_out, ireach_list_out

    def process_node(self,
                     reach,
                     reach_id,
                     reach_idx=None,
                     scalar_max_width=600.,
                     minobs=10,
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
        minobs : int, default 10
            Minimum number of observations for valid node.
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
        # get the number of node inthe reach and only refine if there are
        # enough to do spline
        numNodes = len(np.unique(self.river_obs.index))
        enough_nodes = True if numNodes - 1 > self.river_obs.k else False
        if self.verbose:
            print("numNodes,k:", numNodes, self.river_obs.k)

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
            if self.verbose:
                print('centerline refined')

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

        except:
            segOut = None

        self.write_index_file(self.img_x[self.river_obs.in_channel],
                              self.img_y[self.river_obs.in_channel],
                              self.river_obs.index, self.river_obs.d,
                              self.river_obs.s, self.river_obs.n, reach_idx,
                              segOut, self.h_flg[self.river_obs.in_channel],
                              self.lat_vec[self.river_obs.in_channel],
                              self.lon_vec[self.river_obs.in_channel],
                              self.height_vec[self.river_obs.in_channel])

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

        node_index = np.arange(len(y_prior))
        node_index = node_index[self.river_obs.populated_nodes]
        y_prior = y_prior[self.river_obs.populated_nodes]
        reach_index = np.ones(len(node_index)) * (reach_idx)

        # Add the observations
        self.river_obs.add_obs('h_noise', self.h_noise)
        self.river_obs.add_obs('h_flg', (self.h_flg > 0))
        self.river_obs.add_obs('lon', self.lon)
        self.river_obs.add_obs('lat', self.lat)
        self.river_obs.add_obs('xobs', self.x)
        self.river_obs.add_obs('yobs', self.y)
        if self.xtrack is not None:
            self.river_obs.add_obs('xtrack', self.xtrack)
        if self.sig0 is not None:
            self.river_obs.add_obs('sig0', self.sig0)
        self.river_obs.add_obs('inundated_area', self.inundated_area)

        dsets_to_load = [
            'h_noise', 'h_flg', 'lon', 'lat', 'xobs', 'yobs', 'inundated_area'
        ]

        if self.xtrack is not None:
            dsets_to_load.append('xtrack')

        if self.sig0 is not None:
            dsets_to_load.append('sig0')

        self.river_obs.load_nodes(dsets_to_load)

        if self.verbose:
            print('Observations added to nodes')

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
            self.river_obs.get_node_stat('sincos_median', 'lon'))
        lat_median = np.asarray(self.river_obs.get_node_stat('median', 'lat'))
        ds = np.asarray(self.river_obs.get_node_stat('value', 'ds'))

        # number of good heights
        nobs_h = np.asarray(self.river_obs.get_node_stat('countGood', 'h_flg'))
        # heights using same pixels as widths
        h_noise_ave0 = np.asarray(
            self.river_obs.get_node_stat('median', 'h_noise'))
        h_noise_std0 = np.asarray(
            self.river_obs.get_node_stat('std', 'h_noise'))
        # heights using only "good" heights
        h_noise_ave = np.asarray(
            self.river_obs.get_node_stat(
                'median', 'h_noise', good_flag='h_flg'))
        h_noise_std = np.asarray(
            self.river_obs.get_node_stat('std', 'h_noise', good_flag='h_flg'))

        # The following are estimates of river width
        width_ptp = np.asarray(self.river_obs.get_node_stat('width_ptp', ''))
        width_std = np.asarray(self.river_obs.get_node_stat('width_std', ''))

        # These are area based estimates, the nominal SWOT approach
        area = np.asarray(
            self.river_obs.get_node_stat('sum', 'inundated_area'))

        rdr_sig0 = np.asarray(
            self.river_obs.get_node_stat('median', 'sig0', good_flag='h_flg'))

        # area of pixels used to compute heights
        area_of_ht = np.asarray(
            self.river_obs.get_node_stat('sum', 'inundated_area',
                                         good_flag='h_flg'))

        width_area = np.asarray(
            self.river_obs.get_node_stat('width_area', 'inundated_area'))

        # These are the values from the width database
        width_db = np.ones(
            self.river_obs.n_nodes,
            dtype=np.float64) * self.river_obs.missing_value

        try:
            windex = self.river_obs.centerline_obs['max_width'].populated_nodes
            width_db[windex] = self.river_obs.centerline_obs['max_width'].v
            width_db = width_db[self.river_obs.populated_nodes]

        except:
            width_db = np.ones(
                self.river_obs.n_nodes,
                dtype=np.float64) * self.river_obs.missing_value
            width_db = width_db[self.river_obs.populated_nodes]

        # type cast node outputs and pack it up for RiverReach constructor
        river_reach_kw_args = {
            'lat': lat_median.astype('float64'),
            'lon': lon_median.astype('float64'),
            'x': x_median.astype('float64'),
            'y': y_median.astype('float64'),
            'nobs': nobs.astype('int32'),
            's': s_median.astype('float64'),
            'ds': ds.astype('float64'),
            'w_ptp': width_ptp.astype('float32'),
            'w_std': width_std.astype('float32'),
            'w_area': width_area.astype('float32'),
            'w_db': width_db.astype('float32'),
            'area': area.astype('float32'),
            'area_of_ht': area_of_ht.astype('float32'),
            'h_n_ave': h_noise_ave.astype('float32'),
            'h_n_std': h_noise_std.astype('float32'),
            'h_a_ave': h_noise_ave0.astype('float32'),
            'h_a_std': h_noise_std0.astype('float32'),
            'nobs_h': nobs_h.astype('int32'),
            'x_prior': x_prior.astype('float64'),
            'y_prior': y_prior.astype('float64'),
            'node_indx': node_index.astype('int32'),
            'reach_indx': reach_index.astype('int32'),
            'rdr_sig0': rdr_sig0.astype('float32'),
        }

        if xtrack_median is not None:
            river_reach_kw_args['xtrack'] = xtrack_median
        else:
            river_reach_kw_args['xtrack'] = None

        river_reach = RiverReach(**river_reach_kw_args)

        # Store, if desired
        if self.store_obs:
            self.river_obs_collection[reach_idx] = self.river_obs

        return river_reach

    def process_reach(self,
                      river_reach,
                      reach,
                      reach_id,
                      reach_idx=None,
                      min_fit_points=3,
                      fit_types=['OLS', 'WLS', 'RLM']):
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
        min_fit_points : int, default 3
            Minimum number of populated nodes required for height/slope fit
        fit_types : list
            A list of fit types to perform. Can contain 'OLS','WLS','RLM'.

        Returns
        -------
        Nothing

        Modifies river_reach with the reach quantities (river_reach.metadata)
        """
        # Initialize the fitter
        self.fitter = FitRiver(self.river_obs)

        # Check to see if there are sufficient number of points for fit
        ngood = len(river_reach.s)
        if self.verbose:
            print(('number of fit points: %d' % ngood))

        if ngood < min_fit_points:
            if self.verbose:
                print('not enough good points for fit')

            nresults = None
            return None

        else:
            # Get the start and end
            smin = river_reach.s.min()
            smax = river_reach.s.max()
            # Do the fitting for this reach
            nresults = self.estimate_height_slope(
                smin, smax, fit_types=fit_types, mean_stat='median')
            if self.verbose:
                print('Estimation finished')

        if self.store_fits:
            self.fit_collection[reach_id, 'noise'] = nresults

        # Get the reach statistics for this subreach
        reach_stats = self.get_reach_stats(
            reach_id, reach_idx, river_reach, nresults)

        # trap out of range / missing data
        if reach.metadata['lakeFlag'] < 0 or reach.metadata['lakeFlag'] > 255:
            uint8_flg = 255
        else:
            uint8_flg = reach.metadata['lakeFlag']
        reach_stats['lake_flag'] = uint8_flg

        river_reach.metadata = reach_stats
        return river_reach

    def estimate_height_slope(
        self, smin, smax, fit_types=['OLS', 'WLS', 'RLM'], mean_stat='median'):
        """Get fit results for this subreach."""
        if type(fit_types) == str:
            fit_types = [fit_types]
        nresults = collections.OrderedDict()
        load_inputs = True
        for fit_type in fit_types:
            try:
                nresults[fit_type] = self.fitter.fit_linear(
                    smin, smax, 'h_noise', fit=fit_type, mean_stat=mean_stat,
                    load_inputs=load_inputs)
                load_inputs = False
            except (ZeroDivisionError, ValueError):
                pass
        return nresults

    def get_reach_stats(self, reach_id, reach_idx, river_reach, nresults):
        """Get statistics for a given reach."""

        ds = np.divide(river_reach.area, river_reach.w_area)

        reach_stats = collections.OrderedDict()
        reach_stats['reach_id'] = reach_id
        reach_stats['reach_idx'] = reach_idx
        reach_stats['lon_min'] = np.min(river_reach.lon)
        reach_stats['lon_max'] = np.max(river_reach.lon)
        reach_stats['lat_min'] = np.min(river_reach.lat)
        reach_stats['lat_max'] = np.max(river_reach.lat)
        reach_stats['area'] = np.sum(river_reach.area)
        reach_stats['area_of_ht'] = np.sum(river_reach.area_of_ht)
        reach_stats['length'] = np.sum(ds)
        reach_stats['smin'] = np.min(river_reach.s)
        reach_stats['smax'] = np.max(river_reach.s)
        reach_stats['save'] = np.median(river_reach.s)
        reach_stats['w_ptp_ave'] = np.median(river_reach.w_ptp)
        reach_stats['w_ptp_min'] = np.min(river_reach.w_ptp)
        reach_stats['w_ptp_max'] = np.max(river_reach.w_ptp)
        reach_stats['w_std_ave'] = np.median(river_reach.w_std)
        reach_stats['w_std_min'] = np.min(river_reach.w_std)
        reach_stats['w_std_max'] = np.max(river_reach.w_std)
        reach_stats['w_area_ave'] = np.sum(
            river_reach.area) / reach_stats['length']

        reach_stats['n_good_nod'] = len(river_reach.s)
        reach_stats['w_area_min'] = np.min(river_reach.w_area)
        reach_stats['w_area_max'] = np.max(river_reach.w_area)
        reach_stats['frac_obs'] = (
            len(river_reach.s) / len(self.river_obs.centerline.s))

        if river_reach.xtrack is not None:
            reach_stats['xtrck_ave'] = np.median(river_reach.xtrack)
            reach_stats['xtrck_min'] = np.min(river_reach.xtrack)
            reach_stats['xtrck_max'] = np.max(river_reach.xtrack)

        # fitting results for h_true
        if nresults is not None:
            for fit_type in ['OLS', 'WLS', 'RLM']:

                tag = {'OLS': 'o', 'WLS': 'w', 'RLM': 'r'}[fit_type]
                try:
                    this_result = nresults[fit_type]

                except KeyError:
                    reach_stats['h_n%s' % tag] = self.river_obs.missing_value
                    reach_stats['slp_n%s' % tag] = self.river_obs.missing_value
                    reach_stats['n%s_rsqrd' %
                                tag] = self.river_obs.missing_value
                    reach_stats['n%s_mse' % tag] = self.river_obs.missing_value
                    continue

                reach_stats['h_n%s' % tag] = this_result.params[1]
                reach_stats['slp_n%s' % tag] = this_result.params[0]

                try:
                    reach_stats['n%s_rsqrd' % tag] = this_result.rsquared
                    reach_stats['n%s_mse' % tag] = this_result.mse_resid

                except AttributeError:
                    reach_stats['n%s_rsqrd' %
                                tag] = self.river_obs.missing_value
                    reach_stats['n%s_mse' % tag] = self.river_obs.missing_value
                    pass

        return reach_stats

    def create_index_file(self):
        """
        Create an empty file with appropriate variable to append to
        """
        # if the filename does not exist create it
        with nc.Dataset(self.output_file, 'w') as ofp:
            ofp.createDimension('record', None)
            ofp.createVariable('range_index', 'i4', 'record', fill_value=-128)
            ofp.createVariable(
                'azimuth_index', 'i4', 'record', fill_value=-128)
            ofp.createVariable('node_index', 'i4', 'record', fill_value=-128)
            ofp.createVariable('reach_index', 'i4', 'record', fill_value=-128)
            ofp.createVariable('river_tag', 'i4', 'record', fill_value=-128)
            ofp.createVariable(
                'segmentation_label', 'i4', 'record', fill_value=-128)
            ofp.createVariable(
                'good_height_flag', 'i1', 'record', fill_value=-1)
            ofp.createVariable(
                'distance_to_node', 'f4', 'record', fill_value=-128)
            ofp.createVariable(
                'along_reach', 'f4', 'record', fill_value=-9990000000)
            ofp.createVariable(
                'cross_reach', 'f4', 'record', fill_value=-9990000000)
            ofp.createVariable(
                'latitude_vectorproc', 'f8', 'record', fill_value=-9990000000)
            ofp.createVariable(
                'longitude_vectorproc', 'f8', 'record', fill_value=-9990000000)
            ofp.createVariable(
                'height_vectorproc', 'f8', 'record', fill_value=-9990000000)

            # copy attributes from pixel cloud product
            for att_name in self.nc.__dict__:
                setattr(ofp, att_name, getattr(self.nc, att_name))

    def write_index_file(self, img_x, img_y, node_index, dst, along_reach,
                         cross_reach, reach_index, seg_lbl, h_flg, lat, lon,
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
            ofp.variables['node_index'][curr_len:new_len] = node_index
            ofp.variables['reach_index'][curr_len:new_len] = reach_index
            # set river_tag to reach_index + 1 for now (0 assumes not a reach)
            #  (TODO: figure out what this should be)
            ofp.variables['river_tag'][curr_len:new_len] = reach_index + 1
            ofp.variables['segmentation_label'][curr_len:new_len] = seg_lbl
            ofp.variables['good_height_flag'][curr_len:new_len] = h_flg
            ofp.variables['distance_to_node'][curr_len:new_len] = dst
            ofp.variables['along_reach'][curr_len:new_len] = along_reach
            ofp.variables['cross_reach'][curr_len:new_len] = cross_reach
            # for improved geolocation
            ofp.variables['latitude_vectorproc'][curr_len:new_len] = lat
            ofp.variables['longitude_vectorproc'][curr_len:new_len] = lon
            ofp.variables['height_vectorproc'][curr_len:new_len] = height
        return

    def compute_enhanced_slopes(
        self, river_reach_collection, max_window_size, min_sigma,
        window_size_sigma_ratio, enhanced):
        """
        This function calculate enhanced reach slope from smoothed
        node height using Gaussian moving average.
        For more information, pleasec see Dr. Renato Frasson's paper:
        https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017WR020887

        inputs:
        river_reach_collection: collection of river_reach instance
        max_window_size: the max of Gaussian window, default is 10km
        min_sigma : min sigma for gaussian averaging, default is 1km
        window_size_sigma_ratio : default is 5

        output:
        enhanced_slopes: enhanced reach slopes
        """
        # get list of reach index
        n_reach = len(river_reach_collection)
        ind = [reach.reach_indx[0] for reach in river_reach_collection]

        enhanced_slopes = []
        for river_reach in river_reach_collection:

            this_reach_len = river_reach.s.max() - river_reach.s.min()
            this_reach_id = river_reach.reach_indx[0]

            if enhanced:
                # Build up array of data to be smoothed from downstream to
                # upstream.  Adjust along-reach to be cumulative across
                # reach boundaries.
                first_node = 0
                distances = np.array([])
                heights = np.array([])

                if this_reach_id < n_reach:
                    reach_downstream = river_reach_collection[
                        ind.index(this_reach_id+1)]
                    distances = np.concatenate([
                        reach_downstream.s, distances])
                    heights = np.concatenate([
                        reach_downstream.h_n_ave, heights])

                distances = np.concatenate([
                    river_reach.s, distances+river_reach.s[-1]])
                heights = np.concatenate([river_reach.h_n_ave, heights])

                if this_reach_id > 1:
                    reach_upstream = river_reach_collection[
                        ind.index(this_reach_id-1)]
                    first_node = first_node + len(reach_upstream.h_n_ave)

                    distances = np.concatenate([
                        reach_upstream.s, distances+reach_upstream.s[-1]])
                    heights = np.concatenate([
                        reach_upstream.h_n_ave, heights])

                last_node = first_node + len(river_reach.h_n_ave) - 1

                # window size and sigma for Gaussian averaging
                window_size = np.min([max_window_size, this_reach_len])
                sigma = np.max([
                    min_sigma, window_size/window_size_sigma_ratio])

                # smooth h_n_ave, and get slope
                slope = np.polyfit(distances, heights, 1)[0]
                heights_detrend = heights - slope*distances
                heights_smooth = self.gaussian_averaging(
                    distances, heights_detrend, window_size, sigma)
                heights_smooth = heights_smooth + slope*(
                    distances - distances[0])
                enhanced_slopes.append(
                    (heights_smooth[last_node] - heights_smooth[first_node]
                    )/this_reach_len)

            else:
                enhanced_slopes.append(
                    (river_reach.h_n_ave[0] - river_reach.h_n_ave[-1])
                    /this_reach_len)
        return enhanced_slopes

    @staticmethod
    def gaussian_averaging(distances, heights, window_size, sigma):
        """
        Gaussian smoothing of heights using distances
        distances:   along-flow distance
        heights:     water heights
        window_size: size of data window to use for averaging
        sigma:       STD of Gaussian used for averaging

        outputs:
        smooth_heights : smoothed elevations
        """
        smooth_heights = np.zeros(heights.shape)
        for ii, this_distance in enumerate(distances):

            # get data window
            mask = np.logical_and(
                np.abs(this_distance-distances) <= window_size / 2,
                ~np.isnan(heights))

            weights = scipy.stats.norm.pdf(
                this_distance, distances[mask], sigma)

            smooth_heights[ii] = (
                np.multiply(weights, heights[mask]).sum() /
                weights.sum())

        return smooth_heights
