"""
Module for running RiverObs on a SWOT L2 PixelCloud data product

Author (s): Alex Fore
"""
import sys
import os
import copy
import argparse
import warnings
import netCDF4
import numpy as np
import logging
import datetime

import RDF
import SWOTRiver.EstimateSWOTRiver
from SWOTRiver.products.rivertile import L2HRRiverTile
from SWOTRiver.products.pixcvec import L2PIXCVector
from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9
from SWOTRiver.errors import RiverObsException

LOGGER = logging.getLogger(__name__)

class L2PixcToRiverTile(object):
    """
    Class for running RiverObs on a SWOT L2 PixelCloud data product
    """
    def __init__(self, l2pixc_file, index_file, is_new_pixc=None):
        self.pixc_file = l2pixc_file
        self.index_file = index_file
        self.is_new_pixc = is_new_pixc
        self.node_outputs, self.reach_outputs = None, None

        # if is_new_pixc is not supplied, test pixc file to see if it is true
        if self.is_new_pixc is None:
            self.is_new_pixc = False
            with netCDF4.Dataset(self.pixc_file) as ifp:
                self.is_new_pixc = 'pixel_cloud' in ifp.groups

        # compute day of year
        try:
            with netCDF4.Dataset(self.pixc_file) as ifp:
                t_str_start = ifp.time_coverage_start
            datetime_start = datetime.datetime.strptime(
                t_str_start[0:10], '%Y-%m-%d')
            datetime_ = datetime.datetime.strptime(t_str_start[0:4], '%Y')
            self.day_of_year = (datetime_start-datetime_).days+1
        except (AttributeError, ValueError):
            LOGGER.warn('Unable to parse day of year from PIXC input file,'+
                        'not doing ice flagging!')
            self.day_of_year = None

    def load_config(self, config):
        """Copies config object into self's storage from main"""
        LOGGER.info('load_config')
        self.config = copy.deepcopy(config)
        self.config['subsample_factor'] = 1

    def compute_bounding_box(self, from_attrs=True):
        """Get bounding box of self.pixc_file"""
        LOGGER.info('compute_bounding_box')
        with netCDF4.Dataset(self.pixc_file, 'r') as ifp:
            if from_attrs:
                lat_keys = [a+'_'+b+'_latitude' for a in ('inner', 'outer')
                    for b in ('first', 'last')]
                lon_keys = [a+'_'+b+'_longitude' for a in ('inner', 'outer')
                    for b in ('first', 'last')]
                try:
                    lat = np.array([getattr(ifp, item) for item in lat_keys])
                    lon = np.array([getattr(ifp, item) for item in lon_keys])

                # if older pixc format
                except AttributeError:
                    lat_keys = [a+'_'+b+'_lat' for a in ('inner', 'outer')
                        for b in ('first', 'last')]
                    lon_keys = [a+'_'+b+'_lon' for a in ('inner', 'outer')
                        for b in ('first', 'last')]
                    lat = np.array([getattr(ifp, item) for item in lat_keys])
                    lon = np.array([getattr(ifp, item) for item in lon_keys])

            else:
                data_dict = (ifp.variables if not self.is_new_pixc else
                             ifp.groups['pixel_cloud'])

                lat = np.asarray(data_dict['latitude'][:])
                lon = np.asarray(data_dict['longitude'][:])

        # wrap to [-180, 180) interval
        lon[lon >= 180] -= 360

        mask = ~np.isnan(lat)
        return (lon[mask].min(), lat[mask].min(), lon[mask].max(),
                lat[mask].max())

    def validate_inputs(self):
        """Validates that the input products meet requirements"""
        LOGGER.info('validate_inputs')
        # Check for empty PIXC file (dimension points in /pixel_cloud/ group
        # is zero).
        with netCDF4.Dataset(self.pixc_file, 'r') as ifp:
            if ifp.groups['pixel_cloud'].dimensions['points'].size == 0:
                LOGGER.error('Input L2_HR_PIXC product has zero valid pixels!')
                raise RiverObsException(
                    'Input L2_HR_PIXC product has zero valid pixels!')

        # ...etc for more validation tests

    def do_river_processing(self):
        """Does the river processing"""
        LOGGER.info('do_river_processing')

        qual_words = ("geo_qual_wse_suspect",
                      "geo_qual_wse_degraded",
                      "geo_qual_wse_bad",
                      "class_qual_area_suspect",
                      "class_qual_area_degraded",
                      "class_qual_area_bad",
                      "sig0_qual_suspect",
                      "sig0_qual_bad")
        for word in qual_words:
            if word not in self.config:
                self.config[word] = 0x00000000

        if 'fractional_inundation_kwd' not in self.config:
            self.config['fractional_inundation_kwd'] = 'water_frac'

        if 'height_agg_method' not in self.config:
            self.config['height_agg_method'] = 'weight'

        if 'area_agg_method' not in self.config:
            self.config['area_agg_method'] = 'composite'

        if 'preseg_dilation_iter' not in self.config:
            self.config['preseg_dilation_iter'] = 0

        if 'slope_method' not in self.config:
            self.config['slope_method'] = 'bayes'

        if 'prior_unc_alpha' not in self.config:
            self.config['prior_unc_alpha'] = 1.5

        if 'char_length_tau' not in self.config:
            self.config['char_length_tau'] = 10000

        if 'use_multiple_reaches' not in self.config:
            self.config['use_multiple_reaches'] = False

        if 'use_ext_dist_coef' not in self.config:
            self.config['use_ext_dist_coef'] = True

        if 'pixc_qual_handling' not in self.config:
            self.config['pixc_qual_handling'] = None

        # Set sub-params of iterative_linear outlier flagging method to
        # None if outlier_method is not iterative_linear
        if self.config['outlier_method'] != 'iterative_linear':
            for key in ['outlier_abs_thresh', 'outlier_rel_thresh',
                        'outlier_upr_thresh', 'outlier_iter_num']:
                self.config[key] = None

        # key/value arguments for constructing SWOTRiverEstimator
        kwargs = {
            'bounding_box': self.compute_bounding_box(),
            'lat_kwd': 'latitude', 'lon_kwd': 'longitude',
            'class_kwd': 'classification', 'height_kwd': 'height',
            'rngidx_kwd': 'range_index', 'aziidx_kwd': 'azimuth_index',
            'class_list': self.config['class_list'],
            'xtrack_kwd': 'cross_track',
            'fractional_inundation_kwd': self.config['fractional_inundation_kwd'],
            'use_fractional_inundation': self.config['use_fractional_inundation'],
            'use_segmentation': self.config['use_segmentation'],
            'use_heights': self.config['use_heights'],
            'min_points': self.config['min_points'],
            'trim_ends': False, 'store_obs': False, 'store_reaches': False,
            'output_file': self.index_file,
            'proj': 'laea', 'x_0': 0, 'y_0': 0, 'lat_0': None, 'lon_0': None,
            'subsample_factor': 1,
            'height_agg_method': self.config['height_agg_method'],
            'area_agg_method': self.config['area_agg_method'],
            'preseg_dilation_iter': self.config['preseg_dilation_iter'],
            'slope_method': self.config['slope_method'],
            'prior_unc_alpha': self.config['prior_unc_alpha'],
            'char_length_tau': self.config['char_length_tau'],
            'use_multiple_reaches': self.config['use_multiple_reaches'],
            'use_ext_dist_coef': self.config['use_ext_dist_coef'],
            'outlier_method': self.config['outlier_method'],
            'outlier_abs_thresh': self.config['outlier_abs_thresh'],
            'outlier_rel_thresh': self.config['outlier_rel_thresh'],
            'outlier_upr_thresh': self.config['outlier_upr_thresh'],
            'outlier_iter_num': self.config['outlier_iter_num'],
            'pixc_qual_handling': self.config['pixc_qual_handling'],
            'geo_qual_wse_suspect': self.config['geo_qual_wse_suspect'],
            'geo_qual_wse_degraded': self.config['geo_qual_wse_degraded'],
            'geo_qual_wse_bad': self.config['geo_qual_wse_bad'],
            'class_qual_area_suspect': self.config['class_qual_area_suspect'],
            'class_qual_area_degraded': self.config['class_qual_area_degraded'],
            'class_qual_area_bad': self.config['class_qual_area_bad'],
            'sig0_qual_suspect': self.config['sig0_qual_suspect'],
            'sig0_qual_bad':self.config['sig0_qual_bad']
        }

        river_estimator = SWOTRiver.SWOTRiverEstimator(
            self.pixc_file, **kwargs)

        river_estimator.get_reaches(
            self.config['reach_db_path'], day_of_year=self.day_of_year)

        if len(river_estimator.reaches) == 0:
            LOGGER.info('No valid reaches in PRD for this PIXC data')
        else:
            self.reach_collection = river_estimator.process_reaches(
                minobs=self.config['minobs'],
                min_fit_points=self.config['min_fit_points'],
                enhanced=True)

            if len(self.reach_collection) > 0:
                reach_variables = list(self.reach_collection[0].metadata.keys())
                node_variables = list(self.reach_collection[0].__dict__.keys())
                node_variables.remove('ds')
                node_variables.remove('metadata')

                num_nodes_per_reach = [
                    len(item.lat) for item in self.reach_collection]
                num_nodes = sum(num_nodes_per_reach)

                self.node_outputs = {}
                self.reach_outputs = {}
                for node_variable in node_variables:
                    self.node_outputs[node_variable] = np.concatenate(
                        [getattr(reach, node_variable) for reach in
                         self.reach_collection])

                for reach_variable in reach_variables:
                    self.reach_outputs[reach_variable] = np.squeeze(np.array(
                        [reach.metadata[reach_variable] for reach in
                         self.reach_collection]))

                self.node_outputs['reach_idx'] = np.zeros(
                    self.node_outputs['lat'].shape).astype('int32')
                i_start = 0
                for ireach, num_nodes in enumerate(num_nodes_per_reach):
                    self.node_outputs['reach_idx'][
                        i_start:i_start + num_nodes] = ireach
                    i_start = i_start + num_nodes

            else:
                LOGGER.info('Reach collection has zero entries')

        # save for use later to fill in missing nodes/reaches
        self.prd_reaches = river_estimator.reaches

        # reformat index file to L2PIXCVector format
        pixcvec = L2PIXCVector.from_ncfile(self.index_file)
        pixcvec.update_from_pixc(self.pixc_file)
        pixcvec.to_ncfile(self.index_file)

    def do_improved_geolocation(self):
        """
        Uses output of river processing (nodes) and rare sensor data to
        improve geolocation on lat, lon datasets in index.nc file.
        """
        LOGGER.info('do_improved_geolocation')
        if (self.node_outputs is None or not
            self.config['do_improved_geolocation']):
            return

        if not self.is_new_pixc:
            print("Sensor information not provided, skipping improved ",
                  "geolocation")
            return

        try:
            import cnes.modules.geoloc.scripts.geoloc_river as geoloc_river
        except ModuleNotFoundError:
            print("Cant load CNES improved geolocation, skipping!")
            return

        cnes_sensor = geoloc_river.Sensor.from_pixc(self.pixc_file)

        # compute improved geolocation
        lat_corr, lon_corr, height_corr = geoloc_river.geoloc_river(
            geoloc_river.PixelCloud.from_file(self.pixc_file),
            geoloc_river.PixcvecRiver(self.index_file),
            cnes_sensor,
            geoloc_river.RiverTile.from_node_outputs(self.node_outputs),
            fit_heights_per_reach='fitted_height_from_riverobs',
            interpolate_pixc_between_nodes=True,
            method=self.config['geolocation_method'])

        # update geoloc in place in index file
        with netCDF4.Dataset(self.index_file, 'a') as ifp:
            ifp.variables['latitude_vectorproc'][:] = lat_corr
            ifp.variables['longitude_vectorproc'][:] = lon_corr
            ifp.variables['height_vectorproc'][:] = height_corr

    def match_pixc_idx(self):
        """Matches the pixels from pixcvector to input pixc"""
        LOGGER.info('match_pixc_idx')
        with netCDF4.Dataset(self.pixc_file, 'r') as ifp:

            if self.is_new_pixc:
                nr_pixels = ifp.groups['pixel_cloud'].interferogram_size_range
                azi_index = ifp.groups['pixel_cloud']['azimuth_index'][:]
                rng_index = ifp.groups['pixel_cloud']['range_index'][:]
            else:
                nr_pixels = ifp.nr_pixels
                azi_index = ifp.variables['azimuth_index'][:]
                rng_index = ifp.variables['range_index'][:]

            pixc_idx = np.array(azi_index * int(nr_pixels) + rng_index)

        with netCDF4.Dataset(self.index_file, 'a') as ofp:
            pixcvec_idx = np.array(
                ofp.variables['azimuth_index'][:] * int(nr_pixels) +
                ofp.variables['range_index'][:])

            indx, indx_pv, indx_pixc = np.intersect1d(
                pixcvec_idx, pixc_idx, return_indices=True)

            # re-order PIXCVecRiver datasets to ordering of pixc_index.
            for dset in ofp.variables.keys():
                data = ofp.variables[dset][:]
                ofp.variables[dset][:] = data[indx_pv]

            ofp.variables['pixc_index'][:] = indx_pixc.astype('int32')

    def build_products(self):
        """Constructs the L2HRRiverTile data product / updates the index file"""
        LOGGER.info('build_products')
        # If lake flag is set don't output width, area, or slope.
        try:
            for ireach, reach_id in enumerate(self.reach_outputs['reach_idx']):
                if self.reach_outputs['lake_flag'][ireach] != 0:
                    self.reach_outputs['slope'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['slope2'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['slope_u'][ireach] = MISSING_VALUE_FLT
                    # TODO mask out slope2_u and slope?_r_u datasets when available
                    self.reach_outputs['width'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['width_u'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['area'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['area_u'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['area_det'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['area_det_u'][ireach] = MISSING_VALUE_FLT
                    self.reach_outputs['area_of_ht'][ireach] = MISSING_VALUE_FLT

                    mask = self.node_outputs['reach_indx'] == reach_id
                    self.node_outputs['w_area'][mask] = MISSING_VALUE_FLT
                    self.node_outputs['width_u'][mask] = MISSING_VALUE_FLT
                    self.node_outputs['area'][mask] = MISSING_VALUE_FLT
                    self.node_outputs['area_det'][mask] = MISSING_VALUE_FLT
                    self.node_outputs['area_of_ht'][mask] = MISSING_VALUE_FLT
                    self.node_outputs['area_u'][mask] = MISSING_VALUE_FLT
                    self.node_outputs['area_det_u'][mask] = MISSING_VALUE_FLT

            self.rivertile_product = L2HRRiverTile.from_riverobs(
                self.node_outputs, self.reach_outputs, self.reach_collection,
                self.prd_reaches)

        except TypeError:
            LOGGER.warn('Output products are empty')
            self.rivertile_product = L2HRRiverTile()

        # add in a bunch more stuff from PIXC
        if not os.path.isfile(self.index_file):
            L2PIXCVector().to_ncfile(self.index_file)
        self.rivertile_product.update_from_pixc(
            self.pixc_file, self.index_file)

        history_string = "Created {}".format(
            datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S.%f'))

        pixcvec = L2PIXCVector.from_ncfile(self.index_file)
        pixcvec.update_from_rivertile(self.rivertile_product)
        pixcvec.update_from_pixc(self.pixc_file)
        pixcvec.history = history_string
        pixcvec.to_ncfile(self.index_file)

        self.rivertile_product.nodes.title = \
            "Level 2 KaRIn High Rate River Tile Node Data Product"
        self.rivertile_product.nodes.history = history_string
        self.rivertile_product.nodes.xref_l2_hr_pixc_files = self.pixc_file

        self.rivertile_product.reaches.title = \
            "Level 2 KaRIn High Rate River Tile Reach Data Product"
        self.rivertile_product.reaches.history = history_string
        self.rivertile_product.reaches.xref_l2_hr_pixc_files = self.pixc_file

        # Fixup some other things
        with netCDF4.Dataset(self.pixc_file, 'r') as ifp:
            pol = ifp.polarization

        self.rivertile_product.nodes.rdr_pol = \
            np.ones(self.rivertile_product.nodes.rdr_pol.shape, dtype='S1')
        self.rivertile_product.nodes.rdr_pol[:] = pol[0]
