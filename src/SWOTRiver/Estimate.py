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

import RDF
import SWOTRiver.EstimateSWOTRiver
from SWOTRiver.products.rivertile import L2HRRiverTile
from SWOTRiver.products.pixcvec import L2PIXCVector

LOGGER = logging.getLogger(__name__)

class L2PixcToRiverTile(object):
    """
    Class for running RiverObs on a SWOT L2 PixelCloud data product
    """
    def __init__(self, l2pixc_file, index_file, is_new_pixc=None):
        self.pixc_file = l2pixc_file
        self.index_file = index_file
        self.is_new_pixc = is_new_pixc

        # if is_new_pixc is not supplied, test pixc file to see if it is true
        if self.is_new_pixc is None:
            self.is_new_pixc = False
            with netCDF4.Dataset(self.pixc_file) as ifp:
                self.is_new_pixc = 'pixel_cloud' in ifp.groups

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

    def do_river_processing(self):
        """Does the river processing"""
        LOGGER.info('do_river_processing')
        print(self.config['trim_ends'])

        if 'fractional_inundation_kwd' not in self.config:
            self.config['fractional_inundation_kwd'] = 'water_frac'

        if 'height_agg_method' not in self.config:
            self.config['height_agg_method'] = 'weight'

        if 'area_agg_method' not in self.config:
            self.config['area_agg_method'] = 'composite'

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
            'trim_ends': self.config['trim_ends'],
            'store_obs': False, 'store_reaches': False, 'store_fits': False,
            'output_file': self.index_file,
            'proj': 'laea', 'x_0': 0, 'y_0': 0, 'lat_0': None, 'lon_0': None,
            'subsample_factor': 1, 
            'height_agg_method': self.config['height_agg_method'],
            'area_agg_method': self.config['area_agg_method'],}

        river_estimator = SWOTRiver.SWOTRiverEstimator(
            self.pixc_file, **kwargs)

        river_estimator.get_reaches(
            self.config['shape_file_root'],
            clip_buffer=self.config['clip_buffer'])

        if self.config['use_width_db']:
            river_estimator.get_width_db(self.config['width_db_file'])

        self.reach_collection = river_estimator.process_reaches(
            scalar_max_width=self.config['scalar_max_width'],
            minobs=self.config['minobs'],
            min_fit_points=self.config['min_fit_points'],
            fit_types=self.config['fit_types'],
            use_width_db=self.config['use_width_db'],
            ds=self.config['ds'],
            refine_centerline=self.config['refine_centerline'],
            smooth=self.config['smooth'],
            alpha=self.config['alpha'],
            max_iter=self.config['max_iter'])

        self.node_outputs, self.reach_outputs = None, None
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
                self.reach_outputs[reach_variable] = np.array(
                    [reach.metadata[reach_variable] for reach in
                     self.reach_collection])

            self.node_outputs['reach_idx'] = np.zeros(
                self.node_outputs['lat'].shape).astype('int32')
            i_start = 0
            for ireach, num_nodes in enumerate(num_nodes_per_reach):
                self.node_outputs['reach_idx'][
                    i_start:i_start + num_nodes] = ireach
                i_start = i_start + num_nodes
        else:
            warnings.warn('Reach collection has zero entries')

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
            fit_heights_per_reach=True, interpolate_pixc_between_nodes=True,
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
                ifgram_shape = ifp.groups['pixel_cloud'].interferogram_shape
                splits = ifgram_shape.replace(',', '').split()
                nr_pixels = int(splits[1])
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

            pixc_index = np.array(np.where(
                np.in1d(pixc_idx, pixcvec_idx))).astype('int32')[0]

            var = ofp.createVariable(
                'pixc_index', pixc_index.dtype.str, ('record', ))
            var[:] = pixc_index[:]

    def flag_lakes_pixc(self):
        """
        Adds a lake/river flag from prior database to pixcvector
        """
        LOGGER.info('flag_lakes_pixc')
        with netCDF4.Dataset(self.index_file, 'a') as ofp:
            pixc_reach = ofp.variables['reach_index'][:]

            # make lake_flag dataset and fill with 255
            var_lake_flag = ofp.createVariable(
                'lake_flag', np.dtype('uint8'), ('record', ))
            var_lake_flag[:] = 255*np.ones(pixc_reach.shape)

            if self.reach_outputs is not None:
                for reach, lake_flag in zip(self.reach_outputs['reach_idx'],
                                            self.reach_outputs['lake_flag']):
                    var_lake_flag[pixc_reach == reach] = lake_flag

    def build_products(self):
        """Constructs the L2HRRiverTile data product / updates the index file"""
        LOGGER.info('build_products')
        self.rivertile_product = L2HRRiverTile.from_riverobs(
            self.node_outputs, self.reach_outputs, self.reach_collection)

        # add in a bunch more stuff from PIXC
        self.rivertile_product.update_from_pixc(
            self.pixc_file, self.index_file)

        # copy attributes from pixel cloud to pixel cloud vector
        with netCDF4.Dataset(self.index_file, 'a') as ofp,\
             netCDF4.Dataset(self.pixc_file, 'r') as ifp:

            # put lon in [0, 360)
            lon = ofp.variables['longitude_vectorproc'][:]
            lon[lon < 0] += 360
            ofp.variables['longitude_vectorproc'][:] = lon

            for attr in L2PIXCVector.ATTRIBUTES.keys():
                try:
                    value = getattr(ifp, attr)
                except AttributeError:
                    value = getattr(ifp.groups['pixel_cloud'], attr, 'None')
                setattr(ofp, attr, value)
