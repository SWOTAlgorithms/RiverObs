"""
Module for running RiverObs on a SWOT L2 PixelCloud data product

Author (s): Alex Fore
"""
import sys
import os
import copy
import argparse
import netCDF4
import numpy as np

import RDF
import SWOTRiver.EstimateSWOTRiver
import RiverObs.ShapeWriter
import RiverObs.NetCDFReachWriter

class L2PixcToRiverTile(object):
    """
    Class for running RiverObs on a SWOT L2 PixelCloud data product
    """
    def __init__(self, l2pixc_file, index_file, sensor_file=None):
        self.pixc_file = l2pixc_file
        self.index_file = index_file
        self.sensor_file = sensor_file

        self.is_new_pixc = False
        with netCDF4.Dataset(self.pixc_file) as ifp:
            self.is_new_pixc = 'pixel_cloud' in ifp.groups

    def load_config(self, config):
        """Copies config object into self's storage from main"""
        self.config = copy.deepcopy(config)
        self.config['subsample_factor'] = 1

    def compute_bounding_box(self, from_attrs=True):
        """Get bounding box of self.pixc_file"""
        with netCDF4.Dataset(self.pixc_file, 'r') as ifp:
            if from_attrs:
                lat_keys = [a+'_'+b+'_lat' for a in ('inner', 'outer') for b
                    in ('first', 'last')]
                lon_keys = [a+'_'+b+'_lon' for a in ('inner', 'outer') for b
                    in ('first', 'last')]

                attid = (ifp if not self.is_new_pixc else
                         ifp.groups['pixel_cloud'])

                lat = np.array([getattr(attid, item) for item in lat_keys])
                lon = np.array([getattr(attid, item) for item in lon_keys])
            else:
                data_dict = (ifp.variables if not self.is_new_pixc else
                             ifp.groups['pixel_cloud'])

                lat = np.asarray(data_dict['latitude'][:])
                lon = np.asarray(data_dict['longitude'][:])

        mask = ~np.isnan(lat)
        return (lon[mask].min(), lat[mask].min(), lon[mask].max(),
                lat[mask].max())

    def do_river_processing(self):
        """Does the river processing"""
        print(self.config['trim_ends'])

        # key/value arguments for constructing SWOTRiverEstimator
        kwargs = {
            'bounding_box': self.compute_bounding_box(),
            'lat_kwd': 'latitude', 'lon_kwd': 'longitude',
            'class_kwd': 'classification', 'height_kwd': 'height',
            'rngidx_kwd': 'range_index', 'aziidx_kwd': 'azimuth_index',
            'class_list': self.config['class_list'],
            'xtrack_kwd': 'cross_track',
            'fractional_inundation_kwd': 'continuous_classification',
            'use_fractional_inundation': self.config['use_fractional_inundation'],
            'use_segmentation': self.config['use_segmentation'],
            'use_heights': self.config['use_heights'],
            'min_points': self.config['min_points'],
            'trim_ends': self.config['trim_ends'],
            'verbose': True, 'store_obs': False,
            'store_reaches': False, 'store_fits': False,
            'output_file': self.index_file,
            'proj': 'laea', 'x_0': 0, 'y_0': 0, 'lat_0': None, 'lon_0': None,
            'subsample_factor': 1}

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

        self.node_outputs, self.reach_outputs = \
            RiverObs.NetCDFReachWriter.gather_outputs(self.reach_collection)

    def do_improved_geolocation(self):
        """
        Uses output of river processing (nodes) and rare sensor data to
        improve geolocation on lat, lon datasets in index.nc file.
        """
        if (self.node_outputs is None or not
            self.config['do_improved_geolocation']):
            return

        if self.sensor_file is None and not self.is_new_pixc:
            print("Sensor information not provided, skipping improved ",
                  "geolocation")

        try:
            import cnes.modules.geoloc.scripts.geoloc_river as geoloc_river
        except ModuleNotFoundError:
            print("Cant load CNES improved geolocation, skipping!")
            return

        if self.sensor_file is None:
            cnes_sensor = geoloc_river.Sensor.from_pixc(self.pixc_file)
        else:
            cnes_sensor = geoloc_river.Sensor.from_file(self.sensor_file)

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
        with netCDF4.Dataset(self.pixc_file, 'r') as ifp:

            if self.is_new_pixc:
                nr_pixels = ifp.groups['pixel_cloud'].nr_pixels
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
