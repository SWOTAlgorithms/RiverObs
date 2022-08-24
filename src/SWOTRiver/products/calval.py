'''
Copyright (c) 2022-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore


This module is for various products we want to run through RiverObs for calval
activities.
'''

import os
import datetime
import numpy as np
from collections import OrderedDict as odict

from SWOTWater.products.product import Product

class RiverNCProductMixIn(object):
    """MixIn class implementing some common methods for calval data"""
    def compute_bounding_box(self):
        mask = np.logical_and(
            ~np.isnan(self.latitude), ~np.isnan(self.longitude))
        return (self.longitude[mask].min(), self.latitude[mask].min(),
                self.longitude[mask].max(), self.latitude[mask].max())

class GPSProfile(RiverNCProductMixIn, Product):
    ATTRIBUTES = odict([
        ['Conventions',{}],
        ['title',{}],
        ['institution',{}],
        ['source',{}],
        ['history',{}],
        ['platform',{}],
        ['description',{}],
        ['ellipsoid_semi_major_axis',{}],
        ['ellipsoid_flattening',{}],
        ['xref_input_gnss_files',{}],
        ['xref_input_height_offset_files',{}],
        ])
    GROUPS = odict()
    DIMENSIONS = odict([['record', 0]])
    VARIABLES = odict([
        ['time', {'dtype': 'f8'}],
        ['time_tai', {'dtype': 'f8'}],
        ['latitude', {'dtype': 'f8'}],
        ['longitude', {'dtype': 'f8'}],
        ['height', {'dtype': 'f8'}],
        ['x', {'dtype': 'f8'}],
        ['y', {'dtype': 'f8'}],
        ['z', {'dtype': 'f8'}],
        ['position_3drss_formal_error', {'dtype': 'f8'}],
        ['height_water', {'dtype': 'f8'}],
        ['wse', {'dtype': 'f8'}],
        ['geoid', {'dtype': 'f8'}],
        ['solid_earth_tide', {'dtype': 'f8'}],
        ['load_tide_fes', {'dtype': 'f8'}],
        ['load_tide_got', {'dtype': 'f8'}],
        ['pole_tide', {'dtype': 'f8'}],
        ['height_offset_cor', {'dtype': 'f8'}],
        ['model_dry_tropo_cor', {'dtype': 'f8'}],
        ['model_wet_tropo_cor', {'dtype': 'f8'}],
        ['gnss_wet_tropo_cor', {'dtype': 'f8'}],
        ['classification', {'dtype': 'u1'}],
        ['azimuth_index', {'dtype': 'i4'}],
        ['range_index', {'dtype': 'i4'}]
        ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    @classmethod
    def from_ncfile(cls, gps_profile_file):
        product = super(GPSProfile, cls).from_ncfile(gps_profile_file)
        # Create dummy classification field filled with constant value of 0
        classification = np.zeros(product.latitude.shape)
        # use the position_3drss_formal_error to throw out bad data
        # by setting classification
        classification[product.position_3drss_formal_error > 0.25] = 1
        setattr(product, 'classification', classification)
        # make fake range and azimuth indices so that segmentation algorithms work
        # in riverobs processing
        range_index = np.zeros(product.latitude.shape)
        azimuth_index = np.zeros(product.latitude.shape) + 2
        rind = np.arange(len(classification[classification==0]), dtype=int)
        range_index[classification==0] = rind
        setattr(product, 'azimuth_index', azimuth_index)
        setattr(product, 'range_index', range_index)
        return product

    @classmethod
    def from_native(cls, gps_profile_file):
        """Creates GPSProfile class instance from native format text file"""
        with open(gps_profile_file, 'r') as ifp:
            lines = ifp.readlines()

        data = {}
        keys = lines[0].split()
        for key in keys:
            data[key] = []

        for line in lines[1:]:
            splits = line.split()
            for key, value in zip(keys, splits):
                data[key].append(value)

        latitude = np.array([float(d) for d in data['latitude']])
        longitude = np.array([float(d) for d in data['longitude']])
        height = np.array([float(d) for d in data['WS_83E']])
        swot_tt = np.ones(latitude.shape) * np.nan

        for ii, (date, time) in enumerate(zip(data['Date'], data['Time'])):
            year = int(date[-4:])
            day = int(date[-6:-4])
            month = int(date[0:-6])
            hour = int(time[0:2])
            minute = int(time[2:4])
            second = int(time[4:])

            dt = datetime.datetime(year, month, day, hour, minute, second)
            swot_tt[ii] = (dt-datetime.datetime(2000,1,1)).total_seconds()

        # Create dummy classification field filled with constant value
        classification = np.zeros(latitude.shape)

        # create dummy range and azimuth indices to enable segmentation
        # algorithms in riverobs processing
        range_index = np.arange(latitude.shape, dtype=int)
        azimuth_index = np.zeros(latitude.shape) + 2
        
        klass = cls()
        klass.time = swot_tt
        klass.latitude = latitude
        klass.longitude = longitude
        klass.height_water = height # put it in water_height since we dont have cor fields
        klass.classification = classification
        klass.azimuth_index = azimuth_index
        klass.range_index = range_index
        return klass


