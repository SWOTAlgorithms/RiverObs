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
    ATTRIBUTES = odict()
    GROUPS = odict()
    DIMENSIONS = odict([['record', 0]])
    VARIABLES = odict([
        ['time', {'dtype': 'f8'}],
        ['latitude', {'dtype': 'f8'}],
        ['longitude', {'dtype': 'f8'}],
        ['height', {'dtype': 'f8'}],
        ['classification', {'dtype': 'u1'}]
        ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

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

        # Create dummy classification field filled with contant value
        classification = np.zeros(latitude.shape)

        klass = cls()
        klass.time = swot_tt
        klass.latitude = latitude
        klass.longitude = longitude
        klass.height = height
        klass.classification = classification
        return klass


