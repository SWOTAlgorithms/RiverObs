"""
Classes and methods to interact with the "new" reach database
"""

import netCDF4
import numpy as np
from collections import OrderedDict as odict

from .RiverReach import RiverReach
from SWOTRiver.products.product import Product, FILL_VALUES, textjoin

class ReachExtractor(object):
    """
    Looks / acts / quacks like ReachExtractor.ReachExtractor
    """
    def __init__(
        self, reach_db_file, lat_lon_region, clip=True, clip_buffer=0.1):

        reach_db = ReachDatabase.from_ncfile(reach_db_file)
        self.reach_idx = reach_db.reaches.extract(lat_lon_region.bounding_box)
        self.reach = []

        for ii, reach_idx in enumerate(self.reach_idx):

            this_reach = reach_db(reach_id)
            lon, lat = this_reach['nodes']['x'], this_reach['nodes']['y']

            # TODO: check wrapping here
            if clip:
                inbbox = ((lon >= bbox[0] - clip_buffer) &
                          (lat >= bbox[1] - clip_buffer) &
                          (lon <= bbox[2] + clip_buffer) &
                          (lat <= bbox[3] + clip_buffer))
                lon = lon[inbbox]
                lat = lat[inbbox]

            x, y = lat_lon_region.proj(lon, lat)

            # TODO add stuff from DB herre
            metadata = {}
            self.reach.append(RiverReach(
                lon=lon, lat=lat, x=x, y=y, metadata=metadata,
                reach_index=ii))

        self.idx = 0
        self.nreaches = len(self.reach)

    def __iter__(self):
        """This and the next function define an iterator over reaches."""
        return self

    def __next__(self):  ## Python 3: def __next__(self)
        """This and the previous function define an iterator over reaches."""
        if self.idx >= self.nreaches:
            self.idx = 0
            raise StopIteration

        self.idx += 1
        return self.reach[self.idx - 1]

    next = __next__

    def __len__(self):
        """Number of reaches."""
        return self.nreaches

    def __getitem__(self, index):
        """Get reaches or slices of reaches."""
        return self.reach[index]


class ReachDatabase(Product):
    """Prior Reach database"""
    ATTRIBUTES = ['x_min', 'x_max', 'y_min', 'y_max',]
    GROUPS = odict([
        ['nodes', 'ReachDatabaseNodes'],
        ['reaches', 'ReachDatabaseReaches'],
        ['centerlines', 'ReachDatabaseCenterlines']
    ])

    def __call__(self, reach_id):
        """Returns dict-o-stuff for reach_id"""
        return {"nodes": self.nodes(reach_id),
                "reaches": self.reaches(reach_id)}

class ReachDatabaseNodes(Product):
    """Prior Reach database nodes"""
    ATTRIBUTES = []
    DIMENSIONS = odict([['nodes', 0], ['depth', 2]])
    DIMENSIONS_NODES = 
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['node_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['width',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['segment_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['cls_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['node_length',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['pass',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ])

    def __call__(self, reach_id):
        """Returns dict-o-stuff for reach_id"""
        mask = self.reach_id == reach_id
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

class ReachDatabaseReaches(Product):
    """Prior Reach database reaches"""
    ATTRIBUTES = []
    DIMENSIONS = odict([['reaches', 0], ['depth', 2]])
    DIMENSIONS_REACHES = 
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['x_min',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['x_max',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y_min',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y_max',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['reach_length',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['lakeflag',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['cls_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['pass',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ])

    def __call__(self, reach_id):
        """Returns dict of reach attributes for reach_id"""
        mask = self.reach_id == reach_id
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def extract(self, bounding_box):
        """
        Returns reach_ids for reaches that intersect an input bounding box

        bounding_box = [lonmin, latmin, lonmax, latmax]
        """
        lonmin, latmin, lonmax, latmax = bounding_box
        if lonmax < lonmin:
            lonmax += 360

        # iterate over reaches in self.reaches
        reach_zips = zip(
            self.x_min, self.y_min, self.x_max, self.y_max, self.reach_id)

        overlapping_reach_ids = []
        for reach_zip in reach_zips:

            reach_lonmin, reach_latmin, reach_lonmax, reach_latmax, reach_id =\
                reach_zip

            if reach_lonmax < reach_lonmin:
                reach_lonmax += 360

            # test for overlap (assumption is 2D decomp into two 1D intervals)
            if(latmin < reach_latmax and latmax > reach_latmin and
               lonmin < reach_lonmax and lonmax > reach_lonmin):
                overlapping_reach_ids.append(reach_id)
        return overlapping_reach_ids

class ReachDatabaseCenterlines(Product):
    """Prior Reach database centerlines"""
    ATTRIBUTES = []
    DIMENSIONS = odict([['points', 0]])
    DIMENSIONS_REACHES = 
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['node_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['cl_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ])
