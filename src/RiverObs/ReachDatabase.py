"""
Classes and methods to interact with the "new" reach database
"""

import os
import glob
import netCDF4
import logging
import numpy as np
from collections import OrderedDict as odict

from RiverObs.RiverReach import RiverReach
from SWOTRiver.products.product import Product, FILL_VALUES, textjoin

LOGGER = logging.getLogger(__name__)

class ReachExtractor(object):
    """
    Looks / acts / quacks like ReachExtractor.ReachExtractor
    """
    def __init__(
        self, reach_db_path, lat_lon_region, clip=True, clip_buffer=0.1):

        if os.path.isdir(reach_db_path):
            LOGGER.info('Extracting reaches')
            # figure out which db tiles to use
            reach_db = None
            for db_file in glob.glob(os.path.join(reach_db_path, '*.nc')):
                with netCDF4.Dataset(db_file, 'r') as ifp:
                    bbox = [ifp.x_min, ifp.y_min, ifp.x_max, ifp.y_max]
                    if (bbox[0] < lat_lon_region.bounding_box[2] and
                        bbox[2] > lat_lon_region.bounding_box[0] and
                        bbox[1] < lat_lon_region.bounding_box[3] and
                        bbox[3] > lat_lon_region.bounding_box[1]):

                        LOGGER.info('Using reach db tile {}'.format(db_file))
                        this_db = ReachDatabase.from_ncfile(db_file)
                        if reach_db is None:
                            reach_db = this_db
                        else:
                            reach_db = reach_db + this_db

        else:
            # assume already done
            reach_db = ReachDatabase.from_ncfile(reach_db_path)

        try_reach_idx = reach_db.reaches.extract(lat_lon_region.bounding_box)
        self.reach = []
        self.reach_idx = []
        for ii, reach_idx in enumerate(try_reach_idx):

            this_reach = reach_db(reach_idx)
            lon = this_reach['nodes']['x']
            lat = this_reach['nodes']['y']

            if clip:
                lonmin, latmin, lonmax, latmax = lat_lon_region.bounding_box
                clip_lon = lon.copy()

                # check for wraps
                if lonmax < lonmin: lonmax += 360
                clip_lon[clip_lon < this_reach['reaches']['x_min']] += 360

                inbbox = np.logical_and(
                    np.logical_and(clip_lon >= lonmin - clip_buffer,
                                   clip_lon <= lonmax + clip_buffer),
                    np.logical_and(lat >= latmin - clip_buffer,
                                   lat <= latmax + clip_buffer))

                lon = lon[inbbox]
                lat = lat[inbbox]

            if len(lon) == 0:
                continue

            x, y = lat_lon_region.proj(lon, lat)

            # TODO add stuff from DB herre
            metadata = {'lakeFlag': this_reach['reaches']['lakeflag'][0]}
            self.reach_idx.append(reach_idx)
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
    ATTRIBUTES = odict([
        ['x_min', {'dtype': 'f8' , 'value': None}],
        ['x_max', {'dtype': 'f8' , 'value': None}],
        ['y_min', {'dtype': 'f8' , 'value': None}],
        ['y_max', {'dtype': 'f8' , 'value': None}],
        ])
    GROUPS = odict([
        ['nodes', 'ReachDatabaseNodes'],
        ['reaches', 'ReachDatabaseReaches'],
        ['centerlines', 'ReachDatabaseCenterlines']
    ])

    def __call__(self, reach_id):
        """Returns dict-o-stuff for reach_id"""
        return {"nodes": self.nodes(reach_id),
                "reaches": self.reaches(reach_id)}

    def __add__(self, other):
        # hack it up
        klass = ReachDatabase()
        klass.nodes = self.nodes + other.nodes
        klass.reaches = self.reaches + other.reaches
        klass.centerlines = self.centerlines + other.centerlines
        klass.x_min = np.min([self.x_min, other.x_min])
        klass.x_max = np.max([self.x_max, other.x_max])
        klass.y_min = np.min([self.y_min, other.y_min])
        klass.y_max = np.max([self.y_max, other.y_max])
        return klass

class ReachDatabaseNodes(Product):
    """Prior Reach database nodes"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['depth', 2], ['nodes', 0]])
    DIMENSIONS_NODES = odict([['nodes', 0]])
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
        ['cl_ids',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['node_length',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ])

    def __call__(self, reach_id):
        """Returns dict-o-stuff for reach_id"""
        mask = self.reach_id == reach_id
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseNodes()
        for dset in ['x', 'y', 'node_id', 'reach_id', 'width', 'segment_id',
                     'node_length']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        klass.cl_ids = np.concatenate([self.cl_ids, other.cl_ids], 1)
        return klass

class ReachDatabaseReaches(Product):
    """Prior Reach database reaches"""
    ATTRIBUTES = odict([])
    DIMENSIONS = odict([['depth', 2], ['reach_neighbors', 4], ['reaches', 0]])
    DIMENSIONS_CLIDS = odict([['depth', 2], ['reaches', 0]])
    DIMENSIONS_REACH_UPDOWN = odict([['reach_neighbors', 4], ['reaches', 0]])
    DIMENSIONS_REACHES = odict([['reaches', 0]])
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
        ['segment_id',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['cl_ids',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_CLIDS]])],
        ['rch_id_up',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACH_UPDOWN]])],
        ['rch_id_dn',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACH_UPDOWN]])],
        ])

    def __call__(self, reach_id):
        """Returns dict of reach attributes for reach_id"""
        mask = self.reach_id == reach_id
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseReaches()
        for dset in ['x', 'x_min', 'x_max', 'y', 'y_min', 'y_max',
                     'reach_id', 'reach_length', 'lakeflag', 'segment_id']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)]))

        for dset in ['cl_ids', 'rch_id_up', 'rch_id_dn']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))
        return klass

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
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['depth', 4], ['points', 0]])
    DIMENSIONS_POINTS = odict([['points', 0]])
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_POINTS]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_POINTS]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['node_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['cl_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_POINTS]])],
        ])

    def __add__(self, other):
        klass = ReachDatabaseCenterlines()
        for dset in ['x', 'y', 'cl_id']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)]))

        for dset in ['reach_id', 'node_id']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))
        return klass
