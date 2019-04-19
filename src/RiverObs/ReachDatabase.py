"""
Classes and methods to interact with the "new" reach database
"""

import os
import glob
import netCDF4
import logging
import warnings
import numpy as np
import math
import pyproj

from collections import OrderedDict as odict

from RiverObs.RiverReach import RiverReach
from SWOTRiver.products.product import Product, FILL_VALUES, textjoin

LOGGER = logging.getLogger(__name__)

class LatLonRegion:
    """
    Hacked up class to use for cheap-n-easy use of ReachExtractor

    Look/act/quack like SWOTL2 at least for ReachExtractor's use
    """
    def __init__(self, bounding_box):
        self.bounding_box = bounding_box
        lat_0 = (self.bounding_box[3] + self.bounding_box[1]) / 2.0
        lon_0 = (self.bounding_box[2] + self.bounding_box[0]) / 2.0
        self.proj = pyproj.Proj(
            proj='laea', x_0=0, y_0=0, lat_0=lat_0, lon_0=lon_0, ellps='WGS84')

def get_blocking_widths(x, y):
    """
    Computes the blocking widths for nodes at x, y

    Assumes x, y are ordered from upstream to downstream (i.e. generally from
    highest point to lowest).

    widths are signed: looking from upstream to downstream, bending left is
    positive and bending right is negative.
    """
    widths = np.nan*np.ones(x.shape)

    for ii in range(len(x)):
        # consider a limited set of permutations of nodes that can block
        # this node.
        for joff_m in [1, 2, 3]:
            for joff_p in [1, 2, 3]:
                # skip first/last
                if ii - joff_m < 0 or ii + joff_p >= len(x):
                    continue

                x1, y1 = x[ii-joff_m], y[ii-joff_m]
                x2, y2 = x[ii], y[ii]
                x3, y3 = x[ii+joff_p], y[ii+joff_p]

                this_width = blocking_width(x1, y1, x2, y2, x3, y3)
                if (np.isnan(widths[ii]) or
                    np.abs(this_width) < np.abs(widths[ii])):
                    widths[ii] = this_width
    return widths

def blocking_width(x1, y1, x2, y2, x3, y3):
    """
    Computes width at which point (x2, y2) is blocked by other nodes
    i.e. point in plane eqi-distant to all three.
    """
    # points midway beteen adjacent nodes
    xa = 0.5*(x1+x2);
    ya = 0.5*(y1+y2);

    xb = 0.5*(x2+x3);
    yb = 0.5*(y2+y3);

    # construct normal to connecting vectors, which is the line of
    # equidistant points

    # normal to r2-f1:
    x_normal_21 = y1-y2;
    y_normal_21 = x2-x1;
    norm = math.sqrt(x_normal_21**2 + y_normal_21**2);
    x_normal_21 = x_normal_21 / norm;
    y_normal_21 = y_normal_21 / norm;

    # normal to r3-r2:
    x_normal_32 = y2-y3;
    y_normal_32 = x3-x2;
    norm = math.sqrt(x_normal_32**2 + y_normal_32**2);
    x_normal_32 = x_normal_32 / norm;
    y_normal_32 = y_normal_32 / norm;

    # test if points co-linear
    if x_normal_21 == x_normal_32 and y_normal_21 == y_normal_32:
        return np.Infinity, np.Infinity

    # Solve this set of equations:
    # xa+alpha*x_normal_21 = xb + beta*x_normal_32
    # ya+alpha*y_normal_21 = yb + beta*y_normal_32

    if x_normal_32 != 0:
        # doing it on paper gives me (normal solution):
        alpha = (
            ((yb-ya) + y_normal_32/x_normal_32 * (xa-xb)) /
            (y_normal_21-y_normal_32*x_normal_21/x_normal_32))
        beta = (xa-xb)/x_normal_32 + alpha * x_normal_21/x_normal_32
    else:
        # equations simplify if x_normal_32 == 0
        alpha = (xb-xa)/x_normal_21
        beta = ((ya-yb) + (xa-xb)*y_normal_21 / x_normal_21)/y_normal_32

    xc = xa+alpha*x_normal_21
    yc = ya+alpha*y_normal_21

    dist = math.sqrt((x2-xc)**2 + (y2-yc)**2)
    if alpha < 0:
        dist *= -1
    return dist

class ReachExtractor(object):
    """
    Looks / acts / quacks like ReachExtractor.ReachExtractor
    """
    def __init__(
        self, reach_db_path, lat_lon_region, clip=True, clip_buffer=0.1):

        lonmin, latmin, lonmax, latmax = lat_lon_region.bounding_box
        # check for wraps
        if lonmax < lonmin: lonmax += 360

        if os.path.isdir(reach_db_path):
            LOGGER.info('Extracting reaches')
            # figure out which db tiles to use
            reach_db = None
            for db_file in glob.glob(os.path.join(reach_db_path, '*.nc')):
                with netCDF4.Dataset(db_file, 'r') as ifp:

                    reach_lonmin, reach_lonmax = ifp.x_min, ifp.x_max
                    reach_latmin, reach_latmax = ifp.y_min, ifp.y_max

                    # check for wraps
                    if reach_lonmax < reach_lonmin: reach_lonmax += 360

                    if (reach_lonmin < lonmax and reach_lonmax > lonmin and
                        reach_latmin < latmax and reach_latmax > latmin):

                        LOGGER.info('Using reach db tile {}'.format(db_file))
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            this_db = ReachDatabase.from_ncfile(db_file)

                        if reach_db is None:
                            reach_db = this_db
                        else:
                            reach_db = reach_db + this_db

        else:
            # assume already done
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                reach_db = ReachDatabase.from_ncfile(reach_db_path)

        try_reach_idx = reach_db.reaches.extract(lat_lon_region.bounding_box)
        self.reach = []
        self.reach_idx = []
        for ii, reach_idx in enumerate(try_reach_idx):

            if ii % 100 == 0:
                LOGGER.debug('Appending reach {} of {}'.format(
                    ii, len(try_reach_idx)))
            this_reach = reach_db(reach_idx)
            lon = this_reach['nodes']['x']
            lat = this_reach['nodes']['y']
            node_indx = this_reach['nodes']['node_id']

            if clip:
                clip_lon = lon.copy()

                # check for wraps
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

            blocking_widths = get_blocking_widths(x, y)

            # TODO add stuff from DB here
            metadata = {
                'lakeFlag': this_reach['reaches']['lakeflag'][0],
                'lon': this_reach['reaches']['x'][0],
                'lat': this_reach['reaches']['y'][0],
                'centerline_lon': this_reach['centerlines']['x'],
                'centerline_lat': this_reach['centerlines']['y'],
                'area_fits': this_reach['reaches']['area_fits'],
                'rch_id_up': this_reach['reaches']['rch_id_up'],
                'rch_id_dn': this_reach['reaches']['rch_id_dn']}

            self.reach_idx.append(reach_idx)
            self.reach.append(RiverReach(
                lon=lon, lat=lat, x=x, y=y, metadata=metadata,
                reach_index=ii, node_indx=node_indx,
                blocking_widths=blocking_widths,
                width=this_reach['nodes']['width']))

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
                "reaches": self.reaches(reach_id),
                "centerlines": self.centerlines(reach_id)}

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
    DIMENSIONS = odict([['centerlines', 2], ['nodes', 0]])
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
    GROUPS = odict([['area_fits', 'ReachDatabaseReachAreaFits']])

    DIMENSIONS = odict([['centerlines', 2], ['reach_neighbors', 4], ['reaches', 0]])
    DIMENSIONS_CLIDS = odict([['centerlines', 2], ['reaches', 0]])
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
        outputs['area_fits'] = self.area_fits(mask)
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

        klass.area_fits = self.area_fits + other.area_fits
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

class ReachDatabaseReachAreaFits(Product):
    """class for prior reach database reach area_fits datagroup"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([
        ['nCoeffs', 0], ['nReg', 0], ['hbreak_dim', 0], ['reaches', 0]])
    DIMENSIONS_HW_BREAK = odict([['hbreak_dim', 0], ['reaches', 0]])
    DIMENSIONS_REACHES = odict([['reaches', 0]])
    DIMENSIONS_FIT = odict([['nCoeffs', 0], ['nReg', 0], ['reaches', 0]])
    VARIABLES = odict([
        ['h_break',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_HW_BREAK]])],
        ['w_break',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_HW_BREAK]])],
        ['h_variance',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['w_variance',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['hw_covariance',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['med_flow_area',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['h_err_stdev',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['w_err_stdev',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['h_w_nobs',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['fit_coeffs',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_FIT]])],
        ])

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachAreaFits()
        for dset in ['h_variance', 'w_variance', 'hw_covariance',
                     'med_flow_area', 'h_err_stdev', 'w_err_stdev', 'h_w_nobs']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        for dset in ['h_break', 'w_break']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))
        for dset in ['fit_coeffs',]:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)], 2))
        return klass

class ReachDatabaseCenterlines(Product):
    """Prior Reach database centerlines"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['centerlines', 4], ['points', 0]])
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

    def __call__(self, reach_id):
        """Returns dict of reach attributes for reach_id"""
        mask = np.any(self.reach_id == reach_id, axis=0)
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseCenterlines()
        for dset in ['x', 'y', 'cl_id']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)]))

        for dset in ['reach_id', 'node_id']:
            setattr(klass, dset, np.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))
        return klass
