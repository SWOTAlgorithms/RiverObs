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

from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9
from RiverObs.RiverReach import RiverReach
from SWOTWater.products.product import Product, FILL_VALUES, textjoin

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
        self, reach_db_path, lat_lon_region, clip=False, clip_buffer=0.1,
        day_of_year=None):

        if day_of_year is None:
            LOGGER.warn(
                "Day of year not specified, not extracting ice flag!")

        lonmin, latmin, lonmax, latmax = lat_lon_region.bounding_box
        # check for wraps
        if lonmax < lonmin: lonmax += 360

        if os.path.isdir(reach_db_path):
            LOGGER.info('Extracting reaches')
            # figure out which db tiles to use
            reach_db = ReachDatabase.from_dir(
                reach_db_path, lat_lon_region.bounding_box)

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

            # Remove centerline vertices that joint adjacent reaches and are
            # not close to nodes.
            if np.any(this_reach['centerlines']['is_extra_vertex']):
                is_extra_mask = this_reach['centerlines']['is_extra_vertex']
                not_extra_mask = np.logical_not(is_extra_mask)

                lons = this_reach['centerlines']['x'][not_extra_mask]
                lats = this_reach['centerlines']['y'][not_extra_mask]
                xx, yy = lat_lon_region.proj(lons, lats)

                found_start = 0
                found_stop = 0
                extra_indicies = is_extra_mask.nonzero()[0]
                for ii, extra_index in enumerate(extra_indicies):
                    try_lon = this_reach['centerlines']['x'][extra_index]
                    try_lat = this_reach['centerlines']['y'][extra_index]
                    try_xx, try_yy = lat_lon_region.proj(try_lon, try_lat)

                    dist_start = np.sqrt(
                        (try_xx-xx[0:found_start+1])**2 +
                        (try_yy-yy[0:found_start+1])**2)

                    dist_stop = np.sqrt(
                        (try_xx-xx[-(found_stop+1):])**2 +
                        (try_yy-yy[-(found_stop+1):])**2)

                    if any(dist_start < 50):
                        # put vertex before one it is closest to
                        cut_idx = np.argmin(dist_start)
                        found_start += 1

                    elif any(dist_stop < 50):
                        # put vertex after one it is closest to
                        cut_idx = len(xx)-found_stop+np.argmin(dist_stop)
                        found_stop += 1

                    else:
                        # skip this one
                        continue

                    lons = np.ma.concatenate([
                        lons[:cut_idx], [try_lon], lons[cut_idx:]])
                    lats = np.ma.concatenate([
                        lats[:cut_idx], [try_lat], lats[cut_idx:]])
                    xx = np.ma.concatenate([
                        xx[:cut_idx], [try_xx], xx[cut_idx:]])
                    yy = np.ma.concatenate([
                        yy[:cut_idx], [try_yy], yy[cut_idx:]])

                this_reach['centerlines']['x'] = lons
                this_reach['centerlines']['y'] = lats

            blocking_widths = get_blocking_widths(x, y)
            reach_metadata = {
                'lakeFlag': this_reach['reaches']['lakeflag'][0],
                'low_slope_flag': this_reach['reaches']['low_slope_flag'][0],
                'lon': this_reach['reaches']['x'][0],
                'lat': this_reach['reaches']['y'][0],
                'centerline_lon': this_reach['centerlines']['x'],
                'centerline_lat': this_reach['centerlines']['y'],
                }
            reach_metadata_keys = [
                'area_fits', 'discharge_models', 'reach_length', 'n_nodes',
                'wse', 'wse_var', 'width', 'width_var', 'n_chan_max',
                'n_chan_mod', 'grod_id', 'slope', 'dist_out', 'n_rch_up',
                'n_rch_down', 'rch_id_up', 'rch_id_dn', 'lakeflag', 'iceflag',
                'river_name'
            ]

            for key in reach_metadata_keys:
                if key in ['rch_id_up', 'rch_id_dn', 'area_fits',
                           'discharge_models']:
                    reach_metadata[key] = this_reach['reaches'][key]
                elif key in ['iceflag']:
                    if day_of_year is not None:
                        reach_metadata[key] = this_reach['reaches'][key][
                            day_of_year, 0]
                    else:
                        reach_metadata[key] = MISSING_VALUE_INT4
                else:
                    reach_metadata[key] = this_reach['reaches'][key][0]

            node_metadata_keys = [
                'node_length', 'wse', 'wse_var', 'width', 'width_var',
                'n_chan_max', 'n_chan_mod', 'grod_id', 'dist_out', 'wth_coef',
                'ext_dist_coef', 'river_name']

            node_metadata = {
                key: this_reach['nodes'][key] for key in node_metadata_keys}

            # replace NODATA with no_data
            node_metadata['river_name'][
                node_metadata['river_name']=='NODATA'] = 'no_data'
            if reach_metadata['river_name'] == 'NODATA':
                reach_metadata['river_name'] = 'no_data'

            self.reach_idx.append(reach_idx)
            self.reach.append(RiverReach(
                lon=lon, lat=lat, x=x, y=y, metadata=reach_metadata,
                reach_index=ii, node_indx=node_indx,
                blocking_widths=blocking_widths,
                **node_metadata))

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
        ['Name', {}], ['production_date', {}],
        ['pass_number', {}], ['tile_number', {}], ['swath_side', {}],
        ['version', {'dtype': 'str' , 'value': ''}],
        ])
    GROUPS = odict([
        ['nodes', 'ReachDatabaseNodes'],
        ['reaches', 'ReachDatabaseReaches'],
        ['centerlines', 'ReachDatabaseCenterlines']
    ])

    def subset(self, reach_ids):
        """Subsets the PRD by reach_ids"""
        klass = ReachDatabase()
        klass.nodes = self.nodes.subset(reach_ids)
        klass.reaches = self.reaches.subset(reach_ids)
        klass.centerlines = self.centerlines.subset(reach_ids)
        klass.x_min = np.min(klass.centerlines.x)
        klass.x_max = np.max(klass.centerlines.x)
        klass.y_min = np.min(klass.centerlines.y)
        klass.y_max = np.max(klass.centerlines.y)
        return klass


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

    @classmethod
    def from_dir(cls, reach_db_path, bounding_box):
        """
        Builds a ReachDatabase from a directory of ReachDatabases and a 
        bounding box.

        bounding_box = [min_lon, min_lat, max_lon, max_lat]
        """
        lonmin, latmin, lonmax, latmax = bounding_box

        # wrap to [180, 180) interval
        if lonmin > 180:
            lonmin -= 360
        if lonmax > 180:
            lonmax -= 360

        if lonmax < lonmin: lonmax += 360
        klass = None
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
                        this_db = cls.from_ncfile(db_file)
                    if klass is None:
                        klass = this_db
                    else:
                        klass = klass + this_db

        return klass

class ReachDatabaseNodes(Product):
    """Prior Reach database nodes"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['centerlines', 2], ['nodes', 0]])
    DIMENSIONS_NODES = odict([['nodes', 0]])
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['y',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['node_id',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_NODES]])],
        ['reach_id',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_NODES]])],
        ['cl_ids',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS]])],
        ['node_length',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['wse',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['wse_var',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['width',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['width_var',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['n_chan_max',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['n_chan_mod',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['obstr_type',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['grod_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['hfalls_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['dist_out',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['wth_coef',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['ext_dist_coef',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['facc',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['lakeflag',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['max_width',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['meander_length',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['sinuosity',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_NODES]])],
        ['manual_add',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['river_name',
         odict([['dtype', 'U254'], ['dimensions', DIMENSIONS_NODES]])],
        ])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, reach_ids):
        """Subsets the PRD nodes by reach_ids"""
        klass = ReachDatabaseNodes()
        mask = self.reach_id == reach_ids[0]
        for reach_id in reach_ids[1:]:
            mask = np.logical_or(mask, self.reach_id == reach_id)
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, reach_id):
        """Returns dict-o-stuff for reach_id"""
        mask = self.reach_id == reach_id
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseNodes()
        for dset in [
                'x', 'y', 'node_id', 'reach_id', 'node_length', 'wse',
                'wse_var', 'width', 'width_var', 'n_chan_max', 'n_chan_mod',
                'grod_id', 'dist_out', 'wth_coef', 'ext_dist_coef',
                'river_name', 'obstr_type', 'hfalls_id', 'facc', 'lakeflag',
                'max_width', 'meander_length', 'sinuosity', 'manual_add']:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        klass.cl_ids = np.ma.concatenate([self.cl_ids, other.cl_ids], 1)
        return klass

class ReachDatabaseReaches(Product):
    """Prior Reach database reaches"""
    ATTRIBUTES = odict([])
    GROUPS = odict([
        ['area_fits', 'ReachDatabaseReachAreaFits'],
        ['discharge_models', 'ReachDatabaseReachDischargeModelsGroup']])

    DIMENSIONS = odict([
        ['centerlines', 2], ['reach_neighbors', 4], ['julian_day', 0],
        ['reaches', 0], ['orbits', 0]])
    DIMENSIONS_CLIDS = odict([['centerlines', 2], ['reaches', 0]])
    DIMENSIONS_REACH_UPDOWN = odict([['reach_neighbors', 4], ['reaches', 0]])
    DIMENSIONS_REACHES = odict([['reaches', 0]])
    DIMENSIONS_ICEFLAG = odict([['julian_day', 0], ['reaches', 0]])
    DIMENSIONS_ORBITS = odict([['orbits', 0], ['reaches', 0]])
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['x_min',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['x_max',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y_min',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y_max',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['reach_id',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['reach_length',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['n_nodes',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['wse',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['wse_var',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['width',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['width_var',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['facc',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['n_chan_max',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['n_chan_mod',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['obstr_type',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['grod_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['hfalls_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['slope',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['dist_out',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ['n_rch_up',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['n_rch_down',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['rch_id_up',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_REACH_UPDOWN]])],
        ['rch_id_dn',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_REACH_UPDOWN]])],
        ['lakeflag',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['low_slope_flag',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['iceflag',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_ICEFLAG]])],
        ['swot_obs',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['swot_orbits',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_ORBITS]])],
        ['cl_ids',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_CLIDS]])],
        ['river_name',
         odict([['dtype', 'U254'], ['dimensions', DIMENSIONS_REACHES]])],
        ['max_width',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_REACHES]])],
        ])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, reach_ids):
        """Subsets the PRD reaches by reach_ids"""
        klass = ReachDatabaseReaches()
        mask = self.reach_id == reach_ids[0]
        for reach_id in reach_ids[1:]:
            mask = np.logical_or(mask, self.reach_id == reach_id)
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        outputs['area_fits'] = self.area_fits.subset(mask)
        outputs['discharge_models'] = self.discharge_models.subset(mask)
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, reach_id):
        """Returns dict of reach attributes for reach_id"""
        mask = self.reach_id == reach_id
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        outputs['area_fits'] = self.area_fits(mask)
        outputs['discharge_models'] = self.discharge_models(mask)
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseReaches()
        for dset in ['x', 'x_min', 'x_max', 'y', 'y_min', 'y_max', 'reach_id',
                    'reach_length', 'n_nodes', 'wse', 'wse_var', 'width',
                    'width_var', 'n_chan_max', 'n_chan_mod', 'grod_id',
                    'slope', 'dist_out', 'n_rch_up', 'n_rch_down', 'lakeflag',
                    'river_name', 'facc', 'obstr_type', 'hfalls_id',
                    'swot_obs', 'max_width']:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))

        for dset in [
                'cl_ids', 'rch_id_up', 'rch_id_dn', 'iceflag', 'swot_orbits']:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))

        klass.area_fits = self.area_fits + other.area_fits
        klass.discharge_models = self.discharge_models + other.discharge_models
        return klass

    def extract(self, bounding_box):
        """
        Returns reach_ids for reaches that intersect an input bounding box

        bounding_box = [lonmin, latmin, lonmax, latmax]
        """
        BUFFER = 0.25
        lonmin, latmin, lonmax, latmax = bounding_box
        if lonmax < lonmin:
            lonmax += 360

        # iterate over reaches in self.reaches
        reach_zips = zip(
            self.x_min-BUFFER, self.y_min-BUFFER, self.x_max+BUFFER,
            self.y_max+BUFFER, self.reach_id)

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

class ReachDatabaseReachDischargeModelsGroup(Product):
    ATTRIBUTES = odict()
    DIMENSIONS = odict()
    VARIABLES = odict()
    GROUPS = odict([
        ['unconstrained', 'ReachDatabaseReachDischargeModels'],
        ['constrained', 'ReachDatabaseReachDischargeModels']])

    def subset(self, mask):
        """Subsets ReachDatabaseReachDischargeModelsGroup by reach_ids"""
        klass = ReachDatabaseReachDischargeModelsGroup()
        for group in self.GROUPS:
            klass[group] = self[group].subset(mask)
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {}
        for group in self.GROUPS:
            outputs[group] = self[group](mask)
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseReachDischargeModelsGroup()
        for group in self.GROUPS:
            klass[group] = self[group] + other[group]
        return klass

class ReachDatabaseReachDischargeModels(Product):
    """class for prior reach database reach discharge_models datagroup"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict()
    VARIABLES = odict()
    GROUPS = odict([
        ['MetroMan', 'ReachDatabaseReachMetroMan'],
        ['BAM', 'ReachDatabaseReachBAM'],
        ['HiVDI', 'ReachDatabaseReachHiVDI'],
        ['MOMMA', 'ReachDatabaseReachMOMMA'],
        ['SADS', 'ReachDatabaseReachSADS'],
        ['SIC4DVar', 'ReachDatabaseReachSIC4DVar'],
    ])

    def subset(self, mask):
        """Subsets ReachDatabaseReachDischargeModels by reach_ids"""
        klass = ReachDatabaseReachDischargeModels()
        for group in self.GROUPS:
            klass[group] = self[group].subset(mask)
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {}
        for group in self.GROUPS:
            outputs[group] = self[group](mask)
        return outputs

    def __add__(self, other):
        klass = ReachDatabaseReachDischargeModels()
        for group in self.GROUPS:
            klass[group] = self[group] + other[group]
        return klass

class ReachDatabaseReachMetroMan(Product):
    """class for PRD reach MetroMan discharge model"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['Abar', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['Abar_stdev', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['ninf', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['ninf_stdev', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['p', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['p_stdev', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['ninf_p_cor', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['p_Abar_cor', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['ninf_Abar_cor', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['sbQ_rel', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ])
    GROUPS = odict([])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachMetroMan by reach_ids"""
        klass = ReachDatabaseReachMetroMan()
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachMetroMan()
        for dset in self.VARIABLES.keys():
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        return klass


class ReachDatabaseReachBAM(Product):
    """class for PRD reach BAM discharge model"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['Abar', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['n', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['sbQ_rel', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ])
    GROUPS = odict([])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachBAM by reach_ids"""
        klass = ReachDatabaseReachBAM()
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachBAM()
        for dset in self.VARIABLES.keys():
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        return klass

class ReachDatabaseReachHiVDI(Product):
    """class for PRD reach HiVDI discharge model"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['Abar', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['alpha', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['beta', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['sbQ_rel', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ])
    GROUPS = odict([])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachHiVDI by reach_ids"""
        klass = ReachDatabaseReachHiVDI()
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachHiVDI()
        for dset in self.VARIABLES.keys():
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        return klass

class ReachDatabaseReachMOMMA(Product):
    """class for PRD reach MOMMA discharge model"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['B', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['H', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['Save', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['sbQ_rel', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ])
    GROUPS = odict([])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachMOMMA by reach_ids"""
        klass = ReachDatabaseReachMOMMA()
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachMOMMA()
        for dset in self.VARIABLES.keys():
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        return klass

class ReachDatabaseReachSADS(Product):
    """class for PRD reach SADS discharge model"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['Abar', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['n', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['sbQ_rel', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ])
    GROUPS = odict([])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachSADS by reach_ids"""
        klass = ReachDatabaseReachSADS()
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachSADS()
        for dset in self.VARIABLES.keys():
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        return klass

class ReachDatabaseReachSIC4DVar(Product):
    """class for PRD reach SIC4DVar discharge model"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['Abar', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['n', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ['sbQ_rel', odict([['dtype', 'f8'], ['dimensions', DIMENSIONS]])],
        ])
    GROUPS = odict([])

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachSIC4DVar by reach_ids"""
        klass = ReachDatabaseReachSIC4DVar()
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, mask):
        """Returns dict of reach attributes for reach_id"""
        # HACK alert, mask is injected by ReachDatabaseReaches class
        # which is the parent group of this datagroup
        outputs = {
            key: self[key][mask] for key in self.VARIABLES.keys()}
        return outputs

    def __add__(self, other):
        """Adds other to self"""
        klass = ReachDatabaseReachSIC4DVar()
        for dset in self.VARIABLES.keys():
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        return klass


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

    for var in VARIABLES:
        if VARIABLES[var]['dtype'][0] in ['f', 'i']:
            VARIABLES[var]['_FillValue'] = -9999

    def subset(self, mask):
        """Subsets ReachDatabaseReachAreaFits by reach_ids"""
        klass = ReachDatabaseReachAreaFits()
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

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
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))
        for dset in ['h_break', 'w_break']:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))
        for dset in ['fit_coeffs',]:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)], 2))
        return klass

class ReachDatabaseCenterlines(Product):
    """Prior Reach database centerlines"""
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['centerlines', 4], ['points', 0]])
    DIMENSIONS_POINTS = odict([['points', 0]])
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_POINTS]])],
        ['y',
         odict([['dtype', 'f8'], ['dimensions', DIMENSIONS_POINTS]])],
        ['reach_id',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS]])],
        ['node_id',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS]])],
        ['cl_id',
         odict([['dtype', 'i8'], ['dimensions', DIMENSIONS_POINTS]])],
        ])

    def subset(self, reach_ids):
        """Subsets the PRD ReachDatabaseCenterlines by reach_ids"""
        klass = ReachDatabaseCenterlines()
        mask = np.any(self.reach_id == reach_ids[0], axis=0)
        for reach_id in reach_ids[1:]:
            mask = np.logical_or(
                mask, np.any(self.reach_id == reach_id, axis=0))
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        for key, value in outputs.items():
            klass[key] = value
        return klass

    def __call__(self, reach_id):
        """Returns dict of reach attributes for reach_id"""
        mask_all = self.reach_id == reach_id
        mask = np.any(mask_all, axis=0)
        mask_extra = np.any(mask_all[1:, mask], axis=0)
        outputs = {
            key: self[key][..., mask] for key in self.VARIABLES.keys()}
        outputs['is_extra_vertex'] = mask_extra

        # sort by cl_id
        idx = np.argsort(outputs['cl_id'])
        for key, value in outputs.items():
            outputs[key] = value[..., idx]

        return outputs

    def __add__(self, other):
        klass = ReachDatabaseCenterlines()
        for dset in ['x', 'y', 'cl_id']:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)]))

        for dset in ['reach_id', 'node_id']:
            setattr(klass, dset, np.ma.concatenate([
                getattr(self, dset), getattr(other, dset)], 1))
        return klass

class ReachDatabaseTile(ReachDatabase):
    """PRD reach database for a single RiverTile"""
    ATTRIBUTES = odict([
        ['x_min', {'dtype': 'f8', 'value': None}],
        ['x_max', {'dtype': 'f8', 'value': None}],
        ['y_min', {'dtype': 'f8', 'value': None}],
        ['y_max', {'dtype': 'f8', 'value': None}],
        ['Name', {}], ['production_date', {}],
        ['pass_number', {'dtype': 'i4', 'value': None}],
        ['tile_number', {'dtype': 'i4', 'value': None}],
        ['swath_side', {}],
        ])

    def subset(self, reach_ids):
        """Subsets the PRD by reach_ids"""
        klass = ReachDatabaseTile()
        klass.nodes = self.nodes.subset(reach_ids)
        klass.reaches = self.reaches.subset(reach_ids)
        klass.centerlines = self.centerlines.subset(reach_ids)
        klass.x_min = np.min(klass.centerlines.x)
        klass.x_max = np.max(klass.centerlines.x)
        klass.y_min = np.min(klass.centerlines.y)
        klass.y_max = np.max(klass.centerlines.y)
        return klass
