'''
Copyright (c) 2022-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore, Brent Williams

This module is for various products we want to run through RiverObs for calval
activities.
'''

import os
import warnings
import logging
import datetime
import fiona
import pyproj
import numpy as np
from collections import OrderedDict as odict

from SWOTWater.products.product import Product

LOGGER = logging.getLogger(__name__)

class RiverNCProductMixIn(object):
    """MixIn class implementing some common methods for calval data"""
    def compute_bounding_box(self):
        mask = np.logical_and(
            ~np.isnan(self.latitude), ~np.isnan(self.longitude))
        return (self.longitude[mask].min(), self.latitude[mask].min(),
                self.longitude[mask].max(), self.latitude[mask].max())

class PressureTransducer(Product):
    """Class for pressure transducer data"""
    ATTRIBUTES = odict([
        ['ID',{'dtype':'str', 'value':''}],
        ['latitude',{'dtype':'f8', 'value':None}],
        ['longitude',{'dtype':'f8', 'value':None}],
        ['river_mile',{'dtype':'f8', 'value':None}],
        ])
    GROUPS = odict()
    DIMENSIONS = odict([['record', 0]])
    VARIABLES = odict([
        ['time', {'dtype': 'f8'}],
        ['wse', {'dtype': 'f8'}],
        ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

class PressureTransducers():
    """ 
    Class that is just a list of PressureTransducer objects
    with some methods to read the PT files (that have many PTs in them)
    and create the PT objects
    """
    pts = []
    num_pts = 0
    @classmethod
    def from_native(cls, pt_file):
        """
        read in the PT file and create the list of PT objects
        """
        with open(pt_file, 'r') as ifp:
            lines = ifp.readlines()
        # parse the header with positions of each transducer
        pos = {}
        idline = lines[0].split()
        northline = lines[1].split()
        eastline = lines[2].split()
        mileline = lines[3].split()
        for key, northing, easting, mile in zip(
                idline[3:], northline[2:],eastline[2:], mileline[3:]):
            pos[key] = (northing, easting, mile)
        # parse the data
        data = {}
        keys = lines[5].split()
        for key in keys:
            data[key] = []
        for line in lines[6:]:
            if not line.isspace():
                splits = line.split('\t')
                for key, value in zip(keys, splits):
                    if (value.isspace() or value is ''):
                        value = 'nan'
                    data[key].append(value)
        pts = []
        # now create the PT objects
        for key in pos.keys():
            this_pt = PressureTransducer()
            for d in data[key]:
                float(d)
            wse = np.array([float(d) for d in data[key]])
            # convert time to swot-time
            swot_tt = np.ones(wse.shape) * np.nan
            for ii, (date, time) in enumerate(
                    zip(data['YYYYMMDD'], data['HHMMSS'])):
                year = int(date[0:4])
                month = int(date[4:6])
                day = int(date[6:])
                hour = int(time[0:2])
                minute = int(time[2:4])
                second = int(time[4:])
                dt = datetime.datetime(year, month, day, hour, minute, second)
                swot_tt[ii] = (dt-datetime.datetime(2000,1,1)).total_seconds()
            # convert Northing, easting to lat/lon
            # assuming UTM zone 10
            northing = pos[key][0]
            easting = pos[key][1]
            myProj = pyproj.Proj(
                ("+proj=utm +zone=10 +north +ellps=WGS84 +datum=WGS84 "
                 "+units=m +no_defs"))
            lon, lat = myProj(easting, northing, inverse=True)
            # set the product
            this_pt.wse = wse
            this_pt.time_tai = swot_tt
            this_pt.river_mile = pos[key][2]
            this_pt.longitude = lon
            this_pt.latitude = lat
            this_pt.ID = key
            # append this pt object to list
            pts.append(this_pt)
        # instantiate this object
        PTs = cls()
        PTs.pts = pts
        PTs.num_pts = len(pts)
        return PTs

class Drifter(RiverNCProductMixIn, Product):
    """This is the GPS/GNSS drifter class/object of river profiles."""
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
        ['height_arp', {'dtype': 'f8'}],
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
        ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    @classmethod
    def from_any(cls, drifter_file):
        """try to read the type of drifter data depending on file format"""
        klass = None
        # TODO: maybe there is a better way of doing this...
        try:
            klass = cls.from_ncfile(drifter_file)
        except OSError:
            # not a nc file
            try:
                klass = cls.from_native(drifter_file)
            except (UnicodeDecodeError, KeyError):
                # not a text file
                try:
                    klass = cls.from_shp(drifter_file)
                except fiona.errors.DriverError:
                    # catch the fiona error and raise something more useful
                    raise NotImplementedError(
                        "input file type not supported by Drifter class: {}".format(drifter_file))
        return klass

    @classmethod
    def from_shp(cls, drifter_file):
        """converter for shp file GNSS drift data"""
        # read in the shape file
        with fiona.open(drifter_file) as ifp:
            orig_proj = pyproj.Proj(ifp.crs)
            records = list(ifp)

        # define the output projection to WGS83 to convert to lat/lon
        dest_proj = pyproj.Proj(init='EPSG:4269')

        # get the basic variables
        wse = np.array([r['properties']['wse_83elev'] for r in records])
        wse_flag = np.array([r['properties']['wse_flag'] for r in records])
        xy = np.array([r['geometry']['coordinates'] for r in records])
        lon, lat = pyproj.transform(orig_proj, dest_proj, xy[:,0], xy[:,1])

        # get the time in swot format
        swot_tt = np.ones(len(records)) * np.nan
        for ii, record in enumerate(records):
            this_date = str(record['properties']['date'])
            time_utc = record['properties']['time_utc']

            # covert time_utc to date and time
            date = '{}'.format(int(this_date[4:8])) + this_date[0:4]
            times = time_utc.split(':')

            # put leading zeros in if less than 10
            time = (
                str(int(times[0])).zfill(2) + str(int(times[1])).zfill(2) +
                str(int(times[2])).zfill(2))
            year = int(date[-4:])
            day = int(date[-6:-4])
            month = int(date[0:-6])
            hour = int(time[0:2])
            minute = int(time[2:4])
            second = int(time[4:])

            dt = datetime.datetime(year, month, day, hour, minute, second)
            swot_tt[ii] = (dt-datetime.datetime(2000,1,1)).total_seconds()

        #create Drifter instance
        # only values of 1 and 0, 1 is better quality, 0 is more suspect
        mask = wse_flag >= 0
        klass = cls()
        klass.wse = wse[mask]
        klass.height_water = wse[mask]
        klass.latitude = lat[mask]
        klass.longitude = lon[mask]
        klass.time = swot_tt[mask]

        # TODO: handle quality, for now assume everything is good
        klass.position_3drss_formal_error = np.zeros(np.shape(wse[mask]))
        return klass

    @classmethod
    def from_native(cls, drifter_file):
        """Creates Drifter class instance from native format text file"""
        with open(drifter_file, 'r') as ifp:
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

        klass = cls()
        klass.time = swot_tt
        klass.latitude = latitude
        klass.longitude = longitude

        # put wse in water_height and wse fields since we dont have cor fields
        klass.height_water = height
        klass.wse = height

        # TODO: Handle quality.  For now assume everything is good
        klass.position_3drss_formal_error = np.zeros(latitude.shape)
        return klass

    def split_profiles(outfile=None):
        """
        identify overlapping river profiles in the drifter data and break the
        current instance into multiple instances, optionally writting to 
        output files.
        """
        # TODO: inplement it
        return []

class SimplePixelCloud(RiverNCProductMixIn, Product):
    """
    This is the class/object that all the various cal/val datasets get converted
    to in order to be able to run through RiverObs.
    """
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
        ['xref_input_files',{}],
        ])
    GROUPS = odict()
    DIMENSIONS = odict([['record', 0]])
    VARIABLES = odict([
        ['time', {'dtype': 'f8'}],
        ['time_tai', {'dtype': 'f8'}],
        ['latitude', {'dtype': 'f8'}],
        ['longitude', {'dtype': 'f8'}],
        ['height', {'dtype': 'f8'}],
        ['geoid', {'dtype': 'f8'}],
        ['solid_earth_tide', {'dtype': 'f8'}],
        ['load_tide_fes', {'dtype': 'f8'}],
        ['load_tide_got', {'dtype': 'f8'}],
        ['pole_tide', {'dtype': 'f8'}],
        ['model_dry_tropo_cor', {'dtype': 'f8'}],
        ['model_wet_tropo_cor', {'dtype': 'f8'}],
        ['classification', {'dtype': 'u1'}],
        ['azimuth_index', {'dtype': 'i4'}],
        ['range_index', {'dtype': 'i4'}],
        ['pixel_area', {'dtype': 'f8'}],
        ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    @classmethod
    def from_any(cls, infile):
        """try to read and convert various file types/formats"""
        # TODO: maybe there is better way of doing this...
        try:
            # catch warnings as error to check if input file is a SimplePixelCloud product
            # TODO: there is probablly a better way of doing this...
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                klass = cls.from_ncfile(infile)
        except (OSError, UserWarning):
            try:
                klass = cls.from_drifter(infile)
            except NotImplementedError: 
                try:
                    klass = cls.from_geotif(infile)
                except (KeyError, AttributeError):
                    klass = cls.from_pressure_transducer(infile, None)
        return klass

    @classmethod
    def from_drifter(cls, drifter_file):
        """converter for various formats of drifter data"""
        drft = Drifter.from_any(drifter_file)
        klass = cls()

        # set the needed variables
        for key in klass.VARIABLES.keys():
            try:
                klass[key] = drft[key]
            except AttributeError:
                LOGGER.warning('key {} not in SimplePixelCloud'.format(key))

        # get the mask for good pixels
        good = drft.position_3drss_formal_error < 0.25

        # populate height with height_water
        klass.height = drft.height_water

        # create fake classification variable
        classification = np.zeros(drft.latitude.shape)
        classification[good] = 1
        klass.classification = classification

        # Make fake range and azimuth indices so that segmentation algorithms
        # work in riverobs processing
        range_index = np.zeros(drft.latitude.shape)
        azimuth_index = np.zeros(drft.latitude.shape) + 2
        rind = np.arange(len(classification[good]), dtype=int)
        range_index[good] = rind
        klass.azimuth_index = azimuth_index
        klass.range_index = range_index
        return klass

    @classmethod
    def from_geotif(cls, geotif_file):
        """converter for water height-only geotifs"""

        # try to import the geotif modules from the python repo
        try:
            import util.geotif
        except ModuleNotFoundError:
            LOGGER.warning("Unable to load util.geotif module!")
            raise ModuleNotFoundError

        # read the geotif
        wse, extent, chunk = util.geotif.read_geotif(geotif_file)
        tx = util.geotif.get_lat_lon_tx(geotif_file)

        # mask out the non-water pixels
        good = (wse[:,:,0]>-1000)
        lat, lon = util.geotif.get_lat_lon(tx, extent, chunk, mask=good)
        height = wse[:,:,0]

        # populate the variables
        klass = cls()
        klass.height = height[good]
        klass.wse = height[good]
        klass.latitude = lat[good]
        klass.longitude = lon[good]

        # create the fake classification
        classif = np.zeros_like(height)
        classif[good] = 1
        klass.classification = classif[good]

        #create the fake slant-plane indices
        M,N = np.shape(height)
        az = np.arange(M)
        r = np.arange(N)
        R, Az = np.meshgrid(r, az)
        klass.azimuth_index = Az[good]
        klass.range_index = R[good]

        # compute the pixel_area from the grid spacing
        # so that area/width estimates also work
        M, N = np.shape(height)
        dx = (extent[1]-extent[0])/N
        dy = (extent[3]-extent[2])/M
        klass.pixel_area = np.ones_like(height[good]) * np.abs(dx*dy)
        return klass

    @classmethod
    def from_pressure_transducer(cls, pt_file, swot_time, maxtime=1000):
        """pressure transducer converter for closest time"""
        pts = PressureTransducers.from_native(pt_file)
        wse = np.zeros(len(pts.pts)) * np.nan
        lat = np.zeros(len(pts.pts)) * np.nan
        lon = np.zeros(len(pts.pts)) * np.nan
        ptime = np.zeros(len(pts.pts)) * np.nan
        dtime = np.zeros(len(pts.pts)) * np.nan
        classification = np.ones(len(pts.pts))

        if swot_time is None:
            # use the first time of the first PT
            swot_time = pts.pts[0].time

        for ii, pt in enumerate(pts.pts):
            # compute the index of closest time
            dt = np.abs(pt.time - swot_time)
            ind = np.nanargmin(dt)
            wse[ii] = pt.wse[ind]
            lat[ii] = pt.latitude
            lon[ii] = pt.longitude
            ptime[ii] = pt.time[ind]
            dtime[ii] = dt[ind]

        # create the fake classification
        classification = np.ones(len(pts.pts))
        classification[dtime>maxtime] = 0

        # make fake azimuth and range indices...all disjoint
        azimuth_index = np.arange(len(pts.pts)) * 5
        range_index = np.arange(len(pts.pts)) * 5
        klass = cls()
        klass.height = wse
        klass.longitude = lon
        klass.latitude = lat
        klass.time = ptime
        klass.classification = classification
        klass.range_index = range_index
        klass.azimuth_index = azimuth_index
        return klass

    @classmethod
    def from_lidar(cls, lidar_file):
        """TODO: implement converter"""
        raise NotImplementedError(
            "SimplePixelCloud.from_lidar not implemented!")

    @classmethod
    def from_airborne_imagery(cls, water_mask_file):
        """TODO: implement converter"""
        raise NotImplementedError(
            "SimplePixelCloud.from_airborne_imagery not implemented!")
