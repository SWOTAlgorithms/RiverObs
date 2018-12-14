"""
Access SWOT L2 data conveniently and provides LatLonRegion protocol.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import netCDF4
import pyproj
import logging


class SWOTL2:
    """
    Access SWOT L2 data conveniently. SWOTL2 implements the LatLonRegion
    object interfaces in that it provides a  bounding_box member and a
    proj function that goes from (lon,lat) -> (x,y) and backwards.

    Parameters
    ----------

    bounding_box : tuple or array_like
        should be of the form (lonmin,latmin,lonmax,latmax).
        If the bounding_box is provided, select only the data in the
        bounding box, otherwise, get the bounding box from the data
        itself.
    good_class : list or array_like
        a list of class labels where the data are returned.
    class_list : list
        a list of the class labels for what is considered good data
    lon_kwd, lat_kwd : str
        netcdf names for the longitudes and latitudes to be used
        for georeferencing the data set.
    class_kwd : str
        the netcdf name of the classification layer to use.
    min_points : int
        If the number of good points is less than this, raise an exception.

    Notes
    ------

    The final set of keywords are projection options for pyproj.
    A full list of projection options to set plus explanations of their
    meaning can be found here: https://trac.osgeo.org/proj/wiki/GenParms

    The default projection is Lambert equiarea, which has a proj4 string with
    the following parameters:

    +proj=laea
    +lat_0=Latitude at projection center, set to bounding box center lat
    +lon_0=Longitude at projection center, set to bounding box center lon
    +x_0=False Easting, set to 0
    +y_0=False Northing, set to 0
    """
    def __init__(self,
                 swotL2_file,
                 bounding_box=None,
                 class_list=[1],
                 lat_kwd='latitude',
                 lon_kwd='longitude',
                 class_kwd='classification',
                 rngidx_kwd='range_index',
                 aziidx_kwd='azimuth_index',
                 min_points=100,
                 project_data=True,
                 proj='laea',
                 x_0=0,
                 y_0=0,
                 lat_0=None,
                 lon_0=None,
                 ellps='WGS84',
                 subsample_factor=1,
                 **proj_kwds):
        self.logger = logging.getLogger(__name__+'.'+self.__class__.__name__)
        self.lat_kwd, self.lon_kwd = lat_kwd, lon_kwd
        self.subsample_factor = subsample_factor
        self.nc = netCDF4.Dataset(swotL2_file)
        self.logger.info('Dataset opened')

        self.set_index_and_bounding_box(
            bounding_box, lat_kwd, lon_kwd, class_list, class_kwd=class_kwd)

        self.logger.debug('Good data selected & bounding box calculated.')

        # Get reference locations for these data
        self.lat = self.get(lat_kwd)
        self.lon = self.get(lon_kwd)

        # wrap to [-180, 180) interval
        self.lon[self.lon >= 180] -= 360

        # Put in the radar/image coordinates too
        try:
            self.img_x = self.get(rngidx_kwd)
            self.img_y = self.get(aziidx_kwd)

        except KeyError:
            try:
                print(
                    'Cant Find range_index, or azimuth_index variables,'
                    ' assuming 2D-image image coordinates (like from a gdem)')
                Ny, Nx = np.shape(self.get(lat_kwd, use_index=False))
                ix, iy = np.meshgrid(np.arange(Nx), np.arange(Ny))
                self.img_x = ix[self.index]
                self.img_y = iy[self.index]

            except:
                print(
                    'WARNING: Input file does not contain range/azimuth index. '
                    'Functions relying on radar coordinates WILL break!')
                self.img_x = None
                self.img_y = None

        self.logger.debug('lat/lon read')

        # If not enough good points are found, raise Exception
        if len(self.lat) < min_points:
            raise Exception(
                'number of good points: %d smaller than required: %d' %
                (len(self.lat), min_points))

        # Project to a coordinate system
        if project_data:
            self.x, self.y = self.project(
                proj=proj,
                x_0=x_0,
                y_0=y_0,
                lat_0=lat_0,
                lon_0=lon_0,
                ellps=ellps,
                **proj_kwds)
            self.logger.debug('projection set and x,y calculated')

    def set_index_and_bounding_box(self,
                                   bounding_box,
                                   lat_kwd,
                                   lon_kwd,
                                   class_list,
                                   class_kwd='no_layover_classification'):
        """
        Set the index for the good data and computes the bounding box for
        the good data.

        class_kwd: the netcdf name of the classification layer to use.
        """
        # Read the classification and update the index
        self.klass = self.get(class_kwd, use_index=False)

        self.class_list = class_list

        self.class_index = (self.klass == class_list[0])
        for i in range(1, len(class_list)):
            self.class_index = self.class_index | (self.klass == class_list[i])

        self.logger.debug('Number of points in these classes: %d' %(
                          np.sum(self.class_index)))

        lat = self.get(lat_kwd, use_index=False)
        lon = self.get(lon_kwd, use_index=False)

        # wrap to [-180, 180) interval
        lon[lon >= 180] -= 360

        if bounding_box is not None:
            self.lonmin, self.latmin, self.lonmax, self.latmax = bounding_box

        else:
            self.lonmin = lon.min()
            self.latmin = lat.min()
            self.lonmax = lon.max()
            self.latmax = lat.max()

        self.index = ((lat >= self.latmin) & (lon >= self.lonmin) &
                      (lat <= self.latmax) & (lon <= self.lonmax))

        self.logger.debug('Number of points in bounding box: %d' %
                         (np.sum(self.index)))

        self.index = self.index & self.class_index
        self.logger.debug('Number of good: %d' % (np.sum(self.index)))

        lat = lat[self.index]
        lon = lon[self.index]

        self.bounding_box = (self.lonmin, self.latmin, self.lonmax,
                             self.latmax)

        # Update the classification
        self.klass = self.klass[self.index]

    def get(self, var, use_index=True):
        """
        Get the values of the variable var within the desired index of
        good sites.

        Subsamples input data based on self.subsample_factor (use for GDEMS!)
        """
        try:
            # try old style pixc first
            data = self.nc.variables[var][:]
        except KeyError:
            # then try new one with groups
            try:
                data = self.nc.groups['pixel_cloud'][var][:]
            except (IndexError, KeyError):
                raise KeyError

        if self.subsample_factor > 1:
            if len(data.shape) == 1:
                data = data[::self.subsample_factor]

            elif len(data.shape) == 2:
                # only subsample in azimuth
                data = data[::self.subsample_factor, :]

            else:
                raise Exception('Unexpected size of input data in SWOTL2::get')

        # self.index already subsampled in set_index_and_bounding_box
        if use_index:
            data = data[self.index]

        return data

    def getatt(self, attname):
        try:
            return getattr(self.nc, attname)
        except AttributeError:
            try:
                return getattr(self.nc.groups['pixel_cloud'], attname)
            except (KeyError, AttributeError):
                raise AttributeError

    def project(self,
                proj='laea',
                x_0=0,
                y_0=0,
                lat_0=None,
                lon_0=None,
                ellps='WGS84',
                **proj_kwds):
        """
        Get x and y coordinates for the selected points using a proj4
        projection.

        The proj4 template should replace lat_0, lon_0, x_0, y_0, proj, which
        are passed as keyword parameters (default as below).

        A full list of projection options to set plus explanations of their
        meaning can be found here: https://trac.osgeo.org/proj/wiki/GenParms

        The default projection is Lambert equiarea, which has a proj4
        string with the following parameters:

        +proj=laea
        +lat_0=Latitude at projection center, set to bounding box center lat
        +lon_0=Longitude at projection center, set to bounding box center lon
        +x_0=False Easting, set to 0
        +y_0=False Northing, set to 0
        """
        # Find lat_0 and lon_0 if not specified previously
        # Use center of bounding box instead of data centroid
        # in order to get same nodes for gdem and l2 file
        if lat_0 == None:
            lat_0 = (self.bounding_box[3] + self.bounding_box[1]) / 2.0

        if lon_0 == None:
            lon_0 = (self.bounding_box[2] + self.bounding_box[0]) / 2.0

        self.proj = pyproj.Proj(
            proj=proj,
            lat_0=lat_0,
            lon_0=lon_0,
            x_0=x_0,
            y_0=y_0,
            ellps=ellps,
            **proj_kwds)

        self.x, self.y = self.proj(self.lon, self.lat)
        return self.x, self.y
