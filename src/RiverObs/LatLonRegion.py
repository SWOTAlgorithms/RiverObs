"""
A base class implemention the LatLonRegion protocol.
"""

from __future__ import absolute_import, division, print_function
import pyproj

class LatLonRegion:
    """Access SWOT L2 data conveniently. SWOTL2 implements the LatLonRegion object
    interfaces in that it provides the following members:

    lat_lon_region.bounding_box: (lonmin,latmin,lonmax,latmax)

    lat_lon_region.proj: a pyproj.Proj projection (lon,lat) -> (x,y)
    and (x,y) -> (lon,lat) when called when called with inverse=True

    Parameters
    ----------

    lonmin : float
        Minimum longitude, in degrees.
    latmin : float
        Minimum latitude, in degrees.
    lonmax : float
        Maximum longitude, in degrees.
    latmax : float
        Maximum latitude, in degrees.

    Notes
    ------

    The final set of keywords are projection options for pyproj.
    A full list of projection options to set plus explanations of their
    meaning can be found here: https://trac.osgeo.org/proj/wiki/GenParms

    The default projection is Lambert equiarea, which has a proj4 string with the
    following parameters:

    +proj=laea
    +lat_0=Latitude at projection center, set to bounding box center lat
    +lon_0=Longitude at projection center, set to bounding box center lon
    +x_0=False Easting, set to 0
    +y_0=False Northing, set to 0
    """

    def __init__(self,
                 lonmin,
                 latmin,
                 lonmax,
                 latmax,
                 proj='laea',
                 x_0=0,
                 y_0=0,
                 lat_0=None,
                 lon_0=None,
                 ellps='WGS84',
                 **proj_kwds):

        self.lonmin, self.latmin, self.lonmax, self.latmax = lonmin, latmin, lonmax, latmax
        self.bounding_box = self.lonmin, self.latmin, self.lonmax, self.latmax

        # Find lat_0 and lon_0 if not specified previously

        if lat_0 == None:
            lat_0 = (latmax + latmin) / 2.

        if lon_0 == None:
            lon_0 = (lonmax + lonmin) / 2.

        self.proj = pyproj.Proj(
            proj=proj,
            lat_0=lat_0,
            lon_0=lon_0,
            x_0=x_0,
            y_0=y_0,
            ellps=ellps,
            **proj_kwds)
