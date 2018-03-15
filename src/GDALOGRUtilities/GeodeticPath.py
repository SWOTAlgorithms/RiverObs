"""
Utilities for calculating
"""

from __future__ import absolute_import, division, print_function

from pyproj import Geod
import numpy as np


class GeodeticPath:
    """Calculate the geodetic path between two points using pyproj Geod"""

    def __init__(self, lon0, lat0, lon1, lat1, ellps='WGS84', radians=False):
        """Initialize with the start and stop points."""

        self.geod = Geod(ellps=ellps)

        self.lon0 = lon0
        self.lat0 = lat0
        self.lon1 = lon1
        self.lat1 = lat1

        # Get the forward and backward azimuths and distance between the two points

        self.azimuth0, self.back_azimuth1, self.distance = self.geod.inv(
            lon0, lat0, lon1, lat1)

    def get_path_lonlats(self, separation):
        """Get the latitude and longitude of the path points, given a point separation in meters."""

        self.npts = int(self.distance / separation + 0.5)
        self.lonlats = np.array(
            self.geod.npts(self.lon0, self.lat0, self.lon1, self.lat1,
                           self.npts))
        self.lon = self.lonlats[:, 0]
        self.lat = self.lonlats[:, 1]

    def get_path(self, separation=None):
        """Get the path latitude, longitude, heading, and distance.
        If separation is None, assume get_path_lon_lats has already been called."""

        if separation != None:
            self.get_path_lonlats(separation)

        # Declare heading and distance arrays

        self.heading = np.zeros(len(self.lat), dtype=self.lat.dtype)
        self.along_track_distance = np.zeros(
            len(self.lat), dtype=self.lat.dtype)

        for i in range(len(self.lat) - 1):
            self.heading[i], back_azimuth, d = self.geod.inv(
                self.lon[i], self.lat[i], self.lon[i + 1], self.lat[i + 1])
            self.along_track_distance[i + 1] = self.along_track_distance[i] + d

        if back_azimuth < 0:
            self.heading[-1] = 180 + back_azimuth
        else:
            self.heading[-1] = back_azimuth - 180.


class GeodeticPathFromPegPoint(GeodeticPath):
    """Calulate a geodetic path given a peg-point and a path length."""

    def __init__(self,
                 peg_lat,
                 peg_lon,
                 peg_heading,
                 distance,
                 peg_h=0.,
                 mode='center',
                 ellps='WGS84',
                 radians=False):
        """Initialize with a peg-point and a heading and a distance.

        mode: 'center', 'start', or 'end'. The peg-point location in the path.
        ellps: ellipsoid.
        """

        self.peg_lat = peg_lat
        self.peg_lon = peg_lon
        self.peg_heading = peg_heading
        self.distance = distance
        self.peg_h = peg_h
        self.mode = mode

        self.geod = Geod(ellps=ellps)

        if mode == 'start':
            lat0 = peg_lat
            lon0 = peg_lon
            lon1, lat1, hdg = self.geod.fwd(lon0, lat0, peg_heading, distance)
        elif mode == 'end':
            lat1 = peg_lat
            lon1 = peg_lon
            hdg = (peg_heading + 180.) % 360.
            lon0, lat0, hdg = self.geod.fwd(lon1, lat1, hdg, distance)
        else:
            hdg = (peg_heading + 180.) % 360.
            lon0, lat0, hdg = self.geod.fwd(peg_lon, peg_lat, hdg,
                                            distance / 2.)
            lon1, lat1, hdg = self.geod.fwd(peg_lon, peg_lat, peg_heading,
                                            distance / 2.)

        GeodeticPath.__init__(
            self, lon0, lat0, lon1, lat1, ellps=ellps, radians=False)
