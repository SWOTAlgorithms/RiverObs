"""
Extract all of the reaches overlapping a given bounding box and
computes the coordinates in the same projection used by the data.
The result is contained in a set of RiverReach objects, which can be clipped
to the same bounding box  as the data set.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import pysal

from GeometryDataBase import GeometryDataBase2D
from .RiverReach import RiverReach


class ReachExtractor:
    """Extract all of the reaches overlapping a given bounding box and
    computes the coordinates in the same projection used by the data.
    The result is contained in a set of RiverReach objects, which can be clipped
    to the same bounding box  as the data set.

    Initialize with the shape file database location and a lat_lon_region
    instance.

    Parameters
    ----------

    shape_file_root : str
        path to shapefile database (no suffix)

    lat_lon_region : object
        an object providing the following members:
        lat_lon_region.bounding_box (lonmin,latmin,lonmax,latmax)
        lat_lon_region.proj: a pyproj.Proj projection (lon,lat) -> (x,y)
        and (x,y) -> (lon,lat) when called when called with inverse=True
    clip : bool
        Clip to the bounding box?
    clip_buffer : float
        buffer to add around bounding box.

    Notes
    ------

    If clip is true, the reach is clipped to lie in a bounding box defined
    by the data bounding box plus a buffer given by clip_buffer
    (default is 0.1 deg or ~11 km).

    """

    def __init__(self,
                 shape_file_root,
                 lat_lon_region,
                 clip=True,
                 clip_buffer=0.1):
        # Open the geometry data base and shape files
        #print('shape_file_root:',shape_file_root)
        self.db = GeometryDataBase2D(shape_file_root)

        # Open the shape and dbf files

        self.shp = pysal.open(shape_file_root + '.shp')
        self.dbf = pysal.open(shape_file_root + '.dbf')
        self.dbf_header = self.dbf.header

        # Get the list of applicable reaches and extract them
        self.shape_idx = self.db.intersects_xy_bbox(
            lat_lon_region.bounding_box)
        #print "####### SHAPE_IDX:",self.shape_idx
        self.reach_idx = []

        # Store the reaches in a list of RiverReaches

        self.reach = []
        bbox = lat_lon_region.bounding_box
        #print('bbox:',bbox)
        for i in self.shape_idx:

            # Get the coordinates as arrays

            lon, lat = np.asarray(self.shp[i].vertices).T

            # Clip the data

            if clip:
                inbbox = ((lon >= bbox[0] - clip_buffer) &
                          (lat >= bbox[1] - clip_buffer) &
                          (lon <= bbox[2] + clip_buffer) &
                          (lat <= bbox[3] + clip_buffer))
                lon = lon[inbbox]
                lat = lat[inbbox]

            # Project into the L2 projection

            x, y = lat_lon_region.proj(lon, lat)

            # Get the metadata and reach index
            # Brent Williams, May 2017: Changed a few things here to handle
            # newer river reach database (may have broken ability to read old
            # one though, havent tested)
            metadata = {}
            record = self.dbf[i][0]
            reach_index = i
            max_width = None
            for j, field in enumerate(self.dbf_header):
                metadata[field] = record[j]
                if field == 'reach_idx':  #old grwl way
                    reach_index = record[j]
                if field == 'reachID':  #new database
                    reach_index = record[j]
                if field == 'Reach_ID':  #osu centerline
                    reach_index = record[j]                    
                #if field == 'Wmean':#new database mean width
                #    max_width = record[j]
                #    print "max width:", max_width
            self.reach_idx.append(reach_index)

            #print "reachID:",reach_index
            #print "reach x:",x
            # Append the river reach
            #if max_width==None:
            self.reach.append(
                RiverReach(
                    lon=lon,
                    lat=lat,
                    x=x,
                    y=y,
                    metadata=metadata,
                    reach_index=reach_index))
            #else:
            #    self.reach.append(RiverReach(lon=lon,lat=lat,x=x,y=y,
            #                                 metadata=metadata,
            #                                 reach_index=reach_index,
            #                                 width_max=width_max))

        # Set the iterator indexes
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
