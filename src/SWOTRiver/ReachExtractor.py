"""
Extract all of the reaches overlapping a given SWOTL2 bounding box and 
computes the coordinates in the same projection used by the SWOTL2 data.
The result is contained in a set of RiverReach objects, which can be clipped
to the same bounding box  as the data set.
"""

import numpy as N
import pysal
#from shapely.geometry import asShape, Polygon
from GeometryDataBase import GeometryDataBase2D

class RiverReach:
    """Class containing river reach coordinates and metadata."""

    def __init__(self,**kwds):
        """Initialize with empty information and, optionally,
        a set of keywords that get turned into class members. 
        The metadata is stored as dictionary variable_name:value."""

        self.lat = None
        self.lon = None
        self.x = None
        self.y = None
        self.metadata = {}

        for k in kwds:
            v = kwds[k]
            exec('self.%s = v'%k)

class ReachExtractor:
    """Extract all of the reaches overlapping a given SWOTL2 bounding box and 
    computes the coordinates in the same projection used by the SWOTL2 data.
    The result is contained in a set of RiverReach objects, which can be clipped
    to the same bounding box  as the data set."""

    def __init__(self, shape_file_root, swotL2,clip=True,clip_buffer=0.1):
        """Initialize with the shape file database location and a SWOTL2 instance.

        If clip is true, the reach is clipped to lie in a bounding box defined
        by the data bounding box plus a buffer given by clip_buffer 
        (default is 0.1 deg or ~11 km).
        """

        # Open the geometry data base and shape files

        self.db = GeometryDataBase2D(shape_file_root)

        # Open the shape and dbf files

        self.shp = pysal.open(shape_file_root+'.shp')
        self.dbf = pysal.open(shape_file_root+'.dbf')
        self.dbf_header = self.dbf.header

        # Get the list of applicable reaches and extract them

        self.reach_idx = self.db.intersects_xy_bbox(swotL2.bounding_box)

        # Store the reaches in a list of RiverReaches

        self.reach = []
        bbox = swotL2.bounding_box
        for i in self.reach_idx:

            # Get the coordinates as arrays

            lon, lat = N.asarray(self.shp[i].vertices).T

            # Clip the data

            if clip:
                inbbox = ( (lon >= bbox[0] - clip_buffer) &
                           (lat >= bbox[1] - clip_buffer) &
                           (lon <= bbox[2] + clip_buffer) &
                           (lat <= bbox[3] + clip_buffer) )
                lon = lon[inbbox]
                lat = lat[inbbox]

            # Project into the L2 projection

            x, y = swotL2.proj(lon, lat)

            # Get the metadata

            metadata = {}
            record = self.dbf[i][0]
            for i,field in enumerate(self.dbf_header):
                metadata[field] = record[i]

            # Append the river reach

            self.reach.append(RiverReach(lon=lon,lat=lat,x=x,y=y,metadata=metadata,reach_index=i))

        
        # Set the iterator indexes

        self.idx = 0
        self.nreaches = len(self.reach)
        
    def __iter__(self):
        """This and the next function define an iterator over reaches."""
        return self
    
    def next(self): ## Python 3: def __next__(self)
        """This and the previous function define an iterator over reaches."""

        if self.idx >= self.nreaches:
            self.idx = 0
            raise StopIteration
        
        self.idx += 1
        return self.reach[self.idx - 1]

    def __len__(self):
        """Number of reaches."""

        return self.nreaches

    def __getitem__(self,index):
        """Get reaches or slices of reaches."""

        return self.reach[index]
    

            
            
        
