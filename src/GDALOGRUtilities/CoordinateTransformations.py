"""
Classes to take objects from one coordinate system to another.
"""

from os.path import join
from osgeo import gdal, gdalconst, osr
from pyproj import Proj, transform
import numpy as N

class CoordinateTransformation:
    """Base class for coordinate transformations.

    Initialize ccordinate systems with source and destination projections as
    proj4 strings.
    """

    def __init__(self,source_projection,
                 destination_projection='+units=m +ellps=WGS84 +datum=WGS84 +proj=longlat '):
        self.source_projection = Proj(source_projection)
        self.destination_projection = Proj(destination_projection)


    def transform_xy(self,x,y):
        """Transform sequences of x, y coordinates in the source coordinate system
        to x,y in the target coordinate system."""
        
        return transform(self.source_projection,self.destination_projection,x,y)

