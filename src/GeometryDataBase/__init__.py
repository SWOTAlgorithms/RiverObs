"""
Fast access to a data base containing a set of geometries. 

The data base is created by reading a shapefile, and write an rtree index for it for future use.

All bounding boxes are assumed to be an iterable (xmin,ymin,xmax,ymax)
"""

from GeometryDataBase import GeometryDataBase
from version import __version__
