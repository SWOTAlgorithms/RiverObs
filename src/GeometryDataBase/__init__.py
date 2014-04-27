"""
Fast access to a data base containing a set of geometries. 

The data base is created by reading a shapefile, and write an rtree index for it for future use.

All bounding boxes are assumed to be an iterable (xmin,ymin,xmax,ymax)
"""

from GeometryDataBase import GeometryDataBase2D, GeometryDataBase3D
from GeometryDataBase import bbox_generator_3D, bbox_generator_2D,  shape_bbox_dbf_as_tuple
from GeometryDataBase import shape_bbox_as_tuple, write_shape_rtree_3D, write_shape_rtree_2D


from version import __version__
