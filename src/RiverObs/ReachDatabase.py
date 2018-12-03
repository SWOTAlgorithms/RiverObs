"""
Classes and methods to interact with the "new" reach database
"""

import netCDF4
import numpy as np

from SWOTRiver.products.product import Product, FILL_VALUES, textjoin

class ReachDatabase(Product):
    """Prior Reach database"""
    ATTRIBUTES = ['x_min', 'x_max', 'y_min', 'y_max',]
    GROUPS = odict([
        ['nodes', 'ReachDatabaseNodes'],
        ['reaches', 'ReachDatabaseReaches'],
        ['centerlines', 'ReachDatabaseCenterlines']
    ])

class ReachDatabaseNodes(Product):
    """Prior Reach database nodes"""
    ATTRIBUTES = []
    DIMENSIONS = odict([['nodes', 0], ['depth', 2]])
    DIMENSIONS_NODES = 
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['node_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['width',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['segment_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ['cls_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['node_length',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_NODES]])],
        ['pass',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_NODES]])],
        ])

class ReachDatabaseReaches(Product):
    """Prior Reach database reaches"""
    ATTRIBUTES = []
    DIMENSIONS = odict([['reaches', 0], ['depth', 2]])
    DIMENSIONS_REACHES = 
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['x_min',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['x_max',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y_min',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['y_max',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['reach_length',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['lakeflag',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ['cls_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['pass',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS_REACHES]])],
        ])

    def extract(self, bounding_box):
        """
        Returns reach_ids for reaches that intersect an input bounding box

        bounding_box = [lonmin, latmin, lonmax, latmax]
        """
        lonmin, latmin, lonmax, latmax = bounding_box
        if lonmax < lonmin:
            lonmax += 360

        # iterate over reaches in self.reaches
        reach_zips = zip(
            self.x_min, self.y_min, self.x_max, self.y_max, self.reach_id)

        overlapping_reach_ids = []
        for ireach, reach_bbox in enumerate(reach_zips):

            reach_lonmin, reach_latmin, reach_lonmax, reach_latmax, reach_id =\
                reach_zips

            if reach_lonmax < reach_lonmin:
                reach_lonmax += 360

            # test for overlap (assumption is 2D decomp into two 1D intervals)
            if(latmin < reach_latmax and latmax > reach_latmin and
               lonmin < reach_lonmax and lonmax > reach_lonmin):
                overlapping_reach_ids.append(reach_id)
        return overlapping_reach_ids

class ReachDatabaseCenterlines(Product):
    """Prior Reach database centerlines"""
    ATTRIBUTES = []
    DIMENSIONS = odict([['points', 0]])
    DIMENSIONS_REACHES = 
    VARIABLES = odict([
        ['x',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS]])],
        ['y',
         odict([['dtype', 'f4'], ['dimensions', DIMENSIONS]])],
        ['reach_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['node_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ['cl_id',
         odict([['dtype', 'i4'], ['dimensions', DIMENSIONS]])],
        ])
