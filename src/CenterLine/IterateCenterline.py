"""
Based on an initial set of points defining a centerline, and another set of
points depicting the coordinates of the water, iterate so that the centerline
is on the (possibly edited) centroid of the points. The number of points
in the original Centerline is preserved."""

import numpy as N
from Centerline import Centerline

class IterateCenterline:
    """Based on an initial set of points defining a centerline, and another set of
    points depicting the coordinates of the water, iterate so that the centerline
    is on the (possibly edited) centroid of the points. The number of points
    in the original Centerline is preserved."""

    def __init__(self,xcl,ycl,xdata,ydata,k=3):
        """Initialize with an initial set of centerline coordinates and a set of data points.

        The following are all iterables:
        
        xcl: centerline xcoordinate
        ycl: centerline ycoordinate
        xdata: data x coordinates
        """

        # This is the initial centerline to iterate
        
        self.centerline = Centerline(xcl,ycl,k=k)
        self.xdata = N.asarray(xdata)
        self.ydata = N.asarray(ydata)

    def __call__(self,max_iter=100,max_distance=None):
