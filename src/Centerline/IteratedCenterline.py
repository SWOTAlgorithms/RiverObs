"""
Based on an initial set of points defining a centerline, and another set of
points depicting the coordinates of the water, iterate so that the centerline
is on the (possibly edited) centroid of the points. This is a class derived 
from Centerline. 
"""

import numpy as N
from Centerline import Centerline

class IteratedCenterline(Centerline):
    """Based on an initial set of points defining a centerline, and another set of
    points depicting the coordinates of the water, iterate so that the centerline
    is on the (possibly edited) centroid of the points. The number of points
    in the original Centerline is preserved."""

    def __init__(self,xcl,ycl,xdata,ydata,k=3,ds=None):
        """Initialize with an initial set of centerline coordinates and a set of data points.

        The following are all iterables:
        
        xcl: centerline xcoordinate
        ycl: centerline ycoordinate
        xdata: data x coordinates
        ydata: data y coordinates

        k: order of the spline to use for interpolation (1<= k <=5)
        ds: if not None, resample x,y points to this spacing 
        (approximately to not extend the interval).
        """

        # Initialize the base class and store the data points
        
        Centerline(xcl,ycl,k=k)

        # Store the data points
        
        self.xdata = N.asarray(xdata)
        self.ydata = N.asarray(ydata)

    def iterate(self,max_iter=100,max_distance=None,eps=1.e-6):
        """Iterate the centerline until convergence is reached or the 
        maximum number of iterations are exceeded."""

        pass
