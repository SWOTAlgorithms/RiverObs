"""A class for computing the location of a point or set of points
relative to a curved line defined by a series of two dimentional
points."""

import numpy as N
from scipy.spatial import cKDTree
from scipy.interpolate import splrep, splev

class Centerline:
    """A class for computing the location of a point or set of points
    relative to a curved line defined by a series of two dimentional
    points."""

    def __init__(self,x,y,k=3,ds=None,eps=1.e-6):
        """Initialize with a set of x,y points. A Euclidean
        metric is assumed for distance calculation. Optionally,
        the centerline is resampled to a uniform along-track
        spacing.

        x,y: iterables of the same dimension.
        
        k: order of the spline to use for interpolation (1<= k <=5)
        ds: if not None, resample x,y points to this spacing 
        (approximately to not extend the interval).
        eps: minimum distance for two points to be considered different
        """

        # The following holds the array for cKDTree, which
        # must not be touched during calls

        self.x = N.asarray(x)
        self.y = N.asarray(y)

        # Compute the point separation along the curve

        self.delta = N.zeros((len(x),),dtype=N.float64)
        
        self.delta[1:] = N.sqrt((self.x[1:]-self.x[:-1])**2 + 
                                (self.y[1:]-self.y[:-1])**2 )

        # Remove any duplicate points (point separation < eps)

        self.delta[0] = 1 # to avoid removing the first point
        unique = self.delta > eps
        self.delta[0] = 0.
        self.x = self.x[unique]
        self.y = self.y[unique]
        self.delta = self.delta[unique]

        # Compute the distance along the curve
                
        self.s = N.cumsum(self.delta)

        # Compute the spline for evaluating x and y as a function of s

        self.xtck = splrep(self.s,self.x,k=k)
        self.ytck = splrep(self.s,self.y,k=k)

        # If resampling is desired, find the new s, x, and y

        if ds != None:
            ns = int(self.s[-1]/ds + 1)
            if ns < 2:
                raise Exception('This ds is too large for the data set')

            self.ds =  self.s[-1]/(ns - 1)
            self.s = N.arange(ns)*self.ds
            self.x = splev(self.s,self.xtck)
            self.y = splev(self.s,self.ytck)

        # Initialize the cKDtree

        self.xy = N.zeros((len(self.x),2),dtype=N.float64)
        self.xy[:,0] = self.x
        self.xy[:,1] = self.y

        self.kdtree = cKDTree(self.xy)
                    
        # Calculate the tangent and normal vectors at each point

        self.dx_ds = splev(self.s,self.xtck,der=1)
        self.dy_ds = splev(self.s,self.ytck,der=1)

        self.tangent = N.zeros((len(self.x),2),dtype=N.float64)
        self.tangent[:,0] = self.dx_ds
        self.tangent[:,1] = self.dy_ds
        self.tangent /= N.sqrt(self.tangent[:,0]**2 + self.tangent[:,1]**2)[:,N.newaxis]

        self.normal = N.zeros((len(self.x),2),dtype=N.float64)
        self.normal[:,0] = -self.tangent[:,1]
        self.normal[:,1] = self.tangent[:,0]

    def __call__(self,x0,y0):
        """For each point in (x0,y0), return the nearest point, as well
        as the along and across track coordinates for that point in the
        local coordinate system as the point.

        x0, y0: 1D iterables of the same dimension which can be cast to 
        numpy 1D arrays.

        Returns:

        index: the index of the nearest point
        x,y: The coordiantes of the nearest point
        s,n: The along and across track coordinates of the point
        relative to the nearest point coordinate system.
        """

        xy = N.column_stack([x0,y0])
        d, i = self.kdtree.query(xy)
        x = self.x[i]
        y = self.y[i]
        tx = self.tangent[i][:,0]
        ty = self.tangent[i][:,1]
        s = (x0 - x)*tx + (y0 - y)*ty
        n = -(x0 - x)*ty + (y0 - y)*tx

        return i,d,x,y,s,n
        
