"""
A derived class from RiverObs, which adjusts the initial guess for the
centerline based on the observation locations.
"""

from collections import OrderedDict as odict
import numpy as N
from RiverObs import RiverObs
from Centerline import Centerline

class CenterlineObs:
    """An auxiliary class to hold values of observations associated
    with the centerline."""

    x = None # Array of x locations
    y = None # Array of y locations
    v = None # Value at (x,y)
    populated_nodes = None # Which of the centerline nodes are populated

class IteratedRiverObs(RiverObs):
    """A derived class from RiverObs, which adjusts the initial guess for the
    centerline based on the observation locations.

    The class has the same initialization as RiverObs and does not override
    any of the class functions. It adds functions to iterate the centerline.
    """

    def __init__(self,reach,xobs,yobs,k=3,ds=None,max_width=None,minobs=1):
        """Initialize with a reach variable (e.g., from ReachExtractor),
        and a set of observation coordinates.

        reach: has reach.x,reach.y (and optionally, reach.metadata).
        
        xobs, yobs: iterables with observation coordinates.
        
        k: centerline spline smoothing degree (default 3)
        
        ds: centerline point separation (default None)
        
        max_width: if !=None, exclude all observations more than max_width/2
                   away from the centerline in the normal direction. 
                   On __init__ max_width should be a number. After centerline
                   convergence, an iterable of the same
                   size as reach.x or reach.y can be added using add_centerline_obs.
                   To reflag the data with the new max_width, run reinititialize.
                   
        minobs: minimum number of observations for each node.
        """

        # Right now, things do not work for width observations along the
        # reach, due to ambiguities when the true reac is far from the original

        if N.iterable(max_width):
            raise Exception('maximum_width arrays must be added using add_maximum_width_array')

        # Initialize the base class

        RiverObs.__init__(self,reach,xobs,yobs,k=k,ds=ds,max_width=max_width,
                           minobs=minobs)

        # Unlike the base class, keep a copy of the original data and some
        # additional variables

        self.k = k
        self.ds = ds
        self.reach = reach

        self.xobs = xobs
        self.yobs = yobs

        # Now add the observed coordinates at each node

        self.add_obs('xo',self.xobs)
        self.add_obs('yo',self.yobs)
        self.load_nodes(['xo','yo'])

        # This holds centerline arrays added by add_centerline_obs

        self.centerline_obs = {}

    def refine_centerline(self,alpha=1.,stat='mean'):
        """Calculate coordinate centroids, and change positions by an attenuated step
        to avoid overshooting.

        alpha: smoothing constant. If the centerline is given by (x0,y0) and the
               data centroid is given by (xc,yc), the updated centerline coordinates
               will be (x1, y1) = (x0,y0) + alpha*( (xc,yc) - (x0,y0) ). For alpha=1,
               this is equivalent to straight replacement, but may overshoot.

        stat: statistic to use to get the centroid based on the observation coordinates.
              It should be implemented by RiverNode. Choices could be 'mean', 'median', etc.
              The default is 'mean'."""

        xc = N.asarray(self.get_node_stat(stat,'xo'))
        yc = N.asarray(self.get_node_stat(stat,'yo'))

        # The following is not efficient (due to copies). Refine later.
        
        x0 = N.asarray(self.get_node_stat(stat,'x'))
        y0 = N.asarray(self.get_node_stat(stat,'y'))

        dx = xc - x0
        dy = yc - y0

        # Compute for stopping criterion
         
        eps = max(N.abs(dx).max(), N.abs(dy).max())

        # New coordinates

        x1 = (1-alpha)*x0 + alpha*xc
        y1 = (1-alpha)*y0 + alpha*yc
            
        return x1, y1, eps

    def iterate(self,max_iter=1,alpha=1.,stat='mean',tol=1.):
        """Iterate until the coordinates change by a most tol.


        Right now, experimental evidence shows that this should be used with
        max_iter=1, alpha=1. Attempting to iterate results in a very jagged
        line since there is no penalty for large curvature. A more sophisticated
        method is required.
        """

        for i in range(max_iter):

            x1, y1, eps = self.refine_centerline(alpha=alpha,stat=stat)

            print 'iteration %d maximum coordinate change: %f'%(i,eps)

            if eps < tol:
                break

            # Recalculate the centerline for the next iteration
            
            self.reinit_centerline_nodes(x1,y1)

    def reinit_centerline_nodes(self,x1,y1):
        """Reinitialize the centerline and rebin the measurements."""

        # Calculate the centerline for this reach

        self.centerline = Centerline(x1,y1,k=self.k,ds=self.ds)

        # Calculate the local coordiantes for each observation point
        # index: the index of the nearest point
        # d: distance to the point
        # x,y: The coordiantes of the nearest point
        # s,n: The along and across track coordinates of the point
        # relative to the nearest point coordinate system.

        self.index,self.d,self.x,self.y,self.s,self.n = (
            self.centerline(self.xobs,self.yobs) )

        # Assign to each point the actual along-track distance, not just the delta s

        self.s += self.centerline.s[self.index]

        # Edit, so that only river points appear

        if self.max_width != None:
            self.in_channel = self.flag_out_channel(self.max_width)
        self.nedited_data = len(self.x)

        # Get the mapping from observation to node position (1 -> many); i.e., the inverse
        # of index (many -> 1), which maps node position to observations

        self.populated_nodes, self.obs_to_node_map = (
            self.get_obs_to_node_map(self.index,self.minobs) )

        # Now add the observed coordinates at each node

        self.add_obs('xo',self.xobs)
        self.add_obs('yo',self.yobs)
        self.load_nodes(['xo','yo'])

    def get_centerline_xy(self):
        """Return the centerline coordinates."""

        return self.centerline.x, self.centerline.y

    def get_centerline_xyv(self,name):
        """Return the centerline coordinates."""

        return ( self.centerline_obs[name].x, self.centerline_obs[name].y,
                 self.centerline_obs[name].v )

    def add_centerline_obs(self,xv,yv,v,name,minobs=1):
        """Given an array along a line x,y, interpolate it to the centerline.

        The result is appended to the class member centerline_obs,
        which is a dictionary containing CenterlineObs objects index
        by the name of the observation. Unlike RiverNode objects,
        these objects only store mean values for each populated
        node position.

        xv,yv: x,y coordinates of the vector
        v: array values at each (x,y)
        name: name of the variable to be added as a class variable
        """

        index,d,x,y,s,n = self.centerline(xv,yv)

        # Find mapping to the nodes. Note that this is
        # similar to self.get_obs_to_node_map, but it
        # does not overwrite the observation obs_to_node_map
        
        nodes = N.unique(index)

        obs_to_node_map = odict()
        populated_nodes = []
        for node in nodes:
            obs_index = N.flatnonzero(index == node)
            nobs = len(obs_index)
            if nobs >= minobs:
                populated_nodes.append(node)
                obs_to_node_map[node] = obs_index

        # Add to the observations along the centerline

        self.centerline_obs[name] = CenterlineObs()
        self.centerline_obs[name].populated_nodes = N.asarray(populated_nodes)

        self.centerline_obs[name].x = N.asarray([ N.mean(x[obs_to_node_map[node]])
                                        for node in populated_nodes ] )
        
        self.centerline_obs[name].y = N.asarray([ N.mean(y[obs_to_node_map[node]])
                                        for node in populated_nodes ])

        v = N.asarray(v)
        self.centerline_obs[name].v = N.asarray([ N.mean(v[obs_to_node_map[node]])
                                        for node in populated_nodes ])

    def reinitialize(self):
        """Reinitialize with the current centerline and max_width measurements."""

        class reach:pass
        if 'max_width' in self.centerline_obs:
            reach.x, reach.y, max_width = self.get_centerline_xyv('max_width')
            reach.metadata = self.metadata
        else:
            reach.x, reach.y  = self.get_centerline_xy()
            reach.metadata = self.metadata
            max_width = self.max_width

        # Initialize the base class

        RiverObs.__init__(self,reach,self.xobs,self.yobs,k=self.k,ds=self.ds,
                          max_width=max_width,
                            minobs=self.minobs)

        
            



        
