"""
A derived class from RiverObs, which adjusts the initial guess for the
centerline based on the observation locations.
"""
from __future__ import absolute_import, division, print_function

import collections
import numpy as np
import scipy.interpolate

from RiverObs.RiverObs import RiverObs
from Centerline import Centerline
from RiverObs import RiverNode
from SWOTRiver.errors import RiverObsUseageException

class CenterlineObs:
    """
    An auxiliary class to hold values of observations associated
    with the centerline.
    """
    x = None  # Array of x locations
    y = None  # Array of y locations
    v = None  # Value at (x,y)
    populated_nodes = None  # Which of the centerline nodes are populated


class IteratedRiverObs(RiverObs):
    """
    A derived class from RiverObs, which adjusts the initial guess for the
    centerline based on the observation locations.

    The class has the same initialization as RiverObs and does not override
    any of the class functions. It adds functions to iterate the centerline.

    Parameters
    ----------

    reach : object
        Has reach.x,reach.y (and optionally, reach.metadata).
    xobs, yobs : array_like
        Iterables with observation coordinates.
    k : int
        Centerline spline smoothing degree (default 3).
    ds : float
        Centerline point separation (default None).
    seg_label : array_like
        Array giving a different index to each feature (connected region)
    max_width : float
        If !=None, exclude all observations more than max_width/2
        away from the centerline in the normal direction.
        On __init__ max_width should be a number. After centerline
        convergence, an iterable of the same
        size as reach.x or reach.y can be added using add_centerline_obs.
        To reflag the data with the new max_width, run reinititialize.
    minobs : int
        Minimum number of observations for each node.
    **kwds :
        Additional keywords to pass to RiverObs.__init__
    """

    def __init__(self,
                 reach,
                 xobs,
                 yobs,
                 k=3,
                 ds=None,
                 seg_label=None,
                 max_width=None,
                 minobs=1,
                 **kwds):
        # Right now, things do not work for width observations along the
        # reach due to ambiguities when the true reach is far from the original
        if np.iterable(max_width):
            raise RiverObsUseageException(
                'maximum_width arrays must be added using add_maximum_width_array'
            )

        # Initialize the base class
        self.robs_kwds = kwds
        RiverObs.__init__(
            self,
            reach,
            xobs,
            yobs,
            k=k,
            ds=ds,
            seg_label=seg_label,
            max_width=max_width,
            minobs=minobs,
            **kwds)

        # Unlike the base class, keep a copy of the original data and some
        # additional variables
        self.k = k
        self.ds_init = ds
        self.reach = reach
        self.xobs = xobs
        self.yobs = yobs
        self.seg_label = seg_label

        # Now add the observed coordinates at each node
        self.add_obs('xo', self.xobs)
        self.add_obs('yo', self.yobs)
        self.load_nodes(['xo', 'yo'])

        # This holds centerline arrays added by add_centerline_obs
        self.centerline_obs = {}

    def refine_centerline(self, alpha=1., statfn='mean', std_statfn='std'):
        """
        Calculate coordinate centroids, and change positions by an attenuated
        step to avoid overshooting.

        Parameters
        ----------
        alpha : float
            Smoothing constant. If the centerline is given by (x0,y0) and
            the data centroid is given by (xc,yc), the updated centerline
            coordinates will be (x1, y1) = (x0,y0) + alpha*((xc,yc) - (x0,y0)).
            For alpha=1, this is equivalent to straight replacement, but may
            overshoot.
        statfn : str
            Statistic to use to get the centroid based on the observation
            coordinates. It should be implemented by RiverNode. Choices could
            be 'mean', 'median', etc. The default is 'mean'.
        std_statfn : str
            Statistic to use to get the centroid dispersion based on the
            observation coordinates. It should be implemented by RiverNode.
            Choices could be 'std', 'stderr', etc. The default is 'std'.

        Returns
        -------
        x1, y1 : array_like
            New centerline coordinate proposals.
        eps : float
            Maximum coordinate change in any direction.
        xdelta, ydelta : array_like
            Estimates of data point dispersion, to be used for
            weighting the centerline fit.
        """
        xc = np.asarray(self.get_node_stat(statfn, 'xo'))
        yc = np.asarray(self.get_node_stat(statfn, 'yo'))

        # Compute the scatter among the points for spline weighting
        xdelta = np.asarray(self.get_node_stat(std_statfn, 'xo'))
        ydelta = np.asarray(self.get_node_stat(std_statfn, 'yo'))

        # The following is not efficient (due to copies). Refine later.
        x0 = np.asarray(self.get_node_stat(statfn, 'x'))
        y0 = np.asarray(self.get_node_stat(statfn, 'y'))
        #print "xc:",xc
        #print "x:",x0
        dx = xc - x0
        dy = yc - y0

        # Compute for stopping criterion
        eps = max(np.abs(dx).max(), np.abs(dy).max())

        # New coordinates
        x1 = (1 - alpha) * x0 + alpha * xc
        y1 = (1 - alpha) * y0 + alpha * yc

        return x1, y1, eps, xdelta, ydelta

    def iterate(self,
                max_iter=1,
                alpha=1.,
                statfn='mean',
                std_statfn='std',
                tol=1.,
                weights=True,
                smooth=1.e-2,
                **kwds):
        """
        Iterate until the coordinates change by a most tol.

        Parameters
        ----------

        max_iter : int
            Maximum number of iterations
        alpha : float <= 1
            Smoothing constant. If the centerline is given by (x0,y0) 
            and the data centroid is given by (xc,yc), the updated
            centerline coordinates will be:
            (x1, y1) = (x0,y0) + alpha*( (xc,yc) - (x0,y0) ).
            For alpha=1, this is equivalent to straight replacement, but
            may overshoot.
        statfn : str
            Statistic to use to get the centroid based on the observation
            coordinates.  It should be implemented by RiverNode. Choices
            could be 'mean', 'median', etc. The default is 'mean'.
        std_statfn : str
            Statistic to use to get the centroid dispersion based on the
            observation coordinates. It should be implemented by RiverNode.
            Choices could be 'std', 'stderr', etc. The default is 'std'.
        tol : float
            Maximum coordinate change to stop iteration.
        weights : bool
            If True, compute weights from point dispersion to pass to 
            Centerline.
        smooth : float
            Smoothing parameter to pass to Centerline, if weight = True.
            This is the related to the splrep s parameter by s = smooth*len(x),
            so that smooth = 1. would be the appropriate weighting when
            w = 1/std.
        **kwds :
            Additional keywords o pass to Centerline

        Notes
        -----

        Right now, experimental evidence shows that this should be used with
        max_iter=1 or 2, alpha=1. smooth should be between 1.e-2 and 1.e-1,
        max_iter = 1 when few points, max_iter = 2 for greater point density.
        Attempting to iterate until convergence results in a very jagged
        for small point density.
        """
        tiny = 1.e-6  # avoid overflows

        for i in range(max_iter):
            x1, y1, eps, xdelta, ydelta = self.refine_centerline(
                alpha=alpha, statfn=statfn, std_statfn=std_statfn)

            print('iteration %d maximum coordinate change: %f' % (i, eps))

            if eps < tol:
                break

            # If weights are desired, calculate based on the estimated standard
            # deviation. For points with std = 0, use the maximum non-zero std
            if weights:
                xstdmax = xdelta.max()
                if xstdmax == 0:
                    xstdmax = 1.

                wx = np.where(xdelta > 0, 1. / (xdelta + tiny), 1. / xstdmax)

                ystdmax = ydelta.max()
                if ystdmax == 0:
                    ystdmax = 1.

                wy = np.where(ydelta > 0, 1. / (ydelta + tiny), 1. / ystdmax)

            else:
                wx, wy = None, None

            # Recalculate the centerline for the next iteration
            self.reinit_centerline_nodes(x1, y1, smooth=smooth, wx=wx, wy=wy)

    def reinit_centerline_nodes(self, x1, y1, smooth=None, wx=None, wy=None):
        """Reinitialize the centerline and re-bin the measurements."""

        # If weights and smoothing are given, estimate the centerline points
        # from a smoothed spline
        if smooth is not None and wx is not None and wy is not None:
            centerline = Centerline(
                x1, y1, k=self.k, ds=self.ds_init, smooth=smooth, wx=wx, wy=wy)
            x1 = scipy.interpolate.splev(centerline.s, centerline.xtck)
            y1 = scipy.interpolate.splev(centerline.s, centerline.ytck)

        # Calculate the centerline for this reach
        self.centerline = Centerline(x1, y1, k=self.k, ds=self.ds_init)

        # Associate an along-track dimension to each node
        if self.ds_init is not None:  # Evenly spaced nodes
            self.ds = np.ones(
                len(self.centerline.s),
                dtype=self.centerline.s.dtype) * self.ds_init

        else:
            self.ds = np.ones(
                len(self.centerline.s), dtype=self.centerline.s.dtype)
            self.ds[1:-1] = (
                self.centerline.s[2:] - self.centerline.s[0:-2]) / 2.
            self.ds[0] = self.ds[1]
            self.ds[-1] = self.ds[-2]

        # Calculate the local coordinates for each observation point
        # index: the index of the nearest point
        # d: distance to the point
        # x,y: The coordinates of the nearest point
        # s,n: The along and across track coordinates of the point
        # relative to the nearest point coordinate system.
        self.index, self.d, self.x, self.y, self.s, self.n = (self.centerline(
            self.xobs, self.yobs))

        # Assign to each point the actual along-track distance, not just
        # the delta s
        self.s += self.centerline.s[self.index]

        # Edit, so that only river points appear
        if self.max_width != None:
            self.in_channel = self.flag_out_channel(self.max_width)

        self.nedited_data = len(self.x)

        # Get the mapping from observation to node position (1 -> many);
        # i.e., the inverse of index (many -> 1), which maps node position
        # to observations
        self.populated_nodes, self.obs_to_node_map = (self.get_obs_to_node_map(
            self.index, self.minobs))

        # Now add the observed coordinates at each node
        self.add_obs('xo', self.xobs)
        self.add_obs('yo', self.yobs)
        self.load_nodes(['xo', 'yo'])

    def get_centerline_xy(self):
        """Return the centerline coordinates."""
        return self.centerline.x, self.centerline.y

    def get_centerline_xyv(self, name):
        """Return the centerline coordinates."""
        return (self.centerline_obs[name].x, self.centerline_obs[name].y,
                self.centerline_obs[name].v)

    def get_centerline_xy_max_width(self):
        """Return the centerline coordinates."""
        return self.centerline.x, self.centerline.y, self.max_width

    def add_centerline_obs(self, xv, yv, v, name, minobs=1):
        """
        Given an array along a line x,y, interpolate it to the centerline.

        The result is appended to the class member centerline_obs,
        which is a dictionary containing CenterlineObs objects index
        by the name of the observation. Unlike RiverNode objects,
        these objects only store mean values for each populated
        node position.

        Parameters
        ----------
        xv, yv : array_like
                 x,y coordinates of the vector
        v : array_like
            array values at each (x,y)
        name : str
               name of the variable to be added as a class variable
        """

        index, d, x, y, s, n = self.centerline(xv, yv)

        # Find mapping to the nodes. Note that this is
        # similar to self.get_obs_to_node_map, but it
        # does not overwrite the observation obs_to_node_map
        nodes = np.unique(index)

        obs_to_node_map = collections.OrderedDict()
        populated_nodes = []
        for node in nodes:
            obs_index = np.flatnonzero(index == node)
            nobs = len(obs_index)
            if nobs >= minobs:
                populated_nodes.append(node)
                obs_to_node_map[node] = obs_index

        # Add to the observations along the centerline
        self.centerline_obs[name] = CenterlineObs()
        self.centerline_obs[name].populated_nodes = np.asarray(populated_nodes)

        self.centerline_obs[name].x = np.asarray(
            [np.mean(x[obs_to_node_map[node]]) for node in populated_nodes])

        self.centerline_obs[name].y = np.asarray(
            [np.mean(y[obs_to_node_map[node]]) for node in populated_nodes])

        v = np.asarray(v)
        self.centerline_obs[name].v = np.asarray(
            [np.mean(v[obs_to_node_map[node]]) for node in populated_nodes])

    def reinitialize(self):
        """
        Reinitialize with the current centerline and max_width measurements.
        """

        class reach:
            pass

        if 'max_width' in self.centerline_obs:
            reach.x, reach.y, max_width = self.get_centerline_xyv('max_width')
            reach.metadata = self.metadata

        else:
            reach.x, reach.y = self.get_centerline_xy()
            reach.metadata = self.metadata
            max_width = self.max_width

        reach.node_length = self.reach.node_length
        reach.ext_dist_coef = self.reach.ext_dist_coef

        # Initialize the base class
        RiverObs.__init__(
            self,
            reach,
            self.xobs,
            self.yobs,
            k=self.k,
            ds=self.ds_init,
            seg_label=self.seg_label,
            max_width=max_width,
            minobs=self.minobs,
            **self.robs_kwds)

        # The obs will have to be reprojected
        self.centerline_obs = {}
