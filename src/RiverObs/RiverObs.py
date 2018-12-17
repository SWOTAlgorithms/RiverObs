from __future__ import absolute_import, division, print_function

import copy
import collections
import scipy.stats
import numpy as np

from Centerline import Centerline
from .RiverNode import RiverNode


class RiverObs:
    """
    A class for holding all of the river observations associated with a reach.
    Observations are broken up into RiverNodes, each node associated with
    a center line point.

    The class supports extracting summary observations from each node and
    returning them for analaysis (e.g., fitting).

    Initialize with a reach variable (e.g., from ReachExtractor),
    and a set of observation coordinates.

    Parameters
    ----------

    reach : object
        has reach.x,reach.y (and optionally, reach.metadata).
    xobs, yobs : iterable
        iterables with observation coordinates.
    k : int
        centerline spline smoothing degree (default 3)
    ds : float
        centerline point separation (default None)
    max_width :
        if !=None, exclude all observations more than max_width/2
        away from the centerline in the normal directionp.
        max_width can be a number or an iterable of the same
        size as reach.x or reach.y. If it is an interable,
        it is added to the centerline as a member.
    minobs : int
        minimum number of observations for each node.
    node_class: class
        either RiverNode, or a class derived from it.
    missing_value : float, default -9999
        This value is reported when a node_stat is requested of an empty node.
    verbose : bool, default False
        Output progress to stdout
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
                 node_class=RiverNode,
                 missing_value=-9999,
                 verbose=False):
        self.verbose = verbose
        self.missing_value = missing_value

        # Register the node class
        self.node_class = node_class

        # Copy metadata, in case it is present
        try:
            self.metadata = reach.metadata
        except AttributeError:
            self.metadata = None

        self.ndata = len(xobs)
        # Calculate the centerline for this reach

        if max_width is None or not np.iterable(max_width):
            self.centerline = Centerline(reach.x, reach.y, k=k, ds=ds)
            self.centerline.max_width = max_width
        else:
            self.centerline = Centerline(
                reach.x,
                reach.y,
                k=k,
                ds=ds,
                obs=[max_width],
                obs_names=['max_width'])
        self.max_width = self.centerline.max_width
        if self.verbose: print('Centerline initialized')

        # Associate an along-track dimension to each node
        if ds is not None:  # Evenly spaced nodes
            self.ds = ds * np.ones(
                len(self.centerline.s), dtype=self.centerline.s.dtype)
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
        # s,n: The along and across river coordinates of the point
        # relative to the nearest point coordinate system.
        self.index, self.d, self.x, self.y, self.s, self.n = self.centerline(
            xobs, yobs)
        # squeeze extra dimensions
        self.index = np.squeeze(self.index)
        self.d = np.squeeze(self.d)
        self.x = np.squeeze(self.x)
        self.y = np.squeeze(self.y)
        self.s = np.squeeze(self.s)
        self.n = np.squeeze(self.n)

        if self.verbose: print('Local coordiantes calculated')

        # Assign to each point the along-track distance, not just delta s
        self.s += self.centerline.s[self.index]

        # Flag out pixels not in the dominant segmentation label
        if self.max_width is not None:
            self.in_channel = self.flag_out_channel_and_label(
                self.max_width, seg_label)

        self.nedited_data = len(self.x)
        print("num nodes in reach %d" % len(np.unique(self.index)))
        # Get the mapping from observation to node position (1 -> many);
        # i.e., the inverse of index (many -> 1), which maps node position
        # to observations

        self.minobs = minobs
        self.populated_nodes, self.obs_to_node_map = self.get_obs_to_node_map(
            self.index, self.minobs)

    def flag_out_channel_and_label(self, max_width, seg_label):
        """
        Gets the indexes of all of the points inside a channel of
        max_width, a segmentation label
        and remove the points from the list of observations.
        """
        # Brent Williams, May 2017: added this function to handle
        # segmentation/exclude unconnected-to-river pixels.
        # get dominant label & map centerline observalble to measurements
        if np.iterable(max_width):
            max_distance = max_width[self.index] / 2.
        else:
            max_distance = max_width / 2.

        dst0 = abs(self.s - self.centerline.s[self.index])

        extreme_dist = 20.0 * np.maximum(
            abs(self.ds[self.index]), max_distance)
        self.in_channel = np.logical_and(
            abs(self.n) <= max_distance,
            dst0 <= 3.0 * abs(self.ds[self.index]))

        # apply seg labels
        if seg_label is not None and self.in_channel.any():
            class_mask = np.logical_and(self.in_channel, seg_label > 0)
            if class_mask.any():
                dominant_label = scipy.stats.mode(seg_label[class_mask])[0][0]
                # keep things already in channel as well as things in dominant
                # segmentation label up to the extreme distance
                # (along and cross river)
                self.in_channel = np.logical_or(
                    self.in_channel,
                    np.logical_and(seg_label == dominant_label,
                                   np.logical_and(
                                       dst0 <= extreme_dist,
                                       abs(self.n) <= extreme_dist)))

                if self.verbose:
                    print("Dominant label in reach: %d" % dominant_label)

            else:
                self.in_channel = class_mask
                print("No valid class labels in reach")

        self.index = self.index[self.in_channel]
        self.d = self.d[self.in_channel]
        self.x = self.x[self.in_channel]
        self.y = self.y[self.in_channel]
        self.s = self.s[self.in_channel]
        self.n = self.n[self.in_channel]
        return self.in_channel

    def flag_out_channel(self, max_width):
        """
        Get the indexes of all of the points inside a channel of max_width,
        and remove the points from the list of observations.
        """
        if np.iterable(max_width):
            max_distance = max_width[self.index] / 2.
        else:
            max_distance = max_width / 2.

        self.in_channel = np.abs(self.n) <= max_distance

        self.index = self.index[self.in_channel]
        self.d = self.d[self.in_channel]
        self.x = self.x[self.in_channel]
        self.y = self.y[self.in_channel]
        self.s = self.s[self.in_channel]
        self.n = self.n[self.in_channel]
        return self.in_channel

    def get_obs_to_node_map(self, index, minobs=1):
        """
        Get the mapping from observation to node position (1 -> many);
        i.e., the inverse of index (many -> 1), which maps node position
        to observations.

        In order for a node to appear, it must have at least minobs
        observations.
        """

        # Get the list of potential nodes
        nodes = np.unique(index)

        self.obs_to_node_map = collections.OrderedDict()
        self.nobs = np.zeros(len(self.centerline.x), dtype=np.int32)
        self.populated_nodes = []
        for node in nodes:
            obs_index = np.flatnonzero(index == node)
            nobs = len(obs_index)
            if nobs >= minobs:
                self.populated_nodes.append(node)
                self.obs_to_node_map[node] = obs_index
                self.nobs[node] = nobs
        self.n_populated_nodes = len(self.populated_nodes)

        # Store also a list of all the potential nodes and all the
        # unpopulated nodes
        self.n_nodes = len(self.centerline.s)
        self.all_nodes = np.arange(self.n_nodes, dtype=np.int32)
        self.unpopulated_nodes = []
        for node in self.all_nodes:
            if not node in self.populated_nodes:
                self.unpopulated_nodes.append(node)
        self.n_unpopulated_nodes = len(self.unpopulated_nodes)
        return self.populated_nodes, self.obs_to_node_map

    def add_obs(self, obs_name, obs):
        """
        Add an observation as a class variable self.obs_name.

        The observation is edited to remove measurements outside
        the channel.

        obs is an iterable of length self.ndata or self.nedited_data.
        """

        if len(obs) != self.ndata and len(obs) != self.nedited_data:
            raise Exception(
                'Observation size incompatible with initial observations')

        if self.max_width is not None and len(obs) == self.ndata:
            #obs = obs[self.in_channel]
            obs = np.asarray(obs)[self.in_channel]
        setattr(self, obs_name, obs)

    def obs_to_node(self, obs, node):
        """
        Get all of the observations in an array obs which map to a node.

        Parameters
        ----------
        obs : iterable
            iterable of the same size as the xobs, yobs
            or the same size as self.x, self.y. If the same size
            as xobs, the observations will be limited to in channel
            observations, if this has been computed. If the same
            size as self.x, no editing occurs.

        node : int
            node to match

        Returns
        -------
        The observations for that node, or an empty array if there
        are no observations for that node.
        """

        if not (int(node) in self.populated_nodes):
            return np.array([])

        # If only certain observations have been kept, get the edited vector
        if self.max_width is not None and len(obs) == self.ndata:
            obs = obs[self.in_channel]

        return np.asarray(obs)[self.obs_to_node_map[node]]

    def load_nodes(self, vars=[]):
        """Load the desired variables into each of the populated nodes.

        All of the vars should have been loaded previously with add_obs.
        """

        if type(vars) == str:
            vars = [vars]

        self.river_nodes = collections.OrderedDict()

        for node in self.populated_nodes:
            d = self.obs_to_node(self.d, node)
            x = self.obs_to_node(self.x, node)
            y = self.obs_to_node(self.y, node)
            s = self.obs_to_node(self.s, node)
            n = self.obs_to_node(self.n, node)
            #h_flg = self.obs_to_node(self.h_flg,node)
            self.river_nodes[node] = self.node_class(
                node, d, x, y, s, n, ds=self.ds[node])

            for var in vars:
                obs = self.obs_to_node(getattr(self, var), node)
                self.river_nodes[node].add_obs(var, obs, sort=False)

    def get_node_stat(self, stat, var, all_nodes=False, good_flag=None):
        """
        Get a list of results of applying a given stat to a river node
        variable.

        Both stat and var are strings. var should be the name of an
        instance variable for the river node.

        A stat is a member function of the river node which returns a
        result given the variable name.

        Example statfns are: 'mean', 'std', 'cdf'

        If all_nodes is True, populated and unpopulated nodes are returned.
        Otherwise, only populated nodes are returned.

        The result is a list over desired nodes, with the populated nodes
        holding the result and the unpopulated nodes (when requested) holding
        the missing_value.
        """
        result = []
        for node in self.all_nodes:
            if node in self.populated_nodes:
                river_node = self.river_nodes[node]
                if good_flag is None:
                    result.append(getattr(river_node, stat)(var))

                else:
                    result.append(getattr(river_node, stat)(var, good_flag))

            elif all_nodes:
                result.append(self.missing_value)

        return result
    def get_node_agg(self, 
                     height_method='weight', area_method='composite',
                     all_nodes=False, good_flag='good'):
        """
        Get lists of height, areas, and uncertainties

        If all_nodes is True, populated and unpopulated nodes are returned.
        Otherwise, only populated nodes are returned.

        The result gives arrays over desired nodes, with the populated nodes
        holding the result and the unpopulated nodes (when requested) holding
        the missing_value.
        """
        h_list = []
        h_std_list = []
        h_uncert_list = []
        a_list = []
        w_a_list = []
        a_uncert_list = []
        for node in self.all_nodes:
            if node in self.populated_nodes:
                river_node = self.river_nodes[node]
                h, h_std, h_unc = \
                    river_node.aggregate_height_with_uncert(
                        method=height_method, goodvar=good_flag)
    
                h_list.append(h if h is not None else self.missing_value)
                h_std_list.append(h_std if h_std is not None else self.missing_value)
                h_uncert_list.append(h_unc if h_unc is not None else self.missing_value)
                a, w_a, a_unc = \
                    river_node.aggregate_area_with_uncert(method=area_method)

                a_list.append(a if a is not None else self.missing_value)
                w_a_list.append(w_a if w_a is not None else self.missing_value)
                a_uncert_list.append(a_unc if a_unc is not None else self.missing_value)
            elif all_nodes:
                h_list.append(self.missing_value)
                h_std_list.append(self.missing_value)
                h_uncert_list.append(self.missing_value)
                a_list.append(self.missing_value)
                w_a_list.append(self.missing_value)
                a_uncert_list.append(self.missing_value)
        # cast to arrays to make life easier later
        h = np.asarray(h_list)
        h_std = np.asarray(h_std_list)
        h_uncert = np.asarray(h_uncert_list)
        a = np.asarray(a_list)
        w_a = np.asarray(w_a_list)
        a_uncert = np.asarray(a_uncert_list)
        return h, h_std, h_uncert, a, w_a, a_uncert

    def trim_nodes(self, fraction, mode='both', sort_variable='n'):
        """
        Trim the data in all the nodes.

        fraction: 0 < f < 1. Fraction of the data to remove.
        mode is 'both', 'low', 'high' for which tails of the distribution
        need to be trimmed.

        Prior to trimming, the data are sorted according to the sort variable.
        """
        for node, river_node in self.river_nodes.items():
            river_node.sort(sort_variable=sort_variable)
            river_node.trim(fraction, mode=mode)

    def remove_nodes(self, node_list, reverse=False):
        """
        Move nodes from the populated node list to the unpopulated node list.

        If reverse is True, move in the opposite directionp. No information is
        lost during this process and it is invertible. Both lists are kept
        sorted at each step.
        """
        if not reverse:
            from_list = copy.copy(self.populated_nodes)
            to_list = copy.copy(self.unpopulated_nodes)
        else:
            to_list = copy.copy(self.populated_nodes)
            from_list = copy.copy(self.unpopulated_nodes)

        for node in node_list:
            try:
                index = from_list.index(node)
                from_list.pop(index)
                to_list.append(node)
            except ValueError:
                pass

        from_list.sort()
        to_list.sort()
        if not reverse:
            self.populated_nodes = from_list
            self.unpopulated_nodes = to_list
        else:
            self.populated_nodes = to_list
            self.unpopulated_nodes = from_list
