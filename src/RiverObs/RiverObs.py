from __future__ import absolute_import, division, print_function

import copy
import collections
import scipy.stats
import numpy as np
import logging
import warnings
import contextlib

from Centerline import Centerline
from RiverObs import RiverNode
from SWOTRiver.errors import RiverObsUseageException

LOGGER = logging.getLogger(__name__)

MISSING_VALUE_INT4 = -999
MISSING_VALUE_INT9 = -99999999
MISSING_VALUE_FLT = -999999999999.0

class RiverObs:
    """
    A class for holding all of the river observations associated with a reach.
    Observations are broken up into RiverNodes, each node associated with
    a center line point.

    The class supports extracting summary observations from each node and
    returning them for analysis (e.g., fitting).

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
                 missing_value=MISSING_VALUE_FLT,
                 use_ext_dist_coef=False):

        self.missing_value = missing_value

        # for using with a 2-pass node assignment
        self.use_ext_dist_coef = use_ext_dist_coef

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

        # Associate an along-track dimension to each node
        if ds is not None:  # Evenly spaced nodes
            self.ds = ds * np.ones(
                len(self.centerline.s), dtype=self.centerline.s.dtype)
        else:
            self.ds = reach.node_length

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

        # Assign to each point the along-track distance, not just delta s
        self.s += self.centerline.s[self.index]

        # Flag out pixels not in the dominant segmentation label
        if self.max_width is not None:

            # Use variable ext_dist_coef on second pass
            if self.use_ext_dist_coef:
                self.in_channel = self.flag_out_channel_and_label(
                    self.max_width, seg_label, ext_dist_coef=reach.ext_dist_coef)

            else:
                self.in_channel = self.flag_out_channel_and_label(
                    self.max_width, seg_label, ext_dist_coef=None)

        self.nedited_data = len(self.x)
        # Get the mapping from observation to node position (1 -> many);
        # i.e., the inverse of index (many -> 1), which maps node position
        # to observations

        self.minobs = minobs
        self.populated_nodes, self.obs_to_node_map = self.get_obs_to_node_map(
            self.index, self.minobs)

    def flag_out_channel_and_label(
        self, max_width, seg_label_in, ext_dist_coef=None):
        """
        Gets the indices of all pixels that may be inside a channel given a
        prior max_width, the node-to-node distances, an extreme distance
        threshold, and the pixel segmentation labels. It first uses the
        maximum distance threshold and node spacing to identify all pixels that
        may belong to the reach centerline. It then uses the segmentation
        labels to identify which segment is largest and calls this the
        "dominant label". This is treated as the river channel segment
        hereafter. Finally, it will "pull in" any pixels connected to the
        dominant label up to the extreme distance. It does not prune any pixels
        that are unconnected to the dominant label but within the original
        max_dist threshold.

        If no extreme distance value is given it will use a very wide threshold
        of 20x the input max_width for all nodes.

        Parameters
        ----------
        max_width : The expected maximum width of the channel. May be provided
        as an array with a width value for each node OR as one value for the
        entire reach.
        seg_label_in : An array of length [num_water_pixels_in_tile] where each
        value corresponds to a given pixel's segmentation label. These are used
        to determine which segment corresponds to the river channel.
        ext_dist_coef : Scales the maximum distance threshold based on
        proximity of the river channel to known PLD lakes. Int, typically 1-5

        Outputs
        ----------
        self.in_channel : A mask of all pixels in the tile that may belong to
        the input reach.

        """
        seg_label = seg_label_in.copy()

        # Brent Williams, May 2017: added this function to handle
        # segmentation/exclude unconnected-to-river pixels.
        # get dominant label & map centerline observable to measurements

        # get the cross-reach distance thresholds
        max_distance, extreme_dist = self.get_ext_dist_threshold(
            max_width, ext_dist_coef)
        # determine the along-reach distances for each pixel
        dst0 = abs(self.s - self.centerline.s[self.index])

        # mask pixels outside of along/cross-reach distance bounds
        node_spacing = abs(self.ds[self.index])
        self.in_channel = np.logical_and(
            abs(self.n) <= max_distance,
            dst0 <= 3.0 * node_spacing)

        # Find the dominant label and include pixels up to the extreme dist
        self.dominant_label = None
        if seg_label is not None and self.in_channel.any():
            class_mask = np.logical_and(self.in_channel, seg_label > 0)
            if class_mask.any():
                # find the largest (i.e. dominant) label within the
                # max_distance bounds
                try:
                    dominant_label = scipy.stats.mode(
                        seg_label[class_mask], keepdims=False)[0]
                except TypeError:
                    # Try previous syntax if TypeError raised (older scipy)
                    dominant_label = scipy.stats.mode(
                        seg_label[class_mask])[0][0]

                self.dominant_label = dominant_label

                # If we are in final pass of pixel assignment (where we use
                # varying widths per node), merge segmentation labels for all
                # water bodies within max_distance/2.
                if np.iterable(max_distance):
                    in_center_channel = np.logical_and(
                        self.in_channel, abs(self.n) <= max_distance/2)

                    # loop over all nodes in center channel
                    for this_index in np.unique(self.index[in_center_channel]):
                        this_node_mask = np.logical_and(
                            in_center_channel, self.index == this_index)

                        # search for segment in this node which has smallest
                        # absolute value of n coordinate
                        min_n_seg_label, min_n = -1, 9999999999999
                        for this_seg_label in np.unique(
                                seg_label[this_node_mask]):

                            sub_node_mask = np.logical_and(
                                self.index == this_index,
                                seg_label == this_seg_label)

                            this_min_n = np.abs(self.n[sub_node_mask]).min()
                            if this_min_n < min_n:
                                min_n = this_min_n
                                min_n_seg_label = this_seg_label

                        # merge seg label of segment with min abs n coordinate
                        # with dominant label
                        if min_n_seg_label != -1:
                            seg_label[seg_label == min_n_seg_label] = (
                                dominant_label)

                # keep things already in channel as well as things in dominant
                # segmentation label up to the extreme distance
                # (along and cross river)
                self.in_channel = np.logical_or(
                    self.in_channel,
                    np.logical_and(seg_label == dominant_label,
                                   np.logical_and(
                                       dst0 <= extreme_dist,
                                       abs(self.n) <= extreme_dist)))

            else:
                self.in_channel = class_mask

        self.index = self.index[self.in_channel]
        self.d = self.d[self.in_channel]
        self.x = self.x[self.in_channel]
        self.y = self.y[self.in_channel]
        self.s = self.s[self.in_channel]
        self.n = self.n[self.in_channel]
        return self.in_channel

    def get_ext_dist_threshold(self, max_width, ext_dist_coef):
        """
        Converts a maximum channel width and extreme distance coefficient into
        cross-channel distance thresholds max_distance and extreme_dist.

        Parameters
        ----------
        max_width: The expected maximum width of the channel.
        ext_dist_coef : Scalar factor for the maximum distance threshold based
        on proximity of the river channel to known PLD lakes. Int, typically 1-5

        Outputs
        ----------
        max_distance : The maximum cross-reach distance from the centerline
        that river pixels are expected to have based on prior river width
        knowledge. Assumes centerline transects the middle of the river
        channel.
        extreme_dist : The maximum cross-reach distance that connected water
        pixels are allowed to have based on known nearby waterbodies e.g. PLD
        lakes
        """
        if np.iterable(max_width):
            max_distance = max_width[self.index] / 2.
        else:
            max_distance = max_width / 2.

        node_spacing = abs(self.ds[self.index])
        if ext_dist_coef is None:
            scale_factor = 20.0
        else:
            scale_factor = ext_dist_coef[self.index]
        extreme_dist = scale_factor * np.maximum(node_spacing, max_distance)
        return max_distance, extreme_dist

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
            raise RiverObsUseageException(
                'Observation size incompatible with initial observations')

        if self.max_width is not None and len(obs) == self.ndata:
            obs = obs[self.in_channel]
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

        return obs[self.obs_to_node_map[node]]

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

    def get_node_stat(self, stat, var, all_nodes=False, **kwargs):
        """
        Get a list of results of applying a given stat to a river node
        variable.

        Both stat and var are strings. var should be the name of an
        instance variable for the river node.

        A stat is a member function of the river node which returns a
        result given the variable name. Keyword arguments can vary depending on
        which stat is being called.

        Example statfns include:
            'mean', 'median', 'std', 'cdf', 'height_weighted_mean'

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
                node_stat = getattr(river_node, stat)(var, **kwargs)
                if np.isnan(node_stat):
                    node_stat = self.missing_value
                result.append(node_stat)
            elif all_nodes:
                result.append(self.missing_value)

        return result

    def get_node_agg(
            self, height_method='weight', area_method='composite',
            all_nodes=False, goodvar_wse='good', goodvar_area='good',
            goodvar_sig0='good'):
        """
        Get lists of height, areas, and uncertainties

        If all_nodes is True, populated and unpopulated nodes are returned.
        Otherwise, only populated nodes are returned.

        The result gives arrays over desired nodes, with the populated nodes
        holding the result and the unpopulated nodes (when requested) holding
        the missing_value.
        """
        outputs = {key: [] for key in [
            'h', 'h_std', 'h_u', 'lat_u', 'lon_u', 'area', 'area_u',
            'area_det', 'area_det_u', 'width_area', 'width_area_u', 'sig0',
            'sig0_u', 'sig0_std']}

        for node in self.all_nodes:
            if node in self.populated_nodes:
                river_node = self.river_nodes[node]

                h, h_std, h_u, lat_u, lon_u = river_node.height_with_uncert(
                    method=height_method, goodvar=goodvar_wse)

                sig0, sig0_std, sig0_u = river_node.sig0_with_uncert(
                    goodvar=goodvar_sig0)

                area, width_area, area_u, width_area_u, area_det, area_det_u =\
                    river_node.area_with_uncert(
                        method=area_method, goodvar=goodvar_area)

                local_vars = locals()
                for key in outputs:
                    value = local_vars[key]
                    if value is None: value = self.missing_value
                    outputs[key].append(value)

            elif all_nodes:
                for key in outputs:
                    outputs[key].append(self.missing_value)

        # Cast to arrays with fill values instead of NaNs
        for key in outputs:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                outputs[key] = np.asarray(outputs[key])
            outputs[key][np.isnan(outputs[key])] = MISSING_VALUE_FLT
        return outputs

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

        If reverse is True, move in the opposite direction. No information is
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

    def mask_pixels(self, good_mask):
        """
        Removes pixels from the RiverObs class based on an input mask. The
        input mask should be True for pixel indices we want to keep and False
        for pixels that should be removed.
        Inputs
        :param good_mask: a boolean array of pixel indices that is True for
        pixels we want to keep, and False for pixels we want to remove.
        """

        # reshape good_mask in order to update RiverObs in_channel mask
        indx = np.argwhere(self.in_channel)[:, 0][good_mask]
        new_in_channel = np.zeros(self.in_channel.shape, dtype=bool)
        new_in_channel[indx] = True
        self.in_channel = new_in_channel

        # Recompute things set in RiverObs constructor
        self.index = self.index[good_mask]
        self.d = self.d[good_mask]
        self.x = self.x[good_mask]
        self.y = self.y[good_mask]
        self.s = self.s[good_mask]
        self.n = self.n[good_mask]
        self.nedited_data = len(self.d)
        self.populated_nodes, self.obs_to_node_map = self.get_obs_to_node_map(
            self.index, self.minobs)

        self.add_obs('xo', self.xobs)
        self.add_obs('yo', self.yobs)
        self.load_nodes(['xo', 'yo'])
