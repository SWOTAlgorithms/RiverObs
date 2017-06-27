"""
A class for holding all of the river observations associated with a reach.
Observations are broken up into RiverNodes, each node associated with
a center line point.

The class supports extracting summary observations from each node and
returning them for analaysis (e.g., fitting).
"""

from copy import copy
from collections import OrderedDict as odict
import numpy as N
from Centerline import Centerline
from RiverNode import RiverNode
from scipy.stats import mode

class RiverObs:
    """A class for holding all of the river observations associated with a reach.
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
        away from the centerline in the normal direction. 
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
    
    def __init__(self,reach,xobs,yobs,k=3,ds=None,seg_label=None,max_width=None,minobs=1,
                 node_class=RiverNode,missing_value=-9999,verbose=False):
        
        self.verbose = verbose
        self.missing_value = missing_value

        # Register the node class

        self.node_class = node_class

        # Copy metadata, in case it is present
        
        try:
            self.metadata = reach.metadata
        except:
            self.metadata = None

        self.ndata = len(xobs)
        # Calculate the centerline for this reach

        if type(max_width) == type(None) or not N.iterable(max_width):
            self.centerline = Centerline(reach.x,reach.y,k=k,ds=ds)
            self.centerline.max_width = max_width
        else:
            self.centerline = Centerline(reach.x,reach.y,k=k,ds=ds,
                                         obs=[max_width],obs_names=['max_width'])
        self.max_width = self.centerline.max_width
        if self.verbose: print('Centerline initialized')

        # Associate an along-track dimension to each node

        if ds != None: # Evenly spaced nodes
            self.ds = N.ones(len(self.centerline.s),dtype=self.centerline.s.dtype)*ds
        else:
            self.ds = N.ones(len(self.centerline.s),dtype=self.centerline.s.dtype)
            self.ds[1:-1] = (self.centerline.s[2:] - self.centerline.s[0:-2])/2.
            self.ds[0] = self.ds[1]
            self.ds[-1] = self.ds[-2]
        
        # Calculate the local coordiantes for each observation point
        # index: the index of the nearest point
        # d: distance to the point
        # x,y: The coordiantes of the nearest point
        # s,n: The along and across river coordinates of the point
        # relative to the nearest point coordinate system.
        
        self.index,self.d,self.x,self.y,self.s,self.n = self.centerline(xobs,yobs)
        if self.verbose: print('Local coordiantes calculated')

        # dst0 is so that we dont go farther than a node-length away in s (along river) when assigning to nodes
        # for some reason the nodes on hte ends were being located bad because they were accumulating pixels too far away in s
        dst0 = abs(self.s)-abs(self.ds[self.index]) 
        
        # Assign to each point the actual along-track distance, not just the delta s
        
        self.s += self.centerline.s[self.index]
        
        # Edit, so that only river points appear
        """
        if type(self.max_width) != type(None):
            self.in_channel = self.flag_out_channel(self.max_width)
        self.nedited_data = len(self.x)
        #print "got here,",self.in_channel
        
        # edit so only pixels of one segmentation label for entire reach is used for all nodes
        if type(seg_label) != type(None):
            self.in_label = self.flag_out_label(seg_label[self.in_channel])
            #self.in_channel = self.in_label
        self.nedited_data = len(self.x)
        """
        # added Brent Williams May 2017 to flag out pixels not in the dominant segmentation label
        if (type(self.max_width) != type(None))&(type(seg_label) != type(None)):
            self.in_channel = self.flag_out_channel_and_label(self.max_width,seg_label,dst0)
        elif type(self.max_width) != type(None):
            self.in_channel = self.flag_out_channel(self.max_width)
        self.nedited_data = len(self.x)
        print("num nodes in reach %d"%len(N.unique(self.index)))
        # Get the mapping from observation to node position (1 -> many); i.e., the inverse
        # of index (many -> 1), which maps node position to observations

        self.minobs = minobs
        self.populated_nodes, self.obs_to_node_map = self.get_obs_to_node_map(self.index,self.minobs)

    def flag_out_channel_and_label(self,max_width,seg_label,dst0):
        """Get the indexes of all of the points inside a channel of
        max_width and a segmentation label
        and remove the points from the list of observations."""
        # Brent Williams, May 2017: added this function to handle segmentation/exclude unconnected-to-river pixels
        # get dominant label
        if N.iterable(max_width): # Map centerline observalble to measurements
            max_distance = max_width[self.index]/2.
        else:
            max_distance = max_width/2.
        #
        
        msk = (N.abs(self.n) <= max_distance)&(seg_label>0)&(dst0<=0)
        
        #print "seg_lbl",seg_label[msk]
        dominant_label=mode(seg_label[msk])[0][0]
        print("DOMINANT LABEL in reach: %d"%dominant_label)
        """
        #get dominant label on a node level?
        ui=N.unique(self.index)
        print "ui",ui
        
        dominant_label=N.zeros(N.shape(seg_label))-1
        for ind in ui:
            tmp=seg_label[self.index==ind]
            #print "mode",mode(tmp[tmp>0])[0][0]
            dominant_label[self.index==ind]=mode(tmp[tmp>0])[0][0]
        #
        
        print "dominant label in node :",dominant_label
        """
        
        self.in_channel = (N.abs(self.n) <= max_distance)&(seg_label == dominant_label)&(dst0<=0)
        
        self.index = self.index[self.in_channel]
        self.d = self.d[self.in_channel]
        self.x = self.x[self.in_channel]
        self.y = self.y[self.in_channel]
        self.s = self.s[self.in_channel]
        self.n = self.n[self.in_channel]
        return self.in_channel

    def flag_out_channel(self,max_width):
        """Get the indexes of all of the points inside a channel of max_width,
        and remove the points from the list of observations."""

        if N.iterable(max_width): # Map centerline observalble to measurements
            max_distance = max_width[self.index]/2.
        else:
            max_distance = max_width/2.
         
        self.in_channel = N.abs(self.n) <= max_distance

        self.index = self.index[self.in_channel]
        self.d = self.d[self.in_channel]
        self.x = self.x[self.in_channel]
        self.y = self.y[self.in_channel]
        self.s = self.s[self.in_channel]
        self.n = self.n[self.in_channel]
        return self.in_channel

    def get_obs_to_node_map(self,index,minobs=1):
        """Get the mapping from observation to node position (1 -> many); i.e., the inverse
        of index (many -> 1), which maps node position to observations.

        In order for a node to appear, it must have at least minobs observations.
        """

        # Get the list of potential nodes
        
        nodes = N.unique(index)

        self.obs_to_node_map = odict()
        self.nobs = N.zeros(len(self.centerline.x),dtype=N.int32)
        self.populated_nodes = []
        for node in nodes:
            obs_index = N.flatnonzero(index == node)
            nobs = len(obs_index)
            if nobs >= minobs:
                self.populated_nodes.append(node)
                self.obs_to_node_map[node] = obs_index
                self.nobs[node] = nobs
                
        self.n_populated_nodes = len(self.populated_nodes)

        # Store also a list of all the potential nodes and all the unpopulated nodes

        self.n_nodes = len(self.centerline.s)
        self.all_nodes = N.arange(self.n_nodes,dtype=N.int32)
        self.unpopulated_nodes = []
        for node in self.all_nodes:
            if not node in self.populated_nodes:
                self.unpopulated_nodes.append(node)
        self.n_unpopulated_nodes = len(self.unpopulated_nodes)
        
        return self.populated_nodes, self.obs_to_node_map

    def add_obs(self,obs_name,obs):
        """Add an observation as a class variable self.obs_name.
        
        The observation is edited to remove measurements outside
        the channel.

        obs is an iterable of length self.ndata or self.nedited_data.
        """

        if (len(obs) != self.ndata) and (len(obs) != self.nedited_data):
            raise Exception('Observation size incompatible with initial observations')

        if (type(self.max_width) != type(None)) and (len(obs) == self.ndata):
            obs = obs[self.in_channel]
        
        exec('self.%s = obs'%obs_name)

    def obs_to_node(self,obs,node):
        """Get all of the observations in an array obs which map to a node.

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
            return N.array([])

        # If only certain observations have been kept, get the edited vector
        
        if (type(self.max_width) != type(None)) and (len(obs) == self.ndata):
            obs = obs[self.in_channel] 

        return N.asarray(obs)[self.obs_to_node_map[node]]

    def load_nodes(self,vars=[]):
        """Load the desired variables into each of the populated nodes.

        All of the vars should have been loaded previously with add_obs.
        """

        if type(vars) == str:
            vars = [vars]

        self.river_nodes = odict()

        for node in self.populated_nodes:
            d = self.obs_to_node(self.d,node)
            x = self.obs_to_node(self.x,node)
            y = self.obs_to_node(self.y,node)
            s = self.obs_to_node(self.s,node)
            n = self.obs_to_node(self.n,node)
            #h_flg = self.obs_to_node(self.h_flg,node)
            self.river_nodes[node] = self.node_class(node,d,x,y,s,n,ds=self.ds[node])
            for var in vars:
                obs = None # fake cython compiler
                exec('obs = self.obs_to_node(self.%s,node)'%var)
                self.river_nodes[node].add_obs(var,obs,sort=False)

    def get_node_stat(self,stat,var,all_nodes=False):
        """Get a list of results of applying a given stat to a river node variable.
        
        Both stat and var are strings. var should be the name of an instance variable
        for the river node.
        
        A stat is a member function of the river node which returns a
        result given the variable name.

        Example stats are: 'mean', 'std', 'cdf'

        If all_nodes is True, populated and unpopulated nodes are returned.
        Otherwise, only populated nodes are returned.

        The result is a list over desired nodes, with the populated nodes
        holding the result and the unpopulated nodes (when requested) holding the
        missing_value.
        """

        result = []
        if all_nodes:
            for node in self.all_nodes:
                if node in self.populated_nodes:
                    river_node = self.river_nodes[node]
                    exec('result.append( river_node.%s("%s") )'%(stat,var) )
                else:
                    result.append(self.missing_value)
        else:
            for node, river_node in self.river_nodes.iteritems():
                exec('result.append( river_node.%s("%s") )'%(stat,var) )
                
        return result


    def trim_nodes(self,fraction,mode='both',sort_variable='n'):
        """Trim the data in all the nodes.

        fraction: 0 < f < 1. Fraction of the data to remove.
        mode is 'both', 'low', 'high' for which tails of the distribution need to be
        trimmed.

        Prior to trimming, the data are sorted according to the sort variable.
        """

        for node, river_node in self.river_nodes.iteritems():
            river_node.sort(sort_variable=sort_variable)
            river_node.trim(fraction,mode=mode)

    def remove_nodes(self,node_list,reverse=False):
        """Move nodes from the populated node list to the unpopulated node list.

        If reverse is True, move in the opposite direction. No information is
        lost during this process and it is invertible. Both lists are kept
        sorted at each step."""

        if not reverse:
            from_list = copy(self.populated_nodes)
            to_list = copy(self.unpopulated_nodes)
        else:
            to_list = copy(self.populated_nodes)
            from_list = copy(self.unpopulated_nodes)

        for node in node_list:
            try:
                index = from_list.index(node)
                from_list.pop(node)
                to_list.append(node)
            except:
                pass

        from_list.sort()
        to_list.sort()
        if not reverse:
            self.populated_nodes = from_list
            self.unpopulated_nodes = to_list
        else:
            self.populated_nodes = to_list
            self.unpopulated_nodes = from_list
                
                
