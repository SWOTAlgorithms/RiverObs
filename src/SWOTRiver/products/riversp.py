'''
Copyright (c) 2019-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import copy
import numpy as np
from collections import OrderedDict as odict

from SWOTRiver.products.rivertile import \
    L2HRRiverTile, RiverTileNodes, RiverTileReaches, RIVERTILE_ATTRIBUTES

RIVERSP_ATTRIBUTES = copy.deepcopy(RIVERTILE_ATTRIBUTES)
RIVERSP_ATTRIBUTES['title']['docstr'] = \
        'Level 2 KaRIn High Rate River Single Pass Vector Product'
for key in ['Conventions', 'title', 'platform']:
    RIVERSP_ATTRIBUTES[key]['value'] = RIVERSP_ATTRIBUTES[key]['docstr']


class L2HRRiverSP(L2HRRiverTile):
    """
    Class for River SP data product
    """
    UID = "l2_hr_riversp"
    ATTRIBUTES = odict()
    GROUPS = odict([
        ['nodes', 'RiverSPNodes'],
        ['reaches', 'RiverSPReaches'],
    ])

    @classmethod
    def from_rivertiles(cls, rivertiles):
        """
        Builds a River SP product from a list of rivetile data products
        """
        klass = cls()
        # It does somthing like this
        for rivertile in rivertiles:
           klass += rivertile

        # sort them by increasing reach id
        klass.sort()
        return klass

    def sort(self):
        """sorts self according to the PDD"""
        # sort first by reach_id, then by node_id
        node_sort_idx = np.argsort(self.nodes.node_id)
        for key, values in self.nodes.variables.items():
            self.nodes[key] = values[node_sort_idx]

        reach_sort_idx = np.argsort(self.reaches.reach_id)
        for key, values in self.reaches.variables.items():
            self.reaches[key] = values[reach_sort_idx]

    def __add__(self, other):
        """Adds other to self"""
        klass = L2HRRiverSP()
        klass.nodes = self.nodes + other.nodes
        klass.reaches = self.reaches + other.reaches
        return klass

class RiverSPNodes(RiverTileNodes):
    ATTRIBUTES = RIVERSP_ATTRIBUTES

class RiverSPReaches(RiverTileReaches):
    ATTRIBUTES = RIVERSP_ATTRIBUTES

