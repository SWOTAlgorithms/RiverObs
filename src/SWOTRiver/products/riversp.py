'''
Copyright (c) 2019-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import numpy as np

from SWOTRiver.products.rivertile import L2HRRiverTile

class L2HRRiverSP(L2HRRiverTile):
    """
    Class for River SP data product
    """
    UID = "l2_hr_riversp"
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
        node_sort_idx = np.lexsort((self.nodes.reach_id, self.nodes.node_id))
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
