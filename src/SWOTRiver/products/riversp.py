'''
Copyright (c) 2019-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import copy
import numpy as np
from collections import OrderedDict as odict

from SWOTWater.products.product import Product, FILL_VALUES, textjoin
from SWOTRiver.products.rivertile import \
    L2HRRiverTile, RiverTileNodes, RiverTileReaches, RIVER_PRODUCT_ATTRIBUTES

RIVERSP_ATTRIBUTES = copy.deepcopy(RIVER_PRODUCT_ATTRIBUTES)
RIVERSP_ATTRIBUTES['xref_input_l2_hr_rivertile_files'] = {
    'dtype': 'str',
    'docstr': 'List of RiverTile products used to generate data in product'}
RIVERSP_ATTRIBUTES['xref_l2_hr_river_sp_param_file'] = {
    'dtype': 'str',
    'docstr': textjoin("""
        Name of PGE_L2_HR_RiverSP parameter file used to generate data in
        product""")}

ATTRIBUTE_KEYS2POP = [
    "_".join([a, b, c]) for a in ['inner', 'outer'] for b in ['first', 'last']
    for c in ['latitude', 'longitude']] + ['xref_l2_hr_river_tile_param_file',]

for key in ATTRIBUTE_KEYS2POP:
    RIVERSP_ATTRIBUTES.pop(key, None)

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

    def __add__(self, other):
        """Adds other to self"""
        klass = L2HRRiverSP()
        klass.nodes = self.nodes + other.nodes
        klass.reaches = self.reaches + other.reaches
        return klass

    @staticmethod
    def dump_xmls(node_xml_file, reach_xml_file):
        with open(node_xml_file, 'w') as ofp:
            RiverSPNodes.print_xml(ofp=ofp, is_shapefile=True)
        with open(reach_xml_file, 'w') as ofp:
            RiverSPReaches.print_xml(ofp=ofp, is_shapefile=True)


class RiverSPNodes(RiverTileNodes):
    ATTRIBUTES = RIVERSP_ATTRIBUTES.copy()
    ATTRIBUTES['title'] = {
        'dtype': 'str',
        'value': 'Level 2 KaRIn High Rate River Single Pass Vector Product - Node',
        'docstr': 'Level 2 KaRIn High Rate River Single Pass Vector Product - Node'}

class RiverSPReaches(RiverTileReaches):
    ATTRIBUTES = RIVERSP_ATTRIBUTES.copy()
    ATTRIBUTES['title'] = {
        'dtype': 'str',
        'value': 'Level 2 KaRIn High Rate River Single Pass Vector Product - Reach',
        'docstr': 'Level 2 KaRIn High Rate River Single Pass Vector Product - Reach'}
