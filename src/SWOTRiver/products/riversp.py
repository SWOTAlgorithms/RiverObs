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
from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9

RIVERSP_ATTRIBUTES = copy.deepcopy(RIVER_PRODUCT_ATTRIBUTES)
RIVERSP_ATTRIBUTES['short_name'] = {
    'dtype': 'str', 'docstr': 'L2_HR_RiverSP', 'value': 'L2_HR_RiverSP'}
ATTRIBUTE_KEYS2POP = [
    "_".join([a, b, c]) for a in ['inner', 'outer'] for b in ['first', 'last']
    for c in ['latitude', 'longitude']] + [
        'xref_l2_hr_river_tile_param_file', 'tile_number', 'swath_side',
        'tile_name']

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
    ATTRIBUTES['product_file_id'] = {
        'dtype': 'str', 'value': 'Node', 'docstr': 'Node'}

    def __add__(self, other):
        """Adds other to self"""
        klass = RiverSPNodes()
        for key in klass.VARIABLES:
            setattr(klass, key, np.concatenate((
                getattr(self, key), getattr(other, key))))
        return klass


class RiverSPReaches(RiverTileReaches):
    ATTRIBUTES = RIVERSP_ATTRIBUTES.copy()
    ATTRIBUTES['title'] = {
        'dtype': 'str',
        'value': 'Level 2 KaRIn High Rate River Single Pass Vector Product - Reach',
        'docstr': 'Level 2 KaRIn High Rate River Single Pass Vector Product - Reach'}
    ATTRIBUTES['product_file_id'] = {
        'dtype': 'str', 'value': 'Reach', 'docstr': 'Reach'}
    def __add__(self, other):
        """Adds other to self"""
        self_n_reaches, self_n_centerlines = self.centerline_lon.shape
        other_n_reaches, other_n_centerlines = other.centerline_lon.shape
        cl_len = max([self_n_centerlines, other_n_centerlines])

        cl_lon = MISSING_VALUE_FLT * np.ones(
            [self_n_reaches+other_n_reaches, cl_len])
        cl_lat = MISSING_VALUE_FLT * np.ones(
            [self_n_reaches+other_n_reaches, cl_len])

        cl_lon[0:self_n_reaches, 0:self_n_centerlines] = self.centerline_lon
        cl_lat[0:self_n_reaches, 0:self_n_centerlines] = self.centerline_lat
        cl_lon[self_n_reaches:, 0:other_n_centerlines] = other.centerline_lon
        cl_lat[self_n_reaches:, 0:other_n_centerlines] = other.centerline_lat

        klass = RiverSPReaches()
        for key in klass.VARIABLES:
            if key in ['centerline_lon', 'centerline_lat']:
                value = {
                    'centerline_lon': cl_lon, 'centerline_lat': cl_lat}[key]
                setattr(klass, key, value)
            else:
                setattr(klass, key, np.concatenate((
                    getattr(self, key), getattr(other, key))))
        return klass


