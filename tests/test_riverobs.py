#!/usr/bin/env python
'''
Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore

Regression tests for RiverTile / RiverObs

Execute tests with:
pytest -v test_riverobs.py

To execute tests with debugging:
pytest -v test_riverobs.py --capture=no
'''
import os
import netCDF4
import pytest
import contextlib
import numpy as np

from SWOTRiver.products.rivertile import L2HRRiverTile
from SWOTRiver.products.pixcvec import L2PIXCVector

class RiverObsTester(object):
    def __init__(self, test_rivertile_ncfile, test_pixcvec_ncfile,
                 ref_rivertile_ncfile, ref_pixcvec_ncfile):

        self.test_rivertile = L2HRRiverTile.from_ncfile(test_rivertile_ncfile)
        self.test_pixcvec = L2PIXCVector.from_ncfile(test_pixcvec_ncfile)
        self.ref_rivertile = L2HRRiverTile.from_ncfile(ref_rivertile_ncfile)
        self.ref_pixcvec = L2PIXCVector.from_ncfile(ref_pixcvec_ncfile)

class TestRiverTile():
    @pytest.fixture(scope='class')
    def rivertile_tester(self):
        """Constructs the RiverTileTester class instance to use for testing"""

        test_rivertile_ncfile = os.path.join(os.getcwd(), 'test_case', 'rt.nc')
        test_pixcvec_ncfile = os.path.join(os.getcwd(), 'test_case', 'pv.nc')
        ref_rivertile_ncfile = os.path.join(
            os.getcwd(), 'reference_case', 'rt.nc')
        ref_pixcvec_ncfile = os.path.join(
            os.getcwd(), 'reference_case', 'pv.nc')
        rivertile_tester = RiverObsTester(
            test_rivertile_ncfile, test_pixcvec_ncfile, ref_rivertile_ncfile,
            ref_pixcvec_ncfile)
        return rivertile_tester

    def test_in_valid_range(self, rivertile_tester):
        """Checks that all variables are within the valid range"""
        # This test has been refactored into ProductTesterMixIn
        rivertile_tester.test_rivertile.test_in_valid_range()

    def test_inf_nan(self, rivertile_tester):
        """Checks for non-finite values in floating point rivertile outputs"""
        # This test has been refactored into ProductTesterMixIn
        rivertile_tester.test_rivertile.test_inf_nan()

    def test_node_q_b_valid_max(self, rivertile_tester):
        """Checks that node_q_b valid_max is correct"""
        test_rt = rivertile_tester.test_rivertile
        node_q_b_dict = test_rt.nodes.VARIABLES['node_q_b']
        assert node_q_b_dict['flag_masks'].sum() == node_q_b_dict['valid_max']

    def test_reach_q_b_valid_max(self, rivertile_tester):
        """Checks that reach_q_b valid_max is correct"""
        test_rt = rivertile_tester.test_rivertile
        reach_q_b_dict = test_rt.reaches.VARIABLES['reach_q_b']
        assert reach_q_b_dict['flag_masks'].sum() == reach_q_b_dict['valid_max']

    def test_populated_nodes(self, rivertile_tester):
        """
        Checks that the same nodes are populated in the test and reference
        RiverTiles
        """
        test_rt = rivertile_tester.test_rivertile
        ref_rt = rivertile_tester.ref_rivertile
        # Check that same nodes have valid wse in test and ref
        assert (test_rt.nodes.wse.mask == ref_rt.nodes.wse.mask).all()

    def test_populated_reaches(self, rivertile_tester):
        """
        Checks that the same nodes are populated in the test and reference
        RiverTiles
        """
        test_rt = rivertile_tester.test_rivertile
        ref_rt = rivertile_tester.ref_rivertile
        # Check that same nodes have valid wse in test and ref
        assert (test_rt.reaches.wse.mask == ref_rt.reaches.wse.mask).all()

    def test_sensible_outputs(self, rivertile_tester):
        """
        Checks that there are not any insensible combinations of outputs
        """
        test_rt = rivertile_tester.test_rivertile

        # Test that any node with area has a width
        valid_area_mask = ~test_rt.nodes.area_total.mask
        assert not test_rt.nodes.width[valid_area_mask].mask.any()

        # Test that any reach with valid slope also has valid wse
        valid_slope_mask = ~test_rt.reaches.slope.mask
        assert not test_rt.reaches.wse[valid_slope_mask].mask.any()

    def test_prior_river_database_fields(self, rivertile_tester):
        """
        Checks that various PRD fields in the output product are unchanged
        """
        test_rt = rivertile_tester.test_rivertile
        ref_rt = rivertile_tester.ref_rivertile

        prd_node_vars = [
            'reach_id', 'node_id', 'lat_prior', 'lon_prior', 'p_wse',
            'p_wse_var', 'p_width', 'p_wid_var', 'p_dist_out', 'p_length',
            'p_dam_id', 'p_n_ch_max', 'p_n_ch_mod']

        prd_reach_vars = [
            'reach_id', 'centerline_lat', 'centerline_lon', 'p_lat', 'p_lon',
            'river_name', 'ice_clim_f', 'n_reach_up', 'n_reach_dn',
            'rch_id_up', 'rch_id_dn', 'p_wse', 'p_wse_var', 'p_width',
            'p_wid_var', 'p_n_nodes', 'p_dist_out', 'p_length', 'p_maf',
            'p_dam_id', 'p_n_ch_max', 'p_n_ch_mod', 'p_low_slp']

        for var in prd_node_vars:
            test_values = test_rt.nodes[var]
            test_values = test_values[~test_values.mask].flatten()

            ref_values = ref_rt.nodes[var]
            ref_values = ref_values[~ref_values.mask].flatten()
            print('test_prior_river_database_fields: /nodes/%s'%var)
            for test_value, ref_value in zip(test_values, ref_values):
                assert test_value == ref_value

        for var in prd_reach_vars:
            test_values = test_rt.reaches[var]
            test_values = test_values[~test_values.mask].flatten()

            ref_values = ref_rt.reaches[var]
            ref_values = ref_values[~ref_values.mask].flatten()
            print('test_prior_river_database_fields: /reaches/%s'%var)
            for test_value, ref_value in zip(test_values, ref_values):
                assert test_value == ref_value
