#!/usr/bin/env python
'''
Analyze RiverObs outputs.
'''
import os
import netCDF4
import pytest
import shapefile
import numpy as np

class RiverTilePart(object):
    '''Class for part of RiverTile product (nodes/reaches)'''
    @classmethod
    def from_shapefile(cls, filename):
        '''makes self from a filename'''
        klass = cls()
        ifp = shapefile.Reader(filename)
        records = np.array(ifp.records())
        fields = ifp.fields[1:]
        for ii, field in enumerate(fields):
            setattr(klass, field[0], records[:, ii])
        return klass

class RiverTile(object):
    '''Class for RiverTile shapefile products'''
    @classmethod
    def from_shapefiles(cls, node_filename, reach_filename):
        '''makes self from two shapefiles (nodes.shp/reaches.shp)'''
        klass = cls()
        klass.nodes = RiverTilePart.from_shapefile(node_filename)
        klass.reaches = RiverTilePart.from_shapefile(reach_filename)
        return klass

    @classmethod
    def from_ncfile(cls, filename):
        '''Constructifies self from a netcdf file'''
        klass = cls()
        with netCDF4.Dataset(filename, 'r') as ifp:
            for group in ifp.groups.keys():
                setattr(klass, group, RiverTilePart())
                for key, value in ifp.groups[group].variables.items():
                    setattr(getattr(klass, group), key, value[:])
        return klass

class Index(object):
    """class for index file"""
    @classmethod
    def from_ncfile(cls, filename):
        """makes self from a netcdf file"""
        klass = cls()
        with netCDF4.Dataset(filename, 'r') as ifp:
            for varname, varvalue in ifp.variables.items():
                setattr(klass, varname, varvalue[:])
        return klass


class RiverObsTester(object):
    def __init__(self, nodes_shpbase, reaches_shpbase, index_file,
                 flat_height=100, std_height=0.005):

        self.rivertile = RiverTile.from_shapefiles(
            nodes_shpbase, reaches_shpbase)
        self.index = Index.from_ncfile(index_file)
        self.flat_height = flat_height
        self.std_height = std_height

    @staticmethod
    def _get_stat(values):
        values_ = values[values==values]
        return {'mean': np.mean(values_), 'std': np.std(values_)}

    def get_index_stats(self):
        stats = {}
        for key, value in self.index.__dict__.items():
            stats[key] = self._get_stat(value)
        return stats

    def get_rivertile_stats(self):
        stats = {}
        for part in ['nodes', 'reaches']:
            this_base = getattr(self.rivertile, part)
            stats[part] = {}
            for key, value in this_base.__dict__.items():
                stats[part][key] = self._get_stat(value)
        return stats

    def analyze(self):
        """Analyze rivertile outputs for expected values"""
        self.stats = {
            'rivertile': self.get_rivertile_stats(),
            'index': self.get_index_stats()}


class TestNosiyFlatRiverTile():
    @pytest.fixture(scope='class')
    def rivertile_tester(
        self, nodes_shpbase=None, reaches_shpbase=None, index_file=None):
        """Constructs the RiverTileTester class instance to use for testing"""

        if nodes_shpbase is None:
            nodes_shpbase = os.path.join(os.getcwd(), 'nodes', 'nodes')
            reaches_shpbase = os.path.join(os.getcwd(), 'reaches', 'reaches')
            index_file = os.path.join(os.getcwd(), 'index.nc')

        rivertile_tester = RiverObsTester(
            nodes_shpbase, reaches_shpbase, index_file, flat_height=100.0,
            std_height=0.005)

        rivertile_tester.analyze()
        return rivertile_tester

    def test_reach_height(self, rivertile_tester):
        """tests the reach average heights"""
        for key in ['h_no', 'h_nr', 'h_nw']:
            mean_height = rivertile_tester.stats[
                'rivertile']['reaches'][key]['mean']
            assert rivertile_tester.flat_height == pytest.approx(
                mean_height, abs=rivertile_tester.std_height/10)

    def test_node_height(self, rivertile_tester):
        """tests the node average heights"""
        mean_height = rivertile_tester.stats[
            'rivertile']['nodes']['h_n_ave']['mean']
        assert rivertile_tester.flat_height == pytest.approx(
            mean_height, abs=rivertile_tester.std_height/10)

    def test_pixc_height(self, rivertile_tester):
        """tests the pixc average heights"""
        mean_height = rivertile_tester.stats[
            'index']['height_vectorproc']['mean']
        assert rivertile_tester.flat_height == pytest.approx(
            mean_height, abs=rivertile_tester.std_height/10)

    def test_pixc_height_std(self, rivertile_tester):
        """tests the pixc std of the height"""
        std_height = rivertile_tester.stats[
            'index']['height_vectorproc']['std']
        assert rivertile_tester.std_height == pytest.approx(
            std_height, abs=rivertile_tester.std_height/100)

    def test_areas(self, rivertile_tester):
        """Checks the nodes and reaches have same area"""
        node_area = rivertile_tester.rivertile.nodes.area.sum()
        reach_area = rivertile_tester.rivertile.reaches.area.sum()
        assert node_area == pytest.approx(reach_area, abs=10)

    def test_unique_node_assingment(self, rivertile_tester):
        """Ensures no PIXC pixels are double-assinged to nodes"""
        idx = np.array([
            rivertile_tester.index.range_index,
            rivertile_tester.index.azimuth_index])
        assert idx.shape == np.unique(idx, axis=1).shape

    def test_no_nans(self, rivertile_tester):
        """Checks for nans in node heights"""
        node_hgt = rivertile_tester.rivertile.nodes.h_a_ave
        assert all(node_hgt == node_hgt)
