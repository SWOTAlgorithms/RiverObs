"""
Output the contents of a RiverReach collection into a GIS or hdf5 format.
"""
from __future__ import absolute_import, division, print_function
import os
import shutil
import collections
import numpy as np
import pandas

from GDALOGRUtilities import OGRWriter


class RiverReachWriter:
    """Output the contents of a RiverReach collection into a GIS or hdf5 format.

    Parameters
    ----------

    river_reach_collection : list
        A list of RiverReach instances.
    node_output_variables : list
        A string list of the names of the desired node output variables.
    reach_output_variables : list
        A string list of the desired reach output variables.
    """

    def __init__(self, river_reach_collection, node_output_variables,
                 reach_output_variables):

        self.reaches = river_reach_collection
        self.node_output_variables = node_output_variables
        self.reach_output_variables = reach_output_variables

        # Extract the desired node output data
        self.node_vars = collections.OrderedDict()
        for var in node_output_variables:
            self.node_vars[var] = []

        for reach in self.reaches:
            for var in node_output_variables:
                self.node_vars[var].append(getattr(reach, var))

        # Extract the desired node output data
        self.reach_vars = collections.OrderedDict()
        for var in reach_output_variables:
            self.reach_vars[var] = []

        for reach in self.reaches:
            for var in reach_output_variables:
                #exec('v = reach.metadata["%s"]'%var)
                v = reach.metadata[var]
                self.reach_vars[var].append(v)

    @staticmethod
    def get_type_string(var):
        """
        Returns int, float, or str
        """
        this_string = 'str'

        # If it is a numpy type this will work
        # only needed this try block because float32 doesn't inherit from float
        try:
            if np.issubdtype(var, float):
                this_string = 'float'
            elif np.issubdtype(var, int):
                this_string = 'int'

        # issubdtype raises TypeError on not a numpy array
        except TypeError:
            if isinstance(var, float):
                this_string = 'float'
            elif isinstance(var, int):
                this_string = 'int'

        return this_string

    def write_nodes_ogr(self, output_file, driver='ESRI Shapefile'):
        """Write the nodes as points in a format supporter by OGR."""

        if os.path.isdir(output_file):
            shutil.rmtree(output_file)

        self.fields = collections.OrderedDict()
        self.fields['reach_idx'] = 'int'
        self.fields['node_indx'] = 'int'
        self.fields['reach_indx'] = 'int'

        for var in self.node_output_variables:
            self.fields[var] = self.get_type_string(self.node_vars[var][0][0])

        ogr_writer = OGRWriter(
            output_file,
            layers=[],
            fields=self.fields,
            driver=driver,
            geometry='Point',
            coordinate_system='WGS84')

        for i, reach in enumerate(self.reaches):
            for j in range(len(reach.lat)):
                x = [reach.lon[j]]
                y = [reach.lat[j]]
                field_record = {}
                field_record['reach_idx'] = i
                for var in self.node_vars:
                    if self.fields[var] == 'int':
                        v = int(self.node_vars[var][i][j])
                    elif self.fields[var] == 'float':
                        v = float(self.node_vars[var][i][j])
                    else:
                        v = str(self.node_vars[var][i][j])

                    field_record[var] = v
                ogr_writer.add_xy_feature(x, y, field_record)

        ogr_writer.close()

    def write_reaches_ogr(self, output_file, driver='ESRI Shapefile'):
        """Write the nodes as points in a format supporter by OGR."""

        if os.path.isdir(output_file):
            shutil.rmtree(output_file)

        self.reach_fields = collections.OrderedDict()
        for var in self.reach_output_variables:
            self.reach_fields[var] = self.get_type_string(
                self.reach_vars[var][0])

        ogr_writer = OGRWriter(
            output_file,
            layers=[],
            fields=self.reach_fields,
            driver=driver,
            geometry='LineString',
            coordinate_system='WGS84')

        for i, reach in enumerate(self.reaches):
            x = reach.lon
            y = reach.lat
            field_record = collections.OrderedDict()
            for var in self.reach_vars:
                if self.reach_fields[var] == 'int':
                    v = int(self.reach_vars[var][i])
                elif self.reach_fields[var] == 'float':
                    v = float(self.reach_vars[var][i])
                else:
                    v = str(self.reach_vars[var][i])

                field_record[var] = v
            ogr_writer.add_xy_feature(x, y, field_record)

        ogr_writer.close()

    def write_width_db(self, width_db_file, output_format='h5'):
        """Write a pytables width db file.

        The output_format can be 'h5' or 'csv'.
        """

        for i, reach in enumerate(self.reaches):
            n = len(reach.lat)

            river_data = collections.OrderedDict()
            river_data['width'] = reach.width.astype(np.int16)
            try:
                river_data['nchannels'] = reach.nchannels.astype(np.int8)
            except:
                river_data['nchannels'] = np.ones(n, dtype=np.int8)
            try:
                river_data['reservoir'] = reach.reservoir.astype(np.int8)
            except:
                river_data['reservoir'] = np.zeros(n, dtype=np.int8)
            river_data['long'] = reach.lon.astype(np.float32)
            river_data['lat'] = reach.lat.astype(np.float32)
            river_data['reach_index'] = np.ones(n, dtype=np.int32) * i

            try:
                ds = arc_distance_xy(reach.x, reach.y)
            except:
                ds = arc_distance(reach.lon, reach.lat)
            river_data['reach'] = np.cumsum(ds).astype(np.float32)

            df = DataFrame(river_data)
            if i == 0:
                river_df = df
            else:
                river_df = river_df.append(df, ignore_index=True)

        # Fill in the reach data base

        break_idx = []
        npoints = []
        reach_index = []
        total_reach = []
        lonmin = []
        lonmax = []
        latmin = []
        latmax = []
        width_mean = []
        width_std = []
        width_min = []
        width_max = []
        ibreak = -1
        for i, reach in enumerate(self.reaches):
            n = len(reach.lat)
            ibreak += n
            break_idx.append(ibreak)
            npoints.append(n)
            reach_index.append(i)
            lonmin.append(reach.lon.min())
            lonmax.append(reach.lon.max())
            latmin.append(reach.lat.min())
            latmax.append(reach.lat.max())
            width_mean.append(reach.width.mean())
            width_std.append(reach.width.std())
            width_min.append(reach.width.min())
            width_max.append(reach.width.max())
            try:
                ds = arc_distance_xy(reach.x, reach.y)
            except:
                ds = arc_distance(reach.lon, reach.lat)
            s = np.cumsum(ds)
            total_reach = s[-1]

        reach_df = pandas.DataFrame({
            'break_idx': break_idx,
            'npoints': npoints,
            'reach': total_reach,
            'lonmin': lonmin,
            'lonmax': lonmax,
            'latmin': latmin,
            'latmax': latmax,
            'width_mean': width_mean,
            'width_std': width_std,
            'width_min': width_min,
            'width_max': width_max,
        })

        # Write the width data base

        if output_format == 'h5':
            if not '.h5' in width_db_file:
                width_db_file = width_db_file + '.h5'
            store = pandas.HDFStore(width_db_file, mode='w', complevel=9)
            store['river'] = river_df
            store['reach'] = reach_df
            store.close()
        else:
            river_df.to_csv(width_db_file + '_river_df.csv')
            reach_df.to_csv(width_db_file + '_reach_df.csv')

        return river_df, reach_df


def arc_distance(lon, lat):
    """Calculate distance in meters between subsequent points."""

    R = 6378.e3  # approximate earth radius
    deg2rad = np.pi / 180.

    lon = np.asarray(lon)
    lat = np.asarray(lat)
    d = np.zeros(len(lon), dtype=np.float32)

    dlon = (lon[1:] - lon[0:-1]) * deg2rad
    dlat = (lat[1:] - lat[0:-1]) * deg2rad
    latmean = np.mean(lat) * deg2rad

    dx = dlon * np.cos(latmean) * R
    dy = dlat * R

    d[0:-1] = np.sqrt(dx**2 + dy**2)

    return d


def arc_distance_xy(x, y):
    """Calculate the arc distance, give Cartesian coordinates."""

    d = np.zeros(len(x), dtype=np.float32)
    dx = x[1:] - x[0:-1]
    dy = y[1:] - y[0:-1]

    d[0:-1] = np.sqrt(dx**2 + dy**2)

    return d
