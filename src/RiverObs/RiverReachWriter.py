"""
Output the contents of a RiverReach collection into a GIS or hdf5 format.
"""

from collections import OrderedDict as odict
import numpy as N
from pandas import HDFStore, DataFrame
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
    def __init__(self,river_reach_collection,node_output_variables,reach_output_variables):

        self.reaches = river_reach_collection
        self.node_output_variables = node_output_variables
        self.reach_output_variables = reach_output_variables

        # Extract the desired node output data

        self.node_vars = odict()
        for var in node_output_variables:
            self.node_vars[var] = []
            
        for reach in self.reaches:
            for var in node_output_variables:
                exec('v = reach.%s'%var)
                self.node_vars[var].append(v)

        # Extract the desired node output data

        self.reach_vars = odict()
        for var in reach_output_variables:
            self.reach_vars[var] = []
            
        for reach in self.reaches:
            for var in reach_output_variables:
                exec('v = reach.metadata["%s"]'%var)
                self.reach_vars[var].append(v)

    def write_nodes_ogr(self,output_file,driver='ESRI Shapefile'):
        """Write the nodes as points in a format supporter by OGR."""

        self.fields = {}
        self.fields['reach_index'] = 'int'
        for var in self.node_output_variables:
            if  ((self.node_vars[var][0].dtype == N.float64) or
                 (self.node_vars[var][0].dtype == N.float32) ):
                self.fields[var] = 'float'
            elif ((self.node_vars[var][0].dtype == N.int32) or
                 (self.node_vars[var][0].dtype == N.int16) ):
                self.fields[var] = 'int'
            else:
                self.fields[var] = 'str'

        ogr_writer = OGRWriter(output_file,layers=[],fields=self.fields,
                               driver=driver,
                                geometry='Point',coordinate_system='WGS84')

        for i,reach in enumerate(self.reaches):
            for j in range(len(reach.lat)):
                x = [reach.lon[j]]
                y = [reach.lat[j]]
                field_record = {}
                field_record['reach_index'] = i
                for var in self.node_vars:
                    if self.fields[var] == 'int':
                        v = int(self.node_vars[var][i][j])
                    elif self.fields[var] == 'float':
                        v = float(self.node_vars[var][i][j])
                    else:
                        v = str(self.node_vars[var][i][j])
                        
                    field_record[var] = v
                ogr_writer.add_xy_feature(x,y,field_record)
                 
        ogr_writer.close()

    def write_reaches_ogr(self,output_file,driver='ESRI Shapefile'):
        """Write the nodes as points in a format supporter by OGR."""

        self.reach_fields = {}
        for var in self.reach_output_variables:
            t = type(self.reach_vars[var][0])
            if  (t == float):
                self.reach_fields[var] = 'float'
            elif (t == int ):
                self.reach_fields[var] = 'int'
            else:
                self.reach_fields[var] = 'str'

        ogr_writer = OGRWriter(output_file,layers=[],fields=self.reach_fields,
                               driver=driver,
                                geometry='LineString',coordinate_system='WGS84')

        for i,reach in enumerate(self.reaches):
            x = reach.lon
            y = reach.lat
            field_record = {}
            for var in self.reach_vars:
                if self.reach_fields[var] == 'int':
                    v = int(self.reach_vars[var][i])
                elif self.reach_fields[var] == 'float':
                    v = float(self.reach_vars[var][i])
                else:
                    v = str(self.reach_vars[var][i])
                        
                field_record[var] = v
            ogr_writer.add_xy_feature(x,y,field_record)
                 
        ogr_writer.close()


        
