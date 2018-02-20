import netCDF4
import os
import configparser
import RiverObs
import configparser
import numpy as np

NC_METADATA_FILE = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), 'nc_metadata.cfg')


CASTED_METADATA_NAMES = ['_FillValue', 'valid_max', 'valid_min']
FLOAT_METADATA_NAMES = [
     'geospatial_lat_max', 'geospatial_lat_min', 'geospatial_lon_max',
     'geospatial_lon_min']
INT_METADATA_NAMES = ['NONE',]

def gather_outputs(reach_collection):
    """
    Gathers outputs for river tile data product file
    """
    reach_variables = list(reach_collection[0].metadata.keys())
    node_variables = list(reach_collection[0].__dict__.keys())
    node_variables.remove('ds')
    node_variables.remove('metadata')

    num_nodes_per_reach = [len(item.lat) for item in reach_collection]
    num_nodes = sum(num_nodes_per_reach)

    node_outputs = {}
    reach_outputs = {}
    for node_variable in node_variables:
        node_outputs[node_variable] = np.concatenate(
            [getattr(reach, node_variable) for reach in reach_collection])

    for reach_variable in reach_variables:
        reach_outputs[reach_variable] = np.array([
            reach.metadata[reach_variable] for reach in reach_collection])

    node_outputs['reach_idx'] = np.zeros(node_outputs['lat'].shape)
    i_start = 0
    for ireach, num_nodes in enumerate(num_nodes_per_reach):
        node_outputs['reach_idx'][i_start:i_start+num_nodes] = ireach
        i_start = i_start + num_nodes

    return node_outputs, reach_outputs

def cast_metadata(name, value, base_id):
    if name in INT_METADATA_NAMES:
        value = int(value)
    elif name in FLOAT_METADATA_NAMES:
        value = float(value)
    elif name in CASTED_METADATA_NAMES:
        value = base_id.dtype.type(value)

    return value

def write(nc_file, node_outputs, reach_outputs):
    with netCDF4.Dataset(nc_file, 'w') as ofp:
        ofp.createDimension('node_record', node_outputs['reach_idx'].shape[0])
        ofp.createDimension('reach_record', reach_outputs['area'].shape[0])

        for varname, varvalue in node_outputs.items():
            var = ofp.createVariable(
                "/nodes/%s" % varname, varvalue.dtype.str, ('node_record',))
            var[:] = varvalue

        for varname, varvalue in reach_outputs.items():
            var = ofp.createVariable(
                "/reaches/%s" % varname, varvalue.dtype.str, ('reach_record',))
            var[:] = varvalue


def fixup_metadata(ncfile):
    parser = configparser.ConfigParser()
    parser.read(NC_METADATA_FILE)

    with netCDF4.Dataset(ncfile, 'a') as ifp:
        for section in parser.keys():

            if section == "DEFAULT":
                continue

            # globals would go here
            if section == 'GLOBAL':
                base_id = ifp

            else:
                group = section.split('/')[1]
                varname = section.split('/')[2]
                base_id = ifp.groups[group][varname]

            for attr_name, attr_value in parser[section].items():
                attr_value = cast_metadata(attr_name, attr_value, base_id)
                setattr(base_id, attr_name, attr_value)


