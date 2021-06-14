'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Dustin Lagoy

Tools for working with memory-mapped variables.
'''
from __future__ import print_function
import os
import netCDF4 as nc
import numpy as np

from SWOTWater.products.constants import FILL_VALUES

DEPTH_DIMNAMES = ['depth', 'complex_depth']

def get_fill(dtype):
    """Get the SWOT approved fill value for a data type"""
    dtype_str = np.dtype(dtype).str[1:]
    if (dtype_str[0] == 'S') or (dtype_str[0] == 'U') or (dtype_str[0]=='O'):
        # handle arbitrary-length strings
        dtype_str = 'S1'
    return FILL_VALUES[dtype_str]


def set_variable(dataset, key, array, dimensions, attributes=None):
    '''Set the NetCDF variable, dealing with complex numbers.

    If array is complex, it is stored in the dataset with a third dimension,
    'complex_depth', so that:
        variable.shape == (lines, pixels, 2)
        variable[:, :, 0] == array.real
        variable[:, :, 1] == array.imag
    '''
    # np.ma.mask_array has fill_value attr, else use default fill value
    if attributes is None:
        fill_value = getattr(array, 'fill_value', get_fill(array.dtype))
    else:
        fill_value = attributes.get('_FillValue', get_fill(array.dtype))
        # drop complex part of fill_value
        if 'complex' in array.dtype.name:
            fill_value = np.real(fill_value)

    complevel = attributes.get('complevel', None)

    def _make_variable(key, data, dimensions, attributes=None, complevel=None):
        dtype = data.dtype
        if dtype == object:
            # assume that objects are strings
            dtype = str
        if complevel is not None:
            if complevel not in range(1, 10):
                raise Exception(
                    'Invalid complevel {} specified in _make_variable'.format(
                    complevel))
            dataset.createVariable(
                key, dtype, dimensions, fill_value=fill_value, zlib=True,
                complevel=complevel)
        else:
            dataset.createVariable(
                key, dtype, dimensions, fill_value=fill_value)

        if ((data.dtype.str[1] == 'S') or (data.dtype.str[1] == 'U') or (
                data.dtype.str[1] == 'O')) and np.ma.isMaskedArray(data):
            dataset[key][:] = data.data
        else:
            dataset[key][:] = data
        if attributes is not None:
            for name, value in attributes.items():
                if name in ('dtype', 'dimensions', '_FillValue', 'complevel'):
                    continue
                if np.iscomplexobj(value):
                    value = value.real
                # cast min/max/fill
                if name in ['valid_min', 'valid_max']:
                    value = data.dtype.type(value)
                dataset[key].setncattr(name, value)

    if 'complex' in array.dtype.name:
        # Add the depth dimension
        if 'complex_depth' not in dataset.dimensions:
            dataset.createDimension('complex_depth', 2)
        shape = array.shape
        n_bytes = int(array.dtype.itemsize/2)
        float_type = np.dtype('f'+str(n_bytes))
        if isinstance(array, np.ma.core.MaskedArray):
            # Somehow MaskedArray.view() doesn't work, so convert to a normal
            # array.
            array = array.filled()
        tmp = array.view(dtype=float_type).reshape(shape+(2,))
        _make_variable(
            key, tmp, dimensions+['complex_depth'], attributes,
            complevel=complevel)
    else:
        _make_variable(key, array, dimensions, attributes, complevel=complevel)


def get_variable_keys(dataset):
    '''Return a list of NetCDF variables, dealing with complex numbers.'''
    # Now, complex variables are stored with the same names, so this function
    # is trivial.
    return [key for key in dataset.variables]


def get_variable(dataset, key):
    '''Get the NetCDF variable and dimensions, dealing with complex numbers.

    Key is assumed to be complex if the variable has a first or last dimension
    of 'complex_depth' or 'depth' with length 2.
    '''
    variable = dataset[key]
    if len(variable.dimensions) == 0:
        return variable[0]
    if variable.dimensions[0] in DEPTH_DIMNAMES and variable.shape[0] == 2:
        tmp = np.ma.MaskedArray(variable[0] + 1j*variable[1])
        return tmp
    if variable.dimensions[-1] in DEPTH_DIMNAMES and variable.shape[-1] == 2:
        n_bytes = int(variable.dtype.itemsize*2)
        complex_type = np.dtype('c'+str(n_bytes))
        tmp = variable[:]
        fill_value = get_fill(complex_type)
        mask = None
        if isinstance(tmp, np.ma.MaskedArray):
            # Somehow MaskedArray.view() doesn't work, so convert to a normal
            # array.
            # Keep the current fill value and mask, TODO: use default fill?
            fill_value = tmp.fill_value
            # Assume the real mask matches the imaginary mask
            if isinstance(tmp.mask, np.ndarray):
                mask = tmp.mask[..., 0]
            else:
                mask = tmp.mask
            tmp = tmp.data
        # Make the data complex
        tmp = tmp.view(dtype=complex_type).squeeze()
        # Turn back into a masked array and re-fill
        tmp = tmp.view(np.ma.MaskedArray)
        tmp.set_fill_value(fill_value)
        if mask is not None:
            tmp[mask] = np.ma.masked
        return tmp
    variable = variable[:]
    if isinstance(variable, np.ma.MaskedArray):
        return variable
    variable = variable.view(np.ma.MaskedArray)
    if dataset[key].dtype is str:
        # handle arrays of strings
        fill_value = get_fill('S1')
    else:
        fill_value = get_fill(dataset[key].dtype.str[1:])
    variable.set_fill_value(fill_value)
    return variable


def get_variable_dimensions(dataset, key):
    '''Get the NetCDF variable dimensons, dealing with complex numbers.

    Key is assumed to be complex if the variable has a first or last dimension
    of 'depth' with length 2. This dimension is removed from the returned list.
    '''
    variable = dataset[key]
    if len(variable.dimensions) == 0:
        return variable.dimensions
    if variable.dimensions[0] in DEPTH_DIMNAMES and variable.shape[0] == 2:
        return variable.dimensions[1:]
    if variable.dimensions[-1] in DEPTH_DIMNAMES and variable.shape[-1] == 2:
        return variable.dimensions[:-1]
    return variable.dimensions
