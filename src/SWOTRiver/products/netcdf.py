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


FILL_VALUE = -9.99e9


def set_variable(dataset, key, array, dimensions, attributes=None):
    '''Set the NetCDF variable, dealing with complex numbers.

    If array is complex, it is stored in the dataset with a third dimension,
    'depth', so that:
        variable.shape == (lines, pixels, 2)
        variable[:, :, 0] == array.real
        variable[:, :, 1] == array.imag
    '''
    # np.ma.mask_array has fill_value attr, else use FILL_VALUE
    fill_value = getattr(array, 'fill_value', FILL_VALUE)
    def _make_variable(key, data, dimensions, attributes=None):
        dataset.createVariable(
            key, data.dtype, dimensions, fill_value=fill_value)
        dataset[key][:] = data
        if attributes is not None:
            for name, value in attributes.items():
                if name in ('dtype', 'dimensions', '_FillValue'):
                    continue
                if np.iscomplexobj(value):
                    value = value.real

                # cast min/max/fill
                if name in ['valid_min', 'valid_max']:
                    value = data.dtype.type(value)

                dataset[key].setncattr(name, value)
    if 'complex' in array.dtype.name:
        # Add the depth dimension
        if 'depth' not in dataset.dimensions:
            dataset.createDimension('depth', 2)
        shape = array.shape
        n_bytes = int(array.dtype.itemsize/2)
        float_type = np.dtype('f'+str(n_bytes))
        if isinstance(array, np.ma.core.MaskedArray):
            # Somehow MaskedArray.view() doesn't work, so convert to a normal
            # array.
            array = array.filled()
        tmp = array.view(dtype=float_type).reshape(shape+(2,))
        _make_variable(key, tmp, dimensions+['depth'], attributes)
    else:
        _make_variable(key, array, dimensions, attributes)


def get_variable_keys(dataset):
    '''Return a list of NetCDF variables, dealing with complex numbers.'''
    # Now, complex variables are stored with the same names, so this function
    # is trivial.
    return [key for key in dataset.variables]


def get_variable(dataset, key):
    '''Get the NetCDF variable and dimensions, dealing with complex numbers.

    Key is assumed to be complex if the variable has a first or last dimension
    of 'depth' with length 2.
    '''
    variable = dataset[key]
    if variable.dimensions[0] == 'depth' and variable.shape[0] == 2:
        tmp = np.ma.MaskedArray(variable[0] + 1j*variable[1])
        tmp.set_fill_value(FILL_VALUE + 1j*FILL_VALUE)
        return tmp
    if variable.dimensions[-1] == 'depth' and variable.shape[-1] == 2:
        n_bytes = int(variable.dtype.itemsize*2)
        complex_type = np.dtype('c'+str(n_bytes))
        tmp = variable[:]
        if isinstance(tmp, np.ma.core.MaskedArray):
            # Somehow MaskedArray.view() doesn't work, so convert to a normal
            # array.
            tmp = tmp.filled()
        # Make the data complex
        tmp = tmp.view(dtype=complex_type).squeeze()
        # Turn back into a masked array and re-fill
        tmp = tmp.view(np.ma.MaskedArray)
        tmp.set_fill_value(FILL_VALUE + 1j*FILL_VALUE)
        tmp[tmp == tmp.fill_value] = np.ma.masked
        return tmp
    variable = variable[:]
    variable = variable.view(np.ma.MaskedArray)
    variable.set_fill_value(FILL_VALUE)
    return variable


def get_variable_dimensions(dataset, key):
    '''Get the NetCDF variable dimensons, dealing with complex numbers.

    Key is assumed to be complex if the variable has a first or last dimension
    of 'depth' with length 2. This dimension is removed from the returned list.
    '''
    variable = dataset[key]
    if variable.dimensions[0] == 'depth' and variable.shape[0] == 2:
        return variable.dimensions[1:]
    if variable.dimensions[-1] == 'depth' and variable.shape[-1] == 2:
        return variable.dimensions[:-1]
    return variable.dimensions
