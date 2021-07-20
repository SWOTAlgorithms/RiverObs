"""
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Dustin Lagoy

"""
from collections import OrderedDict as odict
import copy
import os
import sys
import warnings
import textwrap
import numpy as np
import netCDF4 as nc

import SWOTWater.products.netcdf as netcdf
from SWOTWater.products.constants import FILL_VALUES

FIELD_WARNING = "I'm afraid I can't do that. {} is not in {}"

def textjoin(text):
    """Dedent join and strip text"""
    text = textwrap.dedent(text)
    text = text.replace('\n', ' ')
    text = text.strip()
    return text

def get_subclasses(cls):
    """Recursively get dictionary of names/subclasses of class"""
    products = {}
    for product in cls.__subclasses__():
        products[product.__name__] = product
        products.update(get_subclasses(product))
    return products

def sort_variable_attribute_odict(in_odict):
    """Sorts attributes"""
    # These come first in this order
    FIRST_ATTRS = ['dtype', 'dimensions', 'long_name', 'standard_name',
                   'calendar', 'time', 'standard_time', 'tai_utc_difference',
                   'leap_second']
    # Then put in non-standard ones, and finally these ones in this order
    LAST_ATTRS = [
        'units', 'add_offset', 'scale_factor', 'quality_flag', 'flag_meanings',
        'flag_masks', 'flag_values', 'valid_min', 'valid_max', 'coordinates',
        'comment']

    attr_list = []
    for key in FIRST_ATTRS:
        if key in in_odict:
            attr_list.append([key, in_odict[key]])

    for key in in_odict:
        if key not in (FIRST_ATTRS+LAST_ATTRS):
            attr_list.append([key, in_odict[key]])

    for key in LAST_ATTRS:
        if key in in_odict:
            attr_list.append([key, in_odict[key]])
    return odict(attr_list)

class Product(object):
    """Base class for SWOT-like data products.

    These are basically NetCDF files."""
    # These hold forms of what is allowed to exist in this product
    ATTRIBUTES = odict()
    DIMENSIONS = odict()
    VARIABLES = odict()
    GROUPS = odict()
    ATTRS = [
        'ATTRIBUTES', 'DIMENSIONS', 'VARIABLES', 'GROUPS', '_attributes',
        '_variables', '_groups']

    def __init__(self):
        # These hold what actually exists in memory
        self._attributes = odict()
        self._variables = odict()
        self._groups = odict()
        # reorder the attributes of the variables to the blessed order
        for key, attr_odict in self.VARIABLES.items():
            self.VARIABLES[key] = copy.deepcopy(
                sort_variable_attribute_odict(attr_odict))
        # For bug finding
        assert isinstance(self.ATTRIBUTES, odict)
        assert isinstance(self.DIMENSIONS, odict)
        assert isinstance(self.VARIABLES, odict)
        assert isinstance(self.GROUPS, odict)

    @classmethod
    def get_product(cls, name):
        """Get a product (any subclass of Product()) by it's class name"""
        subclasses = get_subclasses(Product)
        return subclasses[name]

    @property
    def attributes(self):
        dictionary = copy.deepcopy(self._attributes)
        return dictionary

    @property
    def variables(self):
        return odict(
            (key, variable) for key, variable in self._variables.items())

    @property
    def dimensions(self):
        """A dictionary of dimensions and their sizes"""
        # Get default values from DIMENSIONS
        dimensions = odict(
            (key, value) for key, value in self.DIMENSIONS.items())
        # Update with any initialied values
        for key, variable in self.VARIABLES.items():
            for i, vdim in enumerate(variable['dimensions']):
                assert vdim in dimensions, vdim+", "+str(list(dimensions.keys()))
                if key in self.variables:
                    dimensions[vdim] = self[key].shape[i]
        return dimensions

    @property
    def full_attributes(self):
        dictionary = self.attributes
        dictionary.update(self.dimensions)
        return dictionary

    def copy_attributes(self, other_product):
        """Copy the attributes of given product into this one"""
        for key, value in other_product.attributes.items():
            self[key] = value

    def copy(self, **kwargs):
        """Return a deep copy of self, optionally without variable data."""
        # Make a new blank object of this type
        new_product = type(self)()
        # Copy it
        self._copy(new_product, **kwargs)
        return new_product

    def append(self, product):
        """Return a copy of self with values in product appended"""
        new_product = self.copy(with_variables=False)
        for key, value in self._groups.items():
            new_product[key] = value.append(product[key])
        for key, variable in self._variables.items():
            if variable is not None:
                new_product[key] = np.append(variable, product[key], axis=0)
        return new_product

    def _copy(self, new_product, with_variables=True):
        # Copy all of self into new_product
        for key, value in self._attributes.items():
            new_product[key] = value
        for key, value in self._groups.items():
            #print('copy group', key)
            new_product[key] = value.copy(with_variables=with_variables)
        if with_variables:
            for key, variable in self._variables.items():
                if variable is not None:
                    new_product[key] = copy.deepcopy(variable)

    def _copy_from(self, product, with_variables=True):
        # Copy all of self into new_product
        for key in self.ATTRIBUTES:
            self[key] = product[key]
        if with_variables:
            for key in self.VARIABLES:
                self[key] = copy.deepcopy(product[key])
        for key in self.GROUPS:
            print('copy group', key)
            self[key] = self[key]._copy_from(with_variables=with_variables)

    # Override built in functions so calling Product.name or Product['name']
    # return the attribute or variable. We assume attributes and variables have
    # different names.
    def __setattr__(self, key, item):
        if key in self.ATTRS:
            # Allow __init__ to set self up properly
            super(Product, self).__setattr__(key, item)
            return
        if isinstance(item, Product):
            if key in self.GROUPS:
                self._groups[key] = item
                return
        if isinstance(item, np.ndarray):
            # Try to set the variable
            if key in self.VARIABLES:
                form = self.VARIABLES[key]
                # Check the incoming item dimensions against the
                # product-variable ones
                assert len(item.shape) == len(form['dimensions'])
                # Check the product-variable dimensions against the
                # product-global ones
                for i, dimension in enumerate(form['dimensions']):
                    if self.dimensions[dimension] != 0:
                        # Only check if this dimension is already initialized
                        assert self.dimensions[dimension] == item.shape[i]
                self._variables[key] = item
                return
        if key in self.ATTRIBUTES:
            self._attributes[key] = item
            return
        # We couldn't set anything!
        warnings.warn(FIELD_WARNING.format(key, type(self).__name__))

    def __getattr__(self, key):
        if key in self._attributes:
            return self._attributes[key]
        if key in self._variables:
            return self._variables[key]
        if key in self._groups:
            return self._groups[key]
        # Nothing in memory, return empty data
        if key in self.ATTRIBUTES:
            return self.ATTRIBUTES[key].get('value', 'None')
        if key in self.VARIABLES:
            shape = []
            for dimension in self.VARIABLES[key]['dimensions']:
                shape.append(self.dimensions[dimension])
            dtype = self.VARIABLES[key].get('dtype', np.float32)
            quantized_fill = self._getfill(key)
            value = np.ones(shape, dtype=dtype)
            if not value.shape==():
                value[:] = quantized_fill
            return np.ma.masked_array(
                data=value, dtype=dtype, fill_value=quantized_fill,
                mask=np.ones(shape))
        if key in self.GROUPS:
            self[key] = self.get_product(self.GROUPS[key])()
            return self[key]
        raise AttributeError('{} not in product'.format(key))

    @classmethod
    def _getfill(cls, name, dtype=None, is_attribute=False):
        """Returns fill value from cls.VARIABLES or FILL_VALUES"""
        data_dict = cls.VARIABLES
        default_dtype = 'f4'
        if is_attribute:
            data_dict = cls.ATTRIBUTES
            default_dtype = 'str'
        try:
            fill = data_dict[name]['_FillValue']
        except KeyError:
            if dtype is None:
                if 'dtype' in data_dict[name]:
                    dtype = data_dict[name]['dtype']
                else:
                    dtype = default_dtype
            dtype_str = np.dtype(dtype).str[1:]
            if (dtype_str[0] == 'S') or (dtype_str[0] == 'U') or (
                    dtype_str[0]=='O'):
                # handle arbirary-size strings
                dtype_str = 'S1'
            fill = FILL_VALUES[dtype_str]
        return fill

    def __setitem__(self, key, item):
        return setattr(self, key, item)

    def __getitem__(self, key):
        return getattr(self, key)

    def _casted_variable(self, key):
        variable = self[key]
        if 'dtype' in self.VARIABLES[key]:
            dtype = self.VARIABLES[key]['dtype']
        else:
            dtype = variable.dtype
        if isinstance(variable, np.ma.MaskedArray):
            mask = variable.mask
        else:
            mask = False
        quantized_fill = self._getfill(key, dtype)
        return np.ma.masked_array(
            data=variable.astype(dtype), dtype=dtype,
            fill_value=quantized_fill, mask=mask)

    def to_ncfile(self, filename):
        """Write self to a netCDF file."""
        dataset = nc.Dataset(filename, 'w')
        self.to_dataset(dataset)
        dataset.sync()
        dataset.close()

    def to_folder(self, folder):
        """Store self in a folder on disk.

        Will recursively store groups in self in subfolders."""
        os.mkdir(folder)
        for key, value in self._groups.items():
            subfolder = os.path.join(folder, key)
            value.to_folder(subfolder)
        with open(os.path.join(folder, 'attributes.txt'), 'w') as pointer:
            for key, value in self._attributes.items():
                try:
                    pointer.write('{} = {:.17g}\n'.format(key, value))
                except ValueError:
                    pointer.write('{} = {}\n'.format(key, value))
        with open(os.path.join(folder, 'variables.txt'), 'w') as pointer:
            for key, value in self.dimensions.items():
                pointer.write('{} = {:.17g}\n'.format(key, value))
            for key, value in self._variables.items():
                pointer.write('{}: {}\n'.format(key, value.dtype.str))
        for key, variable in self._variables.items():
            if variable is None:
                continue
            if variable.ndim == 1:
                np.savetxt(os.path.join(folder, key), variable)
            else:
                variable.data.tofile(os.path.join(folder, key))

    def to_dataset(self, dataset):
        """Store self in a NetCDF dataset/group.

        Will recursively store groups in self under dataset."""
        for key in self.GROUPS:
            netcdf_group = dataset.createGroup(key)
            # Recursively add the group members
            self[key].to_dataset(netcdf_group)
        for key in self.ATTRIBUTES:
            value = self[key]
            if value is None:
                value = ''
            dataset.setncattr(key, value)
        for key, value in self.dimensions.items():
            dataset.createDimension(key, value)
        for key in self.VARIABLES:
            form = self.VARIABLES[key]
            # Cast the variable to the correct type/fill value
            variable = self._casted_variable(key)
            # Use a helper function to deal with complex numbers
            netcdf.set_variable(
                dataset, key, variable, list(form['dimensions']),
                attributes=form)

    def get_biggest_var(self):
        var_names = [key for key, attr in self.VARIABLES.items()
                     if hasattr(self[key], 'nbytes')]
        var_sizes = [self[key].nbytes for key in var_names]
        biggest = max(zip(var_sizes, var_names))
        return biggest[1], biggest[0]

    @classmethod
    def from_ncfile(cls, filename, variables=None):
        """Generate a product directly from a NetCDF file.

        By default, will only read values defined in this Product.

        If 'force' is True, will read any attribute/variable/dimension, even if
        they are not defined in this Product

        If 'variables' is given, will load only those variables with names in
        the list.
        """
        dataset = nc.Dataset(filename, 'r')
        product = cls()
        product.from_dataset(dataset, variables)
        dataset.close()
        return product


    def from_dataset(self, dataset, variables=None):
        """Load self from a NetCDF dataset/group.

        Will recursively load groups in dataset into self._groups."""
        for key in dataset.ncattrs():
            if key in self.ATTRIBUTES:
                setattr(self, key, dataset.getncattr(key))
            else:
                warnings.warn(FIELD_WARNING.format(key, type(self).__name__))
        for key in netcdf.get_variable_keys(dataset):
            if key in self.VARIABLES:
                # Use a helper function to deal with complex numbers
                variable = netcdf.get_variable(dataset, key)
                setattr(self, key, variable)
            else:
                warnings.warn(FIELD_WARNING.format(key, type(self).__name__))
        for key, group in dataset.groups.items():
            if key in self.GROUPS:
                # Recursively load NetCDF groups into self groups
                class_name = self.GROUPS[key]
                cls = self.get_product(class_name)()
                self[key] = cls
                self[key].from_dataset(dataset[key], variables)
            else:
                raise ValueError('{} not in product'.format(key))

    @classmethod
    def print_xml(cls, prefix=None, ofp=sys.stdout, shape_names=[],
                  shape_dims={}, is_shapefile=False, keys_to_drop=[]):
        """Prints the XML for this data product"""
        def remap_type(dtype):
            return

        klass = cls()

        INDENT = 2*' '
        if prefix is None:
            try:
                uid = klass.UID
            except AttributeError:
                uid = klass.__class__.__name__

            ofp.write(INDENT+'<product>\n')
            ofp.write(2*INDENT+'<science uid="%s">\n'% uid)
            ofp.write(3*INDENT+'<nodes>\n')

        for group in klass.GROUPS:
            if prefix is None:
                next_prefix = group
            else:
                next_prefix = '%s/%s' % (prefix, group)

            klass.get_product(klass.GROUPS[group]).print_xml(
                    prefix=next_prefix, ofp=ofp, shape_names=shape_names,
                    shape_dims=shape_dims, is_shapefile=is_shapefile,
                    keys_to_drop=keys_to_drop)

        for dset in klass.VARIABLES:
            if dset in keys_to_drop:
                continue

            attrs = copy.deepcopy(klass.VARIABLES[dset])
            # get dtype str representation from instance of number
            try:
                type_str = np.dtype(attrs.pop('dtype')).str
            except KeyError:
                type_str = '<f4'

            # _FillValue special handling
            attrs["_FillValue"] = klass._getfill(dset)
            attrs.move_to_end("_FillValue", last=False)

            # hacks for river sp shapefile XMLs
            if is_shapefile:

                # skip these ones
                if dset in ['centerline_lon', 'centerline_lat', 'lat_prior',
                            'lon_prior']:
                    continue

                if dset in ['reach_id', 'node_id', 'rch_id_dn',
                            'rch_id_up', 'time_str', 't_str_avg', 't_str_min',
                            't_str_max', 't_str_med']:
                    attrs['type'] = 'text'
                    if dset in ['reach_id', 'node_id']:
                        attrs.pop("_FillValue", None)
                    else:
                        attrs["fill_value"] = "no_data"
                        attrs.move_to_end("fill_value", last=False)
                        attrs.pop("_FillValue", None)
                else:
                    attrs["fill_value"] = klass._getfill(dset)
                    attrs.move_to_end("fill_value", last=False)
                    attrs.pop("_FillValue", None)

                    # decide if 'text', 'float', 'int4', or 'int9'
                    if type_str[1] == 'f':
                        attrs['type'] = 'float'
                    elif (type_str[1] == 'i' or type_str[1] == 'u'):
                        # hack based on fill value!
                        if attrs['fill_value'] == -999:
                            attrs['type'] = 'int4'
                        else:
                            attrs['type'] = 'int9'
                    else:
                        attrs['type'] = 'text'
                attrs.move_to_end('type', last=False)

            if type_str[1] == 'c':
                type_str = type_str[0] + 'f{}'.format(int(int(type_str[2:])/2))
                attrs['dimensions']['complex_depth'] = 2

            # XML wdith value
            width = int(type_str[2]) * 8

            # XML shape value
            shape_name = "_".join(attrs['dimensions'])+"_shape"

            dims = attrs['dimensions']
            if prefix is not None:
                shape_name = "{}_{}".format(prefix, shape_name)
                dims = odict(
                    ('/{}/{}'.format(prefix, key), value) for key, value in
                    attrs['dimensions'].items())

            if shape_name not in shape_names:
                shape_names.append(shape_name)
                shape_dims[shape_name] = dims

            # Don't write out dimensions
            attrs.pop('dimensions', None)

            # XML node name
            dset_name = '/'+dset if prefix is None else '/%s/%s' % (prefix, dset)

            for name, value in attrs.items():
                if np.iscomplexobj(value):
                    attrs[name] = value.real

            things = []
            for key, value in attrs.items():
                # explicitly use floats for add_offset and scale_factor
                if key in ['add_offset', 'scale_factor']:
                    things.append('%s="%f" ' % (key, value))
                elif key in ['flag_values', ]:
                    things.append('{}="{}"'.format(
                        key, ' '.join(['{}'.format(item) for item in value])))
                else:
                    things.append('{}="{}"'.format(key, value))
            annotations = ' '.join(things)

            # for floats
            if type_str[1] == 'f':
                string = '\n'.join([
                    4*INDENT+'<real name="%s" shape="%s" width="%d">' % (
                        dset_name, shape_name, width),
                    5*INDENT+'<annotation app="conformance" %s/>' % annotations,
                    4*INDENT+'</real>\n'])

            # for integers
            elif type_str[1] == 'i' or type_str[1] == 'u':
                signed_str = 'true' if type_str[1] == 'i' else 'false'
                string = '\n'.join([
                    4*INDENT+'<integer name="%s" shape="%s" width="%d" signed="%s">' % (
                        dset_name, shape_name, width, signed_str),
                    5*INDENT+'<annotation app="conformance" %s/>' % annotations,
                    4*INDENT+'</integer>\n'])

            elif type_str[1:3] == 'S1':
                width = int(type_str[2])
                string = '\n'.join([
                    4*INDENT+'<char name="%s" shape="%s" width="%d">' % (
                        dset_name, shape_name, width),
                    5*INDENT+'<annotation app="conformance" %s/>' % annotations,
                    4*INDENT+'</char>\n'])

            elif type_str[1] in ('S', 'U'):
                width = int(type_str[2])
                string = '\n'.join([
                    4*INDENT+'<string name="%s" shape="%s" width="%d">' % (
                        dset_name, shape_name, width),
                    5*INDENT+'<annotation app="conformance" %s/>' % annotations,
                    4*INDENT+'</string>\n'])

            else:
                raise TypeError("Only ints, floats, strings supported so far!")

            ofp.write(string)

        # add attributes
        for attr, attr_value in klass.ATTRIBUTES.items():
            if attr in keys_to_drop:
                continue

            attr_node = (
                '/@'+attr if prefix is None else '/%s/@%s' % (prefix, attr))

            # get the dtype and doc string
            type_str = np.dtype(attr_value['dtype']).str
            width = int(type_str[2]) * 8

            annotation_desc_string = (
                '<annotation description="%s"/>' % attr_value['docstr']
                if 'docstr' in attr_value else None)

            strings = []
            if type_str[1] == 'f':
                strings += [
                    4*INDENT+'<real name="%s" shape="Scalar" width="%d">' %(
                    attr_node, width),]
                if annotation_desc_string is not None:
                    strings += [5*INDENT+annotation_desc_string,]
                strings += [4*INDENT+'</real>',]

            elif type_str[1] == 'i' or type_str[1] == 'u':
                signed_str = 'true' if type_str[1] == 'i' else 'false'
                strings += [
                    4*INDENT+'<integer name="%s" shape="Scalar" width="%d" signed="%s">' %(
                    attr_node, width, signed_str),]
                if annotation_desc_string is not None:
                    strings += [5*INDENT+annotation_desc_string,]
                strings += [4*INDENT+'</integer>',]

            elif type_str == 'str' or type_str[1] == 'U':
                strings += [
                    4*INDENT+'<string name="%s" shape="Scalar" width="%d">' %(
                    attr_node, width/8),]
                if annotation_desc_string is not None:
                    strings += [5*INDENT+annotation_desc_string,]
                strings += [4*INDENT+'</string>',]
            else:
                raise TypeError("Only ints, floats, strings supported so far!")

            ofp.write('\n'.join(strings)+'\n')

        if prefix is None:
            ofp.write(3*INDENT+'</nodes>\n')
            ofp.write(2*INDENT+'</science>\n')
            ofp.write(2*INDENT+'<shape name="Scalar" order="irrelevant"/>')

            # write out all the dimension shapes
            for shape_name in shape_names:
                string = (2*INDENT +
                    '<shape name="%s" order="slowest...fastest">\n'%shape_name)
                ofp.write(string)

                dims = shape_dims[shape_name]
                for key, value in shape_dims[shape_name].items():
                    string = 3*INDENT+'<dimension extent="%s" name="%s"/>\n' % (
                        value, key)
                    ofp.write(string)
                ofp.write(2*INDENT+"</shape>\n")
            ofp.write(INDENT+"</product>\n")

    def quantize_from(self, my_var, value):
        """
        Sets self[my_var] as the quantized data in data according to
        self.VARIABLES[my_var] scale_factor and dtype
        """
        dtype = self.VARIABLES[my_var]['dtype']
        scale_factor = self.VARIABLES[my_var].get('scale_factor', 1)
        add_offset = self.VARIABLES[my_var].get('add_offset', 0)
        quantized_fill = self._getfill(my_var)
        fill = FILL_VALUES[value.dtype.str[1:]]

        valid = np.logical_and(value != fill, ~np.isnan(value))
        quantized_values = quantized_fill * np.ones(value.shape)

        if dtype[0] in ['i', 'u']:
            # round integers to nearest int versus just typecast
            quantized_values[valid] = np.round(
                (value[valid]-add_offset) / scale_factor).astype(dtype)
        else:
            quantized_values[valid] = (
                (value[valid]-add_offset) / scale_factor).astype(dtype)

        out_value = np.ma.masked_array(
            data=quantized_values, dtype=dtype, fill_value=quantized_fill,
            mask=np.logical_not(valid))
        self[my_var] = out_value

    def requantize(self):
        """Calls quantize_from on all variables that are quantized"""
        for var in self.VARIABLES:
            scale_factor = self.VARIABLES[var].get('scale_factor', 1)
            if scale_factor != 1:
                self.quantize_from(var, self[var])

    def cast(self):
        for group in self.GROUPS:
            self[group].cast()
        for key in self.variables:
            self[key] = self._casted_variable(key)

    def reset_valid_minmax(self):
        """open up valid data range to whole range of datatype"""
        print ('resetting valid min/max range')
        # handle variables not in a group
        for key, variable in self.VARIABLES.items():
            if variable is not None:
                print ('variable: ', key)
                try:
                    info = np.iinfo(variable['dtype'])
                except ValueError:
                    info = np.finfo(variable['dtype'])
                variable['valid_min'] = info.min
                variable['valid_max'] = info.max
        # handle variables in groups
        for groupkey, group in self._groups.items():
            print('group: ', groupkey)
            for key, variable in group.VARIABLES.items():
                if variable is not None:
                    print ('\tvariable: ', key)
                    try:
                        info = np.iinfo(variable['dtype'])
                    except ValueError:
                        info = np.finfo(variable['dtype'])
                    variable['valid_min'] = info.min
                    variable['valid_max'] = info.max

class MutableProduct(Product):
    """A product with forms as instance attributes, not class attributes

    Can be initialized with attributes, dimensions, etc. If arguments are just
    lists:
        dimensions are assumed to be 0
        variables are assumed to have all dimensions
        groups are assumed to be MutableProduct()s
    """
    def __init__(
            self, attributes=None, dimensions=None, variables=None,
            groups=None):
        super().__init__()
        # Override these with instance attributes, so they can be specific to
        # this instance
        self.ATTRIBUTES = []
        self.DIMENSIONS = odict()
        self.VARIABLES = odict()
        self.GROUPS = odict()
        # Set whatever has been given
        if attributes is not None:
            self.ATTRIBUTES.extend(attributes)
        if isinstance(dimensions, list):
            self.DIMENSIONS.update([[name, 0] for name in dimensions])
        elif isinstance(dimensions, dict):
            self.DIMENSIONS.update(dimensions)
        if isinstance(variables, list):
            self.VARIABLES.update([[name, odict()] for name in variables])
            for name, form in self.VARIABLES.items():
                form['dimensions'] = self.DIMENSIONS
        elif isinstance(variables, dict):
            self.VARIABLES.update(variables)
        if isinstance(groups, list):
            self.GROUPS.update([[name, 'MutableProduct'] for name in groups])
        elif isinstance(groups, dict):
            self.GROUPS.update(groups)

    def from_dataset(self, dataset, variables=None):
        """Just like Product()'s version, but set self's forms from the file"""
        self._form_from_dataset(dataset)
        super().from_dataset(dataset, variables)

    def _form_from_dataset(self, dataset):
        """Make a MutableProduct with forms set from dataset"""
        for key in dataset.ncattrs():
            self.ATTRIBUTES.append(key)
        for key in netcdf.get_variable_keys(dataset):
            variable = netcdf.get_variable(dataset, key)
            variable_dimensions = netcdf.get_variable_dimensions(
                dataset, key)
            dimensions = odict()
            for i, dimension in enumerate(variable_dimensions):
                dimensions[dimension] = variable.shape[i]
            self.VARIABLES[key] = odict([
                ['dtype', variable.dtype],
                ['dimensions', odict(
                    [[dimension, 0] for dimension in variable_dimensions])],
            ])
            for attribute in dataset[key].ncattrs():
                self.VARIABLES[key][attribute] = dataset[key].getncattr(
                    attribute)
            for dimension in variable_dimensions:
                if dimension not in self.DIMENSIONS:
                    self.DIMENSIONS[dimension] = 0
        for key in dataset.groups:
            self.GROUPS[key] = 'MutableProduct'

    def _copy(self, new_product, with_variables=True):
        new_product.ATTRIBUTES = self.ATTRIBUTES
        new_product.DIMENSIONS = self.DIMENSIONS
        new_product.VARIABLES = self.VARIABLES
        new_product.GROUPS = self.GROUPS
        super()._copy(new_product, with_variables)
