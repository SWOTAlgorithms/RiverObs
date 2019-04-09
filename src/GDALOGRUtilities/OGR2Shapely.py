"""
Read OGR supported vectors from/to shapely arrays.
"""

from __future__ import absolute_import, division, print_function

from osgeo import ogr, osr
import numpy as np
from shapely import wkt


class ShapelyDataSource:
    """Hold all of the information in an ogr data source as shapely arrays.

    Initilize and, optionally, read ogr data from an ogr supported file
    or an ogr DataSource.
    """

    def __init__(self, ogr_file=None, ogr_data_source=None):
        # The layers are stored as a dictionary of ShapelyLayers

        self.nlayers = 0
        self.layer = {}
        self.name = None

        # This sets the layer to use when iterating over the data source

        self.layer_index = 0

        # Read from a file or an OGR data source, if available

        if ogr_file != None:
            self.from_ogr_file(ogr_file)
        elif ogr_data_source != None:
            self.from_ogr_data_source(ogr_data_source)
        else:
            self.ogr_file = None
            self.ogr_data_source = None

    def from_ogr_file(self, ogr_file):
        """Initialize from an ogr_supported data file."""

        self.ogr_file = ogr_file
        ogr_data_source = ogr.Open(ogr_file)
        if ogr_data_source == None:
            raise Exception('Cannot open ogr file: %s' % ogr_file)
        self.from_ogr_data_source(ogr_data_source)

    def from_ogr_data_source(self, ogr_data_source):
        """Initialize from an open ogr DataSource."""

        self.ogr_data_source = ogr_data_source
        self.name = ogr_data_source.GetName()
        self.nlayers = ogr_data_source.GetLayerCount()

        for i in range(self.nlayers):
            ogr_layer = ogr_data_source.GetLayer(i)
            self.layer[i] = ShapelyLayer(ogr_layer=ogr_layer)

    def __getitem__(self, index):
        """Return the shape in the index feature in the active layer. This is useful for iterating over the shapes in the active layer."""

        return self.layer[self.layer_index].feature[index].shape

    def get_shape(self, index, layer_index=0):
        """Return the shape in the index feature."""

        return self.layer[layer_index].feature[index].shape

    def get_numpy_shape(self, index, layer_index=0):
        """Return the shape in the index feature as a numpy array."""

        return np.array(self.layer[layer_index].feature[index].shape)

    def get_field(self, index, layer_index=0):
        """Return the field in the index feature."""

        return self.layer[layer_index].feature[index].field

    def get_shapes(self, layer_index=0):
        """Return the all the shapes in the layer features."""

        return [
            self.layer[layer_index].feature[i].shape
            for i in range(self.layer[layer_index].nfeatures)
        ]

    def get_numpy_shapes(self, layer_index=0):
        """Return the all the shapes in the layer features."""

        return np.array([
            np.array(self.layer[layer_index].feature[i].shape)
            for i in range(self.layer[layer_index].nfeatures)
        ])

    def get_fields(self, layer_index=0):
        """Return the all the fields in the layer features."""

        return [
            self.layer[layer_index].feature[i].field
            for i in range(self.layer[layer_index].nfeatures)
        ]


class ShapelyLayer:
    """Holds layer information in a dictionary of ShapelyFeatures."""

    def __init__(self, ogr_layer=None):
        """Initialize and, optionally copy an ogr layer."""

        self.nfeatures = 0
        self.feature = []
        self.name = ''
        ## self.geometries = [] # a list containing the geometry for each feature
        ## self.fields = [] # a list containing the fields for each feature

        if ogr_layer != None: self.from_ogr(ogr_layer)

    def from_ogr(self, ogr_layer):
        """Initialize from an OGR layer."""

        self.ogr_layer = ogr_layer
        self.spatial_ref = ogr_layer.GetSpatialRef()
        self.nfeatures = ogr_layer.GetFeatureCount()
        self.name = ogr_layer.GetName()

        for i in range(self.nfeatures):  # Notice that the count starts at 1!
            ogr_feature = ogr_layer.GetNextFeature()
            self.feature.append(ShapelyFeature(ogr_feature=ogr_feature))
            ## self.geometries.append(self.feature[i].shape)
            ## self.fields.append(self.feature[i].field)

    def __getitem__(self, index):
        """Return the shape in the index feature. This is useful for iterating over the shapes in the layer."""

        return self.feature[index].shape

    def get_shape(self, index):
        """Return the shape in the index feature."""

        return self.feature[index].shape

    def get_numpy_shape(self, index):
        """Return the shape in the index feature as a numpy array."""

        return np.array(self.feature[index].shape)

    def get_field(self, index):
        """Return the field in the index feature."""

        return self.feature[index].field

    def get_shapes(self):
        """Return the all the shapes in the layer features."""

        return [self.feature[i].shape for i in range(self.nfeatures)]

    def get_numpy_shapes(self):
        """Return the all the shapes in the layer features."""

        return np.array(
            [np.array(self.feature[i].shape) for i in range(self.nfeatures)])

    def get_fields(self):
        """Return the all the fields in the layer features."""

        return [self.feature[i].field for i in range(self.nfeatures)]


class ShapelyFeature:
    """Holds feature information in shapely arrays."""

    def __init__(self, ogr_feature=None):
        """Initilaize from an ogr layer or a shapely shape"""

        self.nfields = 0
        self.field = {}

        if type(ogr_feature) != type(
                None):  # awkward for direct comparison is broken
            self.from_ogr(ogr_feature)

    def from_ogr(self, ogr_feature):
        """Initialize from an OGR feature."""

        self.ogr_feature = ogr_feature
        self.nfields = ogr_feature.GetFieldCount()

        for i in range(self.nfields):
            self.field[list(ogr_feature.keys())[i]] = ogr_feature.GetField(i)

        # Now get the geometry

        self.geometry = ogr_feature.GetGeometryRef()
        self.geometry_name = self.geometry.GetGeometryName()
        self.npoints = self.geometry.GetPointCount()

        self.shape = wkt.loads(self.geometry.ExportToWkt())
