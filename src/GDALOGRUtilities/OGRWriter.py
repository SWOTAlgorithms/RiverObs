"""
A Class to write vector files supported by OGR with inputs from
various sources.
"""

from os.path import split,splitext
from osgeo import ogr, osr
from shapely.geometry import asShape, asPoint, asLineString
from collections import OrderedDict

class OGRWriter:
    """A Class to write vector files supported by OGR with inputs from
    various sources.

    Open the file for writing and intialialize the fields.

    Parameters
    -----------
    
    output_file : str
        Name of the output file.
    layers : a list of layer names.
        If empty, the root name of the file is used.
        If a string is passed, then only a single layer with that name is created.
    fields :
        a list of dictionaries, one per field, containing the options for each layer.
        The dictionaries are of the form

        field_dict[field_name] = type

        or

        field_dict[field_name] = (type,length)

        type is one of the accepted types in field_types.
        The optional length is the field length ('str') or precision ('float').
        If a single dictionary is passed, it will be used for all layers.

    diver : str
        one of the driver names accepted by OGR.

    coordinate_system :
        one of the well known coordinate systems accepted by OGR.
        Currently supported:
        "WGS84": same as "EPSG:4326" but has no dependence on EPSG data files.
        "WGS72": same as "EPSG:4322" but has no dependence on EPSG data files.
        "NAD27": same as "EPSG:4267" but has no dependence on EPSG data files.
        "NAD83": same as "EPSG:4269" but has no dependence on EPSG data files.
        "EPSG:n": same as doing an ImportFromEPSG(n).

        proj4_string: a string understood by proj4 (optional). For instance:
        proj4_string = '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs '
        If not empty, this is applied after setting the well know coordinate system.

    geometry : str
        one of the keys in OGRWriter.geom_type dictionary

    srs :
        a predifined osr spatial reference. Overrides other srs options.
    """

    # This dictionary translates between various ways of specifying geometries
    # and OGR supported geometries
    
    geom_type = {
        'Point':ogr.wkbPoint, 	 # 0-dimensional geometric object, standard WKB
        'point':ogr.wkbPoint, 	 # 0-dimensional geometric object, standard WKB
        'LineString': ogr.wkbLineString, 	# 1-dimensional geometric object with linear interpolation between Points, standard WKB
        'line': ogr.wkbLineString, 	# 1-dimensional geometric object with linear interpolation between Points, standard WKB
        'Polygon': ogr.wkbPolygon, #planar 2-dimensional geometric object defined by 1 exterior boundary and 0 or more interior boundaries, standard WKB
        'polygon': ogr.wkbPolygon, #planar 2-dimensional geometric object defined by 1 exterior boundary and 0 or more interior boundaries, standard WKB
        ## wkbMultiPoint 	 GeometryCollection of Points, standard WKB
        ## wkbMultiLineString 	 GeometryCollection of LineStrings, standard WKB
        ## wkbMultiPolygon 	 GeometryCollection of Polygons, standard WKB
        ## wkbGeometryCollection 	 geometric object that is a collection of 1 or more geometric objects, standard WKB
        ## wkbNone 	 non-standard, for pure attribute records
        ## wkbLinearRing 	 non-standard, just for createGeometry()
        ## wkbPoint25D 	 2.5D extension as per 99-402
        ## wkbLineString25D 	 2.5D extension as per 99-402
        ## wkbPolygon25D 	 2.5D extension as per 99-402
        ## wkbMultiPoint25D 	 2.5D extension as per 99-402
        ## wkbMultiLineString25D 	 2.5D extension as per 99-402
        ## wkbMultiPolygon25D 	 2.5D extension as per 99-402
        ## wkbGeometryCollection25D 	 2.5D extension as per 99-402
        }

    field_type = {
        'int':ogr.OFTInteger, 	# Simple 32bit integer
        ## OFTIntegerList 	List of 32bit integers
        'float':ogr.OFTReal, 	# Double Precision floating point
        ## OFTRealList 	List of doubles
        'str':ogr.OFTString, 	# String of ASCII chars
        ## OFTStringList 	Array of strings
        ## OFTBinary 	Raw Binary data
        ## OFTDate 	Date
        ## OFTTime 	Time
        ## OFTDateTime 	Date and Time
        }

    def __init__(self,output_file,layers=[],fields={},driver='ESRI Shapefile',
                 geometry='Point',coordinate_system='WGS84',proj4_string='',srs=None):

        # Set the coordinate system

        if srs != None:
            self.srs = srs
        else:
            self.srs = osr.SpatialReference()
            self.srs.SetWellKnownGeogCS(coordinate_system)
            if proj4_string != '':
                self.srs.ImportFromProj4(proj4_string)

        # Open the output file
 
        self.driver = ogr.GetDriverByName(driver)
        if self.driver is None:
            raise Exception("%s driver not available" % driver)

        self.data_source = self.driver.CreateDataSource(output_file)
        if self.data_source is None:
            raise Exception("Creation of output file %s failed"%output_file)

        # Create all the layers

        if layers == []: # Create the layer with the file name as the 
            layers = [splitext(split(output_file)[-1])[0]]
        elif type(layers) == type(''):
            layers = [layers]

        self.layer_names = layers
        self.nlayers = len(layers)

        # make sure the fields are consistent with the layers

        if (type(fields) == type({})) or (type(fields) == type(OrderedDict())) :
            fields = [fields]*self.nlayers

        if len(fields) != self.nlayers:
            raise Exception('Field and layer lengths are inconsistent.')

        self.geometry_type = self.geom_type[geometry]
        self.geometry = geometry
        self.layers = []
        for i,layer in enumerate(layers):
            lyr = self.data_source.CreateLayer(layer,srs=self.srs,
                                               geom_type=self.geometry_type)
            if lyr is None:
                raise Exception('Could not create layer: %s'%layer)

            self.init_fields(lyr,fields[i])
            
            self.layers.append(lyr)

    def close(self):
        """Close the data source to finish flushing the data."""
        self.data_source.Destroy()

    def init_fields(self,layer,fields):
        """Add the appropriate fields to an open layer."""

        for field, value in fields.items():

            # Get the field definitions

            if type(value) == type(''):
                field_type = self.field_type[value]
                n = 0
            else:
                field_type = self.field_type[value[0]]
                n = value[1]

            field_defn = ogr.FieldDefn(field,field_type)
            if n > 0:
                if value[0] == 'str':
                    field_defn.SetWidth(n)
                elif value[0] == 'float':
                    field_defn.SetPrecison(n)

            # Create the fields

            if layer.CreateField(field_defn) != 0:
                raise Exception('Cannot create field: %s'%field)

    def add_wkt_feature(self,wkt,field_record,layer_index=0,layer_name=None):
        """Add a feature from a WKT (well known text) string and a field dictionary.

        Parameters
        ----------

        layer_index : int
            index of the layer to which the feature will be added.
        layer_name : str
            overrides layer index 
        """

        if layer_name != None:
            layer_index = self.layer_names.index(layer_name)

        layer = self.layers[layer_index]
        feature = ogr.Feature(layer.GetLayerDefn())

        for name, value in field_record.items():
            feature.SetField(name,value)

        geometry = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(geometry)

        if layer.CreateFeature(feature) != 0:
            raise Exception('Failed to create feature')

        feature.Destroy()
                    
    def add_geo_feature(self,object,field_record,layer_index=0,layer_name=None):
        """Add a feature from an object that supports the Python Geo interface and
        can be converted to a Shape object by shapely.geometry.asShape.

        Parameters
        -----------

        object : object
            This can be a GeoJSON python dictionary; e.g.,

            object = {"type": "Point", "coordinates": (0.0, 0.0)}

        layer_index : int
            index of the layer to which the feature will be added.
        layer_name: str
            overrides layer index 
        """

        shape = asShape(object)
        self.add_wkt_feature(shape.wkt,field_record,layer_index=layer_index,layer_name=layer_name)

    def add_xy_feature(self,x,y,field_record,layer_index=0,layer_name=None):
        """Add a feature from x,y numpy arrays (or anything that can be zipped into coordinate pairs).
        
        layer_index: index of the layer to which the feature will be added.
        layer_name: overrides layer index 
        """

        coords = zip(x,y)
        if self.geometry in ['Point','point']:
            shape = asPoint(coords)
        elif self.geometry in ['line','LineString']:
            shape = asLineString(coords)
        elif self.geometry in ['polygon','Polygon']: # Right now, only simple polygons are supported
            shape = asPolygon(coords)

        self.add_wkt_feature(shape.wkt,field_record,layer_index=layer_index,layer_name=layer_name)
            
