from __future__ import absolute_import, print_function

from .version import __version__
try:
    from .OGR2Shapely import ShapelyDataSource, ShapelyLayer, ShapelyFeature
    from .OGRWriter import OGRWriter
    from .GDALInfo import GDALInfo
    from .GDALutilities import WarpToLayer
except:
    print(
        "Warning: ShapelyDataSource, ShapelyLayer, ShapelyFeature, OGRWriter disabled."
    )
    print("Warning: GDALInfo, WarpToLayer disabled.")
    print("Please install shapely if you would like to use the classes.")
from .GDALLatLonLayer import GDALLatLonLayer, GDALDEMLayer
from .GeodeticPath import GeodeticPath, GeodeticPathFromPegPoint
from .GDALWriter import write_llh_to_gdal, write_numpy_to_gdal
from .CoordinateTransformations import CoordinateTransformation
