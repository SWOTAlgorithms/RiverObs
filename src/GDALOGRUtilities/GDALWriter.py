"""
Write gdal files from various inputs.
"""

from __future__ import absolute_import, division, print_function

from os.path import join
import numpy as np
from numpy.ma import masked_array, is_masked
from osgeo import gdal
import osr

def get_gdal_type(dtype):
    """Given a numpy data type, return the corresponding """

    if (dtype == np.int8) or (dtype == np.uint8):
        gdal_type = gdal.GDT_Byte
    elif dtype == np.int16:
        gdal_type = gdal.GDT_Int16
    elif dtype == np.uint16:
        gdal_type = gdal.GDT_UInt16
    elif dtype == np.int32:
        gdal_type = gdal.GDT_Int32
    elif dtype == np.uint32:
        gdal_type = gdal.GDT_UInt32
    elif dtype == np.float32:
        gdal_type = gdal.GDT_Float32
    elif dtype == np.float64:
        gdal_type = gdal.GDT_Float64
    else:
        raise Exception('Unknown dtype: %s'%dtype)
    return gdal_type

def write_llh_to_gdal(llh_data,lon_min,dlon,lat_min,dlat,
                      gdal_format, dst_filename, origin_up=True,
                      options=None,nodata_value=None,vflip_data=False):
    """Write an LLH layer to a GIS file in a gdal supported format.

    vflip_data: if True llh_data => llh_data[::-1,:]. Use in case the data
    is not aligned with the desired geotranform.
    """

    gdal_type = get_gdal_type(llh_data.dtype)

    # Get the driver and open the output file

    driver = gdal.GetDriverByName( gdal_format )
    if driver == None:
        raise Exception('Unimplented gdal driver: %s'%driver)

    dst_ds = driver.Create( dst_filename, llh_data.shape[1], llh_data.shape[0],
                            bands=1, eType = gdal_type) #, options=options )

    # Flip the data if needed to be consistent with the geotransform

    if vflip_data:
        llh_data = llh_data[::-1,:]

    # Set all of the transform information

    if origin_up:
        nlat = llh_data.shape[0]
        lat_max = lat_min + (nlat -1)*dlat
        dst_ds.SetGeoTransform( [ lon_min, dlon, 0,
                                  lat_max, 0, -dlat ] )
    else:
        dst_ds.SetGeoTransform( [ lon_min, dlon, 0,
                                      lat_min, 0, dlat ] )
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS( 'WGS84' )
    dst_ds.SetProjection( srs.ExportToWkt() )

    # Now write the raster

    band = dst_ds.GetRasterBand(1)

    if nodata_value != None:
        band.SetNoDataValue(nodata_value)

    if is_masked(llh_data):
        if nodata_value != None:
            llh_data.data[llh_data.mask] = nodata_value
        band.WriteArray( llh_data.data )
    else:
        band.WriteArray( llh_data )

    # Clean up by closing the dataset

    dst_ds = None
    src_ds = None

def write_numpy_to_gdal(data,geotransform,wkt_proj,
                        dst_filename, gdal_format='GTiff',origin_up=True,
                        options=None,nodata_value=None):
    """Given numpy data and projection information, write to a gdal file.

    Parameters
    ----------

    data :
        a 2D numpy array
    geotransform :
        a list containing the affine transformation
        (e.g., the result of gdal data_set.GetGeoTransform())
    wkt_proj :
        well known text projection information
        (e.g., the data_set.GetProjection() )
    dst_filename : str
        destination file name
    origin_up : bool
        if origin_up == True, the data is reversed in its first axis
    option :
        options to pass to gdal.
    nodata_value :
        nodata_value value. If None, no nodata_value value is set.
    """

    gdal_type = get_gdal_type(data.dtype)

    # Get the driver and open the output file

    driver = gdal.GetDriverByName( gdal_format )
    if driver == None:
        raise Exception('Unimplented gdal driver: %s'%driver)

    dst_ds = driver.Create( dst_filename, data.shape[1], data.shape[0],
                            bands=1, eType = gdal_type) #, options=options )

    # Set all of the transform information

    if origin_up:
        data = data[::-1,:]

    dst_ds.SetGeoTransform( geotransform )
    dst_ds.SetProjection( wkt_proj )

    # Now write the raster

    band = dst_ds.GetRasterBand(1)

    if nodata_value != None:
        band.SetNoDataValue(nodata_value)

    if is_masked(data):
        if nodata_value != None:
            data.data[data.mask] = nodata_value
        band.WriteArray( data.data )
    else:
        band.WriteArray( data.data )

    # Clean up by closing the dataset

    dst_ds = None
    src_ds = None
