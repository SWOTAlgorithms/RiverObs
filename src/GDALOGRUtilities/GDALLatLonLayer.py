"""
A layer of raster data than can be read with GDAL and interrogated with geographic
coordinates of known projection.
"""

from __future__ import absolute_import, division, print_function

from os.path import join
from osgeo import gdal, gdalconst, osr
from pyproj import Proj, transform
import numpy as np
from numpy.ma import masked_array
from numpy.linalg import inv
from scipy.ndimage.interpolation import map_coordinates


class GDALLatLonLayer:
    """A layer of raster data than can be read with GDAL and interrogated with geographic
    coordinates of known projection.

    Read the data into an array and get the affine transformations.

    Parameters
    -----------

    input_file :
        file containing the raster data
    dtype :
        if None, use the same data type as the input data. Otherwise, convert to dtype.
    destination_projection : str
        a string that can be given to Proj for initialization (default: WGS84 longlat)
        and defines the destination data projection. If None, the destination projection is the same
        as the input projection. This can be used to interpolate the data.
    band: int
        band index to read from the file (index starts at 1)

    Notes
    -------

    If the base layer is lat/lon and all the longitudes > 0, set positive_lon = True
    to enforce that all of the x values lie in [0,360].

    ## The following are GDAL.RasterBand.ReadAsArray keywords:
    ## xoff=0, yoff=0, win_xsize=None, win_ysize=None, buf_xsize=None, buf_ysize=None, buf_obj=None

    """

    def __init__(
            self,
            input_file,
            band=1,
            dtype=np.float32,
            scale=1,
            offset=0,
            destination_projection='+units=m +ellps=WGS84 +datum=WGS84 +proj=longlat ',
            ## xoff=0, yoff=0, win_xsize=None, win_ysize=None, buf_xsize=None, buf_ysize=None, buf_obj=None
            positive_lon=False,
    ):
        self.positive_lon = positive_lon

        gdal.AllRegister()
        self.data_set = gdal.Open(input_file)
        if self.data_set == None:
            raise Exception("Could not open file: %s", input_file)

        # Get the transformation from projection to pixel coordinates

        self.geotransform = self.data_set.GetGeoTransform()
        self.originX = self.geotransform[0]
        self.originY = self.geotransform[3]
        self.pixel_width = self.geotransform[1]
        self.pixel_height = self.geotransform[5]

        # Convert to a true affine matrix that can be inverted, and invert

        self.ij_to_xy_affine = np.array([[
            self.geotransform[1], self.geotransform[2], self.geotransform[0]
        ], [self.geotransform[4], self.geotransform[5], self.geotransform[3]],
                                         [0., 0., 1.]])
        self.xy_to_ij_affine = inv(self.ij_to_xy_affine)

        # Get the projection and the Proj projection

        self.wkt_proj = self.data_set.GetProjection()
        self.spatial_reference = osr.SpatialReference(wkt=self.wkt_proj)
        self.proj4_proj = self.spatial_reference.ExportToProj4()
        self.dataset_proj = Proj(
            self.proj4_proj)  # Projection with this data set

        if destination_projection == None:
            destination_projection = self.proj4_proj

        self.destination_projection = Proj(destination_projection)

        # Read the data into an array

        self.band = self.data_set.GetRasterBand(band)
        self.data = (self.band.ReadAsArray())
        if dtype != None:
            self.data = self.data.astype(dtype)

        self.nodata_value = self.band.GetNoDataValue()
        mask = (self.data == self.nodata_value) | np.isnan(self.data)

        if scale != 1:
            self.data[~mask] *= scale
        if offset != 0:
            self.data[~mask] += offset

        # Fill the masked data with nan, so that any points that use a masked value for interpolation can
        # be filled with nodata_value

        if mask.any():
            self.data = masked_array(self.data, mask=mask)

    def __call__(self,
                 x,
                 y,
                 output=None,
                 nearest=False,
                 order=3,
                 mode='constant',
                 cval=0.0,
                 prefilter=False):
        """Return the value for a set of coordinates stored using scipy.ndimage.interpolation.map_coordinates
        (if nearest=False) or nearest neighbor interpolation (if nearest=True).

        From the scipy documentations:

        The array of coordinates is used to find, for each point in the output,
        the corresponding coordinates in the input. The value of the input at
        those coordinates is determined by spline interpolation of the
        requested order.

        The shape of the output is derived from that of the coordinate
        array by dropping the first axis. The values of the array along
        the first axis are the coordinates in the input array at which the
        output value is found.

        Parameters
        ----------
        x: a tupple,list,ndarray of geographic x coordinates
        y: a tupple,list,ndarray of geographic y coordinates
        output : ndarray or dtype, optional.
        The array in which to place the output, or the dtype of the returned
        array.
        order : int, optional
        The order of the spline interpolation, default is 3.
        The order has to be in the range 0-5.
        mode : str, optional
        Points outside the boundaries of the input are filled according
        to the given mode ('constant', 'nearest', 'reflect' or 'wrap').
        Default is 'constant'.
        cval : scalar, optional
        Value used for points outside the boundaries of the input if
        ``mode='constant'``. Default is 0.0
        prefilter : bool, optional
        The parameter prefilter determines if the input is pre-filtered with
        `spline_filter` before interpolation (necessary for spline
        interpolation of order > 1).  If False, it is assumed that the input is
        already filtered. Default is True.
        """

        if self.positive_lon:
            x = np.where(x < 0, x + 360., x)

        if not nearest:
            self.i, self.j = self.destinationxy_to_ij(
                np.asarray(x), np.asarray(y))

            # Call the interpolation routine

            if type(self.i) != type(np.array([])):  # a scalar

                # map_coordinates accepts only array_like inputs, so convert to array

                return map_coordinates(
                    self.data, ([self.j], [self.i]),
                    output=output,
                    order=order,
                    mode=mode,
                    cval=cval,
                    prefilter=prefilter)[0]
            else:
                return map_coordinates(
                    self.data, (self.j, self.i),
                    output=output,
                    order=order,
                    mode=mode,
                    cval=cval,
                    prefilter=prefilter)

        self.i, self.j = self.destinationxy_to_ij(
            np.asarray(x), np.asarray(y), asint=True)
        return self.data[self.j, self.i]

    def destinationxy_to_datasetxy(self, x, y):
        """Return the data set x,y given the destination x,y"""

        return transform(self.destination_projection, self.dataset_proj, x, y)

    def datasetxy_to_destinationxy(self, dsx, dsy):
        """Return the destination set x,y given the data set x,y"""

        return transform(self.dataset_proj, self.destination_projection, dsx,
                         dsy)

    def datasetxy_to_ij(self, dsx, dsy, asint=False):
        """Get the fractional (or integer if asint=True) pixel values given
        the dataset coordinates."""

        i = self.xy_to_ij_affine[0,
                                 0] * dsx + self.xy_to_ij_affine[0,
                                                                 1] * dsy + self.xy_to_ij_affine[0,
                                                                                                 2]
        j = self.xy_to_ij_affine[1,
                                 0] * dsx + self.xy_to_ij_affine[1,
                                                                 1] * dsy + self.xy_to_ij_affine[1,
                                                                                                 2]

        if asint == True:
            if type(i) != int:
                i = (i + 0.5).astype(np.int32)
                j = (j + 0.5).astype(np.int32)
                i = np.where(i < 0, 0, i)
                j = np.where(i < 0, 0, j)
                i = np.where(i >= self.data.shape[0], self.data.shape[0] - 1,
                             i)
                j = np.where(j >= self.data.shape[0], self.data.shape[1] - 1,
                             j)
            else:
                i = int(i + 0.5)
                j = int(j + 0.5)
                if i < 0: i = 0
                if j < 0: j = 0
                if i >= self.data.shape[0]: i = self.data.shape[0] - 1
                if j >= self.data.shape[1]: j = self.data.shape[1] - 1

        return i, j

    def ij_to_datasetxy(self, i, j):
        """Get the dataset x,y values given the pixel row and column."""

        dsx = self.ij_to_xy_affine[0,
                                   0] * i + self.ij_to_xy_affine[0,
                                                                 1] * j + self.ij_to_xy_affine[0,
                                                                                               2]
        dsy = self.ij_to_xy_affine[1,
                                   0] * i + self.ij_to_xy_affine[1,
                                                                 1] * j + self.ij_to_xy_affine[1,
                                                                                               2]

        return dsx, dsy

    def destinationxy_to_ij(self, x, y, asint=False):
        """Get the fractional (or integer if asint=True) pixel values given
        the destination coordinates."""

        dsx, dsy = self.destinationxy_to_datasetxy(x, y)
        return self.datasetxy_to_ij(dsx, dsy, asint=asint)

    def ij_to_destinationxy(self, i, j):
        """Get the destination x,y values given the pixel row and column."""

        dsx, dsy = self.ij_to_datasetxy(i, j)
        return self.datasetxy_to_destinationxy(dsx, dsy)


class GDALDEMLayer(GDALLatLonLayer):
    """A layer of topography data than can be read with GDAL and interrogated with geographic
    coordinates of known projection. Geoid addition and subtraction is optional."""

    def __init__(
            self,
            input_file,
            band=1,
            dtype=np.float32,
            scale=1,
            offset=0,
            destination_projection='+units=m +ellps=WGS84 +datum=WGS84 +proj=longlat ',
            to_ellipsoid_height=True,
            to_geoid_height=False,
            geoid_file='egm96_15.gtx',
            geoid_dir=None,
            positive_lon=False,
    ):
        """Read the data into an array and get the affine transformations.

        input_file: file containing the raster data

        keyword options:

        destination_projection: a string that can be given to Proj for initialization (default: WGS84 longlat)
        and defines the destination data projection.

        band: band index to read from the file (index starts at 1)

        to_ellipsoid_height: removes the geoid to report height above the reference ellipsoid (default: True)

        to_geoid_height: adds the geoid to report height above the reference geoid (default: False)

        geoid_file: gtx file containing the geoid (needed if ellipsoid_height=True. These files can be obtained from http://download.osgeo.org/proj/vdatum/

        geoid_dir: directory where the gtx geoid file lives.
        """

        GDALLatLonLayer.__init__(
            self,
            input_file,
            band=band,
            dtype=dtype,
            scale=scale,
            offset=offset,
            destination_projection=destination_projection,
            positive_lon=positive_lon)

        # Initialize the geoid from the gtx file

        if geoid_dir != None:
            geoid_file = join(geoid_dir, geoid_file)

        self.geoid = GDALLatLonLayer(
            geoid_file,
            dtype=dtype,
            destination_projection=destination_projection)
        self.to_ellipsoid_height = to_ellipsoid_height
        self.to_geoid_height = to_geoid_height

    def __call__(self,
                 x,
                 y,
                 output=None,
                 nearest=False,
                 order=3,
                 mode='constant',
                 cval=0.0,
                 prefilter=True):
        """Add or subtract geoid to basic GDALLatLonLayer call."""

        h = GDALLatLonLayer.__call__(
            self,
            x,
            y,
            output=output,
            nearest=nearest,
            order=order,
            mode=mode,
            cval=cval,
            prefilter=prefilter)

        geoid = self.geoid(
            x,
            y,
            output=output,
            nearest=nearest,
            order=order,
            mode=mode,
            cval=cval,
            prefilter=prefilter)

        if self.to_ellipsoid_height:
            h += geoid
        elif self.to_geoid_height:
            h -= geoid

        return h
