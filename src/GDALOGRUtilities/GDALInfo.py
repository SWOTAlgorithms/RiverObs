"""
Get relevant information about a GDAL supported raster layer and support various
queries regarding extent via Shapely.
"""

from __future__ import absolute_import, division, print_function

from os.path import join
from osgeo import gdal, gdalconst, osr
from pyproj import Proj, transform
from shapely.geometry import Polygon

# Useful fotmatting functions to print lat/lon in gdalinfo format

def format_lon(lon):
    """Get a longitude string in deg, min, sec."""

    if lon < 0:
        EW = 'W'
    else:
        EW = 'E'
    lon = abs(lon)
    ideg = int(lon)
    fminutes = (lon - ideg)*60
    iminutes = int(fminutes)
    fseconds = (fminutes - iminutes)*60

    return """%dd%2d'%4.2f"%s"""%(ideg,iminutes,fseconds,EW)

def format_lat(lat):
    """Get a latitude string in deg, min, sec."""

    if lat < 0:
        NS = 'S'
    else:
        NS = 'N'
    lat = abs(lat)
    ideg = int(lat)
    fminutes = (lat - ideg)*60
    iminutes = int(fminutes)
    fseconds = (fminutes - iminutes)*60

    return """%dd%2d'%4.2f"%s"""%(ideg,iminutes,fseconds,NS)

class GDALInfo:
    """Get relevant information about a GDAL supported raster layer and support various
    queries regarding extent via Shapely.

    Parameters
    -----------
    input_file : str
        file containing the raster data
    band : int
        band for additional info. Default: 1
    """

    def __init__(self, input_file,band=1):
        gdal.AllRegister()
        self.data_set = gdal.Open(input_file)
        if self.data_set == None:
            raise Exception("Could not open file: %s",input_file)

        # Get the transformation from projection to pixel coordinates

        self.geotransform = self.data_set.GetGeoTransform()
        self.originX = self.geotransform[0]
        self.originY = self.geotransform[3]
        self.pixel_width = self.geotransform[1]
        self.pixel_height = self.geotransform[5]

        self.bands = self.data_set.RasterCount
        self.xsize = self.data_set.RasterXSize
        self.ysize = self.data_set.RasterYSize
        self.band_type = self.data_set.GetRasterBand(1).DataType

        # Get corner locations in native coordinates

        self.ulx = self.geotransform[0]
        self.uly = self.geotransform[3]
        self.lrx = self.ulx + self.geotransform[1] * self.xsize
        self.lry = self.uly + self.geotransform[5] * self.ysize

        self.llx = self.ulx
        self.lly = self.lry
        self.urx = self.lrx
        self.ury = self.uly

        self.centerx = 0.5*(self.ulx + self.lrx)
        self.centery = 0.5*(self.uly + self.lry)

        # Get the projection and the Proj projection

        self.wkt_proj = self.data_set.GetProjection()
        self.spatial_reference = osr.SpatialReference(wkt=self.wkt_proj)
        self.proj4_proj = self.spatial_reference.ExportToProj4()
        self.dataset_proj = Proj(self.proj4_proj) # Projection with this data set

        # This takes you to lon/lat

        destination_projection='+units=m +ellps=WGS84 +datum=WGS84 +proj=longlat '
        self.destination_projection = Proj(destination_projection)

        # Get the lat/lon corners

        self.ullon, self.ullat = transform(self.dataset_proj, self.destination_projection, self.ulx, self.uly)
        self.lrlon, self.lrlat = transform(self.dataset_proj, self.destination_projection, self.lrx, self.lry)
        self.lllon, self.lllat = transform(self.dataset_proj, self.destination_projection, self.llx, self.lly)
        self.urlon, self.urlat = transform(self.dataset_proj, self.destination_projection, self.urx, self.ury)
        self.centerlon, self.centerlat = transform(self.dataset_proj, self.destination_projection, self.centerx, self.centery)

        # Get the lat/lon bounding box

        self.lonmin = min(self.ullon, self.lrlon, self.lllon, self.urlon)
        self.lonmax = max(self.ullon, self.lrlon, self.lllon, self.urlon)
        self.latmin = min(self.ullat, self.lrlat, self.lllat, self.urlat)
        self.latmax = max(self.ullat, self.lrlat, self.lllat, self.urlat)

        # Get the shapely polygons

        self.lonlat_bbox = Polygon([(self.lonmin,self.latmin),
                                    (self.lonmax,self.latmin),
                                    (self.lonmax,self.latmax),
                                    (self.lonmin,self.latmax),
                                    ])

        self.lonlat_poly = Polygon([(self.ullon,self.ullat),
                                    (self.lllon,self.lllat),
                                    (self.lrlon,self.lrlat),
                                    (self.urlon,self.urlat),
                                    ])

        self.xy_poly = Polygon([(self.ulx,self.uly),
                                (self.llx,self.lly),
                                (self.lrx,self.lry),
                                (self.urx,self.ury),
                                ])

        # Get the band nodata value

        self.band = self.data_set.GetRasterBand(band)
        self.nodata_value = self.band.GetNoDataValue()

        # Close the data set and band

        self.band = None
        self.data_set = None

    def print_corners(self,minsec=True):
        """Print the corners information. If minsec, use minutes and seconds format."""

        print("Corner Coordinates:")
        if minsec:
            print("Upper Left  ( %.3f, %.3f) (%s, %s)"%(self.ulx,self.uly,
                                                       format_lon(self.ullon), format_lat(self.ullat)))
            print("Lower Left  ( %.3f, %.3f) (%s, %s)"%(self.llx,self.lly,
                                            format_lon(self.lllon), format_lat(self.lllat)))
            print("Upper Right ( %.3f, %.3f) (%s, %s)"%(self.urx,self.ury,
                                            format_lon(self.urlon), format_lat(self.urlat)))
            print("Lower Right ( %.3f, %.3f) (%s, %s)"%(self.lrx,self.lry,
                                                        format_lon(self.lrlon), format_lat(self.lrlat)))
            print("Center      ( %.3f, %.3f) (%s, %s)"%(self.centerx,self.centery,
                                                    format_lon(self.centerlon), format_lat(self.centerlat)))
        else:
            print("Upper Left  ( %.3f, %.3f) (%.6f, %.6f)"%(self.ulx,self.uly,
                                                       self.ullon, self.ullat))
            print("Lower Left  ( %.3f, %.3f) (%.6f, %.6f)"%(self.llx,self.lly,
                                            self.lllon, self.lllat))
            print("Upper Right ( %.3f, %.3f) (%.6f, %.6f)"%(self.urx,self.ury,
                                            self.urlon, self.urlat))
            print("Lower Right ( %.3f, %.3f) (%.6f, %.6f)"%(self.lrx,self.lry,
                                                        self.lrlon, self.lrlat))
            print("Center      ( %.3f, %.3f) (%.6f, %.6f)"%(self.centerx,self.centery,
                                                    self.centerlon, self.centerlat))

    def llbbox_intersects(self,geom):
        """Does the specified geometry interect the lon/lat bounding box?"""

        return self.lonlat_bbox.intersects(geom)

    def llbbox_contains(self,geom):
        """Does the specified geometry is contained the lon/lat bounding box?"""

        return self.lonlat_bbox.contains(geom)

    def ll_intersects(self,geom):
        """Does the specified geometry interect the lon/lat polygon?"""

        return self.lonlat_poly.intersects(geom)

    def ll_contains(self,geom):
        """Does the specified geometry is contained the lon/lat bounding box?"""

        return self.lonlat_poly.contains(geom)

    def xy_intersects(self,geom):
        """Does the specified geometry interect the lon/lat polygon?"""

        return self.xy_poly.intersects(geom)

    def xy_contains(self,geom):
        """Does the specified geometry is contained the lon/lat bounding box?"""

        return self.xy_poly.contains(geom)
