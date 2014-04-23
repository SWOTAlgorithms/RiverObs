"""
Interfaces to the gdal utility programs, http://www.gdal.org/gdal_utilities.html.
"""

import shlex
from subprocess import Popen, PIPE
from GDALOGRUtilities import GDALInfo

class WarpToLayer:
    """Given a reference gdal layer, extract data from other files,
    reproject, if required, using gdalwarp and extract the data corresponding to
    the reference file at the reference file resolution."""

    def __init__(self,reference_layer_file):
        """Initialize with a reference layer file."""

        self.reference_layer = GDALInfo(reference_layer_file)

    def __call__(self,srcfile,dstfile,dstnodata_value=None,executable='gdalwarp',
                 resampling_method='bilinear',gdal_format='GTiff'):
        """Warp srcfile to dstfile (same projection, window, and resolution).

        resampling_method: Resampling method to use. Available methods are:
        'near': nearest neighbour resampling (default, fastest algorithm, worst interpolation quality).
        'bilinear': bilinear resampling.
        'cubic': cubic resampling.
        'cubicspline': cubic spline resampling.
        'lanczos': Lanczos windowed sinc resampling.
        'average': average resampling, computes the average of all non-NODATA contributing pixels. (GDAL >= 1.10.0)
        'mode': mode resampling, selects the value which appears most often of all the sampled points. (GDAL >= 1.10.0)
        """

        self.init_source_layer(srcfile)
        self.warp(dstfile,dstnodata_value=dstnodata_value,executable=executable,
                  resampling_method=resampling_method,gdal_format=gdal_format)

    def init_source_layer(self,source_layer_file):
        """Get the gdal information from the source layer file, and decide whether 
        reprojection is necessary."""

        self.source_layer_file = source_layer_file
        self.source_layer = GDALInfo(source_layer_file)

        # Check to see if the projections are the same

        self.same_projection = self.reference_layer.spatial_reference.IsSame(self.source_layer.spatial_reference) 

    def warp(self,dstfile,dstnodata_value=None,executable='gdalwarp',
             resampling_method='bilinear',gdal_format='GTiff'):
        """Warp the input layer to the same coordinate system as the reference layer
        by calling gdalwarp.

        resampling_method: Resampling method to use. Available methods are:
        'near': nearest neighbour resampling (default, fastest algorithm, worst interpolation quality).
        'bilinear': bilinear resampling.
        'cubic': cubic resampling.
        'cubicspline': cubic spline resampling.
        'lanczos': Lanczos windowed sinc resampling.
        'average': average resampling, computes the average of all non-NODATA contributing pixels. (GDAL >= 1.10.0)
        'mode': mode resampling, selects the value which appears most often of all the sampled points. (GDAL >= 1.10.0)

        gdal_format: one of the formats supported by gdalwarp (e.g., gdalwarp --formats).
        """

        self.warp_file = dstfile
        if dstnodata_value == None:
            dstnodata_value = self.source_layer.nodata_value

        proj4 = self.reference_layer.proj4_proj
        srcfile = self.source_layer_file
        
        warp_command = "%(executable)s  -t_srs '%(proj4)s' -r %(resampling_method)s "

        # Set the no data values

        if self.source_layer.nodata_value != None:
            srcnodata_value = self.source_layer.nodata_value
            warp_command += " -srcnodata %(srcnodata_value)s "

        if dstnodata_value != None:
            warp_command += " -dstnodata %(dstnodata_value)s "

        # Set the destination pixel size

        xres = abs(self.reference_layer.pixel_width)
        yres = abs(self.reference_layer.pixel_height)        
        warp_command += " -tr %(xres)s %(yres)s "

        # Set the destination output window
        
        xmin = self.reference_layer.llx
        ymin = self.reference_layer.lly
        xmax = self.reference_layer.urx
        ymax = self.reference_layer.ury

        warp_command += " -te %(xmin)s %(ymin)s %(xmax)s %(ymax)s "

        # Set the destination file and format

        warp_command += " -of %(gdal_format)s %(srcfile)s %(dstfile)s"
        warp_command = warp_command%locals()

        print warp_command
        args = shlex.split(warp_command)
        p = Popen(args,
                  shell=False,
                  stdout=PIPE,
                  stderr=PIPE)

        self.warp_stdout, self.warp_stderr = p.communicate()

