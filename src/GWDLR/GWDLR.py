"""
A class to read GWDLR files, make masks, and output to GIS raster data.

The files are in GrADS format, described in http://www.iges.org/grads/gadoc/aboutgriddeddata.html.
A simple GrADS parser is implemented here.
"""

from __future__ import absolute_import, division, print_function

from os.path import join
import numpy as np
from numpy.ma import masked_array
from skimage.morphology import skeletonize

from GDALOGRUtilities import write_llh_to_gdal

class GWDLR:
    """A class to read GWDLR files, make masks, and output to GIS raster data."""

    def __init__(self,rootname,data_dir='.'):
        """Read the bin and ctl files, and parse the ctl."""

        self.ctl_file = join(data_dir,rootname+'.ctl')
        self.parse_ctl(self.ctl_file)

        self.bin_file = join(data_dir,rootname+'.bin')
        self.data = np.reshape(np.fromfile(self.bin_file,dtype=np.float32),
                                (self.y.size,self.x.size))

    def parse_ctl(self,ctl_file):
        """Parse the ctl file to extract variables and store as object variables."""

        # Read the ctl file and close

        with open(ctl_file) as fin:
            header = fin.readlines()

        # parse the header, line by line

        for line in header:
            s = line.split()
            if not len(s) > 0:
                continue
            kwd = s[0]
            if kwd in ['xdef','ydef','tdef','zdef']:
                self.parse_def(line)
            elif kwd in ['dset','title','options']:
                #exec('self.%s = "%s"'%(s[0],s[1]))
                setattr(self,s[0],s[1])
            elif kwd in ['undef','vars']:
                #exec('self.%s = %s'%(s[0],s[1]))
                setattr(self,s[0],s[1])
            elif 'endian' in line:
                for wd in s:
                    if 'endian' in wd:
                        #exec('self.endian = "%s"'%wd)
                        setattr(self,'endian',wd)

    def parse_def(self,line):
        """Store the def into a an object with the following attributes:
        name, size, scale, origin, step."""

        # An empty class to hold the variables

        class Axis:
            pass

        # Parse the line and store in the class

        kwd, size, scale, origin, step = line.split()
        axis = Axis()
        axis.name = kwd[0]
        axis.size = int(size)
        axis.scale = scale
        try:
            axis.origin = float(origin)
        except:
            axis.origin = origin
        try:
            axis.step = float(step)
        except:
            axis.step = step
        if axis.size > 1:
            axis.coordinates = axis.origin + np.arange(axis.size)*axis.step

        # Store the Axis instance as a member

        #exec('self.%s = axis'%(axis.name))
        setattr(self,axis.name,axis)

    def to_mask(self,width,overwrite=True,thin=False):
        """Return a mask with 1's when the data is >= width, 0 otherwise.
        If overwrite == True, the self.data is replaced with the mask.

        If skeletonize: thin the mask
        """

        mask = (self.data >= width).astype(np.uint8)

        if thin:
            mask = skeletonize(mask).astype(np.uint8)

        if overwrite:
            self.data = mask

        return mask

    def to_gdal(self,dst_filename, gdal_format='GTiff',nodata_value=None):
        """Write the data to a georeferenced GIS file compatible with gdal."""

        lon_min,dlon,lat_min,dlat = ( self.x.coordinates.min(), self.x.step,
                                      self.y.coordinates.min(), self.y.step )

        write_llh_to_gdal(self.data,lon_min,dlon,lat_min,dlat,
                      gdal_format, dst_filename, origin_up=True,
                      options=None,nodata_value=nodata_value,vflip_data=False)
