"""
Access SWOT L2 data conveniently and provides LatLonRegion protocol.
"""

import numpy as N
from netCDF4 import Dataset
from pyproj import Proj

class SWOTL2:
    """Access SWOT L2 data conveniently. SWOTL2 implements the LatLonRegion object
    interfaces in that it provides a  bounding_box member and a proj function that goes
    from (lon,lat) -> (x,y) and backwards.

    Parameters
    ----------
    
    bounding_box : tuple or array_like
        should be of the form (lonmin,latmin,lonmax,latmax).
        If the bounding_box is provided, select only the data in the
        bounding box, otherwise, get the bounding box from the data
        itself.
    good_class : list or array_like
        a list of class labels where the data are returned.
    class_list : list
        a list of the class labels for what is considered good data
    lon_kwd, lat_kwd : str
        netcdf names for the longitudes and latitudes to be used
        for georeferencing the data set.
    class_kwd : str
        the netcdf name of the classification layer to use.
    min_points : int
        If the number of good points is less than this, raise an exception.
    verbose : bool, default False
        If True, print to sdtout while processing.

    Notes
    ------

    The final set of keywords are projection options for pyproj.
    A full list of projection options to set plus explanations of their
    meaning can be found here: https://trac.osgeo.org/proj/wiki/GenParms
        
    The default projection is Lambert equiarea, which has a proj4 string with the
    following parameters:
        
    +proj=laea  
    +lat_0=Latitude at projection center, set to bounding box center lat 
    +lon_0=Longitude at projection center, set to bounding box center lon
    +x_0=False Easting, set to 0
    +y_0=False Northing, set to 0
    """

    azimuth_spacing_default = 5.3 # default azimuth spacing

    def __init__(self,swotL2_file,bounding_box=None,class_list=[1],
                 lat_kwd='no_layover_latitude', lon_kwd='no_layover_longitude',
                 class_kwd='no_layover_classification',
                 min_points=100,
                 project_data=True,
                 verbose=False,
                proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,
                ellps='WGS84',**proj_kwds):

        self.verbose = verbose
        
        self.lat_kwd, self.lon_kwd = lat_kwd, lon_kwd
        self.nc = Dataset(swotL2_file)
        if self.verbose: print('Dataset opened')

        # Get some of the metadata

        try:
            self.azimuth_spacing = float(self.nc.azimuth_spacing)
        except:
            self.azimuth_spacing = self.default_azimuth_spacing
            
        self.set_index_and_bounding_box(bounding_box,lat_kwd,lon_kwd,
                                        class_list,class_kwd=class_kwd)
        if self.verbose: print('Good data selected & bounding box calculated.')

        # Get locations for these data (note that these should be the reference values)

        self.lat = self.get(lat_kwd)
        self.lon = self.get(lon_kwd)
        #### added B. Williams Apr 24, 2017:
        # need to put in the radar/image coordinates too
        try:
            self.img_x = self.get('range_index')
            self.img_y = self.get('azimuth_index')
        except:
            print('Cant Find range_index, or azimuth_index variables, assuming 2D-image image coordinates (like from a gdem)')
            Ny,Nx = N.shape(self.nc.variables[lat_kwd])# this should break if it is not an image like a gdem or an L2 pixel cloud
            ix,iy = N.meshgrid(N.arange(Nx),N.arange(Ny))
            self.img_x = ix[self.index]
            self.img_y = iy[self.index]
        ####

        if self.verbose: print('lat/lon read')

        # If not enough good points are found, raise Exception

        if len(self.lat) < min_points:
            raise Exception('number of good points: %d smaller than required: %d'%(
                len(self.lat),min_points) )

        # Project to a coordinate system

        if project_data:
            self.x, self.y = self.project(proj=proj,x_0=x_0,y_0=y_0,lat_0=lat_0,lon_0=lon_0,
                                        ellps=ellps,**proj_kwds)
            if self.verbose: print('projection set and x,y calculated')

    def set_index_and_bounding_box(self,bounding_box,lat_kwd,lon_kwd,
                                   class_list,class_kwd='no_layover_classification'):
        """Set the index for the good data and computes the bounding box for the good data.

        class_kwd: the netcdf name of the classification layer to use.
        """

         # Read the classification and update the index

        self.klass = self.nc.variables[class_kwd][:]

        self.class_list = class_list

        self.class_index = ( self.klass == class_list[0] )
        for i in range(1,len(class_list)):
            self.class_index = self.class_index | (self.klass == class_list[i])

        if self.verbose: print 'Number of points in these classes: %d'%(N.sum(self.class_index))
        
        lat = self.nc.variables[lat_kwd][:]
        lon = self.nc.variables[lon_kwd][:]
            
        if bounding_box != None:
            self.lonmin,self.latmin,self.lonmax,self.latmax = bounding_box
        else:
            self.lonmin = lon.min()
            self.latmin = lat.min()
            self.lonmax = lon.max()
            self.latmax = lat.max()

        self.index = ( (lat >= self.latmin) &
                       (lon >= self.lonmin) &
                       (lat <= self.latmax) &
                       (lon <= self.lonmax) )
        if self.verbose: print 'Number of points in bounding box: %d'%(N.sum(self.index))

        self.index = self.index & self.class_index
        if self.verbose: print 'Number of good: %d'%(N.sum(self.index))
        
        lat = lat[self.index]
        lon = lon[self.index]

        # Brent Williams: commented out next lines to handle setting of bounding box not by data
        # so that you can run with gdem and L2 file and get same nodes
        #self.lonmin = lon.min()
        #self.latmin = lat.min()
        #self.lonmax = lon.max()
        #self.latmax = lat.max()
        #print(lon.min(),lon.max(),lat.min(),lat.max())
        self.bounding_box = (self.lonmin,self.latmin,self.lonmax,self.latmax)

        # Update the classification

        self.klass = self.klass[self.index]

    def get(self,var):
        """Get the values of the variable var within the desired index of good sites."""

        x = self.nc.variables[var][:]
        return x[self.index]

    def project(self,proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,ellps='WGS84',**proj_kwds):
        """Get x and y coordinates for the selected points using a proj4 projection.

        The proj4 template should replace lat_0, lon_0, x_0, y_0, proj, which
        are passed as keyword parameters (default as below).

        A full list of projection options to set plus explanations of their
        meaning can be found here: https://trac.osgeo.org/proj/wiki/GenParms
        
        The default projection is Lambert equiarea, which has a proj4 string with the
        following parameters:
        
        +proj=laea  
        +lat_0=Latitude at projection center, set to bounding box center lat 
        +lon_0=Longitude at projection center, set to bounding box center lon
        +x_0=False Easting, set to 0
        +y_0=False Northing, set to 0

        """

        # Find lat_0 and lon_0 if not specified previously
        # modified by Brent Williams, May 2017 to use center of bounding box instead of data centroid
        # in order to get same nodes for gdem and l2 file
        if lat_0 == None:
            #lat_0 = N.mean(self.lat)
            # use center of bouding box instead of data centroid 
            lat_0=(self.bounding_box[3]+self.bounding_box[1])/2.0
        if lon_0 == None:
            #lon_0 = N.mean(self.lon)
            # use center of bouding box instead of data centroid
            lon_0=(self.bounding_box[2]+self.bounding_box[0])/2.0
        self.proj = Proj(proj=proj,lat_0=lat_0,lon_0=lon_0,x_0=x_0,y_0=y_0,ellps=ellps,**proj_kwds)

        self.x, self.y = self.proj(self.lon,self.lat)

        return self.x, self.y
