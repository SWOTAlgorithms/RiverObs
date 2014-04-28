"""
Access SWOT L2 data conveniently.

The L2 Data format is not fixed. The current version has the 
following header information (an example).

dimensions:
	water_record = 27378 ;
variables:
	int range_index(water_record) ;
	int azimuth_index(water_record) ;

    double no_layover_x(water_record) ;
	double no_layover_y(water_record) ;
	double no_layover_z(water_record) ;
    double no_layover_latitude(water_record) ;
	double no_layover_longitude(water_record) ;
	double no_layover_height(water_record) ;

    float water_height(water_record) ;

	float no_layover_look_angle(water_record) ;

    double look_unit_x(water_record) ;
	double look_unit_y(water_record) ;
	double look_unit_z(water_record) ;
	double range_to_reference(water_record) ;
	float altitude(water_record) ;
	double range(water_record) ;

    double no_layover_ground_resolution(water_record) ;
	double no_layover_cross_track(water_record) ;

    double dlook_dphase_x(water_record) ;
	double dlook_dphase_y(water_record) ;
	double dlook_dphase_z(water_record) ;

    float no_noise_ifgram_real(water_record) ;
	float no_noise_ifgram_imag(water_record) ;
	float no_layover_ifgram_real(water_record) ;
	float no_layover_ifgram_imag(water_record) ;
	double no_noise_dphase(water_record) ;
	float no_noise_delta_x(water_record) ;
	float no_noise_delta_y(water_record) ;
	float no_noise_delta_z(water_record) ;
	byte no_layover_classification(water_record) ;

    byte classification(water_record) ;
	double no_noise_latitude(water_record) ;
	double no_noise_longitude(water_record) ;
	double no_noise_height(water_record) ;
	double no_noise_cross_track(water_record) ;

"""

import numpy as N
from netCDF4 import Dataset
from pyproj import Proj

class SWOTL2:
    """Access SWOT L2 data conveniently."""

    def __init__(self,swotL2_file,bounding_box=None,class_list=[1],
                 lat_kwd='no_layover_latitude', lon_kwd='no_layover_longitude', 
                proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,
                ellps='WGS84',**proj_kwds):
        """Open netcdf SWOT L2 file for reading.

        If the bounding_box is provided, select only the data in the
        bounding box, otherwise, get the bounding box from the data
        itself.

        bounding_box should be of the form (lonmin,latmin,lonmax,latmax) 
        good_class: a list of class labels where the data are returned.

        class_list: a list of the class labels for what is considered good data

        lon_kwd, lat_kwd: netcdf names for the longitudes and latitudes to be used
        for georeferencing the data set.

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

        self.lat_kwd, self.lon_kwd = lat_kwd, lon_kwd
        
        self.nc = Dataset(swotL2_file)
        print('Dataset opened')

        self.set_bounding_box(bounding_box,lat_kwd,lon_kwd)
        print('Bounding box calculated')
        
        self.set_index(class_list)
        print('Good data selected')

        # Get locations for these data (note that these should be the reference values)

        self.lat = self.get(lat_kwd)
        self.lon = self.get(lon_kwd)
        print('lat/lon read')

        # Project to a coordinate system

        self.x, self.y = self.project(proj=proj,x_0=x_0,y_0=y_0,lat_0=lat_0,lon_0=lon_0,
                                      ellps=ellps,**proj_kwds)
        print('projection set and x,y calculated')

    def set_bounding_box(self,bounding_box,lat_kwd,lon_kwd):
        """Set the bounding box and identify pixes in the bounding box."""

        self.no_layover_lat = self.nc.variables[lat_kwd][:]
        self.no_layover_lon = self.nc.variables[lon_kwd][:]

        if bounding_box != None:
            self.lonmin,self.latmin,self.lonmax,self.latmax = bounding_box
        else:
            self.lonmin = self.no_layover_lon.min()
            self.latmin = self.no_layover_lat.min()
            self.lonmax = self.no_layover_lon.max()
            self.latmax = self.no_layover_lat.max()

        self.bounding_box = (self.lonmin,self.latmin,self.lonmax,self.latmax)

    def set_index(self,class_list):
        """Set the index for the good data."""
        
        # Initialize the index of data in the bounding box

        self.index = ( (self.no_layover_lat >= self.latmin) &
                       (self.no_layover_lon >= self.lonmin) &
                       (self.no_layover_lat <= self.latmax) &
                       (self.no_layover_lon <= self.lonmax) )

        # Read the classification and update the index

        self.klass = self.nc.variables['no_layover_classification'][:]

        self.class_list = class_list

        self.class_index = ( self.klass == class_list[0] )
        for i in range(1,len(class_list)):
            self.class_index = self.class_index | (self.klass == class_list[i])

        self.index = self.index & self.class_index

    def get(self,var):
        """Get the values of the variable var within the desired index of good sites."""

        return self.nc.variables[var][self.index]

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

        if lat_0 == None:
            lat_0 = N.mean(self.lat)

        if lon_0 == None:
            lon_0 = N.mean(self.lon)

        self.proj = Proj(proj=proj,lat_0=lat_0,lon_0=lon_0,x_0=x_0,y_0=y_0,ellps=ellps,**proj_kwds)

        self.x, self.y = self.proj(self.lon,self.lat)

        return self.x, self.y
