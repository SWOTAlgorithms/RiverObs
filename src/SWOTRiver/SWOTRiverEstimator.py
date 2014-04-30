"""
Given a SWOTL2 file, fit all of the reaches observed and output results.
"""

from collections import OrderedDict as odict
from SWOTL2 import SWOTL2
from SWOTRiver import ReachExtractor
from SWOTRiver import RiverObs
from FitRiver import FitRiver

class SWOTRiverEstimator(SWOTL2):
    """Given a SWOTL2 file, fit all of the reaches observed and output results.

    This class is derived from the SWOTL2 class.
    """

    def __init__(self,swotL2_file,bounding_box=None,class_list=[1],
                 lat_kwd='no_layover_latitude', lon_kwd='no_layover_longitude',
                 class_kwd='no_layover_classification',
                 height_kwd='height',true_height_kwd='water_height',
                proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,
                ellps='WGS84',**proj_kwds):
        """Initialize as with a SWOTL2 file, but with some added parameters.

        The added keyword inputs are:

        height_kwd: name of netcdf variable to use as a measurement.
        true_height_kwd: name of variable to use as truth for comparison.
                 
        If the bounding_box is provided, select only the data in the
        bounding box, otherwise, get the bounding box from the data
        itself.

        bounding_box should be of the form (lonmin,latmin,lonmax,latmax) 
        good_class: a list of class labels where the data are returned.

        class_list: a list of the class labels for what is considered good data

        lon_kwd, lat_kwd: netcdf names for the longitudes and latitudes to be used
        for georeferencing the data set.

        class_kwd: the netcdf name of the classification layer to use.

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

        # Initialize the base class

        SWOTL2.__init__(self,swotL2_file,bounding_box=bounding_box,
                        class_list=class_list,
                        lat_kwd=lat_kwd, lon_kwd=lon_kwd,
                        class_kwd=class_kwd,
                        proj=proj,x_0=x_0,y_0=y_0,lat_0=lat_0,lon_0=lon_0,
                        ellps=ellps,**proj_kwds)

        self.h = l2.get(height_kwd)
        self.htrue = l2.get(true_height_kwd)

    def get_reaches(self,shape_file_root, clip=True,clip_buffer=0.1):
        """Get all of the reaches using a ReachExtractor."""

        self.reaches = ReachExtractor(shape_file_root, self, clip=clip, clip_buffer=clip_buffer)

        return self.reaches

    def process_reach(self,reach,subreach_size,step=None,fit_types=['OLS','WLS','RLM'],
                      ds=None,smin=0,minobs=10,max_width=None):
        """Estimate the result for one reach.

        reach: one of the reaches from ReachExtractor.
        subreach_size: size of each subreach
        step: change smin with each new estimate for this reach (if None -> subreach_size)
        fit_types: a list of fit types to perform. Can contain 'OLS','WLS','RLM'.
        ds: separation between centerline nodes
        smin: start reach distance.
        minobs: minimum number of observations for valid node.
        max_width: maximum width to use for accepting points.
        """

        if max_width == None:
            max_width = reach.metadata['width_max']

        # Add the abservations
        
        self.river_obs = RiverObs(reach,self.x,self.y,ds=ds,max_width=max_width,minobs=minobs)
        self.river_obs.add_obs('htrue',self.htrue)
        self.river_obs.add_obs('h',self.h)
        self.river_obs.load_nodes(['h','htrue'])

        # Initialize the fitter
        
        self.fitter = FitRiver(self.river_obs)

        sstop = self.obs.centerline.s.max()
        smax = smin + subreach_size

        while smax <= sstop:

            # Do the fitting for this subreach
            
            nresults, tresults = self.process_subreach(smin,smax,fit_types=fit_types)

            # write the results (MISSING)

            smin += step
            smax += step

    def process_subreach(self,smin,smax,fit_types=['OLS','WLS','RLM']):
        """Get fit results for this subreach."""
        
        if type(fit_types) == str:
            fit_types = [fit_types]

        tresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            tresult[fit_type] = self.fitter.fit_linear(smin,smax,'htrue',
                                                       fit=fit_type,
                                                       load_inputs=load_inputs)
            load_inputs = False

        nresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            nresult[fit_type] = self.fitter.fit_linear(smin,smax,'h',
                                                       fit=fit_type,
                                                       load_inputs=load_inputs)
            load_inputs = False

        return nresults, tresults

         

        
