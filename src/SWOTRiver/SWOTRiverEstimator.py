"""
Given a SWOTL2 file, fit all of the reaches observed and output results.
"""

from os.path import splitext, split, exists
from collections import OrderedDict as odict
import numpy as N
from numpy.ma import is_masked
import pandas as pd
from SWOTL2 import SWOTL2
from RiverObs import ReachExtractor
from RiverObs import WidthDataBase
from RiverObs import IteratedRiverObs
from RiverObs import FitRiver
from RiverObs import RiverNode
from RiverObs import RiverReach

class SWOTRiverEstimator(SWOTL2):
    """Given a SWOTL2 file, fit all of the reaches observed and output results.

    This class is derived from the SWOTL2 class.

    This class contains four basic components:

    1. It is a subclass of the SWOTL2, and therefor is a LatLonRegion and also contains SWOT data.
    2. For each reach of interest, it keeps a dictionary RiverObs instances for each reach.
    3. It keeps a dictionary of RiverReach instances, containing the results of height/slope fitting and width estimation.
    4. A dictionary of FitRiver instances, containing the full fit diagnostics.

    Keeping copies for each observation, although useful for visualization, can be a tax on memory
    for a large number of reaches and observation, so that storing can be turned off, if so desired.

    Parameters
    ----------

    height_kwd : str
        name of netcdf variable to use as a measurement.
    true_height_kwd: str
        name of variable to use as truth for comparison.
    no_noise_height_kwd : str
        name for no noise measurement
    xtrack_kwd : str
        name of cross-track distance variable.
    bounding_box : 4-tuple or array_like
        If the bounding_box is provided, select only the data in the
        bounding box, otherwise, get the bounding box from the data
        itself. bounding_box should be of the form
        (lonmin,latmin,lonmax,latmax) 
    lon_kwd, lat_kwd : str
        netcdf names for the longitudes and latitudes to be used
        for georeferencing the data set. The a priori geolocation can be retrieved
        using 'no_layover_longitude' and 'no_layover_latitude'. The estimated
        latitude and longitude can be retrieved using 'longitude' and 'latitude',
        respectively.
    class_kwd : str, default 'classification'
        the netcdf name of the classification layer to use.
        'classification' will return the result of applying a classifier to the simulated data.
        'no_layover_classification' can be used as the thruth
    class_list : list, default [2,3,4]
        a list of the class labels for what is considered good water data.
        This should be piched as any classes that might contain water, even
        if partially. If truth data are desired (class_kwd='no_layover_classification'),
        the this should be set to [1]. The default has interior water pixels(4),
        border water pixels (3), and border land pixels (2). This should be used with
        inundation_fraction turned on.
    fractional_inundation_kwd : str, default 'continuous_classification'
        Netcdf keyword containing the inundation fraction. If None, the no inundation
        fraction is used. The inundation fraction is, in theory, a number between 0 and
        1 estimating the fraction of the pixel covered with water. In practice, because of noise,
        it can be outside the bounds and even be negative! It will produced an ensemble
        unbiased estimate if the class mean cross sections are known.
    use_fractional_inundation : bool list, default [True, True, False]
        For which classes should the inundation fraction be used instead of a binary classification?
        the default is to assume that interior pixels are 100% water, but to use both land
        and water edge pixels partially by using the fractional inundation kwd.
    min_points : int
        If the number of good points is less than this, raise an exception.
    height_kwd : str, default 'height'
        These are the heights to exract from the water file.
    true_height_kwd : str, default 'water_height'
        These are the reference heights to use for performance assessment
    no_noise_height_kwd : str, default 'no_noise_height'
        These are the a priori heights used for phase unwrapping
    xtrack_kwd : str, 'no_layover_cross_track'
        This is the distance to the nadir track.
    xtrack_res_kwd : str, default 'flat_pixel_resolution'
        This is the cross-track dimension of the radar pixel. The netcdf file
        has an estimate 'no_layover_ground_resolution' that takes into account
        a priori surface slopes to derive an incidence angle. However, this runs
        into trouble if water
    trim_ends : bool, default False
        If True, the first and last nodes are not used, to avoid potential edge effects.
    store_obs : bool, default True
        If True, store each RiverObs instance in a dictionary.
    store_reaches : bool, default True
        If True, store each RiverRiver instance in a dictionary.
    store_fits : bool, default True
        If True, store each fit result in a dictionary.        
    verbose : bool, default False
        If True, print to sdtout while processing.

    The final set of keywords are projection options for pyproj. See Notes.

    Notes
    -----
        
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

    xtrack_res_max = 200. # Maximum allowed cross-track range resolution
    platform_height = 970.e3
    earth_radius = 6378e3

    def __init__(self,swotL2_file,bounding_box=None,
                 lat_kwd='no_layover_latitude', lon_kwd='no_layover_longitude',
                 class_list=[2,3,4],class_kwd='classification',
                 fractional_inundation_kwd='continuous_classification',
                 use_fractional_inundation=[True, True, False],
                 min_points=100,
                 height_kwd='height',true_height_kwd='water_height',
                 no_noise_height_kwd='no_noise_height',
                 xtrack_kwd='no_layover_cross_track',
                 xtrack_res_kwd = 'flat_pixel_resolution',
                 trim_ends=False,
                 store_obs=True,
                store_reaches=True,
                store_fits=True,
                verbose=False,
                proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,
                ellps='WGS84',**proj_kwds):

        
        self.trim_ends = trim_ends
        self.store_obs=store_obs
        self.store_reaches=store_reaches
        self.store_fits=store_fits
        self.verbose=verbose
                
        self.input_file = split(swotL2_file)[-1]
        self.output_file = None

        # Classification inputs

        self.class_kwd = class_kwd
        self.class_list = class_list
        self.fractional_inundation_kwd = fractional_inundation_kwd
        self.use_fractional_inundation = use_fractional_inundation

        # Initialize the base class

        SWOTL2.__init__(self,swotL2_file,bounding_box=bounding_box,
                        class_list=class_list,
                        lat_kwd=lat_kwd, lon_kwd=lon_kwd,
                        class_kwd=class_kwd,
                        min_points=min_points,
                        verbose=verbose,
                        proj=proj,x_0=x_0,y_0=y_0,lat_0=lat_0,lon_0=lon_0,
                        ellps=ellps,**proj_kwds)
        if is_masked(self.lat):
            mask = self.lat.mask
        else:
            mask = N.zeros(len(self.lat),dtype=N.bool)
        
        self.h_noise = self.get(height_kwd)
        if is_masked(self.h_noise):
            mask = mask | self.h_noise.mask
            
        self.h_no_noise = self.get(no_noise_height_kwd)
        if is_masked(self.h_no_noise):
            mask = mask | self.h_no_noise.mask
            
        self.h_true = self.get(true_height_kwd)
        if is_masked(self.h_true):
            mask = mask | self.h_true.mask
            
        self.xtrack = self.get(xtrack_kwd)
        if is_masked(self.xtrack):
            mask = mask | self.xtrack.mask
            
        look_angle = self.get('no_layover_look_angle')
        if is_masked(look_angle):
            mask = mask | look_angle.mask
            
        self.incidence_angle = (look_angle)*(1. + self.platform_height/self.earth_radius)
        self.range_resolution = float(self.nc.range_resolution)

        # Get the pixel area, depending on the choice of range resolution

        if xtrack_res_kwd != 'flat_pixel_resolution':
            self.xtrack_res = self.get(xtrack_res_kwd)
            if is_masked(self.xtrack_res):
                mask = mask | self.xtrack_res.mask

            # This is necessary as some pixel values seem flaky (TO BE FIXED)
        
            self.xtrack_res = N.where(self.xtrack_res > self.xtrack_res_max,
                                    self.xtrack_res_max, self.xtrack_res)
        else:
            self.xtrack_res = self.range_resolution/N.sin(N.radians(self.incidence_angle))

        # Make sure all of the arrays are valid over the same points

        good = ~mask
        self.lat = self.lat[good]
        self.lon = self.lon[good]
        self.x = self.x[good]
        self.y = self.y[good]
        self.klass = self.klass[good]
        
        self.h_noise = self.h_noise[good]
        self.h_no_noise = self.h_no_noise[good]
        self.h_true = self.h_true[good] 
        self.xtrack = self.xtrack[good]
        self.incidence_angle = self.incidence_angle[good]
        self.xtrack_res = self.xtrack_res[good]
            
        self.pixel_area = self.azimuth_spacing*self.xtrack_res

        # Calculate the inundated area for each pixel

        if fractional_inundation_kwd == None: # all water pixels are inundated
            self.fractional_inundation = None
            self.inundated_area = self.pixel_area
        else:
            self.fractional_inundation = self.get(fractional_inundation_kwd)
            self.inundated_area = self.pixel_area
            for i,k in enumerate(class_list):
                if use_fractional_inundation[i]:
                    index = self.klass == k
                    self.inundated_area[index] = self.pixel_area[index]*self.fractional_inundation[index]
        
        if self.verbose: print('Data loaded')

        # Initialize the list of observations and reaches

        self.river_obs_collection = odict()
        self.river_reach_collection = odict()
        self.fit_collection = odict()

    def get_reaches(self,shape_file_root, clip=True,clip_buffer=0.1):
        """Get all of the reaches using a ReachExtractor."""

        self.clip = clip
        self.clip_buffer = clip_buffer
        self.reaches = ReachExtractor(shape_file_root, self, clip=clip, clip_buffer=clip_buffer)

        return self.reaches

    def get_width_db(self,width_db_file):
        """Open the width data base for later use."""

        self.width_db = WidthDataBase(width_db_file)

    def set_width_db(self, width_db):
        """Set width data base from an already opened version."""
        
        self.width_db = width_db
        
    def get_max_width_from_db(self,reach_idx):
        """Get the width associated with a given reach from the width data base."""

        return self.width_db.get_river(reach_idx,
                                       columns=['width'],
                                asarray=True,transpose=False,
                                bounding_box=self.bounding_box,
                                clip_buffer=self.clip_buffer)

    def init_output_file(self,output_file,mode='a'):
        """Initialize the output file HDFStore for later writing."""

        self.output_file = output_file
        self.store = pd.HDFStore(output_file,mode=mode)

    def process_reaches(self,
                scalar_max_width=600.,
                minobs=10,min_fit_points=3,
                fit_types=['OLS','WLS','RLM'],
                use_width_db = False,
                ds=None,
                refine_centerline=False,
                smooth=1.e-2,alpha=1.,max_iter=1,
                ):
        """Process all of the reaches in the data bounding box.

        Parameters
        ----------

        scalar_max_width : float, default 600
            How far away to look for points
        minobs : int, default 10
            Minimum number of observations for valid node.
        min_fit_points : int, default 3
            Minimum number of populated nodes required for height/slope fit
        fit_types : list
            A list of fit types to perform. Can contain 'OLS','WLS','RLM'.
        use_width_db : bool, default False
            Use the width data base for setting widths?
        ds : float, optional
            Separation between centerline nodes (in m). If None, uses reach points.
        refine_centerline: bool, default True
            Refine the centerline?
        smooth : float, default 1.e-2
            Centerline smoothing constant (see Centerline)
        alpha : float, default 1
            Centerline refinement update weight
        max_iter : int, default 1
            Maximum number of centerline iterations

        Returns
        -------

        Returns a list containg a RiverReach instance for each reach in the bounding box.
        """

        river_reach_collection = []
        for i, reach_idx in enumerate(self.reaches.reach_idx):
            if self.verbose: print('Reach %d/%d Reach index: %d'%(i+1,self.reaches.nreaches,reach_idx))
            
            if use_width_db:
                max_width = self.get_max_width_from_db(reach_idx)
                if self.verbose: print('max_width read')
            else:
                max_width = None

            if type(max_width) == type(None):
                try:
                    max_width = reach.metadata['width_max']
                except:
                    max_width = scalar_max_width

            river_reach = self.process_reach(self.reaches[i],i,reach_idx,
                                            scalar_max_width=scalar_max_width,
                                            minobs=minobs,min_fit_points=min_fit_points,
                                            fit_types=fit_types,
                                            use_width_db = use_width_db,
                                            max_width=max_width,
                                            ds=ds,
                                            refine_centerline=refine_centerline,
                                            smooth=smooth,alpha=alpha,max_iter=max_iter)

            river_reach_collection.append(river_reach)

            if self.store_reaches:
                self.river_reach_collection[i] = river_reach

            if self.verbose: print('reach pocessed')
            
        return river_reach_collection

    def process_reach(self,reach,reach_id,reach_idx=None,
                scalar_max_width=600.,
                minobs=10,min_fit_points=3,
                fit_types=['OLS','WLS','RLM'],
                use_width_db = False,
                max_width=None,
                ds=None,
                refine_centerline=False,
                smooth=1.e-2,alpha=1.,max_iter=1
                ):
        """Estimate the width, height, and slope for one reach.

        Parameters
        ----------

        reach : Reach instance
            One of the reaches from ReachExtractor.
        reach_id : int
            Index in the list of reaches extracted for this scene.
        reach_idx, int
            Reach index used as pointer to a reach collection that may be different than
            the one used as input (e.g., a global reach collection).
            If using a width database, the retrieved width will be associated with this reach index.
            If None, it will be set equal to the reach_id. 
        scalar_max_width : float, default 600
            How far away to look for points
        minobs : int, default 10
            Minimum number of observations for valid node.
        min_fit_points : int, default 3
            Minimum number of populated nodes required for height/slope fit
        fit_types : list
            A list of fit types to perform. Can contain 'OLS','WLS','RLM'.
        use_width_db : bool, default False
            Use the width data base for setting widths?
        max_width: float or array_like, optional
            Maximum width to use for accepting points. From width database or apriori.
        ds : float, optional
            Separation between centerline nodes (in m). If None, uses reach points.
        refine_centerline: bool, default True
            Refine the centerline?
        smooth : float, default 1.e-2
            Centerline smoothing constant (see Centerline)
        alpha : float, default 1
            Centerline refinement update weight
        max_iter : int, default 1
            Maximum number of centerline iterations

        Returns
        -------

        Return a RiverReach instance with the reach information.
        """

        # Initialize the observations

        self.river_obs = IteratedRiverObs(reach,
                                          self.x,
                                          self.y,
                                          ds=ds,
                                          max_width=scalar_max_width,
                                          minobs=minobs,
                                          verbose=self.verbose)
        if self.verbose: print('river_obs initilized')

        # Refine the centerline, if desired

        if refine_centerline:
            self.river_obs.iterate(max_iter=max_iter,alpha=alpha,
                                        weights=True,smooth=smooth)

            # Associate the width to the new centerline

            if N.iterable(max_width):
                xw = reach.x
                yw = reach.y
                self.river_obs.add_centerline_obs(xw,yw,max_width,'max_width')

            # Reinitialize to the new centerline and max_width

            self.river_obs.reinitialize()
            if self.verbose: print('centerline refined')
        else:
            # Associate the width to the new centerline

            if N.iterable(max_width):
                xw = reach.x
                yw = reach.y
                self.river_obs.add_centerline_obs(xw,yw,max_width,'max_width')

        # Exclude beginning and end nodes, if desired

        if self.trim_ends:
            first_node = self.river_obs.populated_nodes[0]
            last_node = self.river_obs.populated_nodes[-1]
            self.river_obs.remove_nodes([first_node,last_node])

        # Add the observations

        self.river_obs.add_obs('h_true',self.h_true)
        self.river_obs.add_obs('h_noise',self.h_noise)
        self.river_obs.add_obs('h_no_noise',self.h_no_noise)
        self.river_obs.add_obs('lon',self.lon)
        self.river_obs.add_obs('lat',self.lat)
        self.river_obs.add_obs('xobs',self.x)
        self.river_obs.add_obs('yobs',self.y)
        self.river_obs.add_obs('xtrack',self.xtrack)
        self.river_obs.add_obs('inundated_area',self.inundated_area)
        self.river_obs.load_nodes(['h_noise','h_true','h_no_noise',
                                  'lon','lat','xobs','yobs','xtrack','inundated_area'])
        if self.verbose: print('Observations added to nodes')

        # Get various node statistics

        nobs = N.asarray( self.river_obs.get_node_stat('count','') )
        s_median = N.asarray( self.river_obs.get_node_stat('median','s') )
        x_median = N.asarray( self.river_obs.get_node_stat('median','xobs') )
        y_median = N.asarray( self.river_obs.get_node_stat('median','yobs') )
        xtrack_median = N.asarray( self.river_obs.get_node_stat('median','xtrack') )
        lon_median = N.asarray( self.river_obs.get_node_stat('median','lon') )
        lat_median = N.asarray( self.river_obs.get_node_stat('median','lat') )
        ds = N.asarray( self.river_obs.get_node_stat('value','ds') )

        h_true_ave = N.asarray( self.river_obs.get_node_stat('median','h_true') )
        h_true_std = N.asarray( self.river_obs.get_node_stat('std','h_true') )
        h_noise_ave = N.asarray( self.river_obs.get_node_stat('median','h_noise') )
        h_noise_std = N.asarray( self.river_obs.get_node_stat('std','h_noise') )
        h_no_noise_ave = N.asarray( self.river_obs.get_node_stat('median','h_no_noise') )
        h_no_noise_std = N.asarray( self.river_obs.get_node_stat('std','h_no_noise') )

        # The following are estimates of river width

        width_ptp =  N.asarray( self.river_obs.get_node_stat('width_ptp','') )
        width_std =  N.asarray( self.river_obs.get_node_stat('width_std','') )

        # These are area based estimates, the nominal SWOT approach
        
        area = N.asarray( self.river_obs.get_node_stat('sum','inundated_area') )
        width_area = N.asarray( self.river_obs.get_node_stat('width_area','inundated_area') )

        # These are the values from the width database

        width_db = N.ones(self.river_obs.n_nodes,dtype=N.float64)*self.river_obs.missing_value
        if use_width_db:
            windex = self.river_obs.centerline_obs['max_width'].populated_nodes
            width_db[windex] = self.river_obs.centerline_obs['max_width'].v
            width_db = width_db[self.river_obs.populated_nodes]

        # Get the start and end

        smin = s_median.min()
        smax = s_median.max()

        # Initialize the fitter
        
        self.fitter = FitRiver(self.river_obs)

        # Check to see if there are sufficient number of points for fit
            

        ngood = len(s_median)
        if self.verbose: print('number of fit points: %d'%ngood)
            
        if ngood < min_fit_points:
            nresults = None
            nnresults = None
            tresults = None
            if self.verbose: print('not enough good points for fit')
        else:
                
            # Do the fitting for this reach
          
            nresults, nnresults, tresults = self.estimate_height_slope(smin,smax,
                                                                        fit_types=fit_types,
                                                                        mean_stat='median')
            if self.verbose: print('Estimation finished')

        if self.store_fits:
            self.fit_collection[reach_id,'true'] = tresults
            self.fit_collection[reach_id,'noise'] = nresults
            self.fit_collection[reach_id,'no_noise'] = nnresults

        # Get the reach statistics for this subreach
        
        reach_stats = self.get_reach_stats(reach_id,reach_idx,smin,smax,
                                           lon_median, lat_median, 
                                            xtrack_median,
                                            width_ptp,width_std,width_area,area,
                                            tresults,nresults,nnresults)

        # Put the results into a RiverReach object

        river_reach = RiverReach(lat=lat_median,
                                 lon=lon_median,
                                 x = x_median,
                                 y = y_median,
                                 nobs = nobs,
                                 s = s_median,
                                 ds = ds,                                 
                                 xtrack = xtrack_median,
                                 w_ptp = width_ptp,
                                 w_std = width_std,
                                 w_area = width_area,
                                 w_db = width_db,
                                 area = area,
                                 h_t_ave = h_true_ave,
                                 h_t_std = h_true_std,                                 
                                 h_n_ave = h_noise_ave,
                                 h_n_std = h_noise_std,
                                 h_nn_ave = h_no_noise_ave,
                                 h_nn_std = h_no_noise_std,
                                 metadata = reach_stats)

        # Store, if desired
        
        if self.store_obs:
            self.river_obs_collection[reach_idx] = self.river_obs
                        
        return river_reach

    def estimate_height_slope(self,smin,smax,fit_types=['OLS','WLS','RLM'],mean_stat='median'):
        """Get fit results for this subreach."""
        
        if type(fit_types) == str:
            fit_types = [fit_types]

        tresults = odict()
        load_inputs=True
        tresults['OLS'] = self.fitter.fit_linear(smin,smax,'h_true',
                                                    fit='OLS',
                                                    mean_stat=mean_stat,
                                                    load_inputs=load_inputs)
        load_inputs = False

        nresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            nresults[fit_type] = self.fitter.fit_linear(smin,smax,'h_noise',
                                                       fit=fit_type,
                                                       mean_stat=mean_stat,
                                                       load_inputs=load_inputs)
            load_inputs = False

        nnresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            nnresults[fit_type] = self.fitter.fit_linear(smin,smax,'h_no_noise',
                                                       fit=fit_type,
                                                       mean_stat=mean_stat,
                                                       load_inputs=load_inputs)
            load_inputs = False

        return nresults, nnresults, tresults

    def get_reach_stats(self,reach_id,reach_idx,smin,smax,
                        lon_median, lat_median, 
                        xtrack_median,
                        width_ptp,width_std,width_area,area,
                        tresults,nresults,nnresults):
        """Get statistics for a given reach."""

        reach_stats = odict()

        reach_stats['reach_id'] = reach_id
        reach_stats['reach_idx'] = reach_idx
        reach_stats['lon_min'] = N.min(lon_median)
        reach_stats['lon_max'] = N.max(lon_median)
        reach_stats['lat_min'] = N.min(lat_median)
        reach_stats['lat_max'] = N.max(lat_median)
        reach_stats['area'] = N.sum(area)
        reach_stats['length'] = smax - smin
        reach_stats['smin'] = smin
        reach_stats['smax'] = smax
        reach_stats['w_ptp_ave'] = N.median(width_ptp)
        reach_stats['w_ptp_min'] = N.min(width_ptp)
        reach_stats['w_ptp_max'] = N.max(width_ptp)
        reach_stats['w_std_ave'] = N.median(width_std)
        reach_stats['w_std_min'] = N.min(width_std)
        reach_stats['w_std_max'] = N.max(width_std)
        reach_stats['w_area_ave'] = N.median(width_area)
        reach_stats['w_area_min'] = N.min(width_area)
        reach_stats['w_area_max'] = N.max(width_area)
        reach_stats['xtrck_ave'] = N.median(xtrack_median)
        reach_stats['xtrck_min'] = N.min(xtrack_median)
        reach_stats['xtrck_max'] = N.max(xtrack_median)

        # fitting results for h_true

        if type(tresults) != type(None):
            reach_stats['h_t'] = tresults['OLS'].params[1]
            reach_stats['slp_t'] = tresults['OLS'].params[0]
            reach_stats['t_rsqrd'] = tresults['OLS'].rsquared
            reach_stats['t_mse'] = tresults['OLS'].mse_resid
        else:
            reach_stats['h_t'] = self.river_obs.missing_value
            reach_stats['slp_t'] = self.river_obs.missing_value
            reach_stats['t_rsqrd'] = self.river_obs.missing_value
            reach_stats['t_mse'] = self.river_obs.missing_value
            
        tag = {'OLS':'o','WLS':'w','RLM':'r'}

        if type(nresults) != type(None):
            for fit_type in ['OLS','WLS','RLM']:
                if fit_type in nresults:
                    reach_stats['h_n%s'%tag[fit_type] ] = nresults[fit_type].params[1]
                    reach_stats['slp_n%s'%tag[fit_type] ] = nresults[fit_type].params[0]
                    try:
                        reach_stats['n%s_rsqrd'%tag[fit_type] ] = nresults[fit_type].rsquared
                        reach_stats['n%s_mse'%tag[fit_type] ] = nresults[fit_type].mse_resid
                    except:
                        reach_stats['n%s_rsqrd'%tag[fit_type] ] = self.river_obs.missing_value
                        reach_stats['n%s_mse'%tag[fit_type] ] = self.river_obs.missing_value
                else:
                    reach_stats['h_n%s'%tag[fit_type] ] = self.river_obs.missing_value
                    reach_stats['slp_n%s'%tag[fit_type] ] = self.river_obs.missing_value
                    reach_stats['n%s_rsqrd'%tag[fit_type] ] = self.river_obs.missing_value
                    reach_stats['n%s_mse'%tag[fit_type] ] = self.river_obs.missing_value

        if type(nnresults) != type(None):                    
            for fit_type in ['OLS','WLS','RLM']:
                if fit_type in nresults:
                    reach_stats['h_nn%s'%tag[fit_type] ] = nnresults[fit_type].params[1]
                    reach_stats['slp_nn%s'%tag[fit_type] ] = nnresults[fit_type].params[0]
                    try:
                        reach_stats['nn%s_rsqrd'%tag[fit_type] ] = nnresults[fit_type].rsquared
                        reach_stats['nn%s_mse'%tag[fit_type] ] = nnresults[fit_type].mse_resid
                    except:
                        reach_stats['nn%s_rsqrd'%tag[fit_type] ] = self.river_obs.missing_value
                        reach_stats['nn%s_mse'%tag[fit_type] ] = self.river_obs.missing_value
                else:
                    reach_stats['h_nn%s'%tag[fit_type] ] = self.river_obs.missing_value
                    reach_stats['slp_nn%s'%tag[fit_type] ] = self.river_obs.missing_value
                    reach_stats['nn%s_rsqrd'%tag[fit_type] ] = self.river_obs.missing_value
                    reach_stats['nn%s_mse'%tag[fit_type] ] = self.river_obs.missing_value
                            
        return reach_stats


        
