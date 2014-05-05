"""
Given a SWOTL2 file, fit all of the reaches observed and output results.
"""

from os.path import splitext, split, exists
from collections import OrderedDict as odict
import numpy as N
import pandas as pd
from SWOTL2 import SWOTL2
from SWOTRiver import ReachExtractor
from SWOTRiver import IteratedRiverObs
from FitRiver import FitRiver

class SWOTRiverEstimator(SWOTL2):
    """Given a SWOTL2 file, fit all of the reaches observed and output results.

    This class is derived from the SWOTL2 class.
    """

    xtrack_res_max = 100. # Maximum allowed cross-track range resolution 

    def __init__(self,swotL2_file,bounding_box=None,class_list=[1],
                 lat_kwd='no_layover_latitude', lon_kwd='no_layover_longitude',
                 class_kwd='no_layover_classification',
                 min_points=100,
                 height_kwd='height',true_height_kwd='water_height',
                 no_noise_height_kwd='no_noise_height',
                 xtrack_kwd='no_layover_cross_track',
                 xtrack_res_kwd = 'no_layover_ground_resolution',
                proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,
                ellps='WGS84',**proj_kwds):
        """Initialize as with a SWOTL2 file, but with some added parameters.

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
        good_class : list
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

        self.input_file = split(swotL2_file)[-1]
        self.output_file = None

        # Initialize the base class

        SWOTL2.__init__(self,swotL2_file,bounding_box=bounding_box,
                        class_list=class_list,
                        lat_kwd=lat_kwd, lon_kwd=lon_kwd,
                        class_kwd=class_kwd,
                        min_points=min_points,
                        proj=proj,x_0=x_0,y_0=y_0,lat_0=lat_0,lon_0=lon_0,
                        ellps=ellps,**proj_kwds)

        

        self.h_noise = self.get(height_kwd)
        self.h_no_noise = self.get(no_noise_height_kwd)
        self.h_true = self.get(true_height_kwd)
        self.xtrack = self.get(xtrack_kwd)
        self.xtrack_res = self.get(xtrack_res_kwd)

        # This is necessary as some pixel values seem flaky (TO BE FIXED)
        
        self.xtrack_res = N.where(self.xtrack_res > self.xtrack_res_max,
                                  self.xtrack_res_max, self.xtrack_res)

        self.pixel_area = self.azimuth_spacing*self.xtrack_res
        
        print('Data loaded')

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

    def process_reaches(self,use_width_db=True,refine_centerline=True,
                        smooth=1.e-2,alpha=1.,max_iter=1,
                        scalar_max_width=600.,
                        subreach_size=10.e3,
                        min_fit_points=150,
                        step=None,fit_types=['OLS','WLS','RLM'],
                        ds=None,smin=0,minobs=10,max_width=None):
        """Process all of the reaches in the data bounding box.

        Parameters
        ----------

        use_width_db : bool
            Use the width data base for setting widths?
        refine_centerline: bool
            Refine the centerline?
        smooth : float
            Centerline smoothing constant (see Centerline)
        alpha : float
            Centerline refinement update weight
        max_iter : int
            Maximum number of centerline iterations
        scalar_max_width : float
            How far away to look for points
        subreach_size : float
            Size of each subreach (in m)
        min_fit_points : int
            Minimum number of points required for fit                
        step : float
            Change smin (in m) with each new estimate for this reach
            (if None -> subreach_size)
        fit_types : list
            A list of fit types to perform. Can contain 'OLS','WLS','RLM'.
        ds : float
            Separation between centerline nodes (in m)
        smin : float
            Start reach distance (in m).
        minobs : int
            Minimum number of observations for valid node.
        max_width: float or array_like
            Maximum width to use for accepting points.                           
        """

        for i, reach_idx in enumerate(self.reaches.reach_idx):
            print('Reach %d/%d Reach index: %d'%(i+1,self.reaches.nreaches,reach_idx))
            
            if use_width_db:
                max_width = self.get_max_width_from_db(reach_idx)
                print('max_width read')

            results_df = self.process_reach(self.reaches[i],subreach_size,
                            refine_centerline=refine_centerline,
                            smooth=smooth,alpha=alpha,max_iter=max_iter,
                            scalar_max_width=scalar_max_width,
                            min_fit_points=min_fit_points,
                            step=step,fit_types=fit_types,
                            ds=ds,smin=smin,minobs=minobs,max_width=max_width)
            print('reach pocessed')

            if self.output_file != None:
                try:
                    csv_file = splitext(self.input_file)[0]+'_fit_%d.csv'%reach_idx
                    results_df.to_csv(csv_file)
                    self.store.append('fit_results',results_df)
                except:
                    print 'Could not write results'

    def process_reach(self,reach,subreach_size,
                      refine_centerline=True,
                        smooth=1.e-2,alpha=1.,max_iter=1,
                        scalar_max_width=600.,
                        min_fit_points=150,
                        step=None,fit_types=['OLS','WLS','RLM'],
                        ds=None,smin=0,minobs=10,max_width=None):
        """Estimate the result for one reach.

        Parameters
        ----------

        reach : Reach instance
            One of the reaches from ReachExtractor.
        use_width_db : bool
            Use the width data base for setting widths?
        refine_centerline: bool
            Refine the centerline?
        smooth : float
            Centerline smoothing constant (see Centerline)
        alpha : float
            Centerline refinement update weight
        max_iter : int
            Maximum number of centerline iterations
        scalar_max_width : float
            How far away to look for points                
        subreach_size : float
            Size of each subreach (in m)
        min_fit_points : int
            Minimum number of points required for fit
        step : float
            Change smin (in m) with each new estimate for this reach
            (if None -> subreach_size)
        fit_types : list
            A list of fit types to perform. Can contain 'OLS','WLS','RLM'.
        ds : float
            Separation between centerline nodes (in m)
        smin : float
            Start reach distance (in m).
        minobs : int
            Minimum number of observations for valid node.
        max_width: float or array_like
            Maximum width to use for accepting points.
        """

        if max_width == None:
            max_width = reach.metadata['width_max']

        # Initialize the observations

        if refine_centerline:
            self.river_obs = IteratedRiverObs(reach,self.x,self.y,ds=ds,
                                              max_width=scalar_max_width,minobs=minobs)
        else:
            self.river_obs = IteratedRiverObs(reach,self.x,self.y,ds=ds,
                                              max_width=max_width,minobs=minobs)
        print('river_obs initilized')

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
            print('centerline refined')

        # Add the abservations
        
        self.river_obs.add_obs('h_true',self.h_true)
        self.river_obs.add_obs('h_noise',self.h_noise)
        self.river_obs.add_obs('h_no_noise',self.h_no_noise)
        self.river_obs.add_obs('lon',self.lon)
        self.river_obs.add_obs('lat',self.lat)
        self.river_obs.add_obs('xtrack',self.xtrack)
        self.river_obs.add_obs('pixel_area',self.pixel_area)
        self.river_obs.load_nodes(['h_noise','h_true','h_no_noise',
                                   'lon','lat','xtrack','pixel_area'])
        print('Observations added to nodes')

        # Get various  node statistics
        
        s = N.asarray( self.river_obs.get_node_stat('mean','s') )
        xtrack_mean = N.asarray( self.river_obs.get_node_stat('mean','xtrack') )
        xtrack_min = N.asarray( self.river_obs.get_node_stat('min','xtrack') )
        xtrack_max = N.asarray( self.river_obs.get_node_stat('max','xtrack') )
        lon_min = N.asarray( self.river_obs.get_node_stat('min','lon') )
        lon_max = N.asarray( self.river_obs.get_node_stat('max','lon') )
        lat_min = N.asarray( self.river_obs.get_node_stat('min','lat') )
        lat_max = N.asarray( self.river_obs.get_node_stat('max','lat') )

        # The following are estimates of river width

        width_ptp =  N.asarray( self.river_obs.get_node_stat('width_ptp','') )
        width_std =  N.asarray( self.river_obs.get_node_stat('width_std','') )

        # These are area based estimates, the nominal SWOT approach
        
        area = N.asarray( self.river_obs.get_node_stat('sum','pixel_area') )
        width_area = N.asarray( self.river_obs.get_node_stat('width_area','pixel_area') )

        # Initialize the fitter
        
        self.fitter = FitRiver(self.river_obs)

        sstop = self.river_obs.centerline.s.max()
        smax = smin + subreach_size

        first = True
        results_df = None
        while smax <= sstop:

            print('Processing subreach smin: %.0f smax: %.0f'%(smin/1.e3,smax/1.e3))

            # Check to see if there are sufficient number of points for fit
            
            good = ( s >= smin ) & ( s <= smax )
            ngood = N.sum(good)
            print('number of fit points: %d'%ngood)
            
            if ngood < min_fit_points:
                print('not enough good points going to next reach')
                smin += step
                smax += step
                continue
                
            # Do the fitting for this subreach
          
            nresults, nnresults, tresults = self.process_subreach(smin,smax,
                                                                  fit_types=fit_types)
            print('Estimation finished')

            # Get the reach statistics for this subreach

            reach_stats = self.get_reach_stats(good,
                                               lon_min, lat_min, lon_max, lat_max,
                                               xtrack_mean,xtrack_min,xtrack_max,
                                               width_ptp,width_std,width_area,area)
            print('stats calculated')
            print('width_ptp: %.1f'%reach_stats.width_ptp_mean)
            print('width_std: %.1f'%reach_stats.width_std_mean)
            print('width_area: %.1f'%reach_stats.width_area_mean)

            # Put the results into a data frame

            if first:
                results_df = self.get_results_df(smin,smax,reach,
                                            nresults, nnresults, tresults,
                                            reach_stats)
                first = False
            else:
                df = self.get_results_df(smin,smax,reach,
                                        nresults, nnresults, tresults,
                                        reach_stats)
                print df.shape
                #results_df.append(df)#,ignore_index=True)
                results_df = pd.concat([results_df, df])
                print results_df.shape
                
            print('Data frame updated')
                

            smin += step
            smax += step

        return results_df

    def process_subreach(self,smin,smax,fit_types=['OLS','WLS','RLM']):
        """Get fit results for this subreach."""
        
        if type(fit_types) == str:
            fit_types = [fit_types]

        tresults = odict()
        load_inputs=True
        for fit_type in ['OLS']: #fit_types:
            tresults[fit_type] = self.fitter.fit_linear(smin,smax,'h_true',
                                                       fit=fit_type,
                                                       load_inputs=load_inputs)
            load_inputs = False

        nresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            nresults[fit_type] = self.fitter.fit_linear(smin,smax,'h_noise',
                                                       fit=fit_type,
                                                       load_inputs=load_inputs)
            load_inputs = False

        nnresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            nnresults[fit_type] = self.fitter.fit_linear(smin,smax,'h_no_noise',
                                                       fit=fit_type,
                                                       load_inputs=load_inputs)
            load_inputs = False

        return nresults, nnresults, tresults

    def get_reach_stats(self,good,
                        lon_min, lat_min, lon_max, lat_max,
                        xtrack_mean,xtrack_min,xtrack_max,
                        width_ptp,width_std,width_area,area):
        """Get statistics for a given reach."""

        class reach_stats:
            pass

        reach_stats.lon_min = N.min(lon_min[good])
        reach_stats.lon_max = N.max(lon_max[good])
        reach_stats.lat_min = N.min(lat_min[good])
        reach_stats.lat_max = N.max(lat_max[good])
        reach_stats.xtrack_mean = N.mean(xtrack_mean[good])
        reach_stats.xtrack_min = N.min(xtrack_mean[good])
        reach_stats.xtrack_max = N.max(xtrack_mean[good])
        reach_stats.width_ptp_mean = N.mean(width_ptp[good])
        reach_stats.width_ptp_min = N.min(width_ptp[good])
        reach_stats.width_ptp_max = N.max(width_ptp[good])
        reach_stats.width_std_mean = N.mean(width_std[good])
        reach_stats.width_std_min = N.min(width_std[good])
        reach_stats.width_std_max = N.max(width_std[good])
        reach_stats.width_area_mean = N.mean(width_area[good])
        reach_stats.width_area_min = N.min(width_area[good])
        reach_stats.width_area_max = N.max(width_area[good])
        reach_stats.area = N.sum(area[good])

        return reach_stats

    def get_results_df(self,
                       smin,smax,reach,nresults, nnresults, tresults,
                       reach_stats):
        """Stuff the estimation results into a pandas DataFrame."""

        reach_idx = reach.metadata['reach_idx']
        
        results_df = pd.DataFrame(odict([

            # Location results
            
            ('reach_idx' , [reach_idx]),
            ('lonmin' , [reach_stats.lon_min]),
            ('latmin' , [reach_stats.lat_min]),
            ('lonmax' , [reach_stats.lon_max]),
            ('latmax' , [reach_stats.lat_max]),
            ('input_file' , [self.input_file]),
            ('smin' , [smin]),
            ('smax' , [smax]),
            ('xtrack_mean' , [reach_stats.xtrack_mean]),
            ('xtrack_min' , [reach_stats.xtrack_min]),
            ('xtrack_max' , [reach_stats.xtrack_max]),
            
            # width results
            
            ('width_ptp_mean' , [reach_stats.width_ptp_mean]),
            ('width_ptp_min' , [reach_stats.width_ptp_min]),
            ('width_ptp_max' , [reach_stats.width_ptp_max]),
            ('width_std_mean' , [reach_stats.width_std_mean]),
            ('width_std_min' , [reach_stats.width_std_min]),
            ('width_std_max' , [reach_stats.width_std_max]),
            ('width_area_mean' ,[reach_stats.width_area_mean]),
            ('width_area_min' , [reach_stats.width_area_min]),
            ('width_area_max' , [reach_stats.width_area_max]),
            
            # fitting results for h_true
            
            ('ols_h_true_mean' , [tresults['OLS'].params[1]]),
            ('ols_h_true_slope' , [tresults['OLS'].params[0]]),
            ('ols_h_true_rsquared' , [tresults['OLS'].rsquared]),
            ('ols_h_true_mse_resid' , [tresults['OLS'].mse_resid]),
            
            # fitting results for h_noise
            
            ('ols_h_noise_mean' , [nresults['OLS'].params[1]]),
            ('ols_h_noise_slope' , [nresults['OLS'].params[0]]),
            ('ols_h_noise_rsquared' , [nresults['OLS'].rsquared]),
            ('ols_h_noise_mse_resid' , [nresults['OLS'].mse_resid]),
            
            ('wls_h_noise_mean' , [nresults['WLS'].params[1]]),
            ('wls_h_noise_slope' , [nresults['WLS'].params[0]]),
            ('wls_h_noise_rsquared' , [nresults['WLS'].rsquared]),
            ('wls_h_noise_mse_resid' , [nresults['WLS'].mse_resid]),
            
            ('rlm_h_noise_mean' , [nresults['RLM'].params[1]]),
            ('rlm_h_noise_slope' , [nresults['RLM'].params[0]]),
            #('rlm_h_noise_rsquared' , [nresults['RLM'].rsquared]),
            #('rlm_h_noise_mse_resid' , [nresults['RLM'].mse_resid]),                        
            
            # fitting results for h_no_noise
            
            ('ols_h_no_noise_mean' , [nnresults['OLS'].params[1]]),
            ('ols_h_no_noise_slope' , [nnresults['OLS'].params[0]]),
            ('ols_h_no_noise_rsquared' , [nnresults['OLS'].rsquared]),
            ('ols_h_no_noise_mse_resid' , [nnresults['OLS'].mse_resid]),
            
            ('wls_h_no_noise_mean' , [nnresults['WLS'].params[1]]),
            ('wls_h_no_noise_slope' , [nnresults['WLS'].params[0]]),
            ('wls_h_no_noise_rsquared' , [nnresults['WLS'].rsquared]),
            ('wls_h_no_noise_mse_resid' , [nnresults['WLS'].mse_resid]),
            
            ('rlm_h_no_noise_mean' , [nnresults['RLM'].params[1]]),
            ('rlm_h_no_noise_slope' , [nnresults['RLM'].params[0]]),
            #('rlm_h_no_noise_rsquared' , [nnresults['RLM'].rsquared]),
            #('rlm_h_no_noise_mse_resid' , [nnresults['RLM'].mse_resid]),                        
            
            ]))

        return results_df

        
