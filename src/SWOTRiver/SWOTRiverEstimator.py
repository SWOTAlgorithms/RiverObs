"""
Given a SWOTL2 file, fit all of the reaches observed and output results.
"""

from os.path import splitext, split, exists
from collections import OrderedDict as odict
import numpy as N
from numpy.ma import is_masked
import pandas as pd
import netCDF4 as nc
from SWOTL2 import SWOTL2
from RiverObs import ReachExtractor
from RiverObs import WidthDataBase
from RiverObs import IteratedRiverObs
from RiverObs import FitRiver
from RiverObs import RiverNode
from RiverObs import RiverReach
from scipy.ndimage import label, grey_dilation

# Brent Williams, May 2017: changed several things to be more consistent with current pixel cloud
# took out some of the fields that we only have access to in the simulator, added output the index file etc...

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
    class_list : list, default [2,3,4,5]
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
    use_fractional_inundation : bool list, default [True, True, False, False]
        For which classes should the inundation fraction be used instead of a binary classification?
        the default is to assume that interior pixels are 100% water, but to use both land
        and water edge pixels partially by using the fractional inundation kwd.
    use_segmentation :  bool list, default [False, True, True, True]
        Selects which classes to assume as water for segmentation purposes 
    use_heights : bool list, default [False, False, True, False]
        Selects which classes to use for estimating heights
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
    use_segmentation : bool list, default [False, True, True, True]
        Defines which classes should the assumed as water for segmatation
        algorithm to label disjoint features
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
                 class_list=[2,3,4,5],class_kwd='classification',
                 fractional_inundation_kwd='continuous_classification',
                 use_fractional_inundation=[True, True, False, False],
                 use_segmentation=[False, True, True, True],
                 use_heights=[False, False, True, False],
                 min_points=100,
                 height_kwd='height',#true_height_kwd='water_height',no_noise_height_kwd='no_noise_height',xtrack_kwd='no_layover_cross_track',xtrack_res_kwd = 'flat_pixel_resolution',
                 trim_ends=False,
                 store_obs=True,
                store_reaches=True,
                store_fits=True,
                verbose=True,
                proj='laea',x_0=0,y_0=0,lat_0=None,lon_0=None,
                ellps='WGS84',output_file=None,**proj_kwds):

        
        self.trim_ends = trim_ends
        self.store_obs=store_obs
        self.store_reaches=store_reaches
        self.store_fits=store_fits
        self.verbose=verbose
                
        self.input_file = split(swotL2_file)[-1]
        self.output_file = output_file # index file
        self.createIndexFile()
        # Classification inputs

        self.class_kwd = class_kwd
        self.class_list = class_list
        self.fractional_inundation_kwd = fractional_inundation_kwd
        self.use_fractional_inundation = use_fractional_inundation

        self.use_segmentation = use_segmentation
        self.use_heights = use_heights

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
        
        # Make sure all of the arrays are valid over the same points
        
        good = ~mask
        self.lat = self.lat[good]
        self.lon = self.lon[good]
        self.x = self.x[good]
        self.y = self.y[good]
        self.klass = self.klass[good]
        
        self.h_noise = self.h_noise[good]
        
        self.img_x = self.img_x[good] # range or x index
        self.img_y = self.img_y[good] # azimuth or y index
        
        #try to read the pixel area from the L2 file, or compute it from look angle and azimuth spacing, or from azimuth spacing and ground spacing
        try:
            #read directly from file
            self.pixel_area = self.get('pixel_area')
            
            print("got pixel area from l2 file")
        except:
            try:
                # read the old style pixc
                # compute from inc angle
                look_angle = self.get('no_layover_look_angle')
                #if is_masked(look_angle):
                #    mask = mask | look_angle.mask
                    
                self.incidence_angle = (look_angle)*(1. + self.platform_height/self.earth_radius)
                self.incidence_angle = self.incidence_angle[good]
                self.xtrack_res = float(self.nc.range_resolution)/N.sin(N.radians(self.incidence_angle))
                
                self.pixel_area = float(self.nc.azimuth_spacing)*self.xtrack_res
                print("computed pixel area from look angle and azimuth spacing")
            except:
                try:
                    # compute from gdem-like attributes
                    self.pixel_area = float(self.nc.azimuth_spacing)*float(self.nc.ground_spacing)*N.ones(N.shape(self.h_noise))
                    print("computed pixel area from gdem style attributes")
                except:
                    self.pixel_area = 10.0*N.zeros(len(self.h_noise))
                    print("could not find correct pixel area parameters setting to zero")
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
        # Segment the image into disjoint features
        if len(self.use_segmentation) == len(self.class_list):
            self.isWater=N.zeros(N.shape(self.h_noise))
            for i,k in enumerate(class_list):
                if self.use_segmentation[i]:
                    index = self.klass == k
                    self.isWater[index] = 1
            self.segmentWaterClass()
        else:
            self.seg_label=None
        # set which class pixels to use for heights (should be modified by layover flag later)
        if len(self.use_heights) == len(self.class_list):
            self.h_flg=N.zeros(N.shape(self.h_noise))
            for i,k in enumerate(class_list):
                if self.use_heights[i]:
                    index = self.klass == k
                    self.h_flg[index] = 1
        else:
            self.h_flg=None
        #print "shape h_noise",N.shape(self.h_noise)
        #print "shape h_flg",N.shape(self.h_flg)
        if self.verbose: print('Data loaded')
        
        # Initialize the list of observations and reaches

        self.river_obs_collection = odict()
        self.river_reach_collection = odict()
        self.fit_collection = odict()
    def segmentWaterClass(self):
        """do image segmentation algorithm on the water class to label unconnected features"""
        # added by Brent Williams, May 2017
        # put pixc back into radar/image coordinates
        maxX=N.max(self.img_x)
        #minX=N.min(self.img_x)
        maxY=N.max(self.img_y)
        #minY=N.min(self.img_y)
        cls_img=N.zeros((maxY+1,maxX+1))
        cls_img[self.img_y,self.img_x]=self.isWater
        # Do some regularization with morphological operations? (TODO)
        #
        # segment the water class image
        lbl,nlbl=label(cls_img)
        
        #
        # assign land edge segments (label 0) to nearest water feature (label >0)
        """using grey dilation for this, probably need to re-think 
        how to handle land edge pixels that touch two different 
        segments (this will arbitrarily assign it to the one with 
        the largest label index)"""
        lbl2=grey_dilation(lbl,3)
        lbl_out=lbl.copy()
        lbl_out[lbl==0]=lbl2[lbl==0]
        #
        # write out segmentation images to a nc file for testing
        dumpLbl=False
        if dumpLbl:
            M0,M1=N.shape(lbl)
            f=nc.Dataset("segDump.nc",'w')
            f.createDimension('azimuth',M0)
            f.createDimension('range',M1)
            ol=f.createVariable('orig_label','i4',('azimuth','range'),fill_value=-128)
            fl=f.createVariable('final_label','i4',('azimuth','range'),fill_value=-128)
            ol[:]=lbl[:]
            fl[:]=lbl_out[:]
            f.close()
        #
        # create the segmentation label variable 
        self.seg_label=lbl_out[self.img_y,self.img_x]
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

    #def init_output_file(self,output_file,mode='a'):
    #    """Initialize the output file HDFStore for later writing."""
    #
    #    self.output_file = output_file
    #    self.store = pd.HDFStore(output_file,mode=mode)
    
    def process_reaches(self,
                scalar_max_width=600.,
                minobs=10,min_fit_points=3,
                fit_types=['OLS','WLS','RLM'],
                use_width_db = False,
                ds=None,
                refine_centerline=False,
                smooth=1.e-2,alpha=1.,max_iter=1
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
                    # probably should scale this to look some fraction farther than the database width
                    #max_width = self.reaches[i].metadata['width_max']#
                    max_width = self.reaches[i].metadata['Wmean']*N.ones(N.shape(self.reaches[i].x))#*2.0
                except:
                    max_width = scalar_max_width
            #print "reach idx",reach_idx
            #print "reaches[i].x",self.reaches[i].x
            river_reach = self.process_reach(self.reaches[i],i,reach_idx,
                                            scalar_max_width=scalar_max_width,
                                            minobs=minobs,min_fit_points=min_fit_points,
                                            fit_types=fit_types,
                                            use_width_db = use_width_db,
                                            max_width=max_width,
                                            ds=ds,
                                            refine_centerline=refine_centerline,
                                            smooth=smooth,alpha=alpha,max_iter=max_iter)
            if type(river_reach)!=type(None):
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
        #print self.seg_label
        # Initialize the observations
        try:
            self.river_obs = IteratedRiverObs(reach,
                                              self.x,
                                              self.y,
                                              ds=ds,
                                              seg_label=self.seg_label,#h_flg=self.h_flg,
                                              max_width=scalar_max_width,
                                              minobs=minobs,
                                              verbose=self.verbose)
        except:
            return None

        if self.verbose: print('river_obs initilized')
        
        # Refine the centerline, if desired

        if len(self.river_obs.x)==0:
            print('No observations mapped to nodes in this reach')
            return None
        # get the number of node inthe reach and only refine if there are enough to do spline
        enough_nodes=False
        numNodes=len(N.unique(self.river_obs.index))
        print "numNodes,k:",numNodes,self.river_obs.k
        if (numNodes-1>self.river_obs.k):
            enough_nodes=True
        if (refine_centerline)&(enough_nodes):
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
        
        # write out the image coordinates for each node in a netcdf file
        try:
            segOut=self.seg_label[self.river_obs.in_channel]
        except:
            segOut=None
        self.writeIndexFile(self.img_x[self.river_obs.in_channel],self.img_y[self.river_obs.in_channel],self.river_obs.index,self.river_obs.d,reach_idx,segOut,self.h_flg[self.river_obs.in_channel])
        
        
        # get the prior locations and indices of the nodes
        xw=reach.x
        yw=reach.y
        self.river_obs.add_centerline_obs(xw,yw,xw,'x_prior')
        self.river_obs.add_centerline_obs(xw,yw,yw,'y_prior')
        # should probably compute prior lat/lon from prior x,y too
        windex = self.river_obs.centerline_obs['x_prior'].populated_nodes
        x_prior=xw
        x_prior[windex] = self.river_obs.centerline_obs['x_prior'].v
        x_prior = x_prior[self.river_obs.populated_nodes]
        #
        windex = self.river_obs.centerline_obs['y_prior'].populated_nodes
        y_prior=yw
        y_prior[windex] = self.river_obs.centerline_obs['y_prior'].v
        node_index=N.arange(len(y_prior))
        node_index = node_index[self.river_obs.populated_nodes]
        y_prior = y_prior[self.river_obs.populated_nodes]
        reach_index=N.ones(len(node_index))*(reach_idx)
        #
        # Add the observations
        self.river_obs.add_obs('h_noise',self.h_noise)
        self.river_obs.add_obs('h_flg',(self.h_flg>0))
        
        self.river_obs.add_obs('lon',self.lon)
        self.river_obs.add_obs('lat',self.lat)
        
        self.river_obs.add_obs('xobs',self.x)
        self.river_obs.add_obs('yobs',self.y)
        self.river_obs.add_obs('inundated_area',self.inundated_area)
        self.river_obs.load_nodes(['h_noise','h_flg',
                                  'lon','lat','xobs','yobs','inundated_area'])
        if self.verbose: print('Observations added to nodes')

        # Get various node statistics

        nobs = N.asarray( self.river_obs.get_node_stat('count','') )
        
        s_median = N.asarray( self.river_obs.get_node_stat('median','s') )
        
        x_median = N.asarray( self.river_obs.get_node_stat('median','xobs') )
        y_median = N.asarray( self.river_obs.get_node_stat('median','yobs') )
        
        lon_median = N.asarray( self.river_obs.get_node_stat('median','lon') )
        lat_median = N.asarray( self.river_obs.get_node_stat('median','lat') )
        ds = N.asarray( self.river_obs.get_node_stat('value','ds') )
        
        # number of good heights
        nobs_h = N.asarray( self.river_obs.get_node_stat('countGood','h_flg') )
        # heights using same pixels as widths
        h_noise_ave0 = N.asarray( self.river_obs.get_node_stat('median','h_noise') )
        h_noise_std0 = N.asarray( self.river_obs.get_node_stat('std','h_noise') )
        # heights using only "good" heights
        h_noise_ave = N.asarray( self.river_obs.get_node_stat('median','h_noise','h_flg') )
        h_noise_std = N.asarray( self.river_obs.get_node_stat('std','h_noise','h_flg') )
        
        
        
        # The following are estimates of river width
        
        width_ptp =  N.asarray( self.river_obs.get_node_stat('width_ptp','') )
        width_std =  N.asarray( self.river_obs.get_node_stat('width_std','') )
        
        # These are area based estimates, the nominal SWOT approach
        
        area = N.asarray( self.river_obs.get_node_stat('sum','inundated_area') )
        width_area = N.asarray( self.river_obs.get_node_stat('width_area','inundated_area') )
        
        # These are the values from the width database
        
        width_db = N.ones(self.river_obs.n_nodes,dtype=N.float64)*self.river_obs.missing_value
       
        #if use_width_db:
        try:
            windex = self.river_obs.centerline_obs['max_width'].populated_nodes
            width_db[windex] = self.river_obs.centerline_obs['max_width'].v
            width_db = width_db[self.river_obs.populated_nodes]
        except:
            width_db = N.ones(self.river_obs.n_nodes,dtype=N.float64)*self.river_obs.missing_value
            
        # Get the start and end
        #smin = s_median.min()
        #smax = s_median.max()
        
        # Initialize the fitter
        
        self.fitter = FitRiver(self.river_obs)
        
        # Check to see if there are sufficient number of points for fit
        ngood = len(s_median)
        if self.verbose: print('number of fit points: %d'%ngood)
        
        if ngood < min_fit_points:
            nresults = None
            if self.verbose: print('not enough good points for fit')
            return None #only process reaches with min number of nodes
        else:
            # Get the start and end
            smin = s_median.min()
            smax = s_median.max()    
            # Do the fitting for this reach
            
            nresults = self.estimate_height_slope(smin,smax,fit_types=fit_types,mean_stat='median')
            if self.verbose: print('Estimation finished')
        
        if self.store_fits:
            self.fit_collection[reach_id,'noise'] = nresults

        # Get the reach statistics for this subreach
        reach_stats = self.get_reach_stats(reach_id,reach_idx,smin,smax,
                                           lon_median, lat_median, 
                                            width_ptp,width_std,width_area,area,
                                            nresults)
        # Put the results into a RiverReach object
        river_reach = RiverReach(lat=lat_median,
                                 lon=lon_median,
                                 x = x_median,
                                 y = y_median,
                                 nobs = nobs,
                                 s = s_median,
                                 ds = ds,             
                                 w_ptp = width_ptp,
                                 w_std = width_std,
                                 w_area = width_area,
                                 w_db = width_db,
                                 area = area,                              
                                 h_n_ave = h_noise_ave,
                                 h_n_std = h_noise_std,
                                 h_a_ave = h_noise_ave0,
                                 h_a_std = h_noise_std0,
                                 nobs_h = nobs_h,
                                 x_prior=x_prior,
                                 y_prior=y_prior,
                                 node_indx=node_index,
                                 reach_indx=reach_index,
                                 metadata = reach_stats)

        # Store, if desired
        if self.store_obs:
            self.river_obs_collection[reach_idx] = self.river_obs
        return river_reach

    def estimate_height_slope(self,smin,smax,fit_types=['OLS','WLS','RLM'],mean_stat='median'):
        """Get fit results for this subreach."""
        
        if type(fit_types) == str:
            fit_types = [fit_types]

        nresults = odict()
        load_inputs=True
        for fit_type in fit_types:
            nresults[fit_type] = self.fitter.fit_linear(smin,smax,'h_noise',
                                                       fit=fit_type,
                                                       mean_stat=mean_stat,
                                                       load_inputs=load_inputs)
            load_inputs = False
        return nresults#, tresults

    def get_reach_stats(self,reach_id,reach_idx,smin,smax,
                        lon_median, lat_median, 
                        width_ptp,width_std,width_area,area,
                        nresults):
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
        
        # fitting results for h_true
        
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
        
        return reach_stats

    # Brent Williams, May 2017: added these functions to output node/reach/slantrangeimage indices to a file 
    def createIndexFile(self):#fname,img_x,img_y,node_index,reach_index):
        """
        Create an empty file with appropriate variable to append to 
        """
        # if the filename does not exist create it
        f=nc.Dataset(self.output_file,'w')
        f.createDimension('record',None)
        ri=f.createVariable('range_index','i4','record',fill_value=-128)
        ai=f.createVariable('azimuth_index','i4','record',fill_value=-128)
        ni=f.createVariable('node_index','i4','record',fill_value=-128)
        rri=f.createVariable('reach_index','i4','record',fill_value=-128)
        seglbl=f.createVariable('segmentation_label','i4','record',fill_value=-128)
        hflg=f.createVariable('good_height_flag','i1','record',fill_value=-1)
        d=f.createVariable('distance_to_node','f4','record',fill_value=-128)       
        f.close()
        return
    def writeIndexFile(self,img_x,img_y,node_index,dst,reach_index,seg_lbl,h_flg):
        """
        Write out the river obs indices for each pixel that get mapped to a node 
        as well as the pixel cloud coordinates (range and azimuth, or original 
        image coordinate [e.g., gdem along- and cross-track index])
        """
        if self.verbose: print('Writing to Index File: ',self.output_file)
        # if the filename does not exist create it
        f=nc.Dataset(self.output_file,'a')
        # first get the data
        ri=f.variables['range_index']
        ai=f.variables['azimuth_index']
        ni=f.variables['node_index']
        rri=f.variables['reach_index']
        seglbl=f.variables['segmentation_label']
        hflg=f.variables['good_height_flag']
        d=f.variables['distance_to_node']
        # append the new data
        L=len(ri[:])
        K=len(img_x)
        ri[L:L+K]=img_x
        ai[L:L+K]=img_y
        ni[L:L+K]=node_index
        rri[L:L+K]=reach_index
        if not(type(seg_lbl) == type(None)):
            seglbl[L:L+K]=seg_lbl
        hflg[L:L+K]=h_flg
        d[L:L+K]=dst
        f.close()
        return
        
