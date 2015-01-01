"""
Refine an intial set of reaches to iterate the centerline or break the reach into
subreaches of desired characteristics.
"""

from collections import OrderedDict as odict
import numpy as N
from ReachExtractor import ReachExtractor
from WidthDataBase import WidthDataBase
from Centerline import Centerline
from RiverReach import RiverReach


class ReachPreProcessor(ReachExtractor):
    """Refine an intial set of reaches to iterate the centerline or break the reach into
    subreaches of desired characteristics. This is derived class of ReachExtractor.

    On intialization, a candidate set of reaches, which may be much longer than
    desired or be out of date relative to the actual reach, is provided a a shapefile
    of line strings. In addition, a lat_lon_region object is provided 

    Parameters
    ----------
    
    shape_file_root : str
        path to shapefile database (no suffix)
    lat_lon_region : object
        an object satisfying the LatLonRegion protocol providing the following members:
        lat_lon_region.bounding_box: (lonmin,latmin,lonmax,latmax)
        lat_lon_region.proj: a pyproj.Proj projection (lon,lat) -> (x,y)
        and (x,y) -> (lon,lat) when called when called with inverse=True
    clip : bool, optional
        Clip to the bounding box?
    clip_buffer : float, optional
        buffer to add around bounding box.
    width_db_file : str, optional
        The path to a width database file, if one is available.
    ds : float
        If not None, resample x,y points to this spacing 
        (approximately to not extend the interval).
    """
    def __init__(self, shape_file_root, lat_lon_region,clip=True,
                 clip_buffer=0.1,width_db_file=None,ds=None):

        # Initialize the base class
        
        ReachExtractor.__init__(self,shape_file_root, lat_lon_region,clip=clip,
                                clip_buffer=clip_buffer)

        # Remember the lat_lon_region

        self.lat_lon_region = lat_lon_region

        # If a width data base is available, read it and compute maximum width

        if  width_db_file != None:
            self.width_db = WidthDataBase(width_db_file)
        else:
            self.width_db = None

        # Define a center line for each reach

        self.centerline = []
        self.max_width = []
        for i,r in enumerate(self.reach):
            if  self.width_db != None:
                max_width = self.width_db.get_river(self.reach_idx[i],
                                            columns=['width'],
                                            asarray=True,transpose=False,
                                            bounding_box=lat_lon_region.bounding_box,
                                            clip_buffer=clip_buffer).squeeze()
                self.max_width.append(max_width)
                self.centerline.append(Centerline(r.x,r.y,ds=ds,
                                                  obs=[max_width],
                                                  obs_names=['max_width'])
                                                  )
            else:
                self.centerline.append(Centerline(r.x,r.y,ds=ds))

    def split_by_coordinates(self,reach_start_list,reach_end_list,max_distance=None):
        """Split the reaches by a predefined set of reach starts and finishes.

        Parameters
        ----------

        reach_start_list : list
            List of (lon,lat) tuples of reach start coordinates (degrees)
        reach_end_list : list
            List of (lon,lat) tuples of reach start coordinates (degrees)
        max_distance : float
            Maximum distance allowed between start and stop points and the input reach.
            If exceeded, no new reach is generated.

        Returns
        -------

        list of edited reaches.
        """

        self.edited_reach = []
        for i in range(len(reach_start_list)):
            
            # Project to the centerline coordinates
            lon,lat = reach_start_list[i]
            x,y = self.lat_lon_region.proj(lon,lat)

            indexstart, icl, dmin = self.nearest_centerline_node(x,y,max_distance=max_distance)
            
            if icl < 0:
                continue

            lon,lat = reach_end_list[i]
            x,y = self.lat_lon_region.proj(lon,lat)
            indexend, _ , dmin = self.nearest_centerline_node(x,y,max_distance=max_distance,
                                                              cl_index=icl)
            # Extract the reach
            if indexend > -1:
                x = self.centerline[icl].x[indexstart:indexend+1]
                y = self.centerline[icl].y[indexstart:indexend+1]
                lon, lat = self.lat_lon_region.proj(x,y,inverse=True)
                reach_length = self.centerline[icl].s[indexend] - self.centerline[icl].s[indexstart]
                metadata = odict([
                    ('parent_reach_idx',icl),
                    ('parent_start_s',self.centerline[icl].s[indexstart]),
                    ('parent_end_s',self.centerline[icl].s[indexend]),
                    ('indexstart',indexstart),
                    ('indexend',indexend),
                    ('npoints',indexend - indexstart + 1),
                    ('reach_length',reach_length)
                    ])
                if  self.width_db == None:
                    reach = RiverReach(lat=lat,lon=lon,x=x,y=y,metadata=metadata)
                else:
                    width = self.centerline[icl].max_width[indexstart:indexend+1]
                    metadata['width_mean'] = N.mean(width)
                    metadata['width_max'] = N.max(width)
                    metadata['width_min'] = N.min(width)
                    reach = RiverReach(lat=lat,lon=lon,x=x,y=y,metadata=metadata,width=width)
                    
                self.edited_reach.append(reach)
                
        return self.edited_reach


    def nearest_centerline_node(self,x,y,max_distance=None,cl_index=None):
        """Find the nearest centerline and node to a centerline or set of centerlines."""
             
        if cl_index == None:
            CL = self.centerline
        else:
            CL = [self.centerline[cl_index]]
                
        if max_distance == None:
            max_distance = 1.e12
                
        dmin = 1.e12
        icl = -1
        indexmin = -1
        for i,centerline in enumerate(CL):
            index, distance,x,y,s,n = centerline(x,y)
            if (distance < dmin) and (distance < max_distance):
                dmin = distance[0]
                indexmin = index[0]
                icl = i

        return indexmin, icl, dmin

    def split_by_reach_length(self,ds,start_s=0,end_s=None):
        """Split the reaches by a predefined set of reach lengths.

        Parameters
        ----------

        ds : float
            Length of each reach (m)
        start_s : float, optional
            Start the reach split at this length (m).
        end_s : float, optional
            End the reach split at this length or before (m).

        Returns
        -------

        list of edited reaches.
        """

        self.ds,self.start_s,self.end_s = ds,start_s,end_s
        
        self.edited_reach = []
        for icl,reach in enumerate(self.reach):
            s = self.centerline[icl].s

            if end_s != None:
                smax = min(s[-1],end_s)
            else:
                stop_s = s[-1]
            s0 = start_s 
            s1 = s0 + ds

            while s1 <= smax:
                print smax, s0, s1
                i0 = N.flatnonzero(s >= s0 )
                if len(i0) > 0:
                    i0 = i0[0]
                else:
                    break

                i1 = N.flatnonzero(s <= s1 )
                if len(i1) > 0:
                    i1 = i1[-1]
                else:
                    break

                x = self.centerline[icl].x[i0:i1+1]
                y = self.centerline[icl].y[i0:i1+1]
                lon, lat = self.lat_lon_region.proj(x,y,inverse=True)
                reach_length = self.centerline[icl].s[i1] - self.centerline[icl].s[i0]
                metadata = odict([
                    ('parent_reach_idx',icl),
                    ('parent_start_s',self.centerline[icl].s[i0]),
                    ('parent_end_s',self.centerline[icl].s[i1]),
                    ('indexstart',i0),
                    ('indexend',i1),
                    ('npoints',i1 - i0 + 1),
                    ('reach_length',reach_length)
                    ])
                if  self.width_db == None:
                    reach = RiverReach(lat=lat,lon=lon,x=x,y=y,metadata=metadata)
                else:
                    width = self.centerline[icl].max_width[i0:i1+1]
                    metadata['width_mean'] = N.mean(width)
                    metadata['width_max'] = N.max(width)
                    metadata['width_min'] = N.min(width)
                    reach = RiverReach(lat=lat,lon=lon,x=x,y=y,metadata=metadata,width=width)
                    
                self.edited_reach.append(reach)

                s0 = s1
                s1 += ds
        return self.edited_reach
            
                    
        

                
            

        
