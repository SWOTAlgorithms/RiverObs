import numpy as np
import pyproj
#import RiverObs
from Centerline import Centerline
from RiverObs import ReachDatabase

class RiverProjector:
    """
    This object enables direct projection of lat/lon data to equiarea
    and then river coordinates for specific reaches using the various
    RiverObs functions/objects

    An object providing the following members
    bounding_box:  (lonmin,latmin,lonmax,latmax)
    proj:          a pyproj.Proj projection (lon,lat) -> (x,y)
                   and (x,y) -> (lon,lat) when called when called
                   with inverse=True
    sword_file:    filename of tiled netcdf SWORD file
    reaches:       a list of reaches read in from the SWORD netcdf file
    centerine:     a centerline object from SWORD centerline
    reach_id:      the reach_id for the reach to do projections
    setup_reach(): method to update centerline for specific reach
    project():     method to do the equiarea and river-oriented projections
    
    Added Aug. 2025 by Brent Williams to handle these projections outside
    the normal flow of the river processor (e.g., for projecting all PIXCVec
    pixels back into river coordinates for a given node etc). This is useful
    for debugging PIXC assignment from PIXCVec, or for projecting PIXC or
    other data (like the prior water masks) to the river coordinates for
    various purposes (e.g., like creating river channel shape estimates from
    multitemporal PIXC/PIXCVec data).
    """

    def __init__(self, bbox, sword_file, **proj_kwds):
        self.bounding_box = bbox
        self.sword_file = sword_file
        lat_0 = (bbox[3] + bbox[1]) / 2.0
        lon_0 = self.wrap_lons(
                self.wrap_lons(
                    bbox[2] - bbox[0]) / 2.0 + bbox[0])
        self.proj = pyproj.Proj(
            proj='laea',
            lat_0=lat_0,
            lon_0=lon_0,
            x_0=0,
            y_0=0,
            ellps='WGS84',
            **proj_kwds)
        #self.reaches = RiverObs.ReachDatabase.ReachExtractor(
        self.reaches = ReachDatabase.ReachExtractor(
                self.sword_file, self, day_of_year=None)
        self.reach_id = 0 # default for unset reach

    def setup_reach(self, reach_id):
        self.reach_id = reach_id
        ind = np.where(np.array(
            self.reaches.reach_idx) == int(reach_id))[0][0]
        reach = self.reaches[ind]
        self.centerline = Centerline(reach.x, reach.y)

    def project(self, lats, lons):
        """
        method to project lat/lon points to equiarea then to
        river-centric coordinate
        Inputs:
            lats: latitude array (or 2d image) in degrees north
            lons: longitude array (or 2d image) in degrees east
        
        Outputs:
            index: index of nearest centerline point?
            d:     euclidian distance to nearest centerline point?
            x:     quiarea x coordinate (of nearest centeline point)?
            y:     equiarea y coorediate (of nearest centerline point)?
            s:     along-reach distance (m)
            n:     cross-reach distance (m) (signed, negative to left)?
        """
        # TODO: crop input data around bbox first?
        wrapped_lons = self.wrap_lons(lons)
        xobs, yobs = self.proj(wrapped_lons, lats)
        index, d, x, y, s, n = self.centerline(xobs, yobs)
        return index, d, x, y, s, n

    def wrap_lons(self, lons):
        """
        Wraps lons to [-180,180]
        """
        return  lons - 360 * np.floor((lons + 180) / 360)


