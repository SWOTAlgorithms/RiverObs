"""Associate observations with a river."""

from __future__ import absolute_import

from .RiverNode import RiverNode
from .RiverObs import RiverObs
from .WidthDataBase import WidthDataBase
from .IteratedRiverObs import IteratedRiverObs
from .LatLonRegion import LatLonRegion
# from .ReachPreProcessor import ReachPreProcessor
from .RiverReach import RiverReach
try:
    from .RiverReachWriter import RiverReachWriter
except ModuleNotFoundError as e:
    print("Warning: RiverReachWriter disabled.")
    print("please install gdal if you want to use it")
    print("Warning:", e)
from .version import __version__
