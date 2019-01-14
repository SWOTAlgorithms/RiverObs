"""
Base class containing river reach coordinates and metadata.
"""
from __future__ import absolute_import, division, print_function

class RiverReach:
    """Base class containing river reach coordinates and metadata.

    Initialize with empty information and, optionally,
    a set of keywords that get turned into class members.
    The metadata is stored as dictionary variable_name:value.

    This base class can be derived to implement specific reach behavior.
    """

    def __init__(self, **kwds):
        self.lat = None
        self.lon = None
        self.x = None
        self.y = None
        self.node_indx = None
        self.metadata = {}

        for k in kwds:
            v = kwds[k]
            setattr(self, k, v)
