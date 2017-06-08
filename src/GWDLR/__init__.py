"""
This package contains for turning GWDLR data set tiles into vector centerlines.

The format of the GWDLR is (Dai Yamazaki personal communication):

"The width data is prepared as 5degx5deg tiles. It's 4 byte real plain binary data,
and you can find the details in .ctl files. Note that width value is not accurate
around 60N because it's on the domain boundary. The value >0 represents river width
(on river centerline), value=-1 represents non-centerline waterbody,
value=-2 represents islands. value=0 represents land, value=-9999 represents sea.""

The data are in [GrADS](http://www.iges.org/grads/gadoc/aboutgriddeddata.html) format, but a simpleminded GrADS
parser is part of this package.
"""

from __future__ import absolute_import

from .GWDLR import GWDLR
from .GWDLR2shape import GWDLR2shape

from .version import __version__
