from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    name = 'RivWidthHeight',
    author = "Ernesto Rodriguez",
    version = "0.1",
    description = "Get width and height of rivers from interferometric data.",
    ext_modules=[ 
        Extension("RivWidthHeight", ["RivWidthHeight.pyx"], include_dirs = [numpy.get_include()]),
        Extension("label_regions", ["label_regions.pyx"], include_dirs = [numpy.get_include()]),
        ],
    cmdclass = {'build_ext': build_ext}
)

#
# to build: python setup.py build_ext --inplace
#
