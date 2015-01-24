#!/usr/bin/env python
"""
Setup for installing RiverObs as Cython shared libraries.
"""

from distutils.core import setup
from distutils.extension import Extension
import numpy as N
from Cython.Build import cythonize

import os.path
import sys
sys.path.append('./src')
from CythonizeProject import CythonizeProject


##exec(open('src/GDALOGRUtilities/version.py').read())
##setup(name='GDALOGRUtilities',
##      version=__version__,
##      description='Utilities for interacting with GDAL and OGR',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['GDALOGRUtilities',]
##      ## scripts=[script_dir+'binary_to_netcdf.py']
##     )
##
##exec(open('src/Centerline/version.py').read())
##setup(name='Centerline',
##      version=__version__,
##      description='Project coordinates to a curved coordinate system.',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['Centerline',]
##      ## scripts=[script_dir+'binary_to_netcdf.py']
##     )
##
##setup(name='GWDLR',
##      version='0.1',
##      description='Read and process Global Width Database for Large River data.',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['GWDLR',]
##      ## scripts=[script_dir+'binary_to_netcdf.py']
##     )
##
##exec(open('src/GeometryDataBase/version.py').read())
##setup(name='GeometryDataBase',
##      version=__version__,
##      description='Find geometries within bounding boxes.',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['GeometryDataBase',]
##      ## scripts=[script_dir+'binary_to_netcdf.py']
##     )
##
##exec(open('src/SWOTRiver/version.py').read())
##script_dir = 'src/SWOTRiver/scripts/'
##setup(name='SWOTRiver',
##      version=__version__,
##      description='Extract hydrology observables from SWOT data.',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['SWOTRiver',],
##      scripts=[script_dir+'make_simulation_catalog.py',
##               script_dir+'estimate_swot_rivers.py',]
##     )
##
##exec(open('src/RiverObs/version.py').read())
###script_dir = 'src/SWOTRiver/scripts/'
##setup(name='RiverObs',
##      version=__version__,
##      description='Associate observations with river reaches and nodes.',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['RiverObs',],
##      ## scripts=[script_dir+'make_simulation_catalog.py',
##      ##         script_dir+'estimate_swot_rivers.py',]
##     )
##

included_dirs = ['./src/RDF']
build_dir = './cython_build'
ignored_files=['__init__.py','version.py']
excluded_dirs=['scripts']
cyp = CythonizeProject(included_dirs,build_dir,ignored_files,excluded_dirs)
modules = cyp.make_extensions()
package_dir = os.path.join(cyp.build_dir,'src')

version_file = os.path.join(package_dir,'RDF','version.py') 
exec(open(version_file).read())

setup(name='RDF',
      version=__version__,
      description='Read, write various flavours of RDF format files.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': package_dir},
      ## packages = find_packages()
      packages=['RDF',],
      
      ## scripts=[script_dir+'']
     )

##setup(name='toggle_input',
##      version=__version__,
##      description='Toggle notebook input cells.',
##      author='Ernesto Rodriguez',
##      author_email='ernesto.rodriguez@jpl.nasa.gov',
##      ##      url='http://www.python.org/sigs/distutils-sig/',
##      package_dir = {'': 'src'},
##      ## packages = find_packages()
##      packages=['toggle_input',],
##      ## scripts=[script_dir+'']
##     )
