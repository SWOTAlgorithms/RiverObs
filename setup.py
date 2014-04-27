#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import numpy as N

setup(name='RDF',
      version='1.0',
      description='Read RDF files',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['RDF',]
      ## scripts=[script_dir+'binary_to_netcdf.py']
     )

exec(open('src/GDALOGRUtilities/version.py').read())
setup(name='GDALOGRUtilities',
      version=__version__,
      description='Utilities for interacting with GDAL and OGR',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['GDALOGRUtilities',]
      ## scripts=[script_dir+'binary_to_netcdf.py']
     )

exec(open('src/Centerline/version.py').read())
setup(name='Centerline',
      version=__version__,
      description='Project coordinates to a curved coordinate system.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['Centerline',]
      ## scripts=[script_dir+'binary_to_netcdf.py']
     )

setup(name='GWDLR',
      version='0.1',
      description='Read and process Global Width Database for Large River data.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['GWDLR',]
      ## scripts=[script_dir+'binary_to_netcdf.py']
     )

exec(open('src/GeometryDataBase/version.py').read())
setup(name='GeometryDataBase',
      version=__version__,
      description='Find geometries within bounding boxes.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['GeometryDataBase',]
      ## scripts=[script_dir+'binary_to_netcdf.py']
     )

exec(open('src/SWOTRiver/version.py').read())
setup(name='SWOTRiver',
      version=__version__,
      description='Extract hydrology observables from SWOT data.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['SWOTRiver',]
      ## scripts=[script_dir+'binary_to_netcdf.py']
     )
