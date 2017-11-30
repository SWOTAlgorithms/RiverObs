#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

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

exec(open('src/GWDLR/version.py').read())
setup(name='GWDLR',
      version=__version__,
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
script_dir = 'src/bin/'
setup(name='SWOTRiver',
      version=__version__,
      description='Extract hydrology observables from SWOT data.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['SWOTRiver',],
      scripts=[script_dir+'estimate_swot_river.py',]
     )

exec(open('src/RiverObs/version.py').read())
#script_dir = 'src/SWOTRiver/scripts/'
setup(name='RiverObs',
      version=__version__,
      description='Associate observations with river reaches and nodes.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['RiverObs',],
      ## scripts=[script_dir+'make_simulation_catalog.py',
      ##         script_dir+'estimate_swot_rivers.py',]
     )

exec(open('src/RDF/version.py').read())
setup(name='RDF',
      version=__version__,
      description='Read, write various flavours of RDF format files.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['RDF',],
      ## scripts=[script_dir+'']
     )

exec(open('src/toggle_input/version.py').read())
setup(name='toggle_input',
      version=__version__,
      description='Toggle notebook input cells.',
      author='Ernesto Rodriguez',
      author_email='ernesto.rodriguez@jpl.nasa.gov',
      ##      url='http://www.python.org/sigs/distutils-sig/',
      package_dir = {'': 'src'},
      ## packages = find_packages()
      packages=['toggle_input',],
      ## scripts=[script_dir+'']
     )
