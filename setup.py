#!/usr/bin/env python
"""
This setup.py script is used to manage the installation of multiple packages
related to RiverObs. The script consolidates the setup configuration for
several packages into a single file to minimize changes versus its original
design. The version of the meta-package (SWOTRiverAlgorithms) matches the
RiverObs version number. You should be able to import subpackages like RDF,
SWOTRiver, RiverObs etc following pip install.
"""

from setuptools import setup, find_packages

def get_version(package):
    version_file = f'src/{package}/version.py'
    with open(version_file) as f:
        exec(f.read())
    return locals()['__version__']

setup(
    name='SWOTRiverAlgorithms',
    version=get_version("RiverObs"),
    description='Meta-package for SWOT River Algorithms',
    author='Ernesto Rodriguez',
    author_email='ernesto.rodriguez@jpl.nasa.gov',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    scripts=[
        'src/SWOTRiver/scripts/make_simulation_catalog.py',
        'src/SWOTRiver/scripts/estimate_swot_rivers.py',
    ],
    extras_require={
        'RDF': [f'RDF=={get_version("RDF")}'],
        'GDALOGRUtilities': [f'GDALOGRUtilities=={get_version("GDALOGRUtilities")}'],
        'Centerline': [f'Centerline=={get_version("Centerline")}'],
        'GWDLR': [f'GWDLR=={get_version("GWDLR")}'],
        'GeometryDataBase': [f'GeometryDataBase=={get_version("GeometryDataBase")}'],
        'SWOTWater': [f'SWOTWater=={get_version("SWOTWater")}'],
        'SWOTRiver': [f'SWOTRiver=={get_version("SWOTRiver")}'],
        'RiverObs': [f'RiverObs=={get_version("RiverObs")}'],
        'toggle_input': [f'toggle_input=={get_version("toggle_input")}'],
    },
    install_requires=[
        line.strip() for line in open('requirements.txt')
    ],
    
    
)
