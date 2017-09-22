# RiverObs Installation Instructions

    declare -x PYTHONPATH=/path/to/RiverObs/src/:$PYTHONPATH

or do the following


## Preliminaries

These are the instructions for installing the RiverObs package
written by Ernesto Rodriguez in a Unix (linux or Mac) machine with an
[anaconda](https://store.continuum.io/cshop/anaconda) python setup.
The nominal installation instructions have been tested with python3.6,
but should also work with python2.7. Future developments may stop
supporting python2.7, as it is no longer the community standard.

In what follows, it is assumed that the environment variable RIVER_DIR has been
set to point to the root directory of the RiverObs package cloned
by git. For instance, using bash

	export RIVER_DIR=/home/erodrigu/SWOT/RiverObs

## Python virtual environment installation

Note that the dependence on scikit-image is optional and
required only if one wants to vectorize GWDLR data. In
that case, a working grass installation is required (tested
with grass 6.4; grass70 beta has a bug in r.to.vector as of
this writing).

### Setting up an anaconda virtual environment

To make sure that you are retrieving the same version packages as have
been used for testing, make sure that the conda-forge channel is added
to your conda configuration. This can be done by issuing the command

    conda config --add channels conda-forge

or modifying your ~/.condarc file to look something like this:

    channels:
      - conda-forge
      - defaults
    show_channel_urls: true

To create an anaconda virtual environment, execute (Simplest):

    conda create -n RiverObs python=3.6 numpy jupyter notebook matplotlib
    gdal scipy pip scikit-image statsmodels pysal pandas pytables
    shapely netcdf4 sphinx  numpydoc rtree pyproj

or, if you want to keep the code and executables under the RiverObs folder:

    cd $RIVER_DIR
    conda create -p $RIVER_DIR/anaconda python=3.6 numpy jupyter notebook matplotlib
    gdal scipy pip scikit-image statsmodels pysal pandas pytables
    shapely netcdf4 sphinx  numpydoc rtree pyproj

Note: if you must run python 2.7, substitute python=2.7 in the lines above
(not recommended).

To activate this environment, if the first option was used, type

	source activate RiverObs

or, if building in the RiverObs folder,

    source activate $RIVER_DIR/anaconda

if anaconda/bin is in your path. Otherwise, use /path/to/anaconda/bin/source.

To deactivate this environment, type

	source deactivate

If you would like to use jupyter notebooks within the RiverObs environment,
issue the following command while inside the environment:

    python -m ipykernel install --user

## Build the package

Then, to build the RiverObs and associated packages:

	cd $RIVER_DIR
	python setup.py install --force

For an anaconada local virtual environment, this will install the libraries in

	$RIVER_DIR/anaconda/python3.6/site-packages

and the executables in

	$RIVER_DIR/anaconda/bin

Otherwise, they are in similar directories in ~/anaconda/envs/RiverObs
