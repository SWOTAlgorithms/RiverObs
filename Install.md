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

some thrid party packages may have trouble with the newer python 3.6, if you have trouble you can try with 3.5.  Also, it may be necessary to use version 8d version of jpeg.  If so try the following:

     conda create -n RiverObs python=3.5 numpy jupyter notebook matplotlib gdal scipy pip scikit-image statsmodels pysal pandas pytables shapely netcdf4 sphinx  numpydoc rtree pyproj jpeg=8d

Here is what I got working on a linux box with all the versions explicitly stated:

     conda create -n RiverObs python=3.5 numpy=1.13.1 jupyter=1.0.0 notebook=5.0.0 matplotlib=2.0.2 gdal=2.1.0 libgdal=2.1.0 scipy=0.19.1 pip=9.0.1 scikit-image=0.13.0 statsmodels=0.8.0 pysal=1.13.0 pandas=0.20.3 pytables=3.4.2 shapely=1.5.16 netcdf4=1.2.4 sphinx=1.6.3 numpydoc=0.7.0 rtree=0.8.3 pyproj=1.9.5.1 jpeg=8d

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
