# SWOTRiver Installation Instructions


##0- PRELIMINARIES

These are the instructions for installing the SWOTRiver package
written by Ernesto Rodriguez in a Unix (linux or Mac) machine with an
[anaconda](https://store.continuum.io/cshop/anaconda) python setup,
or using a [virtualenv](http://www.virtualenv.org/en/latest) with another python
installation. In both cases, it is assumed that
[numpy](http://scipy.org) is available. In addition, to use the
ipython notebooks, it is assumed that [ipython](http://ipython.org) is also installed.
If you do not have a lot of python experience, it is recommended that
you follow the anaconda install route.

In what follows, it is assumed that the environment variable RIVER_DIR has been 
set to point to the root directory of the AirSWOTAnomalyAnalysis package cloned
by git. For instance, using bash

	export RIVER_DIR=/home/erodrigu/SWOT/AirSWOT/AirSWOTAnomalyAnalysis

##1- PYTHON VIRTUAL ENVIRONMENT INSTALLATION

###1.1- Setting up an anaconda virtual environment

To create an anaconda virtual environment, execute:

	cd $RIVER_DIR
	conda create -p $RIVER_DIR/anaconda numpy ipython ipython-notebook matplotlib gdal scipy pip

or
	conda create -n SWOTRiver numpy ipython ipython-notebook matplotlib gdal scipy pip 
	
To activate this environment, type

	source activate $RIVER_DIR/anaconda

or 

	source activate SWOTRiver
	
if anaconda/bin is in your path. Otherwise, use /path/to/anaconda/bin/source.

To deactivate this environment, type

	source deactivate

Equivalnetly, one can set

	export PATH=$RIVER_DIR/anaconda/bin:$PATH

###1.2- Setting up a python virtual environment

In case you have a prexisting python environment with virtualenv,
numpy, and ipython already installed, create a virtual environment for
this project as follows:

	virtualenv --system-site-packages $RIVER_DIR

To activate this environment, type

	source $RIVER_DIR/bin/activate

and to deactivate

	source deactivate


##2- BUILD THE PACKAGE

In addition to the packages installed by conda, pyproj is required.
Working inside the virtual environment, the following command:

	pip install pyproj

Then, to build the package

	cd $RIVER_DIR
	python setup.py install --force

For an anaconada local virtual environment, this will install the libraries in

	$RIVER_DIR/anaconda/python2.7/site-packages

and the executables in

	$RIVER_DIR/anaconda/bin

Otherwise, they are in similar directories in ~/anaconda/envs/SWOTRiver 

For a virtualenv virtual environment, this will install the libraries in

	$RIVER_DIR/lib/python2.7/site-packages

and the executables in

	$RIVER_DIR/bin





