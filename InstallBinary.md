# RiverObs Installation Instructions


##Preliminaries

These are the instructions for installing the RiverObs **binary** package
written by Ernesto Rodriguez in a Unix (linux or Mac) machine with an
[anaconda](https://store.continuum.io/cshop/anaconda) python setup.

##Python virtual environment installation

Note that the dependence on scikit-image is optional and
required only if one wants to vectorize GWDLR data. In
that case, a working grass installation is required (tested
with grass 6.4; grass70 beta has a bug in r.to.vector as of
this writing).

###Setting up an anaconda virtual environment (Simplest)

To create an anaconda virtual environment in a preselected directory
(here, pointed by the environment variable RIVER_DIR):

	cd $RIVER_DIR
	conda create -p $RIVER_DIR/anaconda numpy ipython ipython-notebook
	matplotlib gdal scipy pip scikit-image statsmodels pysal pandas
	pytables shapely netcdf4 sphinx  numpydoc cython

or, to create an environment in the user's anaconda distribution  (Simplest)

	conda create -n RiverObsBinary numpy ipython ipython-notebook matplotlib
	gdal scipy pip scikit-image statsmodels pysal pandas pytables
	shapely netcdf4 sphinx numpydoc cython
	
To activate this environment, type

	source activate $RIVER_DIR/anaconda

or 

	source activate RiverObsBinary
	
if anaconda/bin is in your path. Otherwise, use /path/to/anaconda/bin/source.

To deactivate this environment, type

	source deactivate

Equivalnetly, one can set

	export PATH=$RIVER_DIR/anaconda/bin:$PATH


##Build additional package requirements

In addition to the packages installed by conda, pyproj, rtree are
required, and they are not directly part of the official anaconda
distribution.  There are two ways to install these packages: either
through conda and [binstar](https://binstar.org/) or by using (inside the conda
virtual environment) pip and compilation.

###Install Using Binstar (Simplest)

There are multiple versions of these packages in binstar for linux or
osx architectures. The user can select which one he wants to install
by searching in the binstar site. The following instructions are valid
as of December 30, 2014.

For osx or linux installation of pyproj:

	conda install -c https://conda.binstar.org/pingucarsti pyproj

For osx or linux installation of rtree and the required library
libspatial index:

	conda install -c https://conda.binstar.org/dougal libspatialindex
	conda install -c https://conda.binstar.org/dougal rtree

###Install Using pip and brew or manual compilation

Working inside the virtual environment, the following command:

	pip install pyproj
	pip install rtree

In addition, [rtree](https://github.com/Toblerity/rtree) requires the
[libspatialindex](http://libspatialindex.github.io) library. On a Mac with
[brew](http://brew.sh) this can be done easily:

	brew install spatialindex

On a generic Unix system, this can be done by downloading the code from
[osgeo](http://download.osgeo.org/libspatialindex), and following the
usual Unix install process. To avoid pat/ownership conflict, one can install into the
anaconda installation:

	tar xvzf spatialindex-src-1.8.1.tar.gz
	cd spatialindex-src-1.8.1
	./configure --prefix=~/anaconda
	make
	make install

###Install numpydoc for sphinx documentation (Optional)

This is only required if you want to build the sphinx documentation:

	pip install numpydoc

##Install the binary package





