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

	conda create -p $RIVER_DIR/RiverObsBinary --file RiverObsSetupPackages.txt
	
To activate this environment, type

	source activate $RIVER_DIR/RiverObsBinary

To deactivate this environment, type

	source deactivate

Equivalnetly, one can set

	export PATH=$RIVER_DIR/RiverObsBinary/bin:$PATH


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

However, it has been reported that, depending of the linux
distribution, these packages may have problems of library
dependency. If these installations have problems, remove them using

	conda remove libspatialindex
	conda remove rtree

and use the manual installation methods below.

###Install Using pip and brew (OSX) or manual compilation

Working inside the virtual environment, the following command:

	pip install pyproj

[Rtree](https://github.com/Toblerity/rtree) requires the
[libspatialindex](http://libspatialindex.github.io) library. On a Mac with
[brew](http://brew.sh) this can be done easily:

	brew install spatialindex

On a generic Unix system, this can be done by downloading the code from
[osgeo](http://download.osgeo.org/libspatialindex), and following the
usual Unix install process. To avoid ownership conflict, one can install into the
RiverObsBinary installation (e.g.):

	tar xvzf spatialindex-src-1.8.5.tar.gz
	cd spatialindex-src-1.8.5
	./configure --prefix=$RIVER_DIR/RiverObsBinary
	make
	make install

This should put the library in $RIVER_DIR/RiverObsBinary/lib and the
includes in $RIVER_DIR/RiverObsBinary/include. To make sure they are
in the load path, you may have to do the following in *nix systems:

	export LD_LIBRARY_PATH=$RIVER_DIR/RiverObsBinary/lib

To install rtree, it should just be possible to issue the following
command from within the virtual environment:

	pip install rtree

Unfortunately, this does not always work because in some systems
ctypes.util.find_library fails to find the libspatialindex
libraries. If you do have problems, you will have to download and
install the riverobs branch of rtree as follows:

	git clone https://github.com/skinkpad/rtree.git
	cd rtree
	git checkout riverobs
	pip install .

If for some reason that still does not work, point to the location of the libspatialindex_c library as follows
(modify according to your system and installation):

	export LIBSPATIALINDEX=$RIVER_DIR/RiverObsBinary/lib/libspatialindex_c.so

before installing with pip as above.
	

##Install the binary package

For release 1.0.0, RiverObs packages for OSX and 64-bit linux
architectures are provided. Once the appropriate package has been
obtained, the RiverObs binary installation is completed by

	cd $RIVER_DIR
	tar xvzf RiverObsBinary_ARCH_version1.0.0.tgz

where ARCH is either linux64 or OSX. The packages can the be accessed
while working in the RiverObsBinary environment. The the actual
packages provided are:

    tar tzf RiverObsBinary_linux64_version1.0.0.tgz 
    RiverObsBinary/bin/estimate_swot_river.py
    RiverObsBinary/lib/python2.7/site-packages/Centerline/
    RiverObsBinary/lib/python2.7/site-packages/Centerline/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/Centerline/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/Centerline/version.py
    RiverObsBinary/lib/python2.7/site-packages/Centerline/Centerline.so
    RiverObsBinary/lib/python2.7/site-packages/Centerline/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/GDALWriter.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/GDALutilities.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/GeodeticPath.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/CoordinateTransformations.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/OGRWriter.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/version.py
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/GDALLatLonLayer.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/OGR2Shapely.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/GDALInfo.so
    RiverObsBinary/lib/python2.7/site-packages/GDALOGRUtilities/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/GeometryDataBase/
    RiverObsBinary/lib/python2.7/site-packages/GeometryDataBase/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/GeometryDataBase/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/GeometryDataBase/version.py
    RiverObsBinary/lib/python2.7/site-packages/GeometryDataBase/GeometryDataBase.so
    RiverObsBinary/lib/python2.7/site-packages/GeometryDataBase/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/GWDLR2shape.so
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/version.py
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/GWDLR.so
    RiverObsBinary/lib/python2.7/site-packages/GWDLR/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/RDF/
    RiverObsBinary/lib/python2.7/site-packages/RDF/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/RDF/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/RDF/RDF_to_class.so
    RiverObsBinary/lib/python2.7/site-packages/RDF/version.py
    RiverObsBinary/lib/python2.7/site-packages/RDF/MRDF.so
    RiverObsBinary/lib/python2.7/site-packages/RDF/ExecuteRDF.so
    RiverObsBinary/lib/python2.7/site-packages/RDF/RDF.so
    RiverObsBinary/lib/python2.7/site-packages/RDF/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/LatLonRegion.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/RiverNode.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/version.py
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/ReachPreProcessor.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/RiverReach.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/IteratedRiverObs.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/ReachExtractor.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/RiverObs.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/WidthDataBase.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/FitRiver.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/RiverReachWriter.so
    RiverObsBinary/lib/python2.7/site-packages/RiverObs/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/SWOTL2.so
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/version.py
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/version.pyc
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/EstimateSWOTRiver.so
    RiverObsBinary/lib/python2.7/site-packages/SWOTRiver/SWOTRiverEstimator.so
    RiverObsBinary/lib/python2.7/site-packages/toggle_input/
    RiverObsBinary/lib/python2.7/site-packages/toggle_input/__init__.py
    RiverObsBinary/lib/python2.7/site-packages/toggle_input/__init__.pyc
    RiverObsBinary/lib/python2.7/site-packages/toggle_input/toggle_input.so
    RiverObsBinary/lib/python2.7/site-packages/toggle_input/version.py
    RiverObsBinary/lib/python2.7/site-packages/toggle_input/version.pyc

The files are owned by erodrigu. To correct this (which should
not matter) you can do

	chown -R user RiverObsBinary

where user is your user name.
