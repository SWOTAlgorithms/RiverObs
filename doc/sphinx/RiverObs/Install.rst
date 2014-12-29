.. _Installation:

RiverObs Installation Instructions
==================================

Preliminaries
-------------

These are the instructions for installing the SWOTRiver package written
by Ernesto Rodriguez in a Unix (linux or Mac) machine with an
`anaconda <https://store.continuum.io/cshop/anaconda>`__ python setup,
or using a `virtualenv <http://www.virtualenv.org/en/latest>`__ with
another python installation. In both cases, it is assumed that
`numpy <http://scipy.org>`__ is available. In addition, to use the
ipython notebooks, it is assumed that `ipython <http://ipython.org>`__
is also installed. If you do not have a lot of python experience, it is
recommended that you follow the anaconda install route.

In what follows, it is assumed that the environment variable RIVER\_DIR
has been set to point to the root directory of the
AirSWOTAnomalyAnalysis package cloned by git. For instance, using bash

::

    export RIVER_DIR=/home/erodrigu/SWOT/RiverObs

Python virtual environment installation
---------------------------------------

Note that the dependence on scikit-image is optional and required only
if one wants to vectorize GWDLR data. In that case, a working grass
installation is required (tested with grass 6.4; grass70 beta has a bug
in r.to.vector as of this writing).

Setting up an anaconda virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create an anaconda virtual environment, execute:

::

    cd $RIVER_DIR
    conda create -p $RIVER_DIR/anaconda numpy ipython ipython-notebook matplotlib gdal scipy pip scikit-image statsmodels pysal pandas pytables

or conda create -n SWOTRiver numpy ipython ipython-notebook matplotlib
gdal scipy pip scikit-image statsmodels pysal pandas pytables

To activate this environment, type

::

    source activate $RIVER_DIR/anaconda

or

::

    source activate SWOTRiver

if anaconda/bin is in your path. Otherwise, use
/path/to/anaconda/bin/source.

To deactivate this environment, type

::

    source deactivate

Equivalnetly, one can set

::

    export PATH=$RIVER_DIR/anaconda/bin:$PATH

Setting up a python virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case you have a prexisting python environment with virtualenv, numpy,
and ipython already installed, create a virtual environment for this
project as follows:

::

    virtualenv --system-site-packages $RIVER_DIR

To activate this environment, type

::

    source $RIVER_DIR/bin/activate

and to deactivate

::

    source deactivate

Build the package
-----------------

In addition to the packages installed by conda, pyproj, rtree, and
shapely are required. Working inside the virtual environment, the
following command:

::

    pip install pyproj
    pip install shapely
    pip install rtree

In addition, `rtree <https://github.com/Toblerity/rtree>`__ requires the
`libspatialindex <http://libspatialindex.github.io>`__ library. On am
Mac with `brew <http://brew.sh>`__ this can be done easily:

::

    brew install spatialindex

On a generic Unix system, this can be done by downloading the code from
`osgeo <http://download.osgeo.org/libspatialindex>`__, and following the
usual Unix install process. To avoid pat/ownership conflict, one can
install into the anaconda installation:

::

    tar xvzf spatialindex-src-1.8.1.tar.gz
    cd spatialindex-src-1.8.1
    ./configure --prefix=~/anaconda
    make
    make install

Then, to build the SWOTRiver and associated packages:

::

    cd $RIVER_DIR
    python setup.py install --force

For an anaconada local virtual environment, this will install the
libraries in

::

    $RIVER_DIR/anaconda/python2.7/site-packages

and the executables in

::

    $RIVER_DIR/anaconda/bin

Otherwise, they are in similar directories in ~/anaconda/envs/SWOTRiver

For a virtualenv virtual environment, this will install the libraries in

::

    $RIVER_DIR/lib/python2.7/site-packages

and the executables in

::

    $RIVER_DIR/bin

