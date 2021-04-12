![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.png?v=103)
# RiverObs

This is a package written initially by
[Ernesto Rodriguez](mailto:ernesto.rodriguez@jpl.nasa.gov) to estimate
various river parameters starting from remote sensing data.
[Alex Fore](mailto:alexander.fore@jpl.nasa.gov) and [Brent Williams](mailto:brent.a.williams@jpl.nasa.gov) have also provided
code to reflect the evolving SWOT project.
The code is currently maintained by the SWOT Algorithm Definition Team.

Detailed installation instructions are in the Install.md file.

# Usage

For generating data products are most similar to the SWOT project's data products, the following script is recommended (found in src/bin):
```
usage: swot_pixc2rivertile.py [-h] [--shpbasedir SHPBASEDIR] [-l LOG_LEVEL]
                              [--gdem-file GDEM_FILE]
                              pixc_file out_riverobs_file out_pixc_vector_file
                              rdf_file
```
where ```pixc_file``` is the SWOT high-resolution pixel cloud data product, ```out_riverobs_file``` is the filename of the output rivertile data product, ```out_pixc_vector_file``` is the filename of the output pixel cloud vector data product, ```rdf_file``` is the configuration file (see [this link](https://github.com/SWOTAlgorithms/RiverObs/blob/develop/src/bin/swot_pixc2rivertile.py#L13)) for the recomended configuration). Additionally there are some optional arguments: ```--shpbasedir SHPBASEDIR``` will write out the nodes and reaches as shapefile format (written as netCDF to ```out_riverobs_file```), ```-l LOG_LEVEL``` controls the verbosity of the logging, and ```--gdem-file GDEM_FILE``` will create a pixc_file from the GDEM file and run RiverObs on that as a type of truth processing.

# Prior Reach Database
RiverObs requires a prior reach and node database. The database contains fixed node locations, reach boundaries, and high-resolution reach centerlines. It is distributed as a set of netcdf files, broken by continent (first two characters in the file name) and "major basins" in the continent (3rd and 4th characters in the file name). Metadata describing the database fields and the current version of the database is available here http://gaia.geosci.unc.edu/SWORD/SWORD_v08_iceflag.zip. 

## Summary of packages provided

**RiverObs**: This is the main package for associating data with river
reaches, and estimating hydrology parameters base on reach
averaging (or not...). In addition to the homegrown packages listed
below, this package requires the following open source packages:

* [scipy](http://www.scipy.org/): Science algorithms swiss army knife.
* [numpy](http://www.scipy.org/): Numerics swiss army knife.
* [netCDF4](code.google.com/p/netcdf4-python): Reading netcdf4 files,
  including SWOT L2 data files.
* [StatsModels](http://statsmodels.sourceforge.net): Fitting and
  estimation tools.
* [pysal](http://pysal.org): nice interface to shapefiles and
      shapely bridge.  
* [pyproj](http://code.google.com/p/pyproj): Cartographic
      projections swiss army knife.
* [pandas](http://pandas.pydata.org): The Python Data Analysis
  Library for DataFrames and HDFStore.
* [pytables](http://www.pytables.org): easy HDF5 support, required for
  pandas HDFStore.

**Centerline**: Provides a class that can be used to project data
   or refine a river center line. Requires the following packages:

* [scipy](http://www.scipy.org/): Science algorithms swiss army knife.
* [numpy](http://www.scipy.org/): Numerics swiss army knife.

**GeometryDataBase**: Find quickly which reach intersects with a
   geometry of interest. The geometries are assumed to be stored in a
   shapefile. Requires the following packages:

* [Rtree](https://github.com/Toblerity/rtree): Fast bounding box queries.
* [libspatialindex](http://libspatialindex.github.io): Required by Rtree.
* [pysal](http://pysal.org): nice interface to shapefiles and
      shapely bridge.
* [shapely](https://github.com/sgillies/shapely): geometry
      calculations.

**SWOTRiver**: This package contains classes that use the RiverObs
capabilities to produce hydrology outputs from SWOT (simulated) data.

* [numpy](http://www.scipy.org/): Numerics swiss army knife.
* [netCDF4](code.google.com/p/netcdf4-python): Reading netcdf4 files,
  including SWOT L2 data files.
* [pyproj](http://code.google.com/p/pyproj): Cartographic
      projections swiss army knife.
* [pandas](http://pandas.pydata.org): The Python Data Analysis
  Library for DataFrames and HDFStore.
* [pytables](http://www.pytables.org): easy HDF5 support, required for
  pandas HDFStore.

**GDALOGRUtilities**: Provides homegrown utilities for reading and writing
   various GIS files. Requires the following packages:

* [gdal](http://www.gdal.org): GIS files swiss army knife.
* [pyproj](http://code.google.com/p/pyproj): Cartographic
      projections swiss army knife.

**GWDLR**: This is an optional package to convert Global Width
   Database-Large Rivers raster data provided by
   [Dai Yamazaki](mailto:bigasmountain1022@gmail.com)  to vectors that can be used as
   centerlines. Requires:

* [grass](http://grass.osgeo.org): for raster to vector program.
* [scikit-image](http://scikit-image.org): for skeletonize.
