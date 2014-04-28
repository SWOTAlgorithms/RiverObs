#SWOTRiver

This is a package written initially by
[Ernesto Rodriguez](mailto:ernesto.rodriguez@jpl.nasa.gov)  to estimate
various river parameters starting from SWOT L2 data.

Detailed installation instructions are in the Install.md file.

##Summary of packages provided

**SWOTRiver**: This is the main package for handling SWOT data,
associating data with river reaches, and estimating hydrology
parameters base on reach averaging (or not...). In addition to the
homegrown packages listed below, this package requires the following
open source packages:

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

**Centerline**: Provides a class that can be used to project data 
   or refine a river center line. Requires the following packages:

* [scipy](http://www.scipy.org/): Science algorithms swiss army knife. 
* [numpy](http://www.scipy.org/): Numerics swiss army knife. 

**GeometryDataBase**: Find quickly which reach intersects with a 
   geometry of interest. The geometries are assumed to be stored in a 
   shafile. Requires the following packages:

* [Rtree](https://github.com/Toblerity/rtree): Fast bounding box queries. 
* [pysal](http://pysal.org): nice interface to shapefiles and 
      shapely bridge. 
* [shapely](https://github.com/sgillies/shapely): geometry 
      calculations. 

**GDALOGRUtilities**: Provides homegrown utilities for reading and writing
   various GIS files. Requires the following packages:

* [gdal](http://www.gdal.org): GIS files swiss army knife.
* [pyproj](http://code.google.com/p/pyproj): Cartographic
      projections swiss army knife.

**GWDLR**: This is an optional package to convert Global Width 
   Database-Large Rivers raster data provided by 
   [Dai Yamazaki](mailto:bigasmountain1022@gmail.com)  to vectors that can be used as 
   centerlines. Requires:

* [grass](grass.osgeo.org): for raster to vector program. 
* [scikit-image](http://scikit-image.org): for skeletonize. 




