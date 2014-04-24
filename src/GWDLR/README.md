#Global Width Database for Large Rivers (GWDLR)

This package contains for turning GWDLR data set tiles into vector centerlines.

The format of the GWDLR is (Dai Yamazaki personal communication):

*The width data is prepared as 5degx5deg tiles. It's 4 byte real plain binary data, and you can find the details in .ctl files. Note that width value is not accurate around 60N because it's on the domain boundary. The value >0 represents river width (on river centerline), value=-1 represents non-centerline waterbody, value=-2 represents islands. value=0 represents land, value=-9999 represents sea.*

The data are in [GrADS](http://www.iges.org/grads/gadoc/aboutgriddeddata.html) format, but a simpleminded GrADS
parser is part of this package.

##Requirements

To go from raster to vectors, [GRASS](http://grass.osgeo.org/) is required. This version was tested
with grass64, Notice that grass70 version 7.0 (beta1) has a broken
version of thin.

##Setting up a Grass GISDBASE

To use the raster to vector tool, the data will have to be imported to grass.
As a first step an empty repository will have to be created. To avoid this,
the SWOTRiver/data directory contains a GISDBASE with a location called TEMP
that can be used as a starting point. Copy this grassdata directory to the
place that will contain all of the grass output.

The following are the steps to create an GISDBASE:

1. Create a directory to hold it:

   mkdir grassdata

2. Start grass in text version

   grass -text

3. The follwing screen pops up:


                                    GRASS 6.4.2
        
        DATABASE: A directory (folder) on disk to contain all GRASS maps and data.
        
        LOCATION: This is the name of a geographic location. It is defined by a
                  co-ordinate system and a rectangular boundary.
        
        MAPSET:   Each GRASS session runs under a particular MAPSET. This consists of
                  a rectangular REGION and a set of maps. Every LOCATION contains at
                  least a MAPSET called PERMANENT, which is readable by all sessions.
        
                 The REGION defaults to the entire area of the chosen LOCATION.
                 You may change it later with the command: g.region
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        LOCATION:   <UNKNOWN>________________  (enter list for a list of locations)
        MAPSET:     <UNKNOWN>________________  (or mapsets within a location)
        
        DATABASE: /home/erodrigu_______________________________________________________
        
        
        
                   AFTER COMPLETING ALL ANSWERS, HIT <ESC><ENTER> TO CONTINUE
                                    (OR <Ctrl-C> TO CANCEL)
        

replace the <UNKNOWN> location with temp and the MASET with PERMANENT. Modify the
DATABASE to where the grassdata directory is going to be.

4. Create the location and follow the instructions to create an x,y dataset with unit
resolution and extent. End by exiting grass.

##Grass steps to create vector data sets

Assuming grass v6.4 and the working directory is grassdata, start
grass on the mac using:

	/Applications/GRASS-6.4.app/Contents/MacOS/grass.sh -text grassdata/temp/PERMANENT/

in linux:

	grass -text grassdata/n35w125/PERMANENT/

Although the following is implemented in Python, for the record, below
are the grass commands required to get to vectors from the sample data:

	# Start grass on the temp dataset

	grass -text grassdata/temp/PERMANENT

	# Start a new location with the mask data

	g.proj -c georef=n35w125_wth_center_line_25.tif location=n35w125 --verbose

	# Switch to the new location

	g.mapset mapset=PERMANENT location=n35w125

	# Read the mask data into the location

	r.in.gdal input=n35w125_wth_center_line_25.tif output=center_line_mask

	# Thin the mask to make sure that r.to.vect does not crash

	r.thin input=center_line_mask output=center_line_mask_thin --verbose

	# Make a vector from the thinned data

	r.to.vect input=center_line_mask_thin output=center_line feature=line --verbose

	# Export to a shapefile

	v.out.ogr input=center_line dsn=n35w125_center_lines
