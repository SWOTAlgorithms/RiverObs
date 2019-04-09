Here are the main scripts / executables for running RiverObs and ananlyzing its outputs.


# swot_pixc2rivertile.py
```
usage: swot_pixc2rivertile.py [-h] [--shpbasedir SHPBASEDIR] [-l LOG_LEVEL]
                              [--gdem-file GDEM_FILE]
                              pixc_file out_riverobs_file out_pixc_vector_file
                              rdf_file
```
The main river processing script.

# fake_pixc_from_gdem.py
```
usage: fake_pixc_from_gdem.py [-h] [--subsample-factor SUBSAMPLE_FACTOR]
                              pixc_file gdem_file fake_pixc_file
```
Creates a pixel cloud file from a gdem suitable for using in RiverObs processing.  The pixc.nc is a pixel cloud file that is used to determine the correct format for RiverObs.  

# plot_riverobs.py
```
usage: plot_riverobs.py [-h] [-t TITLE] [-p] pixc_rivertile gdem_rivertile
```
Compares the pixel cloud rivertile data product to that which was generated using the fake pixc from the GDEM.  Plots height and width errors for nodes.

# plot_reach_stats.py
```
usage: plot_reach_stats.py [-h] [-t TITLE] [-p] pixc_rivertile gdem_rivertile
```
Compares the pixel cloud and gdem rivertile data products similar to plot_riverobs.py but compares reaches instead.

# preproc_gdem.py
```
usage: preproc_gdem.py [-h] [-l LOG_LEVEL] [--plot]
                       [--erosion-iter EROSION_ITER]
                       in_gdem_file out_gdem_file reachdb_path
```
Pre-processes the GDEM to trim out some non-river / main channel water.

Basically, it works in two steps:
* segmenting the water features into disconnected (or close-to-disconnected) regions
* assigning a particular feature label to each river reach using proximity to the river database.
 
There is an optional parameter called ```--erosion-iter``` that you can fiddle with that can be used to disconnect features that are technically touching, but only barely.  What it does is first erode the water mask before the initial segmentation, then figures out how to handle the things that got eroded in a fancy way. 
 