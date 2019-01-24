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
