.. _centerline-usage:

Centerline Usage Example
========================

This example show how to intialize a centerline using a shapefile and to
add point river location measurements.

Below are the Python modules required for this example:

.. code:: python

    from os.path import join
    import numpy as N
    import pysal
    from SWOTRiver import SWOTL2
    from Centerline import Centerline
    
    # For plotting
    
    %pylab inline

.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib


Read the example data
---------------------

In this step, the data are read, a projection to (x,y) coordinates is
defined, and the data bounding box is found.

The data locations are given below:

.. code:: python

    # This is the example data
    
    data_dir = '../../data/examples/'
    l2_file = join(data_dir,'simulated_swot_test_data.nc')
    centerline_file = join(data_dir,'sacramento_centerline.shp')
Read the simulated SWOT data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simplified version of the SWOT simulator data containing only
latitudes, longitudes, water classification, true and simulated noisy
height is provided in the ``data/examples`` directory.

For this example, only data with classification=1 labels are used (pure
water pixels) and the latitude and longitude are assumed to be known
from the reference interferogram. The following lines read the
latitude-longitude data for the water pixels, and projects it to a
Lambert Equiarea projection.

.. code:: python

    class_list=[1]
    lat_kwd='no_layover_latitude'
    lon_kwd='no_layover_longitude'
    
    l2 = SWOTL2(l2_file,class_list=class_list,lat_kwd=lat_kwd,lon_kwd=lon_kwd)

.. parsed-literal::

    Dataset opened
    Bounding box calculated
    Good data selected
    lat/lon read
    projection set and x,y calculated


Read a candidate centerline
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A candidate centerline for part of the Sacramento River is provided in
the ``data/examples`` directory courtesy of `Tamlin
Pavelsky <mailto:pavelsky@unc.edu>`__ and `George
Allen <mailto:georgehenryallen@gmail.com>`__ and modified by `Ernesto
Rodriguez <mailto:ernesto.rodriguez@jpl.nasa.gov>`__ so that reaches
appeared as topologically connected.

The shapefile is read using the pysal package and the latitude and
longitude are extracted.

.. code:: python

    shp = pysal.open(centerline_file)
    cline = shp[0]
    lon, lat = N.asarray(cline.vertices).T
To get true distances, the centerline must be given points so that the
Euclidean distance can be calculated. Below, the SWOTL2 projection
function is used to project to the same projection as the SWOT data.

.. code:: python

    xc, yc = l2.proj(lon,lat)
To see what has been done, plot the centerline points and the measurement locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The centerline is plotted as black points, while the water measurements
are plotted as blue points. As can be seen, the simulated measurements
have some gaps over the river and also include water points outside the
river. Notice also that the measurements show meanders not present in
the candidate centerline: this defect is corrected by the class
IteratedRiverObs, discussed later on.

.. code:: python

    figsize(6,6)
    plot(l2.x/1.e3,l2.y/1.e3,'.b',alpha=0.01)
    plot(xc/1.e3,yc/1.e3,',k',alpha=1)
    xlim(-5,5)
    ylim(-10,10)
    xlabel('X (km)')
    ylabel('Y (km)')



.. parsed-literal::

    <matplotlib.text.Text at 0x10b97c410>




.. image:: CenterlineExample_files/CenterlineExample_11_1.png


Initialize the Centerline instance
----------------------------------

The following step shows how the Centerline is initialized with default
parameters.

.. code:: python

    centerline = Centerline(xc,yc)
The following step shows how to associate the simulated measurements
with centerline node locations and assign along-track and normal
coordinates to each point. *(Note that instead of calling the instance,
the member function to\_centerline could also have been used.)*

.. code:: python

    index, distance,x,y,s,n = centerline(l2.x,l2.y)
Plot the distribution of distances to the input centerline nodes
(*distance*), as well as the normal coordinate (*n*), for each of the
SWOT simulated data points.

.. code:: python

    figsize(10,5)
    subplot(1,2,1)
    hist(distance,bins=100,log=True)
    xlabel('Distance to node (m)')
    ylabel('N observations')
    grid();
    subplot(1,2,2)
    hist(n,bins=100,log=True)
    xlabel('Normal coordinate (m)')
    ylabel('N observations')
    grid();


.. image:: CenterlineExample_files/CenterlineExample_17_0.png


Notice that most of the data points are close to the centerline, but, as
expected from the data picture, some of the points are far away. These
points can be filtered using the RiverObs class.

Below is a zoom around the centerline.

.. code:: python

    figsize(10,5)
    subplot(1,2,1)
    hist(distance,bins=arange(0,500,10),log=False)
    xlabel('Distance to node (m)')
    ylabel('N observations')
    grid();
    subplot(1,2,2)
    hist(n,bins=arange(-500,500,50),log=False)
    xlabel('Normal coordinate (m)')
    ylabel('N observations')
    grid();


.. image:: CenterlineExample_files/CenterlineExample_19_0.png


Below is a plot of the normal coordinate ploted as a function of the
reach distance along the centerline. Where the centerline and the river
measurements agree, one can estimate the river width. The missed
meanders in the centerline are also easily identified.

.. code:: python

    plot(centerline.s[index]/1.e3,n,'.',alpha=0.1)
    xlim(40,90)
    ylim(-500,500)
    xlabel('Centerline Reach (km)')
    ylabel('Normal coordinate (m)');


.. image:: CenterlineExample_files/CenterlineExample_21_0.png

