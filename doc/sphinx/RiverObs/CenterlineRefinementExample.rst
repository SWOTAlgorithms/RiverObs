.. _centerline-refinement:

Centerline Refinement and RiverObs Data Loading Example
==============================================

This section shows an example of how, given an initial guess for the
Centerline, that centerline can be further refined to match the data
locations using the IteratedRiverObs class. It also shows how, once the
Centerline is refined, measurements can be associated with it.

This example illustrates using the output of the SWOT L2 data (i.e.,
water file) without any supplementary classification information, so
that the true water classification, rather than the estimated water
classification is used.

.. code:: python

    import numpy as N
    from SWOTRiver import SWOTL2
    from RiverObs import ReachExtractor
    from RiverObs import RiverObs
    from RiverObs import FitRiver
    from RiverObs import WidthDataBase
    from RiverObs import IteratedRiverObs

.. code:: python

    # This is the file for the width data base
    
    #width_db_dir = '/Volumes/Reservoir/Data/GWD-LR/nAmerica_GRWDL/'
    width_db_dir = '../../data/'
    width_db_file = width_db_dir+'nAmerica_GRWDL.h5'
    
    # This is the SWOT data
    
    #data_dir = '/Volumes/Reservoir/Data/SWOT/'
    data_dir = '../../data/examples/'
    l2_file = (data_dir + 'simulated_swot_test_data.nc')
    
    # This is the file for the reach data base
    
    db_dir = '../../data/nAmerica_GRWDL_river_topo/'
    shape_file_root = db_dir+'nAmerica_GRWDL_river_topo'
For this example, only data with classification=1 labels are used (pure
water pixels) and the latitude and longitude are assumed to be known
from the reference interferogram.

.. code:: python

    class_list=[1]
    lat_kwd='no_layover_latitude'
    lon_kwd='no_layover_longitude'
    
    l2 = SWOTL2(l2_file,class_list=class_list,lat_kwd=lat_kwd,lon_kwd=lon_kwd)

.. parsed-literal::

    Dataset opened
    Good data selected & bounding box calculated.
    lat/lon read
    projection set and x,y calculated


The heights to be used are the noisy height and the water\_height, the
true height weighted by the radar response.

.. code:: python

    h = l2.get('height')
    htrue = l2.get('water_height')
    
Extract the river reaches overlapping the data bounding box
----------------------------

In this case, the reaches are from the North America Global River Width
Database for Large Rivers. The data were made by `Tamlin
Pavelsky <mailto:pavelsky@unc.edu>`__ and `George
Allen <mailto:georgehenryallen@gmail.com>`__ and modified by `Ernesto
Rodriguez <mailto:ernesto.rodriguez@jpl.nasa.gov>`__ so that reaches
appeared as topologically connected. The results are in the shapefile
*nAmerica\_GRWDL\_river\_topo.shp*.

The next step scans through the database and finds all of the
overlapping reaches in the data base (only one is found). A small buffer
of about 2km is placed about the data bounding box to make sure the
river is not cut.

.. code:: python

    clip_buffer = 0.02
    reaches = ReachExtractor(shape_file_root, l2,clip_buffer=clip_buffer)
The following step prints out the reach index and metadata:

.. code:: python

    print 'Reach indexes:',reaches.reach_idx
    print 'Metadata:'
    reaches[0].metadata

.. parsed-literal::

    Reach indexes: [704]
    Metadata:




.. parsed-literal::

    {'break_idx': 1388940,
     'latmax': 40.03200149536133,
     'latmin': 38.16350173950195,
     'lonmax': -121.5009994506836,
     'lonmin': -122.11900329589844,
     'npoints': 10144,
     'reach': 375745.03125,
     'reach_idx': 704,
     'width_max': 532.0,
     'width_mean': 115.74655151367188,
     'width_min': 30.0,
     'width_std': 52.802425384521484}



Extract the width information for this reach from the width database
----------------------------
.. code:: python

    width_db = WidthDataBase(width_db_file)
    
    max_width = width_db.get_river(reaches.reach_idx[0],
                                          columns=['width'],
                                 asarray=True,transpose=False,
                                 bounding_box=l2.bounding_box,
                                 clip_buffer=clip_buffer)
.. code:: python

    print 'max_width length:',len(max_width)
    print 'x length:',len(reaches[0].x)

.. parsed-literal::

    max_width length: 1882
    x length: 1882


Look at where the centerline lies with respect to the data
----------------------------

From this figure, it is clear that the centerline and the SRTM water
mask align well in many places, but can have disagreements that are
significant compared to the width.

.. code:: python

    figsize(8,8)
    plot(l2.x/1.e3,l2.y/1.e3,'k.',alpha=0.05)
    scatter(reaches[0].x/1.e3,reaches[0].y/1.e3,
            c=max_width,
            edgecolor='none',alpha=0.4,
            vmin=0,vmax=300)
    xlim(-5,5 )
    ylim(-10,10)
    colorbar(label='Width');


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_18_0.png


**Caption**: In the observations above, the measurement locations are
black pixels, while the centerline is in color, width the color
representing the estimated width. Notice that the river meanders
significantly (i.e., by more than the estimated width) away from the
centerline. Notice also that the estimated width does not seem to agree
well with the data point distribution.

Refining the centerline and associating a width to it
----------------------------

Use IteratedRiverObs to refine the centerline and resample the width to
the new centerline. Note that the parameter ``scalar_max_with`` is set
large, since a significant number of points are far away from the
guessed centerline.

.. code:: python

    # First step, initialize observations
    
    scalar_max_width = 600.
    
    ds = 50.
    minobs = 10
    river_obs = IteratedRiverObs(reaches[0],l2.x,l2.y,
                             ds=ds,max_width=scalar_max_width,minobs=minobs) 

.. parsed-literal::

    Centerline initialized
    Local coordiantes calculated


Plot the distribution of distances to the input centerline:

.. code:: python

    figsize(10,5)
    subplot(1,2,1)
    hist(river_obs.d,bins=100,log=False)
    xlabel('Distance to node (m)')
    ylabel('N observations')
    grid();
    subplot(1,2,2)
    hist(river_obs.n,bins=100,log=False)
    xlabel('Normal coordinate (m)')
    ylabel('N observations')
    grid()
    tight_layout();


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_23_0.png


**Caption**: Notice that the distribution is skewed and contains points
that outside the nominal river width.

The following step, refines the centerline, but a max\_width vector is
not yet associated with it.

.. code:: python

    weights = True
    smooth = 1.e-2
    river_obs.iterate(weights=weights,smooth=smooth)

.. parsed-literal::

    iteration 0 maximum coordinate change: 279.979018


Now look at the new centerline (no widths associated).

.. code:: python

    # retrieve the centerline coordinates
    
    xc, yc = river_obs.get_centerline_xy()
    
    figsize(8,8)
    plot(l2.x/1.e3,l2.y/1.e3,'k.',alpha=0.05)
    scatter(xc/1.e3,yc/1.e3,c='b',
            #c=max_width,
            edgecolor='none',alpha=0.4)#,
            #vmin=0,vmax=300)
    xlim(-5,5 )
    ylim(-10,10);


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_28_0.png


Now look at the new distribution of points relative to the new
centerline:

.. code:: python

    figsize(10,5)
    subplot(1,2,1)
    hist(river_obs.d,bins=100,log=False)
    xlabel('Distance to node (m)')
    ylabel('N observations')
    grid();
    subplot(1,2,2)
    hist(river_obs.n,bins=100,log=False)
    xlabel('Normal coordinate (m)')
    ylabel('N observations')
    grid()
    tight_layout();


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_30_0.png


**Caption**: Much nicer! The points are close to the new river
centerline.

In the next step, associate the estimated widths with the new
centerline. Notice that this can only be done approximately, since there
is no unique one to one mapping.

.. code:: python

    # These are the old centerline coordinates
    
    xw = reaches[0].x
    yw = reaches[0].y
    
    # This step makes the association
    
    river_obs.add_centerline_obs(xw,yw,max_width,'max_width')
Now get the centerline coordinates that could be associated with a
max\_width and the associated max\_width and plot the results.

.. code:: python

    xi, yi, wi = river_obs.get_centerline_xyv('max_width')
    
    figsize(8,8)
    plot(l2.x/1.e3,l2.y/1.e3,'k.',alpha=0.05)
    scatter(xi/1.e3,yi/1.e3,
            c=wi,
            edgecolor='none',alpha=0.4,
            vmin=0,vmax=300)
    xlim(-5,5 )
    ylim(-10,10);
    colorbar(label='Width');


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_35_0.png


**Caption**: Not too bad, but notice the widths do not always match the
data distribution.

Compute the center line and associate observations with it
----------------------------

In the following step, the new centerline and max\_width vector are used
to exclude the bad data.

.. code:: python

    river_obs.reinitialize()

.. parsed-literal::

    636 636 636
    Centerline initialized
    Local coordiantes calculated


Plot the maximum width along the channel:

.. code:: python

    figsize(5,5)
    scatter(river_obs.centerline.x/1.e3,river_obs.centerline.y/1.e3,
            c=river_obs.max_width,
            edgecolor='none',alpha=0.4,
            vmin=0,vmax=300)
    colorbar(label='Width');


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_40_0.png


In the following step, the height and true heights are associated with
the nodes:

.. code:: python

    river_obs.add_obs('htrue',htrue)
    river_obs.add_obs('h',h)
    river_obs.load_nodes(['h','htrue'])
Get some rough statistics of the differences in heights for the data
associated with each node:

.. code:: python

    hn_mean = N.array(river_obs.get_node_stat('mean','h'))
    hn_median = N.array(river_obs.get_node_stat('median','h'))
    hstdn = N.array(river_obs.get_node_stat('stderr','h'))
    htn = N.array(river_obs.get_node_stat('mean','htrue'))
    sn = N.array(river_obs.get_node_stat('mean','s'))
.. code:: python

    ave = N.mean(hn_mean - htn)*100
    err = N.std(hn_mean - htn)*100
    print 'Mean statitics:   average: %.1f cm std: %.1f cm'%(ave,err)
    
    ave = N.mean(hn_median - htn)*100
    err = N.std(hn_median - htn)*100
    print 'Median statitics: average: %.1f cm std: %.1f cm'%(ave,err)

.. parsed-literal::

    Mean statitics:   average: 4.8 cm std: 66.8 cm
    Median statitics: average: 7.5 cm std: 71.6 cm


Plot the distances from the points to the nodes, and the cross-river
normal coordinate:

.. code:: python

    figsize(10,5)
    subplot(1,2,1)
    hist(river_obs.d,bins=100,log=False)
    xlabel('Distance to node (m)')
    ylabel('N observations')
    grid();
    subplot(1,2,2)
    hist(river_obs.n,bins=100,log=False)
    xlabel('Normal coordinate (m)')
    ylabel('N observations')
    grid()
    tight_layout();


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_47_0.png


Plot the river geometry, number of observations per node, height vs true
height, and a plot of height vs the reach distance downriver.

.. code:: python

    figsize(10,10)
    subplot(2,2,1)
    idx = river_obs.populated_nodes
    plot(river_obs.centerline.x[idx]/1.e3,river_obs.centerline.y[idx]/1.e3,
         '.',alpha=0.5)
    xlabel('x (km)')
    ylabel('y (km)')
    subplot(2,2,2)
    idx = river_obs.populated_nodes
    plot(river_obs.centerline.s[idx]/1.e3,river_obs.nobs[idx],'o',alpha=0.5)
    grid()
    xlabel('Reach (km)')
    ylabel('Number of observations')
    subplot(2,2,3)
    plot(htrue,h,'.',alpha=0.05)
    plot([-10,60],[-10,60],'--k',alpha=0.5)
    xlim(-10,60)
    ylim(-10,60)
    grid()
    xlabel('h true (m)')
    ylabel('h measured (m)')
    subplot(2,2,4)
    plot(river_obs.s/1.e3,river_obs.htrue,'.',alpha=0.1)
    plot(river_obs.s/1.e3,river_obs.h,'.',alpha=0.05)
    grid()
    xlabel('Reach (km)')
    ylabel('Height (m)')
    tight_layout();


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_49_0.png


Plot the true river height averaged over each node, and the measured
height, similarly averaged for mean and median averaging. The gray lines
are 2 standard errors away from the truth.

.. code:: python

    figsize(10,5)
    subplot(1,2,1)
    plot(sn/1.e3,htn,'.',alpha=0.1)
    plot(sn/1.e3,hn_mean,'.',alpha=0.2)
    plot(sn/1.e3,htn+2*hstdn,'-k',alpha=0.1)
    plot(sn/1.e3,htn-2*hstdn,'-k',alpha=0.1)
    xlabel('Reach (km)')
    ylabel('Height (m)')
    title('Mean')
    
    subplot(1,2,2)
    plot(sn/1.e3,htn,'.',alpha=0.1)
    plot(sn/1.e3,hn_median,'.',alpha=0.2)
    plot(sn/1.e3,htn+2*hstdn,'-k',alpha=0.1)
    plot(sn/1.e3,htn-2*hstdn,'-k',alpha=0.1)
    xlabel('Reach (km)')
    title('Median')
    ylabel('Height (m)')
    tight_layout();


.. image:: CenterlineRefinementExample_files/CenterlineRefinementExample_51_0.png

