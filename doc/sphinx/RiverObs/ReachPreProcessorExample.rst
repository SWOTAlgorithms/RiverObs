.. _reach-preprocessing:

ReachPreProcessor Example
=========================

In this example, we demonstrate how to load a set of predefined reaches
and process them so that they are better suited to the data and
applications.

In this case, the reaches are from the North America Global River Width
Database for Large Rivers. The data were made by `Tamlin
Pavelsky <mailto:pavelsky@unc.edu>`__ and `George
Allen <mailto:georgehenryallen@gmail.com>`__ and modified by `Ernesto
Rodriguez <mailto:ernesto.rodriguez@jpl.nasa.gov>`__ so that reaches
appeared as topologically connected. The results are in the shapefile
``nAmerica_GRWDL_river_topo.shp``. A rapidly searchable width database,
segemented by reach has been made into a pytables hdf5 file and will be
used as width data base.

As sample data for this data collection, we use some simulated SWOT data
over the Sacrramento River.

.. code:: python

    # This is the file for the width data base
    
    width_db_dir = '../../data/'
    width_db_file = width_db_dir+'nAmerica_GRWDL.h5'
    
    # This is the SWOT data
    
    data_dir = '../../data/examples/'
    l2_file = (data_dir + 'simulated_swot_test_data.nc')
    
    # This is the file for the reach data base
    
    db_dir = '../../data/nAmerica_GRWDL_river_topo/'
    shape_file_root = db_dir+'nAmerica_GRWDL_river_topo'
Import the Required Modules and Initialize Objects
--------------------------------------------------

.. code:: python

    from SWOTRiver import SWOTL2
    from RiverObs import ReachPreProcessor, RiverReachWriter
Read the SWOT simulated data to find the data bounding box and define
the projection function. For the editing function, use the true data
classification and true water location.

.. code:: python

    class_list=[1]
    lat_kwd='no_layover_latitude'
    lon_kwd='no_layover_longitude'
    
    l2 = SWOTL2(l2_file,class_list=class_list,lat_kwd=lat_kwd,lon_kwd=lon_kwd)


Initialize the pre-processor with the GRWDL data. The next step scans
through the database and finds all of the overlapping reaches in the
data base (only one is found). A small buffer of about 2km is placed
about the data bounding box to make sure the river is not cut.

.. code:: python

    clip_buffer = 0.02
    reaches = ReachPreProcessor(shape_file_root, l2,clip_buffer=clip_buffer,width_db_file=width_db_file)
.. code:: python

    print 'Number of reaches found:',reaches.nreaches
    print 'Reach indexes:',reaches.reach_idx
    for i,reach in enumerate(reaches):
        print('Reach %d Metadata'%i)
        print(reach.metadata)

.. parsed-literal::

    Number of reaches found: 1
    Reach indexes: [704]
    Reach 0 Metadata
    {'latmax': 40.03200149536133, 'width_std': 52.802425384521484, 'width_min': 30.0, 'reach': 375745.03125, 'width_mean': 115.74655151367188, 'break_idx': 1388940, 'reach_idx': 704, 'lonmin': -122.11900329589844, 'width_max': 532.0, 'npoints': 10144, 'lonmax': -121.5009994506836, 'latmin': 38.16350173950195}


As a guide to editing the reaches, plot the data and reach locations for
each reach found.

.. code:: python

    figsize(6,6)
    plot(l2.lon[::10],l2.lat[::10],'o',alpha=0.1,color='aqua')
    for reach in reaches:
        plot(reach.lon,reach.lat,'.',alpha=0.4,label='Reach %d'%i)
    legend(loc='best')
    title('Reaches vs Data')
    if interactive:
        plugins.connect(gcf(),plugins.MousePosition(fmt='.3f'))


.. image:: ReachPreProcessorExample_files/ReachPreProcessorExample_12_0.png


**Figure 1**: Clearly the reach is much longer than the bit of data
imaged. We can cut it down in several way, as demonstrated in the
following.

Break reach by predefined points
--------------------------------

Use interactive graphics in the iPython notebook to pick the location of
the breaks.

.. code:: python

    start_lons = [-122.010,-121.965,-121.962,-121.978]
    start_lats = [39.760,39.735,39.685,39.648]
    end_lons = [-121.965,-121.962,-121.978,-121.998]
    end_lats = [39.735,39.685,39.648,39.601]
    
    reach_start_list = zip(start_lons,start_lats)
    reach_end_list = zip(end_lons,end_lats)
.. code:: python

    figsize(6,6)
    plot(l2.lon[::10],l2.lat[::10],'o',alpha=0.1,color='aqua')
    for reach in reaches:
        plot(reach.lon,reach.lat,'.',alpha=0.4,label='Reach %d'%i)
    scatter(start_lons,start_lats,s=300,c='r',marker='+')
    scatter(end_lons,end_lats,s=300,c='g',marker='x')
    legend(loc='best')
    title('Reaches vs Data')
    if interactive:
        plugins.connect(gcf(),plugins.MousePosition(fmt='.3f'))


.. image:: ReachPreProcessorExample_files/ReachPreProcessorExample_16_0.png


.. code:: python

    edited_reaches = reaches.split_by_coordinates(reach_start_list,reach_end_list)
.. code:: python

    for reach in edited_reaches:
        print(reach.metadata)

.. parsed-literal::

    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 26871.735741185847), ('parent_end_s', 33511.174223440466), ('indexstart', 609), ('indexend', 750), ('npoints', 142), ('reach_length', 6639.4384822546199), ('width_mean', 152.21830985915494), ('width_max', 258), ('width_min', 60)])
    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 33511.174223440466), ('parent_end_s', 42970.82455543863), ('indexstart', 750), ('indexend', 985), ('npoints', 236), ('reach_length', 9459.6503319981639), ('width_mean', 133.30084745762713), ('width_max', 379), ('width_min', 42)])
    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 42970.82455543863), ('parent_end_s', 55972.882740801113), ('indexstart', 985), ('indexend', 1280), ('npoints', 296), ('reach_length', 13002.058185362483), ('width_mean', 142.07432432432432), ('width_max', 381), ('width_min', 60)])
    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 55972.882740801113), ('parent_end_s', 64900.758838929585), ('indexstart', 1280), ('indexend', 1496), ('npoints', 217), ('reach_length', 8927.8760981284722), ('width_mean', 136.76036866359448), ('width_max', 270), ('width_min', 60)])


.. code:: python

    figsize(6,6)
    for i,reach in enumerate(edited_reaches):
        plot(reach.lon,reach.lat,'.',alpha=0.4,label='Reach %d'%i)
    scatter(start_lons,start_lats,s=300,c='r',marker='+')
    scatter(end_lons,end_lats,s=300,c='g',marker='x')
    legend(loc='best')
    title('Edited Reaches')
    if interactive:
        plugins.connect(gcf(),plugins.MousePosition(fmt='.3f'));


.. image:: ReachPreProcessorExample_files/ReachPreProcessorExample_19_0.png


Break reach by constant reach length
------------------------------------

One possible way to break the reaches, which is not particularly
sensible from a hydrologic perspective, is to break the reaches into
sections of equal lengtgh. This process is demonstarted below:

.. code:: python

    start_s = 25.e3
    ds = 10.e3
    end_s = start_s + 4*ds
    edited_reaches = reaches.split_by_reach_length(ds,start_s=start_s,end_s=end_s)

.. code:: python

    for reach in edited_reaches:
        print(reach.metadata)

.. parsed-literal::

    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 25057.046393810597), ('parent_end_s', 34980.136587280504), ('indexstart', 568), ('indexend', 784), ('npoints', 217), ('reach_length', 9923.0901934699068), ('width_mean', 143.56682027649771), ('width_max', 258), ('width_min', 60)])
    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 35013.596509739444), ('parent_end_s', 44979.249180419421), ('indexstart', 785), ('indexend', 1036), ('npoints', 252), ('reach_length', 9965.6526706799777), ('width_mean', 140.88095238095238), ('width_max', 381), ('width_min', 42)])
    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 45065.000473226915), ('parent_end_s', 54990.715219114747), ('indexstart', 1037), ('indexend', 1259), ('npoints', 223), ('reach_length', 9925.7147458878317), ('width_mean', 141.30941704035874), ('width_max', 296), ('width_min', 60)])
    OrderedDict([('parent_reach_idx', 0), ('parent_start_s', 55013.16272915258), ('parent_end_s', 64956.665436932344), ('indexstart', 1260), ('indexend', 1498), ('npoints', 239), ('reach_length', 9943.502707779764), ('width_mean', 132.68200836820083), ('width_max', 270), ('width_min', 60)])


.. code:: python

    figsize(6,6)
    for i,reach in enumerate(edited_reaches):
        plot(reach.lon,reach.lat,'.',alpha=1,label='Reach %d'%i)
    legend(loc='best')
    title('Edited Reaches')
    if interactive:
        plugins.connect(gcf(),plugins.MousePosition(fmt='.3f'));


.. image:: ReachPreProcessorExample_files/ReachPreProcessorExample_23_0.png


Output the edited reaches and width data base
---------------------------------------------

Once the reaches are edited, the results can be written to a new OGR
supported data format. The default format is "ESRI Shapefile", which is
the format expected by the ReachExtractor class. Also, if desired,
output the width data base in hdf5 pytables format.

First, initialize the RiverReachWriter:

.. code:: python

    node_output_name = 'edited_nodes'
    reach_output_name = 'edited_reaches'
    
    !rm -rf $node_output_name
    !rm -rf $reach_output_name
    
    node_output_variables = ['width']
    reach_output_variables = edited_reaches[0].metadata.keys()
    reach_writer = RiverReachWriter(edited_reaches,node_output_variables,reach_output_variables)
First, write as node shape files as Point data with dbf fields for each
point:

.. code:: python

    reach_writer.write_nodes_ogr(node_output_name)
Write the reach data as LineString shapefiles:

.. code:: python

    reach_writer.write_reaches_ogr(reach_output_name)
Just to show that it can be done, and for visualization, write them also
as KML:

.. code:: python

    reach_writer.write_nodes_ogr(node_output_name+'.kml',driver='KML')
    reach_writer.write_reaches_ogr(reach_output_name+'.kml',driver='KML')
    !ls *.kml

.. parsed-literal::

    edited_nodes.kml   edited_reaches.kml


Finally, one can write a new width data base, if one was provided to
begin width.

.. code:: python

    width_db_file = 'edited_width_db'
    river_df, reach_df = reach_writer.write_width_db(width_db_file,output_format='h5')
    !ls *.h5

.. parsed-literal::

    edited_width_db.h5


.. code:: python

    print reach_df.head()

.. parsed-literal::

       break_idx     latmax     latmin      lonmax      lonmin  npoints  \
    0        216  39.767700  39.730701 -121.955002 -122.023003      217   
    1        468  39.730400  39.673801 -121.942001 -121.966003      252   
    2        691  39.676601  39.647301 -121.967003 -121.999001      223   
    3        930  39.648998  39.600498 -121.969002 -121.997002      239   
    
             reach  width_max  width_mean  width_min  width_std  
    0  9950.324219        258  143.566820         60  42.893464  
    1  9950.324219        381  140.880952         42  61.197262  
    2  9950.324219        296  141.309417         60  49.618619  
    3  9950.324219        270  132.682008         60  38.002413  


.. code:: python

    print river_df.head()

.. parsed-literal::

       width  nchannels  reservoir        long        lat  reach_index       reach
    0    150          1          0 -122.023003  39.767700            0   22.081333
    1    108          1          0 -122.023003  39.767502            0   55.627975
    2     84          1          0 -122.023003  39.767200            0  147.514374
    3     84          1          0 -122.022003  39.766899            0  180.636383
    4    127          1          0 -122.022003  39.766602            0  203.142365


.. code:: python

    print river_df.tail()

.. parsed-literal::

         width  nchannels  reservoir        long        lat  reach_index  \
    926    150          1          0 -121.997002  39.601501            3   
    927    150          1          0 -121.997002  39.601299            3   
    928    192          1          0 -121.997002  39.601002            3   
    929    150          1          0 -121.997002  39.600700            3   
    930    192          1          0 -121.997002  39.600498            3   
    
               reach  
    926  9861.149414  
    927  9894.271484  
    928  9927.818359  
    929  9950.324219  
    930  9950.324219  

