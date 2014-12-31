.. _Concepts:

RiverObs Concepts
=============

The major components of the RiverObs concept are illustrated in the 
:ref:`centerline_nodes` figure, which illustrates the three key concepts
underlying RiverObs: the :ref:`centerline-overview`, the 
:ref:`river-node-overview`, and the :ref:`river-obs-overview`.


.. _centerline_nodes:
.. figure:: images/centerline_nodes.png
   :width: 400 px
   :scale: 100 %
   :align: center

   *Illustration of the concepts underlying RiverObs.* 

In the figure, the physical river is shown in blue, the Centerline is the black line running along the middle. 
The centerline consists of a set of RiverNodes (red dots), and every RiverNode 
has an associated coordinate system (in green), one of whose axes (the *s* or 
*along-track* axis) is tangent to the Centerline, and other axis (the *n* or *normal*
axis) is perpendicular to it and defines a right-handed system. Data (for example,
the orange cross) are associated with the closest node (node regions of influence 
are shown as dashed lines) and assigned node index, *s* and *n* coordinates.
A RiverObs object contains a Centerline, a list of nodes (and their associated data), 
and functions to gather information for all the nodes.


.. _centerline-overview:

Centerline
---------
 
The

.. _river-node-overview:

RiverNode
---------

The

.. _river-obs-overview:

RiverObs
--------

The
