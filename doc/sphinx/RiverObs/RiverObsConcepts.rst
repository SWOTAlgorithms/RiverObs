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
 
The Centerline object can be thought of as a curved one-dimensional
coordinate line, with a set of river node locations defined along it, and  with
the capability to provide a mapping between any point in the plane (as
long as it is not too far from the Centerline) to one of its node locations,
and assigning normal and tangential coordinates relative to the node
location. The basic functionality of the Centerline is reviewed in
:ref:`centerline-usage`. Refining the Centerline so that it follows
the data more closely is reviewed in :ref:`centerline-refinement`. 

.. _river-node-overview:

RiverNode
---------

A RiverNode is a data container associated with points on the
Centerline. At a minimum, a RiverNode has the following elements:

    index : int
        index in the center line corresponding to this node
    d : array_like
        distance from the node to each of the data points
    x : array_like
        x coordinate of each measurement associated with the node
    y : array_like
        y coordinate of each measurement associated with the node
    s : array_like
        along-track coordinate (relative to the node center) for each point 
    n : array_like
        across-track (normal) coordinate (relative to the node center) for each point
    ds : float
        along-track dimension for this node. Defaults to 1. Needs to be set
        correctly for width_area to work.

In addition to this basic data, any other object can be stored in a
RiverNode. Once data is stored in a node, it can be queried to produce
a node statistic; e.g., the mean and standard deviations of the data
stored in the node. Several statistic functions are provided in the
RiverNode API.

.. _river-obs-overview:

RiverObs
--------

A RiverObs is an object which contains a Centerline and a set of
RiverNodes associated with that centerline. In addition, it stores the
observation data and can provide statistic lists for each node. A
derived class, IteratedRiverObs, also has the capability to iterate
the centerline to fit the data better. An example of using an
IteratedRiberObs to refine the centerline and load data onto all the
nodes is provided in :ref:`centerline-refinement`.
