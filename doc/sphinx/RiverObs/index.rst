
RiverObs's documentation
====================================

The intent of the RiverObs package is to provide a set of Python
classes that enable associating data, remote sensing or *in situ*,
with rivers, as a whole, and with specific river locations. Although
the package was developed to support the generation of hydrology
observables from the NASA `SWOT <http://swot.jpl.nasa.gov>`__ mission,
it is quite general and can be used for muliple data sets including
point gauge data, irregularly distributed data, such as lidar point
clouds, or line data, such as GPS surveys.

Before taking a look at the detailed :ref:`API`, the novice user
should get familiar with the basic :ref:`Concepts`. The
:ref:`Overview` presents a brief description of the packages that form
part of the distribution, and their dependencies. The :ref:`Installation`
give detailed instructions on how to set up and build the
package. The :ref:`API` contains a detailed listing of all of
the packages, and interface documentation.

Tutorial examples of how to use some of the basic classes using simulated data
are contained in :ref:`centerline-usage` and
:ref:`centerline-refinement`.

Finally, a typical workflow will consist in making a set of input
reaches, reading these reaches and some data, estimating river
widths, heights and slopes, and, finally, writing out the results in
files that can be read by GIS programs. These steps are illustrated in
:ref:`reach-preprocessing` and :ref:`end-to-end-example`. 


Contents:

.. toctree::
   :maxdepth: 1

   RiverObsConcepts
   Overview 
   Installation 
   CenterlineUsageExample
   CenterlineRefinementExample
   ReachPreProcessorExample 
   EndToEndProcessingExample.rst
   API


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

