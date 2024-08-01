'''
Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import os
import netCDF4
import warnings
import numpy as np
from collections import OrderedDict as odict

from SWOTWater.products.product import \
    Product, FILL_VALUES, textjoin, ProductTesterMixIn
from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9

from SWOTRiver.products.rivertile import RIVER_PRODUCT_ATTRIBUTES

class L2PIXCVector(ProductTesterMixIn, Product):
    UID = "l2_hr_pixcvector"

    # copied from L2HRPIXC
    ATTRIBUTES = odict([
        ['Conventions', {
            'dtype': 'str', 'docstr': textjoin("""
                NetCDF-4 conventions adopted in this product. This attribute
                should be set to CF-1.7 to indicate that the group is compliant
                with the Climate and Forecast NetCDF conventions."""),
            'value': 'CF-1.7'}],
        ['title', {
            'dtype': 'str', 'docstr': textjoin("""
                Level 2 KaRIn high rate pixel cloud vector river product."""),
            'value': textjoin("""
                Level 2 KaRIn high rate pixel cloud vector river product.""")}],
        ['short_name', {'dtype': 'str',
            'docstr': 'L2_HR_PIXCVecRiver', 'value': 'L2_HR_PIXCVecRiver'}],
        ['institution', {
            'dtype': 'str', 'docstr': textjoin("""
                Name of producing agency."""), 'value': 'JPL'}],
        ['source', {
            'dtype': 'str', 'docstr': textjoin("""
                The method of production of the original data. If it was
                model-generated, source should name the model and its version,
                as specifically as could be useful. If it is observational,
                source should characterize it (e.g., 'radiometer')."""),
            'value': 'Ka-band radar interferometer'}],
        ['history', {'dtype': 'str',
            'docstr': textjoin("""
                UTC time when file generated. Format is:
                'YYYY-MM-DD hh:mm:ss : Creation' """)}],
        ['platform', {'dtype': 'str' ,'value':'SWOT', 'docstr': 'SWOT'}],
        ['references', {'dtype': 'str',
            'docstr': textjoin("""
                Published or web-based references that describe
                the data or methods used to product it. Provides version number
                of software generating product.""")}],
        ['reference_document', {'dtype': 'str',
            'docstr': textjoin("""
                Name and version of Product Description Document
                to use as reference for product.""")}],
        ['product_version', {'dtype': 'str',
            'docstr': 'Version identifier of this data file'}],
        ['crid', {'dtype': 'str',
            'docstr': textjoin("""
                Composite release identifier (CRID) of the data system used to
                generate this file""")}],
        ['pge_name', {'dtype': 'str',
            'docstr': textjoin("""
                Name of the product generation executable (PGE) that created
                this file""")}],
        ['pge_version', {'dtype': 'str',
            'docstr': textjoin("""
                Version identifier of the product generation executable (PGE)
                that created this file""")}],
        ['contact', {'dtype': 'str',
            'docstr': textjoin("""
                Contact information for producer of product.
                (e.g., 'ops@jpl.nasa.gov').""")}],
        ['cycle_number', {'dtype': 'i2',
            'docstr': 'Cycle number of the product granule.'}],
        ['pass_number', {'dtype': 'i2',
            'docstr': 'Pass number of the product granule.'}],
        ['tile_number', {'dtype': 'i2',
            'docstr': 'Tile number in the pass of the product granule.'}],
        ['swath_side', {'dtype': 'str',
            'docstr': textjoin(
                """'L' or 'R' to indicate left and right swath,respectively.
                """)}],
        ['tile_name', {'dtype': 'str',
            'docstr': textjoin("""
                Tile name using format PPP_TTTS, where PPP is a 3 digit
                pass number with leading zeros, TTT is a 3 digit tile number
                within the pass, and S is a character 'L' or 'R' for the left
                and right swath, respectively.""")}],
        ['continent_id', {'dtype': 'str',
            'docstr': 'Two-letter continent identifier of the product granule.'}],
        ['continent_code', {'dtype': 'str',
            'docstr': 'One-digit (C) continent code of the product granule.'}],
        ['time_granule_start', {'dtype': 'str',
            'docstr': textjoin("""
                Nominal starting UTC time of product granule.
                Format is: YYYY-MM-DDThh:mm:ss.ssssssZ""")}],
        ['time_granule_end', {'dtype': 'str',
            'docstr': textjoin("""
                Nominal ending UTC time of product granule.
                Format is: YYYY-MM-DDThh:mm:ss.ssssssZ""")}],
        ['time_coverage_start', {'dtype': 'str',
            'docstr': textjoin("""
                UTC time of first measurement.
                Format is: YYYY-MM-DD hh:mm:ss.ssssssZ""")}],
        ['time_coverage_end', {'dtype': 'str',
            'docstr': textjoin("""
                UTC time of last measurement.
                Format is: YYYY-MM-DD hh:mm:ss.ssssssZ""")}],
        ['geospatial_lon_min',  {'dtype': 'f8',
            'docstr': "Westernmost longitude (deg) of granule bounding box"}],
        ['geospatial_lon_max',  {'dtype': 'f8',
            'docstr': "Easternmost longitude (deg) of granule bounding box"}],
        ['geospatial_lat_min',  {'dtype': 'f8',
            'docstr': "Southernmost latitude (deg) of granule bounding box"}],
        ['geospatial_lat_max',  {'dtype': 'f8',
            'docstr': "Northernmost latitude (deg) of granule bounding box"}],
        ['inner_first_longitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner longitude (degrees_east) for the first
                 range line and inner part of the swath""")}],
        ['inner_first_latitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner latitude (degrees_north) for the first
                 range line and inner part of the swath""")}],
        ['inner_last_longitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner longitude (degrees_east) for the last
                 range line and inner part of the swath""")}],
        ['inner_last_latitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner latitude (degrees_north)  for the last
                 range line and inner part of the swath""")}],
        ['outer_first_longitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner longitude (degrees_east) for the first
                 range line and outer part of the swath""")}],
        ['outer_first_latitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner latitude (degrees_north) for the first
                 range line and outer part of the swath""")}],
        ['outer_last_longitude',{'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner longitude (degrees_east) for the last
                 range line and outer part of the swath""")}],
        ['outer_last_latitude', {'dtype': 'f8',
            'docstr': textjoin("""
                 Nominal swath corner latitude (degrees_north) for the last
                 range line and outer part of the swath""")}],
        ['near_range', {'dtype': 'f8',
            'docstr': 'The slant range (m) for the first image pixel.'}],
        ['nominal_slant_range_spacing', {'dtype': 'f8',
            'docstr': textjoin("""
                The range spacing (m) corresponding to the 200 MHz
                sampling frequency""")}],
        ['xref_l2_hr_pixc_files',
         RIVER_PRODUCT_ATTRIBUTES['xref_l2_hr_pixc_files']],
        ['xref_param_l2_hr_rivertile_files',
         RIVER_PRODUCT_ATTRIBUTES['xref_param_l2_hr_rivertile_files']],
        ['xref_prior_river_db_files',
         RIVER_PRODUCT_ATTRIBUTES['xref_prior_river_db_files']],
        ['xref_reforbittrack_files',
         RIVER_PRODUCT_ATTRIBUTES['xref_reforbittrack_files']],
        ['ellipsoid_semi_major_axis', {'dtype': 'f8',
            'docstr': 'Semi-major axis of reference ellipsoid in meters.'}],
        ['ellipsoid_flattening', {'dtype': 'f8',
            'docstr': 'Flattening of reference ellipsoid'}],
        ])


    DIMENSIONS = odict([['points', 0]])
    VARIABLES = odict([
        ['azimuth_index',
         odict([['dtype', 'i4'],
                ['long_name', 'rare interferogram azimuth index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 999999],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment',
                 'Rare interferogram azimuth index (indexed from 0).'],
                ])],
        ['range_index',
         odict([['dtype', 'i4'],
                ['long_name', 'rare interferogram range index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 999999],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment',
                 'Rare interferogram range index (indexed from 0).'],
                ])],
        ['latitude_vectorproc',
         odict([['dtype', 'f8'],
                ['long_name', 'height-constrained geolocation latitude'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -80],
                ['valid_max', 80],
                ['comment', textjoin("""
                    Height-constrained geodetic latitude of the pixel.
                    Units are in degrees north of the equator.""")],
                ])],
        ['longitude_vectorproc',
         odict([['dtype', 'f8'],
                ['long_name', 'height-constrained geolocation longitude'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', -180],
                ['valid_max', 180],
                ['comment', textjoin("""
                    Height-constrained geodetic longitude of the pixel.
                    Positive=degrees east of the Greenwich meridian.
                    Negative=degrees west of the Greenwich meridian.""")],
                ])],
        ['height_vectorproc',
         odict([['dtype', 'f4'],
                ['long_name', 'height above reference ellipsoid'],
                ['units', 'm'],
                ['valid_min', -1500.0],
                ['valid_max', 15000.0],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    'Height-constrained height of the pixel above the
                    reference ellipsoid.""")],
                ])],
        ['reach_id',
         odict([['dtype', 'i8'],
                ['long_name', 'identifier of the associated prior river reach'],
                ['valid_min', 0],
                ['valid_max', 9223372036854775807],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Unique reach identifier from the prior river database.
                    The format of the identifier is CBBBBBRRRRT, where
                    C=continent, B=basin, R=reach, T=type.""")],
                ])],
        ['node_id',
         odict([['dtype', 'i8'],
                ['long_name',
                 "identifier of the associated prior river node"],
                ['valid_min', 0],
                ['valid_max', 9223372036854775807],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Unique node identifier from the prior river database.
                    The format of the identifier is CBBBBBRRRRNNNT, where
                    C=continent, B=basin, R=reach, N=node, T=type of water
                    body.""")],
                ])],
        ['ice_clim_f',
         odict([['dtype', 'i1'],
                ['long_name', 'climatological ice cover flag'],
                ['standard_name', 'status_flag'],
                ['institution', 'University of North Carolina'],
                ['flag_meanings', textjoin("""
                    no_ice_cover uncertain_ice_cover full_ice_cover""")],
                ['flag_values', np.array([0, 1, 2]).astype('i1')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['source', 'Yang et al. (2020)'],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Climatological ice cover flag indicating whether the
                    pixel is ice-covered on the day of the observation based
                    on external climatological information (not the SWOT
                    measurement). Values of 0, 1, and 2 indicate that the
                    surface is not ice covered, may or may not be partially
                    or fully ice covered, and fully ice covered, respectively.
                    A value of 127 indicates that this flag is not available.
                    """)],
                ])],
        ['ice_dyn_f',
         odict([['dtype', 'i1'],
                ['long_name', 'dynamical ice cover flag'],
                ['standard_name', 'status_flag'],
                ['institution', 'University of North Carolina'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover""")],
                ['flag_values', np.array([0, 1, 2]).astype('i1')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['source', 'Yang et al. (2020)'],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Dynamic ice cover flag indicating whether the pixel is
                    ice-covered on the day of the observation based on
                    analysis of external satellite optical data. Values of
                    0, 1, and 2 indicate that the surface is not ice covered,
                    partially ice covered, and fully ice covered,
                    respectively. A value of 127 indicates that this flag is
                    not available.""")],
                ])],
        ['pixc_index',
         odict([['dtype', 'i4'],
                ['long_name', 'pixel cloud index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Index of the data in the pixel_cloud group of the
                    L2_HR_PIXC file that is associated with the pixel. This
                    index starts counting from zero.""")],
                ])],
        ['lake_flag',
         odict([['dtype', 'u1'],
                ['long_name', 'lake flag'],
                ['flag_meanings', textjoin("""
                    river lake_or_reservoir canal river_under_tide""")],
                ['flag_values', np.array([0, 1, 2, 3]).astype('u1')],
                ['valid_min', 0],
                ['valid_max', 3],
                ['comment', textjoin("""
                    Flag indicating the reach from the PRD.  0= Reach was
                    flagged as belonging to a river. 1= Reach was flagged as
                    belonging to a lake or reservoir. 2= Reach was flagged as
                    being under the influence of a canal. 3= Reach was flagged
                    as a river under the influence of tides.""")],
                ])],
        ['segmentation_label',
         odict([['dtype', 'i4'],
                ['long_name', 'segmentation label'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    A unique number of identifying which connected water
                    segment the pixel was assigned to.""")],
                ])],
        ['distance_to_node',
         odict([['dtype', 'f4'],
                ['long_name', 'distance to node'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 9999],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Distance from the non-improved pixel location to the PRD
                    node that it is associated with.""")],
                ])],
        ['along_reach',
         odict([['dtype', 'f4'],
                ['long_name', 'along reach distance'],
                ['units', 'm'],
                ['valid_min', -999999],
                ['valid_max', 999999],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Along-reach component of non-improved pixel location
                    relative to PRD node location. Negative=nominally upstream
                    of PRD node. Positive=nominally downstream of PRD node""")],
                ])],
        ['cross_reach',
         odict([['dtype', 'f4'],
                ['long_name', 'across reach distance'],
                ['units', 'm'],
                ['valid_min', -999999],
                ['valid_max', 999999],
                ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
                ['comment', textjoin("""
                    Cross-reach component of non-improved pixel location
                    relative to PRD node location. Negative= left side of
                    centerline. Positive= right side of centerline.""")],
                ])],
        ])

    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    @staticmethod
    def dump_xml(pixcvectorriver_xml_file):
        with open(pixcvectorriver_xml_file, 'w') as ofp:
            L2PIXCVector.print_xml(ofp=ofp)

    def update_from_rivertile(self, rivertile):
        """Updates some stuff in PIXCVecRiver from RiverTile"""
        for node_id, ice_clim_f, ice_dyn_f in zip(
            rivertile.nodes.node_id, rivertile.nodes.ice_clim_f,
            rivertile.nodes.ice_dyn_f):
            mask = self.node_id == node_id
            self.ice_clim_f[mask] = ice_clim_f
            self.ice_dyn_f[mask] = ice_dyn_f

    def update_from_pixc(self, pixc_file):
        """Adds some attributes from PIXC file"""

        ATTRS_2COPY_FROM_PIXC = [
            'cycle_number', 'pass_number', 'tile_number', 'swath_side',
            'tile_name', 'inner_first_latitude', 'inner_first_longitude',
            'inner_last_latitude', 'inner_last_longitude',
            'outer_first_latitude', 'outer_first_longitude',
            'outer_last_latitude', 'outer_last_longitude',
            'ellipsoid_semi_major_axis', 'ellipsoid_flattening',
            'near_range', 'nominal_slant_range_spacing',
            'time_granule_start', 'time_granule_end',
            'time_coverage_start', 'time_coverage_end',
            'geospatial_lon_min', 'geospatial_lon_max',
            'geospatial_lat_min', 'geospatial_lat_max']

        with netCDF4.Dataset(pixc_file, 'r') as ifp:
            for attr in ATTRS_2COPY_FROM_PIXC:
                try:
                    value = getattr(ifp, attr)
                except AttributeError:
                    value = getattr(ifp.groups['pixel_cloud'], attr, 'None')
                self[attr] = value

class L2PIXCVectorPlus(L2PIXCVector):
    """L2PIXCVectorPlus is L2PIXCVector with some additional debug outputs"""
    VARIABLES = L2PIXCVector.VARIABLES.copy()
    DIMENSIONS = L2PIXCVector.DIMENSIONS
    VARIABLES['h_flg'] = odict([
        ['dtype', 'i4'],
        ['long_name', 'height flag'],
        ['units', '1'],
        ['valid_min', 0],
        ['valid_max', 1],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'Height flag from RiverObs.'],
        ['dimensions', DIMENSIONS]])
    VARIABLES['area_flg'] = odict([
        ['dtype', 'i4'],
        ['long_name', 'area flag'],
        ['units', '1'],
        ['valid_min', 0],
        ['valid_max', 1],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'Area flag from RiverObs.'],
        ['dimensions', DIMENSIONS]])
    VARIABLES['sig0_flg'] = odict([
        ['dtype', 'i4'],
        ['long_name', 'sigma0 flag'],
        ['units', '1'],
        ['valid_min', 0],
        ['valid_max', 1],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'sigma0 flag from RiverObs.'],
        ['dimensions', DIMENSIONS]])
    VARIABLES['wse_class_flg'] = odict([
        ['dtype', 'i4'],
        ['long_name', 'wse class flag'],
        ['units', '1'],
        ['valid_min', 0],
        ['valid_max', 1],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'wse class flag from RiverObs.'],
        ['dimensions', DIMENSIONS]])
    VARIABLES['geolocation_qual'] = odict([
        ['dtype', 'u4'],
        ['long_name', 'geolocation_qual'],
        ['units', '1'],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'geolocation_qual from PIXC.'],
        ['dimensions', DIMENSIONS]])
    VARIABLES['classification_qual'] = odict([
        ['dtype', 'u4'],
        ['long_name', 'classification_qual'],
        ['units', '1'],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'classification_qual from PIXC.'],
        ['dimensions', DIMENSIONS]])
    VARIABLES['sig0_qual'] = odict([
        ['dtype', 'u4'],
        ['long_name', 'sig0_qual'],
        ['units', '1'],
        ['coordinates', 'longitude_vectorproc latitude_vectorproc'],
        ['comment', 'sig0_qual from PIXC.'],
        ['dimensions', DIMENSIONS]])
