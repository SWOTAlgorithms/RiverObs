'''
Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import os
import numpy as np
from collections import OrderedDict as odict

from SWOTWater.products.product import Product, FILL_VALUES, textjoin
from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9

class L2PIXCVector(Product):
    UID = "l2_hr_pixcvector"

    # copied from L2HRPIXC
    ATTRIBUTES = odict([
        ['Conventions', {}],
        ['title', {}],
        ['institution', {}],
        ['source', {}],
        ['history', {}],
        ['mission_name', {}],
        ['references', {}],
        ['reference_document', {}],
        ['contact', {}],
        ['cycle_number', {}],
        ['pass_number', {}],
        ['tile_number', {}],
        ['swath_side', {}],
        ['tile_name', {}],
        ['wavelength', {}],
        ['near_range', {}],
        ['nominal_slant_range_spacing', {}],
        ['start_time', {}],
        ['stop_time', {}],
        ['ephemeris', {}],
        ['yaw_flip', {}],
        ['hpa_cold', {}],
        ['processing_beamwidth', {}],
        ['inner_first_latitude', {}],
        ['inner_first_longitude', {}],
        ['inner_last_latitude', {}],
        ['inner_last_longitude', {}],
        ['outer_first_latitude', {}],
        ['outer_first_longitude', {}],
        ['outer_last_latitude', {}],
        ['outer_last_longitude', {}],
        ['slc_first_line_index_in_tvp', {}],
        ['slc_last_line_index_in_tvp', {}],
        ['xref_input_l1b_hr_slc_file', {}],
        ['xref_input_static_karin_cal_file', {}],
        ['xref_input_ref_dem_file', {}],
        ['xref_input_water_mask_file', {}],
        ['xref_input_static_geophys_files', {}],
        ['xref_input_dynamic_geophys_files', {}],
        ['xref_input_int_lr_xover_cal_file', {}],
        ['xref_l2_hr_pixc_config_parameters_file', {}],
        ['ellipsoid_semi_major_axis', {}],
        ['ellipsoid_flattening', {}],
        ['interferogram_size_azimuth', {}],
        ['interferogram_size_range', {}],
        ['looks_to_efflooks', {}]])

    DIMENSIONS = odict([['points', 0]])
    VARIABLES = odict([
        ['azimuth_index',
         odict([['dtype', 'i4'],
                ['long_name', 'rare interferogram azimuth index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 999999],
                ['comment', 'rare interferogram azimuth index'],
                ])],
        ['range_index',
         odict([['dtype', 'i4'],
                ['long_name', 'rare interferogram range index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 999999],
                ['comment', 'rare interferogram range index'],
                ])],
        ['latitude_vectorproc',
         odict([['dtype', 'f8'],
                ['long_name', 'latitude'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -90],
                ['valid_max', 90],
                ['comment', 'geodetic latitude (degrees north of equator)'],
                ])],
        ['longitude_vectorproc',
         odict([['dtype', 'f8'],
                ['long_name', 'longitude'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', -180],
                ['valid_max', 180],
                ['comment', 'longitude (east of the prime meridian)'],
                ])],
        ['height_vectorproc',
         odict([['dtype', 'f4'],
                ['long_name', 'height above reference ellipsoid'],
                ['units', 'm'],
                ['valid_min', -999999],
                ['valid_max', 999999],
                ['comment', 'height above reference ellipsoid'],
                ])],
        ['node_index',
         odict([['dtype', 'i4'],
                ['long_name', 'node index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', 'index of node this pixel was assigned to'],
                ])],
        ['reach_index',
         odict([['dtype', 'i4'],
                ['long_name', 'reach index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', 'index of reach this pixel was assigned to'],
                ])],
        ['segmentation_label',
         odict([['dtype', 'i4'],
                ['long_name', 'segmentation label'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', 'segmentation label'],
                ])],
        ['good_height_flag',
         odict([['dtype', 'u1'],
                ['long_name', 'good height flag'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['comment', 'good height flag'],
                ])],
        ['distance_to_node',
         odict([['dtype', 'f4'],
                ['long_name', 'distance to node'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 9999],
                ['comment', 'distance to node'],
                ])],
        ['along_reach',
         odict([['dtype', 'f4'],
                ['long_name', 'along reach distance'],
                ['units', 'm'],
                ['valid_min', -999999],
                ['valid_max', 999999],
                ['comment', 'along reach distance'],
                ])],
        ['cross_reach',
         odict([['dtype', 'f4'],
                ['long_name', 'across reach distance'],
                ['units', 'm'],
                ['valid_min', -999999],
                ['valid_max', 999999],
                ['comment', 'across reach distance'],
                ])],
        ['pixc_index',
         odict([['dtype', 'i4'],
                ['long_name', 'pixel cloud index'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', 'index in pixel cloud product'],
                ])],
        ['lake_flag',
         odict([['dtype', 'u1'],
                ['long_name', 'lake flag'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['comment', 'lake flag'],
                ])],
        ['ice_clim_f',
         odict([['dtype', 'i2'],
                ['long_name', 'climatological ice cover flag'],
                ['source', 'UNC'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover
                    not_available""")],
                ['flag_values', np.array([0, 1, 2, 255]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 255],
                ['comment', textjoin("""
                    Climatological ice cover flag indicating whether the node
                    is ice-covered on the day of the observation based on
                    external climatological information (not the SWOT
                    measurement).  Values of 0, 1, and 2 indicate that the
                    node is not ice covered, partially ice covered, and fully
                    ice covered, respectively. A value of 255 indicates that
                    this flag is not available.""")],
                ])],
        ['ice_dyn_f',
         odict([['dtype', 'i2'],
                ['long_name', 'dynamical ice cover flag'],
                ['source', 'UNC'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover
                    not_available""")],
                ['flag_values', np.array([0, 1, 2, 255]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 255],
                ['comment', textjoin("""
                    Dynamic ice cover flag indicating whether the surface is
                    ice-covered on the day of the observation based on
                    analysis of external satellite optical data.  Values of
                    0, 1, and 2 indicate that the node is not ice covered,
                    partially ice covered, and fully ice covered, respectively.
                    A value of 255 indicates that this flag is not available.
                    """)],
                ])],
        ])

    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    def update_from_rivertile(self, rivertile):
        """Updates some stuff in PIXCVecRiver from RiverTile"""
        for node_id, ice_clim_f, ice_dyn_f in zip(
            rivertile.nodes.node_id, rivertile.nodes.ice_clim_f,
            rivertile.nodes.ice_dyn_f):
            mask = self.node_index == node_id
            self.ice_clim_f[mask] = ice_clim_f
            self.ice_dyn_f[mask] = ice_dyn_f
