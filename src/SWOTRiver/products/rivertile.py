'''
Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import os
import textwrap
import numpy as np
from collections import OrderedDict as odict

from SWOTRiver.products.product import Product, FILL_VALUES

class L2HRRiverTile(Product):
    UID = "l2_hr_rivertile"
    ATTRIBUTES = []
    GROUPS = odict([
        ['nodes', 'RiverTileNodes'],
        ['reaches', 'RiverTileReaches'],
    ])

def textjoin(text):
    """Dedent join and strip text"""
    text = textwrap.dedent(text)
    text = text.replace('\n', ' ')
    text = text.strip()
    return text

class RiverTileNodes(Product):
    ATTRIBUTES = []
    DIMENSIONS = odict([['nodes', 0]])
    VARIABLES = odict([
        ['reach_id',
         odict([['dtype', 'i4'],
                ['long_name', 'ID of the reach to which the node is associated'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', textjoin("""
                    Mandatory. In/from Prior""")],
                ])],
        ['node_id',
         odict([['dtype', 'i4'],
                ['long_name', textjoin("""
                    Nodes numbered sequentially within a Reach. Increasing in
                    the downstream direction""")],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', textjoin("""
                    Mandatory. In/from Prior""")],
                ])],
        ['time',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    UTC seconds since 2000/1/1 00:00 UTC""")],
                ['units', 's'],
                ['valid_min', 0],
                ['valid_max', 1e10],
                ['comment', textjoin("""TBD""")],
                ])],
        ['time_tai',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    TAI seconds since 2000/1/1 00:00 UTC""")],
                ['units', 's'],
                ['valid_min', 0],
                ['valid_max', 1e10],
                ['comment', textjoin("""TBD""")],
                ])],
        ['latitude',
         odict([['dtype', 'f8'],
                ['long_name', 'latitude'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -90],
                ['valid_max', 90],
                ['comment', 'geodetic latitude (degrees north of equator)'],
                ])],
        ['longitude',
         odict([['dtype', 'f8'],
                ['long_name', 'longitude'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', 0],
                ['valid_max', 360],
                ['comment', 'longitude (east of the prime meridian)'],
                ])],
        ['latitude_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in average Longitude of pixels in node""")],
                ['units', 'degrees'],
                ['valid_min', 0],
                ['valid_max', 10],
                ['comment', 'TBD'],
                ])],
        ['longitud_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in longitude'],
                ['units', 'degrees'],
                ['valid_min', 0],
                ['valid_max', 10],
                ['comment', 'TBD'],
                ])],
        ['height',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Node averaged water surface height with respect to the
                    geoid (m) with all corrections and geophysical fields
                    applied from pixels""")],
                ['standard_name', 'height'],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Current baseline is EGM2008. Geoid value used reported in
                    Geoid_model""")],
                ])],
        ['height_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainy in Node height wrt geoid, including
                    uncertainties of corrections, references""")],
                ['standard_name', 'height'],
                ['units', 'm'],
                ['valid_min', 0.01],
                ['valid_max', 50.0],
                ['comment', textjoin("""TBD""")],
                ])],
        ['width',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""Node average river width""")],
                ['units', 'm'],
                ['valid_min', 50.0],
                ['valid_max', 10000],
                ['comment', textjoin("""TBD""")],
                ])],
        ['width_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in node average river width""")],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""TBD""")],
                ])],
        ['area_detct',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Area of detected water pixels""")],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""TBD""")],
                ])],
        ['area_det_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in area of detected water pixels""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""TBD""")],
                ])],
        ['area_total',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Total water area, including estimate of dark water""")],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    Total estimated area including dark water. Best estimate
                    using water fraction.""")],
                ])],
        ['area_tot_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in total water area""")],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""TBD""")],
                ])],
        ['area_of_ht',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Area of pixels used to compute height""")],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    No uncertainty for this area.""")],
                ])],
        ['layovr_val',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Metric of layover effect (TBD)""")],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Placeholder: Numerical variable to supplement layover
                    flag. Could be layover area.""")],
                ])],
        ['node_dist',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Distance of observed node from a priori node location""")],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Method, nodes included to be defined""")],
                ])],
        ['xtrk_dist',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average distance of pixels in nodes to the satellite
                    ground track""")],
                ['units', 'm'],
                ['valid_min', 10000],
                ['valid_max', 65000],
                ['comment', textjoin("""TBD""")],
                ])],
        ['height2',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Centroid height of pixels in node with respect to the
                    reference ellipsoid.  Fully corrected for instrument and
                    geophysical delays, but not geophysical fields.""")],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 5000],
                ['comment', textjoin("""Nominal centroid is average""")],
                ])],
        ['height2_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in height2 estimate""")],
                ['units', 'm'],
                ['valid_min', 0.1],
                ['valid_max', 10.0],
                ['comment', textjoin("""TBD""")],
                ])],
        ['f_dark',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Indicates low signal to noise ration possible due to rain,
                    dark water, and others""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['f_frozen',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Indicates if the surface is frozen""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['f_layover',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Indicates if significant layover effect in node""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['n_good_pix',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Number of good pixels in node""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['f_quality',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Quality indicator on measurement, other quantities""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['f_partial',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Indicates that Reach is near edge and part may be lost due
                    to orbit variation""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['f_xovr_cal',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Quality of the cross-over calibrations""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD""")],
                ])],
        ['rdr_sig0',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn measured backscatter averaged for node""")],
                ['units', '1'],
                ['valid_min', -1],
                ['valid_max', 100],
                ['comment', textjoin("""TBD""")],
                ])],
        ['rdr_sig0_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn measured backscatter uncertainty for node""")],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""TBD""")],
                ])],
        ['c_sig0_atm',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    sigma0 atmospheric correction within the swath from
                    model data""")],
                ['units', '1'],
                ['valid_min', 1],
                ['valid_max', 10],
                ['comment', textjoin("""TBD""")],
                ])],
        ['geoid_hght',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Geoid model height above the ellipsoid at node location""")],
                ['units', 'm'],
                ['valid_min', -200],
                ['valid_max', 2000],
                ['comment', textjoin("""Current baseline is EGM2008""")],
                ])],
        ['earth_tide',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average height of Earth tide model for node location""")],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 100],
                ['comment', textjoin("""TBD""")],
                ])],
        ['pole_tide',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average height of pole tide model for node location""")],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['comment', textjoin("""TBD""")],
                ])],
        ['load_tide',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average height from loading by water (ocean) tide model for
                    node""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['c_dry_trop',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Numerical weather model dry tropospheric correction to
                    surface height""")],
                ['units', 'm'],
                ['valid_min', -2.5],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Negative as additive range correction.""")],
                ])],
        ['c_wet_trop',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Numerical weather model wet tropospheric correction to
                    surface height""")],
                ['units', 'm'],
                ['valid_min', -0.4],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Negative as additive range correction.""")],
                ])],
        ['c_iono',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Ionospheric model correction to surface height""")],
                ['units', 'm'],
                ['valid_min', -0.4],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Negative as additive range correction.""")],
                ])],
        ['c_xovr_cal',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn correction from crossover cal processing evaluated
                    for reach""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['c_kar_att',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Height correction from KaRIn orientation (attitude)
                    determination""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['c_h_bias',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Overall instrument system height bias""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['c_sys_cg',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn to s/c CG correction to height (m)""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['c_intr_cal',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Corrections on height deduced from instrument internal
                    calibrations if applicable""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['p_height',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Prior height estimate from DEM, first year of SWOT, or ?
                    """)],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['p_height_var',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Prior height variability from DEM, first year of SWOT, or ?
                    """)],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['p_width',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""Width from prior database""")],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 10000],
                ['comment', textjoin("""TBD""")],
                ])],
        ['p_width_var',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Prior width variability from historical data,
                    probability mask, or ?""")],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""TBD""")],
                ])],
        ['p_dist_out',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Distance from outlet of this river branch.""")],
                ['units', 'm'],
                ['valid_min', 10],
                ['valid_max', 10000],
                ['comment', textjoin("""TBD""")],
                ])],
        ['p_class',
         odict([['dtype', 'u2'],
                ['long_name', textjoin("""
                    Planform type from prior database""")],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 65535],
                ['comment', textjoin("""TBD""")],
                ])],
        ['grand_id',
         odict([['dtype', 'u2'],
                ['long_name', textjoin("""
                    Identification number of dam from GranD database""")],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 65535],
                ['comment', textjoin("""TBD""")],
                ])],
            ])

class RiverTileReaches(Product):
    ATTRIBUTES = []
    DIMENSIONS = odict([['reaches', 0]])
