'''
Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Alex Fore
'''
import os
import textwrap
import numpy as np
import fiona
import netCDF4
import datetime
import warnings
from shapely.geometry import Point, mapping
from collections import OrderedDict as odict

from SWOTRiver.products.pixcvec import L2PIXCVector
from SWOTRiver.products.product import Product, FILL_VALUES, textjoin

class L2HRRiverTile(Product):
    UID = "l2_hr_rivertile"
    ATTRIBUTES = odict()
    GROUPS = odict([
        ['nodes', 'RiverTileNodes'],
        ['reaches', 'RiverTileReaches'],
    ])

    @staticmethod
    def dump_xmls(node_xml_file, reach_xml_file):
        with open(node_xml_file, 'w') as ofp:
            RiverTileNodes.print_xml(ofp=ofp)
        with open(reach_xml_file, 'w') as ofp:
            RiverTileReaches.print_xml(ofp=ofp)

    @classmethod
    def from_riverobs(cls, node_outputs, reach_outputs, reach_collection):
        """Constructs self from riverobs outputs"""
        klass = cls()
        klass.nodes = RiverTileNodes.from_riverobs(node_outputs)
        klass.reaches = RiverTileReaches.from_riverobs(
            reach_outputs, reach_collection)
        return klass

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.nodes.update_from_pixc(pixc_file, index_file)
            self.reaches.update_from_pixc(pixc_file, index_file)

class ShapeWriterMixIn(object):
    """MixIn to support shapefile output"""
    def write_shapes(self, shp_fname):
        REMAP_DICT = {
            'i1': 'int', 'i2': 'int', 'i4': 'int',
            'u1': 'int', 'u2': 'int', 'u4': 'int',
            'f4': 'float', 'f8': 'float'}

        properties = odict([
            [key, REMAP_DICT[self.VARIABLES[key]['dtype']]] for key in
            self.VARIABLES])

        try:
            # these are for the geometry part of schema
            properties.pop('latitude')
            properties.pop('longitude')
            is_reach = False

        except KeyError:
            is_reach = True

        # add time-string
        properties_ = properties.copy()
        properties['time_str'] = 'str'

        # mash up the schema
        schema = {'geometry': 'Point', 'properties': properties}
        with fiona.open(shp_fname, 'w', 'ESRI Shapefile', schema) as ofp:
            for ii in range(self.reach_id.shape[0]):
                this_property = odict([[
                    key, np.asscalar(self[key][ii])] for key in
                    properties_])

                if is_reach:
                    point = Point(float(self.p_longitud[ii]),
                                  float(self.p_latitud[ii]))
                else:
                    point = Point(float(self.longitude[ii]),
                                  float(self.latitude[ii]))

                # add time-string
                this_property['time_str'] = (
                    datetime.datetime(2000, 1, 1) + datetime.timedelta(
                        seconds=this_property['time'])
                    ).strftime('%Y-%m-%dT%H:%M%S.%fZ')

                ofp.write({'geometry': mapping(point), 'id': ii,
                           'properties': this_property, 'type': 'Feature'})

class RiverTileNodes(Product, ShapeWriterMixIn):
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['nodes', 0]])
    VARIABLES = odict([
        ['reach_id',
         odict([['dtype', 'i4'],
                ['long_name', 'Reach with which node is associated'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', textjoin("""
                    "Mandatory. In Prior. Format:  CBBBBBRRRNNNT, where
                    C=continent, B=basin,R=reach,N=node, T=type. See PDD for
                    continent,type code details. Nodes number sequentially in
                    reach. Implementation note: Could be 4B integer with
                    current definition with all items as numbers.""")],
                ])],
        ['node_id',
         odict([['dtype', 'i4'],
                ['long_name', "Node Id"],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', textjoin("""
                    Nodes numbered sequentially within a reach. Increasing in
                    the downstream direction. Same format as reach_id.""")],
                ])],
        ['time',
         odict([['dtype', 'f8'],
                ['long_name', "Time in UTC sec"],
                ['units', 's'],
                ['valid_min', 0],
                ['valid_max', 1e10],
                ['comment', textjoin("""
                    Time of measurement in seconds in the UTC time scale since
                    1 Jan 2000 00:00:00 UTC. [tai_utc_difference] is the
                    difference between TAI and UTC reference time (seconds) for
                    the first measurement of the data set. If a leap second
                    occurs within the data set, the attribute leap_second is
                    set to the UTC time at which the leap second occurs.
                    Precision 1 microsecond. """)],
                ])],
        ['time_tai',
         odict([['dtype', 'f8'],
                ['long_name', "Time in TAI sec"],
                ['units', 's'],
                ['valid_min', 0],
                ['valid_max', 1e10],
                ['comment', textjoin("""
                    Time of measurement in seconds in the TAI time scale since
                    1 Jan 2000 00:00:00 TAI. This time scale contains no leap
                    seconds. The difference (in seconds) with time in UTC is
                    given by the attribute [time:tai_utc_difference].
                    time_tai = time + total_leap_sec. Precision 1 microsecond.
                    """)],
                ])],
        ['latitude',
         odict([['dtype', 'f8'],
                ['long_name', 'Latitude of centroid of detected pixels'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -78],
                ['valid_max', 78],
                ['comment', textjoin("""
                    TBD[Average latitude, not necessarily along stream]. 13
                    digits adequate for microdeg gives sub-meter location.""")],
                ])],
        ['longitude',
         odict([['dtype', 'f8'],
                ['long_name', 'Longitude of centroid of detected pixels'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', 0],
                ['valid_max', 360],
                ['comment', textjoin("""
                    East longitude is convention for all products.
                    TBD[Average longitude, not necessarily along stream]. 13
                    digits adequate for microdeg gives sub-meter location.""")],
                ])],
        ['latitude_u',
         odict([['dtype', 'f4'],
                ['long_name', "Uncertainty in latitude of node"],
                ['units', 'degrees'],
                ['valid_min', 0],
                ['valid_max', 10],
                ['comment', textjoin("""
                    TBD additional comment.""")],
                ])],
        ['longitud_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in longitude of node'],
                ['units', 'degrees'],
                ['valid_min', 0],
                ['valid_max', 10],
                ['comment', textjoin("""
                    TBD additional comment.""")],
                ])],
        ['height',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Node average water height with respect to geoid'],
                ['standard_name', 'height'],
                ['units', 'm'],
                ['valid_min', -500],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Node averaged water surface height with respect to the
                    geoid (m) with all corrections and geophysical fields
                    applied from pixels. Current Geoid baseline is EGM2008.
                    Value given in geoid_hght.""")],
                ])],
        ['height_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in node height'],
                ['standard_name', 'height'],
                ['units', 'm'],
                ['valid_min', 0.1],
                ['valid_max', 50.0],
                ['comment', textjoin("""
                    Uncertainy in Node height wrt geoid, including
                    uncertainties of corrections, references.""")],
                ])],
        ['width',
         odict([['dtype', 'f4'],
                ['long_name', "Node average river width"],
                ['units', 'm'],
                ['valid_min', 50.0],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Node average river width based on area_total.""")],
                ])],
        ['width_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in node width'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""TBD additional comment.""")],
                ])],
        ['area_detct',
         odict([['dtype', 'f4'],
                ['long_name', 'Area of detected water pixels'],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    TBC: Aggregation of node areas of pixels used.""")],
                ])],
        ['area_det_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in area of detected water pixels""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    TBD comment on method relative to area_total""")],
                ])],
        ['area_total',
         odict([['dtype', 'f4'],
                ['long_name', 'Total water area with estimate of dark water'],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    Total estimated area including dark water. Best estimate
                    using water fraction.""")],
                ])],
        ['area_tot_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in total water area'],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""TBD additional comment on method.""")],
                ])],
        ['area_of_ht',
         odict([['dtype', 'f4'],
                ['long_name', 'Area of pixels used to compute height'],
                ['units', 'm^2'],
                ['valid_min', 100],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    No uncertainty for this area.""")],
                ])],
        ['layovr_val',
         odict([['dtype', 'f4'],
                ['long_name', 'Metric of layover effect'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Numerical variable to supplement layover flag.  Will be
                    defined later depending on layover algorithm(s). Could be
                    layover area.""")],
                ])],
        ['node_dist',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Distance of observed node from a priori node location'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    TBD additional comment [not necessarily along stream].""")],
                ])],
        ['xtrk_dist',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Average distance of pixels in node to ground track'],
                ['units', 'm'],
                ['valid_min', 10000],
                ['valid_max', 65000],
                ['comment', textjoin("""
                    TBD:Different sign for Left/Right? Precision 1 m OK? If
                    unsigned, could be 2B.""")],
                ])],
        ['height2',
         odict([['dtype', 'f4'],
                ['long_name', 
                textjoin("""
                    Centoid height of pixels in node with respect to the
                    reference ellipsoid""")],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Centoid of height of pixels in node with respect to the
                    reference ellipsoid. Fully corrected for instrument and
                    media delays, but NOT[?] geophysical fields. Nominal
                    "centorid" is average. Nominal in Prior?""")],
                ])],
        ['height2_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in height2 estimate'],
                ['units', 'm'],
                ['valid_min', 0.1],
                ['valid_max', 10.0],
                ['comment', textjoin("""TBD additional comment.""")],
                ])],
        ['n_chan_max',
         odict([['dtype', 'u1'],
                ['long_name', 'Maximum number of channels detected in node'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""Value determined for each node.""")],
                ])],
        ['n_chan_mod',
         odict([['dtype', 'u1'],
                ['long_name', 'Mode of number of channels detected in node'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""Value determined for each node.""")],
                ])],
        ['dark_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Dark water flag'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates low signal to noise ratio possibly due to rain,
                    dark water, and other effects thqt  significantly affect
                    mesurements for this node.""")],
                ])],
        ['frozen_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Frozen flag'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates if the surface is frozen based on TBD.""")],
                ])],
        ['layover_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Layover flag'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates if significant layover effect in node. Area of
                    layover pixels in Expert.""")],
                ])],
        ['n_good_pix',
         odict([['dtype', 'u1'],
                ['long_name', 'Number of good pixels in node'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD additional comment.""")],
                ])],
        ['node_q',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Summary quality indicator on Node measurement""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    May Include instrument, model flags, obs_frac""")],
                ])],
        ['partial_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Indicator that node is partially observed'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates that node is near edge and part may be lost due
                    to orbit variation. TBD[Count lost to frozen, dark,
                    layover].""")],
                ])],
        ['xovr_cal_q',
         odict([['dtype', 'u1'],
                ['long_name', 'Quality indicator of cross-over calibration'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""TBD additional comment.""")],
                ])],
        ['rdr_sigma0',
         odict([['dtype', 'f4'],
                ['long_name', 'Averaged measured sgma0'],
                ['units', '1'],
                ['valid_min', -10],
                ['valid_max', 100],
                ['comment', textjoin("""
                    KaRIn measured backscatter (sigma0) averaged for Node. In
                    linear units, not dB, to allow for negative values.""")],
                ])],
        ['rdr_sig0_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in sigma0'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""TBD additional comment.""")],
                ])],
        ['sig0_atm_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Atmospheric sigma0 calibration'],
                ['units', '1'],
                ['valid_min', 0.5],
                ['valid_max', 5],
                ['comment', textjoin("""
                    Atmospheric sigma0 correction within the swath from model
                    data.""")],
                ])],
        ['geoid_hght',
         odict([['dtype', 'f4'],
                ['long_name', 'Geoid model height above ellipsoid'],
                ['units', 'm'],
                ['valid_min', -200],
                ['valid_max', 2000],
                ['comment', textjoin("""Current baseline is EGM2008""")],
                ])],
        ['solid_tide',
         odict([['dtype', 'f4'],
                ['long_name', 'Height of solid Earth tide model'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 1],
                ['comment', textjoin("""
                    Height of solid Earth tide model at node location.""")],
                ])],
        ['pole_tide',
         odict([['dtype', 'f4'],
                ['long_name', 'Height of Pole tide model'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 1],
                ['comment', textjoin("""
                    Height of Earth Pole tide model at node location.""")],
                ])],
        ['load_tide',
         odict([['dtype', 'f4'],
                ['long_name', 'Height from loading by ocean tide model'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 1],
                ['comment', textjoin("""
                    Height from loading by ocean tide model at node location.
                    Ocean loading extends some distance inland.""")],
                ])],
        ['dry_trop_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Model dry tropospheric correction to height'],
                ['units', 'm'],
                ['valid_min', -2.5],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Numerical weather model dry tropospheric correction to
                    surface height. To replace, subtract from height, add new
                    value with same sign convention.""")],
                ])],
        ['wet_trop_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Model wet tropospheric correction to height'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Numerical weather model dry tropospheric correction to
                    surface height. To replace, subtract from height, add new
                    value with same sign convention.""")],
                ])],
        ['iono_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Model ionospheric correction to height'],
                ['units', 'm'],
                ['valid_min', -0.1],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Ionospheric model correction to surface height. To
                    replace, subtract from height, add new value with same
                    sign convention.""")],
                ])],
        ['xovr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Crossover correction to height'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    KaRIn correction from crossover cal processing evaluated
                    for node.""")],
                ])],
        ['kar_att_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Height correction for KaRIn attitude'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Height correction from KaRIn orientation (attitude)
                    determination.""")],
                ])],
        ['h_bias_c',
         odict([['dtype', 'f4'],
                ['long_name', 'KaRIn height bias'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Overall instrument system height bias.""")],
                ])],
        ['sys_cg_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Center of gravity correction to height'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    KaRIn to s/c CG correction to height.""")],
                ])],
        ['intr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Instrument calibration correction to height'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Total of corrections to height deduced from instrument
                    internal calibration(s).""")],
                ])],
        ['p_height',
         odict([['dtype', 'f4'],
                ['long_name', 'Prior height estimate'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Prior height estimate from DEM, first year of SWOT, or
                    TBD.""")],
                ])],
        ['p_height_var',
         odict([['dtype', 'f4'],
                ['long_name', 'Prior height variability'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Prior height variability from historical data, probability
                    mask, or TBD.""")],
                ])],
        ['p_width',
         odict([['dtype', 'f4'],
                ['long_name', 'Prior width'],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Width from prior database.""")],
                ])],
        ['p_width_var',
         odict([['dtype', 'f4'],
                ['long_name', 'Prior width variability'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Prior width variability from historical data, probability
                    mask, or TBD.""")],
                ])],
        ['p_dist_out',
         odict([['dtype', 'f4'],
                ['long_name', 'Distance from outlet'],
                ['units', 'm'],
                ['valid_min', 1],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Additional comment TBD.""")],
                ])],
        ['p_class',
         odict([['dtype', 'u2'],
                ['long_name', 'Planform type'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 65535],
                ['comment', textjoin("""
                    Planform type from prior database. Type list is TBD.""")],
                ])],
        ['grand_id',
         odict([['dtype', 'u2'],
                ['long_name', 'Dam Id from GranD database'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 65535],
                ['comment', textjoin("""
                    http://www.gwsp.org/products/grand-database.html """)],
                ])],
    ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    @classmethod
    def from_riverobs(cls, node_outputs):
        """
        Constructs self from RiverObs node_outputs
        """
        klass = cls()
        if node_outputs is not None:
            klass['reach_id'] = node_outputs['reach_indx']
            klass['node_id'] = node_outputs['node_indx']
            klass['latitude'] = node_outputs['lat']
            klass['longitude'] = node_outputs['lon']
            klass['latitude_u'] = node_outputs['latitude_u']
            klass['longitud_u'] = node_outputs['longitud_u']
            klass['longitude'][klass['longitude'] < 0] += 360
            klass['height'] = node_outputs['h_n_ave']
            klass['height_u'] = node_outputs['h_n_std']
            klass['height2'] = node_outputs['h_a_ave']
            klass['height2_u'] = node_outputs['h_a_std']
            klass['width'] = node_outputs['w_area']
            klass['width_u'] = node_outputs['width_u']
            klass['area_detct'] = node_outputs['area']
            klass['area_det_u'] = node_outputs['area_u']
            klass['area_total'] = node_outputs['area']
            klass['area_tot_u'] = node_outputs['area_u']
            klass['area_of_ht'] = node_outputs['area_of_ht']
            klass['xtrk_dist'] = node_outputs['xtrack']
            klass['n_good_pix'] = node_outputs['nobs']
            klass['rdr_sigma0'] = node_outputs['rdr_sig0']
            klass['rdr_sig0_u'] = node_outputs['rdr_sig0_u']
            # compute node distance from prior
            klass['node_dist'] = np.sqrt(
                (node_outputs['x']-node_outputs['x_prior'])**2 +
                (node_outputs['y']-node_outputs['y_prior'])**2)
        return klass

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
        pixc_vec = L2PIXCVector.from_ncfile(index_file)

        pixc2rivertile_map = {
            '/pixel_cloud/solid_earth_tide': 'solid_tide',
            '/pixel_cloud/pole_tide': 'pole_tide',
            '/pixel_cloud/load_tide_sol1': 'load_tide',
            '/pixel_cloud/model_dry_tropo_cor': 'dry_trop_c',
            '/pixel_cloud/model_wet_tropo_cor': 'wet_trop_c',
            '/pixel_cloud/iono_cor_gim_ka': 'iono_c',
            '/pixel_cloud/xover_height_cor': 'xovr_cal_c',
            '/tvp/time': 'time',
            '/tvp/time_tai': 'time_tai'}

        pixc_data = {}
        with netCDF4.Dataset(pixc_file, 'r') as ifp:
            for key in pixc2rivertile_map:
                group, dset = key.split('/')[1::]
                pixc_data[key] = ifp.groups[group][dset][:]

        for inkey, outkey in pixc2rivertile_map.items():
            # subset pixel cloud data to look like pixcvec data
            if inkey.split('/')[1] == 'tvp':
                # silly hack
                subdata = pixc_data[inkey][pixc_vec.azimuth_index]
            else:
                subdata = pixc_data[inkey][pixc_vec.pixc_index]

            # index into pixcvec shaped data
            outdata = np.ones(self[outkey].shape)
            for reach_id in np.unique(self.reach_id):
                for node_id in self.node_id[self.reach_id == reach_id]:
                    pixc_mask = np.logical_and(
                        pixc_vec.reach_index == reach_id,
                        pixc_vec.node_index == node_id)

                    out_mask = np.logical_and(
                        self.reach_id == reach_id, self.node_id == node_id)

                    # TBD some other operation than mean (median?)
                    outdata[out_mask] = np.mean(subdata[pixc_mask])

            # stuff in product
            self[outkey] = outdata

class RiverTileReaches(Product, ShapeWriterMixIn):
    ATTRIBUTES = odict()
    DIMENSIONS = odict([['reaches', 0]])
    VARIABLES = odict([
        ['reach_id',
         odict([['dtype', 'i4'],
                ['long_name', 'Reach Id from prior database'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', textjoin("""
                    "Mandatory. In Prior. Format:  CBBBBBRRRNNNT, where
                    C=continent, B=basin,R=reach,N=node, T=type. See PDD for
                    continent, type code details. Nodes number sequentially in
                    reach. Implementation note: Could be 4B integer with
                    current definition with all items as numbers.""")],
                ])],
        ['time', RiverTileNodes.VARIABLES['time'].copy()],
        ['time_tai', RiverTileNodes.VARIABLES['time_tai'].copy()],
        ['p_latitud',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Latitude of the center of the reach from prior database '],
                ['units', 'degrees_north'],
                ['valid_min', -90],
                ['valid_max', 90],
                ['comment', textjoin("""
                    Along-stream center. 13 digits adequate for microdeg gives
                    sub-meter location.""")],
                ])],
        ['p_longitud',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    East Longitude of the center of the reach from prior
                    database""")],
                ['units', 'degrees_east'],
                ['valid_min', 0],
                ['valid_max', 360],
                ['comment', textjoin("""
                    Along-stream center. 13 digits adequate for microdeg gives
                    sub-meter location.""")],
                ])],
        ['height',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Fitted reach surface height with respect to geoid'],
                ['units', 'm'],
                ['valid_min', -500],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Reach water surface height with respect to the geoid with
                    all corrections and geophysical fields applied. Computed
                    via analytical evaluation at nominal center of polynomial
                    fit to Node heights. Geoid value used reported in
                    Geoid_modl.""")],
                ])],
        ['height_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainy in reach height'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 20],
                ['comment', textjoin("""
                    Uncertainy in Reach average height wrt geoid, including
                    uncertainties of corrections and references, and variation
                    about the fit.""")],
                ])],
        ['slope',
         odict([['dtype', 'f4'],
                ['long_name', 'Reach surface slope with respect to the geoid'],
                ['units', '1e-6'],
                ['valid_min', -1000],
                ['valid_max', 100],
                ['comment', textjoin("""
                    Reach water surface slope with respect to the geoid with
                    all corrections and geophysical fields applied. Computed
                    via polynomial  fit to Node heights. Negative slope means
                    downstream: downstream height is lower.""")],
                ])],
        ['slope_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainy in reach slope'],
                ['units', '1e-6'],
                ['valid_min', 0],
                ['valid_max', 50],
                ['comment', textjoin("""
                    Uncertainy in Reach fitted slope, including uncertainties
                    of corrections and references, and variation about the
                    fit.""")],
                ])],
        ['loc_offset',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Location offset between prior and observed reach location'],
                ['units', 'm'],
                ['valid_min', -10000],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Along reach offset distance of observed points from nominal
                    center. Sign convention TBD.""")],
                ])],
        ['width',
         odict([['dtype', 'f4'],
                ['long_name', 'Reach average river width'],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Reach average river width based on area_total.
                    Prior value in Prior section.""")],
                ])],
        ['width_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in reach average river width'],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 100],
                ['comment', textjoin("""
                    TBD additional comment.""")],
                ])],
        ['slope2',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Reach enhanced surface slope with respect to the geoid'],
                ['units', '1e-6'],
                ['valid_min', -1000],
                ['valid_max', 100],
                ['comment', textjoin("""
                    Enhanced reach slope relative to geoid produced using
                    smoothing (window TBD) of node heights. Negative slope
                    means downstream: downstream height is lower.""")],
                ])],
        ['slope2_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in reach enhanced surface slope'],
                ['units', '1e-6'],
                ['valid_min', 0],
                ['valid_max', 50],
                ['comment', textjoin("""
                    Uncertainty in enhanced reach slope.""")],
                ])],
        ['d_x_area',
         odict([['dtype', 'f4'],
                ['long_name', 'Change in cross-sectional area'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""
                    Change in channel cross sectional area from baseline.
                    Determination of baseline is TBD.""")],
                ])],
        ['d_x_area_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in change in cross-section""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 99999],
                ['comment', textjoin("""
                    Uncertainty in reach average cross-sectional area change.
                    """)],
                ])],
        ['area_detct',
         odict([['dtype', 'f4'],
                ['long_name', 'Area of detected water pixels'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    TBC: Aggregation of node areas of pixcels used.""")],
                ])],
        ['area_det_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in area of detected water""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    Uncertainty in area of detected water pixels.""")],
                ])],
        ['area_total',
         odict([['dtype', 'f4'],
                ['long_name', 'Total water area with estimate of dark water'],
                ['units', 'm^2'],
                ['valid_min', 0],
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
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    TBD additional comment on method.""")],
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
                ['long_name', 'Metric of layover effect'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Numerical variable to supplement layover flag. Will be
                    defined later depending on layover algorithm(s). Could be
                    layover area.""")],
                ])],
        ['node_dist',
         odict([['dtype', 'f4'],
                ['long_name',
                    'Mean distance between a priori nodes and observed nodes'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    TBD comment on method, nodes included to be defined.""")],
                ])],
        ['xtrk_dist',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average distance of nodes in reach to satellite ground
                    track""")],
                ['units', 'm'],
                ['valid_min', 10000],
                ['valid_max', 65000],
                ['comment', textjoin("""
                    TBD additional comment on method.""")],
                ])],
        ['n_chan_max',
         odict([['dtype', 'u1'],
                ['long_name', 'Maximum number of channels detected in nodes'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""
                    From values determined for each node.""")],
                ])],
        ['n_chan_mod',
         odict([['dtype', 'u1'],
                ['long_name', 'Mode of number of channels detected in nodes'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""
                    From values determined for each node.""")],
                ])],
        ['discharge',
         odict([['dtype', 'f4'],
                ['long_name', 'Discharge consensus value'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Value from TBD consensus discharge algorithm. No Discharge
                    in distributed product until validated.""")],
                ])],
        ['dischg_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in consensus discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    TBD additional comment on method.""")],
                ])],
        ['discharge1',
         odict([['dtype', 'f4'],
                ['long_name', 'Discharge from model_1'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Added 3 models as placeholders; may be as many as 7. Lake
                    product has only 1 model for Storage Change.""")],
                ])],
        ['dischg1_u',
         odict([['dtype', 'f4'],
                ['long_name', 'Uncertainty in model_1 discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    TBD additional comment on method.""")],
                ])],
        ['discharge2',
         odict([['dtype', 'f4'],
                ['long_name', 'Discharge from model_2'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Second of 3 added models as placeholders.""")],
                ])],
        ['dischg2_u',
         odict([['dtype', 'f4'],
                ['long_name', 'model 2 discharge uncertainty'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    TBD additional comment on method.""")],
                ])],
        ['discharge3',
         odict([['dtype', 'f4'],
                ['long_name', 'Discharge from model_3'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Third of 3 added models as placeholders""")],
                ])],
        ['dischg3_u',
         odict([['dtype', 'f4'],
                ['long_name', 'model 3 discharge uncertainty'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    TBD additional comment on method.""")],
                ])],
        ['dark_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Dark water flag'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates low signal to noise ratio possibly due to rain,
                    dark water, and other effects thqt  significantly affect
                    mesurements for this reach.""")],
                ])],
        ['frozen_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Frozen flag'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates if the surface is frozen based on TBD.""")],
                ])],
        ['layover_f',
         odict([['dtype', 'u1'],
                ['long_name', 'Layover flag'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Indicates if significant layover effect in reach. See
                    layovr_val in Expert.""")],
                ])],
        ['partial_f',
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
        ['n_good_nod',
         odict([['dtype', 'u1'],
                ['long_name', 'Number of good nodes used in fit for height'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    TBD additional comment.""")],
                ])],
        ['obs_frac_n',
         odict([['dtype', 'f4'],
                ['long_name', 'Fraction of nodes observed in reach'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['comment', textjoin("""
                    Fraction based on number of nodes. Indicates that Reach is
                    near edge and part may be lost due to orbit variation. TBD
                    counting of nodes lost to dark water or layover.""")],
                ])],
        ['reach_q',
         odict([['dtype', 'u1'],
                ['long_name', textjoin("""
                    Summary quality indicator for reach measurement""")],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    May Include instrument, model flags, obs_frac.""")],
                ])],
        ['xovr_cal_q',
         odict([['dtype', 'u1'],
                ['long_name', 'Quality of the cross-over calibrations'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    "Method for combining for Reach is TBD. Basic because all
                    flags Basic?""")],
                ])],
        ['geoid_hght',
         odict([['dtype', 'f4'],
                ['long_name', 'Geoid height'],
                ['units', 'm'],
                ['valid_min', -200],
                ['valid_max', 2000],
                ['comment', textjoin("""
                    Geoid model height above the ellipsoid. Current baseline
                    is EGM2008.""")],
                ])],
        ['geoid_slop',
         odict([['dtype', 'f4'],
                ['long_name', 'Geoid slope'],
                ['units', '1e-6'],
                ['valid_min', -1000],
                ['valid_max', 1000],
                ['comment', textjoin("""
                    Geoid model slope in the direction of the Reach.""")],
                ])],
        ['solid_tide',
         odict([['dtype', 'f4'],
                ['long_name', 'Height of solid Earth tide'],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 100],
                ['comment', textjoin("""
                    Avg height of solid Earth tide model for Reach.""")],
                ])],
        ['pole_tide',
         odict([['dtype', 'f4'],
                ['long_name', 'Height of pole tide'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['comment', textjoin("""
                    Avg height of solid Earth Pole tide model for Reach.""")],
                ])],
        ['load_tide',
         odict([['dtype', 'f4'],
                ['long_name', 'Ocean loading tide'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Avg height from loading by water (ocean) tide model for
                    reach. Ocean loading may extend beyond tidally affected
                    reaches.""")],
                ])],
        ['dry_trop_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Model dry tropospheric correction to height'],
                ['units', 'm'],
                ['valid_min', -2.5],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Numerical weather model dry tropospheric correction to
                    surface height. To replace, subtract from height, add new
                    value with same sign convention.""")],
                ])],
        ['wet_trop_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Model wet tropospheric correction to height'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Numerical weather model dry tropospheric correction to
                    surface height. To replace, subtract from height, add new
                    value with same sign convention.""")],
                ])],
        ['iono_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Model ionospheric correction to height'],
                ['units', 'm'],
                ['valid_min', -0.4],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Ionospheric model correction to surface height. To replace,
                    subtract from height, add new value with same sign
                    convention. """)],
                ])],
        ['xovr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Crossover correction to height'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    KaRIn correction from crossover cal processing evaluated
                    for Reach.""")],
                ])],
        ['kar_att_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Height correction for KaRIn attitude'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Height correction from KaRIn orientation (attitude)
                    determination.""")],
                ])],
        ['h_bias_c',
         odict([['dtype', 'f4'],
                ['long_name', 'KaRIn height bias'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Overall instrument system height bias.""")],
                ])],
        ['sys_cg_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Center of gravity correction to height'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    KaRIn to s/c CG correction to height.""")],
                ])],
        ['intr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', 'Instrument calibration correction to height'],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""
                    Total of corrections to height deduced from instrument
                    internal calibration(s).""")],
                ])],
    ])
    for name, reference in VARIABLES.items():
        reference['dimensions'] = DIMENSIONS

    @classmethod
    def from_riverobs(cls, reach_outputs, reach_collection):
        """
        Constructs self from RiverObs node_outputs
        """
        klass = cls()
        if reach_outputs is not None:
            klass['reach_id'] = reach_outputs['reach_idx']
            klass['height'] = reach_outputs['height']
            klass['height_u'] = reach_outputs['height_u']
            klass['slope'] = reach_outputs['slope']
            klass['slope_u'] = reach_outputs['slope_u']
            klass['width'] = reach_outputs['width']
            klass['width_u'] = reach_outputs['width_u']
            klass['area_detct'] = reach_outputs['area']
            klass['area_total'] = reach_outputs['area']
            klass['area_det_u'] = reach_outputs['area_u']
            klass['area_tot_u'] = reach_outputs['area_u']
            klass['area_of_ht'] = reach_outputs['area_of_ht']
            klass['xtrk_dist'] = reach_outputs['xtrk_dist']
            klass['n_good_nod'] = reach_outputs['n_good_nod']
            klass['obs_frac_n'] = reach_outputs['frac_obs']
            klass['node_dist'] = reach_outputs['node_dist']
            klass['loc_offset'] = reach_outputs['loc_offset']

            # may not be populated depending on run config
            for inkey, outkey in {'slp_enhncd': 'slope2'}.items():
                try:
                    klass[outkey] = reach_outputs[inkey]
                except KeyError:
                    pass

            # set quality flag on less than 1/2 reach observed
            klass['partial_f'] = np.zeros(reach_outputs['frac_obs'].shape)
            klass['partial_f'][reach_outputs['frac_obs'] < 0.5] = 1

            # set quality bad if partial flag is set
            klass['quality_f'] = klass['partial_f']

            assert(len(reach_collection) == len(klass['reach_id']))

            plon = np.zeros(klass['reach_id'].shape)
            plat = np.zeros(klass['reach_id'].shape)

            # Add some prior db information from reach_collection
            for ii, reach in enumerate(reach_collection):
                plat[ii] = np.mean(reach.lat)
                plon[ii] = np.rad2deg(np.arctan2(
                    np.mean(np.sin(reach.lon)), np.mean(np.cos(reach.lon))))
                # wrap to [0, 360)
                if(plon[ii] < 0): plon[ii] += 360

            klass['p_latitud'] = plat
            klass['p_longitud'] = plon
        return klass

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
        pass
