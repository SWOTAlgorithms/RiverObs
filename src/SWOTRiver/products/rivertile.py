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
        klass = cls()
        klass.nodes = RiverTileNodes.from_riverobs(node_outputs)
        klass.reaches = RiverTileReaches.from_riverobs(
            reach_outputs, reach_collection)
        return klass

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
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
                    Numerical variable to supplement layover flag""")],
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
        ['dark_f',
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
        ['frozen_f',
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
        ['layover_f',
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
        ['quality_f',
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
        ['xovr_cal_f',
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
        ['sig0_atm_c',
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
        ['dry_trop_c',
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
        ['wet_trop_c',
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
        ['iono_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Ionospheric model correction to surface height""")],
                ['units', 'm'],
                ['valid_min', -0.4],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Negative as additive range correction.""")],
                ])],
        ['xovr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn correction from crossover cal processing evaluated
                    for reach""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['kar_att_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Height correction from KaRIn orientation (attitude)
                    determination""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['h_bias_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Overall instrument system height bias""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['sys_cg_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn to s/c CG correction to height (m)""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['intr_cal_c',
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
            klass['rdr_sig0'] = node_outputs['rdr_sig0']
            # compute node distance from prior
            klass['node_dist'] = np.sqrt(
                (node_outputs['x']-node_outputs['x_prior'])**2 +
                (node_outputs['y']-node_outputs['y_prior'])**2)
        return klass

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
        pixc_vec = L2PIXCVector.from_ncfile(index_file)

        pixc2rivertile_map = {
            '/pixel_cloud/solid_earth_tide': 'earth_tide',
            '/pixel_cloud/pole_tide': 'pole_tide',
            '/pixel_cloud/loading_tide_sol1': 'load_tide',
            '/pixel_cloud/model_dry_tropo_cor': 'dry_trop_c',
            '/pixel_cloud/model_wet_tropo_cor': 'wet_trop_c',
            '/pixel_cloud/model_iono_ka_cor': 'iono_c',
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
                ['long_name', 'ID of the reach to which the node is associated'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 2147483647],
                ['comment', textjoin("""
                    Mandatory. In/from Prior""")],
                ])],
        ['time', RiverTileNodes.VARIABLES['time'].copy()],
        ['time_tai', RiverTileNodes.VARIABLES['time_tai'].copy()],
        ['p_latitud',
         odict([['dtype', 'f4'],
                ['long_name', 'prior latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -90],
                ['valid_max', 90],
                ['comment', textjoin("""
                    Latitude of the along-stream center point of the reach
                    from a priori database""")],
                ])],
        ['p_longitud',
         odict([['dtype', 'f4'],
                ['long_name', 'prior longitude'],
                ['units', 'degrees_east'],
                ['valid_min', 0],
                ['valid_max', 360],
                ['comment', textjoin("""
                    Longitude of the along-stream center point of the reach
                    from a priori database""")],
                ])],
        ['height',
         odict([['dtype', 'f4'],
                ['long_name', 'height above geoid'],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Reach average water surface height with respect to the
                    geoid with all corrections and geophysical fields applied.
                    Computed via analytical average of polynomial fit to node
                    heights.""")],
                ])],
        ['height_u',
         odict([['dtype', 'f4'],
                ['long_name', ''],
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
                ['long_name', 'slope'],
                ['units', '1e-6'],
                ['valid_min', -1000],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Reach average water surface slope with respect to the
                    geoid with all corrections and geophysical fields applied.
                    Computed via polynomial fit to Node heights.""")],
                ])],
        ['slope_u',
         odict([['dtype', 'f4'],
                ['long_name', 'slope uncertainy'],
                ['units', '1e-6'],
                ['valid_min', 0],
                ['valid_max', 50],
                ['comment', textjoin("""
                    Uncertainy in Reach fitted slope, including
                    uncertainties of corrections and references, and variation
                    about the fit.""")],
                ])],
        ['width',
         odict([['dtype', 'f4'],
                ['long_name', 'slope'],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Reach average river width based on area_total.""")],
                ])],
        ['width_u',
         odict([['dtype', 'f4'],
                ['long_name', 'width uncertainy'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['comment', textjoin("""
                    Uncertainty in reach average river width.""")],
                ])],
        ['slope2',
         odict([['dtype', 'f4'],
                ['long_name', 'slope'],
                ['units', '1e-6'],
                ['valid_min', -1000],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Enhanced Reach average slope relative to geoid produced
                    using smoothing of node heights.""")],
                ])],
        ['slope2_u',
         odict([['dtype', 'f4'],
                ['long_name', 'slope uncertainy'],
                ['units', '1e-6'],
                ['valid_min', 0],
                ['valid_max', 50],
                ['comment', textjoin("""
                    Uncertainty in enhanced Reach average slope.""")],
                ])],
        ['d_x_area',
         odict([['dtype', 'f4'],
                ['long_name', 'cross-sectional area change wrt first pass'],
                ['units', 'm^2'],
                ['valid_min', -9999999],
                ['valid_max', 9999999],
                ['comment', textjoin("""
                    Reach average cross-sectional area change with respect to
                    the first overpass.""")],
                ])],
        ['d_x_area_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in cross-sectional area change wrt first
                    pass""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 99999],
                ['comment', textjoin("""
                    Uncertainty in reach average cross-sectional area change.
                    """)],
                ])],
        ['area_detct',
         odict([['dtype', 'f4'],
                ['long_name', 'area of detected water'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    Area of detected water pixels.""")],
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
                ['long_name', 'total area of detected water'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    Total water area, including estimate of dark water.""")],
                ])],
        ['area_tot_u',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Uncertainty in total area of detected water""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['comment', textjoin("""
                    Uncertainty in total water area.""")],
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
                ['long_name', 'layover metric'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 5000],
                ['comment', textjoin("""
                    Numerical variable to supplement layover flag""")],
                ])],
        ['node_dist',
         odict([['dtype', 'f4'],
                ['long_name', 'mean distance of nodes to prior nodes'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['comment', textjoin("""
                    Mean distance of observed node from a priori node
                    location""")],
                ])],
        ['xtrk_dist',
         odict([['dtype', 'f4'],
                ['long_name', 'cross-track distance'],
                ['units', 'm'],
                ['valid_min', 10000],
                ['valid_max', 65000],
                ['comment', textjoin("""
                    Average distance from nodes in reach to the satellite
                    ground track""")],
                ])],
        ['discharge',
         odict([['dtype', 'f4'],
                ['long_name', 'discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Consensus value reach average river discharge.""")],
                ])],
        ['dischg_u',
         odict([['dtype', 'f4'],
                ['long_name', 'discharge uncertainty'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Uncertainty in consensus reach average river discharge.""")],
                ])],
        ['discharge1',
         odict([['dtype', 'f4'],
                ['long_name', 'model 1 discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Model 1 value reach average river discharge.""")],
                ])],
        ['dischg1_u',
         odict([['dtype', 'f4'],
                ['long_name', 'model 1 discharge uncertainty'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Uncertainty in model 1 reach average river discharge.""")],
                ])],
        ['discharge2',
         odict([['dtype', 'f4'],
                ['long_name', 'model 2 discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Model 2 value reach average river discharge.""")],
                ])],
        ['dischg2_u',
         odict([['dtype', 'f4'],
                ['long_name', 'model 2 discharge uncertainty'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Uncertainty in model 2 reach average river discharge.""")],
                ])],
        ['discharge3',
         odict([['dtype', 'f4'],
                ['long_name', 'model 3 discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Model 3 value reach average river discharge.""")],
                ])],
        ['dischg3_u',
         odict([['dtype', 'f4'],
                ['long_name', 'model 3 discharge uncertainty'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 9999999999999],
                ['comment', textjoin("""
                    Uncertainty in model 3 reach average river discharge.""")],
                ])],
        ['dark_f',
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
        ['frozen_f',
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
        ['layover_f',
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
        ['n_good_nod',
         odict([['dtype', 'u1'],
                ['long_name', 'number of good nodes'],
                ['flag_meanings', textjoin("""TBD""")],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 254],
                ['comment', textjoin("""
                    Number of good nodes in reach used for slope fit""")],
                ])],
        ['quality_f',
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
        ['xovr_cal_f',
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
        ['frac_obs',
         odict([['dtype', 'f4'],
                ['long_name', 'fraction of nodes in reach observed'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['comment', textjoin("""
                    TBD""")],
                ])],
        ['loc_offset',
         odict([['dtype', 'f4'],
                ['long_name', 'node mean - prior mean along-reach coordinate'],
                ['units', 'm'],
                ['valid_min', -99999],
                ['valid_max', 99999],
                ['comment', textjoin("""
                    Along reach offset distance of observed points from nominal
                    center.""")],
                ])],
        ['geoid_hght',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Geoid model height above the ellipsoid at reach location""")],
                ['units', 'm'],
                ['valid_min', -200],
                ['valid_max', 2000],
                ['comment', textjoin("""Current baseline is EGM2008""")],
                ])],
        ['earth_tide',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average height of Earth tide model for reach location""")],
                ['units', 'm'],
                ['valid_min', 50],
                ['valid_max', 100],
                ['comment', textjoin("""TBD""")],
                ])],
        ['pole_tide',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average height of pole tide model for reach location""")],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['comment', textjoin("""TBD""")],
                ])],
        ['load_tide',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Average height from loading by water (ocean) tide model for
                    reach""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['dry_trop_c',
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
        ['wet_trop_c',
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
        ['iono_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Ionospheric model correction to surface height""")],
                ['units', 'm'],
                ['valid_min', -0.4],
                ['valid_max', 0],
                ['comment', textjoin("""
                    Negative as additive range correction.""")],
                ])],
        ['xovr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn correction from crossover cal processing evaluated
                    for reach""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['kar_att_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Height correction from KaRIn orientation (attitude)
                    determination""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['h_bias_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Overall instrument system height bias""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['sys_cg_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    KaRIn to s/c CG correction to height (m)""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
                ])],
        ['intr_cal_c',
         odict([['dtype', 'f4'],
                ['long_name', textjoin("""
                    Corrections on height deduced from instrument internal
                    calibrations if applicable""")],
                ['units', 'm'],
                ['valid_min', -9999],
                ['valid_max', 9999],
                ['comment', textjoin("""TBD""")],
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
            klass['height'] = 1/3*(
                reach_outputs['h_nw'] + reach_outputs['h_no'] +
                reach_outputs['h_nr'])
            klass['slope'] = 1/3*(
                reach_outputs['slp_nw'] + reach_outputs['slp_no'] +
                reach_outputs['slp_nr'])
            klass['width'] = reach_outputs['width']
            klass['width_u'] = reach_outputs['width_u']
            klass['slope2'] = reach_outputs['slp_enhncd']
            klass['area_detct'] = reach_outputs['area']
            klass['area_total'] = reach_outputs['area']
            klass['area_det_u'] = reach_outputs['area_u']
            klass['area_tot_u'] = reach_outputs['area_u']
            klass['area_of_ht'] = reach_outputs['area_of_ht']
            klass['xtrk_dist'] = reach_outputs['xtrk_dist']
            klass['n_good_nod'] = reach_outputs['n_good_nod']
            klass['frac_obs'] = reach_outputs['frac_obs']
            klass['loc_offset'] = reach_outputs['loc_offset']

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
