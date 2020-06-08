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

from shapely.geometry import Point, mapping, LineString
from collections import OrderedDict as odict

from SWOTRiver.products.pixcvec import L2PIXCVector
from SWOTWater.products.product import Product, FILL_VALUES, textjoin
from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9


ATTRS_2COPY_FROM_PIXC = [
    'time_coverage_start', 'time_coverage_end', 'cycle_number', 'pass_number']

RIVERTILE_ATTRIBUTES = odict([
    ['Conventions',{'dtype': 'str' ,'value': 'CF-1.7',
        'docstr': textjoin("""
            Esri conventions as given in 'ESRI Shapefile Technical Description,
            an ESRI White Paper, July 1998'
            http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf""") }],
    # title gets overridden for SLC product
    ['title', {'dtype': 'str',
        'value': 'Level 2 KaRIn High Rate River Single Pass Vector Product',
        'docstr': 'Level 2 KaRIn High Rate River Single Pass Vector Product'}],
    ['institution', {'dtype': 'str', 'value': 'JPL',
         'docstr': 'Name of producing agency.'}],
    ['source', {'dtype': 'str', 'value': 'Ka-band radar interferometer',
        'docstr': textjoin("""
            The method of production of the original data.
            If it was model-generated, source should name the model and its
            version, as specifically as could be useful. If it is
            observational, source should characterize it (e.g., 'Ka-band radar
            interferometer').""")}],
    ['history', {'dtype': 'str',
        'docstr': textjoin("""
            UTC time when file generated. Format is:
            'YYYY-MM-DD hh:mm:ss : Creation'""")}],
    ['platform', {'dtype': 'str' ,'value':'SWOT','docstr': 'SWOT'}],
    ['references', {'dtype': 'str',
        'docstr': textjoin("""
            Published or web-based references that describe
            the data or methods used to product it. Provides version number of
            software generating product.""")}],
    ['reference_document', {'dtype': 'str',
        'docstr': textjoin("""
            Name and version of Product Description Document
            to use as reference for product.""")}],
    ['contact', {'dtype': 'str',
        'docstr': textjoin("""
            Contact information for producer of product.
            (e.g., 'ops@jpl.nasa.gov').""")}],
    ['cycle_number', {'dtype': 'i2',
        'docstr': 'Cycle number of the product granule.'}],
    ['pass_number', {'dtype': 'i2',
        'docstr': 'Pass number of the product granule.'}],
    ['continent', {'dtype': 'str',
        'docstr': 'Continent the product belongs to.'}],
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
    ['left_first_longitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner longitude for the first range line and left
            edge of the swath (degrees_east)""")}],
    ['left_first_latitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner latitude for the first range line and left
            edge of the swath (degrees_north)""")}],
    ['left_last_longitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner longitude for the last range line and left
            edge of the swath (degrees_east)""")}],
    ['left_last_latitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner latitude for the last range line and left
            edge of the swath (degrees_north)""")}],
    ['right_first_longitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner longitude for the first range line and right
            edge of the swath (degrees_east)""")}],
    ['right_first_latitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner latitude for the first range line and right
            edge of the swath (degrees_north)""")}],
    ['right_last_longitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner longitude for the last range line and right
            edge of the swath (degrees_east)""")}],
    ['right_last_latitude',  {'dtype': 'float',
        'docstr': textjoin("""
            Nominal swath corner latitude for the last range line and right
            edge of the swath (degrees_north)""")}],
    ['xref_input_l2_hr_pixc_files', {'dtype': 'str',
        'docstr': 'List of L2_HR_PIXC files used to generate data in product'}],
    ['xref_input_l2_hr_rivertile_files', {'dtype': 'str',
        'docstr': 'List of RiverTile products used to generate data in product'}],
    ['xref_static_river_db_file', {'dtype': 'str',
        'docstr': textjoin("""
            Name of static river a-priori database file used to generate data
            in product""")}],
    ['xref_l2_hr_river_sp_param_file', {'dtype': 'str',
        'docstr': textjoin("""
            Name of PGE_L2_HR_RiverSP parameter file used to generate data in
            product""")}],
    ])


for key in ['Conventions', 'title', 'platform']:
    RIVERTILE_ATTRIBUTES[key]['value'] = RIVERTILE_ATTRIBUTES[key]['docstr']


class L2HRRiverTile(Product):
    UID = "l2_hr_rivertile"
    ATTRIBUTES = odict()
    GROUPS = odict([
        ['nodes', 'RiverTileNodes'],
        ['reaches', 'RiverTileReaches'],
    ])

    def sort(self):
        """sorts self according to the PDD"""
        # sort first by reach_id, then by node_id
        node_sort_idx = np.argsort(self.nodes.node_id)
        for key, values in self.nodes.variables.items():
            self.nodes[key] = values[node_sort_idx]

        reach_sort_idx = np.argsort(self.reaches.reach_id)
        for key, values in self.reaches.variables.items():
            self.reaches[key] = values[reach_sort_idx]

    @staticmethod
    def dump_xmls(node_xml_file, reach_xml_file):
        with open(node_xml_file, 'w') as ofp:
            RiverTileNodes.print_xml(ofp=ofp, is_shapefile=True)
        with open(reach_xml_file, 'w') as ofp:
            RiverTileReaches.print_xml(ofp=ofp, is_shapefile=True)

    @classmethod
    def from_riverobs(
            cls, node_outputs, reach_outputs, reach_collection, prd_reaches):
        """Constructs self from riverobs outputs"""
        klass = cls()

        # add missing reaches and nodes
        for reach, reach_id in zip(prd_reaches, prd_reaches.reach_idx):

            # skip ghost reaches
            if reach_id % 10 == 6:
                continue

            # check for missing nodes
            mask_nodes = node_outputs['reach_indx'] == reach_id

            missing_node_ids = np.setdiff1d(
                reach.node_indx, node_outputs['node_indx'][mask_nodes])

            if len(missing_node_ids) > 0:
                for missing_node_id in np.sort(missing_node_ids):

                    rch_idx = np.argwhere(
                        reach.node_indx == missing_node_id)[0][0]

                    try:
                        insert_idx = np.where(
                            missing_node_id > node_outputs['node_indx'])[0][-1]
                    except IndexError:
                        insert_idx = 0

                    node_outputs['x_prior'] = np.insert(
                        node_outputs['x_prior'], insert_idx, reach.x[rch_idx])
                    node_outputs['y_prior'] = np.insert(
                        node_outputs['y_prior'], insert_idx, reach.y[rch_idx])
                    node_outputs['lon_prior'] = np.insert(
                        node_outputs['lon_prior'], insert_idx, reach.lon[rch_idx])
                    node_outputs['lat_prior'] = np.insert(
                        node_outputs['lat_prior'], insert_idx, reach.lat[rch_idx])
                    node_outputs['p_wse'] = np.insert(
                        node_outputs['p_wse'], insert_idx, reach.wse[rch_idx])
                    node_outputs['p_wse_var'] = np.insert(
                        node_outputs['p_wse_var'], insert_idx,
                        reach.wse_var[rch_idx])
                    node_outputs['p_width'] = np.insert(
                        node_outputs['p_width'], insert_idx,
                        reach.width[rch_idx])
                    node_outputs['p_wid_var'] = np.insert(
                        node_outputs['p_wid_var'], insert_idx,
                        reach.width_var[rch_idx])
                    node_outputs['p_dist_out'] = np.insert(
                        node_outputs['p_dist_out'], insert_idx,
                        reach.dist_out[rch_idx])
                    node_outputs['p_length'] = np.insert(
                        node_outputs['p_length'], insert_idx,
                        reach.node_length[rch_idx])
                    node_outputs['grand_id'] = np.insert(
                        node_outputs['grand_id'], insert_idx,
                        reach.grod_id[rch_idx])
                    node_outputs['n_chan_max'] = np.insert(
                        node_outputs['n_chan_max'], insert_idx,
                        reach.n_chan_max[rch_idx])
                    node_outputs['n_chan_mod'] = np.insert(
                        node_outputs['n_chan_mod'], insert_idx,
                        reach.n_chan_mod[rch_idx])
                    node_outputs['node_indx'] = np.insert(
                        node_outputs['node_indx'], insert_idx, missing_node_id)
                    node_outputs['reach_indx'] = np.insert(
                        node_outputs['reach_indx'], insert_idx, reach_id)

                    for key in ['nobs', 'nobs_h', 'node_blocked']:
                        node_outputs[key] = np.insert(
                            node_outputs[key], insert_idx, MISSING_VALUE_INT4)

                    for key in [
                        'lat', 'lon', 'x', 'y', 's', 'w_ptp', 'w_std', 'w_area',
                         'w_db', 'area', 'area_u', 'area_det', 'area_det_u',
                         'area_of_ht', 'wse', 'wse_std', 'wse_u', 'rdr_sig0',
                         'rdr_sig0_u', 'latitude_u', 'longitud_u', 'width_u',
                         'geoid_hght', 'solid_tide', 'load_tidef', 'load_tideg',
                         'pole_tide', 'flow_dir', 'dark_frac', 'xtrack',
                         'h_n_ave', 'fit_height']:
                        node_outputs[key] = np.insert(
                            node_outputs[key], insert_idx, MISSING_VALUE_FLT)

            # for missing reaches
            if reach_id not in reach_outputs['reach_idx']:
                for key in ['centerline_lon', 'centerline_lat']:
                    reach_outputs[key] = np.array(list(
                        reach_outputs[key])+[reach.metadata[key],])

                this_rch_id_up = reach.metadata['rch_id_up'].T
                reach_outputs['rch_id_up'] = np.concatenate(
                    (reach_outputs['rch_id_up'], this_rch_id_up))

                this_rch_id_dn = reach.metadata['rch_id_dn'].T
                reach_outputs['rch_id_dn'] = np.concatenate(
                    (reach_outputs['rch_id_dn'], this_rch_id_up))
                reach_outputs['n_reach_up'] = np.append(
                    reach_outputs['n_reach_up'], (this_rch_id_up>0).sum())
                reach_outputs['n_reach_dn'] = np.append(
                    reach_outputs['n_reach_dn'], (this_rch_id_dn>0).sum())
                reach_outputs['reach_idx'] = np.append(
                        reach_outputs['reach_idx'], reach_id)
                reach_outputs['p_lon'] = np.append(
                        reach_outputs['p_lon'], reach.metadata['lon'])
                reach_outputs['p_lat'] = np.append(
                        reach_outputs['p_lat'], reach.metadata['lat'])
                reach_outputs['p_wse'] = np.append(
                        reach_outputs['p_wse'], reach.metadata['wse'])
                reach_outputs['p_wse_var'] = np.append(
                        reach_outputs['p_wse_var'], reach.metadata['wse_var'])
                reach_outputs['p_width'] = np.append(
                    reach_outputs['p_width'], reach.metadata['width'])
                reach_outputs['p_wid_var'] = np.append(
                    reach_outputs['p_wid_var'], reach.metadata['width_var'])
                reach_outputs['p_n_nodes'] = np.append(
                        reach_outputs['p_n_nodes'], len(reach.x))
                reach_outputs['p_dist_out'] = np.append(
                        reach_outputs['p_dist_out'], reach.metadata['dist_out'])
                reach_outputs['p_length'] = np.append(
                    reach_outputs['p_length'], reach.metadata['reach_length'])
                reach_outputs['grand_id'] = np.append(
                    reach_outputs['grand_id'], reach.metadata['grod_id'])
                reach_outputs['n_chan_max'] = np.append(
                    reach_outputs['n_chan_max'], reach.metadata['n_chan_max'])
                reach_outputs['n_chan_mod'] = np.append(
                    reach_outputs['n_chan_mod'], reach.metadata['n_chan_mod'])
                reach_outputs['reach_id'] = np.append(
                    reach_outputs['reach_id'], MISSING_VALUE_INT9)
                reach_outputs['n_good_nod'] = np.append(
                    reach_outputs['n_good_nod'], MISSING_VALUE_INT4)
                reach_outputs['lake_flag'] = np.append(
                    reach_outputs['lake_flag'], MISSING_VALUE_INT4)

                for key in ['length', 'node_dist', 'area', 'area_u', 'area_det',
                            'area_det_u', 'area_of_ht', 'width', 'width_u',
                            'loc_offset', 'xtrk_dist', 'frac_obs',
                            'slope', 'height', 'slope_u', 'height_u',
                            'geoid_slop', 'geoid_hght', 'prior_node_s',
                            'd_x_area', 'd_x_area_u', 'dark_frac', 'slope2',
                            'metro_q', 'bam_q', 'hivdi_q', 'momma_q', 'sads_q']:
                    reach_outputs[key] = np.append(
                        reach_outputs[key], MISSING_VALUE_FLT)

                # TODO: set discharge flags based on ???

        klass.nodes = RiverTileNodes.from_riverobs(node_outputs)
        klass.reaches = RiverTileReaches.from_riverobs(
            reach_outputs, reach_collection)
        # sort them by increasing reach id
        klass.sort()
        return klass

    @classmethod
    def from_shapes(cls, node_shape_path, reach_shape_path):
        """Constructs self from shapefiles"""
        klass = cls()
        klass.nodes = RiverTileNodes.from_shapes(node_shape_path)
        klass.reaches = RiverTileReaches.from_shapes(reach_shape_path)
        return klass

    def uncorrect_tides(self):
        """Removes geoid, solid earth tide, pole tide, and load tide"""
        self.nodes.uncorrect_tides()
        self.reaches.uncorrect_tides()

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.nodes.update_from_pixc(pixc_file, index_file)
            self.reaches.update_from_pixc(pixc_file, index_file)
            self.reaches.update_from_nodes(self.nodes)

    def enforce_shapes(self, node_shapefile, reach_shapefile):
        """Checks that self contains whats in the shapefile"""
        for id_key, part, shapefile in zip(
            ['node_id', 'reach_id'], [self.nodes, self.reaches],
            [node_shapefile, reach_shapefile]):

            with fiona.open(shapefile) as fc:
                schema = fc.schema
                items = list(fc.items())

            shp_vars = {}
            for varname in schema['properties']:
                shp_vars[varname] = np.array([
                    item[1]['properties'][varname] for item in items])

            shp_ids, shp_idx = np.unique(
                shp_vars[id_key], return_index=True)
            nc_ids, nc_idx = np.unique(
                part[id_key], return_index=True)
            common_ids = np.intersect1d(shp_ids, nc_ids)

            shp_idx = shp_idx[np.isin(nc_ids, shp_ids)]
            nc_idx = nc_idx[np.isin(shp_ids, nc_ids)]

            if len(shp_idx) != len(nc_idx):
                print("Shape and NetCDF not same size!")

            for varname in shp_vars:

                # skip these ones
                if varname in ['rch_id_up', 'rch_id_dn', 'time_str']:
                    continue

                shp_var = shp_vars[varname][shp_idx]
                try:
                    nc_var = part[varname][nc_idx]
                    print(varname, (nc_var-shp_var).mean(),
                          (nc_var-shp_var).std())
                except (KeyError, AttributeError):
                    print(varname+" not found.")

    def __add__(self, other):
        """Adds other to self"""
        klass = L2HRRiverTile()
        klass.nodes = self.nodes + other.nodes
        klass.reaches = self.reaches + other.reaches
        return klass

class ShapeWriterMixIn(object):
    """MixIn to support shapefile output"""
    @staticmethod
    def get_schema(dtype, valid_max=None, valid_min=None):
        """Returns the float:13.XX schema so valid max fits"""
        if dtype in ['i1', 'i2', 'u1', 'u2']:
            schema = 'int:4'
        elif dtype in ['i4', 'u4']:
            schema = 'int:9'
        elif dtype in ['i8', 'u8']:
            schema = 'int:14'
        elif dtype in ['f4', 'f8']:
            if valid_max is None or valid_min is None:
                schema = 'float:13.3'
            else:
                num_digits_left = 1+np.log10(
                    np.max([np.abs(valid_max), np.abs(valid_min)]))
                if valid_min < 0:
                    num_digits_left += 1
                if num_digits_left >= 13:
                    schema = 'float:13'
                else:
                    schema = 'float:13.%d'%(13-num_digits_left-1)
        return schema

    def write_shape_xml(self, filename):
        """Writes the XML metadata to filename"""
        with open(filename, 'w') as ofp:
            ofp.write("<?xml version='1.0' encoding='UTF-8'?>\n")
            ofp.write('<swot_product>\n')
            ofp.write('  <global_attributes>\n')
            for key in self.ATTRIBUTES:
                value = self[key]
                if value is None or value is 'None':
                    value = ''
                ofp.write('    <{}>{}</{}>\n'.format(key, value, key))
            ofp.write('  </global_attributes>\n')
            ofp.write('  <attributes>\n')

            my_vars = odict()
            for key, value in self.VARIABLES.items():
                # Skip the geometry datasets
                if key in ['lat_prior', 'lon_prior', 'centerline_lat',
                            'centerline_lon']:
                    continue

                if key not in ['node_id', 'reach_id']:
                    try:
                        value['fill_value'] = value['_FillValue']
                    except KeyError:
                        pass

                my_vars[key] = value.copy()
                if key in ['time_str',]:
                    my_vars[key]['fill_value'] = 'no_data'

                # remove these fill values
                if key in ['rdr_pol',]:
                    my_vars[key].pop('fill_value')

            for dset, attr_dict in my_vars.items():
                ofp.write('    <{}>\n'.format(dset))
                for attrname, attrvalue in attr_dict.items():
                    # skip these ones
                    if attrname not in ['dtype', 'dimensions', '_FillValue']:
                        ofp.write('      <{}>{}</{}>\n'.format(
                            attrname, attrvalue, attrname))
                ofp.write('    </{}>\n'.format(dset))

            ofp.write('  </attributes>\n')
            ofp.write('</swot_product>\n')

    def write_shapes(self, shp_fname):
        """Writes self to a shapefile"""

        properties = odict()
        for key, var in self.VARIABLES.items():
            if key in ['rdr_pol', 'reach_id', 'node_id']:
                schema = 'str'
            else:
                schema = self.get_schema(
                    var['dtype'], var.get('valid_max'), var.get('valid_min'))
            properties[key] = schema

        try:
            # these are for the geometry part of schema
            properties.pop('lat_prior')
            properties.pop('lon_prior')
            is_reach = False

        except KeyError:
            properties.pop('centerline_lat')
            properties.pop('centerline_lon')
            is_reach = True

        # add time-string
        properties_ = properties.copy()
        properties['time_str'] = 'str'

        # mash up the schema
        schema = {'geometry': 'Point', 'properties': properties}

        # special treatment of these
        if is_reach:
            properties['rch_id_up'] = 'str'
            properties['rch_id_dn'] = 'str'
            schema['geometry'] = 'LineString'

        with fiona.open(shp_fname, 'w', 'ESRI Shapefile', schema) as ofp:
            for ii in range(self.time.shape[0]):

                this_property = odict()
                for key in properties_:
                    if np.ma.isMaskedArray(self[key]):
                        this_item = self[key].data
                    else:
                        this_item = self[key]

                    if key in ['rch_id_up', 'rch_id_dn']:
                        strings = []
                        for item in this_item[ii]:
                            if item == self.VARIABLES[key]['_FillValue']:
                                thisstr = 'no_data'
                            else:
                                thisstr = str(item)
                            strings.append(thisstr)

                        this_property[key] = ' '.join(strings)

                    elif key in ['reach_id', 'node_id']:
                        if this_item[ii] == self.VARIABLES[key]['_FillValue']:
                            this_item[ii] = 'no_data'
                        else:
                            this_property[key] = str(this_item[ii])

                    elif key in ['rdr_pol',]:
                        if (this_item[ii].decode() ==
                            self.VARIABLES[key]['_FillValue']):
                            this_property[key] = 'no_data'
                        else:
                            this_property[key] = this_item[ii].astype(
                                '|S1').decode()

                    else:
                        this_property[key] = np.asscalar(this_item[ii])

                if is_reach:
                    lons = self.centerline_lon[ii]
                    lats = self.centerline_lat[ii]
                    is_valid = np.abs(lats)<90

                    this_geo = LineString([
                        (x, y) for x, y in zip(lons[is_valid], lats[is_valid])])
                else:
                    this_geo = Point(float(self.lon_prior[ii]),
                                     float(self.lat_prior[ii]))

                # add time-string
                try:
                    this_property['time_str'] = (
                        datetime.datetime(2000, 1, 1) + datetime.timedelta(
                            seconds=this_property['time'])
                        ).strftime('%Y-%m-%dT%H:%M%SZ')
                except (OverflowError, ValueError):
                    this_property['time_str'] = 'no_data'

                ofp.write({'geometry': mapping(this_geo), 'id': ii,
                           'properties': this_property, 'type': 'Feature'})

        # write shape XML metadata and prj file
        self.write_shape_xml(shp_fname.replace('.shp', '.xml'))
        with open(shp_fname.replace('.shp', '.prj'), 'w') as ofp:
            ofp.write((
                'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'+
                'SPHEROID["WGS_1984",6378137,298.257223563]],'+
                'PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]\n'))


class RiverTileNodes(Product, ShapeWriterMixIn):

    ATTRIBUTES = RIVERTILE_ATTRIBUTES
    DIMENSIONS = odict([['nodes', 0]])
    VARIABLES = odict([
        ['reach_id',
         odict([['dtype', 'i8'],
                ['long_name', 'reach ID with which node is associated'],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Unique reach identifier from the prior river database. The
                    format of the identifier is CBBBBBRRRRT, where C=continent,
                    B=basin, R=reach, T=type.""")],
                ])],
        ['node_id',
         odict([['dtype', 'i8'],
                ['long_name',
                 "node ID of the node in the prior river database"],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Unique node identifier from the prior river database.
                    The format of the identifier is CBBBBBRRRRNNNT, where
                    C=continent, B=basin, R=reach, N=node, T=type.""")],
                ])],
        ['time',
         odict([['dtype', 'f8'],
                ['long_name', "time (UTC)"],
                ['standard_name', 'time'],
                ['calendar', 'gregorian'],
                ['tai_utc_difference',
                 '[value of TAI-UTC at time of first record]'],
                ['leap_second', 'YYYY-MM-DD hh:mm:ss'],
                ['units', 'seconds since 2000-01-01 00:00:00.000'],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""
                    Time of measurement in seconds in the UTC time scale since
                    1 Jan 2000 00:00:00 UTC. [tai_utc_difference] is the
                    difference between TAI and UTC reference time (seconds) for
                    the first measurement of the data set. If a leap second
                    occurs within the data set, the metadata leap_second is
                    set to the UTC time at which the leap second occurs.""")],
                ])],
        ['time_tai',
         odict([['dtype', 'f8'],
                ['long_name', "time (TAI)"],
                ['standard_name', 'time'],
                ['calendar', 'gregorian'],
                ['units', 'seconds since 2000-01-01 00:00:00.000'],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""
                    Time of measurement in seconds in the TAI time scale since
                    1 Jan 2000 00:00:00 TAI. This time scale contains no leap
                    seconds. The difference (in seconds) with time in UTC is
                    given by the metadata [time:tai_utc_difference].""")],
                ])],
        ['time_str',
         odict([['dtype', 'f8'],
                ['long_name', "time (UTC)"],
                ['standard_name', 'time'],
                ['calendar', 'gregorian'],
                ['tai_utc_difference',
                 '[value of TAI-UTC at time of first record]'],
                ['leap_second', 'YYYY-MM-DD hh:mm:ss'],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""
                    Time string giving UTC time.  The format is
                    YYYY-MM-DDThh:mm:ssZ, where the Z suffix indicates
                    UTC time.""")],
                ])],
        ['lat',
         odict([['dtype', 'f8'],
                ['long_name', 'latitude of centroid of water-detected pixels'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -80],
                ['valid_max', 80],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""
                    Geodetic latitude of the centroid of water-detected pixels
                    assigned to the node.  Positive latitude values increase
                    northward from the equator.""")],
                ])],
        ['lon',
         odict([['dtype', 'f8'],
                ['long_name', 'longitude of centroid of water-detected pixels'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', -180],
                ['valid_max', 180],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""
                    Geodetic longitude of the centroid of water-detected
                    pixels assigned to the node.  The longitude values become
                    more positive to the east and more negative to the west of
                    the Prime Meridian.""")],
                ])],
        ['lat_u',
         odict([['dtype', 'f8'],
                ['long_name', "uncertainty in the node latitude"],
                ['units', 'degrees_north'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty in the latitude of the centroid
                    of water-detected pixels assigned to the node.""")],
                ])],
        ['lon_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in the node longitude'],
                ['units', 'degrees_east'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty in the longitude of the
                    centroid of water-detected pixels assigned to the node.
                    """)],
                ])],
        ['lat_prior',
         odict([['dtype', 'f8'],
                ['long_name', 'Latitude of prior node in database'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -78],
                ['valid_max', 78],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    TBD[Average latitude, not necessarily along stream]. 13
                    digits adequate for microdeg gives sub-meter location.""")],
                ])],
        ['lon_prior',
         odict([['dtype', 'f8'],
                ['long_name', 'Longitude of prior node in database'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', -180],
                ['valid_max', 180],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    East longitude is convention for all products.
                    TBD[Average longitude, not necessarily along stream]. 13
                    digits adequate for microdeg gives sub-meter location.""")],
                ])],
        ['wse',
         odict([['dtype', 'f8'],
                ['long_name',
                 'water surface elevation with respect to the geoid'],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Fitted node water surface elevation, relative to the
                    provided model of the geoid (geoid_hght), with all
                    corrections for media delays (wet and dry troposphere,
                    and ionosphere), crossover correction, and tidal effects
                    (solid_tide, load_tidef, and pole_tide) applied.""")],
                ])],
        ['wse_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'total uncertainty in the water surface elevation'],
                ['units', 'm'],
                ['valid_min', 0.0],
                ['valid_max', 999999],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in the
                    node WSE, including uncertainties of corrections, and
                    variation about the fit.""")],
                ])],
        ['wse_r_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'random-only uncertainty in the water surface elevation'],
                ['units', 'm'],
                ['valid_min', 0.0],
                ['valid_max', 999999],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Random-only uncertainty component in the node WSE,
                    including uncertainties of corrections, and variation about
                    the fit.""")],
                ])],
        ['width',
         odict([['dtype', 'f8'],
                ['long_name', "node width"],
                ['units', 'm'],
                ['valid_min', 0.0],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""Node width.""")],
                ])],
        ['width_u',
         odict([['dtype', 'f8'],
                ['long_name', 'total uncertainty in the node width'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in
                    the node width.""")],
                ])],
        ['area_total',
         odict([['dtype', 'f8'],
                ['long_name', 'total water surface area including dark water'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total estimated water surface area, including dark water
                    that was not detected as water in the SWOT observation but
                    identified through the use of a prior water likelihood
                    map.""")],
                ])],
        ['area_tot_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in the total water surface area'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in the
                    total estimated water surface area area_total.""")],
                ])],
        ['area_detct',
         odict([['dtype', 'f8'],
                ['long_name', 'surface area of detected water pixels'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Surface area of node that was detected as water by the
                    SWOT observation.""")],
                ])],
        ['area_det_u',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    uncertainty in the surface area of detected water""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in
                    the surface area of the detected water pixels.""")],
                ])],
        ['area_wse',
         odict([['dtype', 'f8'],
                ['long_name', 'area used to compute water surface elevation'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Surface area of the node that contributed to the
                    computation of the WSE.""")],
                ])],
        ['layovr_val',
         odict([['dtype', 'f8'],
                ['long_name', 'metric of layover effect'],
                ['units', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Value indicating an estimate of the height error due to
                    layover (TBD).""")],
                ])],
        ['node_dist',
         odict([['dtype', 'f8'],
                ['long_name',textjoin("""
                    distance between observed and prior river database node
                    location""")],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Distance between the observed node location and the node
                    location in the prior river database. """)],
                ])],
        ['xtrk_dist',
         odict([['dtype', 'f8'],
                ['long_name', 'distance to the satellite ground track'],
                ['units', 'm'],
                ['valid_min', -75000],
                ['valid_max', 75000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Distance of the observed node location from the spacecraft
                    nadir track.  A negative value indicates the left side of
                    the swath, relative to the spacecraft velocity vector.  A
                    positive value indicates the right side of the swath.""")],
                ])],
        ['flow_angle',
         odict([['dtype', 'f8'],
                ['long_name',
                 'river flow direction relative to the satellite ground track'],
                ['units', 'degrees'],
                ['valid_min', 0],
                ['valid_max', 360],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    River flow direction for the node relative to the direction
                    of the spacecraft nadir track based on prior data (not the
                    SWOT observation).  A value of zero indicates that the flow
                    is in the same direction as the spacecraft velocity
                    direction.  A value of 90 degrees indicates that the flow
                    is toward the right side of the SWOT swath.""")],
                ])],
        ['node_q',
         odict([['dtype', 'i2'],
                ['long_name', 'summary quality indicator for the node'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""good bad""")],
                ['flag_masks', 'TBD'],
                ['flag_values', np.array([0, 1]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Summary quality indicator for the node measurement.
                    Values of 0 and 1 indicate nominal and off-nominal
                    measurements.""")],
                ])],
        ['dark_frac',
         odict([['dtype', 'f8'],
                ['long_name', 'fractional area of dark water'],
                ['units', 1],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Fraction of node area_total covered by dark water.""")],
                ])],
        ['ice_clim_f',
         odict([['dtype', 'i2'],
                ['long_name', 'climatological ice cover flag'],
                ['standard_name', 'status_flag'],
                ['source', 'Yang et al. (2020)'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Climatological ice cover flag indicating whether the node
                    is ice-covered on the day of the observation based on
                    external climatological information (not the SWOT
                    measurement).  Values of 0, 1, and 2 indicate that the node
                    is likely not ice covered, likely partially ice covered,
                    and likely fully ice covered, respectively""")],
                ])],
        ['ice_dyn_f',
         odict([['dtype', 'i2'],
                ['long_name', 'dynamical ice cover flag'],
                ['standard_name', 'status_flag'],
                ['source', 'Yang et al. (2020)'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Dynamic ice cover flag indicating whether the surface is
                    ice-covered on the day of the observation based on
                    analysis of external satellite optical data.  Values of
                    0, 1, and 2 indicate that the node is not ice covered,
                    partially ice covered, and fully ice covered, respectively.
                    """)],
                ])],
        ['partial_f',
         odict([['dtype', 'i2'],
                ['long_name', 'partial node coverage flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""covered not_covered""")],
                ['flag_values', np.array([0, 1]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Flag that indicates only partial node coverage.  The
                    flag is 0 if at least 10 pixels have a valid WSE
                    measurement; the flag is 1 otherwise and node-level
                    quantities are not computed.""")],
                ])],
        ['n_good_pix',
         odict([['dtype', 'i4'],
                ['long_name', 'number of pixels that have a valid WSE'],
                ['units', 1],
                ['valid_min', 0],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Number of pixels assigned to the node that have a valid
                    node WSE.""")],
                ])],
        ['xovr_cal_q',
         odict([['dtype', 'i2'],
                ['long_name', 'quality of the cross-over calibration'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Quality of the cross-over calibration.""")],
                ])],
        ['rdr_sig0',
         odict([['dtype', 'f8'],
                ['long_name', 'sigma0'],
                ['units', '1'],
                ['valid_min', -1000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Median of the sigma0 from the pixel cloud points assigned
                    to the node in determining the node WSE.  The value is
                    provided as a linear power ratio, not a value in decibels. 
                    A decibel value may be obtained from: rdr_sig0_in_dB =
                    10*log10(rdr_sig0).""")],
                ])],
        ['rdr_sig0_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in sigma0'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 1000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Uncertainty of sig0. The value is provided in linear units.
                    This value is a one-sigma additive (not multiplicative)
                    uncertainty term, which can be added to or subtracted from
                    rdr_sig0.""")],
                ])],
        ['rdr_pol',
         odict([['dtype', 'S1'],
                ['long_name', 'polarization of sigma0'],
                ['_FillValue', 'no_data'],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Flag indicating whether the node is observed with a
                    horizontal (H) or vertical (V) signal polarization.""")],
                ])],
        ['geoid_hght',
         odict([['dtype', 'f8'],
                ['long_name', 'geoid height'],
                ['standard_name','geoid_height_above_reference_ellipsoid'],
                ['source', 'EGM2008 (Pavlis et al., 2012)'],
                ['institution', 'GSFC'],
                ['units', 'm'],
                ['valid_min', -150],
                ['valid_max', 150],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Geoid height above the reference ellipsoid with a
                    correction to refer the value to the mean tide system
                    i.e., includes the permanent tide (zero frequency).""")],
                ])],
        ['solid_tide',
         odict([['dtype', 'f8'],
                ['long_name', 'solid Earth tide height'],
                ['source', textjoin("""
                    Cartwright and Taylor (1971) and Cartwright and Edden
                    (1973)""")],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Solid-Earth (body) tide height. The zero-frequency
                    permanent tide component is not included.""")],
                ])],
        ['load_tidef',
         odict([['dtype', 'f8'],
                ['long_name', 'geocentric load tide height (FES)'],
                ['source', 'FES2014b (Carrere et al., 2016)'],
                ['institution', 'LEGOS/CNES'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Geocentric load tide height. The effect of the ocean tide
                    loading of the Earth's crust. This value is used to
                    compute wse.""")],
                ])],
        ['load_tideg',
         odict([['dtype', 'f8'],
                ['long_name', 'geocentric load tide height (GOT)'],
                ['source', 'GOT4.10c (Ray, 2013)'],
                ['institution', 'GSFC'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Geocentric load tide height. The effect of the ocean tide
                    loading of the Earth's crust.""")],
                ])],
        ['pole_tide',
         odict([['dtype', 'f8'],
                ['long_name', 'geocentric pole tide height'],
                ['source', 'Wahr (1985) and Desai et al. (2015)'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Geocentric pole tide height.  The sum total of the
                    contribution from the solid-Earth (body) pole tide height
                    and the load pole tide height (i.e., the effect of the
                    ocean pole tide loading of the Earth's crust).""")],
                ])],
        ['dry_trop_c',
         odict([['dtype', 'f8'],
                ['long_name', 'dry troposphere vertical correction'],
                ['source', 'European Centre for Medium-Range Weather Forecasts'],
                ['institution', 'ECMWF'],
                ['units', 'm'],
                ['valid_min', -3.0],
                ['valid_max', -1.5],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Equivalent vertical correction due to dry troposphere
                    delay.  Adding the reported correction to the reported
                    reach WSE results in the uncorrected reach WSE.""")],
                ])],
        ['wet_trop_c',
         odict([['dtype', 'f8'],
                ['long_name', 'wet troposphere vertical correction'],
                ['source', 'European Centre for Medium-Range Weather Forecasts'],
                ['institution', 'ECMWF'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 0],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Equivalent vertical correction due to wet troposphere
                    delay.  Adding the reported correction to the reported
                    reach WSE results in the uncorrected reach WSE.""")],
                ])],
        ['iono_c',
         odict([['dtype', 'f8'],
                ['long_name', 'ionosphere vertical correction'],
                ['source', 'Global Ionosphere Maps'],
                ['institution', 'JPL'],
                ['units', 'm'],
                ['valid_min', -0.5],
                ['valid_max', 0],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Equivalent vertical correction due to ionosphere delay. 
                    Adding the reported correction to the reported reach WSE
                    results in the uncorrected reach WSE.""")],
                ])],
        ['xovr_cal_c',
         odict([['dtype', 'f8'],
                ['long_name', 'WSE correction from KaRIn crossovers'],
                ['units', 'm'],
                ['valid_min', -10],
                ['valid_max', 10],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Height correction from KaRIn crossover calibration. 
                    The correction is applied before geolocation but reported
                    as an equivalent height correction.""")],
                ])],
        ['p_wse',
         odict([['dtype', 'f8'],
                ['long_name', 'node water surface elevation'],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Node WSE from the prior river database.""")],
                ])],
        ['p_wse_var',
         odict([['dtype', 'f8'],
                ['long_name', 'node water surface elevation variability'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Node WSE spatial variability from the prior river
                    database.""")],
                ])],
        ['p_width',
         odict([['dtype', 'f8'],
                ['long_name', 'node width'],
                ['units', 'm'],
                ['valid_min', 10],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Node width from prior river database.""")],
                ])],
        ['p_wid_var',
         odict([['dtype', 'f8'],
                ['long_name', 'node width variability'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Node width spatial variability from the prior river
                    database.""")],
                ])],
        ['p_dist_out',
         odict([['dtype', 'f8'],
                ['long_name', 'distance from the node to the outlet'],
                ['units', 'm'],
                ['valid_min', -10000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Along-stream distance from the node to the outlet, from
                    the prior river database.""")],
                ])],
        ['p_length',
         odict([['dtype', 'f8'],
                ['long_name', 'length of node'],
                ['units', 'm'],
                ['valid_min', 10],
                ['valid_max', 1000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Length of the node from the prior river database.
                    This value is used to compute the node width from the
                    water surface area.""")],
                ])],
        ['p_dam_id',
         odict([['dtype', 'i2'],
                ['long_name', 'dam ID from GRanD database'],
                ['source', 'Lehner et al. (2011)'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Dam ID from the Global Reservoir and Dam (GRanD) database.
                    The value is 0 if there is no influence of dams at the
                    node, and a positive value indicates there is an influence
                    of a dam at the node. The value of grand_id identifies the
                    dam ID in the GRanD database.  Nodes influenced by dams
                    are indicated by the type code in node_id.  """)],
                ])],
        ['p_n_ch_max',
         odict([['dtype', 'i2'],
                ['long_name', 'maximum number of channels detected in node'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Maximum number of channels at the node, from the prior
                    river database.""")],
                ])],
        ['p_n_ch_mod',
         odict([['dtype', 'i2'],
                ['long_name', 'mode of the number of channels at the node'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'lon lat'],
                ['comment', textjoin("""
                    Mode of the number of channels at the node, from the prior
                    river database.""")],
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
            klass['lat'] = node_outputs['lat']
            klass['lon'] = node_outputs['lon']
            klass['lat_u'] = node_outputs['latitude_u']
            klass['lon_u'] = node_outputs['longitud_u']
            klass['wse'] = node_outputs['wse']
            klass['wse_r_u'] = node_outputs['wse_std']
            klass['width'] = node_outputs['w_area']
            klass['width_u'] = node_outputs['width_u']
            klass['area_detct'] = node_outputs['area_det']
            klass['area_det_u'] = node_outputs['area_det_u']
            klass['area_total'] = node_outputs['area']
            klass['area_tot_u'] = node_outputs['area_u']
            klass['area_wse'] = node_outputs['area_of_ht']
            klass['xtrk_dist'] = node_outputs['xtrack']
            klass['n_good_pix'] = node_outputs['nobs']
            klass['rdr_sig0'] = node_outputs['rdr_sig0']
            klass['rdr_sig0_u'] = node_outputs['rdr_sig0_u']
            klass['geoid_hght'] = node_outputs['geoid_hght']
            klass['solid_tide'] = node_outputs['solid_tide']
            klass['load_tidef'] = node_outputs['load_tidef']
            klass['load_tideg'] = node_outputs['load_tideg']
            klass['pole_tide'] = node_outputs['pole_tide']
            klass['flow_angle'] = node_outputs['flow_dir']
            # compute node distance from prior
            klass['node_dist'] = np.sqrt(
                (node_outputs['x']-node_outputs['x_prior'])**2 +
                (node_outputs['y']-node_outputs['y_prior'])**2)
            klass['dark_frac'] = node_outputs['dark_frac']

            for key in ['lat_prior', 'lon_prior', 'p_wse', 'p_wse_var',
                        'p_width', 'p_wid_var', 'p_dist_out', 'p_length',
                        'grand_id', 'n_chan_max', 'n_chan_mod']:
                klass[key] = node_outputs[key]

            # set quality flag
            klass['node_q'] = np.zeros(node_outputs['nobs'].shape).astype(
                klass.VARIABLES['node_q']['dtype'])
            klass['node_q'][node_outputs['node_blocked']==1] |= 1

        return klass

    @classmethod
    def from_shapes(cls, shape_path):
        """Constructs self from shapefiles"""
        klass = cls()
        with fiona.open(shape_path) as ifp:
            records = list(ifp)
        data = {}
        data['lon_prior'] = np.array([
            record['geometry']['coordinates'][0] for record in records])
        data['lat_prior'] = np.array([
            record['geometry']['coordinates'][1] for record in records])
        for key, reference in klass.VARIABLES.items():
            if key not in ['lat_prior', 'lon_prior']:
                data[key] = np.array([
                    record['properties'][key] for record in records])
        for key, value in data.items():
            if key in ['reach_id', 'node_id']:
                value = value.astype('int')
            setattr(klass, key, value)
        return klass

    def uncorrect_tides(self):
        """Removes geoid, solid earth tide, pole tide, and load tide"""
        mask = np.logical_and(
            self.geoid_hght > -200, self.geoid_hght < 200)
        self.wse[mask] += (
            self.geoid_hght[mask] + self.solid_tide[mask] +
            self.load_tidef[mask] + self.pole_tide[mask])

    def update_from_pixc(self, pixc_file, index_file):
        """Adds more datasets from pixc_file file using index_file"""
        pixc_vec = L2PIXCVector.from_ncfile(index_file)

        pixc2rivertile_map = {
            '/pixel_cloud/model_dry_tropo_cor': 'dry_trop_c',
            '/pixel_cloud/model_wet_tropo_cor': 'wet_trop_c',
            '/pixel_cloud/iono_cor_gim_ka': 'iono_c',
            '/pixel_cloud/height_cor_xover': 'xovr_cal_c',
            '/pixel_cloud/xover_height_cor': 'xovr_cal_c',# old format
            '/tvp/time': 'time',
            '/tvp/time_tai': 'time_tai'}

        pixc_data = {}
        with netCDF4.Dataset(pixc_file, 'r') as ifp:
            for key in pixc2rivertile_map:
                group, dset = key.split('/')[1::]
                try:
                    pixc_data[key] = ifp.groups[group][dset][:]
                except IndexError:
                    pass
            for attr in ATTRS_2COPY_FROM_PIXC:
                try:
                    value = getattr(ifp, attr)
                except AttributeError:
                    value = getattr(ifp.groups['pixel_cloud'], attr, 'None')
                self[attr] = value

        for inkey, outkey in pixc2rivertile_map.items():
            # subset pixel cloud data to look like pixcvec data
            if inkey.split('/')[1] == 'tvp':
                # silly hack
                subdata = pixc_data[inkey][pixc_vec.azimuth_index]
            else:
                try:
                    subdata = pixc_data[inkey][pixc_vec.pixc_index]
                except KeyError:
                    pass

            # index into pixcvec shaped data
            outdata = np.ones(self[outkey].shape)
            for node_id in self.node_id:
                pixc_mask = pixc_vec.node_id == node_id
                out_mask = self.node_id == node_id

                # TBD some other operation than mean (median?)
                outdata[out_mask] = np.mean(subdata[pixc_mask])

            # replace NaNs with _FillValue
            outdata[np.isnan(outdata)] = self.VARIABLES[outkey]['_FillValue']

            # stuff in product
            self[outkey] = outdata

    def __add__(self, other):
        """Adds other to self"""
        klass = RiverTileNodes()
        for key in klass.VARIABLES:
            setattr(klass, key, np.concatenate((
                getattr(self, key), getattr(other, key))))
        return klass

class RiverTileReaches(Product, ShapeWriterMixIn):

    ATTRIBUTES = RIVERTILE_ATTRIBUTES
    DIMENSIONS = odict([
        ['reaches', 0], ['reach_neighbors', 4], ['centerlines', 1000]])
    DIMENSIONS_REACHES = odict([['reaches', 0]])
    DIMENSIONS_CENTERLINES = odict([['reaches', 0], ['centerlines', 1000]])
    DIMENSIONS_REACH_NEIGHBORS = odict([['reaches', 0], ['reach_neighbors', 4]])
    VARIABLES = odict([
        ['reach_id',
         odict([['dtype', 'i8'],
                ['long_name', 'reach ID from prior river database'],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Unique reach identifier from the prior river database.
                    The format of the identifier is CBBBBBRRRRT, where
                    C=continent, B=basin, R=reach, T=type.""")],
                ])],
        ['time', RiverTileNodes.VARIABLES['time'].copy()],
        ['time_tai', RiverTileNodes.VARIABLES['time_tai'].copy()],
        ['time_str', RiverTileNodes.VARIABLES['time_str'].copy()],
        ['centerline_lat',
         odict([['dtype', 'f8'],
                ['long_name',
                 'Latitude of the centerline of the reach from prior database'],
                ['units', 'degrees_north'],
                ['valid_min', -90],
                ['valid_max', 90],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""""")],
                ])],
        ['centerline_lon',
         odict([['dtype', 'f8'],
                ['long_name',
                 'Longitude of the centerline of the reach from prior database'],
                ['units', 'degrees_east'],
                ['valid_min', 0],
                ['valid_max', 360],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""""")],
                ])],
        ['p_lat',
         odict([['dtype', 'f8'],
                ['long_name', 'latitude of the center of the reach'],
                ['standard_name', 'latitude'],
                ['units', 'degrees_north'],
                ['valid_min', -80],
                ['valid_max', 80],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Basic'],
                ['comment', textjoin("""
                    Geodetic latitude of the reach center from the prior river
                    database.  Positive values increase northward of the
                    equator.""")],
                ])],
        ['p_lon',
         odict([['dtype', 'f8'],
                ['long_name', 'longitude of the center of the reach'],
                ['standard_name', 'longitude'],
                ['units', 'degrees_east'],
                ['valid_min', -180],
                ['valid_max', 180],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['comment', textjoin("""
                    Geodetic longitude of the reach center from the prior
                    river database.  The longitude values become more positive
                    to the east and more negative to the west of the Prime
                    Meridian.""")],
                ])],
        ['wse',
         odict([['dtype', 'f8'],
                ['long_name',
                 'water surface elevation with respect to the geoid'],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Fitted reach water surface elevation, relative to the
                    provided model of the geoid (geoid_hght), with corrections
                    for media delays (wet and dry troposphere, and ionosphere),
                    crossover correction, and tidal effects (solid_tide,
                    load_tidef, and pole_tide) applied.""")],
                ])],
        ['wse_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'total uncertainty in the water surface elevation'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 999999],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in
                    the reach WSE, including uncertainties of corrections,
                    and variation about the fit.""")],
                ])],
        ['wse_r_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'random-only uncertainty in the water surface elevation'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 999999],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Random-only component of the uncertainty in the reach WSE,
                    including uncertainties of corrections, and variation about
                    the fit.""")],
                ])],
        ['slope',
         odict([['dtype', 'f8'],
                ['long_name', 'water surface slope with respect to the geoid'],
                ['units', 'm/m'],
                ['valid_min', -0.001],
                ['valid_max', 0.1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Fitted water surface slope relative to the geoid, and
                    with the same corrections and geophysical fields applied
                    as wse. The units are m/m.  The upstream or downstream
                    direction is defined by the prior river database.  A
                    positive slope means that the downstream WSE is lower.""")],
                ])],
        ['slope_u',
         odict([['dtype', 'f8'],
                ['long_name', 'total uncertainty in the water surface slope'],
                ['units', 'm/m'],
                ['valid_min', 0],
                ['valid_max', 0.1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in the
                    water surface slope, including uncertainties of corrections
                    and variation about the fit.""")],
                ])],
        ['slope_r_u',
         odict([['dtype', 'f8'],
                ['long_name', 'random uncertainty in the water surface slope'],
                ['units', 'm/m'],
                ['valid_min', 0],
                ['valid_max', 0.1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Random-only component of the uncertainty in the water
                    surface slope.""")],
                ])],
        ['slope2',
         odict([['dtype', 'f8'],
                ['long_name',
                    'enhanced water surface slope with respect to the geoid'],
                ['units', 'm/m'],
                ['valid_min', -0.001],
                ['valid_max', 0.1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Enhanced water surface slope relative to the geoid,
                    produced using a smoothing of the node wse. The upstream
                    or downstream direction is defined by the prior river
                    database.  A positive slope means that the downstream WSE
                    is lower.""")],
                ])],
        ['slope2_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in the enhanced water surface slope'],
                ['units', 'm/m'],
                ['valid_min', 0],
                ['valid_max', 0.1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in
                    the enhanced water surface slope, including uncertainties
                    of corrections and variation about the fit.""")],
                ])],
        ['slope2_r_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'random uncertainty in the enhanced water surface slope'],
                ['units', 'm/m'],
                ['valid_min', 0],
                ['valid_max', 0.1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Random-only component of the uncertainty in the enhanced
                    water surface slope.""")],
                ])],
        ['width',
         odict([['dtype', 'f8'],
                ['long_name', 'reach width'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""Reach width.""")],
                ])],
        ['width_u',
         odict([['dtype', 'f8'],
                ['long_name', 'total uncertainty in the reach width'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in
                    the reach width.""")],
                ])],
        ['area_total',
         odict([['dtype', 'f8'],
                ['long_name', 'total water surface area including dark water'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total estimated water surface area, including dark water
                    that was not detected as water in the SWOT observation but
                    identified through the use of a prior water likelihood
                    map.""")],
                ])],
        ['area_tot_u',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    uncertainty in the total water surface area""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000*200],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in the
                    total estimated water surface area area_total.""")],
                ])],
        ['area_detct',
         odict([['dtype', 'f8'],
                ['long_name', 'surface area of detected water pixels'],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Surface area of reach that was detected as water by the
                    SWOT observation.""")],
                ])],
        ['area_det_u',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    uncertainty in the surface area of detected water""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in
                    the surface area of the detected water pixels.""")],
                ])],
        ['area_wse',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    area used to compute water surface elevation""")],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 2000000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Surface area of the reach that contributed to the
                    computation of the WSE.""")],
                ])],
        ['d_x_area',
         odict([['dtype', 'f8'],
                ['long_name', 'change in cross-sectional area'],
                ['units', 'm^2'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Change in channel cross sectional area from the value
                    reported in the prior river database.""")],
                ])],
        ['d_x_area_u',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    total uncertainty of the change in the cross-sectional
                    area """)],
                ['units', 'm^2'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Total one-sigma uncertainty (random and systematic) in the
                    change in the cross-sectional area.""")],
                ])],
        ['layovr_val',
         odict([['dtype', 'f8'],
                ['long_name', 'metric of layover effect'],
                ['units', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Value indicating an estimate of the height error due to
                    layover (TBD). """)],
                ])],
        ['node_dist',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    mean distance between observed and prior river database
                    node locations""")],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Mean distance between the observed node locations and the
                    node locations in the prior river database.""")],
                ])],
        ['loc_offset',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    along-stream location offset between the observed and
                    prior reach location""")],
                ['units', 'm'],
                ['valid_min', -20000],
                ['valid_max', 20000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Location offset between the observed and prior reach
                    locations.  This is defined as the mean of the along-stream
                    node distances of only observed nodes in the reach, minus
                    the mean of all prior river database along-stream node
                    distances in the reach.""")],
                ])],
        ['xtrk_dist',
         odict([['dtype', 'f8'],
                ['long_name', textjoin("""
                    distance to the satellite ground track""")],
                ['units', 'm'],
                ['valid_min', -75000],
                ['valid_max', 75000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Average distance of the observed node locations in the
                    reach from the spacecraft nadir track.  A negative value
                    indicates the left side of the swath, relative to the
                    spacecraft velocity vector.  A positive value indicates
                    the right side of the swath.""")],
                ])],
        ['dschg_c',
         odict([['dtype', 'f8'],
                ['long_name', 'consensus discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the consensus discharge algorithm.""")],
                ])],
        ['dschg_c_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in consensus discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the
                    consensus algorithm.""")],
                ])],
        ['dschg_c_q',
         odict([['dtype', 'i4'],
                ['long_name', 'consensus discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the consensus discharge.
                    Values of 0, 1, and 2 indicate that the consensus discharge
                    is valid, questionable, and invalid, respectively.""")],
                ])],
        ['dschg_gc',
         odict([['dtype', 'f8'],
                ['long_name', 'gauge-constrained consensus discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the gauge-constrained consensus discharge
                    algorithm.""")],
                ])],
        ['dschg_gc_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'uncertainty in gauge-constrained consensus discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the gauge-constrained
                    consensus algorithm.""")],
                ])],
        ['dschg_gc_q',
         odict([['dtype', 'i4'],
                ['long_name',
                 'gauge-constrained consensus discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the gauge-constrained
                    consensus discharge. Values of 0, 1, and 2 indicate that
                    the gauge-constrained consensus discharge is valid,
                    questionable, and invalid, respectively.""")],
                ])],
        ['dschg_m',
         odict([['dtype', 'f8'],
                ['long_name', 'MetroMan discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the MetroMan discharge algorithm.""")],
                ])],
        ['dschg_m_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in MetroMan discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the
                    MetroMan algorithm.""")],
                ])],
        ['dschg_m_q',
         odict([['dtype', 'i4'],
                ['long_name', 'MetroMan discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the MetroMan discharge.
                    Values of 0, 1, and 2 indicate that the MetroMan discharge
                    is valid, questionable, and invalid, respectively.""")],
                ])],
        ['dschg_gm',
         odict([['dtype', 'f8'],
                ['long_name', 'gauge-constrained MetroMan discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the gauge-constrained MetroMan discharge
                    algorithm.""")],
                ])],
        ['dschg_gm_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'uncertainty in gauge-constrained MetroMan discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the gauge-constrained
                    MetroMan algorithm.""")],
                ])],
        ['dschg_gm_q',
         odict([['dtype', 'i4'],
                ['long_name',
                 'gauge-constrained MetroMan discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the gauge-constrained
                    MetroMan discharge. Values of 0, 1, and 2 indicate that
                    the gauge-constrained MetroMan discharge is valid,
                    questionable, and invalid, respectively.""")],
                ])],
        ['dschg_b',
         odict([['dtype', 'f8'],
                ['long_name', 'BAM discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the BAM discharge algorithm.""")],
                ])],
        ['dschg_b_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in BAM discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the
                    BAM algorithm.""")],
                ])],
        ['dschg_b_q',
         odict([['dtype', 'i4'],
                ['long_name', 'BAM discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the BAM discharge.
                    Values of 0, 1, and 2 indicate that the BAM discharge
                    is valid, questionable, and invalid, respectively.""")],
                ])],
        ['dschg_gb',
         odict([['dtype', 'f8'],
                ['long_name', 'gauge-constrained BAM discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the gauge-constrained BAM discharge
                    algorithm.""")],
                ])],
        ['dschg_gb_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'uncertainty in gauge-constrained BAM discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the gauge-constrained
                    BAM algorithm.""")],
                ])],
        ['dschg_gb_q',
         odict([['dtype', 'i4'],
                ['long_name',
                 'gauge-constrained BAM discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the gauge-constrained
                    BAM discharge. Values of 0, 1, and 2 indicate that
                    the gauge-constrained BAM discharge is valid,
                    questionable, and invalid, respectively.""")],
                ])],
        ['dschg_h',
         odict([['dtype', 'f8'],
                ['long_name', 'HiVDI discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the HiVDI discharge algorithm.""")],
                ])],
        ['dschg_h_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in HiVDI discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the
                    HiVDI algorithm.""")],
                ])],
        ['dschg_h_q',
         odict([['dtype', 'i4'],
                ['long_name', 'HiVDI discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the HiVDI discharge.
                    Values of 0, 1, and 2 indicate that the HiVDI discharge
                    is valid, questionable, and invalid, respectively.""")],
                ])],
        ['dschg_gh',
         odict([['dtype', 'f8'],
                ['long_name', 'gauge-constrained HiVDI discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the gauge-constrained HiVDI discharge
                    algorithm.""")],
                ])],
        ['dschg_gh_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'uncertainty in gauge-constrained HiVDI discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the gauge-constrained
                    HiVDI algorithm.""")],
                ])],
        ['dschg_gh_q',
         odict([['dtype', 'i4'],
                ['long_name',
                 'gauge-constrained HiVDI discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the gauge-constrained
                    HiVDI discharge. Values of 0, 1, and 2 indicate that
                    the gauge-constrained HiVDI discharge is valid,
                    questionable, and invalid, respectively.""")],
                ])],
        ['dschg_o',
         odict([['dtype', 'f8'],
                ['long_name', 'MOMMA discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the MOMMA discharge algorithm.""")],
                ])],
        ['dschg_o_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in MOMMA discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the
                    MOMMA algorithm.""")],
                ])],
        ['dschg_o_q',
         odict([['dtype', 'i4'],
                ['long_name', 'MOMMA discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the MOMMA discharge.
                    Values of 0, 1, and 2 indicate that the MOMMA discharge
                    is valid, questionable, and invalid, respectively.""")],
                ])],
        ['dschg_go',
         odict([['dtype', 'f8'],
                ['long_name', 'gauge-constrained MOMMA discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the gauge-constrained MOMMA discharge
                    algorithm.""")],
                ])],
        ['dschg_go_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'uncertainty in gauge-constrained MOMMA discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -999999999998],
                ['valid_max', 9999999999999],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the gauge-constrained
                    MOMMA algorithm.""")],
                ])],
        ['dschg_go_q',
         odict([['dtype', 'i4'],
                ['long_name',
                 'gauge-constrained MOMMA discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the gauge-constrained
                    MOMMA discharge. Values of 0, 1, and 2 indicate that
                    the gauge-constrained MOMMA discharge is valid,
                    questionable, and invalid, respectively.""")],
                ])],
        ['dschg_s',
         odict([['dtype', 'f8'],
                ['long_name', 'SADS discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the SADS discharge algorithm.""")],
                ])],
        ['dschg_s_u',
         odict([['dtype', 'f8'],
                ['long_name', 'uncertainty in SADS discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the
                    SADS algorithm.""")],
                ])],
        ['dschg_s_q',
         odict([['dtype', 'i4'],
                ['long_name', 'SADS discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the SADS discharge.
                    Values of 0, 1, and 2 indicate that the SADS discharge
                    is valid, questionable, and invalid, respectively.""")],
                ])],
        ['dschg_gs',
         odict([['dtype', 'f8'],
                ['long_name', 'gauge-constrained SADS discharge'],
                ['units', 'm^3/s'],
                ['valid_min', -10000000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Discharge from the gauge-constrained SADS discharge
                    algorithm.""")],
                ])],
        ['dschg_gs_u',
         odict([['dtype', 'f8'],
                ['long_name',
                 'uncertainty in gauge-constrained SADS discharge'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Uncertainty in the discharge from the gauge-constrained
                    SADS algorithm.""")],
                ])],
        ['dschg_gs_q',
         odict([['dtype', 'i4'],
                ['long_name',
                 'gauge-constrained SADS discharge quality flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""
                    valid questionable invalid""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates quality of the gauge-constrained
                    SADS discharge. Values of 0, 1, and 2 indicate that
                    the gauge-constrained SADS discharge is valid,
                    questionable, and invalid, respectively.""")],
                ])],
        ['reach_q',
         odict([['dtype', 'i2'],
                ['long_name', 'summary quality indicator for the reach'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""good bad""")],
                ['flag_masks', 'TBD'],
                ['flag_values', np.array([0, 1]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Summary quality indicator for the reach measurement. 
                    Values of 0 and 1 indicate nominal (good) and off-nominal
                    (suspect) measurements.""")],
                ])],
        ['dark_frac',
         odict([['dtype', 'f8'],
                ['long_name', 'fractional area of dark water'],
                ['units', 1],
                ['valid_min', -1000],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Fraction of reach area_total covered by dark water.""")],
                ])],
        ['ice_clim_f',
         odict([['dtype', 'i2'],
                ['long_name', 'climatological ice cover flag'],
                ['standard_name', 'status_flag'],
                ['source', 'Yang et al. (2020)'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Climatological ice cover flag indicating whether the reach
                    is ice-covered on the day of the observation based on
                    external climatological information (not the SWOT
                    measurement).  Values of 0, 1, and 2 indicate that the
                    reach is likely not ice covered, likely partially ice
                    covered, and likely fully ice covered, respectively.""")],
                ])],
        ['ice_dyn_f',
         odict([['dtype', 'i2'],
                ['long_name', 'dynamical ice cover flag'],
                ['standard_name', 'status_flag'],
                ['source', 'Yang et al. (2020)'],
                ['flag_meanings', textjoin("""
                    no_ice_cover partial_ice_cover full_ice_cover""")],
                ['flag_values', np.array([0, 1, 2]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 2],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Dynamic ice cover flag indicating whether the surface is
                    ice-covered on the day of the observation based on
                    analysis of external satellite optical data.  Values of
                    0, 1, and 2 indicate that the reach is not ice covered,
                    partially ice covered, and fully ice covered, respectively.
                    """)],
                ])],
        ['partial_f',
         odict([['dtype', 'i2'],
                ['long_name', 'partial reach coverage flag'],
                ['standard_name', 'status_flag'],
                ['flag_meanings', textjoin("""covered not_covered""")],
                ['flag_values', np.array([0, 1]).astype('i2')],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Flag that indicates only partial reach coverage.  The flag
                    is 0 if at least half the nodes of the reach have valid WSE
                    measurements; the flag is 1 otherwise and reach-level
                    quantities are not computed.""")],
                ])],
        ['n_good_nod',
         odict([['dtype', 'i2'],
                ['long_name',
                 'number of nodes in the reach that have a valid WSE'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Number of nodes in the reach that have a valid node WSE.
                    Note that the total number of nodes from the prior river
                    database is given by p_n_nodes.""")],
                ])],
        ['obs_frac_n',
         odict([['dtype', 'f8'],
                ['long_name', 'fraction of nodes that have a valid WSE'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Fraction of nodes (n_good_nod/p_n_nodes) in the reach that
                    have a valid node WSE.  The value is between 0 and 1.""")],
                ])],
        ['xovr_cal_q',
         odict([['dtype', 'i2'],
                ['long_name', 'quality of the cross-over calibration'],
                ['flag_meanings', textjoin("""TBD""")],
                ['flag_masks', 'TBD'],
                ['flag_values', 'TBD'],
                ['valid_min', 0],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Quality of the cross-over calibration.""")],
                ])],
        ['geoid_hght',
         odict([['dtype', 'f8'],
                ['long_name', 'geoid height'],
                ['standard_name','geoid_height_above_reference_ellipsoid'],
                ['source', 'EGM2008 (Pavlis et al., 2012)'],
                ['institution', 'GSFC'],
                ['units', 'm'],
                ['valid_min', -150],
                ['valid_max', 150],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Geoid height above the reference ellipsoid with a
                    correction to refer the value to the mean tide system
                    i.e., includes the permanent tide (zero frequency).""")],
                ])],
        ['geoid_slop',
         odict([['dtype', 'f8'],
                ['long_name', 'geoid slope'],
                ['source', 'EGM2008 (Pavlis et al., 2012)'],
                ['institution', 'GSFC'],
                ['units', 'm/m'],
                ['valid_min', -0.001],
                ['valid_max', 0.01],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Geoid slope in the along-stream direction, based upon a
                    least-square linear fit along the reach.  A positive slope
                    means that the downstream geoid model height is lower.""")],
                ])],
        ['solid_tide',
         odict([['dtype', 'f8'],
                ['long_name', 'solid Earth tide height'],
                ['source', textjoin("""
                    Cartwright and Taylor (1971) and Cartwright and Edden
                    (1973)""")],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 1],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Solid-Earth (body) tide height. The zero-frequency
                    permanent tide component is not included.""")],
                ])],
        ['load_tidef',
         odict([['dtype', 'f8'],
                ['long_name', 'geocentric load tide height (FES)'],
                ['source', 'FES2014b (Carrere et al., 2016)'],
                ['institution', 'LEGOS/CNES'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Geocentric load tide height. The effect of the ocean tide
                    loading of the Earth's crust. This value is used to
                    compute wse.""")],
                ])],
        ['load_tideg',
         odict([['dtype', 'f8'],
                ['long_name', 'geocentric load tide height (GOT)'],
                ['source', 'GOT4.10c (Ray, 2013)'],
                ['institution', 'GSFC'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Geocentric load tide height. The effect of the ocean tide
                    loading of the Earth's crust.""")],
                ])],
        ['pole_tide',
         odict([['dtype', 'f8'],
                ['long_name', 'geocentric pole tide height'],
                ['source', 'Wahr (1985) and Desai et al. (2015)'],
                ['units', 'm'],
                ['valid_min', -0.2],
                ['valid_max', 0.2],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Geocentric pole tide height.  The sum total of the
                    contribution from the solid-Earth (body) pole tide height
                    and the load pole tide height (i.e., the effect of the
                    ocean pole tide loading of the Earth's crust).""")],
                ])],
        ['dry_trop_c',
         odict([['dtype', 'f8'],
                ['long_name', 'dry troposphere vertical correction'],
                ['source', 'European Centre for Medium-Range Weather Forecasts'],
                ['institution', 'ECMWF'],
                ['units', 'm'],
                ['valid_min', -3.0],
                ['valid_max', -1.5],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Equivalent vertical correction due to dry troposphere
                    delay.  Adding the reported correction to the reported
                    reach WSE results in the uncorrected reach WSE.""")],
                ])],
        ['wet_trop_c',
         odict([['dtype', 'f8'],
                ['long_name', 'wet troposphere vertical correction'],
                ['source', 'European Centre for Medium-Range Weather Forecasts'],
                ['institution', 'ECMWF'],
                ['units', 'm'],
                ['valid_min', -1],
                ['valid_max', 0],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Equivalent vertical correction due to wet troposphere
                    delay.  Adding the reported correction to the reported
                    reach WSE results in the uncorrected reach WSE.""")],
                ])],
        ['iono_c',
         odict([['dtype', 'f8'],
                ['long_name', 'ionosphere vertical correction'],
                ['source', 'Global Ionosphere Maps'],
                ['institution', 'JPL'],
                ['units', 'm'],
                ['valid_min', -0.5],
                ['valid_max', 0],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Equivalent vertical correction due to ionosphere delay. 
                    Adding the reported correction to the reported reach WSE
                    results in the uncorrected reach WSE.""")],
                ])],
        ['xovr_cal_c',
         odict([['dtype', 'f8'],
                ['long_name', 'WSE correction from KaRIn crossovers'],
                ['units', 'm'],
                ['valid_min', -10],
                ['valid_max', 10],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert', 'Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Height correction from KaRIn crossover calibration. 
                    The correction is applied before geolocation but reported
                    as an equivalent height correction.""")],
                ])],
        ['n_reach_up',
         odict([['dtype', 'i2'],
                ['long_name', 'number of upstream reaches'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 4],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert', 'Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Number of upstream reaches, from the prior river database.
                    A value of 4 indicates 4 or more upstream reaches.""")],
                ])],
        ['n_reach_dn',
         odict([['dtype', 'i2'],
                ['long_name', 'number of downstream reaches'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 4],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Number of downstream reaches, from the prior river
                    database.  A value of 4 indicates 4 or more downstream
                    reaches.""")],
                ])],
        ['rch_id_up',
         odict([['dtype', 'i8'],
                ['long_name', 'reach_id  of upstream reaches'],
                ['units', '1'],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Values of reach_id for the upstream reaches, from the
                    prior river database.  The values are strings of
                    comma-separated lists of at most 4 reach identifiers
                    corresponding to the upstream reaches.""")],
                ])],
        ['rch_id_dn',
         odict([['dtype', 'i8'],
                ['long_name', 'reach_id  of downstream reaches'],
                ['units', '1'],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Values of reach_id for the downstream reaches, from the
                    prior river database.  The values are strings of
                    comma-separated lists of at most 4 reach identifiers
                    corresponding to the downstream reaches.""")],
                ])],
        ['p_wse',
         odict([['dtype', 'f8'],
                ['long_name', 'reach water surface elevation'],
                ['units', 'm'],
                ['valid_min', -1000],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Reach WSE from the prior river database.""")],
                ])],
        ['p_wse_var',
         odict([['dtype', 'f8'],
                ['long_name', 'reach water surface elevation variability'],
                ['units', 'm'],
                ['valid_min', 0],
                ['valid_max', 9999],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Reach WSE spatial variability from the prior river
                    database.""")],
                ])],
        ['p_width',
         odict([['dtype', 'f8'],
                ['long_name', 'reach width'],
                ['units', 'm'],
                ['valid_min', 10],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Reach width from the prior river database.""")],
                ])],
        ['p_wid_var',
         odict([['dtype', 'f8'],
                ['long_name', 'reach width variability'],
                ['units', 'm'],
                ['valid_min', 10],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Reach width spatial variability from the prior river
                    database.""")],
                ])],
        ['p_n_nodes',
         odict([['dtype', 'i2'],
                ['long_name', 'number of nodes in the reach'],
                ['units', '1'],
                ['valid_min', 1],
                ['valid_max', 500],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Number of nodes in the reach from the prior river
                    database.""")],
                ])],
        ['p_dist_out',
         odict([['dtype', 'f8'],
                ['long_name', 'distance from the reach to the outlet '],
                ['units', 'm'],
                ['valid_min', -10000],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Along-stream distance from the reach center to the outlet,
                    from the prior river database.""")],
                ])],
        ['p_length',
         odict([['dtype', 'f8'],
                ['long_name', 'length of reach'],
                ['units', 'm'],
                ['valid_min', 100],
                ['valid_max', 100000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Basic'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Length of the reach from the prior river database.  This
                    value is used to compute the reach width from the water
                    surface area.""")],
                ])],
        ['p_maf',
         odict([['dtype', 'f8'],
                ['long_name', 'mean annual flow'],
                ['units', 'm^3/s'],
                ['valid_min', 0],
                ['valid_max', 10000000],
                ['_FillValue', MISSING_VALUE_FLT],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Mean annual flow from the prior river datavase.""")],
                ])],
        ['p_dam_id',
         odict([['dtype', 'i2'],
                ['long_name', 'dam ID from GRanD database'],
                ['source', 'Lehner et al. (2011)'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 10000],
                ['_FillValue', MISSING_VALUE_INT9],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Dam ID from the Global Reservoir and Dam (GRanD) database.
                    The value is 0 if there is no influence of dams along the
                    reach, and a positive value indicates there is an influence
                    of a dam along the reach. The value of grand_id identifies
                    the dam ID in the GRanD database.  Reaches influenced by
                    dams are indicated by the type code in reach_id.""")],
                ])],
        ['p_n_ch_max',
         odict([['dtype', 'i2'],
                ['long_name',
                 'maximum number of channels detected in the reach'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Maximum number of channels in the reach from the prior
                    river database.""")],
                ])],
        ['p_n_ch_mod',
         odict([['dtype', 'i2'],
                ['long_name', 'mode of the number of channels in the reach'],
                ['units', '1'],
                ['valid_min', 0],
                ['valid_max', 100],
                ['_FillValue', MISSING_VALUE_INT4],
                ['tag_basic_expert','Expert'],
                ['coordinates', 'p_lon p_lat'],
                ['comment', textjoin("""
                    Mode of the number of channels in the reach from the
                    prior river database.""")],
                ])],
    ])
    for name, reference in VARIABLES.items():
        if name in ['rch_id_up', 'rch_id_dn']:
            reference['dimensions'] = DIMENSIONS_REACH_NEIGHBORS
        elif name in ['centerline_lat', 'centerline_lon']:
            reference['dimensions'] = DIMENSIONS_CENTERLINES
        else:
            reference['dimensions'] = DIMENSIONS_REACHES

    @classmethod
    def from_riverobs(cls, reach_outputs, reach_collection):
        """
        Constructs self from RiverObs node_outputs
        """
        klass = cls()
        if reach_outputs is not None:
            klass['reach_id'] = reach_outputs['reach_idx']
            klass['wse'] = reach_outputs['height']
            klass['wse_r_u'] = reach_outputs['height_u']
            klass['slope'] = reach_outputs['slope']
            klass['slope_r_u'] = reach_outputs['slope_u']
            klass['slope2'] = reach_outputs['slope2']
            klass['width'] = reach_outputs['width']
            klass['width_u'] = reach_outputs['width_u']
            klass['area_total'] = reach_outputs['area']
            klass['area_tot_u'] = reach_outputs['area_u']
            klass['area_detct'] = reach_outputs['area_det']
            klass['area_det_u'] = reach_outputs['area_det_u']
            klass['area_wse'] = reach_outputs['area_of_ht']
            klass['xtrk_dist'] = reach_outputs['xtrk_dist']
            klass['n_good_nod'] = reach_outputs['n_good_nod']
            klass['obs_frac_n'] = reach_outputs['frac_obs']
            klass['node_dist'] = reach_outputs['node_dist']
            klass['loc_offset'] = reach_outputs['loc_offset']
            klass['geoid_hght'] = reach_outputs['geoid_hght']
            klass['geoid_slop'] = reach_outputs['geoid_slop']
            klass['rch_id_up'] = reach_outputs['rch_id_up']
            klass['rch_id_dn'] = reach_outputs['rch_id_dn']
            klass['n_reach_up'] = reach_outputs['n_reach_up']
            klass['n_reach_dn'] = reach_outputs['n_reach_dn']
            klass['d_x_area'] = reach_outputs['d_x_area']
            klass['d_x_area_u'] = reach_outputs['d_x_area_u']
            klass['dark_frac'] = reach_outputs['dark_frac']
            klass['p_n_ch_max'] = reach_outputs['n_chan_max']
            klass['p_n_ch_mod'] = reach_outputs['n_chan_mod']

            for key in ['p_wse', 'p_wse_var', 'p_width', 'p_wid_var',
                        'p_dist_out', 'p_length', 'grand_id', 'p_n_nodes',
                        'p_lat', 'p_lon']:
                klass[key] = reach_outputs[key]

            cl_lon = klass['centerline_lon'][:]
            cl_lat = klass['centerline_lat'][:]
            for ii in range(klass.dimensions['reaches']):
                this_len = len(reach_outputs['centerline_lon'][ii])
                cl_lon[ii, 0:this_len] = reach_outputs['centerline_lon'][ii]
                cl_lat[ii, 0:this_len] = reach_outputs['centerline_lat'][ii]

            klass['centerline_lon'] = cl_lon
            klass['centerline_lat'] = cl_lat

            # set quality flag on less than 1/2 reach observed
            klass['partial_f'] = np.zeros(reach_outputs['frac_obs'].shape)
            klass['partial_f'][reach_outputs['frac_obs'] < 0.5] = 1

            # set quality bad if partial flag is set
            klass['reach_q'] = klass['partial_f']

        return klass

    @classmethod
    def from_shapes(cls, shape_path):
        """Constructs self from shapefiles"""
        klass = cls()
        with fiona.open(shape_path) as ifp:
            records = list(ifp)

        cl_fill = klass.VARIABLES['centerline_lon']['_FillValue']
        cl_len = klass.DIMENSIONS['centerlines']

        data = {}
        data['centerline_lon'] = np.ones([len(records), cl_len]) * cl_fill
        data['centerline_lat'] = np.ones([len(records), cl_len]) * cl_fill

        for irec, record in enumerate(records):
            this_cl = np.array(records[0]['geometry']['coordinates'])
            data['centerline_lon'][irec, :this_cl.shape[0]] = this_cl[:, 0]
            data['centerline_lat'][irec, :this_cl.shape[0]] = this_cl[:, 1]

        for key, reference in klass.VARIABLES.items():
            if key in ['centerline_lon', 'centerline_lat']:
                pass

            elif key in ['rch_id_up', 'rch_id_dn']:
                fill = klass.VARIABLES[key]['_FillValue']
                n_ids = klass.DIMENSIONS['reach_neighbors']
                data[key] = np.ones([len(records), n_ids])*fill
                for irec, record in enumerate(records):
                    tmp = record['properties'][key].replace(
                        'no_data', str(fill))
                    data[key][irec, :] = np.array([
                        int(item) for item in tmp.split(' ')])

            else:
                data[key] = np.array([
                    record['properties'][key] for record in records])

        for key, value in data.items():
            if key in ['reach_id', 'node_id']:
                value = value.astype('int')
            setattr(klass, key, value)
        return klass

    def uncorrect_tides(self):
        """Removes geoid, solid earth tide, pole tide, and load tide"""
        mask = np.logical_and(
            self.geoid_hght > -200, self.geoid_hght < 200)
        self.wse[mask] += (
            self.geoid_hght[mask] + self.solid_tide[mask] +
            self.load_tidef[mask] + self.pole_tide[mask])

    def update_from_pixc(self, pixc_file, index_file):
        """Copies some attributes from input PIXC file"""
        with netCDF4.Dataset(pixc_file, 'r') as ifp:
            for attr in ATTRS_2COPY_FROM_PIXC:
                try:
                    value = getattr(ifp, attr)
                except AttributeError:
                    value = getattr(ifp.groups['pixel_cloud'], attr, 'None')
                self[attr] = value

    def update_from_nodes(self, nodes):
        """Averages node things to reach things and populates self"""
        keys = ['time', 'time_tai', 'geoid_hght', 'solid_tide',
                'pole_tide', 'load_tidef', 'load_tideg', 'dry_trop_c',
                'wet_trop_c', 'iono_c', 'xovr_cal_c']

        node_reach_type = nodes.node_id % 10
        node_reach_ids = (
            np.floor(nodes.node_id / 10000).astype('int'))*10 + node_reach_type

        for key in keys:
            node_value = getattr(nodes, key)
            reach_value = getattr(self, key)
            for ii, reach_id in enumerate(self.reach_id):
                reach_value[ii] = np.mean(
                    node_value[node_reach_ids == reach_id])
            self[key] = reach_value

    def __add__(self, other):
        """Adds other to self"""
        klass = RiverTileReaches()
        for key in klass.VARIABLES:
            setattr(klass, key, np.concatenate((
                getattr(self, key), getattr(other, key))))
        return klass
