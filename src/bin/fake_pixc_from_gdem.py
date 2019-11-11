#!/usr/bin/env python
"""
fake_pixc_from_gdem.py: Makes a fake pixc from a gdem

Useage: fake_pixc_from_gdem.py pixc.nc gdem.nc fake_pixc.nc

Author (s): Alex Fore

To run it through RiverObs, run swot_pixc2rivertile.py with this config:

width_db_file             (-) = None
use_width_db              (-) = False
reach_db_path             (-) = /u/turner-z0/fore/work/rivertile/new-reach-db/20190415
class_list                (-) = [4,24]
use_fractional_inundation (-) = [False,]
use_segmentation          (-) = [True,]
use_heights               (-) = [True,]
min_points                (-) = 100
clip_buffer               (-) = 0.1
ds                        (-) = None
refine_centerline         (-) = False
smooth                    (-) = 0.01
alpha                     (-) = 1
max_iter                  (-) = 1
scalar_max_width          (-) = 600.0
minobs                    (-) = 10
trim_ends                 (-) = False
min_fit_points            (-) = 3
do_improved_geolocation   (-) = False
geolocation_method        (-) = taylor
"""
import argparse
import os
import netCDF4
import numpy as np
from SWOTWater.constants import PIXC_CLASSES

def fake_pixc_from_gdem(
    gdem_file, pixc_file, fake_pixc_file, subsample_factor=2, dark_water_thresh=1):
    """Fakes a pixel cloud file from a gdem file"""
    with netCDF4.Dataset(gdem_file, 'r') as ifp_gdem,\
         netCDF4.Dataset(pixc_file, 'r') as ifp_pixc,\
         netCDF4.Dataset(fake_pixc_file, 'w') as ofp:

        # copy attributes
        for attr in ifp_pixc.__dict__:
            value = ifp_pixc.__dict__[attr]
            setattr(ofp, attr, value)

        # copy tvp data
        ofp.createGroup('tvp')
        ofp.groups['tvp'].createDimension('nr_tvps', 0)
        for key, value in ifp_pixc.groups['tvp'].variables.items():
            var = ofp.createVariable('/tvp/'+key, value.dtype.str, ('nr_tvps',))
            var[:] = value[:]

        # copy noise data
        try:
            noise_grp = ifp_pixc.groups['noise']
            ofp.createGroup('noise')
            ofp.groups['noise'].createDimension('num_lines', 0)
            for key, value in noise_grp.variables.items():
                var = ofp.createVariable(
                    '/noise/'+key, value.dtype.str, ('num_lines',))
                var[:] = value[:]
        except KeyError:
            pass

        # copy pixel_cloud attributes
        ofp.createGroup('pixel_cloud')
        ofp.groups['pixel_cloud'].createDimension('record', 0)
        ofp.groups['pixel_cloud'].createDimension('complex_depth', 2)

        for attr in ifp_pixc.groups['pixel_cloud'].__dict__:
            value = ifp_pixc.groups['pixel_cloud'].__dict__[attr]
            setattr(ofp.groups['pixel_cloud'], attr, value)

        landtype = ifp_gdem.variables['landtype'][:][::subsample_factor]
        latitude = ifp_gdem.variables['latitude'][:][::subsample_factor]
        longitude = ifp_gdem.variables['longitude'][:][::subsample_factor]
        longitude[longitude>180] -= 360
        elevation = ifp_gdem.variables['elevation'][:][::subsample_factor]
        try:
            media_attenuation = ifp_gdem['media_attenuation'][:][::subsample_factor]
        except IndexError:
            media_attenuation = np.zeros_like(elevation) + dark_water_thresh + 1
        cross_track_ = ifp_gdem.variables['cross_track'][:]
        range_spacing = ifp_gdem.ground_spacing
        azimuth_spacing = ifp_gdem.azimuth_spacing

        gdem_shape = landtype.shape
        range_index, azimuth_index = np.meshgrid(
            np.arange(gdem_shape[1]), np.arange(gdem_shape[0]))
        cross_track, tmp = np.meshgrid(cross_track_, np.arange(gdem_shape[0]))

        # set azimuth and range size
        ofp.groups['pixel_cloud'].interferogram_size_azimuth = gdem_shape[0]
        ofp.groups['pixel_cloud'].interferogram_size_range = gdem_shape[1]

        mask = landtype == 1
        pixc_shape = range_index[mask].shape
        pixel_area = subsample_factor * range_spacing * azimuth_spacing

        tvp_time = ifp_pixc.groups['tvp'].variables['time'][:]

        # set landtypes to pixc classes
        landtype[landtype==1] = PIXC_CLASSES['open_water']
        dark_water_mask = media_attenuation < dark_water_thresh
        landtype[dark_water_mask] = PIXC_CLASSES['dark_water']

        out_pixc_dsets = {}
        out_pixc_dsets['range_index'] = range_index[mask]
        out_pixc_dsets['azimuth_index'] = azimuth_index[mask]
        out_pixc_dsets['classification'] = landtype[mask]
        out_pixc_dsets['water_frac'] = out_pixc_dsets[
            'classification'].copy()*0 + 1.0
        out_pixc_dsets['latitude'] = latitude[mask]
        out_pixc_dsets['longitude'] = longitude[mask]
        out_pixc_dsets['height'] = elevation[mask]
        out_pixc_dsets['cross_track'] = cross_track[mask]
        out_pixc_dsets['illumination_time'] = tvp_time[azimuth_index[mask]]
        out_pixc_dsets['num_rare_looks'] = np.zeros(pixc_shape) + subsample_factor
        out_pixc_dsets['pixel_area'] =  np.zeros(pixc_shape) + pixel_area

        out_shape = out_pixc_dsets['range_index'].shape

        for varname in ifp_pixc.groups['pixel_cloud'].variables:
            if varname not in out_pixc_dsets:
                if varname == 'interferogram':
                    value = np.zeros(list(out_shape)+[2,])
                else:
                    value = np.zeros(out_shape)
                out_pixc_dsets[varname] = value

        for varname, varvalue in out_pixc_dsets.items():
            if varname == 'interferogram':
                var = ofp.createVariable(
                    '/pixel_cloud/'+varname, varvalue.dtype.str,
                    ('record', 'complex_depth'))
            else:
                var = ofp.createVariable(
                    '/pixel_cloud/'+varname, varvalue.dtype.str, ('record',))
            var[:] = varvalue

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_file', help='pixel cloud file')
    parser.add_argument('gdem_file', help='GDEM file')
    parser.add_argument('fake_pixc_file', help='Output fake pixc file')
    parser.add_argument('--subsample-factor', default=2, type=int)
    args = parser.parse_args()

    # Fake it!
    fake_pixc_from_gdem(
        args.gdem_file, args.pixc_file, args.fake_pixc_file,
        args.subsample_factor)

if __name__ == "__main__":
    main()
