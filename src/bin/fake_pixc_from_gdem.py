#!/usr/bin/env python
"""
fake_pixc_from_gdem.py: Makes a fake pixc from a gdem

Useage: fake_pixc_from_gdem.py pixc.nc gdem.nc fake_pixc.nc

Author (s): Alex Fore
"""
import argparse
import os
import netCDF4
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_file', help='pixel cloud file')
    parser.add_argument('gdem_file', help='GDEM file')
    parser.add_argument('fake_pixc_file', help='Output fake pixc file')
    parser.add_argument('--subsample-factor', default=2, type=int)
    args = parser.parse_args()

    subsample_factor = args.subsample_factor

    with netCDF4.Dataset(args.gdem_file, 'r') as ifp_gdem,\
         netCDF4.Dataset(args.pixc_file, 'r') as ifp_pixc,\
         netCDF4.Dataset(args.fake_pixc_file, 'w') as ofp:

        # copy tvp data
        ofp.createGroup('tvp')
        ofp.groups['tvp'].createDimension('nr_tvps', 0)

        for key, value in ifp_pixc.groups['tvp'].variables.items():
            var = ofp.createVariable('/tvp/'+key, value.dtype.str, ('nr_tvps',))
            var[:] = value[:]

        # copy pixel_cloud attributes
        ofp.createGroup('pixel_cloud')
        ofp.groups['pixel_cloud'].createDimension('record', 0)

        for attr in ifp_pixc.groups['pixel_cloud'].__dict__:
            value = ifp_pixc.groups['pixel_cloud'].__dict__[attr]
            setattr(ofp.groups['pixel_cloud'], attr, value)

        landtype = ifp_gdem.variables['landtype'][:][::subsample_factor]
        latitude = ifp_gdem.variables['latitude'][:][::subsample_factor]
        longitude = ifp_gdem.variables['longitude'][:][::subsample_factor]
        elevation = ifp_gdem.variables['elevation'][:][::subsample_factor]
        cross_track_ = ifp_gdem.variables['cross_track'][:]
        range_spacing = ifp_gdem.ground_spacing
        azimuth_spacing = ifp_gdem.azimuth_spacing

        gdem_shape = landtype.shape
        range_index, azimuth_index = np.meshgrid(
            np.arange(gdem_shape[1]), np.arange(gdem_shape[0]))
        cross_track, tmp = np.meshgrid(cross_track_, np.arange(gdem_shape[0]))

        mask = landtype > 0
        pixc_shape = range_index[mask].shape
        pixel_area = subsample_factor * range_spacing * azimuth_spacing

        out_pixc_dsets = {}
        out_pixc_dsets['range_index'] = range_index[mask]
        out_pixc_dsets['azimuth_index'] = azimuth_index[mask]
        out_pixc_dsets['classification'] = landtype[mask]
        out_pixc_dsets['continuous_classification'] = out_pixc_dsets[
            'classification'].copy()*0 + 1.0
        out_pixc_dsets['latitude'] = latitude[mask]
        out_pixc_dsets['longitude'] = longitude[mask]
        out_pixc_dsets['height'] = elevation[mask]
        out_pixc_dsets['cross_track'] = cross_track[mask]
        out_pixc_dsets['illumination_time'] = azimuth_index[mask]
        out_pixc_dsets['num_rare_looks'] = np.zeros(pixc_shape) + subsample_factor
        out_pixc_dsets['pixel_area'] =  np.zeros(pixc_shape) + pixel_area

        for varname, varvalue in out_pixc_dsets.items():
            var = ofp.createVariable(
                '/pixel_cloud/'+varname, varvalue.dtype.str, ('record',))
            var[:] = varvalue


if __name__ == "__main__":
    main()
