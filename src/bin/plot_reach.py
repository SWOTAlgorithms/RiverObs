#!/usr/bin/env python
'''
Copyright (c) 2020-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Plots summary data from the rivertiles in a series of plots for error characterization.

Author(s): Alexander Corben, Cassie Stuurman
'''
import os
import re
import warnings
import argparse
import numpy as np
import matplotlib.axes
import matplotlib.pyplot as plt
import SWOTWater.products.product
import SWOTRiver.analysis.riverobs
import pdb

from netCDF4 import Dataset

from get_reach_metrics import *
from reach_comparison import *

FIGSIZE = (6, 3)
DPI = 200

matplotlib.rcParams.update({'font.size': 6})

CUSTOM_COLORS = {
    'r': '#ff0000',
    'g': '#00ff00',
    'b': '#0000ff',
    'c': '#00ffff',
    'm': '#ff00ff',
    'y': '#ffff00',
}

def plot_wse(data, truth, errors, reach_id, axis, reach_fit=True, title=None, style='.'):
    # plots the water surface elevation (wse) for each node, for the observed and truth data, and the fit for the reach
    node_i = np.logical_and(data.nodes['reach_id']==reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                            np.logical_not(truth.nodes['wse'].mask))
    node_id_truth = truth.nodes['node_id'][node_i_truth]

    wse = data.nodes['wse'][node_i]
    truth_wse = truth.nodes['wse'][node_i_truth]
    #wse_u = data.nodes['wse_u'][node_i] # not currently populated in rivertiles
    wse_r_u = data.nodes['wse_r_u'][node_i]
    truth_wse_r_u = truth.nodes['wse_r_u'][node_i_truth]

    reach_i = data.reaches['reach_id']==reach_id
    reach_i_truth = truth.reaches['reach_id'] == reach_id
    reach_wse = data.reaches['wse'][reach_i]
    truth_reach_wse = truth.reaches['wse'][reach_i_truth]
    reach_slope = data.reaches['slope'][reach_i]
    truth_slope = truth.reaches['slope'][reach_i_truth]
    print('data reaches is', data.reaches)
    reach_width = data.reaches['width'][reach_i]
    truth_width = truth.reaches['width'][reach_i_truth]
    reach_width = str(round(reach_width[0],1))
    reach_xtrk = data.reaches['xtrk_dist'][reach_i]
    reach_xtrk = str(round(reach_xtrk[0]/1000,1))#.encode('utf-8','ignore')
    reach_lat = data.reaches['p_lat'][reach_i]
    reach_long = data.reaches['p_lon'][reach_i]

    axis.errorbar(node_id, wse, wse_r_u, fmt='o', markersize=4, label='pixc', zorder=0)
    axis.plot(node_id_truth, truth_wse, 'kx', markersize=2, label='truth', zorder=10)
    if reach_fit:
        mid_node = int(min(node_id) + (max(node_id) - min(node_id))/2)
        node_spacing = 200 #hard code for now
        data_fit = reach_wse - reach_slope*(node_id-mid_node)*node_spacing
        axis.plot(node_id, data_fit, '--', markersize=10, color='b', label='RiverObs fit')

    # Add summary metrics in text
    left, width = .05, .5
    bottom, height = .02, .8
    top = bottom + height
    quad_coeff, lin_coeff = get_quadratic_fit(node_id, wse)
    lake_proximity = str(get_lake_proximity(reach_id, reach_lat, reach_long))
    if errors:
        str1 = 'Slope Error=' + str(round(errors[0],2)) + 'cm/km\n'
        axis.text(left + 0.1, top, str1,
                  horizontalalignment='left',
                  verticalalignment='bottom',
                  fontsize=5,
                  color=get_passfail_color(errors[0], 'slp e (cm/km)'),
                  transform=axis.transAxes)
    summary_string = 'Lake proximity= ' + lake_proximity + ' m\n' + 'Reach Width= ' + reach_width + ' m\n' + 'Quadratic coeff=' + str(quad_coeff) + '\n' + 'Linear coeff=' + str(lin_coeff) + '\n' + 'Reach X-track distance =' + reach_xtrk + ' km'
    #'Flow angle= \n ' + \
    axis.text(left, bottom, summary_string,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=5,
              transform = axis.transAxes)

    axis.grid()
    axis.set_xlabel('node_id')
    axis.set_ylabel('wse (cm)')
    leg = axis.legend()
    if leg:
        leg.set_draggable(1)
    if title is not None:
        axis.set_title(title)

def plot_area(data, truth, errors, reach_id, axis, title=None, style='.'):
    # plot the truth and observed area, for detected and total
    node_i = np.logical_and(data.nodes['reach_id']==reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                            np.logical_not(truth.nodes['wse'].mask))
    node_id_truth = truth.nodes['node_id'][node_i_truth]

    area_detct = data.nodes['area_detct'][node_i]
    area_total = data.nodes['area_total'][node_i] # includes dark water flag
    area_truth = truth.nodes['area_detct'][node_i_truth]

    reach_i = data.reaches['reach_id']==reach_id
    reach_area_detct = data.reaches['area_detct'][reach_i]
    reach_area_total = data.reaches['area_total'][reach_i]

    axis.plot(node_id, area_detct, style, markersize=4, alpha=.5)
    axis.plot(node_id, area_total, style, markersize=4, alpha=.5)
    axis.plot(node_id_truth, area_truth, 'kx', markersize=2)

    # add text with error summary
    left, width = .05, .5
    bottom, height = .1, .8
    right = left + width
    top = bottom + height
    if errors:
        str1 = 'Area detect e=' + str(round(errors[3], 1)) + '%\n'
        str2 = 'Area total e=' + str(round(errors[2], 1)) + '%'
        str3 = 'Width e=' + str(round(errors[4], 1)) + ' m'
        axis.text(left, top, str1,
                     horizontalalignment='left',
                     verticalalignment='top',
                     fontsize=5,
                     transform=axis.transAxes)
        axis.text(left, top-0.06, str2,
                     horizontalalignment='left',
                     verticalalignment='top',
                     fontsize=5,
                     color=get_passfail_color(errors[2], 'area_tot e (%)'),
                     transform=axis.transAxes)
        axis.text(left, top - 0.12, str3,
                  horizontalalignment='left',
                  verticalalignment='top',
                  fontsize=5,
                  transform=axis.transAxes)

    axis.grid()
    axis.set_xlabel('node_id')
    axis.set_ylabel('area')
    leg = axis.legend(['area detected', 'area total', 'truth'])
    if leg:
        leg.set_draggable(1)
    if title is not None:
        axis.set_title(title)

def plot_pix_assgn(data, reach_id, axis, style='.'):
    # plot the pixel assignment
    # to do: find a way to either loop colourmap or find one with one colour per node
    pix_i = (data['reach_id'] == reach_id) #np.logical_and...np.logical_not(data['wse'].mask))
    node_id = data['node_id'][pix_i]
    reach_i = data['reach_id'] == reach_id

    lat = data['latitude_vectorproc'][pix_i]
    lon = data['longitude_vectorproc'][pix_i]
    
    plot = axis.scatter(lon, lat, cmap=plt.cm.get_cmap('tab20b', len(lon)), s=2, c=node_id, edgecolor='none')
    axis.grid()
    axis.set_aspect('equal', adjustable='box')
    axis.set_xlabel('lon')
    axis.set_ylabel('lat')
    #axis.legend(['pix location'])
    axis.set_title('pixel locations')
    colorbar = plt.colorbar(plot, ax=axis)
    colorbar.set_label('node_id')

def plot_locations(data, truth, reach_id, axis, plot_prior=True, gdem_dem_file=None, title=None):
    # creates the plot with the observation centroids and the prior node locations
    node_i = np.logical_and(data.nodes['reach_id']==reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                            np.logical_not(truth.nodes['wse'].mask))
    node_id_truth = truth.nodes['node_id'][node_i_truth]
    lat = data.nodes['lat'][node_i]
    lon = data.nodes['lon'][node_i]

    if gdem_dem_file is not None:
        with Dataset(gdem_dem_file, 'r') as fin:
            gdem_dem_lat = fin['latitude'][:]
            gdem_dem_lon = fin['longitude'][:]
            gdem_dem_el = fin['elevation'][:]
        lon_bounds = [np.max((np.min(gdem_dem_lon), np.min(lon)-(np.max(lon)-np.min(lon)))),
                      np.min((np.max(gdem_dem_lon), np.max(lon)+(np.max(lon)-np.min(lon))))]
        lat_bounds = [np.max((np.min(gdem_dem_lat), np.min(lat)-(np.max(lat)-np.min(lat)))),
                      np.min((np.max(gdem_dem_lat), np.max(lat)+(np.max(lat)-np.min(lat))))]
        lon_mask = np.logical_and(gdem_dem_lon >= lon_bounds[0],
                                  gdem_dem_lon <= lon_bounds[1])
        lat_mask = np.logical_and(gdem_dem_lat >= lat_bounds[0],
                                  gdem_dem_lat <= lat_bounds[1])
        gdem_dem_el = gdem_dem_el[lat_mask][:,lon_mask]
        cmap = [CUSTOM_COLORS['c'], CUSTOM_COLORS['m'],
                CUSTOM_COLORS['y'], CUSTOM_COLORS['c']]
        color_map = matplotlib.colors.LinearSegmentedColormap.from_list(
            'cmyc', cmap)
        axis.imshow(gdem_dem_el, origin='lower', cmap=color_map,
                    extent=[np.min(gdem_dem_lon[lon_mask]),
                            np.max(gdem_dem_lon[lon_mask]),
                            np.min(gdem_dem_lat[lat_mask]),
                            np.max(gdem_dem_lat[lat_mask])])
    plot = axis.scatter(lon, lat, cmap=plt.cm.get_cmap('tab20b', len(lon)), s=50, c=node_id, edgecolor='none')
    if plot_prior:
        axis.scatter(truth.nodes['lon_prior'][node_i_truth], truth.nodes['lat_prior'][node_i_truth],
                     marker='x', s=5, c='k')
        print('prior location is', truth.nodes['lon_prior'][node_i_truth], truth.nodes['lat_prior'][node_i_truth])
    colorbar = plt.colorbar(plot, ax=axis)
    colorbar.set_label('node_id')
    #if plot_prior:
        #axis.legend(['data node'])#, 'prior node'])
    axis.grid()
    axis.set_xlabel('longitude')
    axis.set_ylabel('latitude')
    if title is not None:
        axis.set_title(title)

def get_passfail_color(error_value, parameter):
    # returns a colour that signifies how a number relates to the scientific requirements for SWOT
    passfail = SWOTRiver.analysis.riverobs.get_passfail()
    if abs(error_value) < passfail[parameter][0] and abs(error_value) < passfail[parameter][1]:
        return 'green'
    elif abs(error_value) > passfail[parameter][0] and abs(error_value) < passfail[parameter][1]:
        return 'orange'
    else:
        return 'red'

def make_plots(rivertile_file, reach_id, gdem_dem_file, errors=None, scene=None):
    # creates the figure window and populates it
    pixc_vec = rivertile_file[0:-12] + '/pixcvec.nc'
    data = SWOTWater.products.product.MutableProduct.from_ncfile(rivertile_file)
    pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc_vec)
    truth_file = rivertile_file[0:-17] + 'gdem/rivertile.nc'
    truth = SWOTWater.products.product.MutableProduct.from_ncfile(truth_file)

    if scene is not None:
        title_str = 'Scene: ' + str(scene) + ' Reach: ' + str(reach_id)
    else:
        title_str = 'Reach: ' + str(reach_id)

    figure, axes = plt.subplots(2, 2, figsize=FIGSIZE, dpi=DPI)
    plot_wse(data, truth, errors, reach_id, axes[0][0], title=title_str + ' - wse')
    plot_area(data, truth, errors, reach_id, axes[1][0], title=title_str + ' - area')
    plot_locations(data, truth, reach_id, axes[0][1], gdem_dem_file=gdem_dem_file,
                   title=title_str + ' - locations')
    plot_pix_assgn(pixc_data, reach_id, axes[1][1])

    plt.tight_layout()
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(0, 0, 1500, 500)

    return figure, axes

def get_reach_error(errors, reach_id):
    # this gets the slope, wse, and area errors for the reach of interest
    scene = 0
    for reach_index, reach in enumerate(errors[0]['reach']):
        if int(reach) == int(reach_id):
            slope_error = errors[0]['slp e (cm/km)'][reach_index]
            wse_error = errors[0]['wse e (cm)'][reach_index]
            area_error = errors[0]['area_tot e (%)'][reach_index]
            area_dtct_error = errors[0]['area_det e (%)'][reach_index]
            width_error = errors[0]['width e (m)'][reach_index]
            scene = errors[0]['scene_pass_tile'][reach_index]
    if not scene:
        raise Exception('Reach ID is not found in this scene/pass/side')
    reach_error = [0, scene, reach_id, slope_error, wse_error, area_error, area_dtct_error, width_error]
    return reach_error

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_rivertile', help='pixc/rivertile.nc')
    parser.add_argument('reach_id', help='reach id', type=int)
    args = parser.parse_args()

    gdem_dem = get_gdem_from_pixc(args.pixc_rivertile)
    gdem_tile = args.pixc_rivertile[0:-17] + 'gdem/rivertile.nc'
    errors = get_errors(args.pixc_rivertile, gdem_tile, test=False, verbose=False)
    reach_error = get_reach_error(errors, args.reach_id)
    make_plots(args.pixc_rivertile, args.reach_id, gdem_dem, reach_error)
    plt.show()

if __name__ == "__main__":
    main()
