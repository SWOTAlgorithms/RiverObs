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
import math
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
cmap_custom = [CUSTOM_COLORS['c'], CUSTOM_COLORS['m'],
                CUSTOM_COLORS['y'], CUSTOM_COLORS['c']]
cmaph = matplotlib.colors.LinearSegmentedColormap.from_list(
            'cmyc', cmap_custom)

def plot_wse(data, truth, errors, reach_id, axis, reach_fit=True, title=None, style='.'):
    # plots the water surface elevation (wse) for each node, for the observed and truth data, and the fit for the reach
    reach_id = int(reach_id)
    node_i = np.logical_and(data.nodes['reach_id'] == reach_id, np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id, np.logical_not(truth.nodes['wse'].mask))
    node_id_truth = truth.nodes['node_id'][node_i_truth]

    node_dist = data.nodes['p_dist_out'][node_i]
    node_dist_truth = truth.nodes['p_dist_out'][node_i_truth]

    wse = data.nodes['wse'][node_i]
    truth_wse = truth.nodes['wse'][node_i_truth]
    avg_wse = np.mean(wse)
    avg_truth_wse = np.mean(truth_wse)
    # wse_u = data.nodes['wse_u'][node_i] # not currently populated in rivertiles
    wse_r_u = data.nodes['wse_r_u'][node_i]
    truth_wse_r_u = truth.nodes['wse_r_u'][node_i_truth]
    node_dist = data.nodes['p_dist_out'][node_i]
    truth_node_dist = truth.nodes['p_dist_out'][node_i_truth]

    reach_i = data.reaches['reach_id'] == reach_id
    reach_i_truth = truth.reaches['reach_id'] == reach_id
    reach_wse = data.reaches['wse'][reach_i]
    truth_reach_wse = truth.reaches['wse'][reach_i_truth]
    reach_slope = data.reaches['slope'][reach_i]
    truth_slope = truth.reaches['slope'][reach_i_truth]
    reach_width = data.reaches['width'][reach_i]
    truth_width = truth.reaches['width'][reach_i_truth]
    reach_width = str(round(reach_width[0], 1))
    reach_xtrk = data.reaches['xtrk_dist'][reach_i]
    reach_xtrk = str(round(reach_xtrk[0] / 1000, 1))  # .encode('utf-8','ignore')
    reach_lat = data.reaches['p_lat'][reach_i]
    reach_long = data.reaches['p_lon'][reach_i]
    reach_wse_r_u = data.reaches['wse_r_u'][reach_i]
    reach_wse_t_r_u = truth.reaches['wse_r_u'][reach_i]


    axis.errorbar(node_dist, wse, wse_r_u, fmt='o', markersize=4, label='pixc', zorder=0)
    axis.plot(node_dist_truth, truth_wse, 'kx', markersize=2, label='truth', zorder=10)
    axis2 = axis.twiny()
    axis2.plot(node_id, avg_wse*np.ones(len(node_id)))
    axis2.cla()
    axis2.xaxis.get_offset_text().set_visible(False)
    axis2.set_xlabel('node id')

    print('reach wse is', reach_wse,
          'truth wse is', truth_reach_wse,
          'reach slope is', reach_slope,
          'truth slope is', truth_slope)
    print('Difference of node-level obs and truth means', avg_wse - avg_truth_wse)
    print('calculated slope error is', reach_slope - truth_slope, 'calculated wse error is',
          reach_wse - truth_reach_wse)
    print('normalized error is',
          (reach_wse - truth_reach_wse) * 1e2 / math.sqrt(reach_wse_r_u ** 2 + reach_wse_t_r_u ** 2))
    print('input slope error is', str(round(errors[0], 2)), 'cm/km')
    if reach_fit:
        mid_node = int(min(node_id) + (max(node_id) - min(node_id)) / 2)
        mid_node_truth = int(min(node_id_truth) + (max(node_id_truth) - min(node_id_truth)) / 2)
        node_spacing = abs(np.max(node_dist) - np.min(node_dist))/(len(node_id)-1)
        node_spacing_truth = abs(np.max(node_dist_truth) - np.min(node_dist_truth))/(len(node_id_truth)-1)
        print('average node spacing is', node_spacing)
        print('average truth node spacing is', node_spacing_truth)
        data_fit = reach_wse - reach_slope/10 * (mid_node-node_id) * node_spacing
        truth_fit = truth_reach_wse - truth_slope/10 * (mid_node_truth-node_id_truth) * node_spacing_truth
        axis.plot(truth_node_dist, truth_fit, '--', markersize=10, color='r', label='Truth fit')
        axis.plot(node_dist, data_fit, '--', markersize=10, color='b', label='RiverObs fit')

    # Add summary metrics in text
    left, width = .05, .5
    bottom, height = .02, .85
    top = bottom + height
    #quad_coeff, lin_coeff = get_quadratic_fit(node_id_truth, truth_wse)
    try: # fix this later, only important for slope assessments
        r_sq = str(get_r_squared(node_id_truth, truth_wse))
    except:
        print('couldnt get r squared')
        r_sq = '0'
    #lake_proximity = str(get_lake_proximity(reach_id, reach_lat, reach_long))
    if errors:
        str1 = 'Slope Error=' + str(round(errors[0], 2)) + 'cm/km\n' #+ 'linear r-sq=' + r_sq
        axis.text(left + 0.1, top, str1,
                  horizontalalignment='left',
                  verticalalignment='bottom',
                  fontsize=5,
                  color=get_passfail_color(errors[0], 'slp e (cm/km)'),
                  transform=axis.transAxes)
        str2 = 'wse error=' + str(round(errors[1], 2)) + ' cm\n'
        axis.text(left + 0.1, top - 0.1, str2,
                  horizontalalignment='left',
                  verticalalignment='bottom',
                  fontsize=5,
                  color=get_passfail_color(errors[1], 'wse e (cm)'),
                  transform=axis.transAxes)
    summary_string = 'Reach Width= ' + reach_width + ' m\n' + 'Reach X-track distance =' + reach_xtrk + ' km'
    # 'Lake proximity= ' + lake_proximity + ' m\n' +
    # 'Flow angle= \n ' + \
    # + 'Truth quad coeff=' + str(
    #         quad_coeff) + '\n' + 'Truth Lin coeff=' + str(
    #         lin_coeff) +
    axis.text(left, bottom, summary_string,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=5,
              transform=axis.transAxes)

    axis.grid()
    axis.set_xlabel('dist from outlet (m)')
    axis.set_ylabel('wse (m)')
    leg = axis.legend()
    if leg:
        leg.set_draggable(1)
    if title is not None:
        axis.set_title(title)

def plot_area(data, truth, errors, reach_id, axis, title=None, style='.'):
    # plot the truth and observed area, for detected and total
    node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                                  np.logical_not(truth.nodes['wse'].mask))
    node_id_truth = truth.nodes['node_id'][node_i_truth]

    area_detct = data.nodes['area_detct'][node_i]
    area_total = data.nodes['area_total'][node_i]  # includes dark water flag
    area_truth = truth.nodes['area_detct'][node_i_truth]

    reach_i = data.reaches['reach_id'] == reach_id
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
        axis.text(left, top - 0.06, str2,
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
    pix_i = (data['reach_id'] == reach_id)  # np.logical_and...np.logical_not(data['wse'].mask))
    node_id = data['node_id'][pix_i]
    reach_i = data['reach_id'] == reach_id

    lat = data['latitude_vectorproc'][pix_i]
    lon = data['longitude_vectorproc'][pix_i]

    plot = axis.scatter(lon, lat, cmap=plt.cm.get_cmap('tab20b', len(lon)), s=2, c=node_id, edgecolor='none')
    axis.grid()
    axis.set_aspect('equal', adjustable='box')
    axis.set_xlabel('lon')
    axis.set_ylabel('lat')
    # axis.legend(['pix location'])
    axis.set_title('pixel locations')
    colorbar = plt.colorbar(plot, ax=axis)
    colorbar.set_label('node_id')


def plot_locations(data, truth, reach_id, axis, plot_prior=True, gdem_dem_file=None, title=None):
    # creates the plot with the observation centroids and the prior node locations
    reach_id = int(reach_id)
    node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                                  np.logical_not(truth.nodes['wse'].mask))
    node_id_truth = truth.nodes['node_id'][node_i_truth]
    lat = data.nodes['lat'][node_i]
    lon = data.nodes['lon'][node_i]

    if gdem_dem_file is not None:
        try:
            with Dataset(gdem_dem_file, 'r') as fin:
                gdem_dem_lat = fin['latitude'][:]
                gdem_dem_lon = fin['longitude'][:]
                gdem_dem_el = fin['elevation'][:]
        except FileNotFoundError:
            gdem_dem_file = gdem_dem_file[0:-8] + '_sim.nc'
            if '__' in gdem_dem_file:
                gdem_dem_file = gdem_dem_file.replace('__', '_')
            print('gdem_file is', gdem_dem_file)
            with Dataset(gdem_dem_file, 'r') as fin:
                gdem_dem_lat = fin['latitude'][:]
                gdem_dem_lon = fin['longitude'][:]
                gdem_dem_el = fin['elevation'][:]
        lon_bounds = [np.max((np.min(gdem_dem_lon), np.min(lon) - (np.max(lon) - np.min(lon)))),
                      np.min((np.max(gdem_dem_lon), np.max(lon) + (np.max(lon) - np.min(lon))))]
        lat_bounds = [np.max((np.min(gdem_dem_lat), np.min(lat) - (np.max(lat) - np.min(lat)))),
                      np.min((np.max(gdem_dem_lat), np.max(lat) + (np.max(lat) - np.min(lat))))]
        lon_mask = np.logical_and(gdem_dem_lon >= lon_bounds[0],
                                  gdem_dem_lon <= lon_bounds[1])
        lat_mask = np.logical_and(gdem_dem_lat >= lat_bounds[0],
                                  gdem_dem_lat <= lat_bounds[1])
        gdem_dem_el = gdem_dem_el[lat_mask][:, lon_mask]
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
    colorbar = plt.colorbar(plot, ax=axis)
    colorbar.set_label('node_id')
    # if plot_prior:
    # axis.legend(['data node'])#, 'prior node'])
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


def make_plots(rivertile_file, truth_file, pixcvec, pixc,
        truth_pixcvec, truth_pixc, reach_id,
        gdem_dem_file, errors=None, scene=None, nodes=None):
    reach_id = int(reach_id)
    # creates the figure window and populates it
    #pixc_vec = rivertile_file[0:-12] + '/pixcvec.nc'
    data = SWOTWater.products.product.MutableProduct.from_ncfile(rivertile_file)
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
    if pixcvec is not None:
        pixcvec_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixcvec)
        plot_pix_assgn(pixcvec_data, reach_id, axes[1][1])
    else:
        pixcvec_data = None

    plt.tight_layout()
    mngr = plt.get_current_fig_manager()
    # mngr.window.setGeometry(0, 0, 1500, 500)
    if pixc and pixcvec_data:
        pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc)
        plot_pixcs(pixcvec_data, pixc_data, reach_id, nodes)
    else:
        print('Missing pixc or pixcvec file, skipping pixel assignment plot')

    if pixc and truth_pixc:# only plot these if pixc was also given
        truth_pixcvec_data = SWOTWater.products.product.MutableProduct.from_ncfile(truth_pixcvec)
        truth_pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(truth_pixc)
        plot_pixcs(truth_pixcvec_data, truth_pixc_data, reach_id, nodes, title_tag='(truth)')

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
    reach_error = [slope_error, wse_error, area_error, area_dtct_error, width_error]
    return reach_error

def plot_pixcs(pixc_vec, pixc, reach_id, nodes=None, title_tag='(slant-plane)'):
    reach_id = int(reach_id)
    # get only the reach_id for pixels in pixc_vec
    pix_i = (pixc_vec['reach_id'] == reach_id)
    node_id0 = pixc_vec['node_id'][pix_i]
    reach_i = pixc_vec['reach_id'] == reach_id
    #print('node0: ', node_id0[0])
    #print('reach_id: ', reach_id)
    node_id = node_id0.astype('int') - (reach_id-1)*1000
    #print('node[0]: ', node_id[0])

    aziv = pixc_vec['azimuth_index'][pix_i]
    riv = pixc_vec['range_index'][pix_i]

    latv = pixc_vec['latitude_vectorproc'][pix_i]
    lonv = pixc_vec['longitude_vectorproc'][pix_i]
    heightv = pixc_vec['height_vectorproc'][pix_i]

    #azi = pixc['azimuth_index']
    #ri = pixc['range_index']

    # map to slant_plane
    M1 = np.max(aziv)+1
    N1 = np.max(riv)+1
    M0 = np.min(aziv)
    N0 = np.min(riv)
    M = M1-M0
    N = N1-N0
    Node_id = np.zeros((M,N)) + np.nan
    Node_id[aziv-M0,riv-N0] = node_id[:]
    Heightv = np.zeros((M,N)) + np.nan
    Heightv[aziv-M0,riv-N0] = heightv[:]


    # now get PIXC in slant-plane
    azi = pixc.pixel_cloud['azimuth_index']
    ri = pixc.pixel_cloud['range_index']
    height = pixc.pixel_cloud['height']
    geoid = pixc.pixel_cloud['geoid']
    cls = pixc.pixel_cloud['classification']
    m = np.max(azi)+1
    n = np.max(ri)+1
    Height = np.zeros((m,n)) + np.nan
    Geoid = np.zeros((m,n)) + np.nan
    Cls = np.zeros((m,n)) + np.nan
    Height[azi,ri] = height[:]
    Geoid[azi,ri] = geoid[:]
    Cls[azi,ri] = cls[:]

    # now crop it to pixcvec size
    Height1 = Height[M0:M1,N0:N1]
    Geoid1 = Geoid[M0:M1,N0:N1]
    Cls1 = Cls[M0:M1,N0:N1]
    # exclude non-pixcvec things in theis reach
    Height1[np.isnan(Heightv)] = np.nan
    Geoid1[np.isnan(Heightv)] = np.nan
    Cls1[np.isnan(Node_id)] = np.nan
    # now plot them
    c1 = np.nanpercentile(Height1,80)
    c0 = np.nanpercentile(Height1,20)

    plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax1 = plt.subplot(2,3,1)
    pt1 =ax1.imshow(Node_id, interpolation='none',aspect='auto',
        cmap=plt.cm.get_cmap('tab20b'))
    plt.colorbar(pt1,ax=ax1)
    ax1.set_title('node_id '+title_tag)

    # TODO: make a better cmap for classification, also make font bigger 
    ax2 = plt.subplot(2,3,2, sharex=ax1, sharey=ax1)
    pt2 = ax2.imshow(Cls1, interpolation='none', aspect='auto',
        cmap='tab10', clim=(0,5))
    ax2.set_title('classification '+title_tag)
    plt.colorbar(pt2,ax=ax2)

    ax3 = plt.subplot(2,3,4, sharex=ax1, sharey=ax1)
    pt3 = ax3.imshow(Heightv, interpolation='none', aspect='auto',
        cmap=cmaph, clim=(c0,c1))
    ax3.set_title('height_vectorproc (m) '+title_tag)
    plt.colorbar(pt3,ax=ax3)

    ax4 = plt.subplot(2,3,5, sharex=ax1, sharey=ax1)
    pt4 = ax4.imshow(Height1, interpolation='none', aspect='auto',
        cmap=cmaph, clim=(c0,c1))
    ax4.set_title('height (m) '+title_tag)
    plt.colorbar(pt4,ax=ax4)

    ax5 = plt.subplot(2,3,6, sharex=ax1, sharey=ax1)
    pt5 = ax5.imshow(Geoid1, interpolation='none', aspect='auto',
        cmap=cmaph)#, clim=(c0,c1))
    ax5.set_title('geoid height (m) '+title_tag)
    plt.colorbar(pt5,ax=ax5)

    if nodes:
        for node in nodes:
            # plot node-level pixc height histograms
            idx = (Node_id==int(node))
            hgt = Height1[idx]
            hgtv = Heightv[idx]
            hgtv = Heightv[idx]
            klass = Cls1[idx]
            #print('hgt:',hgt)
            #print('hgtv:',hgtv)
            hgt_both = np.concatenate((hgt, hgtv))
            b1 = np.nanpercentile(hgt_both,99)
            b0 = np.nanpercentile(hgt_both,1)
            num = 200
            if len(hgt) < 100:
                num = len(hgt)/2 + 1
            bins = np.linspace(b0,b1, int(num))
            h, bins0 = np.histogram(hgt, bins)
            hv, bins0 = np.histogram(hgtv, bins)
            h4, bins0 = np.histogram(hgt[klass==4], bins)
            h3, bins0 = np.histogram(hgt[klass==3], bins)
            h2, bins0 = np.histogram(hgt[klass==2], bins)
            hd, bins0 = np.histogram(hgt[klass>4], bins)
            binc = bins[0:-1] + (bins[1]-bins[2])/2.0
            mn = np.mean(hgt)
            sd = np.std(hgt)
            plt.figure(figsize=(3,2), dpi=DPI)
            plt.plot(binc, h)#, linewidth=2)
            plt.plot(binc, hv)#, linewidth=2)
            plt.plot(binc, h4)#, linewidth=2)
            plt.plot(binc, h3)#, linewidth=2)
            plt.plot(binc, h2,'--')#, linewidth=2)
            plt.plot(binc, hd, ':')#, linewidth=2)
            plt.title('node %d, mean=%3.2f, std=%3.2f'%(int(node), mn, sd))
            plt.xlabel('height (m)')
            plt.grid()
            plt.legend(['pixc', 'pixc_vec',
                'pixc interior water', 'pixc edge water',
                'pixc edge land', 'pixc dark water'],
                loc='best')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_tile', help='river_data/rivertile.nc')
    parser.add_argument('truth_tile', help='river_data/rivertile.nc')
    parser.add_argument('reach_id', help='reach id', type=int)
    parser.add_argument('--pixcvec', help='pixcvec.nc, defaults to river_data/pixcvec.nc', default=None)
    parser.add_argument('--pixc', help='pixel_cloud.nc', default=None)
    parser.add_argument('--truth_pixcvec', default=None,
        help='river_truth*/river_data/pixcvec.nc, defaults to river_truth*/river_data/pixcvec.nc')
    parser.add_argument('--truth_pixc', help='gdem_pixc.nc', default=None)
    parser.add_argument('--nodes', nargs='*',
        help='list of nodes for which to plot height  histograms', default=None)
    args = parser.parse_args()

    proc_tile = os.path.abspath(args.proc_tile)
    truth_tile = os.path.abspath(args.truth_tile)
    gdem_dem = get_gdem_from_pixc(args.proc_tile)
    gdem_tile = args.truth_tile
    pixcvec = args.pixcvec
    truth_pixcvec = args.truth_pixcvec
    #proc_tile = os.path.abspath(args.proc_tile)
    #truth_tile = os.path.abspath(args.truth_tile)
    if pixcvec is None:
        pixcvec = proc_tile[0:-12] + '/pixcvec.nc'
    if truth_pixcvec is None:
        truth_pixcvec = truth_tile[0:-12] + '/pixcvec.nc'
    if args.pixc is None:
        pixc = None
    else:
        pixc = os.path.abspath(args.pixc)
    if args.truth_pixc is None:
        truth_pixc = truth_tile[0:-12] + '/gdem_pixc.nc'
    else:
        truth_pixc = os.path.abspath(args.truth_pixc)
    errors = get_errors(proc_tile, truth_tile, test=False, verbose=False)
    reach_error = get_reach_error(errors, args.reach_id)
    #make_plots(args.proc_tile, args.truth_tile, args.reach_id, gdem_dem, reach_error)
    make_plots(proc_tile, truth_tile, pixcvec, pixc,
        truth_pixcvec, truth_pixc, args.reach_id,
        gdem_dem, reach_error, nodes=args.nodes)
    plt.show()


if __name__ == "__main__":
    main()
