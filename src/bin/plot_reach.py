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
import pandas as pd
import matplotlib.axes
import matplotlib.pyplot as plt
import SWOTWater.products.product
import SWOTRiver.analysis.riverobs
import pdb
import statsmodels.api as sm

from netCDF4 import Dataset

from reach_comparison import *
from SWOTRiver.products.rivertile import RiverTileNodes

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
    'w': '#ffffff'
}

cmap_custom = [CUSTOM_COLORS['b'], CUSTOM_COLORS['w'],
               CUSTOM_COLORS['r']]
cmaph = matplotlib.colors.LinearSegmentedColormap.from_list(
    'bwr', cmap_custom)


def get_simple_node_id(node_id, reach_id):
    return np.floor((node_id.astype(int) - (reach_id - 1) * 1000) / 10).astype(
        int)


def plot_wse(data, truth, errors, reach_id, axis, plot_slope2=True,
             title=None, prd_heights=False, plot_bit_qual=False):
    # plots the water surface elevation (wse) for each node, for the observed
    # and truth data, and the fit for the reach
    reach_id = int(reach_id)

    data_df = pd.DataFrame.from_dict(data['nodes'].variables)
    if truth is not None:
        truth_df = pd.DataFrame.from_dict(truth['nodes'].variables)
        node_i_truth = truth.nodes['reach_id'] == reach_id
        truth_df = truth_df[node_i_truth]
        truth_df.set_index('node_id')
        node_p_dist_truth = truth.nodes['p_dist_out'][node_i_truth]
        truth_wse = truth.nodes['wse'][node_i_truth]
        reach_i_truth = truth.reaches['reach_id'] == reach_id
        truth_reach_wse = truth.reaches['wse'][reach_i_truth]
        truth_slope = truth.reaches['slope'][reach_i_truth]

    # node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
    #                         np.logical_not(data.nodes['wse'].mask))
    node_i = data.nodes['reach_id'] == reach_id
    node_id = data.nodes['node_id'][node_i]
    node_q = data.nodes['node_q'][node_i]
    node_q_b = data.nodes['node_q_b'][node_i]
    data_df = data_df[node_i]
    data_df.set_index('node_id')

    node_p_dist = data.nodes['p_dist_out'][node_i]

    wse = data.nodes['wse'][node_i]
    avg_wse = np.mean(wse)
    wse_r_u = data.nodes['wse_r_u'][node_i]

    reach_i = data.reaches['reach_id'] == reach_id
    reach_wse = data.reaches['wse'][reach_i]
    reach_slope = data.reaches['slope'][reach_i]
    reach_slope2 = data.reaches['slope2'][reach_i]
    reach_width = data.reaches['width'][reach_i]
    reach_xtrk = data.reaches['xtrk_dist'][reach_i]
    reach_xtrk = str(round(reach_xtrk[0] / 1000, 1))

    axis.errorbar(node_p_dist, wse, wse_r_u, fmt='o',
                  markersize=2, label='node wse', zorder=0)
    # plot the bad quality nodes in different colour
    sus_qual_mask = node_q == 1
    deg_qual_mask = node_q == 2
    bad_qual_mask = node_q == 3
    axis.errorbar(node_p_dist[bad_qual_mask], wse[bad_qual_mask],
                  wse_r_u[bad_qual_mask],
                  fmt='o',
                  color='red',
                  markersize=2,
                  markerfacecolor='red',
                  markeredgecolor='red',
                  markeredgewidth=1,
                  label='bad qual',
                  zorder=1)
    axis.errorbar(node_p_dist[deg_qual_mask], wse[deg_qual_mask],
                  wse_r_u[deg_qual_mask],
                  fmt='o',
                  color='orange',
                  markersize=2,
                  markerfacecolor='orange',
                  markeredgecolor='orange',
                  markeredgewidth=1,
                  label='degraded qual',
                  zorder=1)
    axis.errorbar(node_p_dist[sus_qual_mask], wse[sus_qual_mask],
                  wse_r_u[sus_qual_mask],
                  fmt='o',
                  color='yellow',
                  markersize=2,
                  markerfacecolor='yellow',
                  markeredgecolor='yellow',
                  markeredgewidth=1,
                  label='suspect qual',
                  zorder=1)
    if plot_bit_qual:
        for node_dist, wse, node_q in zip(node_p_dist, wse, node_q_b):
            axis.text(
                node_dist, wse + 0.5, node_q, fontsize=3, color='lightgrey')
    if truth is not None:
        axis.plot(node_p_dist_truth, truth_wse, 'kx',
                  markersize=2, label='truth', zorder=10)

    # mark outlier nodes
    wse_outlier_ind = RiverTileNodes. \
        VARIABLES['node_q_b']['flag_meanings'].split(' ').index('wse_outlier')
    wse_outlier_q_b = RiverTileNodes. \
        VARIABLES['node_q_b']['flag_masks'][wse_outlier_ind]
    outlier_qual_mask = (node_q_b & wse_outlier_q_b) == wse_outlier_q_b
    axis.plot(node_p_dist[outlier_qual_mask],
              wse[outlier_qual_mask], 'mo',
              markersize=3, label='outlier', zorder=11)

    axis2 = axis.twiny()
    node_id = node_id - node_id[0] + 11  # no reach in node_id, for readability
    axis2.plot(node_id, avg_wse * np.ones(len(node_id)))
    axis2.cla()
    axis2.xaxis.get_offset_text().set_visible(False)
    axis2.set_xlabel('node id')

    # plot the reach slope
    # reset around PRD center
    reach_center_dist = np.mean(node_p_dist)
    ss = node_p_dist - reach_center_dist
    ss_min = min(ss)
    ss_max = max(ss)
    fit_x = [ss_min, 0, ss_max] + reach_center_dist
    # get slope end-points using slope and PRD center height
    obs_fit_y = [reach_wse + ss_min * reach_slope,
                 reach_wse,
                 reach_wse + ss_max * reach_slope]
    if truth is not None:
        truth_fit_y = [truth_reach_wse + ss_min * truth_slope,
                       truth_reach_wse,
                       truth_reach_wse + ss_max * truth_slope]
        axis.plot(fit_x, truth_fit_y, '--', markersize=10,
                  color='r', label='truth fit')
        axis.plot(np.mean(node_p_dist_truth), truth_reach_wse,
                  'r*', markersize=5, label='truth wse', zorder=0)

    axis.plot(fit_x, obs_fit_y, '--', markersize=10,
              color='b', label='obs fit')
    axis.plot(np.mean(node_p_dist), reach_wse, 'b*', markersize=5,
              color='g', label='obs wse', zorder=1)
    axis.axvline(x=reach_center_dist, ls='--', lw=0.2)
    # plot the wse_r_u shading
    axis.fill_between(node_p_dist, wse + 3 * wse_r_u, wse - 3 * wse_r_u,
                      facecolor='gray', alpha=0.3, interpolate=True)
    if plot_slope2:
        # plot enhanced reach slope
        obs_fit_y2 = [reach_wse + ss_min * reach_slope2, reach_wse,
                      reach_wse + ss_max * reach_slope2]
        axis.plot(fit_x, obs_fit_y2, '--', markersize=10, color='g',
                  label='slp2 fit')
    if prd_heights:
        axis.plot(node_p_dist, data.nodes['p_wse'][node_i],
                  'D', markersize=2, label='PRD wse')

    # Add summary metrics in text
    left, width = .05, .5
    bottom, height = .02, .82
    top = bottom + height

    if truth is not None:
        str1 = 'slope_e=' + str(round(errors['slp e (cm/km)'], 2)) + 'cm/km\n'
        axis.text(left + 0.1, top, str1,
                  horizontalalignment='left',
                  verticalalignment='bottom',
                  fontsize=5,
                  color=get_passfail_color(errors['slp e (cm/km)'],
                                           'slp e (cm/km)'),
                  transform=axis.transAxes)
        str2 = 'wse_e=' + str(round(errors['wse e (cm)'], 2)) + ' cm\n'
        axis.text(left + 0.1, top - 0.1, str2,
                  horizontalalignment='left',
                  verticalalignment='bottom',
                  fontsize=5,
                  color=get_passfail_color(errors['wse e (cm)'], 'wse e (cm)'),
                  transform=axis.transAxes)
    summary_string = 'w = ' + str(reach_width[0]) + ' m\n' \
                     + 'x-trk =' + str(reach_xtrk) + ' km'
    axis.text(left, bottom, summary_string,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=5,
              transform=axis.transAxes)

    axis.grid()
    axis.set_xlabel('dist from outlet (m)')
    axis.set_ylabel('wse (m)')
    leg = axis.legend(
        bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize=5, ncol=1)
    leg.set_draggable(1)

    if title is not None:
        axis.set_title(title)


def plot_area(data, truth, errors, reach_id, axis, title=None, style='.'):
    # plot the truth and observed area, for detected and total
    node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
                            np.logical_not(data.nodes['area_total'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_id = get_simple_node_id(node_id, reach_id)
    if truth is not None:
        node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                                      np.logical_not(truth.nodes['wse'].mask))
        node_id_truth = truth.nodes['node_id'][node_i_truth]
        node_id_truth = get_simple_node_id(node_id_truth, reach_id)
        area_truth = truth.nodes['area_total'][node_i_truth]
        axis.plot(node_id_truth, area_truth, 'kx', markersize=2)

    area_detct = data.nodes['area_detct'][node_i]
    area_total = data.nodes['area_total'][node_i]  # includes dark water pixels

    axis.plot(node_id, area_detct, style, markersize=4, alpha=.5)
    axis.plot(node_id, area_total, style, markersize=4, alpha=.5)

    # add text with error summary
    left, width = .05, .5
    bottom, height = .1, .8
    top = bottom + height
    if truth is not None:
        str1 = 'Area detect e=' + str(round(errors['area_det e (%)'], 1)) + '%\n'
        str2 = 'Area total e=' + str(round(errors['area_tot e (%)'], 1)) + '%'
        str3 = 'Width e=' + str(round(errors['width e (m)'], 1)) + ' m'
        axis.text(left, top, str1,
                  horizontalalignment='left',
                  verticalalignment='top',
                  fontsize=5,
                  transform=axis.transAxes)
        axis.text(left, top - 0.06, str2,
                  horizontalalignment='left',
                  verticalalignment='top',
                  fontsize=5,
                  color=get_passfail_color(errors['area_tot e (%)'],
                                           'area_tot e (%)'),
                  transform=axis.transAxes)
        axis.text(left, top - 0.12, str3,
                  horizontalalignment='left',
                  verticalalignment='top',
                  fontsize=5,
                  transform=axis.transAxes)

    axis.grid()
    axis.set_xlabel('node_id')
    axis.set_ylabel('area (m^2)')
    leg = axis.legend(['area detected', 'area total', 'truth'], fontsize=5)
    leg.set_draggable(1)
    if title is not None:
        axis.set_title(title)


def plot_pix_assgn(data, reach_id, axis, style='.'):
    # plot the pixel assignment

    pix_i = (data['reach_id'] == reach_id)
    node_id = data['node_id'][pix_i]
    reach_i = data['reach_id'] == reach_id

    lat = data['latitude_vectorproc'][pix_i]
    lon = data['longitude_vectorproc'][pix_i]

    plot = axis.scatter(lon, lat, cmap=plt.cm.get_cmap(
        'tab20b', len(lon)), s=2, c=node_id, edgecolor='none')
    axis.grid()
    axis.set_aspect('equal', adjustable='box')
    axis.set_xlabel('lon')
    axis.set_ylabel('lat')
    # axis.legend(['pix location'])
    axis.set_title('pixel locations')
    colorbar = plt.colorbar(plot, ax=axis)
    colorbar.set_label('node_id')


def plot_locations(data, truth, reach_id, axis, plot_prior=True,
                   gdem_dem_file=None, title=None):
    # creates plot with the observation centroids and the prior node locations
    reach_id = int(reach_id)
    node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    if truth is not None:
        node_i_truth = np.logical_and(truth.nodes['reach_id'] == reach_id,
                                      np.logical_not(truth.nodes['wse'].mask))
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
        lon_bounds = [np.max((np.min(gdem_dem_lon),
                              np.min(lon) - (np.max(lon) - np.min(lon)))),
                      np.min((np.max(gdem_dem_lon),
                              np.max(lon) + (np.max(lon) - np.min(lon))))]
        lat_bounds = [np.max((np.min(gdem_dem_lat),
                              np.min(lat) - (np.max(lat) - np.min(lat)))),
                      np.min((np.max(gdem_dem_lat),
                              np.max(lat) + (np.max(lat) - np.min(lat))))]
        lon_mask = np.logical_and(gdem_dem_lon >= lon_bounds[0],
                                  gdem_dem_lon <= lon_bounds[1])
        lat_mask = np.logical_and(gdem_dem_lat >= lat_bounds[0],
                                  gdem_dem_lat <= lat_bounds[1])
        gdem_dem_el = gdem_dem_el[lat_mask][:, lon_mask]
        cmap = [CUSTOM_COLORS['c'], CUSTOM_COLORS['m'],
                CUSTOM_COLORS['y'], CUSTOM_COLORS['c']]
        color_map = matplotlib.colors.LinearSegmentedColormap.from_list(
            'bwr', cmap)
        axis.imshow(gdem_dem_el, origin='lower', cmap=color_map,
                    extent=[np.min(gdem_dem_lon[lon_mask]),
                            np.max(gdem_dem_lon[lon_mask]),
                            np.min(gdem_dem_lat[lat_mask]),
                            np.max(gdem_dem_lat[lat_mask])])
    plot = axis.scatter(lon, lat, cmap=plt.cm.get_cmap('tab20b', len(lon)),
                        s=50, c=node_id, edgecolor='none')
    if plot_prior and truth is not None:
        axis.scatter(truth.nodes['lon_prior'][node_i_truth],
                     truth.nodes['lat_prior'][node_i_truth],
                     marker='x', s=5, c='k')
    colorbar = plt.colorbar(plot, ax=axis)
    colorbar.set_label('node_id')
    axis.grid()
    axis.set_xlabel('longitude')
    axis.set_ylabel('latitude')
    if title is not None:
        axis.set_title(title)


def get_passfail_color(error_value, parameter):
    # returns a colour that signifies how a number relates to the scientific
    # requirements for SWOT
    passfail = SWOTRiver.analysis.riverobs.get_passfail()
    if abs(error_value) < passfail[parameter][0] \
            and abs(error_value) < passfail[parameter][1]:
        return 'green'
    elif abs(error_value) > passfail[parameter][0] \
            and abs(error_value) < passfail[parameter][1]:
        return 'orange'
    else:
        return 'red'


def make_plots(rivertile_file, truth_file, pixcvec, pixc,
               truth_pixcvec, truth_pixc, reach_id,
               gdem_dem_file, errors=None, scene=None, nodes=None,
               pixc_truth=None):
    reach_id = int(reach_id)

    # create MutableProduct objects using input files
    # contains node group and reach group for each input netcdf

    data = SWOTWater.products.product.MutableProduct.from_ncfile(
        rivertile_file)
    if truth_file is not None:
        truth = SWOTWater.products.product.MutableProduct.from_ncfile(
            truth_file)
    else:
        truth = None

    if scene is not None:
        title_str = 'Scene: ' + str(scene) + ' Reach: ' + str(reach_id)
    else:
        title_str = 'Reach: ' + str(reach_id)

    figure, axes = plt.subplots(2, 2, figsize=FIGSIZE, dpi=DPI)
    plot_wse(data, truth, errors, reach_id, axes[0][0],
             title=title_str + ' - wse')
    plot_area(data, truth, errors, reach_id, axes[1][0],
              title=title_str + ' - area')
    plot_locations(data, truth, reach_id, axes[0][1],
                   gdem_dem_file=gdem_dem_file,
                   title=title_str + ' - locations')
    if pixcvec is not None:
        pixcvec_data = SWOTWater.products.product.MutableProduct.from_ncfile(
            pixcvec)
        plot_pix_assgn(pixcvec_data, reach_id, axes[1][1])
    else:
        pixcvec_data = None

    plt.tight_layout()
    mngr = plt.get_current_fig_manager()
    # mngr.window.setGeometry(0, 0, 1500, 500)
    if pixc and pixcvec_data:
        pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc)
        pixc_truth_data = None
        if pixc_truth is not None:
            pixc_truth_data = SWOTWater.products.product.MutableProduct.from_ncfile(
                pixc_truth)
        plot_pixcs(pixcvec_data, pixc_data, reach_id, nodes,
                   reach_data=data, pixc_truth=pixc_truth_data)
    else:
        print('Missing pixc or pixcvec file, skipping pixel assignment plot')

    if pixc and truth_pixc:  # only plot these if pixc was also given
        truth_pixcvec_data = SWOTWater.products.product.MutableProduct.from_ncfile(
            truth_pixcvec)
        truth_pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(
            truth_pixc)
        plot_pixcs(truth_pixcvec_data, truth_pixc_data, reach_id, nodes,
                   title_tag='(truth)', reach_data=truth)
    return figure, axes


def get_reach_error(errors, reach_id):
    # this gets the slope, wse, and area errors for the reach of interest
    reach_error = {}
    index = errors[0]['reach'].index(str(reach_id))
    for key in errors[0].keys():
        reach_error[key] = errors[0][key][index]

    return reach_error


def plot_pixcs(pixc_vec, pixc, reach_id, nodes=None,
               title_tag='(slant-plane)', reach_data=None, pixc_truth=None,
               apply_corr=True):
    reach_id = int(reach_id)
    # get only the reach_id for pixels in pixc_vec
    pix_i = (pixc_vec['reach_id'] == reach_id)
    node_id0 = pixc_vec['node_id'][pix_i]
    node_id = get_simple_node_id(node_id0, reach_id)
    aziv = pixc_vec['azimuth_index'][pix_i]
    riv = pixc_vec['range_index'][pix_i]
    heightv = pixc_vec['height_vectorproc'][pix_i]

    # map to slant_plane
    M1 = np.max(aziv) + 1
    N1 = np.max(riv) + 1
    M0 = np.min(aziv)
    N0 = np.min(riv)
    M = M1 - M0
    N = N1 - N0
    Node_id = np.zeros((M, N)) + np.nan
    Node_id[aziv - M0, riv - N0] = node_id[:]
    Heightv = np.zeros((M, N)) + np.nan
    Heightv[aziv - M0, riv - N0] = heightv[:]

    # now get PIXC in slant-plane
    azi = pixc.pixel_cloud['azimuth_index']
    ri = pixc.pixel_cloud['range_index']
    height = pixc.pixel_cloud['height']
    geoid = pixc.pixel_cloud['geoid']
    solid_tide = pixc.pixel_cloud['solid_earth_tide']
    load_tide_fes = pixc.pixel_cloud['load_tide_fes']
    load_tide_got = pixc.pixel_cloud['load_tide_got']
    pole_tide = pixc.pixel_cloud['pole_tide']
    cls = pixc.pixel_cloud['classification']
    wf = pixc.pixel_cloud['water_frac']
    pxarea = pixc.pixel_cloud['pixel_area']
    cls_t = None
    if pixc_truth is not None:
        cls_t = pixc_truth['classification'][M0:M1, N0:N1]
    m = np.max(azi) + 1
    n = np.max(ri) + 1
    Height = np.zeros((m, n)) + np.nan
    Geoid = np.zeros((m, n)) + np.nan
    Solid_tide = np.zeros((m, n)) + np.nan
    Load_tide_fes = np.zeros((m, n)) + np.nan
    Load_tide_got = np.zeros((m, n)) + np.nan
    Pole_tide = np.zeros((m, n)) + np.nan
    Cls = np.zeros((m, n)) + np.nan
    Wf = np.zeros((m, n)) + np.nan
    Pxarea = np.zeros((m, n)) + np.nan

    Height[azi, ri] = height[:]
    Geoid[azi, ri] = geoid[:]
    Solid_tide[azi, ri] = solid_tide[:]
    Load_tide_fes[azi, ri] = load_tide_fes[:]
    Load_tide_got[azi, ri] = load_tide_got[:]
    Pole_tide[azi, ri] = pole_tide[:]
    Cls[azi, ri] = cls[:]
    Wf[azi, ri] = wf[:]
    Pxarea[azi, ri] = pxarea[:]

    # now crop it to pixcvec size
    Height1 = Height[M0:M1, N0:N1]
    Geoid1 = Geoid[M0:M1, N0:N1]
    Solid_tide1 = Solid_tide[M0:M1, N0:N1]
    Load_tide_fes1 = Load_tide_fes[M0:M1, N0:N1]
    Load_tide_got1 = Load_tide_got[M0:M1, N0:N1]
    Pole_tide1 = Pole_tide[M0:M1, N0:N1]
    Cls1 = Cls[M0:M1, N0:N1]
    Wf1 = Wf[M0:M1, N0:N1]
    Pxarea1 = Pxarea[M0:M1, N0:N1]

    # exclude non-pixcvec things in this reach
    Height1[np.isnan(Heightv)] = np.nan
    Geoid1[np.isnan(Heightv)] = np.nan
    Solid_tide1[np.isnan(Heightv)] = np.nan
    Load_tide_fes1[np.isnan(Heightv)] = np.nan
    Load_tide_got1[np.isnan(Heightv)] = np.nan
    Pole_tide1[np.isnan(Heightv)] = np.nan
    Cls1[np.isnan(Node_id)] = np.nan
    Wf1[np.isnan(Node_id)] = np.nan
    Pxarea1[np.isnan(Node_id)] = np.nan
    Warea1 = Pxarea1.copy()
    Warea1[Cls1 == 2] = Pxarea1[Cls1 == 2] * 0
    Warea1[Cls1 == 1] = Pxarea1[Cls1 == 1] * 0
    Warea1[Cls1 == 3] = Pxarea1[Cls1 == 3] * 0

    if apply_corr:
        Height1 -= (Geoid1 + Solid_tide1 + Load_tide_fes1 + Pole_tide1)
        Heightv -= (Geoid1 + Solid_tide1 + Load_tide_fes1 + Pole_tide1)
        # Height1 -= (Geoid1 + Load_tide_fes1 + Pole_tide1)
        # Heightv -= (Geoid1 + Load_tide_fes1 + Pole_tide1)
    # now plot them
    cmap_max = np.nanpercentile(Height1, 80)
    cmap_min = np.nanpercentile(Height1, 20)

    plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax1 = plt.subplot(2, 3, 1)
    pt1 = ax1.imshow(Node_id, interpolation='none', aspect='auto',
                     cmap=plt.cm.get_cmap('tab20b'))
    plt.colorbar(pt1, ax=ax1)
    ax1.set_title('node_id ' + title_tag)

    # TODO: make a better cmap for classification, also make font bigger 
    ax2 = plt.subplot(2, 3, 2, sharex=ax1, sharey=ax1)
    pt2 = ax2.imshow(Cls1, interpolation='none', aspect='auto',
                     cmap='tab10', clim=(0, 5))
    ax2.set_title('classification ' + title_tag)
    plt.colorbar(pt2, ax=ax2)
    ha = None
    hat = None
    a_cmax = 500
    a_cmin = 50
    if reach_data is not None:
        NodeArea = Node_id.copy()
        NodeArea_det = Node_id.copy()
        node_i = np.logical_and(
            reach_data.nodes['reach_id'] == reach_id,
            np.logical_not(reach_data.nodes['area_total'].mask)
        )
        node_id = reach_data.nodes['node_id'][node_i]
        area_tot = reach_data.nodes['area_total'][node_i]
        area_det = reach_data.nodes['area_detct'][node_i]
        p_width = reach_data.nodes['p_width'][node_i]
        a_cmax = np.max(p_width) * 1.1
        a_cmin = np.min(p_width) * 0.9
        for node_id1 in node_id0:
            id0 = get_simple_node_id(node_id1, reach_id)
            # print(id0, node_id1)
            NodeArea[Node_id == id0] = area_tot[node_id == node_id1]
            NodeArea_det[Node_id == id0] = area_det[node_id == node_id1]
        NodeArea = NodeArea / 200
        NodeArea_det = NodeArea_det / 200
        ax3 = plt.subplot(2, 3, 3, sharex=ax1, sharey=ax1)
        pt3 = ax3.imshow(NodeArea, interpolation='none', aspect='auto',
                         cmap='jet', clim=(a_cmin,
                                           a_cmax))  # clim=(np.nanpercentile(NodeArea,10),np.nanpercentile(NodeArea,90)))
        ax3.set_title('Node Area (m^2)' + title_tag)
        plt.colorbar(pt3, ax=ax3)
        abins0 = np.linspace(50, 450, 100)
        amsk = NodeArea > -100
        ha, abins = np.histogram(NodeArea[amsk], abins0)
        amn = np.mean(NodeArea[amsk])
        amed = np.median(NodeArea[amsk])
        asd = np.std(NodeArea[amsk])
        #
        ax6 = plt.subplot(2, 3, 6, sharex=ax1, sharey=ax1)
        pt6 = ax6.imshow(NodeArea_det, interpolation='none', aspect='auto',
                         cmap='jet', clim=(a_cmin, a_cmax))  # , clim=(c0,c1))
        ax6.set_title('Node Area det. ' + title_tag)
        plt.colorbar(pt6, ax=ax6)
    else:
        ax3 = plt.subplot(2, 3, 3, sharex=ax1, sharey=ax1)
        pt3 = ax3.imshow(Warea1, interpolation='none', aspect='auto',
                         cmap='jet', clim=(0, np.nanpercentile(Warea1, 90)))
        ax3.set_title('water area (pixel-level)' + title_tag)
        plt.colorbar(pt3, ax=ax3)
        #
        ax6 = plt.subplot(2, 3, 6, sharex=ax1, sharey=ax1)
        pt6 = ax6.imshow(Geoid1, interpolation='none', aspect='auto',
                         cmap=cmaph)  # , clim=(c0,c1))
        ax6.set_title('geoid height (m) ' + title_tag)
        plt.colorbar(pt6, ax=ax6)

    ax4 = plt.subplot(2, 3, 4, sharex=ax1, sharey=ax1)
    pt4 = ax4.imshow(Heightv, interpolation='none', aspect='auto',
                     cmap=cmaph, clim=(cmap_min, cmap_max))
    ax4.set_title('height_vectorproc (m) ' + title_tag)
    plt.colorbar(pt4, ax=ax4)

    ax5 = plt.subplot(2, 3, 5, sharex=ax1, sharey=ax1)
    pt5 = ax5.imshow(Height1, interpolation='none', aspect='auto',
                     cmap=cmaph, clim=(cmap_min, cmap_max))
    ax5.set_title('height (m) ' + title_tag)
    plt.colorbar(pt5, ax=ax5)

    if cls_t is not None:
        # plot an extra set of figures for truth classification
        plt.figure(figsize=FIGSIZE, dpi=DPI)
        ax_1 = plt.subplot(2, 3, 1)
        pt_1 = ax_1.imshow(cls_t, interpolation='none', aspect='auto',
                           cmap=plt.cm.get_cmap('tab10'), clim=(0, 5))
        plt.colorbar(pt_1, ax=ax_1)
        ax_1.set_title('classification pixc_true' + title_tag)
        # map the classification to the pixcvec
        Cls_t = np.zeros_like(cls_t) + np.nan
        Area_t = np.zeros_like(cls_t) + np.nan
        Nid = np.unique(Node_id[Node_id > -1])
        print(Nid)
        clsw_t = np.zeros_like(cls_t)
        clsw_t[cls_t == 4] = 1
        clsw_t[cls_t == 3] = 1
        clsw_t[cls_t == 5] = 1

        for nid in Nid:
            # print(cls_t[Node_id==nid])
            Cls_t[Node_id == nid] = cls_t[Node_id == nid]
            Area_t[Node_id == nid] = np.nansum(
                Pxarea1[Node_id == nid] * clsw_t[Node_id == nid])
        Cls_t[Cls_t == 0] = np.nan
        Area_t = Area_t / 200
        ax_2 = plt.subplot(2, 3, 2, sharex=ax_1, sharey=ax_1)
        pt_2 = ax_2.imshow(Cls_t, interpolation='none', aspect='auto',
                           cmap=plt.cm.get_cmap('tab10'), clim=(0, 5))
        plt.colorbar(pt_2, ax=ax_2)
        ax_2.set_title('classification pixc_true' + title_tag)
        #
        Cls_t2 = np.zeros(np.shape(Cls_t))
        Cls_t2[np.logical_and(Cls_t > 0, Cls1 > 0)] = 3
        Cls_t2[np.logical_and(Cls_t > 0, np.isnan(Cls1))] = 2
        Cls_t2[np.logical_and(Cls1 > 0, np.isnan(Cls_t))] = 1
        ax_3 = plt.subplot(2, 3, 3, sharex=ax_1, sharey=ax_1)
        pt_3 = ax_3.imshow(Cls_t2, interpolation='none', aspect='auto',
                           cmap=plt.cm.get_cmap('tab10'), clim=(0, 5))
        plt.colorbar(pt_3, ax=ax_3)
        ax_3.set_title('classification ' + title_tag)
        #
        ax_4 = plt.subplot(2, 3, 4, sharex=ax_1, sharey=ax_1)
        pt_4 = ax_4.imshow(Cls1, interpolation='none', aspect='auto',
                           cmap=plt.cm.get_cmap('tab10'), clim=(0, 5))
        plt.colorbar(pt_4, ax=ax_4)
        ax_4.set_title('classification diff ' + title_tag)
        #
        ax_5 = plt.subplot(2, 3, 5, sharex=ax_1, sharey=ax_1)
        pt_5 = ax_5.imshow(Area_t, interpolation='none', aspect='auto',
                           cmap='jet', clim=(a_cmin, a_cmax))
        plt.colorbar(pt_5, ax=ax_5)
        ax_5.set_title('NodeArea pixc_true ' + title_tag)
        #
        ax_6 = plt.subplot(2, 3, 6, sharex=ax_1, sharey=ax_1)
        pt_6 = ax_6.imshow((NodeArea - Area_t) / Area_t * 100,
                           interpolation='none',
                           aspect='auto', cmap='jet')
        plt.colorbar(pt_6, ax=ax_6)
        ax_6.set_title('node area % error ' + title_tag)
        if ha is not None:
            amsk = Area_t > -100
            hat, abinst = np.histogram(Area_t[amsk], abins0)
            amnt = np.mean(Area_t[amsk])
            amedt = np.median(Area_t[amsk])
            asdt = np.std(Area_t[amsk])

    # plot area histograms
    if ha is not None:
        plt.figure()
        plt.plot(abins[:-1] + (abins[1] - abins[0]) / 2, ha)
        if hat is not None:
            plt.plot(abinst[:-1] + (abinst[1] - abinst[0]) / 2, hat)
            plt.title(
                'mean=%3.2f:%3.2f, med=%3.2f:%3.2f, std=%3.2f:%3.2f' % (
                    amn, amnt, amed, amedt, asd, asdt))
        else:
            plt.title('mean=%3.2f, med=%3.2f, std=%3.2f' % (amn, amed, asd))
        plt.xlabel('Node Area (m^2)')

    if nodes:
        for node in nodes:
            # plot node-level pixc height histograms
            idx = (Node_id == int(node))
            hgt = Height1[idx]
            hgtv = Heightv[idx]
            klass = Cls1[idx]
            # print('hgt:',hgt)
            # print('hgtv:',hgtv)
            hgt_both = np.concatenate((hgt, hgtv))
            b1 = np.nanpercentile(hgt_both, 99)
            b0 = np.nanpercentile(hgt_both, 1)
            num = 200
            if len(hgt) < 100:
                num = len(hgt) / 2 + 1
            bins = np.linspace(b0, b1, int(num))
            h, bins0 = np.histogram(hgt, bins)
            hv, bins0 = np.histogram(hgtv, bins)
            h4, bins0 = np.histogram(hgt[klass == 4], bins)
            h3, bins0 = np.histogram(hgt[klass == 3], bins)
            h2, bins0 = np.histogram(hgt[klass == 2], bins)
            hd, bins0 = np.histogram(hgt[klass > 4], bins)
            binc = bins[0:-1] + (bins[1] - bins[2]) / 2.0
            mn = np.mean(hgt)
            sd = np.std(hgt)
            plt.figure(figsize=(3, 2), dpi=DPI)
            plt.plot(binc, h)  # , linewidth=2)
            plt.plot(binc, hv)  # , linewidth=2)
            plt.plot(binc, h4)  # , linewidth=2)
            plt.plot(binc, h3)  # , linewidth=2)
            plt.plot(binc, h2, '--')  # , linewidth=2)
            plt.plot(binc, hd, ':')  # , linewidth=2)
            if reach_data is not None:
                ar = np.nanmedian(NodeArea[idx])
            else:
                ar = np.nansum(Warea1[idx])
            plt.title('reach %d, node %d, mean=%3.2f, std=%3.2f, area=%3.2f' %
                      (int(reach_id), int(node), mn, sd, ar))
            plt.xlabel('height (m)')
            plt.grid()
            plt.legend(['pixc', 'pixc_vec', 'pixc interior water',
                        'pixc edge water', 'pixc edge land',
                        'pixc dark water'], loc='best')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_tile', help='river_data/rivertile.nc')
    parser.add_argument('reach_id', help='reach id', type=int)
    parser.add_argument('--truth_tile', help='river_data/rivertile.nc',
                        default=None)
    parser.add_argument('--pixcvec',
                        help='pixcvec.nc, defaults to river_data/pixcvec.nc',
                        default=None)
    parser.add_argument('--pixc', help='pixel_cloud.nc', default=None)
    parser.add_argument('--pixc_truth', help='pixel_cloud.nc (2D-grid truth)',
                        default=None)
    parser.add_argument('--truth_pixcvec', default=None,
                        help='river_truth*/river_data/pixcvec.nc, defaults to river_truth*/river_data/pixcvec.nc')
    parser.add_argument('--truth_pixc', help='gdem_pixc.nc', default=None)
    parser.add_argument('--nodes', nargs='*',
                        help='list of nodes for which to plot height  histograms',
                        default=None)
    args = parser.parse_args()

    proc_tile = os.path.abspath(args.proc_tile)
    if args.truth_tile is not None:  # user wants to plot truth as well
        if os.path.isfile(args.truth_tile):
            truth_tile = os.path.abspath(args.truth_tile)
            gdem_dem = get_gdem_from_rivertile(args.proc_tile)
            truth_pixcvec = args.truth_pixcvec
            errors = get_errors(proc_tile, truth_tile, test=False,
                                truth_filter=None)
            reach_error = get_reach_error(errors, args.reach_id)
        else:
            print('Input truth file is not a file. Check directory names.')
    else:
        truth_tile = None
        truth_pixcvec = None
        truth_pixc = None
        gdem_dem = None
        reach_error = None
    pixcvec = args.pixcvec
    if args.pixc is None:
        pixc = None
    else:
        pixc = os.path.abspath(args.pixc)
    if args.truth_pixc is not None:
        truth_pixc = os.path.abspath(args.truth_pixc)
    pixc_truth = None
    if args.pixc_truth is not None:
        pixc_truth = os.path.abspath(args.pixc_truth)
    if os.path.isfile(proc_tile):
        make_plots(proc_tile, truth_tile, pixcvec, pixc,
                   truth_pixcvec, truth_pixc, args.reach_id,
                   gdem_dem, reach_error, nodes=args.nodes,
                   pixc_truth=pixc_truth)
        plt.show()
    else:
        print('Input file', proc_tile, 'does not exist')


if __name__ == "__main__":
    main()
