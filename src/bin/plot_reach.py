#!/usr/bin/env python
'''
Copyright (c) 2020-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Plots summary data from the rivertiles in a series of plots for error
characterization.

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
from matplotlib import colors
from scipy import stats
import SWOTWater.products.product
import SWOTRiver.analysis.riverobs
import statsmodels.api as sm
import seaborn as sns
import matplotlib.colors as mcolors
import geopandas as gpd

from netCDF4 import Dataset

from reach_comparison import *
from SWOTRiver.products.rivertile import RiverTileNodes

FIGSIZE = (16, 9)
DPI = 200
LEFT, WIDTH = .05, .75
RIGHT = LEFT + WIDTH
BOTTOM, HEIGHT = .02, .85
TOP = BOTTOM + HEIGHT

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
    return np.floor(
        (node_id.astype(int)
         - (reach_id - 1) * 1000) / 10).astype(int)


def assign_value_or_none(data, key):
    if np.isnan(data[key]).any():
        return None
    else:
        return get_first_element(data[key])


def get_first_element(data):
    # Returns the first element of a list or pandas series. Mainly used to grab
    # river names for plot titles/filenames.
    if isinstance(data, pd.Series) or isinstance(data, pd.DataFrame):
        # Handle Pandas Series and DataFrame
        return data.iloc[0] if isinstance(data, pd.Series) else data.iloc[0, 0]
    elif isinstance(data, np.ma.MaskedArray):
        # Handle NumPy Masked Array
        return data[0]
    else:
        # Optional: Handle other types, or raise an error if type is unexpected
        raise TypeError("Unsupported data type")


def plot_wse(data, truth, errors, reach_id, axis, figure,
             title=None, prd_heights=False, plot_bit_qual=True, cycle=None,
             tile=None, pass_no=None, annotate_metrics=True,
             multi_reach=False):
    # plots the water surface elevation (wse) for each node, for the observed
    # and truth data, and the fit for the reach. Input "data" should match
    # the format of the mutable netcdf product from SWOTWater.products.product.
    # MutableProduct.from_ncfile()

    reach_id = int(reach_id)
    if multi_reach:
        # grab adjacent reaches too if they exist
        node_i = np.logical_or.reduce((
            data.nodes['reach_id'] == reach_id,
            data.nodes['reach_id'] == reach_id - 10,
            data.nodes['reach_id'] == reach_id + 10
        ))
    else:
        node_i = data.nodes['reach_id'] == reach_id
    node_id = data.nodes['node_id'][node_i]
    node_q = data.nodes['node_q'][node_i]
    node_q_b = data.nodes['node_q_b'][node_i]
    node_p_dist = data.nodes['p_dist_out'][node_i]

    wse = data.nodes['wse'][node_i]
    if np.sum(wse > -999) == 0:
        # fill value heights only; can't plot
        has_truth = False
        return has_truth
    avg_wse = np.mean(wse)
    wse_r_u = data.nodes['wse_r_u'][node_i]
    if multi_reach:
        # grab adjacent reaches too if they exist
        reach_i = np.logical_or.reduce((
            data.reaches['reach_id'] == reach_id,
            data.reaches['reach_id'] == reach_id - 10,
            data.reaches['reach_id'] == reach_id + 10
        ))
    else:
        reach_i = data.reaches['reach_id'] == reach_id

    reach_wse = data.reaches['wse'][reach_i]
    reach_slope = data.reaches['slope'][reach_i]
    reach_slope2 = data.reaches['slope2'][reach_i]
    river_name = get_first_element(data.reaches['river_name'][reach_i])
    fit_x, ss_min, ss_max = plot_wse_and_qual(
        node_p_dist, wse, wse_r_u, node_q, node_q_b, axis, plot_bit_qual,
        reach_wse, reach_slope, reach_slope2
    )
    # plot the reconstructed WSE, if present
    try:
        w_opt = data.nodes['w_opt'][node_i]
        w_opt_r_u = data.nodes['w_opt_r_u'][node_i]
        opt_mask = w_opt > -999
        (_, caps, bars) = axis.errorbar(
            node_p_dist[opt_mask], w_opt[opt_mask], w_opt_r_u[opt_mask],
            label='opt WSE', linestyle=':', alpha=0.5, color='orange'
        )
    except AttributeError:
        # w_opt not in product
        print('No reconstructed WSE available, skipping...')

    # set grid, title, and labels
    axis.grid()
    axis.set_xlabel('dist from outlet (m)', fontsize=9)
    axis.set_ylabel('WSE (m)', fontsize=9)
    # Increase the fontsize of the tick labels for the primary axis
    axis.tick_params(axis='both', which='major',
                     labelsize=9)  # You can adjust the fontsize as needed
    if title is not None:
        axis.set_title(title[0:20])
    figure.suptitle(river_name + ' ' + str(reach_id)
                    + '_' + str(cycle))

    # create second top axis showing node ID
    axis2 = axis.twiny()
    node_id = node_id # get_simple_node_id(node_id, reach_id)
    axis2.plot(node_id, avg_wse * np.ones(len(node_id)))
    axis2.cla()
    axis2.xaxis.get_offset_text().set_visible(False)
    axis2.set_xlim(max(node_id), min(node_id))
    # Increase the fontsize of the tick labels for the secondary axis
    axis2.tick_params(axis='both', which='major',
                      labelsize=9)  # You can adjust the fontsize as needed

    if prd_heights:
        axis.plot(node_p_dist, data.nodes['p_wse'][node_i],
                  'D', markersize=2, label='PRD wse')

    # Add reach summary metrics in text to bottom left of plot
    # if annotate_metrics:
    plot_summary_metrics(data, axis, reach_i)
    if truth is not None:
        # match and plot data. If there are no matches, capture that in
        # has_truth flag
        has_truth, pt_wse_e, pt_stdev, pt_slp_e, axis2 = get_truth_matches(
            truth, reach_id, pass_no, tile, cycle, axis, axis2, fit_x, ss_min,
            ss_max, annotate_metrics, multi_reach)
        if not has_truth:
            plot_swot_only(axis, reach_id, cycle, pass_no, river_name)
            return has_truth
    else:
        # no truth was input to WSE plotter
        has_truth = False
        pt_wse_e, pt_stdev, pt_slp_e = None, None, None

    if errors is not None:
        slp_e = 'slope_e=' + str(round(errors['slp e (cm/km)'], 2)) + 'cm/km\n'
        wse_e = 'wse_e=' + str(round(errors['wse e (cm)'], 2)) + ' cm\n'
    else:
        slp_e, wse_e = None, None

    if errors is not None:
        slp_e = 'slope_e=' + str(round(errors['slp e (cm/km)'], 2)) + 'cm/km\n'
        wse_e = 'wse_e=' + str(round(errors['wse e (cm)'], 2)) + ' cm\n'
    if annotate_metrics:
        annotate_errors(axis, slp_e, wse_e, pt_wse_e, pt_stdev, pt_slp_e)

    leg = axis.legend(
        bbox_to_anchor=(1.01, 1.0), loc='upper left', fontsize=5, ncol=1)
    leg2 = axis2.legend(fontsize=5, ncol=1)
    leg.set_draggable(1)
    # flip axis so highest part of river is on the left
    axis.set_xlim(axis.get_xlim()[::-1])
    # call `draw` to re-render the graph
    plt.draw()
    plt.tight_layout()

    return has_truth


def get_truth_matches(truth, reach_id, pass_no, tile, cycle, axis, axis2,
                      fit_x, ss_min, ss_max, annotate_metrics, multi_reach=False):
    # Matches PT, drift, or truth rivertile data to the SWOT data and plots it.

    # initialize variables for storing matches (if they exist)
    has_truth = False
    pt_wse_e = None
    pt_stdev = None
    pt_slp_e = None
    pt_node_match = None
    drift_node_match = None
    coarse_node_match = None

    # check for truth type. Could be truth tile or calval dataframe
    if isinstance(truth, dict):
        print('You input a dict of field dataframes! Plotting...')
        if truth['pt_node'] is not None:
            # match and plot PT data for this reach-cycle
            if multi_reach:
                truth_i = np.logical_or.reduce((
                    truth['pt_node']['reach_id'] == str(reach_id),
                    truth['pt_node']['reach_id'] == str(reach_id - 10),
                    truth['pt_node']['reach_id'] == str(reach_id + 10)
                ))
                truth_pt = truth['pt_node'][truth_i]
            else:
                truth_pt = truth['pt_node'][
                    truth['pt_node']['reach_id'] == str(reach_id)]
            if not truth_pt.empty:
                pt_node_match = get_pt_match(truth_pt, pass_no, tile, cycle,
                                             reach_id)
                pt_error = truth['pt_error'][
                    truth['pt_error']['reach_id'] == reach_id]
                pt_err_match = get_pt_match(pt_error, pass_no, tile,
                    int(cycle), reach_id)
                if pt_node_match is not None:
                    print('Truth pt available for reach', reach_id,
                          'cycle', cycle)
                    pt_wse_e, pt_stdev, pt_slp_e = plot_pt_data(
                        pt_node_match, pt_err_match, axis, axis2, fit_x,
                        ss_min, ss_max, reach_id, annotate_metrics
                    )
        if truth['drift_node_match'] is not None:
            # match and plot toolbox drift for this reach-cycle
            if multi_reach:
                drift_i = np.logical_or.reduce((
                    truth['drift_node_match']['reach_id'] == str(reach_id),
                    truth['drift_node_match']['reach_id'] == str(reach_id - 10),
                    truth['drift_node_match']['reach_id'] == str(reach_id + 10)
                ))
                truth_drift = truth['drift_node_match'][drift_i]
            else:
                truth_drift = truth['drift_node_match'][
                    truth['drift_node_match']['reach_id'] == str(reach_id)]
            if not truth_drift.empty:
                drift_node_match = get_drift_match(truth_drift, pass_no, tile,
                                                   cycle, reach_id)
                if drift_node_match.empty:
                    drift_node_match = None
                else:
                    print('Truth toolbox drift available for reach', reach_id,
                          'cycle', cycle)
                    # TODO fix plotting here, i think y var is wrong
                    axis2.plot(truth_drift['p_dist_out'],
                               truth_drift['ellipsoid_height_m'], 'm+',
                               label='drift nodes')
        if truth['coarse_node'] is not None:
            # match and plot coarse drift for this reach-cycle
            if 'reach_id' not in truth['coarse_node'].columns:
                # get reach ID from node ID
                truth['coarse_node']['reach_id'] = truth[
                    'coarse_node']['node_id'].apply(
                    transform_reach_id_from_node_id)
            if multi_reach:
                coarse_i = np.logical_or.reduce((
                    truth['coarse_node']['reach_id'] == reach_id,
                    truth['coarse_node']['reach_id'] == reach_id - 10,
                    truth['coarse_node']['reach_id'] == reach_id + 10
                ))
                truth_coarse = truth['coarse_node'][coarse_i]
            else:
                truth_coarse = truth['coarse_node'][
                    truth['coarse_node']['reach_id'] == reach_id]
            if not truth_coarse.empty:
                coarse_node_match = get_drift_match(truth_coarse, pass_no,
                                                    tile, cycle, reach_id)
                if coarse_node_match.empty:
                    coarse_node_match = None
                else:
                    print('Truth coarse drift available for reach', reach_id,
                          'cycle', cycle)
                    try:
                        axis.plot(coarse_node_match['p_dist_out'],
                                   coarse_node_match['mean_node_drift_wse_m'], 'r+',
                                   label='coarse drift nodes')
                    except KeyError:
                        # try shapefile column name
                        axis.plot(coarse_node_match['p_dist_out_drift'],
                                   coarse_node_match['mean_node_drift_wse_m'], 'r+',
                                   label='coarse drift nodes')

        # Check if any match is not None to set has_truth
        if pt_node_match is not None or drift_node_match is not None or \
                coarse_node_match is not None:
            has_truth = True
    else:  # Assuming it is a rivertile netcdf object
        print('You input a truth rivertile! Plotting...')
        truth_df = pd.DataFrame.from_dict(truth['nodes'].variables)
        node_i_truth = truth.nodes['reach_id'] == reach_id
        truth_df = truth_df[node_i_truth]
        truth_df.set_index('node_id')
        node_p_dist_truth = truth.nodes['p_dist_out'][node_i_truth]
        truth_wse = truth.nodes['wse'][node_i_truth]
        reach_i_truth = truth.reaches['reach_id'] == reach_id
        truth_reach_wse = truth.reaches['wse'][reach_i_truth]
        truth_slope = truth.reaches['slope'][reach_i_truth]

    return has_truth, pt_wse_e, pt_stdev, pt_slp_e, axis2

def transform_reach_id_from_node_id(number):
    # changes node Id's to reach IDs. Assumes all reaches are river reaches
    # (i.e. type "1"). This isn't always true, but should be the case for all
    # calval rivers.
    num_str = str(number)
    # Remove the last three digits
    shortened_str = num_str[:-3]
    # Replace the last remaining digit with "1"
    if len(shortened_str) > 0:
        transformed_str = shortened_str[:-1] + "1"
    else:
        transformed_str = "1"  # In case slicing results in an empty string
    # Convert back to integer
    return int(transformed_str)

def annotate_errors(axis, slp_e, wse_e, pt_wse_e, pt_stdev, pt_slp_e):
    # puts the input errors on the input axis
    def add_text(y_offset, error_value, label, color_key):
        if error_value is not None:
            if 'wse' in label or 'stdev' in label:
                er_str = f'{label} = {round(error_value, 2)} cm\n'
            else:
                er_str = f'{label} = {round(error_value, 2)} cm/km\n'
            axis.text(RIGHT - 0.1, TOP - y_offset, er_str,
                horizontalalignment='left', verticalalignment='bottom',
                fontsize=5, color=get_passfail_color(error_value, color_key),
                transform=axis.transAxes)

    errors = [(0.08, slp_e, 'slp_e', 'slp e (cm/km)'),
        (0.16, wse_e, 'wse_e', 'wse e (cm)'),
        (0.22, pt_wse_e, 'pt_wse_e', 'wse e (cm)'),
        (0.28, pt_stdev, 'pt_stdev', 'wse e (cm)'),
        (0.34, pt_slp_e, 'pt_slp_e', 'slp e (cm/km)')]

    for offset, value, label, color_key in errors:
        add_text(offset, value, label, color_key)

def plot_summary_metrics(data, axis, reach_i):
    reach_width = get_first_element(data.reaches['width'][reach_i])
    reach_xtrk = get_first_element(data.reaches['xtrk_dist'][reach_i])
    reach_xtrk = str(round(np.mean(reach_xtrk) / 1000, 1))
    reach_dark_frac = get_first_element(data.reaches['dark_frac'][reach_i])
    reach_obs_frac = get_first_element(data.reaches['obs_frac_n'][reach_i])
    reach_xovr_cal_q = get_first_element(data.reaches['xovr_cal_q'][reach_i])
    summary_string = 'w = ' + str(round(reach_width, 2)) + ' m\n' \
                     + 'x-trk =' + reach_xtrk + ' km\n' \
                     + 'dark_frac = ' + str(round(reach_dark_frac, 2))  \
                     + '\n' + 'obs_frac = ' + str(round(reach_obs_frac, 2))\
                     + '\n' + 'xovr_q = ' + str(round(reach_xovr_cal_q, 2))
    axis.text(LEFT, BOTTOM, summary_string,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=8,
              transform=axis.transAxes,
              bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))


def plot_wse_and_qual(node_p_dist, wse, wse_r_u, node_q, node_q_b, axis,
                      plot_bit_qual, reach_wse, reach_slope, reach_slope2,
                      plot_slp2=False):
    # HACK to account for very bad wse_r_u values
    wse_r_u[wse_r_u < 0] = 0
    # Remove fill value nodes
    mask = wse >= -999
    node_p_dist, wse, wse_r_u, node_q, node_q_b = (arr[mask] for arr in
    [node_p_dist, wse, wse_r_u, node_q, node_q_b])

    # Plot all nodes
    axis.errorbar(node_p_dist, wse, yerr=wse_r_u, fmt='o', markersize=2,
                  label='node wse', zorder=0)

    # Quality masks and colors
    quality_masks = {
        'bad': (node_q == 3, 'red', 'bad qual'),
        'degraded': (node_q == 2, 'orange', 'degraded qual'),
        'suspect': (node_q == 1, 'yellow', 'suspect qual'),
        'good': (node_q == 0, '#ADD8E6', 'good qual')
    }

    # Plot nodes by quality
    for mask, color, label in quality_masks.values():
        axis.errorbar(node_p_dist[mask], wse[mask], wse_r_u[mask], fmt='o',
                      color=color, markersize=4, markerfacecolor=color,
                      markeredgecolor='black', markeredgewidth=1, label=label,
                      zorder=1)

    if plot_bit_qual:
        # Annotate with node_q_b
        for x, y, z in zip(node_p_dist, wse, node_q_b):
            axis.text(x, y + 0.02, z, fontsize=4, color='black')

    # Mark outlier nodes
    outlier_mask = (node_q_b &
                    RiverTileNodes.VARIABLES['node_q_b']['flag_masks'][
                        RiverTileNodes.VARIABLES['node_q_b'][
                            'flag_meanings'].split(' ').index(
                            'wse_outlier')]) > 0
    axis.plot(node_p_dist[outlier_mask], wse[outlier_mask], 'mo', markersize=4,
              label='outlier')

    # Adjust y-axis limits
    if len(wse) > 0:
        axis.set_ylim(min(wse) - 1.5, max(wse) + 1.5)

    # Plot the SWOT observed reach slope
    reach_center_dist = np.mean(node_p_dist)
    ss = node_p_dist - reach_center_dist
    fit_x = np.array([min(ss), 0, max(ss)]) + reach_center_dist

    # Observed fit
    obs_fit_y = [reach_wse + min(ss) * reach_slope, reach_wse,
                 reach_wse + max(ss) * reach_slope]
    # axis.plot(fit_x, obs_fit_y, '--', markersize=10, color='b',
    #           label='obs fit')
    # axis.plot(reach_center_dist, reach_wse, 'b*', markersize=9, color='g',
    #           label='obs wse', zorder=1)
    # axis.axvline(x=reach_center_dist, ls='--', lw=0.5)

    # WSE RU shading
    axis.fill_between(node_p_dist, wse + 3 * wse_r_u, wse - 3 * wse_r_u,
                      facecolor='gray', alpha=0.3, interpolate=True)

    # Enhanced reach slope
    if plot_slp2:
        obs_fit_y2 = [reach_wse + min(ss) * reach_slope2, reach_wse,
                      reach_wse + max(ss) * reach_slope2]
        axis.plot(fit_x, obs_fit_y2, '--', markersize=10, color='g',
                  label='slp2 fit')

    return fit_x, min(ss), max(ss)


def plot_pt_data(pt_node_match, pt_err_match, axis, axis2, fit_x, ss_min,
                 ss_max, reach_id, annotate_metrics):
    # plot the pt data on the WSE axis
    pt_wse = pt_node_match.mean_node_pt_wse_m
    pt_node_id = pt_node_match['node_id'] # get_simple_node_id(pt_node_match['node_id'], reach_id)
    pt_q = pt_node_match['pt_qual']

    # Define a custom colormap for PT node qual
    color_map = {0: 'green', 1: 'orange', 2: 'red'}

    colors = pt_q.map(color_map)
    # Plot the points with colors based on pt_qual
    scatter = axis2.scatter(pt_node_id, pt_wse, c=colors, marker='x',
                            label='PT node WSE', zorder=10, s=30)
    if annotate_metrics:
        for node_id, wse, q in zip(pt_node_id, pt_wse, pt_node_match['flag']):
            axis2.text(node_id, wse, str(q), fontsize=6, color='blue')

    if not pt_err_match.empty:
        # Usage of the function to set variables
        pt_wse_e = assign_value_or_none(pt_err_match, 'wse_error_cm')
        try:
            pt_stdev = assign_value_or_none(pt_err_match, 'rel_wse_error_cm')
        except KeyError:
            pt_stdev = assign_value_or_none(pt_err_match, 'relative_wse_error_cm')
        try:
            pt_reach_wse = assign_value_or_none(pt_err_match,
                                                'mean_reach_pt_wse_m')
            pt_slope = assign_value_or_none(pt_err_match, 'slope_m_m')
            pt_slp_e = assign_value_or_none(pt_err_match, 'slp_error_cmkm')
        except KeyError:
            # not a reach dataframe
            pt_reach_wse = None
            pt_slp_e = None
            pt_slope = None
        if pt_reach_wse is not None:
            axis2.plot(np.mean(pt_node_id), pt_reach_wse, 'r*', markersize=8,
                       label='PT reach WSE', zorder=0)
            if pt_slope is not None:
                pt_fit_y = [pt_reach_wse + ss_min * pt_slope, pt_reach_wse,
                            pt_reach_wse + ss_max * pt_slope]
                axis.plot(fit_x, pt_fit_y, '--', markersize=10, color='r',
                          label='PT fit')
    else:
        pt_wse_e = None
        pt_stdev = None
        pt_slp_e = None
    return pt_wse_e, pt_stdev, pt_slp_e


def plot_swot_only(axis, reach_id, cycle, pass_no, river_name):
    # make the WSE/slope plot for SWOT data only. Should be on the axis already
    leg = axis.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize=5,
                      ncol=1)
    leg.set_draggable(1)
    print('No PT nor drift node dataframe for Reach:', reach_id,
          '\n                         Cycle:', cycle,
          '\n                         Pass:', pass_no,
          '\n                         River:', river_name)
    axis.set_xlim(axis.get_xlim()[::-1])
    # call `draw` to re-render the graph
    plt.draw()


def get_pt_match(truth_pt, pass_no, tile, cycle, reach_id):
    # gets matching rows of a PT truth node dataframe
    if not truth_pt.empty:
        print(
            'Matching PT dataframes to pass/tile', pass_no, tile,
            'cycle', cycle,
            'reach_id', reach_id
        )
        truth_pt['pass'] = truth_pt['pass'].astype(str).str.zfill(3)
        pt_match = truth_pt[truth_pt['pass'] == pass_no]
        pt_match = pt_match[pt_match['cycle'] == cycle]
        return pt_match
    else:
        return None


def get_drift_match(drift, pass_no, tile, cycle, reach_id):
    # gets matching rows of a PT truth node dataframe
    if not drift.empty:
        print(
            'Matching drift node dataframes to pass/tile', pass_no, tile,
            'cycle', cycle,
            'reach_id', reach_id
        )
        drift_match = drift[
            drift['pass'].astype(str).str.zfill(3) == pass_no]
        drift_match = drift_match[
            drift_match['cycle'].astype(str).str.zfill(3) == cycle]
        return drift_match
    else:
        return None


def plot_area(data, truth, errors, reach_id, axis, title=None, style='.',
              multi_reach=False):
    # plot the truth and observed area, for detected and total
    if multi_reach:
        # grab adjacent reaches too if they exist
        node_i = np.logical_or.reduce((
            data.nodes['reach_id'] == reach_id,
            data.nodes['reach_id'] == reach_id - 10,
            data.nodes['reach_id'] == reach_id + 10,
            np.logical_not(data.nodes['area_total'].mask)
        ))
    else:
        node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
                            np.logical_not(data.nodes['area_total'].mask))
    node_id = data.nodes['node_id'][node_i]
    node_id = get_simple_node_id(node_id, reach_id)

    area_detct = data.nodes['area_detct'][node_i]
    area_total = data.nodes['area_total'][node_i]  # includes dark water pixels

    axis.plot(node_id, area_detct, style, markersize=4, alpha=.5)
    axis.plot(node_id, area_total, style, markersize=4, alpha=.5)

    if truth is not None:
        if isinstance(truth, dict):
            truth = None  # TODO: bring area dataframes in from fineval dataset
        else:
            node_i_truth = np.logical_and(
                truth.nodes['reach_id'] == reach_id,
                np.logical_not(truth.nodes['wse'].mask)
            )
            node_id_truth = truth.nodes['node_id'][node_i_truth]
            node_id_truth = get_simple_node_id(node_id_truth, reach_id)
            area_truth = truth.nodes['area_total'][node_i_truth]
            axis.plot(node_id_truth, area_truth, 'kx', markersize=2)

    # add text with error summary
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
    axis.set_xlim(axis.get_xlim()[::-1])  # flip axis to align with WSE plot
    if title is not None:
        axis.set_title(title)


def plot_pix_assgn(data, reach_id, axis, h_flg=False, area_flg=False,
                   multi_reach=False):
    # Filter data for the specified reach_id
    if multi_reach:
        # grab adjacent reaches too if they exist
        pix_i = np.logical_or.reduce((
            data['reach_id'] == reach_id,
            data['reach_id'] == reach_id - 10,
            data['reach_id'] == reach_id + 10
        ))
    else:
        pix_i = (data['reach_id'] == reach_id)
    node_id = data['node_id'][pix_i]
    lat = data['latitude_vectorproc'][pix_i]
    lon = data['longitude_vectorproc'][pix_i]

    # Check if node_id has valid data and isn't empty
    if node_id.size > 0 and np.all(np.isfinite(node_id)):
        vmin = node_id.min()
        vmax = node_id.max()

        # Check if all values are the same
        if vmin == vmax:
            # Adjust vmax slightly to avoid zero division
            vmax += 1  # or some small epsilon value specific to your data scale

        # Create a normalization object
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        # Fallback normalization (or handle the case specifically as needed)
        norm = mcolors.Normalize(vmin=0,
                                 vmax=1)  # Default or dummy normalization

    # Create a scatter plot manually handling colors
    scatter = axis.scatter(lon, lat, c=node_id, cmap='tab20b', norm=norm,
                           linewidth=0, alpha=0.7, s=5)

    # Set plot properties
    axis.grid(True)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('lon')
    plt.ylabel('lat')

    # get min and max for colorbar
    if not np.ma.is_masked(node_id) and np.ma.count(node_id) > 0:
        vmin = node_id.min()
        vmax = node_id.max()
    else:
        # Set a default minimum value if the array is empty or masked
        vmin = 0
        vmax = 1
    # Manually create ScalarMappable object
    norm = mcolors.Normalize(vmin=vmin,
                             vmax=vmax)  # Define normalization
    sm = plt.cm.ScalarMappable(cmap='tab20b',
                               norm=norm)
    sm.set_array([])  # Dummy array to satisfy ScalarMappable

    # Add colorbar
    colorbar = plt.colorbar(sm, ax=axis)
    colorbar.set_label('Node ID')

    if h_flg:
        # plot bad h_flg over existing plot
        h_flg = data['h_flg'][pix_i]
        bad_wse_lat = lat[h_flg == 0]
        bad_wse_lon = lon[h_flg == 0]
        axis.scatter(bad_wse_lon, bad_wse_lat, color='k', s=0.1)
        axis.set_title('Pixel Locations, h_flg/node_ID')
        axis.plot([], [], 'ko', label='h_flg=False')
        axis.legend()
    elif area_flg:
        area_flg = data['area_flg'][pix_i]
        bad_area_lat = lat[area_flg == 0]
        bad_area_lon = lon[area_flg == 0]
        axis.scatter(bad_area_lon, bad_area_lat, color='k', s=0.1)
        axis.set_title('Pixel Locations, area_flg')
        axis.plot([], [], 'ko', label='area_flg=False')
        axis.legend()


def plot_locations(data, truth, reach_id, axis, plot_prior=True, title=None):
    # creates plot with the observation centroids and the prior node locations
    reach_id = int(reach_id)
    node_i = np.logical_and(data.nodes['reach_id'] == reach_id,
                            np.logical_not(data.nodes['wse'].mask))
    node_id = data.nodes['node_id'][node_i]
    if truth is not None:
        if isinstance(truth, dict):
            truth = None  # TODO: bring more truth dataframes from fineval data
        else:
            node_i_truth = np.logical_and(
                truth.nodes['reach_id'] == reach_id,
                np.logical_not(truth.nodes['wse'].mask)
            )
    lat = data.nodes['lat'][node_i]
    lon = data.nodes['lon'][node_i]

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
    elif passfail[parameter][0] < abs(error_value) < passfail[parameter][1]:
        return 'orange'
    else:
        return 'red'


def decode_rivertile_filename(fname):
    """parse the filename"""
    _, tail = os.path.split(fname)
    parts = tail.split('_')
    if 'Reach' in fname:
        cycle = parts[5]
        pas = parts[6]
        tile = parts[7]
    else:
        cycle = parts[4]
        pas = parts[5]
        tile = parts[6]
    return cycle, pas, tile


def check_for_existing_river_plots(out_dir, title, river_code):
    # returns whether the river plots exist, and whether they have truth or not
    filename1 = out_dir + '/truth/' + title + '_' + river_code + '.png'
    filename2 = out_dir + '/no_truth/' + title + '_' + river_code + '.png'
    if os.path.exists(filename1):
        print('River plot', filename1, 'already exists, continuing...')
        return True, True  # the plots already exist in /truth/ folder
    elif os.path.exists(filename2):
        print('River plot', filename2, 'already exists, continuing...')
        return True, False  # the plots already exist in /no_truth/ folder
    else:
        return False, None  # the plots do not already exist


def check_for_existing_pixc_plots(out_dir, title):
    filename1 = out_dir + '/truth/' + title + '_PIXC.png'
    filename2 = out_dir + '/no_truth/' + title + '_PIXC.png'
    if os.path.exists(filename1) or os.path.exists(filename2):
        print(title + '_PIXC' + ' already exists, continuing...')
        return True  # the plots already exist
    else:
        return False  # the plots do not already exist


def make_pixc_plots(
        pixcvec_data, river_data, pixc, truth_pixcvec, truth_pixc, reach_id,
        nodes=None, pixc_truth=None, out_dir=None, title=None, has_truth=None):

    plt.tight_layout()
    mngr = plt.get_current_fig_manager()
    # mngr.window.setGeometry(0, 0, 1500, 500)
    if pixc and pixcvec_data:
        pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc)
        pixc_truth_data = None
        if pixc_truth is not None:
            pixc_truth_data = SWOTWater.products.product.MutableProduct.from_ncfile(
                pixc_truth)
        slant_plane_fig = plot_pixcs(
            pixcvec_data, pixc_data, reach_id, nodes, reach_data=river_data,
            pixc_truth=pixc_truth_data
        )
        if out_dir is not None:
            # save current figure to file
            plt.title(title, backgroundcolor='white')
            if has_truth:
                filename = out_dir + '/truth/' + title + '_PIXC'
            else:
                filename = out_dir + '/no_truth/' + title + '_PIXC'
            this_river_name = river_data['reaches']['river_name'][
                river_data['reaches']['reach_id'] == reach_id][0]
            fig_title = this_river_name + ' ' + title
            slant_plane_fig.suptitle(fig_title)
            slant_plane_fig.savefig(filename)
            plt.close()
        else:
            plt.title(title, backgroundcolor='white')
            plt.show()

    else:
        print('Missing pixc or pixcvec file, skipping pixel assignment plot')

    if pixc and truth_pixc:  # only plot these if pixc was also given
        truth_pixcvec_data = SWOTWater.products.product.MutableProduct.from_ncfile(
            truth_pixcvec)
        truth_pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(
            truth_pixc)
        plot_pixcs(truth_pixcvec_data, truth_pixc_data, reach_id, nodes,
                   title_tag='(truth)', reach_data=truth)


def get_river_code(river_data, reach_id):
    river_name = get_first_element(river_data.reaches['river_name'][
        river_data.reaches['reach_id'] == reach_id])
    # Dictionary mapping river names to their corresponding codes
    river_codes = {'Willamette River': 'WM',
                   'Waimakariri River': 'WK',
                   'Porcupine River': 'PY',
                   'Yukon River': 'YR',
                   'Grand Canyon': 'GC',
                   'Colorado River': 'GC',
                   'Tanana River': 'TN',
                   'South Santiam River': 'WM',
                   'Santiam River': 'WM',
                   'Santiam River; South Santiam River': 'WM',
                   'Connecticut River': 'CR',
                   'North Saskatchewan River': 'NS',
                   'Peace River': 'PD',
                   'Slave River': 'PD',
                   'Sagavanirktok River': 'SG',
                   'Connecticut River; Westfield River': 'CR',
                   'Mania': 'MNR',
                   'Lawa': 'LW',
                   "L'Aussonnelle": 'LAS',
                   '181601': '181',
                   'La Garonne': 'GR',
                   'Tsiribihina': 'TS',
                   'Merrimack River': 'MRK',
                   'Maroni': 'MRN'}
    # Get the river code from the dictionary, use a default if river not found
    river_code = river_codes.get(river_name, 'N-A')
    return river_code


def make_river_plots(rivertile_file, river_data, truth_data, pixcvec, reach_id,
                     errors=None, out_dir=None, title=None, multi_reach=False):
    # contains node group and reach group for each input netcdf
    cycle, pass_no, tile = decode_rivertile_filename(rivertile_file)
    if truth_data is not None:
        if isinstance(truth_data, str):
            truth = SWOTWater.products.product.MutableProduct.from_ncfile(
                truth_file)
        elif isinstance(truth_data, dict):
            # truth input is a dictionary of field dataframes
            truth = truth_data
    else:
        truth = None
    if truth is None:
        has_truth = False
    title_str = str(cycle) + '_' + str(reach_id)
    # remap river name
    def remap_river_names(river_data):
        # Define a mapping of old names to new names
        river_name_mapping = {
            '822820': 'PAD',
            '822810': 'PAD',
            '572058': 'Waimakariri River',
            'Ribdon River': 'Sagavanirktok River',
            'L\'Aussonelle': 'La Garonne',
            '232141': 'La Garonne',
            'Lawa': 'Maroni'}

        # Create a mask for each remapping and apply the new name
        for old_name, new_name in river_name_mapping.items():
            mask = river_data.reaches['river_name'] == old_name
            river_data.reaches['river_name'][mask] = new_name

        return river_data
    river_data = remap_river_names(river_data)

    figure, axes = plt.subplots(2, 2, figsize=FIGSIZE, dpi=DPI)
    has_truth = plot_wse(
        river_data, truth, errors, reach_id, axes[0][0], figure,
        title=title_str + ' - wse', cycle=cycle, tile=tile, pass_no=pass_no,
        multi_reach=multi_reach
    )
    plot_area(river_data, truth, errors, reach_id, axes[1][0],
              title=title_str + ' - area', multi_reach=multi_reach)
    # uncomment the below block if you'd prefer to plot the centroids rather
    # than the area flag
    # plot_locations(river_data, truth, reach_id, axes[0][1],
    #                title=title_str + ' - locations')
    if pixcvec is not None:
        plot_pix_assgn(pixcvec, reach_id, axes[0][1], area_flg=True,
                       multi_reach=multi_reach)
        plot_pix_assgn(pixcvec, reach_id, axes[1][1], h_flg=True,
                       multi_reach=multi_reach)

    plt.tight_layout()
    if out_dir is not None:
        # save current figure to file
        river_code = get_river_code(river_data, reach_id)
        plt.title(title, backgroundcolor='white')
        if has_truth:
            filename = out_dir + '/truth/' + title + '_' + river_code
        else:
            filename = out_dir + '/no_truth/' + title + '_' + river_code
        plt.savefig(filename)
        plt.close()
    else:
        plt.title(title, backgroundcolor='white')
        plt.show()
    return figure, axes, has_truth


def make_plots(rivertile_file, rivertile_df, truth_data, pixcvec, pixc,
               truth_pixcvec, truth_pixc, reach_id, errors=None,
               nodes=None, pixc_truth=None, out_dir=None, title=None,
               overwrite=False):

    # handle overwriting if user says not to
    make_rivers = True
    make_pixc = True
    if ~overwrite:
        # Don't write files if they exist already!
        river_code = get_river_code(rivertile_df, reach_id)
        river_exists, has_truth = check_for_existing_river_plots(
            out_dir, title, river_code)
        pixc_exists = check_for_existing_pixc_plots(out_dir, title)
        if river_exists:
            make_rivers = False
        if pixc_exists:
            make_pixc = False

    reach_id = int(reach_id)  # ensure reach ID inputs are integer type
    if make_rivers or make_pixc:
        # import the RiverTile and PIXCVecRiver data
        if pixcvec is not None:
            pixcvec_data = SWOTWater.products.product.MutableProduct.from_ncfile(
                pixcvec)
        else:
            pixcvec_data = None

    if make_rivers:
        fig, ax, has_truth = make_river_plots(
            rivertile_file, rivertile_df, truth_data, pixcvec_data, reach_id,
            errors=errors, out_dir=out_dir, title=title)
    if make_pixc:
        make_pixc_plots(
            pixcvec_data, rivertile_df, pixc, truth_pixcvec, truth_pixc, reach_id,
            nodes=nodes, pixc_truth=pixc_truth, out_dir=out_dir, title=title,
            has_truth=has_truth)
    plt.close()


def get_reach_error(errors, reach_id):
    # this gets the slope, wse, and area errors for the reach of interest
    reach_error = {}
    index = errors[0]['reach'].index(str(reach_id))
    for key in errors[0].keys():
        reach_error[key] = errors[0][key][index]

    return reach_error


def plot_pixcs(pixc_vec, pixc, reach_id, nodes=None,
               title_tag='(slant-plane)', reach_data=None, pixc_truth=None,
               apply_corr=True, plot_qual=True):
    reach_id = int(reach_id)
    # get only the reach_id for pixels in pixc_vec
    pix_i = (pixc_vec['reach_id'] == reach_id)
    slant_plane_fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    if len(pix_i) > 0:
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

        ax1 = plt.subplot(2, 3, 1)
        pt1 = ax1.imshow(Node_id, interpolation='none', aspect='auto',
                         cmap=plt.cm.get_cmap('tab20b'))
        plt.colorbar(pt1, ax=ax1)
        ax1.set_title('node_id ' + title_tag)

        # TODO: make a better cmap for classification, also make font bigger
        ax2 = plt.subplot(2, 3, 2, sharex=ax1, sharey=ax1)
        class_cmap = colors.ListedColormap(
            ['pink', 'darkgreen', 'lightgreen', 'aquamarine', 'blue', 'black',
             'yellow', 'red']
        )
        pt2 = ax2.imshow(Cls1, interpolation='none', aspect='auto',
                         cmap=class_cmap, clim=(0, 7))
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
                if len(area_tot[node_id == node_id1]) > 0:
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
        else:
            print('No reach ID matching', reach_id, 'in this pixel cloud')
    return slant_plane_fig


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
                        help='list of nodes for which to plot height histograms',
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
