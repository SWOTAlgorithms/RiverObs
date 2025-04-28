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

from plot_tile_reaches import load_pkl_dataframes

# TODO: maybe should do a try-catch so it is not required
#       to install these next ones?
# import things for Google Maps in plots
from matplotlib_scalebar.scalebar import ScaleBar
import cartopy.crs as ccrs
from cartopy.io.img_tiles import GoogleTiles
# for multitemporal stuff in plots
import rivscale.products.along_stretch


FIGSIZE = (16, 9)
DPI = 200
LEFT, WIDTH = .04, .75
RIGHT = LEFT + WIDTH
BOTTOM, HEIGHT = .02, .85
TOP = BOTTOM + HEIGHT

#matplotlib.rcParams.update({'font.size': 6})
matplotlib.rcParams.update({'font.size': 9})

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

def populate_rivertile(proc_df, node_df, cycle_id, pass_id):
    # populate node data from node_df to proc_df
    common_keys = set(proc_df.nodes.VARIABLES.keys()).intersection(
            set(node_df.keys()))
    for key in common_keys:
        proc_df.nodes[key] = np.ma.masked_array(node_df[key])
    proc_df.nodes.cycle_number = cycle_id
    proc_df.nodes.pass_number = pass_id
    proc_df.nodes.wse[np.abs(proc_df.nodes.wse)>1e6] = np.nan
    proc_df.nodes.width[np.abs(proc_df.nodes.width)>1e6] = np.nan

    # populate some reach variables from nodes
    ignore_keys = [
                    'time_str',
                    'centerline_lat',
                    'centerline_lon',
                    'rch_id_up',
                    'rch_id_dn']
    keys = set(proc_df.reaches.VARIABLES.keys()) - set(ignore_keys)
    #print(keys)
    for key in keys:
        #print(key)
        if key == 'river_name':
            proc_df.reaches[key] = np.ma.masked_array(
                np.unique(node_df[key]))
        elif key == 'reach_id':
            #rid = np.array([int(i) for i in node_df[key]])
            proc_df.reaches[key] = np.ma.masked_array(
                np.unique(node_df[key]))
        elif key in proc_df.nodes.variables.keys():
            proc_df.reaches[key] = np.ma.masked_array(
                [np.mean(proc_df.nodes[key]),])# just average the others
        else:
            #print(key)
            proc_df.reaches[key] = np.ma.masked_array([np.nan,])
    return proc_df

def load_wse_data(pixc_data):
    p_height = toslant(pixc_data.pixel_cloud, 'height')
    p_geoid = toslant(pixc_data.pixel_cloud, 'geoid')
    p_solid = toslant(pixc_data.pixel_cloud, 'solid_earth_tide')
    p_load = toslant(pixc_data.pixel_cloud, 'load_tide_fes')
    p_pole = toslant(pixc_data.pixel_cloud, 'pole_tide')
    pixc_var = p_height - (p_geoid + p_solid + p_load + p_pole)
    return pixc_var

def get_pixc_var(pixc_data, pixcvec_data, reach_id, var='height'):
    # handle special cases
    if var=='height_u':
        p_var = toslant(pixc_data.pixel_cloud, 'phase_noise_std') * np.abs(
                toslant(pixc_data.pixel_cloud, 'dheight_dphase'))
    elif var=='wse':
        p_var = load_wse_data(pixc_data)
    else:
        p_var = toslant(pixc_data.pixel_cloud, var)

    pix_i = (pixcvec_data['reach_id'] == reach_id)
    azimuth_index_vec = pixcvec_data['azimuth_index'][pix_i]
    range_index_vec = pixcvec_data['range_index'][pix_i]
    #
    var1 = p_var[azimuth_index_vec, range_index_vec]
    az0 = np.min(azimuth_index_vec)
    r0 = np.min(range_index_vec)
    # get 2D cropped around reach
    var_pixc = np.zeros((
        np.max(azimuth_index_vec-az0)+1,
        np.max(range_index_vec-r0)+1)) + np.nan
    var_pixc[azimuth_index_vec-az0, range_index_vec-r0] = var1
    return var_pixc, p_var, pix_i, azimuth_index_vec, range_index_vec, az0, r0

def sandbox_run(rivertile_df, pixcvec_data, pixc_data, reach_id):
    """
    do some experimental things
    """
    # get PIXC variables for pixcvec pixels around this reach
    height, p_height, pix_i, azimuth_index_vec, range_index_vec, az0, r0 = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='wse')
    height_u, p_height_u, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='height_u')
    geo_qual, p_geo_qual, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='geolocation_qual')
    klass, p_klass, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='classification')
    klass_qual, p_klass_qual, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='classification_qual')
    prob, p_prob, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='prior_water_prob')
    change, p_change, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='prior_water_change')
    water_frac, p_water_frac, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='water_frac')
    pixc_area, p_pixc_area, _, _, _, _, _ = get_pixc_var(
            pixc_data, pixcvec_data, reach_id, var='pixel_area')
    area_frac = pixc_area.copy()
    edges = np.where(np.logical_or(klass==2, klass==3))
    area_frac[edges] = area_frac[edges] * water_frac[edges]
    # get some pixcvec pixel-wise variables
    node_id = pixcvec_data['node_id'][pix_i]
    outlier = np.zeros((
        np.max(azimuth_index_vec-az0)+1,
        np.max(range_index_vec-r0)+1))
    IQR = []
    ref = []
    ref_sm = []
    #breakpoint()
    for node in np.unique(node_id):
        #breakpoint()
        print(node)
        # compute IQR
        pix_n = (node_id == node)
        az_index_vec = azimuth_index_vec[pix_n]
        rng_index_vec = range_index_vec[pix_n]
        #p_heightn = toslant(pixc_data.pixel_cloud, 'height')
        #p_classn = toslant(pixc_data.pixel_cloud, 'classification')
        heightn = p_height[az_index_vec, rng_index_vec]
        klassn = p_klass[az_index_vec, rng_index_vec]
        probn =p_prob[az_index_vec, rng_index_vec]
        good = np.logical_and(klassn==4, probn>0.05)
        if np.sum(good)<10:
            good = probn>0.05
        p25 = np.nanpercentile(heightn[good], 25)
        p75 = np.nanpercentile(heightn[good], 75)
        p50 = np.nanpercentile(heightn[good], 50)
        IQR.append(p75 - p25)
        ref.append(p50)
        #
        wse_sm = rivertile_df.nodes.wse_sm[
                rivertile_df.nodes.node_id==node]
        ref_sm.append(wse_sm)
    flag = []
    sig = []
    ref_a = []
    ref_a2 = []
    ref_5 = []
    ref_25 = []
    ref_50 = []
    ref_75 = []
    for node,ref0 in zip(np.unique(node_id),ref):
        #breakpoint()
        print(node)
        # compute IQR
        pix_n = (node_id == node)
        az_index_vec = azimuth_index_vec[pix_n]
        rng_index_vec = range_index_vec[pix_n]
        #p_heightn = toslant(pixc_data.pixel_cloud, 'height')
        heightn = p_height[az_index_vec, rng_index_vec]
        heightn_u = p_height_u[az_index_vec, rng_index_vec]
        #p25 = np.nanpercentile(var_pixcn, 25)
        #p75 = np.nanpercentile(var_pixcn, 75)
        #p50 = np.nanpercentile(hein, 50)
        #IQR = p75 - p25
        med_hu = np.median(heightn_u)
        sig.append(med_hu)
        upper = ref0 + 5 * med_hu#heightn_u#np.median(height1_u)#IQR)
        lower = ref0 - 5 * med_hu#heightn_u#np.median(height1_u)#IQR)
        msk = np.logical_or(heightn>upper, heightn<lower)
        outlier[az_index_vec[msk]-az0, rng_index_vec[msk]-r0] = 1
        frac = np.sum(msk) / len(heightn)
        flag.append(frac)
        arean =p_pixc_area[az_index_vec, rng_index_vec]
        wfn =p_water_frac[az_index_vec, rng_index_vec]
        klassn = p_klass[az_index_vec, rng_index_vec]
        probn =p_prob[az_index_vec, rng_index_vec]
        area_fracn = arean.copy()
        inds = np.where(np.logical_or(klassn==2, klassn==3))
        area_fracn[inds] = arean[inds] * wfn[inds]
        area = np.nansum(area_fracn)
        area2 = np.nansum(area_fracn[outlier[
            az_index_vec-az0, rng_index_vec-r0]==0])
        ref_a.append(area)
        ref_a2.append(area2)
        area_5 = np.sum(arean[probn>0.05])
        area_25 = np.sum(arean[probn>0.25])
        area_75 = np.sum(arean[probn>0.75])
        area_50 = np.sum(arean[probn>0.5])
        ref_5.append(area_5)
        ref_25.append(area_25)
        ref_50.append(area_50)
        ref_75.append(area_75)
    ref = np.array(ref)
    ref_sm = np.array(ref_sm)
    IQR = np.array(IQR)
    sig = np.array(sig)
    ref_a = np.array(ref_a)
    ref_a2 = np.array(ref_a2)
    ref_5 = np.array(ref_5)
    ref_25 = np.array(ref_25)
    ref_50 = np.array(ref_50)
    ref_75 = np.array(ref_75)
    # zero out widths with no data
    ref_5[~np.isfinite(ref_5)] = 0
    ref_25[~np.isfinite(ref_25)] = 0
    ref_50[~np.isfinite(ref_50)] = 0
    ref_75[~np.isfinite(ref_75)] = 0
    kwargs = {'interpolation':'none', 'aspect':'auto'}
    plt.figure()
    plt.imshow(outlier, cmap='gray_r', **kwargs)
    plt.imshow(height, cmap='jet',alpha=0.5, **kwargs)
    plt.colorbar()
    plt.title('height')
    #
    plt.figure()
    plt.imshow(height_u, clim=(0,1), cmap='jet', **kwargs)
    plt.colorbar()
    plt.title('height_u')
    #
    plt.figure()
    plt.imshow(klass_qual, cmap='tab20', **kwargs)
    plt.colorbar()
    plt.title('classification qual')
    #
    plt.figure()
    plt.imshow(prob, cmap='jet', **kwargs)
    plt.colorbar()
    plt.title('prior water prob')
    #
    plt.figure()
    plt.imshow(change, clim=(-1,1), cmap='jet', **kwargs)
    plt.colorbar()
    plt.title('prior water change')
    #
    plt.figure()
    plt.imshow(water_frac, cmap='jet', **kwargs)
    plt.colorbar()
    plt.title('water_frac')
    #
    plt.figure()
    plt.imshow(area_frac, cmap='jet', **kwargs)
    plt.colorbar()
    plt.title('area_frac')
    #
    plt.figure()
    plt.imshow(klass, cmap='jet', **kwargs)
    plt.colorbar()
    plt.title('classification')
    #
    #plt.figure()
    #plt.imshow(bright, cmap='jet', **kwargs)
    #plt.colorbar()
    #
    h_qual = IQR / sig
    x = np.arange(len(ref))
    plt.figure()
    plt.plot(x,ref, label='med')
    plt.plot(x,ref_sm, label='bayes')
    plt.plot(x[h_qual>3], ref[h_qual>3],'o', label='flag 2')
    plt.plot(x[h_qual>5], ref[h_qual>5],'x', label='flag 5')
    plt.title('height')
    plt.legend()

    plt.figure()
    plt.plot(x,ref_a, label='data')
    plt.plot(x,ref_a2, label='data no outlier pixels')
    plt.plot(x,ref_5, label='5%')
    plt.plot(x,ref_25, label='25%')
    plt.plot(x,ref_50, label='50%')
    plt.plot(x,ref_75, label='75%')
    #plt.plot(x[h_qual>3], ref[h_qual>3],'o')
    #plt.plot(x[h_qual>5], ref[h_qual>5],'x')
    plt.title('area')
    plt.legend()

    """
    #
    plt.figure()
    plt.plot(IQR)
    plt.plot(sig)
    #
    plt.figure()
    plt.plot(IQR / sig)
    plt.grid()
    #
    """
    plt.show()

    breakpoint()

def get_IQR_range(y, ptiles, ptile_list,
        scale_shade=3, scale_lim=5,
        p_low=None, p_high=None):
    p25 = ptiles[:,ptile_list==25].squeeze()
    p75 = ptiles[:,ptile_list==75].squeeze()
    IQR = p75 - p25
    y_low = y - scale_shade * IQR
    y_high = y + scale_shade * IQR
    if p_low is not None:
        #just return the p_low %ile
        y_low = ptiles[:,ptile_list==p_low].squeeze()
    if p_high is not None:
        #just return the p_low %ile
        y_high = ptiles[:,ptile_list==p_high].squeeze()
    y_min = np.min(y - scale_lim * IQR)
    y_max = np.max(y + scale_lim * IQR)
    return y_low, y_high, y_min, y_max

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
             multi_reach=False, mt_wse=None):
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
        if mt_wse is not None:
            mt_wse = mt_wse.crop_to_reach()
    node_id = data.nodes['node_id'][node_i]
    node_q = data.nodes['node_q'][node_i]
    node_q_b = data.nodes['node_q_b'][node_i]
    node_p_dist = data.nodes['p_dist_out'][node_i]
    p_dist_out = data.nodes['p_dist_out'][node_i]
    along_dist = np.cumsum(data.nodes['p_length'][node_i])
    along_dist = np.max(along_dist) - along_dist # make it downstream
    node_id = get_simple_node_id(node_id, reach_id)
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
        along_dist,#node_p_dist,
        wse, wse_r_u, node_q, node_q_b, axis, plot_bit_qual,
        reach_wse, reach_slope, reach_slope2, mt_wse=mt_wse
    )
    # plot the reconstructed WSE, if present
    try:
        w_opt = data.nodes['wse_sm'][node_i]
        w_opt_u = data.nodes['wse_sm_u'][node_i]
        opt_mask = w_opt > -999
        (_, caps, bars) = axis.errorbar(
            along_dist[opt_mask],#node_p_dist[opt_mask],
            w_opt[opt_mask], w_opt_u[opt_mask],
            label='opt WSE', linestyle=':', alpha=0.5, color='orange'
        )
    except AttributeError:
        # w_opt not in product
        print('No reconstructed WSE available, skipping...')

    # set grid, title, and labels
    axis.grid()
    #axis.set_xlabel('dist from outlet (m)')#, fontsize=9)
    #axis.set_xlabel('along-river distance (m)')#, fontsize=9)
    axis.set_ylabel('WSE (m)')#, fontsize=9)
    # Increase the fontsize of the tick labels for the primary axis
    #axis.tick_params(axis='both', which='major',
    #                 labelsize=9)  # You can adjust the fontsize as needed
    if title is not None:
        axis.set_title(title[0:20])

    if prd_heights:
        axis.plot(along_dist,#node_p_dist,
                data.nodes['p_wse'][node_i],
                'D', markersize=2, label='PRD wse')

    # Add reach summary metrics in text to bottom left of plot
    # if annotate_metrics:
    plot_summary_metrics(data, axis, reach_i)
    if truth is not None:
        # match and plot data. If there are no matches, capture that in
        # has_truth flag
        has_truth, pt_wse_e, pt_stdev, pt_slp_e = get_truth_matches(
            truth, reach_id, pass_no, tile, cycle, axis, fit_x, ss_min,
            ss_max, annotate_metrics, along_dist, node_id,
            p_dist_out, multi_reach)
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
    #leg = axis.legend(loc='best', fontsize=5, ncol=1)
    #leg2 = axis2.legend(fontsize=5, ncol=1)
    leg.set_draggable(1)
    # flip axis so highest part of river is on the left
    #axis.set_xlim(axis.get_xlim()[::-1])
    # call `draw` to re-render the graph
    plt.draw()
    #plt.tight_layout()

    return has_truth


def get_truth_matches(truth, reach_id, pass_no, tile, cycle, axis,
                      fit_x, ss_min, ss_max, annotate_metrics,
                      swot_along_dist, swot_node_id, p_dist_out,
                      multi_reach=False):
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
                        pt_node_match, pt_err_match, axis, fit_x,
                        ss_min, ss_max, reach_id, annotate_metrics,
                        swot_along_dist, swot_node_id
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
                    # using p_dist_out instead of cumcsum(p_length) here
                    # should be ok for within a reach, but may break for
                    # multireaches
                    d_dist_out = truth_drift['p_dist_out']
                    d_along_dist = np.max(p_dist_out) - d_dist_out
                    # TODO fix plotting here, i think y var is wrong
                    #      need to at least add the geoid so we plot wse
                    d_wse = truth_drift['ellipsoid_height_m']
                    axis.plot(d_along_dist,
                            d_wse, 'm+',
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
                    # using p_dist_out instead of cumcsum(p_length) here
                    # should be ok for within a reach, but may break for
                    # multireaches
                    d_along_dist = coarse_node_match['p_dist_out']
                    d_along_dist = np.max(p_dist_out) - d_along_dist
                    try:
                        axis.plot(along_dist,#coarse_node_match['p_dist_out'],
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

    return has_truth, pt_wse_e, pt_stdev, pt_slp_e

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
    reach_xtrk = str(round(np.mean(reach_xtrk) / 1000, 1))  # km
    reach_dark_frac = get_first_element(data.reaches['dark_frac'][reach_i])
    reach_obs_frac = get_first_element(data.reaches['obs_frac_n'][reach_i])
    reach_xovr_cal_q = get_first_element(data.reaches['xovr_cal_q'][reach_i])
    reach_ice_clim_f = get_first_element(data.reaches['ice_clim_f'][reach_i])
    reach_pwidth = get_first_element(data.reaches['p_width'][reach_i])
    reach_plength = get_first_element(data.reaches['p_length'][reach_i])
    # The ice flag, pwidth, and plength from SWORD can be fill valued, which
    # looks like a masked array here
    ice_clim_f = (str(round(reach_ice_clim_f, 2)) if not np.ma.is_masked(
        reach_ice_clim_f) else '-999')
    pwidth = (str(round(reach_pwidth, 2)) if not np.ma.is_masked(
        reach_pwidth) else '-999')
    plength = (str(round(reach_plength, 2)/1000) if not np.ma.is_masked(  # km
        reach_plength) else '-999')
    summary_string = 'p_width =' + pwidth + ' m\n' \
                     + 'p_length = ' + plength + ' km\n' \
                     + 'x-trk =' + reach_xtrk + ' km\n' \
                     + 'w = ' + str(round(reach_width, 2)) + ' m\n' \
                     + 'dark_frac = ' + str(round(reach_dark_frac, 2))  \
                     + '\nobs_frac = ' + str(round(reach_obs_frac, 2))\
                     + '\nxovr_q = ' + str(round(reach_xovr_cal_q, 2)) \
                     + '\n ice_clim_f = ' + ice_clim_f
    axis.text(LEFT, BOTTOM, summary_string,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=7,
              transform=axis.transAxes,
              bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

def plot_wse_and_qual(node_p_dist, wse, wse_r_u, node_q, node_q_b, axis,
                      plot_bit_qual, reach_wse, reach_slope, reach_slope2,
                      plot_slp2=False, mt_wse=None, ylim=None):
    # HACK to account for very bad wse_r_u values
    wse_r_u[wse_r_u < 0] = 0
    # Remove fill value nodes
    mask = wse >= -999
    node_p_dist, wse, wse_r_u, node_q, node_q_b = (arr[mask] for arr in
    [node_p_dist, wse, wse_r_u, node_q, node_q_b])

    # first plot the mt_reference and spread under everything else
    if mt_wse is not None:
        x = mt_wse.along_dist - np.min(mt_wse.along_dist)
        y = mt_wse.reference
        ptiles = mt_wse.percentiles
        ptile_list = mt_wse.percentile_list
        y_low, y_high, y_min, y_max = get_IQR_range(
            y, ptiles, ptile_list, p_low=5, p_high=95)
        if ylim is None:# dont overwrite input ylim
            ylim = (y_min, y_max)
        axis.plot(x, y, color='gray', label='mt reference')
        axis.fill_between(x,y_low, y_high, color='gray',
                label='mt 5-95%ile spread', alpha=0.2)
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
                      facecolor='cornflowerblue', alpha=0.3, interpolate=True,
                      label='3 x wse_r_u spread')

    # Enhanced reach slope
    if plot_slp2:
        obs_fit_y2 = [reach_wse + min(ss) * reach_slope2, reach_wse,
                      reach_wse + max(ss) * reach_slope2]
        axis.plot(fit_x, obs_fit_y2, '--', markersize=10, color='g',
                  label='slp2 fit')
    # 
    if ylim is not None:
        axis.set_ylim(ylim)
    #
    return fit_x, min(ss), max(ss)


def plot_pt_data(pt_node_match, pt_err_match, axis, fit_x, ss_min,
                 ss_max, reach_id, annotate_metrics, along_dist, node_id):
    # plot the pt data on the WSE axis
    pt_wse = pt_node_match.mean_node_pt_wse_m
    pt_node_id = get_simple_node_id(pt_node_match['node_id'], reach_id)
    pt_along_dist = []
    for nid in pt_node_id:
        pt_along_dist.append(along_dist[node_id==nid])
    pt_along_dist = np.array(pt_along_dist).squeeze()
    #breakpoint()
    pt_q = pt_node_match['pt_qual']

    # Define a custom colormap for PT node qual
    color_map = {0: 'green', 1: 'orange', 2: 'red'}

    colors = pt_q.map(color_map)
    # Plot the points with colors based on pt_qual
    scatter = axis.scatter(pt_along_dist, pt_wse, c=colors, marker='x',
                            label='PT node WSE', zorder=10, s=30)
    if annotate_metrics:
        #for node_id, wse, q in zip(pt_node_id, pt_wse, pt_node_match['flag']):
        #    axis.text(node_id, wse, str(q), fontsize=6, color='blue')
        for alng_d, wse, q in zip(pt_along_dist, pt_wse, pt_node_match['flag']):
            axis.text(alng_d, wse, str(q), fontsize=6, color='blue')
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
            axis.plot(np.mean(along_dist), pt_reach_wse, 'r*', markersize=8,
                       label='PT reach WSE', zorder=0)
            if pt_slope is not None:
                pt_fit_y = [pt_reach_wse - ss_min * pt_slope, pt_reach_wse,
                            pt_reach_wse - ss_max * pt_slope]
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
              multi_reach=False, mt_width=None, plot_width=True, ylim=None):
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
        #
        if mt_width is not None:
            mt_width = mt_width.crop_to_reach()
    node_id = data.nodes['node_id'][node_i]
    node_id = get_simple_node_id(node_id, reach_id)
    along_dist = np.cumsum(data.nodes['p_length'][node_i])
    along_dist = np.max(along_dist) - along_dist # make it go downstream
    if plot_width:
        #x_detct = data.nodes['width_detct'][node_i]
        x_total = data.nodes['width'][node_i]
        area_total = data.nodes['area_total'][node_i]
        area_detct = data.nodes['area_detct'][node_i]
        x_detct = area_detct / area_total * x_total # compute detected width
    else:
        x_detct = data.nodes['area_detct'][node_i]
        x_total = data.nodes['area_total'][node_i]  # includes dark water pixels

    #ylim = None
    if mt_width is not None:
        x = mt_width.along_dist - np.min(mt_width.along_dist)
        y = mt_width.reference
        ptiles = mt_width.percentiles
        ptile_list = mt_width.percentile_list
        y_low, y_high, y_min, y_max = get_IQR_range(
                y, ptiles, ptile_list, p_low=5, p_high=95)
        if ylim is None: # dont overwrite input ylim
            ylim = (y_min, y_max)
        axis.plot(x, y, color='gray', label='mt reference')
        axis.fill_between(x,y_low, y_high, color='gray', 
                label='mt 5-95%ile spread', alpha=0.2)
        if y_min < -5.0: # clip width min axis lim
            y_min = -5.0
    axis.scatter(along_dist, x_detct, c='k', marker='x', alpha=.7,
            label='detected')
    scat = axis.scatter(along_dist, x_total, c=np.array(node_id),
            alpha=.7, cmap='tab20b', label = 'total (color=node_id)')
    #cbar = plt.colorbar(scat, ax=axis)
    #cbar.set_label('node_id (local)')
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
            along_truth = np.cumsum(truth.nodes['p_length'][node_i_truth])
            along_truth = np.max(along_truth) - along_truth
            if plot_width:
                x_truth = truth.nodes['width'][node_i_truth]
            else:
                x_truth = truth.nodes['area_total'][node_i_truth]
            #axis.plot(node_id_truth, x_truth, 'kx', markersize=2)
            axis.plot(along_truth, x_truth, 'kx', markersize=2)
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
    axis.set_xlabel('along-river distance (m)')
    if plot_width:
        axis.set_ylabel('width (m)')
    else:
        axis.set_ylabel('area (m^2)')
    #leg = axis.legend(['detected', 'total (color=node_id)', 'truth'], fontsize=5)
    #leg = axis.legend(fontsize=5)
    leg = axis.legend(
        bbox_to_anchor=(1.01, 1.0), loc='upper left', fontsize=5, ncol=1)
    leg.set_draggable(1)
    #axis.set_xlim(axis.get_xlim()[::-1])  # flip axis to align with WSE plot
    # 
    if ylim is not None:
        axis.set_ylim(ylim)
    #
    if title is not None:
        axis.set_title(title)

def toslant(pixc, varname):
    data = pixc[varname]
    var = np.ma.zeros((
        pixc.interferogram_size_azimuth,
        pixc.interferogram_size_range),
        dtype=data.dtype)
    var[var==0] = np.ma.masked
    var[pixc.azimuth_index, pixc.range_index] = data
    return var

def plot_pix_assgn(data, reach_id, axis, h_flg=False, area_flg=False,
                   pixc_data=None, var='node_id', multi_reach=False,
                   plot_map=True, mt_data=None, clim=None):
    # Filter data for the specified reach_id
    if multi_reach:
        # grab adjacent reaches too if they exist
        pix_i = np.logical_or.reduce((
            data['reach_id'] == reach_id,
            data['reach_id'] == reach_id - 10,
            data['reach_id'] == reach_id + 10
        ))
    else:
        #breakpoint()
        pix_i = (data['reach_id'] == reach_id)
        if mt_data is not None:
            mt_data = mt_data.crop_to_reach()
    node_id = data['node_id'][pix_i]
    node_id = get_simple_node_id(node_id, reach_id)
    #breakpoint()
    lat = data['latitude_vectorproc'][pix_i]
    lon = data['longitude_vectorproc'][pix_i]
    height_vec = data['height_vectorproc'][pix_i]
    range_index_vec = data['range_index'][pix_i]
    azimuth_index_vec = data['azimuth_index'][pix_i]
    if pixc_data is not None:
        # get the corresponding pixc pixels
        #pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc)
        #def toslant(pixc, varname):
        #    data = pixc[varname]
        #    var = np.ma.zeros((
        #        pixc.interferogram_size_azimuth,
        #        pixc.interferogram_size_range),
        #        dtype=data.dtype)
        #    var[var==0] = np.ma.masked
        #    var[pixc.azimuth_index, pixc.range_index] = data
        #    return var
        if var=='wse':
            pixc_var = load_wse_data(pixc_data)
            """
            p_height = toslant(pixc_data.pixel_cloud, 'height')
            p_geoid = toslant(pixc_data.pixel_cloud, 'geoid')
            p_solid = toslant(pixc_data.pixel_cloud, 'solid_earth_tide')
            p_load = toslant(pixc_data.pixel_cloud, 'load_tide_fes')
            p_pole = toslant(pixc_data.pixel_cloud, 'pole_tide')
            pixc_var = p_height - (p_geoid + p_solid + p_load + p_pole)
            """
        else:
            pixc_var = toslant(pixc_data.pixel_cloud, var)
        var_pixc = pixc_var[azimuth_index_vec, range_index_vec]
    else:
        var_pixc = height_vec
        if var=='wse':
            var = 'height'
    #breakpoint()
    # Create a scatter plot manually handling colors
    c_var = node_id
    cmap = 'tab20b'
    c_label = 'node_id (local)'
    clim_ptile=False
    if var=='height' or var=='wse':
        c_var = var_pixc
        cmap = 'jet'
        if pixc_data is None:
            c_label = 'pixcvec {} (m)'.format(var)
        else:
            c_label = 'pixc {} (m)'.format(var)
        clim_ptile = True
    #if var=='wse':
    #    c_var = var_pixc
    #    cmap = 'jet'
    #    c_label = 'wse (m)'
    #    clim_ptile = True
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
    # get min and max for colorbar
    if not np.ma.is_masked(c_var) and np.ma.count(c_var) > 0:
        vmin = c_var.min()
        vmax = c_var.max()
        if clim_ptile:
            vmin = np.nanpercentile(c_var, 5)
            vmax = np.nanpercentile(c_var, 95)
    else:
        # Set a default minimum value if the array is empty or masked
        vmin = 0
        vmax = 1
    if clim is None:
        clim = (vmin, vmax)
        if (mt_data is not None) and var=='wse':
            # overwrite using multitemporal IQR
            y = mt_data.reference
            ptiles = mt_data.percentiles
            ptile_list = mt_data.percentile_list
            y_low, y_high, y_min, y_max = get_IQR_range(
                y, ptiles, ptile_list)
            clim = (y_min, y_max)
    """
    # Create a scatter plot manually handling colors
    c_var = node_id
    cmap = 'tab20b'
    c_label = 'node_id (local)'
    clim_ptile=False
    if var=='height':
        c_var = var_pixc
        cmap = 'jet'
        c_label = 'height (m)'
        clim_ptile = True
    """
    lt = lat.copy()
    ln = lon.copy()
    if mt_data is not None:
        # replace bounding box using the mt lat/lon so whole stack of cases
        # get cropped the same
        lt = mt_data.p_lat
        ln = mt_data.p_lon
    # define the bounding box
    buff = 0.01 #
    bbox = [
            np.nanmin(ln) - buff,
            np.nanmax(ln) + buff,
            np.nanmin(lt) - buff,
            np.nanmax(lt) + buff]
    # explicitly set default transform
    transform = axis.transData
    default_tx = axis.transData
    if plot_map:
        """
        # define the bounding box
        buff = 0.01 #
        bbox = [
            np.nanmin(lon) - buff,
            np.nanmax(lon) + buff,
            np.nanmin(lat) - buff,
            np.nanmax(lat) + buff]
        """
        # set zoom based on bbox
        zoom = 14
        sz = np.sqrt((bbox[1] - bbox[0])**2 + (bbox[2] - bbox[3])**2)
        print(sz)
        # decrement zooom at 1/8 degree increments
        zoom = int(16 - np.ceil(sz*8))
        if zoom >14:
            zoom = 14
        if zoom < 6:
            zoom = 6
        print("bbox size: ", sz, ", zoom: ", zoom)
        # plot the google image
        tiler = GoogleTiles(style="satellite")
        mercator = tiler.crs
        transform=ccrs.Geodetic()
        #axis.set_projection(mercator)
        axis.set_extent(bbox, crs=ccrs.Geodetic())#ccrs.PlateCarree())
        axis.add_image(tiler, zoom, alpha=0.7)
        axis.add_artist(ScaleBar(1, location='lower right'))

    # plot the data
    scatter = axis.scatter(lon, lat, c=c_var, cmap=cmap, #norm=norm,
                           linewidth=0, alpha=0.7, s=5, clim=clim,
                           transform=transform)
    # plot swath orientation arrows
    #breakpoint()
    # inner along
    d_lat_ai = (data.inner_first_latitude - data.inner_last_latitude)
    d_lon_ai = (data.inner_first_longitude - data.inner_last_longitude)
    # outer along
    d_lat_ao = (data.outer_first_latitude - data.outer_last_latitude)
    d_lon_ao = (data.outer_first_longitude - data.outer_last_longitude)
    # first cross
    d_lat_xf = (data.inner_first_latitude - data.outer_first_latitude)
    d_lon_xf = (data.inner_first_longitude - data.outer_first_longitude)
    # last cross
    d_lat_xl = (data.inner_last_latitude - data.outer_last_latitude)
    d_lon_xl = (data.inner_last_longitude - data.outer_last_longitude)
    # inner/outer average
    d_lat_a = (d_lat_ai + d_lat_ao) / 2.0
    d_lon_a = (d_lon_ai + d_lon_ao) / 2.0
    # first/last average
    d_lat_x = (d_lat_xf + d_lat_xl) / 2.0
    d_lon_x = (d_lon_xf + d_lon_xl) / 2.0
    # compute angles
    ang_a = np.arctan2(d_lat_a, d_lon_a)
    ang_x = np.arctan2(d_lat_x, d_lon_x)
    # define where to put on plot
    lat_dist = (bbox[3] - bbox[2])
    lon_dist = (bbox[1] - bbox[0])
    scale = 0.2
    origin_lon = bbox[0] + scale * lon_dist
    origin_lat = bbox[3] - scale * lat_dist
    arrow_len = 0.5 * scale * np.sqrt(lat_dist**2 + lon_dist**2)
    #
    xy_origin = np.array((origin_lon, origin_lat))
    xy_head_a = np.array((origin_lon + arrow_len * np.cos(ang_a),
        origin_lat + arrow_len * np.sin(ang_a)))
    xy_head_x = np.array((origin_lon + arrow_len * np.cos(ang_x),
        origin_lat + arrow_len * np.sin(ang_x)))

    # put the arrows on the plot
    axis.annotate("",#"along-track",
            xy=xy_origin, #xycoords=coords,
            xytext=xy_head_a, #textcoords=transform,
            arrowprops=dict(arrowstyle="<-"),
            ha="center", va="center",#rotation=ang_a*180/np.pi,
            transform=transform,
            size=6
            )
    """
    axis.text(
            xy_origin[0], xy_origin[1], #xycoords=coords,
            "along-track",
            rotation=np.rad2deg(ang_a),
            transform=transform,
            size=6
            )
    """
    axis.annotate("",#"cross-track",
            xy=xy_origin, #xycoords=default_tx,
            xytext=xy_head_x, #textcoords=default_tx,
            arrowprops=dict(arrowstyle="<-"),
            ha="center", va="center",#rotation=ang_a*180/np.pi,
            transform=transform,
            size=6
            )

    if not plot_map:
        # Set plot properties
        axis.grid(True)
        # plt.gca().set_aspect('equal', adjustable='box')
        axis.set_xlabel('lon')
        axis.set_ylabel('lat')

    # Add colorbar
    #colorbar = plt.colorbar(sm, ax=axis)
    colorbar = plt.colorbar(scatter,ax=axis)
    colorbar.set_label(c_label)

    if h_flg:
        # plot bad h_flg over existing plot, if it is available
        try:
            h_flg = data['h_flg'][pix_i]
            bad_wse_lat = lat[h_flg == 0]
            bad_wse_lon = lon[h_flg == 0]
            axis.scatter(bad_wse_lon, bad_wse_lat, color='k', s=0.1,
                    #label='h_flg=False',
                    transform=transform)
            #axis.set_title('Pixel Locations, h_flg/node_ID')
            axis.plot([], [], 'ko', label='h_flg=False')
            #axis.legend()
        except AttributeError:
            print('h_flg not in output pixcvecriver.')
    elif area_flg:
        # plot bad area_flg over existing plot, if it is available
        try:
            area_flg = data['area_flg'][pix_i]
            bad_area_lat = lat[area_flg == 0]
            bad_area_lon = lon[area_flg == 0]
            axis.scatter(bad_area_lon, bad_area_lat, color='k', s=0.1,
                    #label='area_flg=False',
                    transform=transform)
            axis.set_title('Pixel Locations, area_flg')
            #axis.plot([], [], 'ko', label='area_flg=False')
        except AttributeError:
            print('area_flg not in output pixcvecriver.')
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
        pixcvec_data, river_data, pixc_data,
        truth_pixcvec, truth_pixc, reach_id,
        nodes=None, pixc_truth=None, out_dir=None, title=None, has_truth=None):

    plt.tight_layout()
    mngr = plt.get_current_fig_manager()
    # mngr.window.setGeometry(0, 0, 1500, 500)
    if pixc_data and pixcvec_data:
        #pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc)
        pixc_truth_data = None
        if pixc_truth is not None:
            pixc_truth_data = SWOTWater.products.product.MutableProduct.from_ncfile(
                pixc_truth)
        slant_plane_fig = plot_pixcs(
            pixcvec_data, pixc_data, reach_id, nodes, reach_data=river_data,
            pixc_truth=pixc_truth_data
        )
        if title is None:
            title=''
        this_river_name = river_data['reaches']['river_name'][
            river_data['reaches']['reach_id'] == reach_id][0]
        date = pixc_data.attributes['time_granule_start'].split('T')[0]
        fig_title = '{} {} {} '.format(
            this_river_name, reach_id, date) + title
        #fig_title = this_river_name + ' ' + date + ' ' + title
        slant_plane_fig.suptitle(fig_title)
        if out_dir is not None:
            # save current figure to file
            plt.title(title, backgroundcolor='white')
            base_fname = date + '_' + title + '_PIXC'
            base_fname = base_fname.replace(':','').replace('.','p').replace(
                '-','_').replace(' ','_')
            if has_truth:
                filename = out_dir + '/truth/' + base_fname
            else:
                filename = out_dir + '/no_truth/' + base_fname
            #this_river_name = river_data['reaches']['river_name'][
            #    river_data['reaches']['reach_id'] == reach_id][0]
            #fig_title = this_river_name + ' ' + date + ' ' + title
            #slant_plane_fig.suptitle(fig_title)
            slant_plane_fig.savefig(filename)
            plt.close()
        else:
            plt.title(title, backgroundcolor='white')
            plt.show()

    else:
        print('Missing pixc or pixcvec file, skipping pixel assignment plot')

    if pixc_data and truth_pixc:  # only plot these if pixc was also given
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
                     errors=None, out_dir=None, title=None, multi_reach=False,
                     pixc_data=None, mt_wse=None, mt_width=None, plot_map=True):
    # contains node group and reach group for each input netcdf
    #cycle, pass_no, tile = decode_rivertile_filename(rivertile_file)
    cycle = '{}'.format(river_data.nodes.cycle_number)
    pass_no = '{}'.format(river_data.nodes.pass_number)
    tile = '{}'.format(river_data.nodes.tile_number)
    #breakpoint()
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

    #figure, axes = plt.subplots(2, 2, figsize=FIGSIZE, dpi=DPI)
    figure = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax1 = plt.subplot(221)
    ax3 = plt.subplot(223, sharex = ax1)
    if plot_map:
        tiler = GoogleTiles(style="satellite")
        mercator = tiler.crs
        ax2 = plt.subplot(222, projection=mercator)
        ax4 = plt.subplot(224, projection=mercator)
    else:
        ax2 = plt.subplot(222)
        ax4 = plt.subplot(224)
    axes = np.array([[ax1, ax2],[ax3, ax4]])
    #breakpoint()
    date = river_data.reaches.time_granule_start.split('T')[0] 
    #breakpoint()
    river_name = river_data.reaches['river_name'][
            river_data.reaches['reach_id']==reach_id][0]
    ttl = '{} {} {},  pass: {}, cycle: {}'.format(
            river_name, reach_id, date, pass_no, cycle)
    #river_name + ' ' + str(reach_id) + ' {}'.format(
    #        date) + ', pass: '', cycle: ' + str(cycle)
    figure.suptitle(ttl)
    has_truth = plot_wse(
        river_data, truth, errors, reach_id, axes[0][0], figure,
        #title=title_str + ' - wse',
        cycle=cycle, tile=tile, pass_no=pass_no,
        multi_reach=multi_reach, mt_wse=mt_wse
    )
    plot_area(river_data, truth, errors, reach_id, axes[1][0],
              #title=title_str + ' - width',
              multi_reach=multi_reach, mt_width=mt_width)
    # uncomment the below block if you'd prefer to plot the centroids rather
    # than the area flag
    # plot_locations(river_data, truth, reach_id, axes[0][1],
    #                title=title_str + ' - locations')
    if pixcvec is not None:
        plot_pix_assgn(pixcvec, reach_id, axes[0][1], h_flg=True,
                       pixc_data=pixc_data, var='wse',#var='height',
                       multi_reach=multi_reach, plot_map=plot_map,
                       mt_data=mt_wse)
        #plot_pix_assgn(pixcvec, reach_id, axes[0][1], area_flg=True,
        #               multi_reach=multi_reach)
        plot_pix_assgn(pixcvec, reach_id, axes[1][1], h_flg=True,
                       multi_reach=multi_reach, plot_map=plot_map)
    plt.tight_layout()
    #breakpoint()
    if out_dir is not None:
        # save current figure to file
        river_code = get_river_code(river_data, reach_id)
        #plt.title(title, backgroundcolor='white')
        base_fname = date + '_' + title + '_' + river_code
        base_fname = base_fname.replace(':','').replace('.','p').replace(
                '-','_').replace(' ','_')
        if has_truth:
            filename = out_dir + '/truth/' + base_fname
        else:
            filename = out_dir + '/no_truth/' + base_fname
        #breakpoint()
        plt.savefig(filename)
        plt.close()
    else:
        plt.title(title, backgroundcolor='white')
        plt.show()
    return figure, axes, has_truth


def make_plots(rivertile_file, rivertile_df, truth_data, pixcvec, pixc,
               truth_pixcvec, truth_pixc, reach_id, errors=None,
               nodes=None, pixc_truth=None, mt_wse=None, mt_width=None,
               out_dir=None, title=None, overwrite=False, sandbox=False):

    # handle overwriting if user says not to
    make_rivers = True
    make_pixc = True
    if (~overwrite) and (out_dir is not None):
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
            # handle pixcvec reach id that is different than pixcvec river
            #breakpoint()
            if (pixcvec_data['reach_id'].dtype == '|S1'):
                rid = pixcvec_data['reach_id']
                rid_int = np.array([
                    np.nan if np.sum(k.mask)>0 else int(k.tostring()) for k in rid])
                nid = pixcvec_data['node_id']
                nid_int = np.array([
                    np.nan if np.sum(k.mask)>0 else int(k.tostring()) for k in nid])
                #pixcvec_data['reach_id'] = rid_int
                ATTRS = pixcvec_data.ATTRIBUTES
                VARS = pixcvec_data.VARIABLES
                DIMS = pixcvec_data.DIMENSIONS
                from collections import OrderedDict as odict
                VARS['reach_id'] = odict([
                    ('dtype', 'i8'),
                    ('dimensions', odict([('points', 0),])),
                    ('coordinates', 'longitude_vectorproc latitude_vectorproc'),
                    ('comment', '')])
                VARS['node_id'] = odict([('dtype', 'i8'),
                    ('dimensions', odict([('points', 0),])),
                    ('coordinates', 'longitude_vectorproc latitude_vectorproc'),
                    ('comment', '')])
                #breakpoint()
                #make a copy of the pixcvec and replace the reach_id
                pixcvec_data2 = SWOTWater.products.product.MutableProduct(
                        attributes=ATTRS,
                        variables=VARS,
                        dimensions = DIMS)
                
                for key in pixcvec_data.attributes.keys():
                    pixcvec_data2[key] = pixcvec_data[key]
                var_keys = list(set(pixcvec_data.variables.keys()) - 
                        set(['reach_id',]) - set(['node_id',]))
                for key in var_keys:
                    pixcvec_data2[key] = pixcvec_data[key]
                pixcvec_data2['reach_id'] = rid_int
                pixcvec_data2['node_id'] = nid_int
                #breakpoint()
                
                # now overwrite the original
                pixcvec_data = pixcvec_data2
        else:
            pixcvec_data = None
        if pixc is not None:
            pixc_data = SWOTWater.products.product.MutableProduct.from_ncfile(pixc)
        else:
            pixc_data = None
    if sandbox:
        sandbox_run(rivertile_df, pixcvec_data, pixc_data, reach_id)
    if make_rivers:
        fig, ax, has_truth = make_river_plots(
            rivertile_file, rivertile_df, truth_data, pixcvec_data, reach_id,
            errors=errors, out_dir=out_dir, title=title, pixc_data=pixc_data,
            mt_wse=mt_wse, mt_width=mt_width)
    if make_pixc:
        make_pixc_plots(
            pixcvec_data, rivertile_df, pixc_data, truth_pixcvec, truth_pixc, reach_id,
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
               apply_corr=True):
    # Plots six pixel cloud plots in a grid, each showing different information.
    # Windows the pixel cloud to the pixels in PIXCVecRiver (the river assigned
    # ones) for the input reach ID.
    reach_id = int(reach_id)
    pix_i = (pixc_vec['reach_id'] == reach_id)


    # If no data for this reach, bail out
    if np.sum(pix_i) == 0:
        print(f'No reach ID matching {reach_id} in this pixel cloud')
        return None

    slant_plane_fig = plt.figure(figsize=FIGSIZE, dpi=DPI)

    # Get pixc_vec arrays for just this reach
    node_id0 = pixc_vec['node_id'][pix_i]
    aziv = pixc_vec['azimuth_index'][pix_i]
    riv = pixc_vec['range_index'][pix_i]
    heightv = pixc_vec['height_vectorproc'][pix_i]

    # Grab entire pixc arrays
    azi_full = pixc.pixel_cloud['azimuth_index']
    ri_full = pixc.pixel_cloud['range_index']

    # Determine array sizes to hold *all* pixc AND pixc_vec indices
    #    (no immediate bounding box clamp)
    msize = max(np.max(azi_full), np.max(aziv)) + 1
    nsize = max(np.max(ri_full), np.max(riv)) + 1

    Node_id_full = np.full((msize, nsize), np.nan, dtype=float)
    Heightv_full = np.full((msize, nsize), np.nan, dtype=float)
    node_id_simpl = np.array([get_simple_node_id(nid, reach_id) for nid in node_id0])
    in_bounds = (
        (aziv >= 0) & (aziv < msize) &
        (riv >= 0)  & (riv < nsize)
    )
    Node_id_full[aziv[in_bounds], riv[in_bounds]] = node_id_simpl[in_bounds]
    Heightv_full[aziv[in_bounds], riv[in_bounds]] = heightv[in_bounds]

    # now get PIXC in slant-plane
    pixc_arrays = {
        'height': pixc.pixel_cloud['height'],
        'geoid': pixc.pixel_cloud['geoid'],
        'solid_earth_tide': pixc.pixel_cloud['solid_earth_tide'],
        'load_tide_fes': pixc.pixel_cloud['load_tide_fes'],
        'load_tide_got': pixc.pixel_cloud['load_tide_got'],
        'pole_tide': pixc.pixel_cloud['pole_tide'],
        'classification': pixc.pixel_cloud['classification'],
        'water_frac': pixc.pixel_cloud['water_frac'],
        'pixel_area': pixc.pixel_cloud['pixel_area'],
    }
    # Initialize 2D grids for each variable
    arrays_2d = {}
    for key, values in pixc_arrays.items():
        arr2d = np.full((msize, nsize), np.nan, dtype=values.dtype)
        arr2d[azi_full, ri_full] = values
        arrays_2d[key] = arr2d

    # Determine a bounding box to "zoom in" on the region of this reach
    M0 = np.min(aziv[in_bounds])
    M1 = min(np.max(aziv[in_bounds]), np.max(azi_full)) + 1
    N0 = np.min(riv[in_bounds])
    N1 = min(np.max(riv[in_bounds]), np.max(ri_full)) + 1

    Node_id = Node_id_full[M0:M1, N0:N1]
    Heightv = Heightv_full[M0:M1, N0:N1]

    # Crop everything to match the pixc_vec bounding box (M0:M1, N0:N1)
    crop_arrays = {}
    for key, arr2d in arrays_2d.items():
        crop_arrays[key] = arr2d[M0:M1, N0:N1]

    # If we have truth classifications, slice that too
    cls_t = None
    if pixc_truth is not None:
        cls_t = pixc_truth['classification'][M0:M1, N0:N1]
    # Mask out anything not valid for this reach
    # i.e., places where heightv or Node_id are NaN
    isnan_heightv = np.isnan(Heightv)
    isnan_nodeid = np.isnan(Node_id)
    for k in ('height', 'geoid', 'solid_earth_tide', 'load_tide_fes',
              'load_tide_got', 'pole_tide'):
        crop_arrays[k][isnan_heightv] = np.nan
    crop_arrays['classification'] = crop_arrays['classification'].astype(float)
    crop_arrays['classification'][isnan_nodeid] = np.nan
    crop_arrays['water_frac'][isnan_nodeid] = np.nan
    crop_arrays['pixel_area'][isnan_nodeid] = np.nan
    # Create a water-area array for display
    Warea1 = crop_arrays['pixel_area'].copy()
    # Remove classes 2, 1, and 3 from area display
    for badclass in [2, 1, 3]:
        Warea1[crop_arrays['classification'] == badclass] *= 0
    # Apply geophysical corrections if specified by user to do so
    if apply_corr:
        correction = (crop_arrays['geoid'] +
                      crop_arrays['solid_earth_tide'] +
                      crop_arrays['load_tide_fes'] +
                      crop_arrays['pole_tide'])
        crop_arrays['height'] -= correction
        Heightv -= correction
    # Compute color ranges for height
    cmap_min = np.nanpercentile(crop_arrays['height'], 20)
    cmap_max = np.nanpercentile(crop_arrays['height'], 80)

    # Top left: Plot node assignments
    ax1 = plt.subplot(2, 3, 1)
    pt1 = ax1.imshow(Node_id, interpolation='none', aspect='auto',
                     cmap=plt.cm.get_cmap('tab20b'))
    colorbar = plt.colorbar(pt1, ax=ax1)
    ax1.set_title('node_id ' + title_tag)
    colorbar.set_label('Node ID')

    # Top center: PIXC water classifications
    ax2 = plt.subplot(2, 3, 2, sharex=ax1, sharey=ax1)
    class_cmap = colors.ListedColormap(
        ['pink', 'darkgreen', 'lightgreen', 'aquamarine', 'blue', 'black',
         'yellow', 'red']
    )
    pt2 = ax2.imshow(crop_arrays['classification'], interpolation='none',
                     aspect='auto', cmap=class_cmap, clim=(0, 7))
    ax2.set_title('classification ' + title_tag)
    colorbar = plt.colorbar(pt2, ax=ax2)
    colorbar.set_label('PIXC water class')


    # Top and bottom right: Area plots. Have to handle detected and total.
    ax3 = plt.subplot(2, 3, 3, sharex=ax1, sharey=ax1)
    if reach_data is not None:
        NodeArea = Node_id.copy()      # total area
        NodeArea_det = Node_id.copy()  # detected area
        node_i = np.logical_and(
            reach_data.nodes['reach_id'] == reach_id,
            ~reach_data.nodes['area_total'].mask
        )
        node_ids = reach_data.nodes['node_id'][node_i]
        area_tot = reach_data.nodes['area_total'][node_i]
        area_det = reach_data.nodes['area_detct'][node_i]
        p_width = reach_data.nodes['p_width'][node_i]

        a_cmax = np.max(p_width) * 1.1
        a_cmin = np.min(p_width) * 0.9
        for raw_nid in node_id0:
            sid = get_simple_node_id(raw_nid, reach_id)
            # if area data exist
            if len(area_tot[node_ids == raw_nid]) > 0:
                NodeArea[Node_id == sid] = area_tot[node_ids == raw_nid]
                NodeArea_det[Node_id == sid] = area_det[node_ids == raw_nid]

        NodeArea = NodeArea / 200
        NodeArea_det = NodeArea_det / 200

        pt3 = ax3.imshow(NodeArea, interpolation='none', aspect='auto',
                         cmap='jet', clim=(a_cmin, a_cmax))
        ax3.set_title('Node Area (m^2)' + title_tag)
        colorbar = plt.colorbar(pt3, ax=ax3)
        colorbar.set_label('Node Area (m^2)')


        abins0 = np.linspace(50, 450, 100)
        amsk = NodeArea > -100
        ha, abins = np.histogram(NodeArea[amsk], abins0)
        amn = np.mean(NodeArea[amsk])
        amed = np.median(NodeArea[amsk])
        asd = np.std(NodeArea[amsk])

        ax6 = plt.subplot(2, 3, 6, sharex=ax1, sharey=ax1)
        pt6 = ax6.imshow(NodeArea_det, interpolation='none', aspect='auto',
                         cmap='jet', clim=(a_cmin, a_cmax))  # , clim=(c0,c1))
        ax6.set_title('Node Area det. ' + title_tag)
        colorbar = plt.colorbar(pt6, ax=ax6)
        colorbar.set_label('Detected Node Area (m^2)')

    else:
        # plot pixel level water area
        wmax = np.nanpercentile(Warea1, 90)
        pt3 = ax3.imshow(Warea1, interpolation='none', aspect='auto',
                         cmap='jet', clim=(0, wmax))
        ax3.set_title('water area (pixel-level) ' + title_tag)
        colorbar = plt.colorbar(pt3, ax=ax3)
        colorbar.set_label('Pixel water area (m^2)')

        ax6 = plt.subplot(2, 3, 6, sharex=ax1, sharey=ax1)
        pt6 = ax6.imshow(Geoid1, interpolation='none', aspect='auto',
                         cmap=cmaph)  # , clim=(c0,c1))
        ax6.set_title('geoid height (m) ' + title_tag)
        colorbar = plt.colorbar(pt6, ax=ax6)
        colorbar.set_label('geoid height (m)')

    # Bottom left: PIXCVecRiver height. Aggregated node heights.
    ax4 = plt.subplot(2, 3, 4, sharex=ax1, sharey=ax1)
    pt4 = ax4.imshow(Heightv, interpolation='none', aspect='auto',
                     cmap=cmaph, clim=(cmap_min, cmap_max))
    ax4.set_title('height_vectorproc (m) ' + title_tag)
    colorbar = plt.colorbar(pt4, ax=ax4)
    colorbar.set_label('PIXCVecRiver Height, (m)')

# Bottom center: Pixel cloud heights
    ax5 = plt.subplot(2, 3, 5, sharex=ax1, sharey=ax1)
    pt5 = ax5.imshow(crop_arrays['height'], interpolation='none', aspect='auto',
                     cmap=cmaph, clim=(cmap_min, cmap_max))
    ax5.set_title('height (m) ' + title_tag)
    colorbar = plt.colorbar(pt5, ax=ax5)
    colorbar.set_label('PIXC height, (m)')

    # Plot an extra set of figures for PIXC truth classifications, if we have
    # them.
    if cls_t is not None:
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

        # Classification differences
        Cls_t2 = np.zeros(np.shape(Cls_t))
        Cls_t2[np.logical_and(Cls_t > 0, Cls1 > 0)] = 3
        Cls_t2[np.logical_and(Cls_t > 0, np.isnan(Cls1))] = 2
        Cls_t2[np.logical_and(Cls1 > 0, np.isnan(Cls_t))] = 1
        ax_3 = plt.subplot(2, 3, 3, sharex=ax_1, sharey=ax_1)
        pt_3 = ax_3.imshow(Cls_t2, interpolation='none', aspect='auto',
                           cmap=plt.cm.get_cmap('tab10'), clim=(0, 5))
        plt.colorbar(pt_3, ax=ax_3)
        ax_3.set_title('classification ' + title_tag)


        ax_4 = plt.subplot(2, 3, 4, sharex=ax_1, sharey=ax_1)
        pt_4 = ax_4.imshow(Cls1, interpolation='none', aspect='auto',
                           cmap=plt.cm.get_cmap('tab10'), clim=(0, 5))
        plt.colorbar(pt_4, ax=ax_4)
        ax_4.set_title('classification diff ' + title_tag)

        ax_5 = plt.subplot(2, 3, 5, sharex=ax_1, sharey=ax_1)
        pt_5 = ax_5.imshow(Area_t, interpolation='none', aspect='auto',
                           cmap='jet', clim=(a_cmin, a_cmax))
        plt.colorbar(pt_5, ax=ax_5)
        ax_5.set_title('NodeArea pixc_true ' + title_tag)

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

        # Plot node level histograms, if specified by user to do so
        if nodes:
            for node in nodes:
                # plot node-level pixc height histograms
                idx = (Node_id == int(node))
                hgt = Height1[idx]
                hgtv = Heightv[idx]
                klass = Cls1[idx]

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
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.set_xlabel("Range index")
        ax.set_ylabel("Azimuth index")
    return slant_plane_fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_tile', help='river_data/rivertile.nc')
    parser.add_argument('reach_id', help='reach id', type=int)
    parser.add_argument('--truth_tile', help='river_data/rivertile.nc or pickle file dir',
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
    parser.add_argument('--mt_wse',
            help='multitemporal <reach_id>_wse_satas.nc file', default=None)
    parser.add_argument('--mt_width',
            help='multitemporal <reach_id>_width_satas.nc file', default=None)
    parser.add_argument('--sandbox', help='apply experimental algorithms',
            action='store_true', default=False)
    parser.add_argument('--cycle_id',
                        help='which cycle to plot (if using hydrochron nodes)',
                        default=None)
    parser.add_argument('--pass_id',
                        help='which pass to plot (if using hydrochron nodes)',
                        default=None)
    args = parser.parse_args()

    proc_tile = os.path.abspath(args.proc_tile)
    truth_tile = None
    truth_pixcvec = None
    truth_pixc = None
    gdem_dem = None
    reach_error = None
    if args.truth_tile is not None:  # user wants to plot truth as well
        if os.path.isfile(args.truth_tile):
            truth_tile = os.path.abspath(args.truth_tile)
            gdem_dem = get_gdem_from_rivertile(args.proc_tile)
            truth_pixcvec = args.truth_pixcvec
            errors = get_errors(proc_tile, truth_tile, test=False,
                                truth_filter=None)
            reach_error = get_reach_error(errors, args.reach_id)
        else:
            # assume it is a pkl file
            print("trying to load pkl file")
            truth_tile = load_pkl_dataframes(args.truth_tile)
            print("finished reading pkl file")
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
    mt_wse = None
    if args.mt_wse is not None:
        mt_wse = rivscale.products.along_stretch.AlongStretchStats.from_ncfile(
                args.mt_wse)
    mt_width = None
    if args.mt_width is not None:
        mt_width = rivscale.products.along_stretch.AlongStretchStats.from_ncfile(
                args.mt_width)
    if os.path.isfile(proc_tile):
        # read the dataframe
        if proc_tile.endswith('.nc'):
            proc_df = SWOTRiver.products.rivertile.L2HRRiverTile.from_ncfile(proc_tile)
        elif proc_tile.endswith('.shp'):
            # assume it is RiverSP shape-file
            proc_df = SWOTRiver.products.rivertile.L2HRRiverTile()
            print("reading RiverSP file:", proc_tile)
            node_df = gpd.read_file(proc_tile, ignore_geometry=True)
            node_df = node_df[node_df['reach_id']=='{}'.format(args.reach_id)]
            # make sure node_id and reach_id are ints
            node_df['reach_id'] = np.array(
                    [int(rid) for rid in node_df['reach_id']])
            node_df['node_id'] = np.array(
                    [int(nid) for nid in node_df['node_id']])
            # sort node_id
            node_df = node_df.sort_values(['node_id'])
            #i get cycle and pass from filename
            parts = os.path.split(proc_tile)[1].split('_')
            cycle_id = parts[5]
            pass_id = parts[6]
            # populate the rivertile object
            proc_df = populate_rivertile(proc_df, node_df, cycle_id, pass_id)
            #breakpoint()
            #
        else:
            proc_df = SWOTRiver.products.rivertile.L2HRRiverTile()
            # assume it is a csv dataframe output from hydrochron
            node_df = pd.read_csv(proc_tile)
            # filter out the desired reach
            node_df = node_df[node_df['reach_id']==args.reach_id]
            # filter out the time/cycle obs based on the pixcvec input
            cycle_id = args.cycle_id
            pass_id = args.pass_id
            if cycle_id is None:
                # try to get from the opixc of pixcvec
                breakpoint()
                if pixcvec is not None:
                    cycle_id = os.path.split(pixcvec)[1].split('_')[4]
                elif pixc is not None:
                    cycle_id = os.path.split(pixc)[1].split('_')[4]
                else:
                    print("no cycle provided")
                    return
            node_df = node_df[node_df['cycle_id']==int(cycle_id)]
            # filter on pass_id too
            pass_id = args.pass_id
            node_df = node_df[node_df['pass_id']==int(pass_id)]
            node_df = node_df.sort_values(['node_id'])
            # populate node data
            proc_df = populate_rivertile(proc_df, node_df, cycle_id, pass_id)
        # call the make plots routine
        make_plots(proc_tile, proc_df, truth_tile, pixcvec, pixc,
                   truth_pixcvec, truth_pixc, args.reach_id,
                   reach_error, nodes=args.nodes,
                   pixc_truth=pixc_truth, mt_wse=mt_wse, mt_width=mt_width,
                   sandbox=args.sandbox)
        plt.show()
    else:
        print('Input file', proc_tile, 'does not exist')


if __name__ == "__main__":
    main()
