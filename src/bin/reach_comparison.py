#!/usr/bin/env python

# This goes through all reaches in all rivertiles, finds the reaches with the
# greatest error, and then plots then on a dashboard with height, area, and
# pixel assignment results.
#
# The reaches are plotted one-by-one, ranked from highest error to lowest
# error. The user may specify which type of error to sort by, e.g. wse, area,
# slope, or slope2 (enhanced slope). Written for performance assessment and
# troubleshooting of SWOT rivertiles generated using simulated data.

# Copyright (c) 2020-, California Institute of Technology ("Caltech"). U.S.
# Government sponsorship acknowledged.
# All rights reserved.

# Author(s): Cassie Stuurman


import argparse
import os
import glob
import math
import numpy as np
import SWOTRiver.analysis.riverobs
import SWOTWater
import plot_reach
import pdb
import matplotlib.pyplot as plt
import warnings
import pandas as pd
from pathlib import Path

from os import path
from plot_reach_stats import load_and_accumulate
from plot_riverobs import NodePlots, ParamPlots, HeightPlots
from collections import namedtuple


def get_input_files(basedir, slc_dir, pixc_dir, proc_rivertile,
                    truth_rivertile, pge, truth_basedir=None, truth_only=False,
                    scene_filter=None):
    # get all pixc files 'rivertile.nc' and find associated truth file
    missing_truth_count = 0
    if truth_only:
        print('doing truth-to-truth file aggregation...')
        proc_rivertile_list = glob.glob(os.path.join(basedir, '*', '*',
                                                     slc_dir,
                                                     proc_rivertile,
                                                     'river_data',
                                                     'rivertile.nc'))
    else:
        if pge:
            # Filenames have different format due to using PGE stand in
            proc_rivertile_list = glob.glob(
                os.path.join(basedir, '*', '*', slc_dir, pixc_dir,
                             proc_rivertile, 'river_data', '*RiverTile*.nc'))
        else:
            proc_rivertile_list = glob.glob(os.path.join(basedir, '*', '*',
                                                         slc_dir,
                                                         pixc_dir,
                                                         proc_rivertile,
                                                         'river_data',
                                                         'rivertile.nc'))
    if scene_filter:  # only keep tiles matching the input scene
        proc_rivertile_list = [k for k in proc_rivertile_list
                               if scene_filter in k]
    if len(proc_rivertile_list) == 0:
        raise Exception('No rivertiles found, check input directory names')
    truth_rivertile_list = []
    for index, rivertile in enumerate(proc_rivertile_list):
        truth_file = get_truth_file(
            proc_rivertile, pixc_dir, rivertile, truth_rivertile, basedir,
            truth_basedir, truth_only)
        if os.path.exists(truth_file):
            truth_rivertile_list.append(truth_file)
        else:
            warn_str = 'Truth rivertile file ' + truth_file + ' does not exist.'
            warnings.warn(warn_str)
            truth_rivertile_list.append(None)
            missing_truth_count += 1
    print('total missing truths = ', missing_truth_count, 'out of',
          len(proc_rivertile_list), 'files.')
    return proc_rivertile_list, truth_rivertile_list


def get_truth_file(proc_dir, pixc_dir, proc_rivertile, truth_rivertile,
                   basedir, truth_basedir=None, truth_only=False):
    proc_file_parts = proc_rivertile.split('/')
    if truth_only:
        n_char = len(proc_dir) + 1 + len(proc_file_parts[-1] + '/'
                                         + proc_file_parts[-2])
    else:
        n_char = len(proc_dir) + len(pixc_dir) + 1 \
                 + len(proc_file_parts[-1] + '/' + proc_file_parts[-2]) + 1
    path = proc_rivertile[0:-n_char]
    if truth_basedir:
        truth_path = proc_rivertile.replace(basedir, truth_basedir)[0:-n_char]
        truth_file = glob.glob(os.path.join(
            truth_path, truth_rivertile, '*', '*iver*ile*.nc'))
        if truth_file:
            truth_file = truth_file[0]
        else:
            truth_file = ''
            print('Missing truth:', truth_path + truth_rivertile +
                  '/river_data/rivertile.nc')
    else:
        truth_file = glob.glob(os.path.join(
            path, truth_rivertile, '*', '*iver*ile*.nc'))
        if truth_file:
            truth_file = truth_file[0]
        else:
            truth_file = ''
            print('Missing truth:', path + truth_rivertile +
                  '/river_data/rivertile.nc')
    return truth_file


def get_pixc_file(proc_dir, proc_rivertile):
    n_char = len(proc_dir) + 1 + len('river_data/rivertile.nc') + 1
    proc_path = proc_rivertile[0:-n_char]
    pixc_file = proc_path + '/pixc_data/pixel_cloud.nc'
    if path.exists(pixc_file):
        return pixc_file
    else:
        print('Missing pixel cloud file, continuing...')
        return None


def get_errors(rivertile_list, truth_list, test, truth_filter):
    # Use existing analysis tools to obtain the node and reach level error
    # metrics and mask for the scientific requirement bounds.
    bad_scene = []

    metrics = None
    truth = None
    data = None
    scene = None
    scene_nodes = None
    sig0 = None
    good_rivertile_list = []
    good_truth_list = []

    test_count = 0  # this only increases if we're in script 'test' mode
    print('Retrieving rivertile errors...')
    if type(rivertile_list) is not list:  # only one input file. Make it a list
        rivertile_list = [rivertile_list]
        truth_list = [truth_list]
    for index, filename in enumerate(rivertile_list):
        if test_count <= 10:
            # get the error of that scene
            if rivertile_list[index] and truth_list[index]:
                metrics, truth, data, scene, scene_nodes, \
                sig0, has_reach_data = load_and_accumulate(rivertile_list[index],
                                                           truth_list[index],
                                                           metrics,
                                                           truth,
                                                           data,
                                                           scene,
                                                           scene_nodes,
                                                           sig0,
                                                           bad_scene,
                                                           truth_filter)
                if has_reach_data:
                    good_rivertile_list.append(rivertile_list[index])
                    good_truth_list.append(truth_list[index])
                    if test:
                        test_count += 1

    if len(metrics['area_total']) > 0:
        passfail = SWOTRiver.analysis.riverobs.get_passfail()
        msk, bounds, dark_frac, reach_len, reach_width, qual_flag, \
        rch_count = SWOTRiver.analysis.riverobs.mask_for_sci_req(
                truth, data, scene, scene_nodes
        )
        preamble = "\nFor " + str(bounds['min_xtrk']) + " km<xtrk_dist<" \
                   + str(bounds['max_xtrk']) + " km and width>" \
                   + str(bounds['min_width']) + " m and area>" \
                   + str(bounds['min_area']) + " m^2 and reach len>=" \
                   + str(bounds['min_length']) + " m\n"
        print(preamble)
        if any(msk):
            table = SWOTRiver.analysis.riverobs.print_errors(metrics,
                                                             with_node_avg=True)
            metrics_table = SWOTRiver.analysis.riverobs.print_metrics(
                metrics, truth, scene,
                msk,
                with_node_avg=True,
                passfail=passfail,
                reach_len=reach_len,
                preamble=preamble)
        else:
            print('No reaches in tile', filename,
                  'meet the science requirements.')
    else:
        print('No metrics obtainable.')
        metrics_table = []

    return metrics_table, good_rivertile_list, good_truth_list


def plot_worst_reaches(metrics, first_reach_index, rivertile_files,
                       truth_files, sort_param, proc_dir, out_dir):
    # calls plot_reach for all reaches sorted by worst error
    sorted_metrics = sort_errors(metrics, sort_param,
                                 rivertile_files, truth_files)
    # slice metrics values according to input first percentile of plotting
    sorted_metrics = sorted_metrics.iloc[first_reach_index:]
    ifig = 1
    for index, reach in enumerate(sorted_metrics['reach']):
        this_reach = sorted_metrics.iloc[index]
        rivertile_file = this_reach['rivertile_file']
        truth_file = this_reach['truth_rivertile']
        pixc_file = get_pixc_file(proc_dir, rivertile_file)
        scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(
            rivertile_file)
        truth_pixc = get_pixc_file(proc_dir, truth_file)
        truth_pixcvec = get_pixc_file(proc_dir, truth_file)
        gdem_file = get_gdem_from_rivertile(rivertile_file)
        pixc_vec = get_pixcvec_from_rivertile(rivertile_file)
        try:
            figure, axes = plot_reach.make_plots(rivertile_file, truth_file,
                                                 pixc_vec, pixc_file,
                                                 truth_pixcvec, truth_pixc,
                                                 reach, gdem_file, this_reach,
                                                 scene)

            if out_dir is not None and os.path.isdir(out_dir):
                ifig += 1
                plt.savefig(out_dir + str(ifig) + '_' + str(scene) + '_' +
                            str(reach))
            else:
                plt.show()


        except TypeError:
            print('couldn\'t make plot for', rivertile_file)


def sort_errors(metrics, sort_param, rivertile_files, truth_files):
    # ranks reach-level errors from largest to smallest absolute value
    errors = []
    sort_dict = {'wse': 'wse e (cm)',
                 'area': 'area_tot e (%)',
                 'width': 'width e (m)',
                 'slope': 'slp e (cm/km)',
                 'slope2': 'slp2 e (cm/km)',
                 'xtrk': 'xtrk (km)',
                 'area_det': 'area_det e (%)',
                 'wse_node_e': 'wse node e (cm)'}

    # match rivertile files to metrics table using scene/pass/side
    matched_rivertiles = []
    matched_truths = []
    for tile_index, tile_id in enumerate(metrics['scene_pass_tile']):
        for file_index, filename in enumerate(rivertile_files):
            if tile_id[-8:] in filename:
                matched_rivertiles.append(filename)
                matched_truths.append(truth_files[file_index])
                break

    # convert to dataframe and combine with input files for sorting
    metrics_df = pd.DataFrame.from_dict(metrics)
    metrics_df['rivertile_file'] = matched_rivertiles
    metrics_df['truth_rivertile'] = matched_truths
    metrics_df['sort_param_abs'] = metrics_df[sort_dict[sort_param]].abs()
    metrics_df.sort_values('sort_param_abs', ascending=False, inplace=True)
    return metrics_df


def get_gdem_from_rivertile(rivertile_file):
    # gets the input gdem file from an output pixc rivertile.nc. Hard coded.
    file_parts = rivertile_file.split('/')

    # assume the lidar scene name is the part before the cycle_* part
    ind = -6
    for n, part in enumerate(file_parts):
        if file_parts[n].startswith('cycle_'):
            ind = n-1
    lidar_scene = file_parts[ind]
    gdem_file = '/u/swot-fn-r0/swot/sim_proc_inputs/gdem-dem-truth-v9-nowet/' +  lidar_scene + '_lidar.nc'
    if os.path.exists(gdem_file):
        return gdem_file
    else:
        return None


def get_pixcvec_from_rivertile(rivertile_file):
    # gets a pixcvec file from an input rivertile. Hard coded.
    # TODO: make more general
    pixcvec_file = rivertile_file[0:-12] + '/pixcvec.nc'
    if path.exists(pixcvec_file):
        return pixcvec_file
    else:
        # pixcvec not found, try PGE format instead
        pixcvec_file = rivertile_file[0:-65] + 'PIXCVecRiver' \
                       + rivertile_file[-56:]
    if path.exists(pixcvec_file):
        return pixcvec_file
    else:
        print('Could not find pixcvec file.')
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_rivertile', type=str, default=None,
                        help='processed rivertile file (or basename)')
    parser.add_argument('truth_rivertile', type=str, default=None,
                        help='truth rivertile file (or basename)')
    parser.add_argument('--basedir', help='base directory of processing')
    parser.add_argument('--truth_basedir', type=str, default=None,
                        help='base directory of truth processing. Only use if'
                             'different than nominal processing.')
    parser.add_argument('-sb', '--slc_basename', type=str, default=None,
                        help='slc directory basename')
    parser.add_argument('-pb', '--pixc_basename', type=str, default=None,
                        help='pixc directory basename')
    parser.add_argument('--test_boolean',
                        help='set to "True" if testing script', default=False,
                        required=False)
    parser.add_argument('--percentile', type=int, default=100, required=False,
                        help='%%ile along the distribution of errors where you'
                             ' want to begin the analysis, 0-100')
    parser.add_argument('--sort_by', type=str, default='wse',
                        help='Which error class to sort by: wse area or slope')
    parser.add_argument('-t', '--truth_only', type=bool, default=False,
                        help='Compare truth rivertiles to truth rivertile, '
                             'True or False')
    parser.add_argument('-sf', '--scene_filter', type=str,
                        help='keep scenes matching this name e.g 3607, tanana')
    parser.add_argument('-trf', '--truth_filter', type=str, default=None,
                        nargs='+', help='filter out reaches by truth class. '
                                        'Options are tribs, non_linear, '
                                        'edge_node, partial_truth, multi_chn,'
                                        'bad_reach, wrong_dir, linear.')
    parser.add_argument('-pge', '--pge', action='store_true', default=None,
                        help='Flag that signifies we are looking for '
                             'pge-generated files. These have different names '
                             'and require a flag to correctly be found.')
    parser.add_argument('--out_dir', help='output directory for reach plots',
                        type=str, default=None)

    args = parser.parse_args()

    # check validity of input sort parameter
    sort_strs = ['wse', 'area', 'slope', 'slope2', 'xtrk', 'area_det',
                 'wse_node_e', 'width']
    if not any(args.sort_by in sort_strs for sort_str in sort_strs):
        raise Exception('Input sort string must be wse, area, or slope.')

    print('base directory is', args.basedir)
    proc_list, truth_list = get_input_files(args.basedir,
                                            args.slc_basename,
                                            args.pixc_basename,
                                            args.proc_rivertile,
                                            args.truth_rivertile,
                                            args.pge,
                                            args.truth_basedir,
                                            args.truth_only,
                                            args.scene_filter)

    metrics, good_rivertile_list, good_truth_list = get_errors(
        proc_list, truth_list, args.test_boolean, args.truth_filter)

    # start plotting at the error percentile of interest
    n_reaches = len(metrics['reach'])
    first_reach_index = int(np.floor(n_reaches*(100-args.percentile)/100))
    plot_worst_reaches(metrics,
                       first_reach_index,
                       good_rivertile_list,
                       good_truth_list,
                       args.sort_by,
                       args.proc_rivertile,
                       args.out_dir)


if __name__ == "__main__":
    main()
