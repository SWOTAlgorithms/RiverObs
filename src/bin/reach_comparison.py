#!/usr/bin/env python

# This goes through all reaches in all lidar scenes and finds the best and worst scene/pass/side scenes and reaches
# It then plots all reaches, reach by reach, ranked from highest slope error to lowest slope error

# To do: clean this up to be an actual functions dump rather than a railroad, and create a separate railroad script

# Copyright (c) 2020-, California Institute of Technology ("Caltech"). U.S.
# Government sponsorship acknowledged.
# All rights reserved.

# Author(s): Cassie Stuurman


import argparse
import os
import glob
import math
from plot_reach_stats import load_and_accumulate
from plot_riverobs import NodePlots, ParamPlots, HeightPlots
import SWOTRiver.analysis.riverobs
import SWOTWater
import plot_reach
import pdb
import matplotlib.pyplot as plt

def get_input_files(basedir, slc_dir, pixc_dir, proc_rivertile, truth_dir):
    # get all pixc files 'rivertile.nc' and find associated truth file
    #file_str = basedir + '**/cycle_**/' + slc_dir + '/' + proc_rivertile + '/' + 'river_data/rivertile.nc'
    proc_rivertile_list = glob.glob(os.path.join(basedir, '*', '*', slc_dir, pixc_dir, proc_rivertile, 'river_data', 'rivertile.nc'))
    truth_rivertile_list = []
    n_char = len(pixc_dir) + 1 + len(proc_rivertile) + len('river_data/rivertile.nc') + 1
    for file in proc_rivertile_list:
        truth_file =  get_truth_file(proc_rivertile, pixc_dir, file, truth_dir)
        truth_rivertile_list.append(truth_file)
    return proc_rivertile_list, truth_rivertile_list

def get_truth_file(proc_dir, pixc_dir, proc_rivertile, truth_dir):
    n_char = len(proc_dir) + len(pixc_dir) +1 + len('river_data/rivertile.nc') + 1
    path = proc_rivertile[0:-n_char]
    #truth_file = path + truth_dir + '/river_data/rivertile.nc'
    truth_file = proc_rivertile[0:43] + proc_rivertile[44:120] + truth_dir + '/river_data/rivertile.nc' # hardcoded for new set of rivertiles to avoid truth re-run
    #truth_file = path + truth_dir + '/river_data/rivertile.nc'
    return truth_file

def get_errors(pixc_list, gdem_list, test, verbose=True):
    worst_area = 0
    worst_detected_area = 0
    worst_wse = 0
    worst_slope = 0
    most_bad_reaches = 0
    worst_reach = 0 # need to populate!
    
    best_area = 10000
    best_detected_area = 10000
    best_wse = 10000
    best_slope = 10000

    scene_error_list = []
    reach_error_list = []

    metrics = None
    truth = None
    data = None
    scene = None
    sig0 = None
    table = None
    bad_scenes = []         #['3356',]# these scenes will be excluded from analysis
    good_pixc_list = []

    test_count = 0                      # this only increases if we're in script testing mode
    if type(pixc_list) is not list:
        # function was called for a single file
        print('Retrieving errors for single file...')
        try:
            metrics, truth, data, scene, sig0 = load_and_accumulate(pixc_list, gdem_list)
        except FileNotFoundError:
            print('Pixc rivertile has no matching gdem rivertile')
        passfail = SWOTRiver.analysis.riverobs.get_passfail()
        if truth:
            msk, fit_error, dark_frac, reach_len = SWOTRiver.analysis.riverobs.mask_for_sci_req(
                metrics, truth, data, scene, sig0=sig0)
            print("\nFor 10km<xtrk_dist<60km and width>100m and area>(1km)^2 and reach len>=10km")
            metrics_table = SWOTRiver.analysis.riverobs.print_metrics(
                metrics, truth, scene, msk, fit_error,
                dark_frac, with_node_avg=True, passfail=passfail, reach_len=reach_len)
            table = SWOTRiver.analysis.riverobs.print_errors(metrics, msk, with_node_avg=True)
            reach_error_list.append(metrics_table)
            return reach_error_list

    else: # function was called for a list of files
        print('Retrieving errors for all rivertiles...')
        for index, filename in enumerate(pixc_list):
            if test_count<=5:
                # get the error of that scene
                try:
                    metrics, truth, data, scene, scene_nodes, sig0 = load_and_accumulate(pixc_list[index], gdem_list[index])
                except FileNotFoundError:
                    print('Pixc rivertile has no matching gdem rivertile')

                passfail = SWOTRiver.analysis.riverobs.get_passfail()

                if truth:
                    msk, fit_error, dark_frac, reach_len = SWOTRiver.analysis.riverobs.mask_for_sci_req(
                        metrics, truth, data, scene, sig0=sig0)
                    if not any(msk):
                        print('No reaches in scene are within sci req bounds')
                        table=None
                    else:
                        print("\nFor 10km<xtrk_dist<60km and width>100m and area>(1km)^2 and reach len>=10km")
                        table = SWOTRiver.analysis.riverobs.print_errors(metrics, msk, with_node_avg=True)
                        metrics_table = SWOTRiver.analysis.riverobs.print_metrics(
                                            metrics, truth, scene, msk, fit_error,
                                            dark_frac, with_node_avg=True, passfail=passfail, reach_len=reach_len)
                        scene_error_list.append(table)
                        reach_error_list.append(metrics_table)
                        good_pixc_list.append(pixc_list[index])
                else:
                    print('Load and accumulation failure')
            if table:
                # compare error max/min to collection (olympics)
                if abs(table['mean'][0]) > worst_area:
                    worst_area = table['mean'][0]
                    worst_area_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][1]) > worst_detected_area:
                    worst_detected_area = table['mean'][1]
                    worst_d_area_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][2]) > worst_wse:
                    worst_wse = table['mean'][2]
                    worst_wse_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][3]) > worst_slope:
                    worst_slope = table['mean'][3]
                    worst_slope_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][0]) < best_area:
                    best_area = table['mean'][0]
                    best_area_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][1]) < best_detected_area:
                    best_detected_area = table['mean'][1]
                    best_d_area_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][2]) < best_wse:
                    best_wse = table['mean'][2]
                    best_wse_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if abs(table['mean'][3]) < best_slope:
                    best_slope = table['mean'][3]
                    best_slope_scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_list[index])
                if test:
                    test_count += 1
        # end file loop

        # if verbose:
        #     print('Worst scenes are', worst_area_scene, worst_d_area_scene, worst_wse_scene, worst_slope_scene)
        #     print('\nWorst Mean Area Error ', worst_area,
        #           '\nWorst mean detected area error', worst_detected_area,
        #           '\nWorst mean wse error', worst_wse,
        #           '\nWorst mean slope error   ', worst_slope)
        #     print('\nBest scenes are', best_area_scene, best_d_area_scene, best_wse_scene, best_slope_scene)
        #     print('\nBest Mean Area Error ', best_area,
        #           '\nBest mean detected area error', best_detected_area,
        #           '\nBest mean wse error', best_wse,
        #           '\nBest mean slope error', best_slope)

        return scene_error_list, reach_error_list, good_pixc_list

def plot_worst_reaches(reach_errors, rivertile_files, proc_dir, pixc_dir, truth_basedir):
    # calls plot_reach for all reaches, from worst to best slope (or area, hardcoded) error
    reach_wse_errors = sort_wse_errors(reach_errors, rivertile_files)
    for index, reach in enumerate(reach_wse_errors):
        rivertile_file = reach[0]
        reach_id = reach[2]
        errors = [reach[3], reach[4], reach[5], reach[6], reach[7]] # slope e, wse_e, area total e, area detected e
        scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(rivertile_file)
        truth_file = get_truth_file(proc_dir, pixc_dir, rivertile_file, truth_basedir)
        gdem_file = get_gdem_from_pixc(rivertile_file)
        figure, axes = plot_reach.make_plots(rivertile_file, truth_file, reach_id, gdem_file, errors, scene)

        plt.show()

def sort_slope_errors(reach_error_list, pixc_list):
    # ranks reach-level slope errors from largest to smallest absolute value
    slope_errors = []
    for scene_index, scene in enumerate(reach_error_list):
        for reach_index, reach in enumerate(reach_error_list[scene_index]['reach']):
            slope_error = reach_error_list[scene_index]['slp e (cm/km)'][reach_index]
            wse_error = reach_error_list[scene_index]['wse e (cm)'][reach_index]
            wse_range = reach_error_list[scene_index]['wse range (cm)'][reach_index]
            area_error = reach_error_list[scene_index]['area_tot e (%)'][reach_index]
            area_dtct_error = reach_error_list[scene_index]['area_det e (%)'][reach_index]
            width_error = reach_error_list[scene_index]['width e (m)'][reach_index]
            scene = reach_error_list[scene_index]['scene_pass_tile'][reach_index]
            slope_errors.append([pixc_list[scene_index], scene, reach,
                                 slope_error, wse_error, wse_range, area_error, area_dtct_error, width_error])
    slope_errors.sort(key=takeSlope, reverse=True)
    return slope_errors

def sort_wse_errors(reach_error_list, pixc_list):
    # ranks reach-level slope errors from largest to smallest absolute value
    wse_errors = []
    for scene_index, scene in enumerate(reach_error_list):
        for reach_index, reach in enumerate(reach_error_list[scene_index]['reach']):
            slope_error = reach_error_list[scene_index]['slp e (cm/km)'][reach_index]
            wse_error = reach_error_list[scene_index]['wse e (cm)'][reach_index]
            wse_r_u = reach_error_list[scene_index]['wse r u (cm)'][reach_index]
            wse_t_r_u = reach_error_list[scene_index]['wse t r u (cm)'][reach_index]
            norm_wse_e = wse_error/math.sqrt(wse_r_u**2 + wse_t_r_u**2)
            area_error = reach_error_list[scene_index]['area_tot e (%)'][reach_index]
            area_dtct_error = reach_error_list[scene_index]['area_det e (%)'][reach_index]
            width_error = reach_error_list[scene_index]['width e (m)'][reach_index]
            scene = reach_error_list[scene_index]['scene_pass_tile'][reach_index]
            wse_errors.append([pixc_list[scene_index], scene, reach,
                                 slope_error, wse_error, norm_wse_e, area_error, area_dtct_error, width_error])
    wse_errors.sort(key=takeWse, reverse=True)
    return wse_errors

def sort_area_errors(reach_error_list, pixc_list):
    # ranks reach-level area errors from largest to smallest absolute value
    area_errors = []
    for scene_index, scene in enumerate(reach_error_list):
        for reach_index, reach in enumerate(reach_error_list[scene_index]['reach']):
            slope_error = reach_error_list[scene_index]['slp e (cm/km)'][reach_index]
            wse_error = reach_error_list[scene_index]['wse e (cm)'][reach_index]
            wse_r_u = reach_error_list[scene_index]['wse r u (cm)'][reach_index]
            wse_t_r_u = reach_error_list[scene_index]['wse t r u (cm)'][reach_index]
            norm_wse_e = wse_error/math.sqrt(wse_r_u**2 + wse_t_r_u**2)
            area_error = reach_error_list[scene_index]['area_tot e (%)'][reach_index]
            area_dtct_error = reach_error_list[scene_index]['area_det e (%)'][reach_index]
            width_error = reach_error_list[scene_index]['width e (m)'][reach_index]
            scene = reach_error_list[scene_index]['scene_pass_tile'][reach_index]
            area_errors.append([pixc_list[scene_index], scene, reach,
                                 slope_error, wse_error, norm_wse_e, area_error, area_dtct_error, width_error])
    area_errors.sort(key=takeArea, reverse=True)
    return area_errors

def sort_scene_errors(scene_error_list, pixc_list, error_var):
    # error var is slope, wse, area detected, or area total error
    # returns scenes sorted by input error parameter
    errors = []
    index_dict = {'slope_error': 4,
                  'wse_error': 2,
                  'area_dtct_error': 1,
                  'area_error': 0}
    for scene_index, scene in enumerate(scene_error_list):
        error = scene_error_list[scene_index]['mean'][index_dict[error_var]]
        errors.append([pixc_list[scene_index], error])
    sorted_errors = errors.sort(key=takeSceneError, reverse=True)
    return sorted_errors

def print_best_worst_scenes(param_errors, param_str):
    print('Worst ' + param_str + ' scenes', param_errors[0:10])
    print('\nBest ' + param_str + ' scenes', param_errors[-10:])

def takeSlope(elem):
    print(elem[3])
    return abs(elem[3])

def takeArea(elem):
    return abs(elem[6])

def takeWse(elem):
    print('wse norm error is', elem[5])
    return abs(elem[5])

def takeSceneError(elem):
    return abs(elem[1])

def get_gdem_from_pixc(pixc_file):
    # gets the input gdem file from an output pixc rivertile.nc. Hard coded.
    file_parts = pixc_file.split('/')
    lidar_scene = file_parts[6]
    gdem_file = '/u/swot-fn-r0/swot/sim_proc_inputs/gdem-dem-truth-v9-nowet/' +  lidar_scene + '_lidar.nc'
    return gdem_file
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_rivertile', type=str, default=None,
                        help='processed rivertile file (or basename)')
    parser.add_argument('truth_rivertile', type=str, default=None,
                        help='truth rivertile file (or basename)')
    parser.add_argument('basedir', help='base directory of processing')
    parser.add_argument('slc_basename', type=str, default=None,
                        help='slc directory basename')
    parser.add_argument('pixc_basename', type=str, default=None,
                        help='pixc directory basename')
    parser.add_argument('--test_boolean', help='set to "True" if testing script', default=False, required=False)
    args = parser.parse_args()
    error_table = []
    
    print('base directory is', args.basedir)
    pixc_list, gdem_list = get_input_files(args.basedir, args.slc_basename, args.pixc_basename, args.proc_rivertile, args.truth_rivertile)
    scene_errors, reach_errors, pixc_list = get_errors(pixc_list, gdem_list, args.test_boolean)
    scene_slope_error = sort_scene_errors(scene_errors, pixc_list, 'slope_error')
    scene_wse_error = sort_scene_errors(scene_errors, pixc_list, 'wse_error')
    scene_area_error = sort_scene_errors(scene_errors, pixc_list, 'area_error')
    scene_area_dtct_error = sort_scene_errors(scene_errors, pixc_list, 'area_dtct_error')
    # print_best_worst_scenes(scene_slope_error, 'Slope')
    # print_best_worst_scenes(scene_wse_error, 'wse')
    # print_best_worst_scenes(scene_area_error, 'area total')
    # print_best_worst_scenes(scene_area_dtct_error, 'area detected')

    plot_worst_reaches(reach_errors, pixc_list, args.proc_rivertile, args.pixc_basename, args.truth_rivertile)

if __name__ == "__main__":
    main()
