#!/usr/bin/env python
'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s): Brent Williams

'''
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

import SWOTRiver.analysis.riverobs
import SWOTRiver.analysis.tabley
import glob
from pathlib import Path
import warnings
import pdb

def handle_bad_reaches(truth_tmp, data_tmp, truth_filter):
    """
    only keep reaches where pertinent data is good for both truth and measured
    """
    bad_reaches = np.array([])
    if truth_filter:
        truth_classes = SWOTRiver.analysis.riverobs.get_truth_classes()
        for key in truth_filter:
            bad_reaches = np.concatenate([bad_reaches, truth_classes.get(key)])
    main_keys = ['area_total','wse','slope','width','slope2']
    for key in main_keys:
        # if any of these are masked, throw out the entire
        # reach by setting all elements to nan
        if (isinstance(truth_tmp.reaches[key], np.ma.MaskedArray)):
            bad = truth_tmp.reaches.reach_id[truth_tmp.reaches[key].mask]
            if len(bad)>0:
                bad_reaches = np.concatenate((bad_reaches, bad.data))
        if (isinstance(data_tmp.reaches[key], np.ma.MaskedArray)):
            bad = data_tmp.reaches.reach_id[data_tmp.reaches[key].mask]
            if len(bad)>0:
                bad_reaches = np.concatenate((bad_reaches, bad.data))
    #print("bad_data_reaches",bad_reaches)
    msk_t = [False for reach in truth_tmp.reaches.reach_id]
    msk_d = [False for reach in data_tmp.reaches.reach_id]
    for i,reach in enumerate(truth_tmp.reaches.reach_id):
        if reach in bad_reaches:
            msk_t[i] = True
    for i,reach in enumerate(data_tmp.reaches.reach_id):
        if reach in bad_reaches:
            msk_d[i] = True
    for key in truth_tmp.reaches.variables:
        # setting all variables to nans for bad reaches makes them be excluded
        # when matching reach ids later between the data and truth.
        # Must check if truth variables also in data for backwards SWORD
        # compatibility.
        try:
            if key in data_tmp.reaches.variables:
                if isinstance(truth_tmp.reaches[key], np.ma.MaskedArray):
                    tmp = truth_tmp.reaches[key].data.astype(np.double)
                    tmp[msk_t]=np.nan
                    truth_tmp.reaches[key] = tmp
                if isinstance(data_tmp.reaches[key], np.ma.MaskedArray):
                    tmp = data_tmp.reaches[key].data.astype(np.double)
                    tmp[msk_d]=np.nan
                    data_tmp.reaches[key] = tmp
        except ValueError:  # occurs if data is not a number
            continue
    return truth_tmp, data_tmp

def load_and_accumulate(
        pixc_rivertile, gdem_rivertile, metrics=None,
        truth=None, data=None, scene=None, scene_nodes=None, sig0=None,
        bad_scenes=None, yukon_good_tiles=None, truth_filter=None):
    '''
    load reaches from a particular scene/tile, compute metrics,
    and accumulate the data, truth and metrics (if input)
    '''

    if bad_scenes is None:
        bad_scenes = []
    if yukon_good_tiles is None:
        yukon_good_tiles = []
    truth_tmp, data_tmp = SWOTRiver.analysis.riverobs.load_rivertiles(
        gdem_rivertile, pixc_rivertile)
    # do nothing if truth or data file have no reach data
    if len(truth_tmp.reaches.reach_id)<=0:
        print('File', gdem_rivertile, 'has no reach data.')
        has_reach_data = False
        return metrics, truth, data, scene, scene_nodes, sig0, has_reach_data
    if len(data_tmp.reaches.reach_id) <= 0:
        print('File', pixc_rivertile, 'has no reach data.')
        has_reach_data = False
        return metrics, truth, data, scene, scene_nodes, sig0, has_reach_data
    # handle masked arrays here
    truth_tmp, data_tmp = handle_bad_reaches(truth_tmp, data_tmp, truth_filter)
    #
    SWOTRiver.analysis.riverobs.match_reaches(truth_tmp, data_tmp)
    wse_node_avg, sig0_avg = SWOTRiver.analysis.riverobs.\
        compute_average_node_error(data_tmp, truth_tmp)
    tmp = SWOTRiver.analysis.riverobs.get_metrics(
        truth_tmp.reaches, data_tmp.reaches, wse_node_avg=wse_node_avg)
    # get the scene
    scene1 = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_rivertile)
    scene_tmp = [scene1 for item in data_tmp.reaches.reach_id]
    scene_splits = scene1.split('_')
    if scene_splits[0] == 'yukonflats' and scene1 not in yukon_good_tiles:
        print('File', pixc_rivertile, 'not in Yukon good list.')
        has_reach_data = False
        return metrics, truth, data, scene, scene_nodes, sig0, has_reach_data
    #print("data_tmp:",data_tmp.nodes.reach_id)
    scene_tmp2 = [scene1 for item in data_tmp.nodes.reach_id]
    # accumulate if needed
    if metrics is None:
        metrics = tmp
        truth = truth_tmp
        data = data_tmp
        scene = scene_tmp
        scene_nodes = scene_tmp2
        sig0 = sig0_avg
    else:
        # don't accumulate scenes designated as bad
        scene_parts = scene1.split('_')
        scene2 = scene_parts[0]
        if len(scene_parts)>3:
            scene2 = scene_parts[0]+'_'+scene_parts[1]+'_'+scene_parts[2]
        if scene2 in bad_scenes:
            return metrics, truth, data, scene, scene_nodes, sig0
        for name, value in tmp.items():
            metrics[name] = np.append(metrics[name], value)
        truth = truth.append(truth_tmp)
        data = data.append(data_tmp)
        scene = np.append(scene, scene_tmp)
        scene_nodes = np.append(scene_nodes, scene_tmp2)
        sig0 = np.append(sig0, sig0_avg)
    has_reach_data = True
    return metrics, truth, data, scene, scene_nodes, sig0, has_reach_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_rivertile', type=str, default=None,
                        help='processed rivertile file (or basename)')
    parser.add_argument('truth_rivertile', type=str, default=None,
                        help='truth rivertile file (or basename)')
    parser.add_argument('--basedir', type=str, default=None,
                        help='base directory of processing')
    parser.add_argument('--truth_basedir', type=str, default=None,
                        help='base directory of truth files (only use if '
                             'different than processing basedir)')
    parser.add_argument('-sb', '--slc_basename', type=str, default=None,
                        help='slc directory basename')
    parser.add_argument('-pb', '--pixc_basename', type=str, default=None,
                        help='pixc directory basename')
    parser.add_argument('-eb','--pixc_errors_basename', type=str, default=None,
                        help = "pixc systematic errors basename")
    parser.add_argument('-sf', '--scene_filter',
                        type=str, help='keep scenes matching this name '
                                       '(e.g. 3607, tanana)')
    parser.add_argument('-trf', '--truth_filter', type=str, default=None,
                        nargs='+',
                        help='filter out reaches by truth class. Options are '
                             'tribs, non_linear, edge_node, partial_truth,'
                             'multi_chn, bad_reach, wrong_dir, linear.')
    parser.add_argument('-o', '--outdir', type=str,
                        help='output directory for print files')
    parser.add_argument('-t', '--title', type=str, default=None,
                        help='output filename for figure/table files')
    parser.add_argument('-p', '--print_me', action='store_true', default=None,
                        help='Add -p if you want to save the output files.'
                             'If you flag -p, it is recommended that you also'
                             ' specify --outdir and --title')
    parser.add_argument('-pge', '--pge', action='store_true', default=None,
                        help='Flag that signifies we are looking for '
                             'pge-generated files. These have different names '
                             'and require a flag to correctly be found.')
    parser.add_argument('-tm', '--test_mode', action='store_true', default=None,
                        help='Flag that triggers testing mode. Only a small'
                             'number of files will be read in.')
    args = parser.parse_args()

    metrics = None
    truth = None
    data = None
    scene = None
    scene_nodes = None
    sig0 = None

    bad_scenes = []
    yukon_good_tiles = ['yukonflats_0358_035R',
                        'yukonflats_0358_035L',
                        'yukonflats_0515_274R']
    truth_filter = args.truth_filter

    print("args.basedir: ", args.basedir)
    if args.basedir:
        if args.slc_basename is None or args.pixc_basename is None:
            print('Must specify at least slc_basename and pixc_basename '
                  + 'if aggregating stats')
            return

        if args.outdir or args.title:
            if args.print_me is None:
                args.print_me = True

        # TODO: Right now it's hardcoded that the truth data lives under the slc
        # base directory, and the proc data lives under the pixc base directory
        if args.pixc_errors_basename is not None:
            proc_rivertile_list = glob.glob(os.path.join(
                args.basedir, '*', '*', args.slc_basename, args.pixc_basename,
                args.pixc_errors_basename, args.proc_rivertile))
        else:
            proc_rivertile_list = glob.glob(os.path.join(
                args.basedir, '*', '*', args.slc_basename, args.pixc_basename,
                args.proc_rivertile))
        if args.scene_filter:  # only keep tiles matching the input scene
            proc_rivertile_list = [k for k in proc_rivertile_list
                                   if args.scene_filter in k]

        # If proc_rivertile input is a basename, get the actual rivertile
        if args.pge:
            proc_rivertile_list = [glob.glob(
                proc_rivertile + '/river_data/' + '*RiverTile*.nc'
            ) if os.path.isdir(proc_rivertile) else proc_rivertile
                                   for proc_rivertile in proc_rivertile_list]
            proc_rivertile_list = [
                ''.join(tile) for tile in proc_rivertile_list
            ]
        else:
            proc_rivertile_list = [os.path.join(
                proc_rivertile, 'river_data', 'rivertile.nc'
            ) if os.path.isdir(proc_rivertile) else proc_rivertile
                                   for proc_rivertile in proc_rivertile_list]
        if args.truth_basedir is not None:
            len_basedir = len(Path(args.basedir).parts)
            if args.pixc_errors_basename is not None:
                truth_rivertile_list = [
                    os.path.join(
                        args.truth_basedir,
                        *Path(proc_rivertile).parts[len_basedir:-5],
                        args.truth_rivertile
                    ) for proc_rivertile in proc_rivertile_list
                ]
            else:
                truth_rivertile_list = [
                    os.path.join(
                        args.truth_basedir,
                        *Path(proc_rivertile).parts[len_basedir:-4],
                        args.truth_rivertile
                    ) for proc_rivertile in proc_rivertile_list
                ]
        else:  # truth basedir same as nominal basedir
            if args.pixc_errors_basename is not None:
                truth_rivertile_list = [
                    os.path.join(
                        *Path(proc_rivertile).parts[:-5], args.truth_rivertile
                    ) for proc_rivertile in proc_rivertile_list
                ]
            else:
                truth_rivertile_list = [
                    os.path.join(
                        *Path(proc_rivertile).parts[:-4], args.truth_rivertile
                    ) for proc_rivertile in proc_rivertile_list
                ]

        # If truth_rivertile input is a basename, get the actual rivertile
        truth_rivertile_list = [
            os.path.join(
                truth_rivertile, 'river_data', 'rivertile.nc'
            ) if os.path.isdir(truth_rivertile) else truth_rivertile
            for truth_rivertile in truth_rivertile_list
        ]
        # error checking
        truth_sum = 0
        nominal_sum = 0
        for truth_rivertile in truth_rivertile_list:
            if os.path.isfile(truth_rivertile):
                truth_sum+=1
        if truth_sum==0:
            print('No truth files found. Check input variable names.')
        for nominal_rivertile in proc_rivertile_list:
            if os.path.isfile(nominal_rivertile):
                nominal_sum+=1
        if nominal_sum==0:
            print('No nominal tiles found. Check input variable names.')

        if args.test_mode:
            proc_rivertile_list = proc_rivertile_list[0:10]

        for proc_rivertile, truth_rivertile in zip(proc_rivertile_list,
                                                   truth_rivertile_list):
            if os.path.isfile(proc_rivertile) and os.path.isfile(
                    truth_rivertile):
                metrics, truth, data, scene, scene_nodes, \
                sig0, has_reach_data = load_and_accumulate(proc_rivertile,
                                                           truth_rivertile,
                                                           metrics, truth,
                                                           data, scene,
                                                           scene_nodes,
                                                           sig0, bad_scenes,
                                                           yukon_good_tiles,
                                                           truth_filter)

    else:
        # Inputs can be either rivertile files, or basenames
        proc_rivertile = args.proc_rivertile
        truth_rivertile = args.truth_rivertile
        if os.path.isdir(proc_rivertile):
            proc_rivertile = os.path.join(
                proc_rivertile, 'river_data', 'rivertile.nc'
            )
        if os.path.isdir(truth_rivertile):
            truth_rivertile = os.path.join(
                truth_rivertile, 'river_data', 'rivertile.nc'
            )

        metrics, truth, data, scene, \
        scene_nodes, sig0, has_reach_data = load_and_accumulate(
            proc_rivertile, truth_rivertile)

    # get pass/fail and mask for science requirements
    passfail = SWOTRiver.analysis.riverobs.get_passfail()
    msk, bounds, dark_frac, reach_len, reach_width, qual_flag, \
    rch_count = SWOTRiver.analysis.riverobs.mask_for_sci_req(
            truth, data, scene, scene_nodes)

    # generate output filenames for metric tables and images
    if args.print_me:
        filenames = ['reach-area.png', 'reach-wse.png',
                     'reach-slope.png','reach-wse-vs-area.png']
        if args.title:
            filenames = [args.title + '-' + name for name in filenames]
            table_fname = args.title
        else:
            table_fname = 'table'
        if args.outdir:
            filenames = [args.outdir + '/' + name for name in filenames]
            table_fname = args.outdir + '/' + table_fname
            # if the directory doesnt exist create it
            if not os.path.exists(args.outdir):
                os.makedirs(args.outdir)
    else:
        filenames = [None, None, None, None]
        table_fname = 'temp'

    # printing all data to table
    SWOTRiver.analysis.riverobs.print_metrics(
        metrics, truth, scene,
        with_node_avg=True,
        passfail=passfail,
        preamble='For All Data',
        dark_frac=dark_frac,
        reach_len=reach_len,
        qual_flag=qual_flag,
        fname=table_fname + '_all.txt')
    print("\nFor All Data")
    SWOTRiver.analysis.riverobs.print_errors(metrics,
                                             fname=table_fname + '_errors_all.txt',
                                             preamble = 'For All Data',
                                             with_node_avg=True)
    # printing masked data to table
    preamble = "Total reach num: " + str(rch_count['rch_num_total']) \
               + "\nnum of reaches filtered due to xtrk: " + str(
                   rch_count['bad_xtrk']) \
               + "\nnum of reaches filtered due to width: " + str(
                   rch_count['bad_width']) \
               + "\nnum of reaches filtered due to area: " + str(
                   rch_count['bad_area']) \
               + "\nnum of reaches filtered due to length: " + str(
                   rch_count['bad_length']) \
               + "\nnum of reaches filtered due to obs frac: " + str(
                   rch_count['bad_obs_frac']) \
               + "\nnum of reaches filtered due to dark frac: " + str(
                   rch_count['bad_dark_frac']) \
               + "\nnum of reaches filtered due to qual flag: " + str(
                   rch_count['bad_qual']) \
               + "\nnum of reaches filtered due to xtrk ratio: " + str(
                   rch_count['bad_xtrk_ratio']) \
               + "\n\nFor " + str(bounds['min_xtrk']) + " km<xtrk_dist<" \
               + str(bounds['max_xtrk']) + " km and width>" \
               + str(bounds['min_width']) + " m and area>" \
               + str(bounds['min_area']) + " m^2 \n and reach len>=" \
               + str(bounds['min_length']) + " m and obs frac >=" \
               + str(bounds['min_obs_frac']) + " and truth ratio >= "\
               + str(bounds['min_truth_obs_frac']) + " and xtrk proportion >= "\
               + str(bounds['min_xtrk_ratio']) + " and qual flag <= "\
               + str(bounds['max_qual_flag'])

    print(preamble)
    SWOTRiver.analysis.riverobs.print_metrics(
        metrics, truth, scene, msk, dark_frac,
        with_node_avg=True,
        passfail=passfail,
        reach_len=reach_len,
        qual_flag=qual_flag,
        preamble=preamble,
        fname=table_fname + '.txt')
    SWOTRiver.analysis.riverobs.print_errors(metrics, msk,
                                             fname=table_fname + '_errors.txt',
                                             preamble=preamble,
                                             with_node_avg=True)

    # plot slope error vs reach length
    fig, (ax1, ax2, ax3) = plt.subplots(1,3)
    ax1.scatter(metrics['slope_t'][msk], metrics['slope'][msk],
                c=reach_len[msk]/1e3)
    ax1.set(xlabel='slope (cm/km)', ylabel='slope error (cm/km)')
    ax1.set_title('reach length (km)')
    ax1.grid()

    ax2.scatter(reach_len[msk]/1e3, metrics['slope_t'][msk],
                c=metrics['slope'][msk])
    ax2.set(xlabel='reach length (km)', ylabel='slope error (cm/km)')
    ax2.set_title('slope (cm/km)')
    #fig.colorbar(metrics['slope'][msk], ax=ax2)
    #ax2.colorbar()
    ax2.grid()

    ax3.scatter(metrics['slope'][msk], metrics['wse'][msk],
                c=reach_len[msk]/1e3)
    ax3.set(xlabel='slope error (cm/km)', ylabel='height error (cm)')
    ax3.set_title('reach length (km)')
    c = reach_len[msk]/1e3
    #fig.colorbar(c, ax=ax3)
    ax3.grid()
    plt.show()

    # plot slope error vs reach length
    # plt.figure()
    # plt.scatter(metrics['slope_t'][msk], metrics['slope'][msk], c=reach_len[msk]/1e3)
    # plt.colorbar()
    # plt.xlabel('slope (cm/km)')
    # plt.ylabel('slope error (cm/km)')
    # plt.title('reach length (km)')
    # plt.grid()
    # #
    # plt.figure()
    # plt.scatter(reach_len[msk]/1e3, metrics['slope_t'][msk], c=metrics['slope'][msk])
    # plt.colorbar()
    # plt.title('slope (cm/km)')
    # plt.ylabel('slope error (cm/km)')
    # plt.xlabel('reach length (km)')
    # plt.grid()
    # #
    # plt.figure()
    # plt.scatter(metrics['slope'][msk], metrics['wse'][msk], c=reach_len[msk]/1e3)
    # plt.colorbar()
    # plt.ylabel('height error (cm)')
    # plt.xlabel('slope error (cm/km)')
    # plt.title('reach length (km)')
    # plt.grid()
    # plot the fit error histogram
    #plt.figure()
    #plt.scatter(fit_error[msk], metrics['height'][msk], c=dark_frac[msk])
    #plt.colorbar()
    #plt.xlabel('fit_error')
    #plt.ylabel('metrics')
    plt.figure()
    ecdf = sm.distributions.ECDF(reach_len/1e3)
    x = np.linspace(min(reach_len)/1e3, max(reach_len)/1e3)
    y = ecdf(x)
    plt.step(x, y, label='before')
    ecdf2 = sm.distributions.ECDF(reach_len[msk]/1e3)
    x2 = np.linspace(min(reach_len[msk])/1e3, max(reach_len[msk])/1e3)
    y2 = ecdf2(x2)
    plt.step(x2, y2, label='masked')
    plt.ylabel('cdf')
    plt.xlabel('reach length (km)')
    plt.title('reach length cdfs')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    ecdf = sm.distributions.ECDF(reach_width)
    x = np.linspace(min(reach_width), max(reach_width))
    y = ecdf(x)
    plt.step(x, y, label='before')
    ecdf2 = sm.distributions.ECDF(reach_width[msk])
    x2 = np.linspace(min(reach_width[msk]), max(reach_width[msk]))
    y2 = ecdf2(x2)
    plt.step(x2, y2, label='masked')
    plt.ylabel('cdf')
    plt.xlabel('reach width (m)')
    plt.title('reach width cdfs')
    plt.legend()
    plt.grid()
    plt.show()

    plt.figure()
    xtrk_dist = abs(truth['reaches']['xtrk_dist']/1e3)
    ecdf = sm.distributions.ECDF(xtrk_dist)
    x = np.linspace(min(xtrk_dist), max(xtrk_dist))
    y = ecdf(x)
    plt.step(x, y, label='before')
    ecdf2 = sm.distributions.ECDF(xtrk_dist[msk])
    x2 = np.linspace(min(xtrk_dist[msk]), max(xtrk_dist[msk]))
    y2 = ecdf2(x2)
    plt.step(x2, y2, label='masked')
    plt.ylabel('cdf')
    plt.xlabel('xtrk dist (km)')
    plt.title('xtrk dist cdfs')
    plt.legend()
    plt.grid()
    plt.show()

    SWOTRiver.analysis.riverobs.AreaPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[0], msk=msk)
    SWOTRiver.analysis.riverobs.HeightPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[1], msk=msk)
    SWOTRiver.analysis.riverobs.SlopePlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[2], msk=msk)
    SWOTRiver.analysis.riverobs.HeightVsAreaPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[3], msk=msk)
    if args.print_me:
        file = table_fname
        print('showing stats and saving to file...')
    else:
        file = 'temp'
        print('showing stats...')
    plt.show()
    f = open(file + '.txt', 'r')
    file_contents = f.read()
    print(file_contents)
    f = open(file + '_errors.txt', 'r')
    file_contents = f.read()
    print(file_contents)
    if args.print_me is None:
        os.remove('temp.txt')
        os.remove('temp_all.txt')
        os.remove('temp_errors.txt')
        os.remove('temp_errors_all.txt')


if __name__ == '__main__':
    main()
