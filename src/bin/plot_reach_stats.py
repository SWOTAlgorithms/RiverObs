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

import SWOTRiver.analysis.riverobs
import SWOTRiver.analysis.tabley
import glob
from pathlib import Path
import pdb

def handle_bad_reaches(truth_tmp, data_tmp):
    """
    only keep reaches that all the pertinent data is good for both truth and measured
    """
    bad_reaches = np.array([])
    # bad_reaches = np.array([73150600041, 73150600551, 73150600581, 73160300011, 73216000261, 73218000071, 73218000251,
    #                73220700271, 73220900221, 73240100201, 73240200041, 74230900181, 74230900191, 74230900251,
    #                74262700251, 74266300011, 74269800121, 74291700071, 74291900011, 74292100271])  # from Rui
    main_keys = ['area_total','wse','slope','width']
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
        # setting all variables to nans for bad reaches
        # makes them be excluded when matching
        # reach ids later between the data and truth
        if isinstance(truth_tmp.reaches[key], np.ma.MaskedArray):
            tmp = truth_tmp.reaches[key].data.astype(np.double)
            tmp[msk_t]=np.nan
            truth_tmp.reaches[key] = tmp
        if isinstance(data_tmp.reaches[key], np.ma.MaskedArray):
            tmp = data_tmp.reaches[key].data.astype(np.double)
            tmp[msk_d]=np.nan
            data_tmp.reaches[key] = tmp
    return truth_tmp, data_tmp

def load_and_accumulate(
        pixc_rivertile, gdem_rivertile, metrics=None,
        truth=None, data=None, scene=None, scene_nodes=None, sig0=None, bad_scenes=[]):
    '''
    load reaches from a particular scene/tile, compute metrics,
    and accumulate the data, truth and metrics (if input)
    '''
    truth_tmp, data_tmp = SWOTRiver.analysis.riverobs.load_rivertiles(
        gdem_rivertile, pixc_rivertile)
    if (len(truth_tmp.reaches.reach_id)<=0) or (
       len(data_tmp.reaches.reach_id)<=0):
        # do nothing if truth or data file have no reach data
        print('File', gdem_rivertile, 'has no reach data')
        return metrics, truth, data, scene, scene_nodes, sig0
    # handle masked arrays here
    truth_tmp, data_tmp = handle_bad_reaches(truth_tmp, data_tmp)
    #
    SWOTRiver.analysis.riverobs.match_reaches(truth_tmp, data_tmp)
    wse_node_avg, sig0_avg = SWOTRiver.analysis.riverobs.compute_average_node_error(data_tmp, truth_tmp)
    tmp = SWOTRiver.analysis.riverobs.get_metrics(
        truth_tmp.reaches, data_tmp.reaches, wse_node_avg=wse_node_avg)
    # get the scene
    scene1 = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(pixc_rivertile)
    scene_tmp = [scene1 for item in data_tmp.reaches.reach_id]
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
    return metrics, truth, data, scene, scene_nodes, sig0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proc_rivertile', type=str, default=None,
                        help='processed rivertile file (or basename)')
    parser.add_argument('truth_rivertile', type=str, default=None,
                        help='truth rivertile file (or basename)')
    parser.add_argument('--basedir', type=str, default=None,
                        help='base directory of processing')
    parser.add_argument('-sb', '--slc_basename', type=str, default=None,
                        help='slc directory basename')
    parser.add_argument('-pb', '--pixc_basename', type=str, default=None,
                        help='pixc directory basename')
    parser.add_argument('-eb','--pixc_errors_basename', type=str, default=None,
                        help = "pixc systematic errors basename")
    parser.add_argument('-o', '--outdir', type=str,help='output directory for print files')
    parser.add_argument('-tf','--table_fname', type=str, default=None, help='output filename for summary table')
    parser.add_argument('-t', '--title')
    parser.add_argument('-p', '--print', action='store_true')
    args = parser.parse_args()

    metrics = None
    truth = None
    data = None
    scene = None
    scene_nodes = None
    sig0 = None

    bad_scenes = [] # e.g. ['3356',...] these scenes will be excluded from analysis
    # bad_scenes = ['1782', '1793', '1830', '1847',
    #               '2367', '2447', '2800', '2801', '2819', '298', '3007', '3024', '316',
    #               '3356', '3376', '3379', '3382', '3396',
    #               '3413', '3421', '3424', '3442', '3455', '3487', '3574', '3575', '3579', '3586', '3596', '3607',
    #               '3769', '3770', '3781', '3783', '3787', '3815',  '3821',
    #               '730', 'platte_0036_0012', 'platte_0120_0290',
    #               'sacramento_0221_0264', 'sacramento_0314_0264', 'sacramento_0220_0249',
    #               'tanana_0434_0002', 'tanana_0156_0001']
    # 2543, 3454, 3583, 3824, 3836,
    # subset [2543 3454 3583 3824 3836]
    # bad_scenes = [316, 730, 1782, 3376]
    print("args.basedir: ", args.basedir)
    if args.basedir is not None:
        if args.slc_basename is None or args.pixc_basename is None:
            print('Must specify at least slc_basename and pixc_basename '
                  + 'if aggregating stats')
            return

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

        # If proc_rivertile input is a basename, get the actual rivertile
        proc_rivertile_list = [os.path.join(proc_rivertile, 'river_data', 'rivertile.nc')
                               if os.path.isdir(proc_rivertile) else proc_rivertile
                               for proc_rivertile in proc_rivertile_list]

        if args.pixc_errors_basename is not None:
            truth_rivertile_list = \
                [os.path.join(*Path(proc_rivertile).parts[:-5], args.truth_rivertile)
                 for proc_rivertile in proc_rivertile_list]
        else:
            truth_rivertile_list = \
                [os.path.join(*Path(proc_rivertile).parts[:-4], args.truth_rivertile)
                 for proc_rivertile in proc_rivertile_list]

        # If truth_rivertile input is a basename, get the actual rivertile
        truth_rivertile_list = [os.path.join(truth_rivertile, 'river_data', 'rivertile.nc')
                                if os.path.isdir(truth_rivertile) else truth_rivertile
                                for truth_rivertile in truth_rivertile_list]

        for proc_rivertile, truth_rivertile in zip(proc_rivertile_list, truth_rivertile_list):
            if os.path.isfile(proc_rivertile) and os.path.isfile(truth_rivertile):
                metrics, truth, data, scene, scene_nodes, sig0 = load_and_accumulate(
                    proc_rivertile, truth_rivertile,
                    metrics, truth, data, scene, scene_nodes, sig0, bad_scenes)

    else:
        # Inputs can be either rivertile files, or basenames
        proc_rivertile = args.proc_rivertile
        truth_rivertile = args.truth_rivertile
        if os.path.isdir(proc_rivertile):
            proc_rivertile = os.path.join(proc_rivertile, 'river_data', 'rivertile.nc')
        if os.path.isdir(truth_rivertile):
            truth_rivertile = os.path.join(truth_rivertile, 'river_data', 'rivertile.nc')

        metrics, truth, data, scene, scene_nodes, sig0 = load_and_accumulate(
            proc_rivertile, truth_rivertile)

    # get pass/fail and mask for science requirements
    passfail = SWOTRiver.analysis.riverobs.get_passfail()
    msk, fit_error, bounds, dark_frac, reach_len = SWOTRiver.analysis.riverobs.mask_for_sci_req(
        metrics, truth, data, scene, scene_nodes, sig0=sig0)

    # generate output filenames for metric tables and images
    if args.print:
        filenames = ['reach-area.png', 'reach-wse.png',
                     'reach-slope.png','reach-wse-vs-area.png']
        if args.table_fname is None:
            args.table_fname = 'table'
        table_fname = args.table_fname
        if args.title is not None:
            filenames = [args.title + '-' + name for name in filenames]
            table_fname = args.title + '_' + table_fname
        if args.outdir is not None:
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
        fname=table_fname + '_all.txt')
    print("\nFor All Data")
    SWOTRiver.analysis.riverobs.print_errors(metrics,
                                             fname=table_fname + '_errors_all.txt',
                                             preamble = 'For All Data',
                                             with_node_avg=True)
    # printing masked data to table
    preamble = "\nFor " + str(bounds['min_xtrk']) + " km<xtrk_dist<" + str(bounds['max_xtrk']) + " km and width>" \
               + str(bounds['min_width']) + " m and area>" + str(bounds['min_area']) + " m^2 and reach len>=" \
               + str(bounds['min_length']) + " m"
    print(preamble)
    SWOTRiver.analysis.riverobs.print_metrics(
        metrics, truth, scene, msk, fit_error,
        dark_frac,
        with_node_avg=True,
        passfail=passfail,
        reach_len=reach_len,
        preamble=preamble,
        fname=table_fname + '.txt')
    SWOTRiver.analysis.riverobs.print_errors(metrics, msk,
                                             fname=table_fname + '_errors.txt',
                                             preamble=preamble,
                                             with_node_avg=True)

    # plot slope error vs reach length
    fig, (ax1, ax2, ax3) = plt.subplots(1,3)
    ax1.scatter(metrics['slope_t'][msk], metrics['slope'][msk], c=reach_len[msk]/1e3)
    c = reach_len[msk]/1e3
    ax1.set(xlabel='slope (cm/km)', ylabel='slope error (cm/km)')
    ax1.set_title('reach length (km)')
    #ax1.colorbar()
    #fig.colorbar(reach_len[msk]/1e3, ax=ax1)
    ax1.grid()

    ax2.scatter(reach_len[msk]/1e3, metrics['slope_t'][msk], c=metrics['slope'][msk])
    ax2.set(xlabel='reach length (km)', ylabel='slope error (cm/km)')
    c = metrics['slope'][msk]
    ax2.set_title('slope (cm/km)')
    #fig.colorbar(metrics['slope'][msk], ax=ax2)
    #ax2.colorbar()
    ax2.grid()

    ax3.scatter(metrics['slope'][msk], metrics['wse'][msk], c=reach_len[msk]/1e3)
    ax3.set(xlabel='slope error (cm/km)', ylabel='height error (cm)')
    ax3.set_title('reach length (km)')
    c = reach_len[msk]/1e3
    #fig.colorbar(c, ax=ax3)
    ax3.grid()

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
    SWOTRiver.analysis.riverobs.AreaPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[0], msk=msk)
    SWOTRiver.analysis.riverobs.HeightPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[1], msk=msk)
    SWOTRiver.analysis.riverobs.SlopePlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[2], msk=msk)
    SWOTRiver.analysis.riverobs.HeightVsAreaPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[3], msk=msk)
    if not args.print:
        print('showing stats...')
        f = open('temp.txt', 'r')
        file_contents = f.read()
        print(file_contents)
        f = open('temp_errors.txt', 'r')
        file_contents = f.read()
        print(file_contents)
        plt.show()
        os.remove('temp.txt')
        os.remove('temp_all.txt')
        os.remove('temp_errors.txt')
        os.remove('temp_errors_all.txt')



if __name__ == '__main__':
    main()
