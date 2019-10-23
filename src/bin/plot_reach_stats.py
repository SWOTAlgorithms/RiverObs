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




def handle_bad_reaches(truth_tmp, data_tmp):
    """
    only keep reaches that all the pertinent data is good for both truth and meaured 
    """
    bad_reaches = np.array([])
    #main_keys = ['area_detct','height','slope','width']
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
    print("bad_data_reaches",bad_reaches)
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
        truth=None, data=None, scene=None, sig0=None, bad_scenes=[]):
    '''
    load reaches from a particular scene/tile, compute metrics,
    and accumulate the data, truth and metrics (if input)
    '''
    truth_tmp, data_tmp = SWOTRiver.analysis.riverobs.load_rivertiles(
        gdem_rivertile, pixc_rivertile)
    if (len(truth_tmp.reaches.reach_id)<=0) or (
       len(data_tmp.reaches.reach_id)<=0):
        # do nothing if truth or data file have no reach data
        return metrics, truth, data, scene, sig0
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
    # accumulate if needed
    if metrics is None:
        metrics = tmp
        truth = truth_tmp
        data = data_tmp
        scene = scene_tmp
        sig0 = sig0_avg
    else:
        # don't accumulate scenes designated as bad
        scene2 = scene1.split('_')[0]
        if scene2 in bad_scenes:
            return metrics, truth, data, scene, sig0
        for name, value in tmp.items():
            metrics[name] = np.append(metrics[name], value)
        truth = truth.append(truth_tmp)
        data = data.append(data_tmp)
        scene = np.append(scene, scene_tmp)
        sig0 = np.append(sig0, sig0_avg)
    return metrics, truth, data, scene, sig0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_rivertile', help='e.g., pixc/rivertile.nc')
    parser.add_argument('gdem_rivertile', help='e.g., gdem/rivertile.nc')
    parser.add_argument('-t', '--title')
    parser.add_argument('-p', '--print', action='store_true')
    parser.add_argument('-d', '--dir',
                        help='base directory to search for rivertile files')
    parser.add_argument('-f', '--flavor',
                        help='flavor of pixel cloud processing')
    args = parser.parse_args()

    metrics = None
    truth = None
    data = None
    scene = None
    sig0 = None
    bad_scenes = []#['3356',]# these scenes will be excluded from analysis
    print("args.dir: ",args.dir)
    if args.dir is not None:
        # TODO: rethink this glob call to be more general
        #   i.e., not depend so much on assumed directory structure
        pixc_rivertile_list = glob.glob(
            args.dir+"/*/*"+args.flavor+'/'+args.pixc_rivertile)
        for pixc_rivertile in pixc_rivertile_list:
            basedir = pixc_rivertile[0:-len(args.pixc_rivertile)]
            gdem_rivertile = basedir + args.gdem_rivertile
            print (pixc_rivertile)
            print (gdem_rivertile)
            if os.path.isfile(gdem_rivertile):
                # check if the gdem and the pixc are the same
                metrics, truth, data, scene, sig0 = load_and_accumulate(
                    pixc_rivertile, gdem_rivertile,
                    metrics, truth, data, scene, sig0, bad_scenes)
                
    else:
        metrics, truth, data, scene, sig0 = load_and_accumulate(
            args.pixc_rivertile, args.gdem_rivertile)

   
    #SWOTRiver.analysis.tabley.print_table(metrics)
    passfail = SWOTRiver.analysis.riverobs.get_passfail()
    msk, fit_error, dark_frac, reach_len = SWOTRiver.analysis.riverobs.mask_for_sci_req(
        metrics, truth, data, scene, sig0=sig0)
    SWOTRiver.analysis.riverobs.print_metrics(
        metrics, truth, scene, with_node_avg=True,
        passfail=passfail, dark_frac=dark_frac, reach_len=reach_len)
    print("\nFor All Data")
    SWOTRiver.analysis.riverobs.print_errors(metrics, with_node_avg=True)
    print("\nFor 10km<xtrk_dist<60km and width>100m and area>(1km)^2 and reach len>=10km")
    SWOTRiver.analysis.riverobs.print_metrics(
        metrics, truth, scene, msk, fit_error,
        dark_frac, with_node_avg=True, passfail=passfail, reach_len=reach_len)
    print("\nFor all Data with 10km<xtrk_dist<60km and width>100m and area>(1km)^2 and reach len>=10km")
    SWOTRiver.analysis.riverobs.print_errors(metrics, msk, with_node_avg=True)

    # plot slope error vs reach length
    plt.figure()
    plt.scatter(metrics['slope_t'][msk], metrics['slope'][msk], c=reach_len[msk]/1e3)
    plt.colorbar()
    plt.xlabel('slope (cm/km)')
    plt.ylabel('slope error (cm/km)')
    plt.title('reach length (km)')
    plt.grid()
    #
    plt.figure()
    plt.scatter(reach_len[msk]/1e3, metrics['slope_t'][msk], c=metrics['slope'][msk])
    plt.colorbar()
    plt.title('slope (cm/km)')
    plt.ylabel('slope error (cm/km)')
    plt.xlabel('reach length (km)')
    plt.grid()
    #
    plt.figure()
    plt.scatter(metrics['slope'][msk], metrics['wse'][msk], c=reach_len[msk]/1e3)
    plt.colorbar()
    plt.ylabel('height error (cm)')
    plt.xlabel('slope error (cm/km)')
    plt.title('reach length (km)')
    plt.grid()
    # plot the fit error histogram
    #plt.figure()
    #plt.scatter(fit_error[msk], metrics['height'][msk], c=dark_frac[msk])
    #plt.colorbar()
    #plt.xlabel('fit_error')
    #plt.ylabel('metrics')
    
    if args.print:
        filenames = ['reach-area.png', 'reach-wse.png',
                     'reach-slope.png','reach-wse-vs-area.png']
        if args.title is not None:
            filenames = [args.title + '-' + name for name in filenames]
    else:
        filenames = [None, None, None, None]
    SWOTRiver.analysis.riverobs.AreaPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[0], msk=msk)
    SWOTRiver.analysis.riverobs.HeightPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[1], msk=msk)
    SWOTRiver.analysis.riverobs.SlopePlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[2], msk=msk)
    SWOTRiver.analysis.riverobs.HeightVsAreaPlot(
        truth.reaches, data.reaches, metrics, args.title, filenames[3], msk=msk)
    if not args.print:
        print('show')
        plt.show()


if __name__ == '__main__':
    main()
