#!/usr/bin/env python
'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s):

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
    pixc_rivertile, gdem_rivertile, metrics=None, truth=None, data=None, scene=None):
    '''
    load reaches from a particular scene/tile, compute metrics,
    and accumulate the data, truth and metrics (if input)
    '''
    truth_tmp, data_tmp = SWOTRiver.analysis.riverobs.load_rivertiles(
        gdem_rivertile, pixc_rivertile)
    if (len(truth_tmp.reaches.reach_id)<=0) or (
       len(data_tmp.reaches.reach_id)<=0):
        # do nothing if truth or data file have no reach data
        return metrics, truth, data, scene
    # handle masked arrays here
    truth_tmp, data_tmp = handle_bad_reaches(truth_tmp, data_tmp)
    #
    SWOTRiver.analysis.riverobs.match_reaches(truth_tmp, data_tmp)
    tmp = SWOTRiver.analysis.riverobs.get_metrics(
        truth_tmp.reaches, data_tmp.reaches)
    # get the scene
    path_parts = os.path.abspath(pixc_rivertile).split('/')
    scene0 = path_parts[-4] # assumes particular directory structure...
    # put in the pass and tile too
    cycle_pass_tile_flavor = path_parts[-3].split('_')
    scene1 = scene0+"_"+cycle_pass_tile_flavor[3]+"_"+cycle_pass_tile_flavor[4]
    scene_tmp = [scene1 for item in data_tmp.reaches.reach_id] 
    #for k in range(len(data_tmp.reaches.reach_id)):
    #    scene_tmp.append(scene0)
    # accumulate if needed
    if metrics is None:
        metrics = tmp
        truth = truth_tmp
        data = data_tmp
        scene = scene_tmp
    else:
        for name, value in tmp.items():
            metrics[name] = np.append(metrics[name], value)
        truth = truth.append(truth_tmp)
        data = data.append(data_tmp)
        scene = np.append(scene, scene_tmp)
    return metrics, truth, data, scene

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
                metrics, truth, data, scene = load_and_accumulate(
                    pixc_rivertile, gdem_rivertile, metrics, truth, data, scene)
                
    else:
        metrics, truth, data, scene = load_and_accumulate(
            args.pixc_rivertile, args.gdem_rivertile)

   
    #SWOTRiver.analysis.tabley.print_table(metrics)
    SWOTRiver.analysis.riverobs.print_metrics(metrics, truth, scene)
    print("\nFor All Data")
    SWOTRiver.analysis.riverobs.print_errors(metrics)
    print("\nFor 10km<xtrk_dist<60km and width>100m and area>(1km)^2")
    msk, fit_error, dark_frac = SWOTRiver.analysis.riverobs.mask_for_sci_req(
        metrics, truth, data, scene)
    SWOTRiver.analysis.riverobs.print_metrics(metrics, truth, scene, msk, fit_error, dark_frac)
    print("\nFor all Data with 10km<xtrk_dist<60km and width>100m and area>(1km)^2")
    SWOTRiver.analysis.riverobs.print_errors(metrics, msk)

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
