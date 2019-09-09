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

FIGSIZE = (8, 4) #(12, 8)
DPI = 200
CMAP = 'plasma'


class ReachPlot():
    def __init__(self, truth, data, metrics, title=None, filename=None, msk=True):
        self.truth = truth
        self.data = data
        self.metrics = metrics
        self.msk = msk
        self.title = title
        self.filename = filename
        self.figure, self.axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        self.plot()

    def plot(self, independent, dependent, color_key, color_abs=True, outlierlim=None):
        if color_abs:
            col = np.abs(self.data.reaches[color_key][self.msk])
        else:
            col = self.data.reaches[color_key][self.msk]
        # clip outliers and plot them with different marker
        outliers = np.ones(np.shape(independent[self.msk]),dtype=bool)
        indep = independent[self.msk].copy()
        dep = dependent[self.msk].copy()
        if outlierlim is not None:
            if len(outlierlim)==4:
                outliers = np.logical_or(
                    np.logical_or(
                        np.logical_or(dep<outlierlim[0],
                                      dep>outlierlim[1]),
                                      indep<outlierlim[2]),
                                      indep>outlierlim[3])
                dep[dep<outlierlim[0]] = outlierlim[0]
                dep[dep>outlierlim[1]] = outlierlim[1]
                indep[indep<outlierlim[2]] = outlierlim[2]
                indep[indep>outlierlim[3]] = outlierlim[3]
            else:
                outliers = np.logical_or(
                    (dep<outlierlim[0]),(dep>outlierlim[1]))
                #indep = independent.copy()
                dep[dep<outlierlim[0]] = outlierlim[0]
                dep[dep>outlierlim[1]] = outlierlim[1]
        not_outliers = np.logical_not(outliers)
        scatter = self.axis.scatter(
            indep[not_outliers], dep[not_outliers], c=col[not_outliers],
            cmap=CMAP)
        scatter_out = self.axis.scatter(
            indep[outliers], dep[outliers], c=col[outliers],
            cmap=CMAP, marker = '+')
        colorbar = plt.colorbar(scatter)
        colorbar.set_label(color_key)
        self.plot_requirements()
        self.finalize()

    def finalize(self):
        self.axis.grid()
        if self.title is not None:
            self.axis.set_title(self.title)
        if self.filename is not None:
            self.figure.savefig(self.filename)

class HeightVsAreaPlot(ReachPlot):
    def plot(self):
        super().plot(self.metrics['area'], self.metrics['height'], 'xtrk_dist',
                     outlierlim=(-50,50,-50,50))

    def plot_requirements(self):
        self.axis.axvline(x=15, color='g')
        self.axis.axvline(x=-15, color='g')
        self.axis.axhline(y=10, color='r')
        self.axis.axhline(y=-10, color='r')

    def finalize(self):
        self.axis.set_xlabel('area % error')
        self.axis.set_ylabel('wse error (cm)')
        super().finalize()

class AreaPlot(ReachPlot):
    def plot(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.reaches.area_total)
        super().plot(true_area, self.metrics['area'], 'xtrk_dist', outlierlim=(-50,50))

    def plot_requirements(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.reaches.area_total[self.msk])
        buff = 100
        i = 1
        self.axis.plot([np.sqrt(50*10e3), 1000], [i*25, i*25], '--g')
        self.axis.plot([1000, np.sqrt(170*10e3)], [i*15, i*15], '--y')
        self.axis.plot(
            [np.sqrt(170*10e3), np.amax(true_area)+buff], [i*15, i*15], '--r')
        self.axis.set_xlim((0, np.amax(true_area)+buff))
        #self.axis.set_ylim((-50,50))
        self.axis.legend(
            ['SG for $A>0.7 km^2$', 'BSM for $A>1 km^2$',
             'TSM for $A>1.3 km^2$', 'data','outlier clipped'],
            loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15))
        i = -1
        self.axis.plot([np.sqrt(50*10e3), 1000], [i*25, i*25], '--g')
        self.axis.plot([1000, np.sqrt(170*10e3)], [i*15, i*15], '--y')
        self.axis.plot(
            [np.sqrt(170*10e3), np.amax(true_area)+buff], [i*15, i*15], '--r')

    def finalize(self):
        self.axis.set_ylabel('area % error')
        self.axis.set_xlabel('sqrt reach area (m)')
        super().finalize()


class HeightPlot(ReachPlot):
    def plot(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.reaches.area_total)
        super().plot(true_area, self.metrics['height'], 'xtrk_dist', outlierlim=(-50,50))

    def plot_requirements(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.reaches.area_total[self.msk])
        buff = 100
        i = 1
        self.axis.plot([250, 1000], [i*25, i*25], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*10, i*10], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*11, i*11], '--r')
        self.axis.set_xlim((0,np.amax(true_area)+buff))
        #self.axis.set_ylim((-50,50))
        self.axis.legend(
            ['BSM for $A>(250m)^2$', 'BSM for $A>1 km^2$',
             'TSM for $A>1 km^2$', 'data','outlier clipped'],
            loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15))
        i = -1
        self.axis.plot([250, 1000], [i*25, i*25], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*10, i*10], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*11, i*11], '--r')

    def finalize(self):
        self.axis.set_ylabel('wse error (cm)')
        self.axis.set_xlabel('sqrt reach area (m)')
        super().finalize()


class SlopePlot(ReachPlot):
    def plot(self):
        true_width = self.truth.reaches.width
        super().plot(true_width, self.metrics['slope'], 'xtrk_dist', outlierlim=(-5,5))

    def plot_requirements(self):
        true_width = self.truth.reaches.width
        buff = 100
        i = 1
        self.axis.plot([100, np.amax(true_width)+buff], [i*1.7, i*1.7], '--y')
        self.axis.plot([100, np.amax(true_width)+buff], [i*3, i*3], '--r')
        self.axis.set_xlim((0, np.amax(true_width)+buff))
        self.axis.legend(
            ['BSM for $A>1 km^2$', 'TSM for $A>1 km^2$',
             'data','outlier clipped'],
            loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.1))
        i = -1
        self.axis.plot([100, np.amax(true_width)+buff], [i*1.7, i*1.7], '--y')
        self.axis.plot([100, np.amax(true_width)+buff], [i*3, i*3], '--r')

    def finalize(self):
        self.axis.set_ylabel('slope error (cm/km)')
        self.axis.set_xlabel('river width (m)')
        super().finalize()

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
    plt.figure()
    plt.scatter(fit_error[msk], metrics['height'][msk], c=dark_frac[msk])
    plt.colorbar()
    plt.xlabel('fit_error')
    plt.ylabel('metrics')
    
    if args.print:
        filenames = ['reach-area.png', 'reach-wse.png',
                     'reach-slope.png','reach-wse-vs-area.png']
        if args.title is not None:
            filenames = [args.title + '-' + name for name in filenames]
    else:
        filenames = [None, None, None, None]
    AreaPlot(truth, data, metrics, args.title, filenames[0], msk=msk)
    HeightPlot(truth, data, metrics, args.title, filenames[1], msk=msk)
    SlopePlot(truth, data, metrics, args.title, filenames[2], msk=msk)
    HeightVsAreaPlot(truth, data, metrics, args.title, filenames[3], msk=msk)
    if not args.print:
        print('show')
        plt.show()


if __name__ == '__main__':
    main()
