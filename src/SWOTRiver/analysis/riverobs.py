"""
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s): Dustin Lagoy

"""
import warnings

import numpy as np
import os.path

import SWOTWater.products.product
import SWOTRiver.analysis.tabley

import matplotlib.pyplot as plt

FIGSIZE = (8, 4) #(12, 8)
DPI = 200
CMAP = 'plasma'


class ReachPlot():
    def __init__(self, truth, data, metrics, title=None, filename=None, msk=True, is_lake=False):
        self.truth = truth
        self.data = data
        self.metrics = metrics
        self.msk = msk
        self.title = title
        self.filename = filename
        self.is_lake = is_lake
        self.figure, self.axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        self.plot()

    def plot(self, independent, dependent, color_key, color_abs=True, outlierlim=None):
        if color_abs:
            col = np.abs(self.data[color_key][self.msk])
        else:
            col = self.data[color_key][self.msk]
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
        super().plot(self.metrics['area_total'], self.metrics['wse'], 'xtrk_dist',
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
        true_area = np.sqrt(self.truth.area_total)
        super().plot(true_area, self.metrics['area_total'], 'xtrk_dist', outlierlim=(-50,50))

    def plot_requirements(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.area_total[self.msk])
        buff = 100
        # set up the default legend and y-axis for req. for rivers
        legnd = ['Goal for $A>0.7 km^2$', 'Req. for $A>1 km^2$',
             'TSM for $A>1.3 km^2$', 'data','outlier clipped']
        goal_x = [np.sqrt(50*10e3), 1000]
        req_x = [1000, np.sqrt(170*10e3)]
        tsm_x = [np.sqrt(170*10e3), np.amax(true_area)+buff]
        if self.is_lake:
            # set up the default legend and y-axis for req. for lakes
            legnd = ['Goal for $(100 m)^2<A<(250m)^2$', 'Req. for $A>(250 m)^2$',
                'TSM for $A>1 km^2$', 'data','outlier clipped']
            goal_x = [100, 250]
            req_x = [250, np.amax(true_area)+buff]
            tsm_x = [1000, np.amax(true_area)+buff]
        
        
        i = 1
        self.axis.plot(goal_x, [i*25, i*25], '--g')
        self.axis.plot(req_x, [i*15, i*15], '--y')
        self.axis.plot(tsm_x, [i*15, i*15], '--r')
        self.axis.set_xlim((0, np.amax(true_area)+buff))
        #self.axis.set_ylim((-50,50))
        self.axis.legend(legnd,loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15))
        i = -1
        self.axis.plot(goal_x, [i*25, i*25], '--g')
        self.axis.plot(req_x, [i*15, i*15], '--y')
        self.axis.plot(tsm_x, [i*15, i*15], '--r')

    def finalize(self):
        self.axis.set_ylabel('area % error')
        if self.is_lake:
            self.axis.set_xlabel('sqrt lake area (m)')
        else:
            self.axis.set_xlabel('sqrt reach area (m)')
        super().finalize()


class HeightPlot(ReachPlot):
    def plot(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.area_total)
        super().plot(true_area, self.metrics['wse'], 'xtrk_dist', outlierlim=(-50,50))

    def plot_requirements(self):
        #true_area = np.sqrt(self.truth.reaches.area_detct)
        true_area = np.sqrt(self.truth.area_total[self.msk])
        buff = 100
        legnd = ['BSM for $A>(250m)^2$', 'BSM for $A>1 km^2$',
             'TSM for $A>1 km^2$', 'data','outlier clipped']
        goal_x =[250, 1000]
        req_x = [1000, np.amax(true_area)+buff]
        tsm_x = [1000, np.amax(true_area)+buff]
        # same numbers for lake and rivers...
        i = 1
        self.axis.plot(goal_x, [i*25, i*25], '--y')
        self.axis.plot(req_x, [i*10, i*10], '--y')
        self.axis.plot(tsm_x, [i*11, i*11], '--r')
        self.axis.set_xlim((0,np.amax(true_area)+buff))
        #self.axis.set_ylim((-50,50))
        self.axis.legend(legnd,loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15))
        i = -1
        self.axis.plot(goal_x, [i*25, i*25], '--y')
        self.axis.plot(req_x, [i*10, i*10], '--y')
        self.axis.plot(tsm_x, [i*11, i*11], '--r')

    def finalize(self):
        self.axis.set_ylabel('wse error (cm)')
        if self.is_lake:
            self.axis.set_xlabel('sqrt lake area (m)')
        else:
            self.axis.set_xlabel('sqrt reach area (m)')
        super().finalize()


class SlopePlot(ReachPlot):
    def plot(self):
        true_width = self.truth.width
        super().plot(true_width, self.metrics['slope'], 'xtrk_dist', outlierlim=(-5,5))

    def plot_requirements(self):
        true_width = self.truth.width
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

def load_rivertiles(truth_file, data_file):
    truth = SWOTWater.products.product.MutableProduct.from_ncfile(truth_file)
    data = SWOTWater.products.product.MutableProduct.from_ncfile(data_file)
    return truth, data


def match_rivertiles(truth, data):
    match_reaches(truth, data)
    match_nodes(truth, data)


def match_reaches(truth, data):
    # get common reach ids from intersection
    common_ids = np.intersect1d(data.reaches.reach_id, truth.reaches.reach_id)
    data_mapping = [
        np.where(data.reaches.reach_id == i)[0][0] for i in common_ids]
    true_mapping = [
        np.where(truth.reaches.reach_id == i)[0][0] for i in common_ids]

    new_reaches = data.reaches.copy(with_variables=False)
    for name in data.reaches.variables:
        new_reaches[name] = data.reaches[name][data_mapping]
    data.reaches = new_reaches

    new_reaches = truth.reaches.copy(with_variables=False)
    for name in truth.reaches.variables:
        new_reaches[name] = truth.reaches[name][true_mapping]
    truth.reaches = new_reaches


def match_nodes(truth, data):
    common_ids = np.intersect1d(data.nodes.node_id, truth.nodes.node_id)
    data_mapping = [
        np.where(data.nodes.node_id == node)[0][0] for node in common_ids]
    true_mapping = [
        np.where(truth.nodes.node_id == node)[0][0] for node in common_ids]

    new_nodes = data.nodes.copy(with_variables=False)
    for name in data.nodes.variables:
        new_nodes[name] = data.nodes[name][data_mapping]
    data.nodes = new_nodes

    new_nodes = truth.nodes.copy(with_variables=False)
    for name in truth.nodes.variables:
        new_nodes[name] = truth.nodes[name][true_mapping]
    truth.nodes = new_nodes


def get_metrics(truth, data, with_slope=True, with_width=True):
    metrics = {
        'area_total': (
            (data.area_total - truth.area_total) / truth.area_total) * 100.0,
        'area_detct':(
            (data.area_detct - truth.area_detct) / truth.area_detct) * 100.0,
        'wse': (data.wse - truth.wse) * 1e2,#convert m to cm
    }
    if with_slope:
        metrics['slope'] = (data.slope - truth.slope) * 1e5#convert from m/m to cm/km
    if with_width:
        metrics['width'] = data.width - truth.width
    return metrics

def compute_reach_fit_error(truth):
    fit_error = []
    for reach in truth.reaches['reach_id']:
        inds = (truth.nodes['reach_id'] == reach)
        ind = (truth.reaches['reach_id'] == reach)
        try:
            y0 = truth.nodes['wse'][inds]
            x0 = truth.nodes['node_id'][inds]
            # exclude nans and infs
            y = y0[np.isfinite(y0)]
            x = x0[np.isfinite(y0)]
            z = np.polyfit( x, y, 1)
            p = np.poly1d(z)
            err = np.nanmean(np.sqrt((y - p(x))**2))*100 #in cm
        except np.linalg.LinAlgError:
            err = 1000000000
            print("linAlgError caught. truth.nodes[wse]:",truth.nodes['wse'][inds])
        fit_error.append(err)
    return np.array(fit_error)
        
def mask_for_sci_req(metrics, truth, data, scene):
    # find reaches where the height profile linear fit is not that good
    # so we can filter out bogus/non-realistic reaches from the analysis
    fit_error = compute_reach_fit_error(truth)
    
    # now make the mask
    msk = np.logical_and((np.abs(truth.reaches['xtrk_dist'])>10000),
          np.logical_and((np.abs(truth.reaches['xtrk_dist'])<60000), 
          np.logical_and((truth.reaches['width']>100),
          np.logical_and((truth.reaches['area_total']>1e6),
          np.logical_and(fit_error < 150.0, data.reaches['dark_frac'] < 0.35)))))
    return msk, fit_error, data.reaches['dark_frac']
#
def get_scene_from_fnamedir(fnamedir):
    path_parts = os.path.abspath(fnamedir).split('/')
    scene0 = path_parts[-4] # assumes particular directory structure...
    # put in the pass and tile too
    cycle_pass_tile_flavor = path_parts[-3].split('_')
    if len(cycle_pass_tile_flavor)<5:
        scene='unknown'
    else:
        scene = scene0+"_"+cycle_pass_tile_flavor[3]+"_"+cycle_pass_tile_flavor[4]
        #scene = [scene1 for item in data_tmp.reaches.reach_id]
    return scene
#
def print_errors(metrics, msk=True, with_slope=True):
    # get statistics of area error
    area_68 = np.nanpercentile(abs(metrics['area_total'][msk]), 68)
    area_50 = np.nanpercentile(metrics['area_total'][msk], 50)
    area_mean = np.nanmean(metrics['area_total'][msk])
    area_num = np.count_nonzero(~np.isnan(metrics['area_total'][msk]))
    # get statistics of area_detct error
    area_d_68 = np.nanpercentile(abs(metrics['area_detct'][msk]), 68)
    area_d_50 = np.nanpercentile(metrics['area_detct'][msk], 50)
    area_d_mean = np.nanmean(metrics['area_detct'][msk])
    area_d_num = np.count_nonzero(~np.isnan(metrics['area_detct'][msk]))
    # get statistics of height error
    height_68 = np.nanpercentile(abs(metrics['wse'][msk]), 68)
    height_50 = np.nanpercentile(metrics['wse'][msk], 50)
    height_mean = np.nanmean(metrics['wse'][msk])
    height_num = np.count_nonzero(~np.isnan(metrics['wse'][msk]))
    table = {
        'metric': ['area_total (%)','area_detct (%)', 'wse (cm)'],
        '|68%ile|': [area_68, area_d_68, height_68],
        '50%ile': [area_50, area_d_50, height_50],
        'mean': [area_mean, area_d_mean, height_mean],
        'count': [area_num, area_d_num, height_num],
    }
    if with_slope:
        # get statistics of slope error
        slope_68 = np.nanpercentile(abs(metrics['slope'][msk]), 68)
        slope_50 = np.nanpercentile(metrics['slope'][msk], 50)
        slope_mean = np.nanmean(metrics['slope'][msk])
        slope_num = np.count_nonzero(~np.isnan(metrics['slope'][msk]))

        table['metric'].append('slope (cm/km)')
        table['|68%ile|'].append(slope_68)
        table['50%ile'].append(slope_50)
        table['mean'].append(slope_mean)
        table['count'].append(slope_num)
    
    SWOTRiver.analysis.tabley.print_table(table, precision=8)


def print_metrics(
        metrics, truth, scene=None, msk=None, fit_error=None,
        dark_frac=None, with_slope=True, with_width=True):
    table = {}
    if msk is None:
        msk = np.ones(np.shape(metrics['wse']),dtype = bool)
    table['wse e (cm)'] = metrics['wse'][msk]
    if with_slope:
        print(metrics)
        print(metrics['slope'])
        print(metrics['slope'][msk])
        print(msk)
        table['slp e (cm/km)'] = metrics['slope'][msk]
    table['area_tot e (%)'] = metrics['area_total'][msk]
    table['area_det e (%)'] = metrics['area_detct'][msk]
    if with_width:
        table['width e (m)'] = metrics['width'][msk]
        table['width (m)'] = truth.reaches['width'][msk]
    else:
        table['sqrt(area) (m)'] = np.sqrt(truth['area_total'][msk])
    try:
        table['reach'] = truth.reaches['reach_id'][msk]
        table['xtrk (km)'] = truth.reaches['xtrk_dist'][msk]/1e3
    except AttributeError as e:
        table['lake_id'] = truth['lakeobs_id'][msk]
        table['xtrk (km)'] = truth['xtrk_dist'][msk]/1e3
    if scene is not None:
        table['scene_pass_tile'] = np.array(scene)[msk]
    if fit_error is not None:
        table['fit_error (cm)'] = np.array(fit_error)[msk]
    if dark_frac is not None:
        table['dark_frac'] = np.array(dark_frac)[msk]
    
    SWOTRiver.analysis.tabley.print_table(table, precision=5)
