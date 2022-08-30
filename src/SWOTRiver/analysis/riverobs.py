"""
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s): Dustin Lagoy, Brent Williams

"""
import warnings
import pdb

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
        self.axis.legend(legnd, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.15))
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
        self.axis.legend(legnd, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.15)) # ,
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
             ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.1))
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

def compute_average_node_error(data, truth):
    err_out = np.zeros(np.shape(data.reaches['reach_id']))+np.nan
    sig0_out = np.zeros(np.shape(data.reaches['reach_id']))+np.nan
    for reach in data.reaches['reach_id']:
        inds = (data.nodes['reach_id'] == reach)
        inds_t = (truth.nodes['reach_id'] == reach)
        ind = (data.reaches['reach_id'] == reach)
        diff_wse = []
        weight = []
        sig0 = []
        nodes = data.nodes['node_id'][inds]
        nodes_t = truth.nodes['node_id'][inds_t]
        wse = data.nodes['wse'][inds]
        s0 = data.nodes['rdr_sig0'][inds].filled(fill_value=-10)
        #print(s0)
        wgt = data.nodes['wse_u'][inds]**2
        wse_t = truth.nodes['wse'][inds_t]
        for node in nodes:
            if (node in nodes and node in nodes_t):
                t_wse = wse_t[nodes_t==node]
                d_wse = wse[nodes==node]
                w = 1.0 / wgt[nodes==node]
                w[wgt[nodes==node]<0.00001] = 10000000
                diff_wse.append((d_wse - t_wse))
                weight.append(w)
                #print("wse: ",wse[nodes==node])
                #print("s0: ",s0[nodes==node])
                sig0.append(s0[nodes==node])
        # exclude outliers
        p_68 = np.nanpercentile(np.abs(diff_wse),68)
        diff_wse = np.array(diff_wse)
        diff_wse = diff_wse.flatten()
        diff_wse[abs(diff_wse) > 1.0e+11] = np.nan # remove fill values
        # 3-sigma is outlier
        outlier_msk = np.array([np.abs(d) > p_68 * 3 if ~np.isnan(d)
                                else False for d in diff_wse])
        weight = np.array(weight)
        diff_wse[outlier_msk] = np.nan
        weight[outlier_msk] = np.nan
        err_out[ind] = np.nanmedian(diff_wse)#np.nanmean((diff_wse * weight)) / np.nanmean(weight)
        sig0_out[ind] = np.nanmean(np.array(sig0))
    return err_out, sig0_out

def get_metrics(truth, data, msk=None, with_slope=True, with_width=True,
                with_wse_r_u=True, with_slope2=True, wse_node_avg=None):
    if msk is None:
        msk = np.ones(np.shape(data.wse),dtype=bool)
    metrics = {
        'area_total': (
            (data.area_total[msk] - truth.area_total[msk]) / truth.area_total[msk]) * 100.0,
        'area_detct':(
            #(data.area_detct[msk] - truth.area_detct[msk]) / truth.area_detct[msk]) * 100.0,
            (data.area_detct[msk] - truth.area_total[msk]) / truth.area_total[msk]) * 100.0,
        'wse': (data.wse[msk] - truth.wse[msk]) * 1e2,#convert m to cm

    }
    if wse_node_avg is not None:
        metrics['wse_node_avg'] = wse_node_avg[msk] * 1e2#convert m to cm
    if with_slope:
        metrics['slope'] = (data.slope[msk] - truth.slope[msk]) * 1e5#convert from m/m to cm/km
        metrics['slope_t'] = (truth.slope[msk]) * 1e5#convert from m/m to cm/km
    if with_slope2:
        metrics['slope2'] = (data.slope2[msk] - truth.slope2[msk]) * 1e5#convert from m/m to cm/km
        metrics['slope2_t'] = (truth.slope2[msk]) * 1e5#convert from m/m to cm/km

    if with_width:
        metrics['width'] = data.width[msk] - truth.width[msk]
    if with_wse_r_u:
        metrics['wse_r_u'] = data.wse_r_u[msk] * 1e2 #convert m to cm
        metrics['wse_t_r_u'] = truth.wse_r_u[msk] * 1e2 #convert m to cm
    return metrics

def get_passfail(is_lake = False):
    if not is_lake:
        passfail = {
            'area_tot e (%)': [15, 30],
            'wse e (cm)': [10, 20],
            'slp e (cm/km)':[1.7, 3.4],
            'slp2 e (cm/km)': [1.7, 3.4]
        }
    else:
        passfail = {
            'area_tot e (%)': [15, 30],
            'wse e (cm)': [10, 20],
        }
    return passfail

def get_truth_classes():
    truth_classes = {
        'bad_reach': [74292100251, 74292200011, 74291800011, 74100600051,
                      74100600061, 74100600071, 74100600091, 74100600551,
                      74100600571, 73260300021, 73270200031, 73270200051,
                      81130400071, 74267600131, 74267600121, 74265000141,
                      74267100031, 81140300021, 81140300031, 81140300041,
                      81140300061, 81140200011, 73240300011, 73240300021,
                      73240200201, 73218000031, 73218000471, 73218000271,
                      73218000441, 74262700351, 74262700231, 73240500151,
                      74269900361, 74269900781, 74269900491, 74269800011,
                      74269900011, 73160200101, 73216000211, 73214000011,
                      73150600541, 73150600171, 73150600031, 73150600011,
                      73150600021, 73150600151, 73150600161, 73150600951,
                      73150600111, 21602100101, 21602600131, 21602600191,
                      21602600311, 21602600371, 21602600871, 21602600891,
                      21602601451, 21602700051, 21602700071, 23221000021,
                      23267000111, 23267000131, 23267000201, 23267000271,
                      23267000371, 21602600241, 81130400011, 74292300011,
                      74292100211],
        'tribs': [74230900151, 74291800111, 74291700051, 74291800081,
                  74284300051, 74284300061, 74100600051, 74100600061,
                  74100600071, 74100600081, 74100600551, 74100600561,
                  74100600571, 73260300061, 73260100045, 73260200021,
                  73260300021, 73270200021, 73270200031, 73270200041,
                  73270200051, 81130400111, 74267600131, 74267600151,
                  74266300051, 74266300061, 74266300071, 74266300081,
                  74267100041, 81140300031, 81140300041, 81140300051,
                  81140300071, 73240200091, 73220700291, 73220700281,
                  73220900211, 73220900321, 73218000611, 73218000321,
                  73240900141, 73240500151, 73240200301, 74269900281,
                  74269900291, 74269900311, 74269900351, 74269700041,
                  74269800031, 74269800171, 74269800241, 74269800271,
                  74269800291, 73160200101, 73160200061, 74269600061,
                  73213000021, 73214000021, 73214000031, 73150600211,
                  73150600221, 73150600241, 73150600251, 73150600541,
                  73150600561, 73150600181, 73160100161, 73160100181,
                  73160100101, 73150600031, 73150600061, 73150600011,
                  73150600161, 73160100071, 78210000261, 23221000061,
                  21602100565, 21602700081, 21602600191],
        'non_linear': [74292300011, 74292100221, 74230900141, 74230900161,
                       74230900261, 74291900051, 74291700061, 74291800051,
                       74100600081, 74100600091, 73260300071, 73260200021,
                       73260300041, 73270200021, 81130400021, 81130400031,
                       81130400051, 81130400061, 81130400111, 81130400071,
                       74265000141, 74266300071, 74267100061, 81140300021,
                       81140300041, 81140300061, 81140300071, 81140300011,
                       73240200081, 73220700291, 73220700281, 73240200101,
                       73218000371, 73218000471, 73220900211, 73218000321,
                       73240500151, 73240500161, 73240200301, 74269900281,
                       74269900291, 74269900301, 74269900311, 74269900341,
                       74269900351, 74269900271, 74269900491, 74269800021,
                       74269800031, 74269800061, 74269800081, 74269800101,
                       74269800221, 74269800261, 74269900031, 74269700031,
                       73160200101, 73160200091, 73160200081, 74269600051,
                       73216000201, 73214000171, 73214000201, 73150600231,
                       73150600251, 73150600541, 73160100181, 73160100191,
                       73160300021, 73160100101, 73150600031, 73150600061,
                       73150600081, 73150600011, 73150600021, 73150600101,
                       73160100071, 73160100091],
        'edge_node': [74291800111, 81140300051, 81140300061, 81140300081,
                      73240900141, 74262700231, 73240200301, 74269900281,
                      74269900781, 74269700041, 73160200071, 73214000011,
                      73150600101, 73160100091, 23230200041, 21602100565],
        'partial_truth': [74292100221, 74292100251, 74292100261],
        'wrong_dir': [74292200061, 74230900241, 74230900261, 74230900271,
                      74291800011, 73260100045, 73270200031, 74267600131,
                      74265000141, 74266300031, 73240300021, 73240200201,
                      73218000011, 73218000041, 73218000051, 73218000351,
                      73218000361, 73218000381, 73218000391, 73218000401,
                      73218000411, 73218000081, 74262700211, 74262700231,
                      74269900361, 74269800041, 74269800171, 74269800181,
                      73216000021, 73216000211, 73216000051, 73213000011,
                      73214000021, 73150600021, 73150600111],
        'multi_chn': [74230900151, 74230900241, 74230900221, 81130400051,
                      74267100071, 74269700031, 73214000021, 21602600191,
                      21602300011],
        'linear': [73150600011, 73150600061, 73150600251, 73160100091,
                   73160100191, 73160200101, 73160300021, 73218000371,
                   73218000471, 73220700281, 73220700291, 73220900211,
                   73240200301, 74269800021, 74269800061, 74269800081,
                   74269800101, 74269800261, 74269900301, 74269900341,
                   74291900051, 74291900071, 81130400021, 81130400061,
                   81130400071, 81130400111, 81140300011, 81140300071,
                   21602600171, 21602100565, 23221000021, 23221000061,
                   21602600121, 21602600111]}
    return truth_classes


def compute_reach_fit_error(truth, scene, scene_nodes):
    fit_error = []
    if truth:
        for reach, scene_id in zip(truth.reaches['reach_id'], np.array(scene)):
            inds = (truth.nodes['reach_id'] == reach)
            inds = np.logical_and(truth.nodes['reach_id'] == reach, 
                np.array(scene_nodes) == scene_id)
            print("scene_nodes",np.array(scene_nodes)[inds], reach, truth.nodes['reach_id'][inds])
            print("reach, wse, node_id", reach, truth.nodes['wse'][inds], truth.nodes['node_id'][inds])
            ind = (truth.reaches['reach_id'] == reach)
            try:
                y0 = truth.nodes['wse'][inds]
                x0 = truth.nodes['node_id'][inds]
                # check masked array
                x1 = x0[y0.mask == False]
                y1 = y0[y0.mask == False]
                # handle bad values set to large fill value or anomolously large wse
                x2 = x1[np.abs(y1) < 1e11]
                y2 = y1[np.abs(y1) < 1e11]
                #print("reach, x1, y1",reach, x1, y1)
                # exclude nans and infs
                y = y2[np.isfinite(y2)]
                x = x2[np.isfinite(y2)]
                #print("reach, x, y",reach, x, y)
                z = np.polyfit( x, y, 1)
                p = np.poly1d(z)
                #err = np.nanmean(np.sqrt((y - p(x))**2))*100 #in cm
                err = np.nanmax(np.sqrt((y - p(x))**2))*100 #in cm
                #print("reach, z, err",reach, z, err)
                #if reach == 73216000741:#81140200151:#81130400121:
                #    print("reach,reach,x,y,z,err", reach, int(truth.reaches['reach_id'][ind][0]) ,x,y,z,err)
            except np.linalg.LinAlgError:
                err = 1000000000
                print("linAlgError caught. truth.nodes[wse]:",truth.nodes['wse'][inds])
            fit_error.append(err)
    return np.array(fit_error)
        
def mask_for_sci_req(metrics, truth, data, scene, scene_nodes=None, sig0=None,
                     print_table=False):
    # find reaches where the height profile linear fit is not that good
    # so we can filter out bogus/non-realistic reaches from the analysis
    fit_error = []#compute_reach_fit_error(truth, scene, scene_nodes)
    #print("p_length",truth.reaches['p_length'][truth.reaches['p_length']>0])
    #print("p_n_nodes",truth.reaches['p_n_nodes'][truth.reaches['p_n_nodes']>0]*200)
    # now make the mask
    # msk = np.logical_and(truth.reaches['width']>0, truth.reaches['area_total']>0)
    bounds = {
        'min_xtrk': 10000,
        'max_xtrk': 60000,
        'min_width': 80,
        'min_area': 800000,
        'min_length': 8000,
        'min_obs_frac': 0.5,
        'max_dark_frac': 1,
        'min_area_obs_frac': 0.2,
        'min_truth_ratio': 0.2,
        'min_xtrk_ratio': 1.0
    }
    msk=[]

    # define some extra masking criteria for each reach based on node values
    obs_area_frac = np.empty(np.size(data.reaches['reach_id']))
    # truth_ratio = np.empty(np.size(data.reaches['reach_id']))
    xtrk_ratio = np.empty(np.size(data.reaches['reach_id']))
    for index, reach in enumerate(data.reaches['reach_id']):
        reach_scene = scene[index]
        # can have multiple reaches with same reach id in dataset
        scene_mask = [s == reach_scene for s in scene_nodes]
        node_mask = np.logical_and(data.nodes['reach_id'] == reach,
                                   scene_mask)
        # node_truth_mask = np.logical_and(truth.nodes['reach_id'] == reach,
        #                                  scene_mask)
        n_good_data = np.sum(data.nodes['area_total'][node_mask] > 0)
        # n_good_truth = np.sum(truth.nodes['area_total'][node_truth_mask] > 0)
        n_prd = len(data.nodes['node_id'][node_mask])
        n_good_xtrk = np.sum(
            np.logical_and(np.abs(data.nodes['xtrk_dist'][node_mask]) > 10000,
                           np.abs(data.nodes['xtrk_dist'][node_mask]) < 60000))

        obs_area_frac[index] = n_good_data / n_prd
        # truth_ratio[index] = n_good_data / n_good_truth
        xtrk_ratio[index] = n_good_xtrk / n_prd

    # some quick plots
    # plt.hist(obs_area_frac)
    # plt.xlabel('Observed fraction for area')
    # plt.show()
    # plt.hist(truth_ratio)
    # plt.xlabel('n_good_nom / n_good_truth')
    # plt.show()
    # plt.hist(xtrk_ratio)
    # plt.xlabel('Proportion of nodes within xtrk bounds')
    # plt.show()
    # sc = plt.scatter(obs_area_frac, truth_ratio, c=metrics['area_total'],
    #                  cmap=plt.cm.coolwarm)
    # cbar = plt.colorbar(sc)
    # plt.clim(-40, 40)
    # cbar.set_label('area total error', rotation=270)
    # plt.xlabel('Observed fraction for area')
    # plt.ylabel('n_good_nom / n_good_truth')
    # plt.title('Truth ratio vs obs frac with area total error')
    # plt.show()
    # plt.scatter(truth_ratio, metrics['area_total'], cmap=plt.cm.coolwarm)
    # plt.clim(-40, 40)
    # plt.xlabel('n_good_nom / n_good_truth')
    # plt.ylabel('area total error (%)')
    # plt.title('Truth ratio vs area total error')
    # plt.ylim(-40, 40)
    # plt.show()
    # plt.scatter(obs_area_frac, metrics['area_total'], cmap=plt.cm.coolwarm)
    # plt.clim(-40, 40)
    # plt.xlabel('Area observed fraction')
    # plt.ylabel('area total error (%)')
    # plt.title('Area observed frac vs area total error')
    # plt.ylim(-40,40)
    # plt.show()

    if truth:
        msk = np.logical_and.reduce([
            np.abs(truth.reaches['xtrk_dist']) > bounds['min_xtrk'],
            np.abs(truth.reaches['xtrk_dist']) < bounds['max_xtrk'],
            truth.reaches['width'] > bounds['min_width'],
            truth.reaches['area_total'] > bounds['min_area'],
            truth.reaches['p_length'] >= bounds['min_length'],
            data.reaches['obs_frac_n'] >= bounds['min_obs_frac'],
            truth.reaches['dark_frac'] <= bounds['max_dark_frac']
        ])
        # add the node-level filters to the mask
        msk = np.logical_and.reduce([
            msk,
            obs_area_frac >= bounds['min_area_obs_frac'], #truth_ratio >= bounds['min_truth_ratio'],
            xtrk_ratio >= bounds['min_xtrk_ratio']
        ])
        if print_table:
            passfail = {
                'Truth width (m)': [bounds['min_width'], 'flip'],
                'Truth area_t (m^2)': [bounds['min_area'], 'flip'],
                'Prior length (m)': [bounds['min_length'], 'flip'],
                'Obs fraction (%)': [bounds['min_obs_frac'] * 100, 'flip'],
                'Xtrk dist (km)': [bounds['min_xtrk'], 'flip'],
                'Dark frac (%)':[bounds['max_dark_frac'] * 100,
                                 bounds['max_dark_frac'] * 100],
                'obs frac area (%)': [bounds['min_area_obs_frac']*100, 'flip'],
                'xtrk_ratio (%)': [bounds['min_xtrk_ratio']*100, 'flip']
            }
            preamble = "\nFor " + str(bounds['min_xtrk']) + " km<xtrk_dist<" \
                   + str(bounds['max_xtrk']) + " km and width>" \
                   + str(bounds['min_width']) + " m and area>" \
                   + str(bounds['min_area']) + " m^2 \n and reach len>=" \
                   + str(bounds['min_length']) + " m and obs frac >" \
                   + str(bounds['min_obs_frac']) + " and truth ratio > "\
                   + str(bounds['min_truth_ratio']) + " and xtrk proportion > "\
                   + str(bounds['min_xtrk_ratio'])
            table = {
                'Reach ID': data.reaches['reach_id'].astype(str),
                'Truth width (m)': truth.reaches['width'],
                'Truth area_t (m^2)': truth.reaches['area_total'],
                'Prior length (m)': truth.reaches['p_length'],
                'Obs fraction (%)': data.reaches['obs_frac_n']*100,
                'Xtrk dist (km)': truth.reaches['xtrk_dist'],
                'Dark frac (%)': data.reaches['dark_frac']*100,
                'River Name': data.reaches['river_name'],
                'sci req msk': msk,
                'obs frac area (%)': obs_area_frac*100,
                'xtrk_ratio (%)': xtrk_ratio*100
            }

            SWOTRiver.analysis.tabley.print_table(table, preamble=preamble,
                                                  fname=None, precision=6,
                                                  passfail=passfail)

        return msk, fit_error, bounds, truth.reaches['dark_frac'],\
               truth.reaches['p_length']
    else:
        return msk, fit_error, bounds, [], []

#
def get_scene_from_fnamedir(fnamedir):
    path_parts = os.path.abspath(fnamedir).split('/')

    # assumes particular directory structure
    scene0 = path_parts[-7]
    cycle_pass_tile = path_parts[-6].split('_')
    if scene0.startswith('cycle'):
        scene0 = path_parts[-8]
        cycle_pass_tile = path_parts[-7].split('_')

    if len(cycle_pass_tile) < 5:
        scene='unknown'
    else:
        scene = scene0 + "_" + cycle_pass_tile[3] + "_" \
                + cycle_pass_tile[4]
    return scene
#
def print_errors(metrics, msk=True, fname=None, preamble=None, with_slope=True,
                 with_slope2=True, with_node_avg=False):
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
    if with_node_avg:
        # get statistics of slope error
        wse_n_68 = np.nanpercentile(abs(metrics['wse_node_avg'][msk]), 68)
        wse_n_50 = np.nanpercentile(metrics['wse_node_avg'][msk], 50)
        wse_n_mean = np.nanmean(metrics['wse_node_avg'][msk])
        wse_n_num = np.count_nonzero(~np.isnan(metrics['wse_node_avg'][msk]))

        table['metric'].append('wse_node_avg (cm)')
        table['|68%ile|'].append(wse_n_68)
        table['50%ile'].append(wse_n_50)
        table['mean'].append(wse_n_mean)
        table['count'].append(wse_n_num)
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
    if with_slope2:
        slope2_68 = np.nanpercentile(abs(metrics['slope2'][msk]), 68)
        slope2_50 = np.nanpercentile(metrics['slope2'][msk], 50)
        slope2_mean = np.nanmean(metrics['slope2'][msk])
        slope2_num = np.count_nonzero(~np.isnan(metrics['slope2'][msk]))

        table['metric'].append('slope2 (cm/km)')
        table['|68%ile|'].append(slope2_68)
        table['50%ile'].append(slope2_50)
        table['mean'].append(slope2_mean)
        table['count'].append(slope2_num)
    SWOTRiver.analysis.tabley.print_table(table, preamble=preamble, fname=fname, precision=8)
    return table


def print_metrics(
        metrics, truth, scene=None, msk=None, fit_error=None,
        dark_frac=None, preamble=None, with_slope=True, with_slope2=True,
        with_width=True, with_node_avg=False, reach_len=None,
        with_wse_r_u=True, fname=None, passfail={}):
    table = {}
    if msk is None:
        msk = np.ones(np.shape(metrics['wse'][:]),dtype = bool)
    table['wse e (cm)'] = metrics['wse'][msk]
    if with_node_avg:
        table['wse node e (cm)'] = metrics['wse_node_avg'][msk]
    if with_wse_r_u:
        table['wse r u (cm)'] = metrics['wse_r_u'][msk]
        table['wse t r u (cm)'] = metrics['wse_t_r_u'][msk]
    if with_slope:
        table['slp e (cm/km)'] = metrics['slope'][msk]
        table['slope (cm/km)'] = metrics['slope_t'][msk]
    if with_slope2:
        table['slp2 e (cm/km)'] = metrics['slope2'][msk]
        table['slope2 (cm/km)'] = metrics['slope2_t'][msk]
    table['area_tot e (%)'] = metrics['area_total'][msk]
    table['area_det e (%)'] = metrics['area_detct'][msk]
    if with_width:
        table['width e (m)'] = metrics['width'][msk]
        table['width (m)'] = truth.reaches['width'][msk]
    else:
        table['sqrt(area) (m)'] = np.sqrt(truth['area_total'][msk])
    try:
        table['reach'] = [str(int(rid)) for rid in truth.reaches['reach_id'][msk]]
        table['xtrk (km)'] = truth.reaches['xtrk_dist'][msk]/1e3
    except AttributeError as e:
        table['lake_id'] = truth['obs_id'][msk]
        table['xtrk (km)'] = truth['xtrk_dist'][msk]/1e3
    if scene is not None:
        table['scene_pass_tile'] = np.array(scene)[msk]
    #if fit_error is not None:
    #    table['fit_error (cm)'] = np.array(fit_error)[msk]
    if dark_frac is not None:
        table['dark_frac'] = np.array(dark_frac)[msk]
    if reach_len is not None:
        table['reach_len (km)'] = np.array(reach_len/1e3)[msk]


    SWOTRiver.analysis.tabley.print_table(table, preamble=preamble, precision=5, passfail=passfail, fname=fname)
    return table
