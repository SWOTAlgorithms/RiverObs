"""
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s): Dustin Lagoy

"""
import warnings

import numpy as np

import SWOTRiver.products.product
import SWOTRiver.analysis.tabley


def load_rivertiles(truth_file, data_file):
    truth = SWOTRiver.products.product.MutableProduct.from_ncfile(truth_file)
    data = SWOTRiver.products.product.MutableProduct.from_ncfile(data_file)
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


def get_metrics(truth, data):
    metrics = {
        'area': (
            (data.area_total - truth.area_total) / truth.area_total) * 100.0,
        'area_detct':(
            (data.area_detct - truth.area_detct) / truth.area_detct) * 100.0,
        'height': (data.height - truth.height) * 1e2,
        'slope': (data.slope - truth.slope) / 10,
        'width': data.width - truth.width,
    }
    return metrics

def mask_for_sci_req(metrics, truth):
    metrics_msk={}
    for key in metrics:
        msk = np.logical_and((np.abs(truth.reaches['xtrk_dist'])>10000),
              np.logical_and((np.abs(truth.reaches['xtrk_dist'])<60000), 
                             (truth.reaches['width']>100)))
        metrics_msk[key] = metrics[key][msk]
    return metrics_msk

def print_errors(metrics):
    # get statistics of area error
    area_68 = np.nanpercentile(abs(metrics['area']), 68)
    area_50 = np.nanpercentile(metrics['area'], 50)
    area_mean = np.nanmean(metrics['area'])
    area_num = np.count_nonzero(~np.isnan(metrics['area']))
    # get statistics of height error
    height_68 = np.nanpercentile(abs(metrics['height']), 68)
    height_50 = np.nanpercentile(metrics['height'], 50)
    height_mean = np.nanmean(metrics['height'])
    height_num = np.count_nonzero(~np.isnan(metrics['height']))
    # get statistics of slope error
    slope_68 = np.nanpercentile(abs(metrics['slope']), 68)
    slope_50 = np.nanpercentile(metrics['slope'], 50)
    slope_mean = np.nanmean(metrics['slope'])
    slope_num = np.count_nonzero(~np.isnan(metrics['slope']))

    table = {
        'metric': ['area (%)', 'height (cm)', 'slope (cm/km)'],
        '|68%ile|': [area_68, height_68, slope_68],
        '50%ile': [area_50, height_50, slope_50],
        'mean': [area_mean, height_mean, slope_mean],
        'count': [area_num, height_num, slope_num],
    }
    SWOTRiver.analysis.tabley.print_table(table, precision=8)


def print_metrics(metrics, truth, scene):
    table = {}
    table['hgt e (cm)'] = metrics['height']
    table['slp e (cm/km)'] = metrics['slope']
    table['area e (%)'] = metrics['area']
    table['wid e (m)'] = metrics['width']
    table['width (m)'] = truth.reaches['width']
    table['reach'] = truth.reaches['reach_id']
    table['xtrk (km)'] = truth.reaches['xtrk_dist']/1e3
    table['scene'] = scene
    
    SWOTRiver.analysis.tabley.print_table(table, precision=5)
