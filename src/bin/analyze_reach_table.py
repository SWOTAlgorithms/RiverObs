#!/usr/bin/env python
'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s): Brent Williams

Read in the river tables and make some plots

'''
import pdb

import numpy as np
import matplotlib.pyplot as plt
import argparse
import SWOTRiver.analysis.tabley

plt.rcParams.update({'font.size': 12})


def get_percentiles(data):
    return np.nanpercentile(data, 16), np.nanpercentile(data, 50), \
           np.nanpercentile(data, 84)


def plot_hist(data1, xlab='', legend_label=['',], bins=100, data2=None):
    h1, bins1 = np.histogram(data1, bins)
    b = bins1[0:-1] + (bins1[1]-bins1[0])/2
    pcntle1 = get_percentiles(data1)
    legend_label[0] = legend_label[0] + '[16, 50, 84]%-ile = '\
                + '[%2.1f, %2.1f, %2.1f]'%pcntle1
    plt.figure()
    plt.plot(b, h1)
    if data2 is not None:
        pcntle2 = get_percentiles(data2)
        legend_label[1] = legend_label[1] + '[16, 50, 84]%-ile = '\
                   + '[%2.1f, %2.1f, %2.1f]'%pcntle2
        h2, _ = np.histogram(data2, bins1)
        plt.plot(b, h2)
    plt.grid()
    legend = plt.legend(legend_label)
    legend.set_draggable(state=True)
    plt.xlabel(xlab)
    plt.ylabel('count')


def plot_cdf(data, requirement, xlabel='', bins=250):
    plt.figure()
    abs_data = abs(data)
    plt.hist(abs_data, cumulative=True, label=xlabel+' CDF', histtype='step',
             bins=bins, density=True, linewidth=1.5)
    plt.xlabel('|' + xlabel + '|')
    plt.ylabel('Likelihood of occurrence')
    plt.axhline(y=0.68, color='r', linestyle='--', label='|68%ile|')
    plt.axvline(x=requirement, color='b', linestyle='--', label='requirement')
    x_min = -1
    x_max = np.percentile(abs_data, 98)
    plt.xlim([x_min, x_max])
    plt.grid()
    legend = plt.legend(loc='lower right')
    legend.set_draggable(state=True)
    plt.plot()


def plot_2d_hist(data, xdata, bins=100, xlab='',ylab=''):
    h, ye, xe = np.histogram2d(data, xdata, bins)
    xb = xe[0:-1] + (xe[1]-xe[0])/2
    #print(xb)
    #yb = ye[0:-1] + (ye[1]-ye[0])/2
    extent = (xe[0], xe[-1], ye[0], ye[-1])
    #print(extent)
    hp = h.copy()
    hp[h==0]=np.nan
    plt.figure()
    plt.imshow(hp, aspect='auto', interpolation='none', origin='lower',
               extent=extent)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.grid()
    # compute percentiles vs x-dim
    p_negsig = []
    med = []
    p_possig = []
    for i in range(len(xe)-1):
        msk = np.where(np.logical_and(xdata >= xe[i], xdata < xe[i+1]))
        #print(msk)
        #breakpoint()
        #dat = data[msk]
        pt = get_percentiles(data[msk])
        p_negsig.append(pt[0])
        med.append(pt[1])
        p_possig.append(pt[2])
    plt.plot(xb,p_negsig, '--r')
    plt.plot(xb,med,'r')
    plt.plot(xb,p_possig,'--r')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="""river errors table file with each
     observed reach itemized (not the summary stats table) as produced by  
     plot_reach_stats.py""", type=str)
    args = parser.parse_args()
    table = SWOTRiver.analysis.tabley.Table.from_file(args.input_file)
    wse_reach_e = np.array([line[0] for line in table.data])
    wse_node_e = np.array([line[1] for line in table.data])
    slp_e = np.array([line[4] for line in table.data])
    area_tot_e = np.array([line[8] for line in table.data])
    area_det_e = np.array([line[9] for line in table.data])
    width_e = np.array([line[10] for line in table.data])
    wse_bins = np.linspace(-100, 100, 50)
    slp_bins = np.linspace(-10, 10, 50)
    area_bins = np.linspace(-70, 70, 50)
    width_bins = np.linspace(-150, 150, 50)

    reach_ids = np.array([line[12] for line in table.data])
    scene_pass_tiles = np.array([line[14] for line in table.data])
    tiles = scene_pass_tiles
    plot_hist(
        wse_reach_e, 'wse error (cm)', ['reach, ', 'node, '], wse_bins,
        wse_node_e
    )
    plot_hist(slp_e, 'slope error (cm/km)', ['',], slp_bins)
    plot_hist(
        area_tot_e, 'area error (%)', ['tot, ', 'det, '], area_bins,
        area_det_e
    )
    plot_hist(width_e, 'river width (m)', ['',], width_bins)

    plot_cdf(wse_reach_e, 10,  'wse error (cm)')
    plot_cdf(slp_e, 1.7, 'slope error (cm/km)')
    plot_cdf(area_tot_e, 15, 'area error (%)')

    # plot 2d histograms vs river width
    plot_2d_hist(area_tot_e, width_e, [area_bins, width_bins],
                 'river width (m)', 'area_tot error (%)')
    plot_2d_hist(area_det_e, width_e, [area_bins, width_bins],
                 'river width (m)', 'area_det error (%)')
    plt.show()
    """
    h_wse, b_wse = np.histogram(wse, 100)
    h_wse_n, _ = np.histogram(wse_n, 100)
    bins_wse = b_wse[0:-1] + (b_wse[1]-b_wse[0])/2
    plt.figure()
    plt.plot(bins_wse, h_wse)
    plt.plot(bins_wse, h_wse_n)
    plt.xlabel('wse error (cm)')

    h_slp, b_slp = np.histogram(slp, 100)
    bins_slp = b_wse[0:-1] + (b_slp[1]-b_slp[0])/2
    plt.figure()
    plt.plot(bins_slp,h_slp)
    plt.xlabel('slp error (cm)')

    h_area_tot, b_area = np.histogram(area_tot, 100)
    h_area_det, _ = np.histogram(area_det, b_area)
    bins_area = b_area[0:-1] + (b_area[1]-b_area[0])/2
    plt.figure()
    plt.plot(bins_area,h_area_tot)
    plt.plot(bins_area,h_area_det)
    plt.xlabel('area % error (%)')
    plt.legend(['area_tot', 'area_det'])
    plt.show()
    """


if __name__ == '__main__':
    main()

