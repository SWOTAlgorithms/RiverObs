#!/usr/bin/env python
import argparse
import pdb
import matplotlib
import numpy as np
from reach_comparison import *


def plot_area_vs_slope(errors):
    plt.figure()
    slopes = [el[3] for el in errors]
    areas = [el[5] for el in errors]
    plt.plot(abs(slopes), abs(areas), 'o')
    plt.xlabel('slope error (%)')
    plt.ylabel('area total error (%)')
    print('n is', len(areas))
    plt.show()


def mask_for_sci_req_old(metrics, truth, data, scene, scene_nodes=None, sig0=None):
    # find reaches where the height profile linear fit is not that good
    # so we can filter out bogus/non-realistic reaches from the analysis
    fit_error = []  # SWOTRiver.analysis.riverobs.compute_reach_fit_error(truth, scene, scene_nodes)
    bad_reaches = [73150600261, 73150600551, 73150601011, 73160300011, 73216000231, 73218000071, 73218000251,
                   73220700271, 73220900221, 73220900231, 73240100211, 73240200041, 73240200231, 74230900191,
                   74230900201, 74230900271, 74262700251, 74265000131, 74266300011, 74269800121, 74291700051,
                   74291700071, 74291900011, 74292100271, 81130400071]
    msk = np.logical_and((np.abs(truth.reaches['xtrk_dist']) > 15000),
          np.logical_and((np.abs(truth.reaches['xtrk_dist']) < 55000),
          np.logical_and((truth.reaches['width'] > 100),
          np.logical_and((truth.reaches['area_total'] > 8e5),
          np.logical_and(np.isin(truth.reaches['reach_id'], bad_reaches, invert=True),
          np.logical_and((truth.reaches['p_length'] >= 8000),
          np.logical_and(truth.reaches['obs_frac_n'] >= 1.0, truth.reaches['dark_frac'] < 0.9)))))))
    return msk, fit_error, truth.reaches['dark_frac'], truth.reaches['p_n_nodes'] * 200.0

def load_data(pixc_files, gdem_files, test_bool):
    data_collection = []
    truth_collection = []
    msk_collection = []
    test_counter = 0
    truth = 0
    for index, filename in enumerate(pixc_files):
        if test_counter < 5:
            # get the error of that scene
            try:
                metrics, truth, data, scene, scene_nodes, sig0 = load_and_accumulate(pixc_files[index], gdem_files[index])
            except FileNotFoundError:
                print('Pixc rivertile has no matching gdem rivertile')
            if truth:
                msk, fit_error, dark_frac, reach_len = mask_for_sci_req_old(metrics, truth, data, scene, sig0=sig0) #SWOTRiver.analysis.riverobs.
                data_collection.append(data)
                truth_collection.append(truth)
                msk_collection.append(msk)
                if not any(msk):
                    print('No reaches in file', pixc_files[index], 'are within sci req bounds')
                else:
                    print('load and accumulate successful for', pixc_files[index])
            if test_bool: # only increase if in test mode
                test_counter=index

    return data_collection, truth_collection, msk_collection


def dark_frac_hist(pixc_data, gdem_data, msks, title_str):
    data_fracs = np.array([])
    truth_fracs = np.array([])
    for index, scene in enumerate(pixc_data):
        msk = msks[index]
        data_fracs = np.append(data_fracs, pixc_data[index].reaches.dark_frac[msk])
        truth_fracs = np.append(truth_fracs, gdem_data[index].reaches.dark_frac[msk])
    plt.figure()
    plt.hist(data_fracs, bins=50, alpha=0.5, label='data dark frac')
    plt.hist(truth_fracs, bins=50, alpha=0.5, label='truth dark frac')
    plt.legend(loc='upper right')
    plt.xlabel('Dark fraction')
    plt.title(title_str + 'Collection Dark Fractions')
    plt.show()


def slope_hist(pixc_data, gdem_data, msks, title_str):
    data_slopes = np.array([])
    truth_slopes = np.array([])
    for index, scene in enumerate(pixc_data):
        msk = msks[index]
        data_slopes = np.append(data_slopes[data_slopes<0.002], pixc_data[index].reaches.slope[msk])
        truth_slopes = np.append(truth_slopes[truth_slopes<0.002], gdem_data[index].reaches.slope[msk])
    plt.figure()
    plt.hist(data_slopes, bins=50, alpha=0.5, label='data slopes')
    plt.hist(truth_slopes, bins=50, alpha=0.5, label='truth slopes')
    plt.legend(loc='upper right')
    plt.xlabel('slope, m/m')
    plt.title(title_str + ' Slopes')
    plt.show()


def wse_hist(pixc_data, gdem_data, msks, title_str):
    data_wse = np.array([])
    truth_wse = np.array([])
    for index, scene in enumerate(pixc_data):
        msk = msks[index]
        data_wse = np.append(data_wse, pixc_data[index].reaches.wse[msk])
        truth_wse = np.append(truth_wse, gdem_data[index].reaches.wse[msk])
    plt.figure()
    plt.hist(data_wse, bins=50, alpha=0.5, label='data reach wse')
    plt.hist(truth_wse, bins=50, alpha=0.5, label='truth reach wse')
    plt.legend(loc='upper right')
    plt.xlabel('wse, cm')
    plt.title(title_str + ' wse')
    plt.show()


def area_hist(pixc_data, gdem_data, msks, title_str):
    data_area = np.array([])
    truth_area = np.array([])
    for index, scene in enumerate(pixc_data):
        msk = msks[index]
        data_area = np.append(data_area, pixc_data[index].reaches.area_total[msk])
        truth_area = np.append(truth_area, gdem_data[index].reaches.area_total[msk])
    plt.figure()
    plt.hist(data_area, bins=50, alpha=0.5, label='data reach area total')
    plt.hist(truth_area, bins=50, alpha=0.5, label='truth reach area total')
    plt.legend(loc='upper right')
    plt.xlabel('area total (incl dark water), m^2')
    plt.title(title_str + ' Area Total')
    plt.show()


def area_dtct_hist(pixc_data, gdem_data, msks, title_str):
    data_area = np.array([])
    truth_area = np.array([])
    for index, scene in enumerate(pixc_data):
        msk = msks[index]
        data_area = np.append(data_area, pixc_data[index].reaches.area_detct[msk])
        truth_area = np.append(truth_area, gdem_data[index].reaches.area_detct[msk])
    plt.figure()
    plt.hist(data_area, bins=50, alpha=0.5, label='data reach area detected')
    plt.hist(truth_area, bins=50, alpha=0.5, label='truth reach area detected')
    plt.legend(loc='upper right')
    plt.xlabel('area detected, m^2')
    plt.title(title_str + ' Area Detected')
    plt.show()


def get_collection_node_error(datas, truths, title_str):
    scene_node_errors = []
    for index, scene in enumerate(datas):
        data = datas[index]
        truth = truths[index]
        scene_node_errors = np.append(scene_node_errors, SWOTRiver.analysis.riverobs.compute_average_node_error(data, truth))
    collection_avg_node_error = np.mean(abs(scene_node_errors))
    collection_med_node_error = np.median(scene_node_errors)
    plt.figure()
    plt.hist(scene_node_errors[scene_node_errors<100], bins=100, alpha=0.5, label='Average node errors for each scene in collection')
    plt.legend(loc='upper right')
    plt.xlabel('Average node wse error of scene')
    plt.title(title_str + ' Node WSE')
    plt.show()
    print('average (abs) wse node error of the collection is', collection_avg_node_error)
    print('median wse node error of the collection is', collection_med_node_error)


def get_collection_errors(datas, truths, msks):
    metrics = []
    wse_errors = []
    slope_errors = []
    area_total_errors = []
    area_detct_errors = []
    for index, scene in enumerate(datas):
        msk = msks[index]
        truth = truths[index]
        data = datas[index]
        metric = SWOTRiver.analysis.riverobs.get_metrics(truth.reaches, data.reaches, msk)
        metrics.append(metric)
        wse_errors = np.append(wse_errors, metric['wse'])
        slope_errors = np.append(slope_errors, metric['slope'])
        area_total_errors = np.append(area_total_errors, metric['area_total'])
        area_detct_errors = np.append(area_detct_errors, metric['area_detct'])

    # compute average errors for whole collection of rivertiles
    print('mean absolute wse reach error is', str(round(np.mean(abs(wse_errors)), 2)))
    print('median wse reach error is', str(round(np.median(wse_errors), 2)))
    print('mean absolute slope reach error is', str(round(np.mean(abs(slope_errors)), 3)))
    print('median slope reach error is', str(round(np.median(slope_errors), 3)))
    print('mean absolute area total reach error is', str(round(np.mean(abs(area_total_errors)), 2)))
    print('median area total reach error is', str(round(np.median(area_total_errors), 2)))
    print('mean absolute area detected reach error is', str(round(np.mean(abs(area_detct_errors)), 2)))
    print('median area detected reach error is', str(round(np.median(area_detct_errors), 2)))

    # compute proportion that do not meet science requirements
    passfail = {
        'area_tot e (%)': [15, 30],
        'wse e (cm)': [10, 20],
        'slp e (cm/km)': [1.7, 3.4],
    }
    # plot histograms of errors
    wse_max = 1000
    slope_max = 100
    area_t_max = 100
    area_d_max = 100
    left, width = .05, .5
    bottom, height = .1, .85
    top = bottom + height

    fig, [(ax0, ax1), (ax2, ax3)] = plt.subplots(2, 2, sharey=False, tight_layout=True)
    ax0.hist(wse_errors[abs(wse_errors) < wse_max], 100)
    ax0.axvline(passfail['wse e (cm)'][0], color='k', linestyle='dashed', linewidth=1)
    ax0.axvline(-1*passfail['wse e (cm)'][0], color='k', linestyle='dashed', linewidth=1)
    ax0.axvline(passfail['wse e (cm)'][1], color='r', linestyle='dashed', linewidth=1)
    ax0.axvline(-1*passfail['wse e (cm)'][1], color='r', linestyle='dashed', linewidth=1)
    ax0.set_title('Reach wse errors', fontsize=12)
    ax0.set_xlabel('wse error, cm')
    wse_percent_good = 100*len(wse_errors[abs(wse_errors)<passfail['wse e (cm)'][0]])/len(wse_errors)
    wsestr = '% of reaches that meet scientific requirements = ' + str(round(wse_percent_good,2))
    num_wse = 'n reaches=' + str(len(wse_errors))
    ax0.text(left, top, wsestr,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=8,
              color='k',
              transform=ax0.transAxes)
    ax0.text(left, top-0.1, num_wse,
             horizontalalignment='left',
             verticalalignment='bottom',
             fontsize=8,
             color='k',
             transform=ax0.transAxes)

    ax1.hist(slope_errors[abs(slope_errors)<slope_max], 100)
    ax1.set_title('Reach Slope Errors', fontsize=12)
    ax1.set_xlabel('Slope error, cm/km')
    ax1.axvline(passfail['slp e (cm/km)'][0], color='k', linestyle='dashed', linewidth=1)
    ax1.axvline(-1*passfail['slp e (cm/km)'][0], color='k', linestyle='dashed', linewidth=1)
    ax1.axvline(passfail['slp e (cm/km)'][1], color='r', linestyle='dashed', linewidth=1)
    ax1.axvline(-1*passfail['slp e (cm/km)'][1], color='r', linestyle='dashed', linewidth=1)
    slope_percent_good = 100*len(slope_errors[abs(slope_errors)<passfail['slp e (cm/km)'][0]])/len(slope_errors)
    slopestr = '% of reaches that meet scientific requirements = ' + str(round(slope_percent_good,2))
    num_slope = 'num reaches=' + str(len(slope_errors))
    ax1.text(left, top, slopestr,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=8,
              color='k',
              transform=ax1.transAxes)

    ax2.hist(area_total_errors[abs(area_total_errors)<area_t_max], 100)
    ax2.set_title('Reach Area total errors', fontsize=12)
    ax2.set_xlabel('Area total errors, %')
    ax2.axvline(passfail['area_tot e (%)'][0], color='k', linestyle='dashed', linewidth=1)
    ax2.axvline(-1*passfail['area_tot e (%)'][0], color='k', linestyle='dashed', linewidth=1)
    ax2.axvline(passfail['area_tot e (%)'][1], color='r', linestyle='dashed', linewidth=1)
    ax2.axvline(-1*passfail['area_tot e (%)'][1], color='r', linestyle='dashed', linewidth=1)
    area_t_percent_good = 100*len(area_total_errors[abs(area_total_errors)<passfail['area_tot e (%)'][0]])\
                          /len(area_total_errors)
    num_a_t = 'num reaches=' + str(len(area_total_errors))
    areastr = '% of reaches that meet scientific requirements = ' + str(round(area_t_percent_good,2))

    # getting 68%iles
    wse_med = np.median(wse_errors)
    wse_ordered = np.sort(abs(wse_errors))
    index = int(np.floor(0.68 * len(wse_errors)))
    print('68%ile wse error (cm) is', str(round(wse_ordered[index],2)), 'median is', str(round(wse_med,2)))

    area_med = np.median(area_total_errors)
    area_ordered = np.sort(abs(area_total_errors))
    index = int(np.floor(0.68 * len(area_ordered)))
    print('68%ile total area error (%) is', str(round(area_ordered[index], 2)), 'median is', str(round(area_med, 2)))

    slope_med = np.median(slope_errors)
    slope_ordered = np.sort(abs(slope_errors))
    index = int(np.floor(0.68 * len(slope_ordered)))
    print('68%ile total slope error (cm/km) is', str(round(slope_ordered[index],2)), 'median is', str(round(slope_med,2)))

    ax2.text(left, top, areastr,
              horizontalalignment='left',
              verticalalignment='bottom',
              fontsize=8,
              color='k',
              transform=ax2.transAxes)

    ax3.hist(area_detct_errors[abs(area_detct_errors)<area_d_max], 100)
    ax3.set_title('Reach Area detected errors', fontsize=12)
    ax3.set_xlabel('Area detected errors, %')

    plt.show()


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
    parser.add_argument('--title', help='Title of Rivertile Set', required=False)
    args = parser.parse_args()

    # get or create title for rivertile set
    if args.title is None:
        title_str = args.basedir + '_' + args.slc_basename + '_' + args.pixc_basename + '_' + args.proc_rivertile
        print('Title is', title_str)
    else:
        title_str = args.title

    # get all rivertiles
    pixc_files, gdem_files = get_input_files(args.basedir, args.slc_basename, args.pixc_basename, args.proc_rivertile,
                                           args.truth_rivertile)
    # load and accumulate data
    datas, truths, msks = load_data(pixc_files, gdem_files, args.test_boolean)

    # get histogram of slopes
    dark_frac_hist(datas, truths, msks, title_str)
    slope_hist(datas, truths, msks, title_str)
    wse_hist(datas, truths, msks, title_str)
    area_hist(datas, truths, msks, title_str)
    area_dtct_hist(datas, truths, msks, title_str)

    # print general stats
    get_collection_node_error(datas, truths, title_str)
    get_collection_errors(datas, truths, msks)

    # get highest slope error cases
    # scene_errors, reach_errors, pixc_list = get_errors(pixc_files, gdem_files, args.test_boolean)
    # reach_slp_errors = sort_slope_errors(reach_errors, pixc_list)
    # plot_area_vs_slope(reach_slp_errors)

    # for index, reach in enumerate(reach_slp_errors):
    #     rivertile_file = reach[0]
    #     reach_id = reach[2]
    #     errors = [reach[3], reach[4], reach[5], reach[6], reach[7]] # slope e, wse_e, area total e, area detected e
    #     scene = SWOTRiver.analysis.riverobs.get_scene_from_fnamedir(rivertile_file)
    #     gdem_file = get_gdem_from_pixc(rivertile_file)
    #     figure, axes = plot_reach.make_plots(rivertile_file, reach_id, gdem_file, errors, scene)
    #
    #     plt.show()


if __name__ == "__main__":
    main()