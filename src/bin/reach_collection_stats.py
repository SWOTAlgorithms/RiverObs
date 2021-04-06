#!/usr/bin/env python
import argparse
import pdb
import matplotlib
import pandas as pd
import xarray as xr
import numpy as np
import seaborn as sns

from plotnine import *
from plot_reach_stats import load_and_accumulate
from reach_comparison import *


def mask_for_sci_req_local(metrics, truth, data, scene, scene_nodes=None, sig0=None):
    # This function masks for the scientific requirements locally (i.e. reach_collection_stats.py ONLY),
    # so we can tweak the histograms independent of the defaults set in SWOTRiver/analysis/riverobs.py

    # uncomment below line to find reaches where the height profile linear fit is not that good
    # so we can filter out bogus/non-realistic reaches from the analysis
    bad_reaches = []
    fit_error = []  # SWOTRiver.analysis.riverobs.compute_reach_fit_error(truth, scene, scene_nodes)
    msk = np.logical_and((np.abs(truth.reaches['xtrk_dist']) > 10000),
          np.logical_and((np.abs(truth.reaches['xtrk_dist']) < 60000),
          np.logical_and((truth.reaches['width'] > 80),
          np.logical_and((truth.reaches['area_total'] > 8e5),
          np.logical_and(np.isin(truth.reaches['reach_id'], bad_reaches, invert=True),
          np.logical_and((truth.reaches['p_length'] >= 8000),
          np.logical_and(truth.reaches['obs_frac_n'] >= 0.5,
                         truth.reaches['dark_frac'] < 1.0)))))))
    return msk, fit_error, truth.reaches['dark_frac']


def load_data_df(data_files, truth_files=None, test_bool=None):
    # takes all input filenames and forms some dataframes out of their data and metrics
    # TO DO: make an independent to_dataframe or dictionary method for the mutable products and use that instead
    data_nodes_df = pd.DataFrame()
    data_reaches_df = pd.DataFrame()
    truth_nodes_df = pd.DataFrame()
    truth_reaches_df = pd.DataFrame()
    metrics_df = pd.DataFrame()
    test_counter = 0
    truth = 0
    if type(truth_files) is not list:  # single file case
        data_files = [data_files]
        truth_files = [truth_files]
    for index, filename in enumerate(data_files):
        if test_counter < 5:
            # get the error of that scene
            try:
                metrics, truth, data, scene, scene_nodes, sig0 = load_and_accumulate(filename, truth_files[index])
            except FileNotFoundError:
                print('\032[93mData rivertile', filename, 'has no matching truth rivertile',
                      truth_files[index], '\032[0m')
            if truth:
                msk, fit_error, dark_frac = mask_for_sci_req_local(metrics, truth, data, scene, sig0=sig0)
                if not any(msk):
                    print('\n\033[93mNo reaches in file', filename, 'are within sci req bounds\n')
                    try:
                        print('For reach_ids', int(truth.reaches['reach_id']),
                              '\nTruth xtrks', truth.reaches['xtrk_dist'],
                              '\nWidths', truth.reaches['width'], 'm',
                              '\nArea tots', truth.reaches['area_total'], 'm^2',
                              '\nReach lengths', truth.reaches['p_length'], 'm',
                              '\nObserved fracs', truth.reaches['obs_frac_n'],
                              '\nDark fracs', truth.reaches['dark_frac'], '\033[0m')
                    except TypeError:
                        print('No reaches in file.\033[0m')
                else:
                    print('\nload and accumulate successful for', filename)
                    # struggled with the product object here, so I restructured into a pandas dataframe  -- Cassie
                    print('populating variables for dataframes...')
                    temp_node_data = pd.DataFrame()
                    temp_node_truth = pd.DataFrame()
                    for var in data.nodes.variables:
                        # get data for nodes
                        temp_node_truth[var] = truth.nodes[var]
                        temp_node_data[var] = data.nodes[var]
                    temp_node_data['scene'] = scene_nodes
                    temp_node_data['msk'] = data.nodes['wse'].mask
                    temp_node_truth['scene'] = scene_nodes
                    temp_node_truth['msk'] = truth.nodes['wse'].mask
                    data_nodes_df = data_nodes_df.append(temp_node_data, ignore_index=True)
                    truth_nodes_df = truth_nodes_df.append(temp_node_truth, ignore_index=True)

                    temp_reach_data = pd.DataFrame()
                    temp_reach_truth = pd.DataFrame()
                    for var in data.reaches.variables:
                        # get data for reaches
                        if var not in ['centerline_lat', 'centerline_lon', 'rch_id_up', 'rch_id_dn']:
                            temp_reach_truth[var] = truth.reaches[var]
                            temp_reach_data[var] = data.reaches[var]

                    temp_reach_data['msk'] = msk
                    temp_reach_data['scene'] = scene
                    temp_reach_truth['msk'] = msk
                    temp_reach_truth['scene'] = scene
                    temp_reach_data['sig0'] = sig0
                    data_reaches_df = data_reaches_df.append(temp_reach_data, ignore_index=True)
                    truth_reaches_df = truth_reaches_df.append(temp_reach_truth, ignore_index=True)

                    metrics_tmp_df = pd.DataFrame(metrics)
                    metrics_tmp_df['reach_id'] = data.reaches['reach_id']
                    metrics_tmp_df['scene'] = scene
                    metrics_df = metrics_df.append(metrics_tmp_df, ignore_index=True)
            if test_bool:  # only increment this counter if user is in 'test' mode
                test_counter = index

    return data_nodes_df, truth_nodes_df, data_reaches_df, truth_reaches_df, metrics_df


def combine_truth_and_data(data_df, truth_df):
    data_df['type'] = 'data'
    truth_df['type'] = 'truth'
    df = data_df.append(truth_df, ignore_index=True)
    # clean up fill values
    df = df[(df['wse'] > -999999999999) &
            (df['area_detct'] > -999999999999) &
            (df['area_total'] > -999999999999)]
    return df


def make_hist(node_df, node_df_truth, reach_df, reach_df_truth, reach_metrics_df, node_metrics_df, variables, title_str):
    # creates histograms for the distribution of data at the node and the reach level
    # use sci req mask only
    reach_df = reach_df[reach_df['msk'] == True]
    reach_df_truth = reach_df_truth[reach_df_truth['msk'] == True]
    node_df = node_df[node_df['msk'] == False]
    node_df_truth = node_df_truth[node_df_truth['msk'] == False]

    reach_df_comb = combine_truth_and_data(reach_df, reach_df_truth)
    node_df_comb = combine_truth_and_data(node_df, node_df_truth)
    for var in variables:
        print('Creating', var, 'histogram...')
        if var in reach_df.columns:
            anno_text = "data median is " + str(reach_df[var].quantile(0.5)) + \
                "\ntruth median is " + str(reach_df_truth[var].quantile(0.5))
            print(anno_text)
            g = (
                    ggplot(reach_df_comb)
                    + aes(x=var, fill='type')
                    + geom_histogram(alpha=0.5, bins=50)
                    + geom_vline(xintercept=[reach_df[var].quantile(0.5), reach_df_truth[var].quantile(0.5)],
                                 colour=['red', 'green'],
                                 linetype='dotted')
                    + labs(title=title_str + " Reach-level " + var)
                    + annotate("text", label=anno_text)
            )
            print(g)
        if var in node_df.columns:
            anno_text = "data median is " + str(node_df[var].quantile(0.5)) + \
                        "\ntruth median is " + str(node_df_truth[var].quantile(0.5))
            x_min = min(node_df_comb[var])
            x_max = node_df_comb[var].quantile(0.95)
            g = (
                    ggplot(node_df_comb)
                    + aes(x=var, fill='type')
                    + geom_histogram(alpha=0.5, bins=100)
                    + labs(title=title_str + " Node-level " + var)
                    + annotate("text", label=anno_text)
                    + coord_cartesian(xlim=(x_min, x_max))
            )
            print(g)
        if var in reach_metrics_df.columns:
            g = (
                    ggplot(reach_metrics_df[reach_metrics_df[var] < reach_metrics_df[var].quantile(0.90)])
                    + aes(x=var, fill='scene')
                    + geom_histogram(alpha=0.5, bins=100)
                    + labs(title=title_str + " Reach-level " + var + " error by scene")
            )
            print(g)
        if var in node_metrics_df.columns:
            var_metrics_df = node_metrics_df[node_metrics_df[var].abs() > 0]
            # x_bound = var_metrics_df[var].abs().quantile(0.9)
            g = (
                    ggplot(var_metrics_df[var_metrics_df[var].abs() < var_metrics_df[var].quantile(0.90)])
                    + aes(x=var, fill='scene')
                    + geom_histogram(alpha=0.5, bins=100)
                    + labs(title=title_str + " Node-level " + var + " error") # + coord_cartesian(xlim=(-1*x_bound, x_bound))
            )
            print(g)


def get_node_errors(node_df, node_df_truth):
    print('getting node-level errors...')
    node_metrics = pd.DataFrame()
    node_metrics['node_id'] = node_df['node_id']
    node_metrics['reach_id'] = node_df['reach_id']
    node_metrics['scene'] = node_df['scene']
    node_metrics['area_total'] = ((node_df['area_total'] - node_df_truth['area_total'])/node_df['area_total']) * 100.0
    node_metrics['area_detct'] = ((node_df['area_detct'] - node_df_truth['area_detct'])/node_df['area_detct']) * 100.0
    node_metrics['wse'] = (node_df['wse'] - node_df_truth['wse']) * 1e2   # convert m to cm
    node_metrics['width'] = node_df['width'] - node_df_truth['width']
    node_metrics['lat'] = node_df['lat'] - node_df_truth['lat']
    node_metrics['lon'] = node_df['lon'] - node_df_truth['lon']
    # remove very large error values
    node_metrics = node_metrics[(abs(node_metrics['wse']) < 3000) &
                                (abs(node_metrics['area_total']) < 999) &
                                (abs(node_metrics['area_detct']) < 999)]
    return node_metrics


def combine_metrics(node_df, node_df_truth, reach_df, reach_df_truth, node_metrics, reach_metrics):
    print('combining all river dataframes...')
    # remove masked reaches/nodes
    reach_df = reach_df[reach_df['msk'] == True]
    reach_df_truth = reach_df_truth[reach_df_truth['msk'] == True]
    node_df = node_df[node_df['msk'] == False]
    node_df_truth = node_df_truth[node_df_truth['msk'] == False]
    # rename error columns for uniqueness
    node_metrics = node_metrics.rename(columns={'wse': 'node_wse_e',
                                                'area_total': 'node_area_t_e',
                                                'area_detct': 'node_area_d_e',
                                                'width': 'node_width_e',
                                                'lat': 'node_lat_e',
                                                'lon': 'node_lon_e'})
    reach_metrics = reach_metrics.rename(columns={'wse': 'reach_wse_e',
                                                  'slope': 'slope_e',
                                                  'area_total': 'reach_area_t_e',
                                                  'area_detct': 'reach_area_d_e',
                                                  'width': 'reach_width_e',
                                                  'lat': 'reach_lat_e',
                                                  'lon': 'reach_lon_e'})
    # combine all dataframes
    all_metrics = pd.merge(node_metrics, reach_metrics, on=['reach_id', 'scene'])
    all_data = pd.merge(left=node_df, right=reach_df, on=['reach_id', 'scene'], suffixes=('_node', '_reach'))
    all_truth = pd.merge(left=node_df_truth, right=reach_df_truth,
                         on=['reach_id', 'scene'], suffixes=('_node', '_reach'))
    data_and_truth = pd.merge(left=all_data, right=all_truth, on=['reach_id', 'node_id', 'scene'],
                              suffixes=('_data', '_truth'))
    river_metrics = pd.merge(left=all_metrics, right=data_and_truth,
                             on=['reach_id', 'scene', 'node_id'],
                             suffixes=('_met', '_d'))
    return river_metrics


def plot_correlation_matrix(river_metrics):
    # get correlations and plot matrix
    river_matrix = river_metrics.corr()
    river_matrix = river_matrix.dropna(axis=0, how='all')
    river_matrix = river_matrix.dropna(axis=1, how='all')
    print('error correlations...')
    print_coeffs(river_matrix)

    mask = np.triu(np.ones_like(river_matrix, dtype=bool))
    cmap = sns.diverging_palette(250, 15, s=75, l=40, n=9, center="light", as_cmap=True)
    plt.figure(figsize=(16, 12))

    sns.heatmap(river_matrix, mask=mask, center=0, annot=False,
                    fmt='.2f', square=True, cmap=cmap)
    plt.show()

    # plot smaller matrix for error types of interest
    error_matrix = river_matrix[['reach_wse_e', 'node_wse_e', 'slope_e', 'dark_frac_reach_truth', 'dark_frac_node_truth',
                                'node_area_d_e', 'reach_area_d_e', 'node_area_t_e', 'reach_area_t_e',
                                'node_dist_node_data', 'node_dist_reach_data', 'node_lat_e', 'node_lon_e'
                                 ]].sort_values(by=['reach_wse_e'], ascending=False).transpose()
    h = sns.heatmap(error_matrix, cmap=cmap, xticklabels=True, yticklabels=True, square=True)
    h.set_yticklabels(h.get_yticklabels(), fontsize=12)
    plt.show()

    return river_matrix


def plot_xy(dataframe, var1, var2, title_str):
    print('Plotting x y', var1, var2)
    try:
        g = (
                ggplot(dataframe)
                + aes(x=var1, y=var2, fill='scene')
                + geom_point(alpha=0.5)
                + labs(title=title_str + ' ' + var1 + ' vs ' + var2 + ' by scene')
        )
        print(g)
    except KeyError:
        print('Input plot_xy variables not found in dataframe. Check variable names and try again.')


def print_coeffs(corr_matrix):
    error_vars = ['reach_wse_e', 'node_wse_e', 'slope_e',
                  'reach_area_t_e', 'node_area_t_e',
                  'reach_area_d_e', 'node_area_d_e',
                  'reach_width_e', 'node_width_e',
                  'node_lat_e', 'node_lon_e']
    for var in error_vars:
        print('\033[93mHighest', var, 'correlations: \n')
        print(corr_matrix[var].nlargest(10), '\n\033[0m')
        print('\033[92mLowest', var, 'correlations: \n')
        print(corr_matrix[var].nsmallest(10), '\n\033[0m')


def get_collection_node_error(datas, truths, title_str):
    scene_node_errors = []
    for index, scene in enumerate(datas):
        data = datas[index]
        truth = truths[index]
        scene_node_errors = np.append(scene_node_errors,
                                      SWOTRiver.analysis.riverobs.compute_average_node_error(data, truth))
    collection_avg_node_error = np.mean(abs(scene_node_errors))
    collection_med_node_error = np.median(scene_node_errors)
    plt.figure()
    plt.hist(scene_node_errors[scene_node_errors<100], bins=100, alpha=0.5,
             label='Average node errors for each scene in collection')
    plt.legend(loc='upper right')
    plt.xlabel('Average node wse error of scene')
    plt.title(title_str + ' Node WSE')
    plt.show()
    print('average (abs) wse node error of the collection is', collection_avg_node_error)
    print('median wse node error of the collection is', collection_med_node_error)


def get_collection_errors(river_metrics):
    # compute average errors for whole collection of rivertiles
    print('mean absolute wse reach error is', round(river_metrics['reach_wse_e'].abs().mean(), 2))
    print('median wse reach error is', round(river_metrics['reach_wse_e'].median(), 2))
    print('mean absolute slope reach error is', round(river_metrics['slope_e'].abs().mean(), 2))
    print('median slope reach error is', round(river_metrics['slope_e'].median(), 2))
    print('mean absolute area total reach error is', round(river_metrics['reach_area_t_e'].abs().mean(), 2))
    print('median area total reach error is', round(river_metrics['reach_area_t_e'].median(), 2))
    print('mean absolute area detected reach error is', round(river_metrics['reach_area_d_e'].abs().mean(), 2))
    print('median area detected reach error is', round(river_metrics['reach_area_d_e'].median(), 2))
    river_metrics = river_metrics.drop_duplicates('reach_id')
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
    reach_wse_e = river_metrics['reach_wse_e']
    slope_e = river_metrics['slope_e']
    area_t_e = river_metrics['reach_area_t_e']
    area_d_e = river_metrics['reach_area_d_e']

    fig, [(ax0, ax1), (ax2, ax3)] = plt.subplots(2, 2, sharey=False, tight_layout=True)
    ax0.hist(reach_wse_e, 100)
    ax0.axvline(passfail['wse e (cm)'][0], color='k', linestyle='dashed', linewidth=1)
    ax0.axvline(-1*passfail['wse e (cm)'][0], color='k', linestyle='dashed', linewidth=1)
    ax0.axvline(passfail['wse e (cm)'][1], color='r', linestyle='dashed', linewidth=1)
    ax0.axvline(-1*passfail['wse e (cm)'][1], color='r', linestyle='dashed', linewidth=1)
    ax0.set_title('Reach wse errors', fontsize=12)
    ax0.set_xlabel('wse error, cm')
    wse_percent_good = 100*(sum(reach_wse_e.abs() < passfail['wse e (cm)'][0]))/len(reach_wse_e)
    wsestr = '% of reaches that meet scientific requirements = ' + str(round(wse_percent_good, 2))
    # get number of reaches
    num_wse = 'n reaches=' + str(sum(reach_wse_e.abs() < passfail['wse e (cm)'][0]))
    ax0.text(left, top, wsestr,
             horizontalalignment='left',
             verticalalignment='bottom',
             fontsize=8,
             color='k',
             transform=ax0.transAxes)
    # ax0.text(left, top-0.1, num_wse,
    #          horizontalalignment='left',
    #          verticalalignment='bottom',
    #          fontsize=8,
    #          color='k',
    #          transform=ax0.transAxes)

    ax1.hist(slope_e, 100)
    ax1.set_title('Reach Slope Errors', fontsize=12)
    ax1.set_xlabel('Slope error, cm/km')
    ax1.axvline(passfail['slp e (cm/km)'][0], color='k', linestyle='dashed', linewidth=1)
    ax1.axvline(-1*passfail['slp e (cm/km)'][0], color='k', linestyle='dashed', linewidth=1)
    ax1.axvline(passfail['slp e (cm/km)'][1], color='r', linestyle='dashed', linewidth=1)
    ax1.axvline(-1*passfail['slp e (cm/km)'][1], color='r', linestyle='dashed', linewidth=1)
    slope_percent_good = 100*sum(slope_e.abs() < passfail['slp e (cm/km)'][0])/len(slope_e)
    slopestr = '% of reaches that meet scientific requirements = ' + str(round(slope_percent_good, 2))
    ax1.text(left, top, slopestr,
             horizontalalignment='left',
             verticalalignment='bottom',
             fontsize=8,
             color='k',
             transform=ax1.transAxes)

    ax2.hist(area_t_e[abs(area_t_e) < area_t_max], 100)
    ax2.set_title('Reach Area total errors', fontsize=12)
    ax2.set_xlabel('Area total errors, %')
    ax2.axvline(passfail['area_tot e (%)'][0], color='k', linestyle='dashed', linewidth=1)
    ax2.axvline(-1*passfail['area_tot e (%)'][0], color='k', linestyle='dashed', linewidth=1)
    ax2.axvline(passfail['area_tot e (%)'][1], color='r', linestyle='dashed', linewidth=1)
    ax2.axvline(-1*passfail['area_tot e (%)'][1], color='r', linestyle='dashed', linewidth=1)
    area_t_percent_good = 100*sum(area_t_e.abs() < passfail['area_tot e (%)'][0])/len(area_t_e)
    areastr = '% of reaches that meet scientific requirements = ' + str(round(area_t_percent_good, 2))

    # getting 68%iles
    print('68%ile wse error (cm) is', river_metrics['reach_wse_e'].abs().quantile(0.68),
          'median is', river_metrics['reach_wse_e'].median())
    print('68%ile total slope error (cm/km) is', river_metrics['slope_e'].abs().quantile(0.68),
          'median is', river_metrics['slope_e'].median())
    print('68%ile total area error (%) is', river_metrics['reach_area_t_e'].abs().quantile(0.68),
          'median is', river_metrics['reach_area_t_e'].median())
    print('68%ile detected area error (%) is', river_metrics['reach_area_d_e'].abs().quantile(0.68),
          'median is', river_metrics['reach_area_d_e'].median())

    ax2.text(left, top, areastr,
             horizontalalignment='left',
             verticalalignment='bottom',
             fontsize=8,
             color='k',
             transform=ax2.transAxes)

    ax3.hist(area_d_e[abs(area_d_e) < area_d_max], 100)
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
    parser.add_argument('-v', '--show_vars', type=str, nargs='+',
                        default=['wse', 'area_total', 'area_detct', 'slope'],
                        help='Variables to plot. Default is wse, area, slope')
    parser.add_argument('-xy', '--plot_xy', type=str, nargs='+', default=None,
                        help='Plot x vs y of any two variables. e.g. \'-xy p_length slope_e\'')
    parser.add_argument('--test_boolean', help='set to "True" if testing script', default=False, required=False)
    parser.add_argument('--title', help='Title of Rivertile Set', required=False)
    parser.add_argument('-c', '--corr', type=bool, help='Plot correlation matrix (True or False)', default=False)
    args = parser.parse_args()

    # get or create title for rivertile set
    if args.title is None:
        title_str = args.basedir + '_' + args.slc_basename + '_' + args.pixc_basename + '_' + args.proc_rivertile
        print('Title is', title_str)
    else:
        title_str = args.title

    # get all rivertiles
    data_files, truth_files = get_input_files(args.basedir, args.slc_basename, args.pixc_basename, args.proc_rivertile,
                                           args.truth_rivertile)

    # load and accumulate data
    node_df, node_df_truth, reach_df, reach_df_truth, reach_metrics_df = load_data_df(data_files, truth_files,
                                                                                      args.test_boolean)

    # get node-level errors
    node_metrics_df = get_node_errors(node_df, node_df_truth)

    # combine data and error dataframes
    river_metrics = combine_metrics(node_df, node_df_truth, reach_df, reach_df_truth, node_metrics_df, reach_metrics_df)

    # plot correlation matrix
    if args.corr:
        corr_matrix = plot_correlation_matrix(river_metrics)

    # make xy plots (uncomment and populate if you want to see individual x, y relationships)
    if args.plot_xy:
        for index, var1 in enumerate(args.plot_xy):
            if index % 2 == 0:
                var2 = args.plot_xy[index+1]
                plot_xy(river_metrics, var1, var2, title_str)

    # get distribution of each result
    make_hist(node_df, node_df_truth, reach_df, reach_df_truth, reach_metrics_df,
              node_metrics_df, args.show_vars, title_str)

    # print general stats
    #get_collection_node_error(node_df, node_df_truth, title_str)
    get_collection_errors(river_metrics)

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
