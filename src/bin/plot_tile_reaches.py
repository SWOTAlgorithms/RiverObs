#!/usr/bin/env python

"""

This script runs the reach dashboard tool given an input rivertile and PIXC
RUN_ID from the standard swot-adt-data directory structure. It will generate as
many reach dashboard plots as there are unique reach/pass/tile/cycle
combinations. It will save them to the specified output directory in PNG
format.
"""

import os
import argparse
import pdb
import glob
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import pandas as pd
import plot_reach


def get_input_files(basedir, pixc_run_id, river_run_id,
                    cycles=None, passes=None, tiles=None, pkl=None):
    if 'fwd' in pixc_run_id:
        rivertiles = glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                         'SWOT_L2_HR_PIXC_*/'
                                         '/SWOT_L2_HR_RiverTile*/'
                                         'SWOT_L2_HR_RiverTile*' +
                                         river_run_id +
                                         '/SWOT_L2_HR_RiverTile*.nc',
                               recursive=True)
    else:
        rivertiles = glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                         'SWOT_L2_HR_PIXC_*/' +
                                          pixc_run_id +
                                         '/SWOT_L2_HR_RiverTile*/'
                                         'SWOT_L2_HR_RiverTile*' +
                                         river_run_id +
                                         '/SWOT_L2_HR_RiverTile*.nc',
                               recursive=True)
    # get the associated pixcvecs and pixc files for each rivertile
    pixcvecs = np.empty(len(rivertiles), dtype=object)
    pixcs = np.empty(len(rivertiles), dtype=object)
    for index, rivertile in enumerate(rivertiles):
        pixcvecs[index] = glob.glob(
                '/' + os.path.join(*rivertile.split('/')[:-1])
                + '/SWOT_L2_HR_PIXCVecRiver_*.nc'
        )[0]
        pixcs[index] = glob.glob(
            '/' + os.path.join(*rivertile.split('/')[:-3])
            + '/SWOT_L2_HR_PIXC_*.nc'
        )[0]
    if len(rivertiles) == 0:
        raise Exception('No rivertile found, check input directory names')
    if len(rivertiles) != len(pixcvecs):
        raise Exception('The number of rivertiles found doesnt match with '
                        'the number of pixcvecs found, some will be missing')

    if pkl is not None:
        pt_node_file = pkl + '/dataframe/matched_pt_node_df.pkl'
        pt_reach_file = pkl + '/dataframe/pt_reach_wse_df.pkl'
        pt_reach_matched_wse = pkl + '/dataframe/matched_pt_reach_wse_df.pkl'
        pt_matched_slope = pkl + '/dataframe/matched_pt_reach_slope_df.pkl'
        drift_matched_nodes = pkl + '/dataframe/matched_drift_node_df.pkl'
        error_dataframe = pkl + '/table/reach_metric_table_pt.csv'
        field_dataframes = {'pt_node': pt_node_file,
                            'pt_reach': pt_reach_file,
                            'pt_match_wse': pt_reach_matched_wse,
                            'pt_match_slope': pt_matched_slope,
                            'drift_node_match': drift_matched_nodes}
        for key in field_dataframes.keys():
            if os.path.isfile(field_dataframes[key]):
                field_dataframes[key] = pd.read_pickle(field_dataframes[key])
            else:
                print('pkl input', field_dataframes[key],
                      'does not exist; check filenames!')
                field_dataframes[key] = None
        # add PT errors separately due to different file format
        field_dataframes['pt_error'] = pd.read_csv(error_dataframe)
    else:
        field_dataframes = None
    return rivertiles, pixcvecs, pixcs, field_dataframes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rivertile_dir', help='rivertile directory')
    parser.add_argument('--out_dir', help='output dir', default=None)
    parser.add_argument('--pixc_run_id', help='pixc_run_id', default=None)
    parser.add_argument('--river_run_id', help='river_run_id', default=None)
    parser.add_argument('-p', '--passes', help='list of passes', nargs='+',
                        type=str.lower, default=None)
    parser.add_argument('-c', '--cycles', help='list of cycles', nargs='+',
                        type=str.lower, default=None)
    parser.add_argument('-t', '--tiles', help='list of tiles, e.g. 038R',
                        nargs='+', type=str.upper, default=None)
    parser.add_argument('-o', '--overwrite', action='store_true', default=False,
                        help='Overwrite existing output files')
    parser.add_argument('-pkl', '--pkl_input_dir',  default=None, type=str,
                        help='Input directory for pkl truth files. If included,'
                             'truth WSE/slope data are plotted alongside SWOT '
                             'profiles.')

    args = parser.parse_args()
    out_dir = args.out_dir
    truth_pixcvec = None
    truth_pixc = None
    pixc_truth = None
    truth = None
    reach_error = None
    nodes = None

    # get input files
    print('PASSES: ', args.passes)
    print('TILES: ', args.tiles)
    print('CYCLES: ', args.cycles)
    rivertiles, pixcvecs, pixcs, field_dataframes = get_input_files(
        args.rivertile_dir,
        args.pixc_run_id,
        args.river_run_id,
        args.cycles,
        args.passes,
        args.tiles,
        args.pkl_input_dir
    )
    for rivertile, pixcvec, pixc in zip(rivertiles, pixcvecs, pixcs):
        if os.path.isfile(rivertile):
            ncf = nc.Dataset(rivertile, 'r')
            # TODO: rewrite the below to avoid deprecation warnings
            reach_ids = ncf['reaches']['reach_id'][:].filled(np.nan)
            reach_wse = ncf['reaches']['wse'][:].filled(np.nan)
            reach_width = ncf['reaches']['width'][:].filled(np.nan)
            reach_ids = reach_ids[~np.isnan(reach_wse) & ~np.isnan(reach_width)]

            ncf.close()

            dir_parts = rivertile.split('/')[-1]
            file_parts = dir_parts.split('_')
            if out_dir is not None:
                this_out_dir = f'{out_dir}/{args.pixc_run_id}/' \
                               f'{args.river_run_id}'
                if not os.path.isdir(this_out_dir):
                    os.umask(0)
                    os.makedirs(this_out_dir, 0o777)
                    os.makedirs(this_out_dir + '/truth/', 0o777)
                    os.makedirs(this_out_dir + '/no_truth/', 0o777)
            else:
                this_out_dir = None

            for reach_id in reach_ids:
                title = file_parts[4] + '_' + file_parts[5] + '_' + \
                        file_parts[6] + '_' + str(reach_id)
                if args.pkl_input_dir is not None:
                    plot_reach.make_plots(
                        rivertile, field_dataframes, pixcvec, pixc,
                        truth_pixcvec, truth_pixc, reach_id,
                        reach_error, nodes, pixc_truth, out_dir=this_out_dir,
                        title=title, overwrite=args.overwrite
                    )
                else:
                    # deprecated; may delete later
                    plot_reach.make_plots(
                        rivertile, truth, pixcvec, pixc, truth_pixcvec,
                        truth_pixc, reach_id, reach_error, nodes,
                        pixc_truth, out_dir=this_out_dir, title=title,
                        overwrite=args.overwrite
                    )
        else:
            print('Input file', rivertile, 'does not exist')


if __name__ == "__main__":
    main()
