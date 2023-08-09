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
import plot_reach


def get_input_files(basedir, pixc_run_id, river_run_id,
                    cycles=None, passes=None, tiles=None):
    if 'fwd' in pixc_run_id:
        rivertiles = glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                         'SWOT_L2_HR_PIXC_*/'
                                         '/SWOT_L2_HR_RiverTile*/'
                                         'SWOT_L2_HR_RiverTile*' +
                                         river_run_id +
                                         '/SWOT_L2_HR_RiverTile*.nc',
                               recursive=True)
        pixcvecs = glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                       'SWOT_L2_HR_PIXC*/' 
                                       '/SWOT_L2_HR_RiverTile*/'
                                       'SWOT_L2_HR_RiverTile*' + river_run_id +
                                       '/SWOT_L2_HR_PIXCVecRiver_*.nc',
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
        pixcvecs = glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                       'SWOT_L2_HR_PIXC*/' +
                                       pixc_run_id +
                                       '/SWOT_L2_HR_RiverTile*/'
                                       'SWOT_L2_HR_RiverTile*' + river_run_id +
                                       '/SWOT_L2_HR_PIXCVecRiver_*.nc',
                             recursive=True)

    if len(rivertiles) == 0:
        raise Exception('No rivertile found, check input directory names')
    if len(rivertiles) != len(pixcvecs):
        raise Exception('The number of rivertiles found doesnt match with '
                        'the number of pixcvecs found')

    pixcvecs_final = []
    rivertiles_final = []
    for pixcvec, rivertile in zip(pixcvecs, rivertiles):
        file_parts = pixcvec.split('/')[-1].split('_')
        cycle, pass_id, tile = file_parts[4], file_parts[5], file_parts[6]
        if cycles is not None and cycle not in cycles:
            continue
        if passes is not None and pass_id not in passes:
            continue
        if tiles is not None and tile not in tiles:
            continue
        rivertiles_final.append(rivertile)
        pixcvecs_final.append(pixcvec)
    if len(pixcvecs_final) == 0:
        raise Exception('No match files found')
    return rivertiles_final, pixcvecs_final


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

    args = parser.parse_args()
    out_dir = args.out_dir
    truth_tile = None
    truth_pixcvec = None
    truth_pixc = None
    pixc_truth = None
    gdem_dem = None
    reach_error = None
    nodes = None
    pixc = None

    # get input files
    print('PASSES: ', args.passes)
    print('TILES: ', args.tiles)
    print('CYCLES: ', args.cycles)
    rivertiles, pixcvecs = get_input_files(args.rivertile_dir,
                                           args.pixc_run_id,
                                           args.river_run_id,
                                           args.cycles,
                                           args.passes,
                                           args.tiles)
    for rivertile, pixcvec in zip(rivertiles, pixcvecs):
        if os.path.isfile(rivertile):
            ncf = nc.Dataset(rivertile, 'r')
            reach_ids = ncf['reaches']['reach_id'][:].filled(np.nan)
            reach_wse = ncf['reaches']['wse'][:].filled(np.nan)
            reach_width = ncf['reaches']['width'][:].filled(np.nan)
            reach_ids = reach_ids[~np.isnan(reach_wse) & ~np.isnan(reach_width)]
            ncf.close()

            dir_parts = rivertile.split('/')[-1]
            file_parts = dir_parts.split('_')
            for reach_id in reach_ids:
                title = file_parts[4] + '_' + file_parts[5] + '_' + \
                        file_parts[6] + '_' + str(reach_id)
                plot_reach.make_plots(rivertile, truth_tile, pixcvec, pixc,
                                      truth_pixcvec, truth_pixc, reach_id,
                                      gdem_dem, reach_error, nodes, pixc_truth)
                if out_dir is not None:
                    this_out_dir = f'{out_dir}/{args.pixc_run_id}/' \
                                   f'{args.river_run_id}'
                    if not os.path.isdir(this_out_dir):
                        os.umask(0)
                        os.makedirs(this_out_dir, 0o777)
                    plt.title(title, backgroundcolor='white')
                    plt.savefig(this_out_dir + '/' + title)
                    plt.close()
                else:
                    plt.title(title, backgroundcolor='white')
                    plt.show()
        else:
            print('Input file', rivertile, 'does not exist')


if __name__ == "__main__":
    main()
