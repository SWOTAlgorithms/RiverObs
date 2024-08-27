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
import geopandas as gpd
import SWOTWater.products.product
import itertools


CALVAL_RIVERS = [
    'Connecticut River', 'Connecticut River; Westfield River',
    'North Saskatchewan River', 'Peace River', 'Ribdon River',
    'Sagavanirktok River', 'Slave River', 'Tanana River',
    'Willamette River', 'Yukon River', 'Colorado River', 'Porcupine River',
    'no_data', 'Santiam River', 'South Santiam River', 'Waimakariri River',
    'Rakaia River', 'Merrimack River', 'Quinebaug River', 'Mania', 'Lawa',
    "L'Aussonnelle", '181601', 'La Garonne', 'Tsiribihina', 'Maroni'
]


def get_input_files(basedir, pixc_run_id, river_run_id,
                    cycles=None, passes=None, tiles=None, pkl=None):
    print('Getting input files....')
    # Ensure cycles, passes, and tiles are lists and replace None with empty list
    cycles = [str(cycle).zfill(3) if cycle is not None else '' for cycle in
              (cycles if isinstance(cycles, list) else [cycles])]
    passes = [str(p).zfill(3) if p is not None else '' for p in
              (passes if isinstance(passes, list) else [passes])]
    tiles = [str(tile) if tile is not None else '' for tile in
             (tiles if isinstance(tiles, list) else [tiles])]    # Generate all possible combinations of cycle/pass/tile search strings
    search_strings = ['*'.join(map(str, combo)) for combo in
                      itertools.product(cycles, passes, tiles)]
    print('Cycle/pass/tile search strings are:', search_strings)

    rivertiles = []
    # uncomment below & comment out the globs for quick .nc test mode
    # rivertiles = ['/u/franz-r0/swot/calval-site-data/garonne/swot-data/016_076R/522/SWOT_L1B_HR_SLC_522_016_076R/SWOT_L2_HR_PIXC_522_016_076R/fineval_L3PGC0c_dev/SWOT_L2_HR_RiverTile_522_016_076R/SWOT_L2_HR_RiverTile_522_016_076R_fineval_L3PGC0c_dev_prd16_KCV-418_20240516/SWOT_L2_HR_RiverTile_522_016_076R_20230516T033113_20230516T033124_PGA1_01.nc']
    # river_fwd = False
    # uncomment below & comment out the globs for quick .shp test mode
    # rivertiles = [(
    # '/u/franz-r0/swot/swot-adt-data/009_228R/509/SWOT_L1B_HR_SLC_509_009_228R/SWOT_L2_HR_PIXC_509_009_228R/bulkreprocPGC0/SWOT_L2_HR_RiverTile_509_009_228R/SWOT_L2_HR_RiverTile_509_009_228R_bulkreprocPGC0/SWOT_L2_HR_RiverTile_Reach_509_009_228R_20230503T000035_20230503T000046_PGC0_01.shp',
    # '/u/franz-r0/swot/swot-adt-data/009_228R/509/SWOT_L1B_HR_SLC_509_009_228R/SWOT_L2_HR_PIXC_509_009_228R/bulkreprocPGC0/SWOT_L2_HR_RiverTile_509_009_228R/SWOT_L2_HR_RiverTile_509_009_228R_bulkreprocPGC0/SWOT_L2_HR_RiverTile_Node_509_009_228R_20230503T000035_20230503T000046_PGC0_01.shp')]
    # river_fwd = True

    for search_str in search_strings:
        if 'fwd' in pixc_run_id:
            rivertiles.extend(glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                                 'SWOT_L2_HR_PIXC_*/'
                                                 '/SWOT_L2_HR_RiverTile*/'
                                                 'SWOT_L2_HR_RiverTile*' +
                                                 river_run_id +
                                                 '/SWOT_L2_HR_RiverTile*' +
                                                 search_str + '*.nc',
                                       recursive=True))
        elif 'bulkreproc' in river_run_id:
            # we use shapefiles and there are no PIXCVecRiver
            rivertiles.extend(glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                                 'SWOT_L2_HR_PIXC_*/' +
                                                  pixc_run_id +
                                                 '/SWOT_L2_HR_RiverTile*/'
                                                 'SWOT_L2_HR_RiverTile*' +
                                                 river_run_id +
                                                 '/SWOT_L2_HR_RiverTile*' +
                                                 search_str + '*.shp',
                                       recursive=True))
            river_fwd = True
        else:
            rivertiles.extend(glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'
                                                 'SWOT_L2_HR_PIXC_*/' +
                                                  pixc_run_id +
                                                 '/SWOT_L2_HR_RiverTile*/'
                                                 'SWOT_L2_HR_RiverTile*' +
                                                 river_run_id +
                                                 '/SWOT_L2_HR_RiverTile*' +
                                                 search_str + '*.nc',
                                       recursive=True))
            river_fwd=False

    if len(rivertiles) == 0:
        raise Exception('No rivertile found, check input directory names')

    pixcvecs = np.empty(len(rivertiles), dtype=object)
    pixcs = np.empty(len(rivertiles), dtype=object)
    if not river_fwd:
        # get the associated pixcvecs and pixc files for each rivertile
        for index, rivertile in enumerate(rivertiles):
            print('RiverTile', rivertile)
            pixcvecs[index] = glob.glob(
                    '/' + os.path.join(*rivertile.split('/')[:-1])
                    + '/SWOT_L2_HR_PIXCVecRiver_*.nc'
            )[0]
            pixcs[index] = glob.glob(
                '/' + os.path.join(*rivertile.split('/')[:-3])
                + '/SWOT_L2_HR_PIXC_*.nc'
            )[0]
        if len(rivertiles) != len(pixcvecs):
            raise Exception('The number of rivertiles found doesnt match with '
                            'the number of pixcvecs found, some will be missing')
    else:
        # pair the reach and node files from a list to a list of tuples
        rivertiles = pair_files(rivertiles)
        # print('test mode') # comment above and uncomment this for shp test
    if pkl is not None:
        pt_node_file = pkl + '/dataframe/matched_pt_node_df.pkl'  # '/dataframe/matched_pt_node_df.pkl'
        pt_reach_file = pkl + '/dataframe/pt_reach_wse_df.pkl'
        pt_reach_matched_wse = pkl + '/dataframe/matched_pt_reach_wse_df.pkl'
        pt_matched_slope = pkl + '/dataframe/matched_pt_reach_slope_df.pkl'
        drift_matched_nodes = pkl + '/dataframe/matched_drift_node_df.pkl'
        drift_matched_reaches = pkl + 'dataframe/matched_drift_reach_df.pkl'
        error_dataframe = pkl + '/stats/ALL/pt_matched_reach_error_table_ALL.csv'
        coarse_nodes = pkl + '/stats/ALL/coarse_matched_node_error_table_ALL.csv'
        coarse_reaches = pkl + '/stats/ALL/offline_coarse_matched_reach_error_table_ALL.csv'
        field_dataframes = {'pt_node': pt_node_file,
                            'pt_reach': pt_reach_file,
                            'pt_match_wse': pt_reach_matched_wse,
                            'pt_match_slope': pt_matched_slope,
                            'drift_node_match': drift_matched_nodes,
                            'drift_reach_match': drift_matched_reaches}
        for key in field_dataframes.keys():
            if os.path.isfile(field_dataframes[key]):
                field_dataframes[key] = pd.read_pickle(field_dataframes[key])
            else:
                print('pkl input', field_dataframes[key],
                      'does not exist; check filenames!')
                field_dataframes[key] = None
        # add PT errors separately due to different file format
        if os.path.isfile(error_dataframe):
            field_dataframes['pt_error'] = pd.read_csv(error_dataframe)
        else:
            # try one other file nesting path; maybe reach doesn't exist
            if os.path.isfile(pkl + '/stats/node/ALL/pt_matched_node_error_table_ALL.csv'):
                field_dataframes['pt_error'] = pd.read_csv(pkl + '/stats/node/ALL/pt_matched_node_error_table_ALL.csv')
        if os.path.isfile(coarse_nodes):
            field_dataframes['coarse_node'] = pd.read_csv(coarse_nodes)
        elif os.path.isfile(pkl + '/dataframe/drift_node_df.pkl'):
            # try the cnes dataframe version
            field_dataframes['coarse_node'] = pd.read_csv(pkl + '/dataframe/drift_node_df.pkl')
        else:
            field_dataframes['coarse_node'] = None
        if os.path.isfile(coarse_reaches):
            field_dataframes['coarse_reach'] = pd.read_csv(coarse_reaches)
        else:
            field_dataframes['coarse_reach'] = None
    else:
        field_dataframes = None
    return rivertiles, pixcvecs, pixcs, field_dataframes

# Function to pair reach and node shapefiles together
def pair_files(files):
    paired_list = []
    # Create a dictionary to track pairs
    file_dict = {f.replace('Reach', 'Node') if 'Reach' in f else f.replace(
        'Node', 'Reach'): f for f in files}
    for file in files:
        # Check if the current file is a 'reach' file and has a corresponding 'node' file
        if 'Reach' in file:
            node_version = file.replace('Reach', 'Node')
            if node_version in file_dict:
                paired_list.append((file, node_version))
    return paired_list

# Function to read NetCDF
def read_netcdf(file_path):
    with nc.Dataset(file_path, 'r') as ncf:
        reach_ids = ncf['reaches']['reach_id'][:].filled(np.nan)
        reach_wse = ncf['reaches']['wse'][:].filled(np.nan)
        reach_width = ncf['reaches']['width'][:].filled(np.nan)
        reach_ids = reach_ids[~np.isnan(reach_wse) & ~np.isnan(reach_width)]
        river_names = ncf['reaches']['river_name'][:]
    river_dict = SWOTWater.products.product.MutableProduct.from_ncfile(
        file_path)
    return reach_ids, reach_wse, reach_width, river_names, river_dict

# Function to read Shapefile
def read_shapefile(file_path_tuple):
    reach_file, node_file = file_path_tuple
    reach_df = gpd.read_file(reach_file)
    node_df = gpd.read_file(node_file)
    node_df['node_id'] = node_df['node_id'].astype(int)
    node_df['reach_id'] = node_df['reach_id'].astype(int)
    reach_df['reach_id'] = reach_df['reach_id'].astype(int)
    reach_ids = reach_df['reach_id'].to_numpy()
    reach_wse = reach_df['wse'].to_numpy()
    reach_width = reach_df['width'].to_numpy()
    river_names = reach_df['river_name'].to_numpy()
    # Create a data container so it later reads like the netcdf data
    river_df = create_data_container(node_df, reach_df)

    return reach_ids, reach_wse, reach_width, river_names, river_df

def create_data_container(df_nodes, df_reaches):
    class DataContainer:
        def __init__(self, df_nodes, df_reaches):
            self.nodes = df_nodes
            self.reaches = df_reaches

        def __repr__(self):
            return f"DataContainer with nodes and reaches data"

    return DataContainer(df_nodes, df_reaches)

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
    parser.add_argument('-tag', '--tag_output', default='', type=str,
                        help='Tag to append to the directory for writing '
                             'outputs.')
    parser.add_argument('-pkl', '--pkl_input_dir',  default=None, type=str,
                        help='Input directory for pkl truth files. If included,'
                             'truth WSE/slope data are plotted alongside SWOT '
                             'profiles.')
    parser.add_argument('-cv', '--calval', action='store_true', default=False,
                        help='If flagged, only run reaches in the SWOT calval '
                             'set (NS, CR, WM, WK, PD, YR, TN, SG, PY, GC).')


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
        try:
            file_extension = os.path.splitext(rivertile)[1].lower()
            if file_extension == '.nc':
                reach_ids, reach_wse, reach_width, river_names, \
                river_df = read_netcdf(rivertile)
                river_file = rivertile
        except TypeError:
            if isinstance(rivertile, tuple):
                # input a list of tuples containing reach & node shapefiles
                reach_ids, reach_wse, reach_width, river_names, \
                river_df = read_shapefile(rivertile)
                river_file = rivertile[0]
            else:
                print('Unsupported file format')
        if args.calval:
            print(
                'Calval arg set by user; only processing calval rivers...'
            )
            if set(river_names) & set(CALVAL_RIVERS):
                print(
                    'Tile', rivertile, 'has a calval river. Continuing...'
                )
            else:
                print('Tile', rivertile, 'does not have a calval river.')
                continue
        dir_parts = river_file.split('/')[-1]
        file_parts = dir_parts.split('_')
        if out_dir is not None:
            last_dir = args.river_run_id + args.tag_output
            this_out_dir = f'{out_dir}/{args.pixc_run_id}/' \
                           f'{last_dir}'
            if not os.path.isdir(this_out_dir):
                os.umask(0)
                os.makedirs(this_out_dir, 0o777)
                os.makedirs(this_out_dir + '/truth/', 0o777)
                os.makedirs(this_out_dir + '/no_truth/', 0o777)
        else:
            this_out_dir = None

        for reach_id in reach_ids:
            if 'Reach' in file_parts:
                # input was a shapefile
                title = file_parts[5] + '_' + file_parts[6] + '_' + file_parts[
                    7] + '_' + str(reach_id) + '_' + args.river_run_id
            else:
                # input was a netcdf
                title = file_parts[4] + '_' + file_parts[5] + '_' + \
                        file_parts[6] + '_' + str(reach_id) + '_' \
                        + args.river_run_id
            if args.pkl_input_dir is not None:
                plot_reach.make_plots(river_file,
                    river_df, field_dataframes, pixcvec, pixc,
                    truth_pixcvec, truth_pixc, reach_id,
                    reach_error, nodes, pixc_truth, out_dir=this_out_dir,
                    title=title, overwrite=args.overwrite
                )
            else:
                # deprecated; may delete later
                plot_reach.make_plots(rivertile,
                    river_df, truth, pixcvec, pixc, truth_pixcvec,
                    truth_pixc, reach_id, reach_error, nodes,
                    pixc_truth, out_dir=this_out_dir, title=title,
                    overwrite=args.overwrite
                )

if __name__ == "__main__":
    main()
