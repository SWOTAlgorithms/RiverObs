#!/usr/bin/env python

"""

This script runs the reach dashboard tool given an input rivertile and PIXC
RUN_ID from the standard swot-adt-data directory structure. It will generate as
many reach dashboard plots as there are unique reach/pass/tile/cycle
combinations. It will save them to the specified output directory in PNG
format.
"""

import os
import re
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
import warnings

# for multitemporal stuff in plots
import rivscale.products.along_stretch

CALVAL_RIVERS = [
    'Connecticut River', 'Connecticut River; Westfield River',
    'North Saskatchewan River', 'Peace River', 'Ribdon River',
    'Sagavanirktok River', 'Slave River', 'Tanana River',
    'Willamette River', 'Yukon River', 'Colorado River', 'Porcupine River',
    'no_data', 'Santiam River', 'South Santiam River', 'Waimakariri River',
    'Rakaia River', 'Merrimack River', 'Quinebaug River', 'Mania', 'Lawa',
    "L'Aussonnelle", '181601', 'La Garonne', 'Tsiribihina', 'Maroni'
]

def load_mt_stats_files(mt_basedir, reach_id, mt_flavor='v1'):
    """
    search for specific reach 
    """
    basename = os.path.join(
        mt_basedir,
        '{}'.format(reach_id),
        '{}'.format(mt_flavor),
        '{}'.format(reach_id))
    #
    wse_file = basename + '_wse_stats.nc'
    mt_wse = None
    if os.path.isfile(wse_file):
        # try to read it
        mt_wse = rivscale.products.along_stretch.AlongStretchStats.from_ncfile(
                wse_file)
    #
    width_file = basename + '_width_stats.nc'
    mt_width = None
    if os.path.isfile(width_file):
        # try to read it
        mt_width = rivscale.products.along_stretch.AlongStretchStats.from_ncfile(
                width_file)
    return mt_wse, mt_width

def get_input_files(basedir, pixc_run_id, river_run_id, slc_run_id=None,
                    cycles=None, passes=None, tiles=None):
    print('Getting input files....')
    # Ensure cycles, passes, and tiles are lists. Replace None with empty list
    cycles = [str(cycle).zfill(3) if cycle is not None else '' for cycle in
              (cycles if isinstance(cycles, list) else [cycles])]
    passes = [str(p).zfill(3) if p is not None else '' for p in
              (passes if isinstance(passes, list) else [passes])]
    tiles = [str(tile) if tile is not None else '' for tile in
             (tiles if isinstance(tiles, list) else [tiles])]
    search_strings = ['*'.join(map(str, combo)) for combo in
                      itertools.product(cycles, passes, tiles)]
    print('Cycle/pass/tile search strings are:', search_strings)
    rivertiles = []

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
        elif 'local' in river_run_id:
            # Files were output to a local directory, outside the usual
            # directory structure. We walk until we find all RiverTiles (since
            # there should only be one version).
            for root, dirs, files in os.walk(basedir, followlinks=True):
                for file in files:
                    if file.startswith("SWOT_L2_HR_RiverTile_") and file.endswith(
                        ".nc"):
                        full_path = os.path.join(root, file)
                        if os.path.isfile(full_path):
                            rivertiles.append(full_path)
            river_fwd = False
        else:
            if slc_run_id is None:
                slc_str = ''
            else:
                slc_str = '{}/'.format(slc_run_id)
            rivertiles.extend(glob.glob(basedir + '/**/SWOT_L1B_HR_SLC*/'+
                                                 slc_str +
                                                 'SWOT_L2_HR_PIXC_*/' +
                                                  pixc_run_id +
                                                 '/SWOT_L2_HR_RiverTile*/'
                                                 'SWOT_L2_HR_RiverTile*' +
                                                 river_run_id +
                                                 '/SWOT_L2_HR_RiverTile*' +
                                                 search_str + '*.nc',
                                       recursive=True))
            river_fwd = False
    if len(rivertiles) == 0:
        raise Exception('No rivertile found, check input directory names')
    pixcvecs = np.empty(len(rivertiles), dtype=object)
    pixcs = np.empty(len(rivertiles), dtype=object)
    if not river_fwd:
        for index, rivertile in enumerate(rivertiles):
            print('RiverTile:', rivertile)
            pixcvecs[index] = get_pixcvecriver_from_rivertile(rivertile, river_run_id)
            pixcs[index] = get_pixc_from_rivertile(rivertile, river_run_id)
    else:
        rivertiles = pair_files(rivertiles)
    return rivertiles, pixcvecs, pixcs


def get_input_files_from_csv(csv_file, pixc_run_id, river_run_id, pkl=None):
    df = pd.read_csv(csv_file)
    if (('river' in df.keys()) and
            ('pixc' in df.keys()) and
            ('pixcvec' in df.keys())):
        rivertiles = df['river']
        pixcs = df['pixc']
        pixcvecs = df['pixcvec']
        reach_ids = df['reach_id']
    else:
        # Read CSV using the header row so it is not interpreted as data
        df = pd.read_csv(csv_file, header=0)
        # Use the first column by its header name
        col_name = df.columns[0]
        rivertiles = df[col_name].tolist()
        print(f"Loaded {len(rivertiles)} rivertiles from CSV: {csv_file}")
        pixcvecs = np.empty(len(rivertiles), dtype=object)
        pixcs = np.empty(len(rivertiles), dtype=object)
        for index, rivertile in enumerate(rivertiles):
            print('RiverTile:', rivertile)
            # Call our helper functions to fetch the corresponding files.
            pixcvecs[index] = get_pixcvecriver_from_rivertile(rivertile, river_run_id)
            pixcs[index] = get_pixc_from_rivertile(rivertile, river_run_id)
        reach_ids = [None,] * len(rivertiles)
    return rivertiles, pixcvecs, pixcs, reach_ids


def get_pixcvecriver_from_rivertile(rivertile, river_run_id):
    """
    Extract cycle/pass/tile info from the rivertile file name and return
    the first matching PIXCVecRiver file found in the rivertile directory.
    """
    cycle_tile_re = r'_\d{3}_\d{3}_\d{3}[LR]_'
    m = re.search(cycle_tile_re, rivertile)
    if not m:
        raise ValueError(
            f"Could not extract cycle/pass/tile info from filename: {rivertile}")
    cycle_pass_tile = m.group()
    rivertile_dir = os.path.dirname(os.path.abspath(rivertile))
    pixcvec_tag = 'SWOT_L2_HR_PIXCVecRiver' + cycle_pass_tile
    pixcvecriver = find_file(rivertile_dir, pixcvec_tag)
    if pixcvecriver is None:
        raise FileNotFoundError(f"No PIXCVecRiver file found for {rivertile}")
    return pixcvecriver


def get_pixc_from_rivertile(rivertile, river_run_id):
    """
    Extract cycle/pass/tile info from the rivertile file name and return
    the associated PIXC file found using the appropriate parent directory.
    """
    cycle_tile_re = r'_\d{3}_\d{3}_\d{3}[LR]_'
    m = re.search(cycle_tile_re, rivertile)
    if not m:
        raise ValueError(f"Could not extract cycle/pass/tile info from filename: {rivertile}")
    cycle_pass_tile = m.group()
    rivertile_dir = os.path.dirname(os.path.abspath(rivertile))
    parent_dir = os.path.dirname(rivertile_dir)
    pixc_tag = 'SWOT_L2_HR_PIXC' + cycle_pass_tile
    if "local" in river_run_id:
        pixc_dir = os.path.join(parent_dir, 'pixc')
    else:
        pixc_dir = os.path.dirname(parent_dir)
    pixc = find_file(pixc_dir, pixc_tag)
    if pixc is None:
        print(f"Warning: No PIXC file found for {rivertile}")
    return pixc


def load_pkl_dataframes(pkl_dir):
    """Load pkl and CSV dataframes from the specified directory."""
    field_dataframes = {}

    pt_node_file = os.path.join(pkl_dir, 'dataframe', 'matched_pt_node_df.pkl')
    pt_reach_file = os.path.join(pkl_dir, 'dataframe', 'pt_reach_wse_df.pkl')
    pt_reach_matched_wse = os.path.join(pkl_dir, 'dataframe',
                                        'matched_pt_reach_wse_df.pkl')
    pt_matched_slope = os.path.join(pkl_dir, 'dataframe',
                                    'matched_pt_reach_slope_df.pkl')
    drift_matched_nodes = os.path.join(pkl_dir, 'dataframe',
                                       'matched_drift_node_df.pkl')
    drift_matched_reaches = os.path.join(pkl_dir, 'dataframe',
                                         'matched_drift_reach_df.pkl')
    error_dataframe = os.path.join(pkl_dir, 'stats', 'ALL',
                                   'pt_matched_reach_error_table_ALL.csv')
    coarse_nodes = os.path.join(pkl_dir, 'stats', 'ALL',
                                'coarse_matched_node_error_table_ALL.csv')
    coarse_reaches = os.path.join(pkl_dir, 'stats', 'ALL',
                                  'offline_coarse_matched_reach_error_table_ALL.csv')

    field_dataframes = {'pt_node': pt_node_file, 'pt_reach': pt_reach_file,
        'pt_match_wse': pt_reach_matched_wse,
        'pt_match_slope': pt_matched_slope,
        'drift_node_match': drift_matched_nodes,
        'drift_reach_match': drift_matched_reaches, }
    for key in field_dataframes.keys():
        if os.path.isfile(field_dataframes[key]):
            field_dataframes[key] = pd.read_pickle(field_dataframes[key])
        else:
            print(
                f"pkl input {field_dataframes[key]} does not exist; check filenames!")
            field_dataframes[key] = None

    if os.path.isfile(error_dataframe):
        field_dataframes['pt_error'] = pd.read_csv(error_dataframe)
    else:
        alt_error = os.path.join(pkl_dir, 'stats', 'node', 'ALL',
                                 'pt_matched_node_error_table_ALL.csv')
        if os.path.isfile(alt_error):
            field_dataframes['pt_error'] = pd.read_csv(alt_error)

    if os.path.isfile(coarse_nodes):
        field_dataframes['coarse_node'] = pd.read_csv(coarse_nodes)
    else:
        alt_coarse = os.path.join(pkl_dir, 'dataframe', 'drift_node_df.pkl')
        if os.path.isfile(alt_coarse):
            field_dataframes['coarse_node'] = pd.read_csv(alt_coarse)
        else:
            field_dataframes['coarse_node'] = None

    if os.path.isfile(coarse_reaches):
        field_dataframes['coarse_reach'] = pd.read_csv(coarse_reaches)
    else:
        field_dataframes['coarse_reach'] = None

    return field_dataframes


def find_file(directory, pattern, recursive=False):
    """Finds the first file matching a pattern in a given directory."""
    if recursive:
        for root, _, files in os.walk(directory):
            for file in files:
                if file.startswith(pattern):
                    return os.path.join(root, file)
    else:
        files = [f for f in os.listdir(directory) if f.startswith(pattern)]
        if files:
            return os.path.join(directory, files[0])
    return None


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

def get_out_dir(out_dir, river_run_id, pixc_run_id, tag_output):
    if out_dir is not None:
        if river_run_id is not None:
            last_dir = river_run_id + tag_output
            if pixc_run_id is not None:
                this_out_dir = f'{out_dir}/{pixc_run_id}/' \
                    f'{last_dir}'
            else:
                this_out_dir = f'{out_dir}/{last_dir}/'
        else:
            this_out_dir = f'{out_dir}'
    else:
        this_out_dir = None
    return this_out_dir

def get_out_basename(river_file, reach_id, river_run_id, date):
    dir_parts = river_file.split('/')[-1]
    file_parts = dir_parts.split('_')
    if ('Reach' in file_parts) or ('Node' in file_parts):
        # input was a shapefile
        out_basename = date+'_'+file_parts[
            5] + '_' + file_parts[6] + '_' + file_parts[
            7] + '_' + str(reach_id) + '_' + river_run_id
    else:
        # input was a netcdf
        out_basename = date+'_'+file_parts[4] + '_' + file_parts[5] + '_' + \
            file_parts[6] + '_' + str(reach_id) + '_' \
            + river_run_id
    #
    out_basename = out_basename.replace(':','').replace(
                    '.','p').replace('-','_').replace(' ','_')
    return out_basename 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rivertile_dir', help='rivertile directory')
    parser.add_argument('--out_dir', help='output dir', default=None)
    parser.add_argument('--rivertile_csv',
                        help='CSV file with rivertile filenames',
                        default=None, type=str)
    parser.add_argument('--mt_basedir',
                        help='multitemporal stats base directory',
                        default=None, type=str)
    parser.add_argument('--slc_run_id', help='slc_run_id', default=None)
    parser.add_argument('--pixc_run_id', help='pixc_run_id', default=None)
    parser.add_argument('--river_run_id', help='river_run_id', default=None)
    parser.add_argument('--mt_flavor',help='multitemproal flavor',
                        default='v1', type=str)
    parser.add_argument('-p', '--passes', help='list of passes', nargs='+',
                        type=str.lower, default=None)
    parser.add_argument('-c', '--cycles', help='list of cycles', nargs='+',
                        type=str.lower, default=None)
    parser.add_argument('-t', '--tiles', help='list of tiles, e.g. 038R',
                        nargs='+', type=str.upper, default=None)
    parser.add_argument('-r', '--reaches', help='list of specific reaches',
                        nargs='+', type=int, default=None)
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
    print('REACHES: ', args.reaches)
    if args.rivertile_csv is not None:
        rivertiles, pixcvecs, pixcs, reach_ids = get_input_files_from_csv(
            args.rivertile_csv, args.pixc_run_id, args.river_run_id)
    else:
        rivertiles, pixcvecs, pixcs = get_input_files(
            args.rivertile_dir,
            args.pixc_run_id,
            args.river_run_id,
            args.slc_run_id,
            args.cycles,
            args.passes,
            args.tiles
        )
    # Process PKL/CSV dataframes if provided.
    if args.pkl_input_dir is not None:
        field_dataframes = load_pkl_dataframes(args.pkl_input_dir)
    else:
        field_dataframes = None
    for rivertile, pixcvec, pixc, this_reach_id in zip(
            rivertiles, pixcvecs, pixcs, reach_ids):
        print('\nloading river file:',rivertile,'\n')
        # check files up front before loading anything
        this_out_dir = get_out_dir(
            out_dir, args.river_run_id, args.pixc_run_id, args.tag_output)
        # create the output dir if needed
        if this_out_dir is not None:
            if not os.path.isdir(this_out_dir):
                os.umask(0)
                os.makedirs(this_out_dir, 0o777)
                os.makedirs(this_out_dir + '/truth/', 0o777)
                os.makedirs(this_out_dir + '/no_truth/', 0o777)
        # now get the outfile basename
        river_file = rivertile
        if isinstance(rivertile, tuple):
            river_file=rivertile[0]
        if this_reach_id is None:
            this_reach_id0 = '*'
        else:
            this_reach_id0 = this_reach_id
        out_basename = get_out_basename(
            river_file, this_reach_id0, args.river_run_id, date='*')
        # now skipcases that already have plots generated
        # unless commanded to overwrite
        if (~args.overwrite):
            glob_str = os.path.join(this_out_dir,'*truth', out_basename+'*.png')
            lst = glob.glob(glob_str)
            #breakpoint()
            if len(lst)>=2:
                print('output pngs already exist, not rerunning...')
                continue
        # now read the input data
        try:
            file_extension = os.path.splitext(rivertile)[1].lower()
            if file_extension == '.nc':
                reach_ids0, reach_wse, reach_width, river_names, \
                river_df = read_netcdf(rivertile)
                #river_file = rivertile
            elif file_extension=='.shp':
                # handle inputting only the node file assuming reach one is in
                # same place with identical filenames except the Reach/Node tag
                # assume it is RiverSP shape-file
                import SWOTRiver.products
                # assume it is RiverSP shape-file
                path, fname = os.path.split(rivertile)
                node_file = fname
                reach_file = node_file.replace('Node', 'Reach')
                if 'Reach' in node_file:
                    reach_file = fname
                    node_file = reach_file.replace('Reach', 'Node')
                # stick back on the path
                node_file = os.path.join(path, node_file)
                reach_file = os.path.join(path, reach_file)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    river_df = SWOTRiver.products.rivertile.L2HRRiverTile.from_shapes(
                        node_file, reach_file)
                # make area_total a masked array with fill-values filled
                river_df.nodes['area_total'] = np.ma.masked_array(
                    river_df.nodes['area_total'])
                #river_file = rivertile
                reach_ids0 = np.unique(river_df.reaches['reach_id'])
        except TypeError:
            if isinstance(rivertile, tuple):
                # input a list of tuples containing reach & node shapefiles
                reach_ids0, reach_wse, reach_width, river_names, \
                river_df = read_shapefile(rivertile)
                #river_file = rivertile[0]
            else:
                print('Unsupported file format')
        if args.calval:
            print('Calval arg set; processing calval rivers...')
            if set(river_names) & set(CALVAL_RIVERS):
                print(
                    'Tile', rivertile, 'has a calval river. Continuing...'
                )
            else:
                print('Tile', rivertile, 'does not have a calval river.')
                continue
        # check if the input is empty here
        if len(river_df.nodes['river_name'])==0:
            continue
        if len(river_df.reaches['river_name'])==0:
            continue
        if this_reach_id is not None:
            this_reach_ids = [this_reach_id,]
        else:
            this_reach_ids = reach_ids0
        # limit to desired reach list
        if args.reaches is not None:
            this_reach_ids = list(
                set(this_reach_ids).intersection(set(args.reaches)))
        # get date for filenames
        date = river_df.reaches.time_granule_start.split('T')[0]
        for reach_id in this_reach_ids:
            # load multitemporal stats
            if args.mt_basedir is not None:
                mt_wse, mt_width = load_mt_stats_files(
                    args.mt_basedir, reach_id, args.mt_flavor)
            #
            out_basename = get_out_basename(
                river_file, reach_id, args.river_run_id, date)
            if args.pkl_input_dir is not None:
                plot_reach.make_plots(river_file,
                    river_df, field_dataframes, pixcvec, pixc,
                    truth_pixcvec, truth_pixc, reach_id,
                    reach_error, nodes, pixc_truth, out_dir=this_out_dir,
                    out_basename=out_basename, overwrite=args.overwrite,
                    mt_wse=mt_wse, mt_width=mt_width
                )
            else:
                plot_reach.make_plots(rivertile,
                    river_df, truth, pixcvec, pixc, truth_pixcvec,
                    truth_pixc, reach_id, reach_error, nodes,
                    pixc_truth, out_dir=this_out_dir,
                    out_basename=out_basename,
                    overwrite=args.overwrite,
                    mt_wse=mt_wse, mt_width=mt_width
                )


if __name__ == "__main__":
    main()
