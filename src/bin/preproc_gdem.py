#!/usr/bin/env python
"""
Pre-processes a GDEM for use with RiverObs validation studies

Author (s): Alex Fore
"""
import shutil
import netCDF4
import argparse
import logging
import scipy.stats
import scipy.ndimage
import numpy as np
import scipy.ndimage.morphology

import RiverObs.ReachDatabase

LOGGER = logging.getLogger(__name__)

def wrap_label_for_plots(arr_in, wrap_value=20):
    # set 0 values to nan and wrap every
    arr_out=np.zeros(np.shape(arr_in)) + np.nan
    arr_out[arr_in>0] = arr_in[arr_in>0]
    return np.mod(arr_out, wrap_value)

def erosion_segmentation(ltype_in, erosion_iter, plot=False):
    LOGGER.info('erosion_segmentation')
    # First erode the water mask
    ltypeb = scipy.ndimage.morphology.binary_erosion(
        ltype_in,iterations = erosion_iter)
    ltype = np.zeros_like(ltype_in)
    ltype[ltypeb] = 1

    # Now segment the image
    type_label, num_labels = scipy.ndimage.label(ltype)
    
    # Assign closest label to orphaned pixels 
    # (i.e., pixels that were originally eroded out of mask)
    # first diate out to get back the edges and assign to nearest feature
    type_label_tmp = type_label.copy()
    for i in range(erosion_iter):
        #print("grey_dilation iter:",i)
        LOGGER.debug(
                "grey dilation iter: {} of {}".format(i,erosion_iter))
        tmp1 = scipy.ndimage.morphology.grey_dilation(
            type_label_tmp,size=5)
        msk = np.logical_and(np.logical_and(type_label_tmp==0, tmp1>0), ltype_in==1)
        type_label_tmp[msk] = tmp1[msk]
    type_label_tmp = type_label_tmp * ltype_in
    # Now handle remaining orphaned regions 
    # (e.g., skinny things that completely eroded)
    orphan_mask = np.zeros_like(ltype_in)
    orphan_mask[np.logical_and(ltype_in==1, type_label_tmp<=0)] = 1
    # Dilate the orphan mask to connect regions that touch
    orphan_dilated = scipy.ndimage.morphology.binary_dilation(
        orphan_mask,iterations = 1)
    # segment the orphans
    orphan_label, num_orphan_labels = scipy.ndimage.label(orphan_dilated)
    # Combine labels with orphans as highest labels
    combo_label = type_label_tmp.copy()
    combo_label[orphan_label>0] = orphan_label[orphan_label>0] + np.max(type_label_tmp)
    combo_label[ltype_in==0] = 0
    combo_label0 = combo_label.copy()

    # See if any orphans are actually touching other features
    # if so, reassign label to those features in the combo_label
    objects = scipy.ndimage.find_objects(orphan_label)
    labels = np.unique(orphan_label)
    for i, slyce in enumerate(objects):
        label_intersect = type_label_tmp[slyce]*orphan_dilated[slyce]
        if np.sum(label_intersect)>0:
            msk = np.logical_and(
                orphan_label[slyce]==labels[i+1], ltype_in[slyce]==1)
            combo_label[slyce][msk] = scipy.stats.mode(
                label_intersect[label_intersect>0])[0]
    if (False):
        import matplotlib.pyplot as plt
        # for setting the zoom
        a0=7700
        a1=8200
        c0=2700
        c1=3100
        
        plt.figure()
        plt.imshow(wrap_label_for_plots(type_label_tmp),cmap='jet')
        plt.colorbar()
        plt.title('type_label_tmp')
        plt.ylim([a0,a1])
        plt.xlim(c0,c1)

        plt.figure()
        plt.imshow(wrap_label_for_plots(type_label),cmap='jet')
        plt.colorbar()
        plt.title('type_label')
        plt.ylim([a0,a1])
        plt.xlim(c0,c1)

        plt.figure()
        plt.imshow(ltype_in,cmap='jet')
        plt.colorbar()
        plt.title('ltype_in')
        plt.ylim([a0,a1])
        plt.xlim(c0,c1)

        
        plt.figure()
        plt.imshow(wrap_label_for_plots(combo_label),cmap='jet')
        plt.colorbar()
        plt.title('combo_label')
        plt.ylim([a0,a1])
        plt.xlim(c0,c1)

       
        plt.figure()
        plt.imshow(wrap_label_for_plots(combo_label0),cmap='jet')
        plt.colorbar()
        plt.title('combo_label0')
        plt.ylim([a0,a1])
        plt.xlim(c0,c1)

       
        plt.figure()
        plt.imshow(wrap_label_for_plots(orphan_label),cmap='jet')
        plt.colorbar()
        plt.title('orphan_label')
        plt.ylim([a0,a1])
        plt.xlim(c0,c1)
        plt.show()

    return combo_label

def select_river_labels(type, type_label, gdem_x, gdem_y, reaches):
    """
    Picks the labels that are nearest to the reach lons/lats.  For each reach
    finds the closest label.
    """
    LOGGER.info('select_river_label')

    min_compare = 3200 # 10m * 3.125m * 3200 == 2*200m*250m (2 nodes)
    max_compare = 100000

    uniq_labels = np.unique(type_label[type == 1])
    cnts, _ = np.histogram(type_label[type == 1], uniq_labels)

    idxsort = np.argsort(cnts)[::-1]
    uniq_labels = uniq_labels[idxsort]
    cnts = cnts[idxsort]

    labels = []
    for ireach, reach in enumerate(reaches):
        type_dist = 9999999*np.ones(uniq_labels.shape)
        for ilabel, uniq_label in enumerate(uniq_labels):

            # probably not it if only 1000 pixels in feature
            if cnts[ilabel] < min_compare:
                continue
            
            # get the pixels in the feature
            these_x = gdem_x[type_label == uniq_label]
            these_y = gdem_y[type_label == uniq_label]

            # subsample randomly if its big 
            if cnts[ilabel] > max_compare:
                these_x, these_y = np.random.permutation(
                    np.array([these_x, these_y]))[:, :max_compare]

            # compute distance between the pixels in the feature and the reach
            delta2 = (
                (these_x[:, np.newaxis] - reach.x)**2 +
                (these_y[:, np.newaxis] - reach.y)**2)

            # use the min distance as the feature distance
            min_d2 = delta2.min(axis=0)
            type_dist[ilabel] = np.mean(np.sqrt(min_d2))

            LOGGER.debug(
                "reach, label, dist, num: {} {} {} {} {}".format(
                    ireach, uniq_label, type_dist[ilabel], len(these_x), cnts[ilabel]))

        this_label = uniq_labels[type_dist.argmin()]
        print ('this_label: ',this_label)
        if this_label not in labels:
            labels.append(this_label)

    return labels

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_gdem_file', help='Input GDEM file')
    parser.add_argument('out_gdem_file', help='Output GDEM file')
    parser.add_argument('reachdb_path', help='reach DB path/file')
    parser.add_argument(
        '-l', '--log-level', type=str, default="warning",
        help="logging level, one of: debug info warning error")
    parser.add_argument('--plot', default=False, action='store_true')
    parser.add_argument('--erosion-iter', type=int, default=0,
                        help='erosion iterations')
    args = parser.parse_args()

    level = {'debug': logging.DEBUG, 'info': logging.INFO,
             'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=format)

    with netCDF4.Dataset(args.in_gdem_file, 'r') as ifp:
        in_type = ifp.variables['landtype'][:]
        lat = ifp.variables['latitude'][:]
        lon = ifp.variables['longitude'][:]

    type = np.zeros(in_type.shape)
    type[in_type == 1] = 1

    # Extract prior node locations, first make a LatLonRegion
    llbox = RiverObs.ReachDatabase.LatLonRegion(
        [lon.min(), lat.min(), lon.max(), lat.max()])

    # project GDEM coordinates
    gdem_x, gdem_y = llbox.proj(lon, lat)

    # Extract Reaches
    reaches = RiverObs.ReachDatabase.ReachExtractor(args.reachdb_path, llbox)

    # Optionally erode before segmentation
    if args.erosion_iter > 0:
        type_label = erosion_segmentation(type, args.erosion_iter, plot=args.plot)
    else:
        # Labeled in descending order of counts
        type_label, num_labels = scipy.ndimage.label(type)

    # Get land and river labels.
    land_label = 0
    river_labels = select_river_labels(
        type, type_label, gdem_x, gdem_y, reaches)

    # Water that is not in the largest water features
    water_not_main_label = np.logical_and(
        np.isin(type_label, river_labels, invert=True),
        type_label != land_label)

    out_type = type.copy()
    out_type[water_not_main_label] = 0

    if args.plot:
        import matplotlib.pyplot as plt

        reach_lon, reach_lat = np.array([]), np.array([])
        for item in reaches:
            reach_lon = np.append(reach_lon, item.lon)
            reach_lat = np.append(reach_lat, item.lat)

        figure, axis = plt.subplots()
        axis.plot(lon[type==1], lat[type==1], 'g.')
        axis.plot(lon[out_type==1], lat[out_type==1], 'k.')
        #axis.scatter(lon[type==1], lat[type==1],
        #             c=wrap_label_for_plots(type_label[type==1]),
        #             edgecolors='none')
        axis.plot(reach_lon, reach_lat, 'rs')
        axis.set_xlabel('longitude')
        axis.set_ylabel('latitude')

        figure, axis = plt.subplots()
        axis.imshow(type+out_type)
        axis.set_title('input+output landtype')

        figure, axis = plt.subplots()
        axis.imshow(wrap_label_for_plots(type_label))
        axis.set_title('segmentation label')

        plt.show()

    shutil.copy(args.in_gdem_file, args.out_gdem_file)
    with netCDF4.Dataset(args.out_gdem_file, 'a') as ofp:
        ifp.variables['landtype'][:] = out_type

if __name__ == "__main__":
    main()
