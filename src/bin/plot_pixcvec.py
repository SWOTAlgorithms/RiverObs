#!/usr/bin/env python
"""
Plots pixels from pixcvec and node locations
"""
import numpy as np
import argparse
import matplotlib.axes
import matplotlib.pyplot as plt

import RiverObs.ReachDatabase
from SWOTRiver.products.rivertile import L2HRRiverTile
from SWOTRiver.products.pixcvec import L2PIXCVector

FIGSIZE = (12, 8)
DPI = 200

class LLBOX(object):
    def __init__(self, bbox):
        self.bounding_box = bbox

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rivertile', help='rivertile.nc')
    parser.add_argument('pixcvec', help='pixel cloud vector.nc')
    args = parser.parse_args()

    rivertile = L2HRRiverTile.from_ncfile(args.rivertile)
    pixc_vector = L2PIXCVector.from_ncfile(args.pixcvec)


    reach_ids = rivertile.reaches.reach_id

    scatter_colors = ['k', 'r', 'b', 'g']

    figure, axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    for ii, node_id in enumerate(rivertile.nodes.node_id):

        this_color = scatter_colors[ii%len(scatter_colors)]

        mask = pixc_vector.node_index == node_id
        axis.scatter(
            pixc_vector.longitude_vectorproc[mask],
            pixc_vector.latitude_vectorproc[mask],
            s=50, c=this_color, edgecolor='none')

    axis.scatter(
        rivertile.nodes.longitude, rivertile.nodes.latitude,
        marker='+', c='k')

    axis.set_xlabel('longitude')
    axis.set_ylabel('latitude')
    plt.show()

#     figure, axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
#     plot = axis.scatter(
#         data['longitude'], data['latitude'], s=50, c=data[color_field],
#         edgecolor='none')

if __name__ == "__main__":
    main()
