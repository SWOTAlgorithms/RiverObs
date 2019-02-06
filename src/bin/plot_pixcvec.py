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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rivertile', help='rivertile.nc')
    parser.add_argument('pixcvec', help='pixel cloud vector.nc')
    args = parser.parse_args()

    rivertile = L2HRRiverTile.from_ncfile(args.rivertile)
    pixc_vector = L2PIXCVector.from_ncfile(args.pixcvec)

    scatter_colors = ['k', 'r', 'b', 'g']
    figure, axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    for ii, node_id in enumerate(rivertile.nodes.node_id):

        mask = pixc_vector.node_index == node_id
        this_color = scatter_colors[ii%len(scatter_colors)]
        axis.scatter(
            pixc_vector.longitude_vectorproc[mask],
            pixc_vector.latitude_vectorproc[mask],
            s=50, c=this_color, edgecolor='none')

    is_bad_width = rivertile.nodes.node_q & 0x1 > 0

    axis.scatter(
        rivertile.nodes.longitude[is_bad_width],
        rivertile.nodes.latitude[is_bad_width],
        marker='s', c='m')

    axis.scatter(
        rivertile.nodes.longitude[~is_bad_width],
        rivertile.nodes.latitude[~is_bad_width],
        marker='x', c='m')

    axis.set_xlabel('longitude')
    axis.set_ylabel('latitude')
    plt.show()

if __name__ == "__main__":
    main()
