#!/usr/bin/env python
'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s):
'''
import warnings
import argparse
import os
import numpy as np
import matplotlib.axes
import matplotlib.pyplot as plt

import SWOTRiver.analysis.riverobs

FIGSIZE = (12, 8)
DPI = 200
CMAP = 'plasma'

class NodePlots():
    def __init__(self, truth, data, title=None, tofile=False):
        self.truth = truth
        self.data = data
        self.title = title
        self.tofile = tofile
        self.figures = []
        self.axes = []
        self.plot()

    def add_figure(self):
        figure, axes = plt.subplots(3, 1, figsize=FIGSIZE, dpi=DPI)
        self.figures.append(figure)
        self.axes.append(axes)
        return figure, axes

    def plot(self):
        reaches = np.unique(self.data['reach_id'])
        for reach in reaches:
            i_truth = (self.truth['reach_id'] == reach)
            i_data = (self.data['reach_id'] == reach)
            reach_id = self.data['reach_id'][i_data][0]
            self.plot_reach(reach, reach_id, i_truth, i_data)
        self.finalize()

    def finalize(self):
        reaches = np.unique(self.data['reach_id'])
        for reach, figure, reach_axes in zip(reaches, self.figures, self.axes):
            i_data = (self.data['reach_id'] == reach)
            reach_id = self.data['reach_id'][i_data][0]
            for axis in reach_axes:
                axis.grid()
            if self.title is not None:
                title = '{}-{:03d}-{:03d}'.format(self.title, reach, reach_id)
                reach_axes[0].set_title(title)
                if self.tofile:
                    figure.savefig(title + '.png')


class ParamPlots(NodePlots):
    def __init__(
            self, truth, data, fields=['area_total'], truth_field='area_total',
            percent=False, **kwargs):
        self.fields = fields
        self.truth_field = truth_field
        self.percent = percent
        super().__init__(truth, data, **kwargs)

    def plot_reach(self, reach, reach_id, i_truth, i_data):
        figure, axes = self.add_figure()
        axes[0].set_title('reach %s %s' % (reach, reach_id))
        self.plot_data(i_truth, i_data, axes[0])
        self.plot_data(i_truth, i_data, axes[1], style='-')
        self.plot_error(i_truth, i_data, axes[2], style='.-')
        # TODO: compute the bulk errors for each reach

    def plot_data(self, i_truth, i_data, axis, style='.'):
        lgnd = []
        for field in self.fields:
            if not field.endswith('_u'):
                lgnd.append('data '+field)
                axis.plot(
                    self.data['node_id'][i_data], self.data[field][i_data],
                    style, markersize=10)
                
        axis.plot(
            self.truth['node_id'][i_truth],
            self.truth[self.truth_field][i_truth], style, markersize=10,
            marker='+', color='g')
        axis.set_ylabel(self.truth_field)
        axis.legend(
            lgnd + ['truth ' + self.truth_field], loc='best')

    def plot_error(self, i_truth, i_data, axis, style='.'):
        for field in self.fields:
            if field.endswith('_u'):
                error = self.data[field][i_data]
                print(error)
                label = field
            else:
                Ldata = len(self.data[field][i_data])
                Ltruth = len(self.truth[field][i_truth])
                if Ldata == Ltruth:
                    error = (
                        self.data[field][i_data]
                        - self.truth[field][i_truth])
                else:
                    warnings.warn(
                        'Different number of nodes truth/data {:d} {:d}'.format(
                        Ltruth, Ldata))
                    return
                label = self.truth_field + ' error'
                if self.percent:
                    error = error / self.truth[self.truth_field][i_truth] * 100.0
                    label = self.truth_field + '% error'
            axis.plot(self.data['node_id'][i_data], error, label=label)
            
        # axis.plot(
        #     self.truth['node_id'][i_truth], self.truth['height'][i_truth],
        #     marker='+')
        axis.set_xlabel('node_id')
        if self.percent:
            axis.set_ylabel(self.truth_field + '% error, data-truth')
            axis.plot(
                [-9e9, 9e9], [15, 15], '--y', scalex=False,
                label='Req. w>100m')
            axis.plot(
                [-9e9, 9e9], [25, 25], ':y', scalex=False,
                label='Goal, w>100m')
            axis.plot(
                [-9e9, 9e9], [-15, -15], '--y', scalex=False)
            axis.plot(
                [-9e9, 9e9], [-25, -25], ':y', scalex=False)
            axis.legend()
        else:
            axis.set_ylabel(self.truth_field + ' error, data-truth')
    def finalize(self):
        for axes in self.axes:
            ylims = axes[2].get_ylim()
            edge = np.max(np.abs(ylims))
            axes[2].set_ylim(-1*edge, edge)
        super().finalize()


class HeightPlots(ParamPlots):
    def finalize(self):
        for axes in self.axes:
            axes[0].set_ylabel('height (m)')
            axes[1].set_ylabel('height (m)')
            axes[2].set_ylabel('height error (m), data-truth')
            axes[2].set_ylim(-0.5, 0.5)
            axes[2].set_yticks(np.linspace(-0.5, 0.5, 11))
            for i in [-1, 1]:
                labels = [
                    'BSM for $A>(250m)^2$', 'BSM for $A>1 km^2$',
                    'TSM for $A>1 km^2$']
                if i == 1:
                    labels = [None, None, None]
                axes[2].plot(
                    [-9e9, 9e9], [i*0.25, i*0.25], ':y', scalex=False,
                    label=labels[0])
                axes[2].plot(
                    [-9e9, 9e9], [i*0.10, i*0.10], '--y', scalex=False,
                    label=labels[1])
                axes[2].plot(
                    [-9e9, 9e9], [i*0.11, i*0.11], '--r', scalex=False,
                    label=labels[2])
            axes[2].legend()
        super().finalize()


def plot_locations(
        truth, data, color_field='reach_id', title=None, tofile=False):
    figure, axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    plot = axis.scatter(
        data['longitude'], data['latitude'], s=50, c=data[color_field],
        edgecolor='none')
    colorbar = plt.colorbar(plot)
    colorbar.set_label(color_field)
    axis.scatter(truth['longitude'], truth['latitude'], marker='+', c='k')
    axis.grid()
    axis.set_xlabel('longitude')
    axis.set_ylabel('latitude')
    axis.legend(['data node', 'truth node'])
    if title is not None:
        axis.set_title(title)
        if tofile:
            figure.savefig(title + '.png')

def plot_prior_locations(truth, data, color_field='reach_id'):
    figure, axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    plot = axis.scatter(
        data['x_prior'], data['y_prior'], s=50, c=data[color_field],
        edgecolor='none')
    colorbar = plt.colorbar(plot)
    colorbar.set_label(color_field)
    axis.scatter(truth['x_prior'], truth['y_prior'], marker='+', c='k')
    axis.set_xlabel('x_prior')
    axis.set_ylabel('y_prior')
    axis.legend(['data node', 'truth node'])
    axis.grid()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_rivertile', help='pixel cloud rivertile.nc')
    parser.add_argument('gdem_rivertile', help='GDEM rivertile.nc')
    parser.add_argument('-t', '--title', default='node')
    parser.add_argument('-p', '--print', action='store_true')
    args = parser.parse_args()

    truth, data = SWOTRiver.analysis.riverobs.load_rivertiles(
        args.gdem_rivertile, args.pixc_rivertile)
    SWOTRiver.analysis.riverobs.match_rivertiles(truth, data)
    HeightPlots(
        truth.nodes, data.nodes, ['height', 'height_u'], 'height',
        title=args.title+'-height', tofile=args.print)
    ParamPlots(
        truth.nodes, data.nodes, ['area_detct',], 'area_detct',
        percent=True, title=args.title+'-area', tofile=args.print)
    plot_locations(
        truth.nodes, data.nodes, title=args.title+'-locations',
        tofile=args.print)
    # plot_prior_locations(data.nodes, truth.nodes)
    if not args.print:
        plt.show()

if __name__ == "__main__":
    main()
