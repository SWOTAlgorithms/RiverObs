#!/usr/bin/env python
'''
Copyright (c) 2017-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author(s):

'''
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

import SWOTRiver.analysis.riverobs
import SWOTRiver.analysis.tabley


FIGSIZE = (12, 8)
DPI = 200
CMAP = 'plasma'


class ReachPlot():
    def __init__(self, truth, metrics, title=None, filename=None):
        self.truth = truth
        self.metrics = metrics
        self.title = title
        self.filename = filename
        self.figure, self.axis = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        self.plot()

    def plot(self, independent, dependent, color_key):
        scatter = self.axis.scatter(
            independent, dependent, c=self.truth.reaches[color_key],
            cmap=CMAP)
        colorbar = plt.colorbar(scatter)
        colorbar.set_label(color_key)
        self.plot_requirements()
        self.finalize()

    def finalize(self):
        self.axis.grid()
        if self.title is not None:
            self.axis.set_title(self.title)
        if self.filename is not None:
            self.figure.savefig(self.filename)


class AreaPlot(ReachPlot):
    def plot(self):
        true_area = np.sqrt(self.truth.reaches.area_detct)
        super().plot(true_area, self.metrics['area'], 'xtrk_dist')

    def plot_requirements(self):
        true_area = np.sqrt(self.truth.reaches.area_detct)
        buff = 100
        i = 1
        self.axis.plot([np.sqrt(50*10e3), 1000], [i*25, i*25], '--g')
        self.axis.plot([1000, np.sqrt(170*10e3)], [i*15, i*15], '--y')
        self.axis.plot(
            [np.sqrt(170*10e3), np.amax(true_area)+buff], [i*15, i*15], '--r')
        self.axis.set_xlim((0, np.amax(true_area)+buff))
        self.axis.legend(
            ['SG for $A>0.7 km^2$', 'BSM for $A>1 km^2$',
             'TSM for $A>1.3 km^2$', '$<10km$ reach','$\geq 10km$ reach'],
            loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15))
        i = -1
        self.axis.plot([np.sqrt(50*10e3), 1000], [i*25, i*25], '--g')
        self.axis.plot([1000, np.sqrt(170*10e3)], [i*15, i*15], '--y')
        self.axis.plot(
            [np.sqrt(170*10e3), np.amax(true_area)+buff], [i*15, i*15], '--r')

    def finalize(self):
        self.axis.set_ylabel('area % error')
        self.axis.set_xlabel('sqrt reach area (m)')
        super().finalize()


class HeightPlot(ReachPlot):
    def plot(self):
        true_area = np.sqrt(self.truth.reaches.area_detct)
        super().plot(true_area, self.metrics['height'], 'xtrk_dist')

    def plot_requirements(self):
        true_area = np.sqrt(self.truth.reaches.area_detct)
        buff = 100
        i = 1
        self.axis.plot([250, 1000], [i*25, i*25], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*10, i*10], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*11, i*11], '--r')
        self.axis.set_xlim((0,np.amax(true_area)+buff))
        self.axis.legend(
            ['BSM for $A>(250m)^2$', 'BSM for $A>1 km^2$',
             'TSM for $A>1 km^2$', '$<10km$ reach','$\geq 10km$ reach'],
            loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15))
        i = -1
        self.axis.plot([250, 1000], [i*25, i*25], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*10, i*10], '--y')
        self.axis.plot([1000, np.amax(true_area)+buff], [i*11, i*11], '--r')

    def finalize(self):
        self.axis.set_ylabel('height error (cm)')
        self.axis.set_xlabel('sqrt reach area (m)')
        super().finalize()


class SlopePlot(ReachPlot):
    def plot(self):
        true_width = self.truth.reaches.width
        super().plot(true_width, self.metrics['slope'], 'xtrk_dist')

    def plot_requirements(self):
        true_width = self.truth.reaches.width
        buff = 100
        i = 1
        self.axis.plot([100, np.amax(true_width)+buff], [i*1.7, i*1.7], '--y')
        self.axis.plot([100, np.amax(true_width)+buff], [i*3, i*3], '--r')
        self.axis.set_xlim((0, np.amax(true_width)+buff))
        self.axis.legend(
            ['BSM for $A>1 km^2$', 'TSM for $A>1 km^2$',
             '$<10km$ reach','$\geq 10km$ reach'],
            loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.1))
        i = -1
        self.axis.plot([100, np.amax(true_width)+buff], [i*1.7, i*1.7], '--y')
        self.axis.plot([100, np.amax(true_width)+buff], [i*3, i*3], '--r')

    def finalize(self):
        self.axis.set_ylabel('slope error (cm/km)')
        self.axis.set_xlabel('river width (m)')
        super().finalize()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pixc_rivertile', help='pixel cloud rivertile.nc')
    parser.add_argument('gdem_rivertile', help='GDEM rivertile.nc')
    parser.add_argument('-t', '--title')
    parser.add_argument('-p', '--print', action='store_true')
    args = parser.parse_args()

    metrics = None
    truth = None
    truth_tmp, data = SWOTRiver.analysis.riverobs.load_rivertiles(
        args.gdem_rivertile, args.pixc_rivertile)
    SWOTRiver.analysis.riverobs.match_reaches(truth_tmp, data)
    tmp = SWOTRiver.analysis.riverobs.get_metrics(
        truth_tmp.reaches, data.reaches)
    if metrics is None:
        metrics = tmp
        truth = truth_tmp
    else:
        for name, value in tmp.items():
            metrics[name] = np.append(metrics[name], value)
        truth = truth.append(truth_tmp)
    print(truth.reaches.area_total.shape)
    SWOTRiver.analysis.tabley.print_table(metrics)
    SWOTRiver.analysis.riverobs.print_errors(metrics)

    if args.print:
        filenames = ['reach-area.png', 'reach-height.png', 'reach-slope.png']
        if args.title is not None:
            filenames = [args.title + '-' + name for name in filenames]
    else:
        filenames = [None, None, None]
    AreaPlot(truth, metrics, args.title, filenames[0])
    HeightPlot(truth, metrics, args.title, filenames[1])
    SlopePlot(truth, metrics, args.title, filenames[2])
    if not args.print:
        print('show')
        plt.show()


if __name__ == '__main__':
    main()
