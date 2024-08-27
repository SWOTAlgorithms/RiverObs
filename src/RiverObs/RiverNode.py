"""
A river node holds all of the measurements associated with a node in a centerline.
It returns various data characteristics when queried.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import SWOTWater.aggregate as aggregate
from SWOTWater.constants import AGG_CLASSES
from SWOTRiver.errors import RiverObsUseageException

class RiverNode:
    """
    A river node holds all of the measurements associated with a node in a
    centerline.  It returns various data characteristics when queried.

    Initialize with the outputs of Centerline plus height and potentially
    other keywords that are turned into class variables.

    Parameters
    ----------

    index : int
        index in the center line corresponding to this node
    d : array_like
        distance from the node to each of the data points
    x : array_like
        x coordinate of each measurement associated with the node
    y : array_like
        y coordinate of each measurement associated with the node
    s : array_like
        along-track coordinate (relative to the node center) for each point
    n : array_like
        across-track (normal) coordinate (relative to the node center) for
        each point
    h_flg: array_like
        flag indicating which pixels to use for height stats
    ds : float
        along-track dimension for this node. Defaults to 1. Needs to be set
        correctly for width_area to work.

    Notes
    -----

    To edit the data, it is useful to sort it by order in one of the variables.
    To sort the data, set sort=True and set sort_variable one of
    ['d','s','n','h'].  If the data are not sorted initially, they can be
    sorted later, as desired by calling RiverNode.sort.
    """

    # These variables are sorted simultaneously when sort is called. Append
    # to this list after class instantiation if you want to add an additional
    # sort variable

    sort_vars = ['d', 'x', 'y', 's', 'n']

    def __init__(self, index, d, x, y, s, n, ds=1, **kwds):
        """
        """

        self.index = index
        self.ds = ds
        self.d = np.asarray(d)
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.s = np.asarray(s)
        self.n = np.asarray(n)

        for key, value in kwds.items():
            setattr(self, key, value)

        self.ndata = len(self.x)
        self.good = np.ones(self.ndata, dtype=bool)
        self.sorted = False
        self.height_std_pix = None

    def add_obs(self, obs_name, obs, sort=True):
        """
        Add a new observations, and, if desired, sort them in accordance to
        the current sort index.

        obs_name: name of the observation, will be added as self.obs_name
        obs: iterable with observations of length self.ndata
        """

        if len(obs) != self.ndata:
            raise RiverObsUseageException(
                'length of observations not consistent with number of node points'
            )

        self.sort_vars.append(obs_name)

        if sort and self.sorted:
            setattr(self, obs_name, obs[self.sort_index])
        else:
            setattr(self, obs_name, obs)

    def count(self, *pars, **kwds):
        """Return the number of points in the node."""
        return self.ndata

    def count_good(self, goodvar):
        """Return the number of good points in the node."""
        tmp = np.zeros(self.ndata)
        #exec('tmp[self.%s]=1'%goodvar)
        tmp[getattr(self, goodvar)] = 1
        return sum(tmp)

    def value(self, var):
        """
        Return self.value. Useful when getting stats for all nodes and a
        scalar variable is desired.
        """
        return getattr(self, var)

    def mean(self, var, goodvar='good'):
        """
        Return mean of a variable.
        """
        good = getattr(self, goodvar)
        mean = np.mean(getattr(self, var)[good])
        return mean

    def height_weighted_mean(self, var, goodvar='good', method='weight'):
        """
        Return weighted mean of a variable, using same weights as pixel heights
        """
        good = getattr(self, goodvar)
        if self.height_std_pix is None:
            getattr(self, 'pix_height_std')()

        var_out, weight_norm = aggregate.height_only(getattr(self, var), good,
                                                     self.height_std_pix,
                                                     method=method)

        return var_out

    def median(self, var, goodvar='good'):
        """Return median of a variable"""
        var = getattr(self, var)[getattr(self, goodvar)]
        if len(var) > 0:
            var_med = np.ma.median(var)
        else:
            var_med = np.nan
        return var_med

    def sincos_median(self, var, goodvar='good'):
        """Return atan2(median(sin(variable)), median(cos(variable)))*180/pi"""
        var = getattr(self, var)[getattr(self, goodvar)]
        if len(var) > 0:
            sincos_med = np.rad2deg(np.ma.arctan2(
                np.ma.median(np.ma.sin(np.deg2rad(var))),
                np.ma.median(np.ma.cos(np.deg2rad(var)))))
        else:
            sincos_med = np.nan
        return sincos_med

    def std(self, var, goodvar='good'):
        """Return std of a variable"""
        return np.std(getattr(self, var)[getattr(self, goodvar)])

    def stderr(self, var, goodvar='good'):
        """Return standard error of a variable"""
        good = getattr(self, goodvar)
        scale = 1. / np.sqrt(np.sum(good).astype(np.float))
        stderr = scale * np.std(getattr(self, var)[good])
        return stderr

    def min(self, var, goodvar='good'):
        """Return min of a variable"""
        good = getattr(self, goodvar)
        try:
            vmin = np.min(getattr(self, var)[good])
        except ValueError:
            vmax = np.nan
        return vmin

    def max(self, var, goodvar='good'):
        """Return max of a variable"""
        good = getattr(self, goodvar)
        try:
            vmax = np.max(getattr(self, var)[good])
        except ValueError:
            vmax = np.nan
        return vmax

    def ptp(self, var, goodvar='good'):
        """Return peak to peak variation of a variable"""
        good = getattr(self, goodvar)
        ptp = np.ptp(getattr(self, var)[good])
        return ptp

    def percentile(self, var, q, goodvar='good'):
        """Return qth percentile of a variable"""
        good = getattr(self, goodvar)
        try:
            percentile = np.percentile(getattr(self, var)[good], q)
        except IndexError:
            return np.nan
        return percentile

    def sum(self, var, goodvar='good'):
        """Return the sum of all variable values (e.g., for area)."""
        good = getattr(self, goodvar)
        if good.any():
            Sum = np.sum(getattr(self, var)[good])
            return Sum
        else:
            return 0

    def cdf(self, var, goodvar='good'):
        """Get the cdf for a variable."""
        good = getattr(self, goodvar)
        ngood = np.sum(good.astype(np.int32))
        cdf = np.cumsum(np.ones(ngood, dtype=np.float64))
        cdf /= cdf[-1]
        x = np.sort(getattr(self, var)[good])
        return x, cdf

    def bitwise_or(self, var, goodvar='good'):
        """Bitwise or of inputs"""
        good = getattr(self, goodvar)
        return np.bitwise_or.reduce(getattr(self, var)[good])

    def bitwise_and(self, var, goodvar='good'):
        """Bitwise and of inputs"""
        good = getattr(self, goodvar)
        return np.bitwise_and.reduce(getattr(self, var)[good])

    def flag_extent(self, var_name, var_min, var_max, goodvar='good'):
        """
        Update the good flag to places where the variable var_name
        is within a range of values.
        """
        good_init = getattr(self, goodvar)
        x = getattr(self, var_name)
        good = (x >= var_min) & (x <= var_max)
        self.good = good_init & good

    def add_sort_variables(self, vars):
        """
        Add sort variables to the sort list.
        vars is a string variable name, or an iterable of variable names.
        """
        if type(vars) == str:
            self.sort_vars.append(vars)
        else:
            for var in vars:
                self.sort_vars.append(var)

    def sort(self, sort_variable='n'):
        """
        Sort the data according to the variable sort_variable.
        self.sort_variable must exist.
        All of the variables in self.sort_vars are sorted in the same way.
        """
        # Get the index for sorting
        self.sort_index = np.argsort(getattr(self, sort_variable))

        # Sort all of the desired variables
        for var in self.sort_vars:
            setattr(self, var, getattr(self, var)[self.sort_index])

        # Make sure the good flag is also sorted
        self.good = self.good[self.sort_index]

        # This flag is checked when calling
        self.sorted = True

    def trim(self, fraction, mode='both'):
        """
        Calculate a index flag that trims the data, after sorting.

        fraction: 0 < f < 1. Fraction of the data to remove.
        mode is 'both', 'low', 'high' for which tails of the distribution
        need to be trimmed.
        """
        if not self.sorted:
            raise RiverObsUseageException('Run sort before calling trim')

        if mode == 'both':
            fraction = fraction / 2.

        ntrim = int(self.ndata * fraction + 0.5)

        self.trim_index = np.ones(self.ndata, dtype=bool)

        if mode == 'both':
            self.trim_index[:ntrim] = False
            self.trim_index[:-ntrim] = False
        if mode == 'low':
            self.trim_index[:ntrim] = False
        else:
            self.trim_index[:ntrim] = False

        self.good = self.good & self.trim_index

        self.trimmed = True

    def width_ptp(self, *pars, **kwds):
        """
        Return the river width at this node, estimated from the max/min
        extent of the normal coordinate.

        Call signature made compatible with other stat calls, but requires
        no inputs.
        """
        return self.max('n') - self.min('n')

    def width_std(self, *pars, **kwds):
        """
        Return the river width at this node, estimated as sqrt(12)*std(n).
        This would be appropriate for uniformly distributed point measurements.

        Call signature made compatible with other stat calls, but requires
        no inputs.
        """
        nstd = self.std('n')
        width_std = np.sqrt(12.) * nstd
        return width_std

    def width_area(self, area_var='pixel_area'):
        """
        Return the river width at this node, estimated as Area/ds.
        This would be appropriate for pixels representing areas.

        area_var: str
            name of the variable associated with the pixel area
            (default 'pixel_area')

        Call signature made compatible with other stat calls, but requires
        no inputs.
        """
        area = self.sum(area_var)
        width_area = area / self.ds
        return width_area

    def sig0_with_uncert(self, goodvar='good'):
        """
        Returns the aggregate sigma0 and sigma0_unc
        """
        good = getattr(self, goodvar)
        return aggregate.sig0_with_uncerts(self.sig0, good, self.sig0_uncert)

    def pix_height_std(self):
        """
        Returns the pixel-level height uncertainties
        """
        height_std_pix = np.abs(self.phase_noise_std * self.dh_dphi)
        # set bad pix height std to high number to de-weight
        # instead of giving infs/nans
        bad_num = 1.0e5
        height_std_pix[height_std_pix <= 0] = bad_num
        height_std_pix[np.isinf(height_std_pix)] = bad_num
        height_std_pix[np.isnan(height_std_pix)] = bad_num
        self.height_std_pix = height_std_pix

    def height_with_uncert(self, goodvar='good', method='weight'):
        """
        Return the aggregate height with corresponding uncertainty 
        """
        good = getattr(self, goodvar)
        if self.height_std_pix is None:
            getattr(self, 'pix_height_std')()

        # call the general function
        return aggregate.height_with_uncerts(
            self.h_noise,  good, self.num_rare_looks, self.num_med_looks,
            self.ifgram, self.power1, self.power2, self.looks_to_efflooks,
            self.dh_dphi, self.dlat_dphi, self.dlon_dphi, self.height_std_pix,
            method=method)

    def area_with_uncert(self, goodvar='good', method='composite'):
        """
        Return the aggregate width_area with corresponding uncertainty 
        """
        # compute the pixel assignment error?
        # call the general function, TODO

        # should normally just use all the data 
        # (not just the use_heights pixels), but could use goodvar 
        # to filter out outliers
        good = getattr(self, goodvar)

        # assumes we dont give land pixels in the class_list if method=simple
        # and that we set the use_fractional_inundation to true for land edge
        # pixels if method=water_fraction or method=composite

        # call the external function to aggregate areas and uncertainties
        # using dark water classes.
        area, area_unc, area_pcnt_uncert = aggregate.area_with_uncert(
            self.pixel_area, self.water_frac, self.water_frac_uncert,
            self.darea_dheight, self.klass, self.false_detection_rate,
            self.missed_detection_rate, good,
            Pca=0.9, Pw=0.5, Ptf=0.5, ref_dem_std=10,
            interior_water_klasses=AGG_CLASSES['interior_water_klasses'],
            water_edge_klasses=AGG_CLASSES['water_edge_klasses'],
            land_edge_klasses=AGG_CLASSES['land_edge_klasses'],
            dark_water_klasses=AGG_CLASSES['dark_water_klasses'],
            method=method)

        width_area = area/self.ds
        width_area_unc = area_unc/self.ds

        # Recompute areas and uncertainties without dark water
        area_det, area_det_unc, area_det_pcnt_uncert = \
            aggregate.area_with_uncert(
                self.pixel_area, self.water_frac, self.water_frac_uncert,
                self.darea_dheight, self.klass, self.false_detection_rate,
                self.missed_detection_rate, good,
                Pca=0.9, Pw=0.5, Ptf=0.5, ref_dem_std=10,
                interior_water_klasses=AGG_CLASSES['interior_water_klasses'],
                water_edge_klasses=AGG_CLASSES['water_edge_klasses'],
                land_edge_klasses=AGG_CLASSES['land_edge_klasses'],
                dark_water_klasses=[],
                method=method)

        return (area, width_area, area_unc, width_area_unc, area_det,
                area_det_unc)
