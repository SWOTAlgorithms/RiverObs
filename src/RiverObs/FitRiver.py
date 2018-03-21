"""
A class to fit various river parameters over a set of reaches.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import statsmodels.api


class FitRiver:
    """A class to fit various river parameters over a set of reaches.

    Initialize with a RiverObs object where all of the relevant variables
    have been loaded into which all of the observables have been loaded into
    the river nodes.
    """

    def __init__(self, river_obs):
        self.river_obs = river_obs
        self.inputs_computed = False

    def get_linear_fit_inputs(self,
                              smin,
                              smax,
                              fit_var,
                              mean_stat='mean',
                              err_stat='stderr'):
        """Get the design matrix, the observations, and the fitting weights
        for all the nodes in the interval [smin,smax].

        Parameters
        ----------

        fit_var : str
            the variable to fit
        mean_stat : str
            how the data are averaged in the node; e.g., 'mean', 'median', etc.
        err_stat: str
            formal errors for weighting the node; e.g., 'stderr', 'std', etc.
        """

        s = np.asarray(self.river_obs.get_node_stat(mean_stat, 's'))

        good = (s >= smin) & (s <= smax)

        s = s[good]
        nsample = len(s)

        # To get the values around the center value, subtract the mean

        x = s - np.mean(s)

        # Observations

        y = np.asarray(self.river_obs.get_node_stat(mean_stat, fit_var))[good]

        # Weights
        tmp = np.asarray(self.river_obs.get_node_stat(err_stat, fit_var))[good]
        w = 1.0 / tmp
        w[tmp <= 0] = 1.0
        #w = 1./np.asarray( self.river_obs.get_node_stat(err_stat,fit_var) )[good]

        # Fitting matrix for linear fit

        X = np.c_[x, np.ones(nsample, dtype=x.dtype)]

        self.inputs_computed = True

        return s, y, X, w

    def fit_linear(self,
                   smin,
                   smax,
                   fit_var,
                   fit='OLS',
                   mean_stat='mean',
                   err_stat='stderr',
                   load_inputs=True):
        """
        Fit the variable in the interval [smin,smax] using various fitting
        methods.

        Parameters
        ----------

        fit : str
            'OLS': ordinary least squares, 'WLS': weighted least squares,
            'RLM': robust linear model (in which case norm is used).
        fit_var : str
            the variable to fit
        mean_stat : str
            how the data are averaged in the node; e.g., 'mean', 'median', etc.
        err_stat: str
            formal errors for weighting the node; e.g., 'stderr', 'std', etc.
        load_inputs : bool
            if True, calculates the fitting inputs. If False, assumes that
            they are available from a previous fit. Useful if using multiple
            fit methods using the same data.
        """

        # Get the fit inputs, if desired

        if load_inputs or (not self.inputs_computed):
            self.s, self.y, self.X, self.w = self.get_linear_fit_inputs(
                smin, smax, fit_var, mean_stat=mean_stat, err_stat=err_stat)

        # Initialize the model depending on the fit

        if fit == 'OLS':
            self.model = statsmodels.api.OLS(self.y, self.X)
        elif fit == 'WLS':
            self.model = statsmodels.api.WLS(self.y, self.X, weights=self.w)
        else:
            self.model = statsmodels.api.RLM(self.y, self.X)

        # Perform the fit and return

        self.results = self.model.fit()

        return self.results
