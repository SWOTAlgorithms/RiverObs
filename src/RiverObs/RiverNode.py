"""
A river node holds all of the measurements associated with a node in a centerline.
It returns various data characteristics when queried.
"""

import numpy as N
from scipy.stats import trimboth

class RiverNode:
    """A river node holds all of the measurements associated with a node in a centerline.
    It returns various data characteristics when queried.

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
        across-track (normal) coordinate (relative to the node center) for each point
    h_flg: array_like
        flag indicating which pixels to use for height stats
    ds : float
        along-track dimension for this node. Defaults to 1. Needs to be set
        correctly for width_area to work.

    Notes
    -----

    To edit the data, it is useful to sort it by order in one of the variables.
    To sort the data, set sort=True and set sort_variable one of ['d','s','n','h'].
    If the data are not sorted intially, they can be sorted later, as desired by calling 
    RiverNode.sort.
    
    """

    # These variables are sorted simultaneously when sort is called. Append
    # to this list after class instantiation if you want to add an additional 
    # sort variable
    
    sort_vars = ['d','x','y','s','n']

    def __init__(self,index,d,x,y,s,n,ds=1,**kwds):
        """
        """

        self.index = index
        self.ds = ds
        self.d = N.asarray(d)
        self.x = N.asarray(x)
        self.y = N.asarray(y)
        self.s = N.asarray(s)
        self.n = N.asarray(n)

        for k in kwds:
            v = kwds[k]
            exec('self.%s = v'%k)

        self.ndata = len(self.x)
        self.good = N.ones(self.ndata,dtype=N.bool)
        self.sorted = False

    def add_obs(self,obs_name,obs,sort=True):
        """Add a new observations, and, if desired, sort them in accordance to the
        current sort index.

        obs_name: name of the observation, will be added as self.obs_name
        obs: iterable with observations of length self.ndata
        """

        if len(obs) != self.ndata:
            raise Exception('length of observations not consistent with number of node points')

        self.sort_vars.append(obs_name)
        
        if sort and self.sorted:
            exec('self.%s = N.asarray(obs[self.sort_index])'%obs_name)
        else:
            exec('self.%s = N.asarray(obs)'%obs_name)

    def count(self,*pars,**kwds):
        """Return the number of points in the node."""

        return self.ndata
    
    def countGood(self,goodvar):
        """Return the number of good points in the node."""
        tmp=N.zeros(self.ndata)
        exec('tmp[self.%s]=1'%goodvar)
        return sum(tmp)
    
    def value(self,var):
        """Return self.value. Useful when getting stats for all nodes and a scalar variable
        is desired."""

        value = 0 # fake cython compiler
        exec('value = self.%s'%var)
        
        return value

    def mean(self,var):
        """Return mean of a variable"""

        mean = 0 # fake cython compiler
        exec('mean = N.mean(self.%s[self.good])'%var)        
                 
        return mean
    # modified by Brent Williams, May 2017:
    # added goodvar to allow height aggregation to use different pixels than rest
    def median(self,var,goodvar='good'):
        """Return mean of a variable"""

        median = 0 # fake cython compiler
        #exec('median = N.median(self.%s[self.good])'%var)
        exec('median = N.median(self.%s[self.%s])'%(var,goodvar))  
                 
        return median

    def std(self,var,goodvar='good'):
        """Return std of a variable"""
        
        std = 0 # fake cython compiler
        #exec('std = N.std(self.%s[self.good])'%var)
        exec('std = N.std(self.%s[self.%s])'%(var,goodvar))
                 
        return std
    
    def stderr(self,var):
        """Return standrad error of a variable"""

        stderr = 0 # fake cython compiler
        scale = 1./N.sqrt(N.sum(self.good).astype(N.float))
        exec('stderr = scale*N.std(self.%s[self.good])'%var)        
                 
        return stderr

    def min(self,var):
        """Return min of a variable"""

        min = 0 # fake cython compiler
        exec('min = N.min(self.%s[self.good])'%var)        
                 
        return min

    def max(self,var):
        """Return max of a variable"""

        max = 0 # fake cython compiler
        exec('max = N.max(self.%s[self.good])'%var)        
                 
        return max

    def ptp(self,var):
        """Return peak to peak variation of a variable"""

        ptp = 0 # fake cython compiler
        exec('ptp = N.ptp(self.%s[self.good])'%var)        
                 
        return ptp

    def percentile(self,var,q):
        """Return qth percentile of a variable"""

        percentile = 0 # fake cython compiler
        exec('percentile = N.percentile(self.%s[self.good],q)'%var)        
                 
        return percentile

    def sum(self,var):
        """Return the sum of all variable values (e.g., for area)."""

        sum = 0 # fake cython compiler
        exec('sum = N.sum(self.%s[self.good])'%var)        
                 
        return sum

    def cdf(self,var):
        """Get the cdf for a variable."""

        x = 0 # fake cython compiler
        ngood = N.sum(self.good.astype(N.int32))
        cdf = N.cumsum(N.ones(ngood,dtype=N.float64))
        cdf /= cdf[-1]
        exec('x = N.sort(self.%s[self.good])'%var)

        return x, cdf

    def flag_extent(self,var_name,var_min,var_max):
        """Update the good flag to places where the variable var_name
        is within a range of values."""

        x = 0 # fake cython compiler
        exec('x = self.%s'%var_name)

        good = ( x >= var_min ) & ( x <= var_max)

        self.good = self.good & good

    def add_sort_variables(self,vars):
        """Add sort variables to the sort list. 

        vars is a string variable name, or an iterable of variable names."""

        if type(vars) == str:
            self.sort_vars.append(vars)
        else:
            for var in vars:
                self.sort_vars.append(var)

    def sort(self,sort_variable='n'):
        """Sort the data according to the variable sort_variable. self.sort_variable must exist.

        All of the variables in self.sort_vars are sorted in the same way.
        """

        # Get the index for sorting
        
        exec('self.sort_index = N.argsort(self.%s)'%sort_variable)

        # Sort all of the desired variables

        for var in self.sort_vars:
            exec('self.%s = self.%s[self.sort_index]'%(var,var))

        # Make sure the good flag is also sorted
            
        self.good = self.good[self.sort_index]

        # This flag is checked when calling 
        
        self.sorted = True

    def trim(self,fraction,mode='both'):
        """Calculate a index flag that trims the data, after sorting.

        fraction: 0 < f < 1. Fraction of the data to remove.
        mode is 'both', 'low', 'high' for which tails of the distribution need to be
        trimmed.
        """

        if not self.sorted:
            raise Exception('Run sort before calling trim')

        if mode == 'both':
            fraction = fraction/2.

        ntrim = int(self.ndata*fraction + 0.5)

        self.trim_index = N.ones(self.ndata, dtype=N.bool)

        if mode == 'both':
            self.trim_index[:ntrim] = False
            self.trim_index[:-ntrim] = False
        if mode == 'low':
            self.trim_index[:ntrim] = False
        else:
            self.trim_index[:ntrim] = False

        self.good = self.good & self.trim_index

        self.trimmed = True

    def width_ptp(self,*pars,**kwds):
        """Return the river width at this node, estimated from the max/min
        extent of the normal coordinate.

        Call signature made compatible with other stat calls, but requires
        no inputs."""

        nmin = self.min('n')
        nmax = self.max('n')
        width_ptp = nmax - nmin

        return width_ptp

    def width_std(self,*pars,**kwds):
        """Return the river width at this node, estimated as sqrt(12)*std(n).
        This would be appropriate for uniformly distributed point measurements.

        Call signature made compatible with other stat calls, but requires
        no inputs."""

        nstd = self.std('n')
        width_std = N.sqrt(12.)*nstd # uniformly distributed (~ 3.5 sigma Gaussian)

        return width_std

    def width_area(self,area_var='pixel_area'):
        """Return the river width at this node, estimated as Area/ds.
        This would be appropriate for pixels representing areas.

        area_var: str
            name of the variable associated with the pixel area (default 'pixel_area')

        Call signature made compatible with other stat calls, but requires
        no inputs."""

        area = self.sum(area_var)
        width_area = area/self.ds

        return width_area
