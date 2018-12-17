'''
Description:
These are general functions to be used in the river node, raster, 
and lake processors for aggregating heights and areas with their 
corresponding uncertainties


Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Brent Williams
'''

import numpy as np

def simple(in_var, metric='mean'):
    """
    Aggregate the input variable according to desired metric/accumulator.
    
    INPUT:
    in_var = 1d list or array of a given variable for all pixels over a feature 
             (river node, raster bin, lake)
    """
    if metric == 'mean':
        out_var = np.nanmean(in_var)
    elif metric == 'median':
        out_var = np.nanmedian(in_var)
    elif metric == 'sum':
        out_var = np.nansum(in_var)
    elif metric == 'std':
        out_var = np.nanstd(in_var)
    elif metric == 'count':
        out_var = np.nansum(np.ones(np.shape(in_var)))
    return out_var
#
def height_only(height, good, height_std=1.0, method='weight'):
    """
    Return the aggregate height
    implements methods: weight (default), median, uniform 
    good is a mask used to filter out heights that are expected to be bad 
    or outliers, if desired
    
    INPUTS:
    height      = 1d array of heights over the water feature
    good        = mask for filtering out some pixels
    height_std  = pixel-wise height error std (phase_std * dh_dphase)
    method      = type of aggregator ('weight', 'uniform', or 'median')

    OUTPUTS:
    height_out  = aggregated height
    weight_norm = normalized weighting array

    Reference: implements Eq. (1) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    """
    if method=='median':
        # for median, use uniform weighting for uncertainty estimation later
        height_agg = simple(height[good], metric = 'median')
        num_pixels = simple(height[good], metric = 'count')
        weight = np.ones(np.shape(height))
        return height_agg, weight/num_pixels
    elif method=='uniform':
        weight = np.ones(np.shape(height)) 
    else:# method=='weight':
        # inverse height variance weighting
        weight = np.ones(np.shape(height))/(height_std)**2
    #
    height_agg = simple(weight[good]*height[good], metric = 'sum')
    weight_sum = simple(weight[good], metric = 'sum')
    num_pixels = simple(height[good], metric = 'count')
    
    weight_sum_pixc = np.ones(np.shape(weight))
    weight_sum_pixc[good] = weight_sum
    
    height_out = height_agg/weight_sum
    weight_norm = weight/weight_sum_pixc
    return height_out, weight_norm
#
def height_uncert_std(height, good, num_rare_looks, num_med_looks):
    """
    Compute the sample standard devieation of the heights and scale by the
    appropriate factor instead of 1/sqrt(N), since the medium pixels are
     correlated
    
    INPUTS:
    height         = 1d array of heights over the water feature
    good           = mask for filtering out some pixels
    num_rare_looks = rare number of looks (either actual looks, or effective)
    num_med_looks  = medium number of looks (either actual looks, or effective)

    Reference: implements Eq. (14) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    """
    h_std  = simple(height[good], metric = 'std')
    num_pixels = simple(height[good], metric = 'count')
    num_ind_pixels = simple(num_med_looks[good]/num_rare_looks[good],'mean')
    height_std_out = h_std * np.sqrt(num_ind_pixels/num_pixels)
    return height_std_out
#
def height_uncert_multilook(
        ifgram, power1, power2, weight_norm, good,
        num_rare_looks, looks_to_efflooks, dh_dphi):
    """
    compute height uncertainty bound by multilooking 
    everything over feature to compute average coh
    then computing phase noise std using CRB,
    then projecting through sensitivities

    INPUTS (all arrays are 1d lists for a given feature) 
    ifgram            = rare complex flattened interferogram
    power(1,2)        = the rare two channel interferpgram powers
    good              = mask for filtering out some pixels if desired
    weight_norm       = the normalized weighting function
    num_rare_looks    = rare looks array
    looks_to_efflooks = scale factor to get effective looks from looks taken
    dh_dphi           = height sensitivity to phase array
    
    OUTPUTS:
    height_uncert_out = scalar height uncertainty for this feature using this method

    Reference: implements square root of Eq. (7) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    """
    # multilook the rare interferogram over the raster bin
    #  by averaging cerain fields
    agg_real = simple(np.real(ifgram[good])*weight_norm[good])
    agg_imag = simple(np.imag(ifgram[good])*weight_norm[good])
    agg_p1 = simple(power1[good]*weight_norm[good])
    agg_p2 = simple(power2[good]*weight_norm[good])
    num_pixels = simple(power1[good],'count')
    # compute coherence
    coh = abs(agg_real + 1j *agg_imag)/np.sqrt(agg_p1*agg_p2)
    # get total num_eff_looks
    rare_looks = num_rare_looks/looks_to_efflooks
    agg_looks = simple(rare_looks[good])
    
    num_looks = agg_looks * num_pixels
    # get phase noise std using CRB
    phase_std = (0.5 / num_looks) * (1.0-coh**2)/(coh**2)
    agg_dh_dphi = simple(dh_dphi[good]*weight_norm[good])
    agg_dh_dphi2 = simple(dh_dphi[good]**2*weight_norm[good])
    
    height_uncert_out = np.sqrt(phase_std) * np.abs(np.array(agg_dh_dphi2)/np.array(agg_dh_dphi))
    return height_uncert_out
#
def height_with_uncerts(
        height,  good, num_rare_looks, num_med_looks,
        ifgram, power1, power2, look_to_efflooks, dh_dphi,
        height_std=1.0, method='weight'):
    """
    Return the aggregate height with corresponding uncertainty
    implements methods: weight(default), median, uniform
    good is a mask used to filter out heights that are expected
    to be bad or outliers, if desired
    """
    # first aggregate the heights
    height_out, weight_norm = height_only(
        height,  good, height_std=height_std, method=method)

    # now compute uncertainties
    height_std_out = height_uncert_std(
        height, good, num_rare_looks, num_med_looks)

    height_uncert_out = height_uncert_multilook(
        ifgram, power1, power2, weight_norm, good,
        num_rare_looks, look_to_efflooks, dh_dphi)
    
    return height_out, height_std_out, height_uncert_out
#


    

def area_only(pixel_area, klass, good, method='composite', 
              interior_water_klass=4, water_edge_klass=3, land_edge_klass=2):
    """
    Return the aggregate height
    implements methods: weight (default), median, uniform 
    good is a mask used to filter out heights that are expected to be bad 
    or outliers, if desired
    
    INPUTS:
    pixel_area  = 1d array of pixel_area
    klass       = classification, with edges
    good        = mask for filtering out some pixels
    method      = type of aggregator ('simple', 'water_fraction', 'composite')

    OUTPUTS:
    area_out  = aggregated height

    TODO: handle dark water

    Reference: implements Eq.s (15), (16), and (17) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    """
    Idw_in = np.zeros(np.shape(pixel_area))
    Idw_in[klass==interior_water_klass] = 1.0
    #
    Idw = np.zeros(np.shape(pixel_area))
    Idw[klass==interior_water_klass] = 1.0
    Idw[klass==water_edge_klass] = 1.0
    #
    Ide = np.zeros(np.shape(pixel_area))
    Ide[klass==water_edge_klass] = 1.0
    Ide[klass==land_edge_klass] = 1.0
    #
    I = np.zeros(np.shape(pixel_area))
    I[(Idw + Idw_in+ Ide) > 0] = 1.0 #all pixels near water
    
    if method=='simple':
        area_agg = simple(pixel_area[good]*Idw[good], metric = 'sum')
        num_pixels = simple(Idw[good], metric = 'sum')
    elif method=='water_fraction':
        area_agg = simple(pixel_area[good]*I[good], metric = 'sum')
        num_pixels = simple(I[good], metric = 'sum')
    else:# method=='composite':
        area_agg_in = simple(pixel_area[good]*Idw_in[good], metric = 'sum')
        area_agg_edge = simple(pixel_area[good]*Ide[good], metric = 'sum')
        area_agg = area_agg_in + area_agg_edge
        num_pixels = simple(Idw_in[good]+Ide[good], metric = 'sum')
    return area_agg, num_pixels
#

def area_with_uncert(pixel_area, klass, good, method='composite', 
              interior_water_klass=4, water_edge_klass=3, land_edge_klass=2):
    
    area_agg, num_pixels = area_only(pixel_area, klass, good, method=method, 
                                     interior_water_klass=interior_water_klass, 
                                     water_edge_klass=interior_water_klass, 
                                     land_edge_klass=interior_water_klass)
    # TODO call the uncertainty stuff
    return area_agg, None
