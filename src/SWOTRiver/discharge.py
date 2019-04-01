"""
Module for computing discharge for river reaches
"""
import numpy as np

def area(observed_height, observed_width, area_fits):
    """
    Provides a nicer interface for _area wrapping up the unpacking of prior
    db area_fits into CalculatedAEIV.m inputs.

    observed_height - swot observed height for this reach
    observed_width - swot observed width for this reach
    area_fits - dictionary of things extracted from prior DB
    """
    height_breakpoints = np.squeeze(area_fits['h_break'])
    poly_fits = np.squeeze(area_fits['fit_coeffs'])
    poly_fits = np.squeeze(area_fits['fit_coeffs'])
    poly_fits = np.squeeze(area_fits['fit_coeffs'])

    area_median_flow = np.squeeze(area_fits['med_flow_area'])

    fit_width_std = np.squeeze(area_fits['w_err_stdev'])
    fit_height_std = np.squeeze(area_fits['h_err_stdev'])

    cov_height_width = np.zeros([2, 2])
    cov_height_width[0, 0] = np.squeeze(area_fits['w_variance'])
    cov_height_width[0, 1] = np.squeeze(area_fits['hw_covariance'])
    cov_height_width[1, 0] = cov_height_width[0, 1]
    cov_height_width[1, 1] = np.squeeze(area_fits['h_variance'])
    num_obs = np.squeeze(area_fits['h_w_nobs'])

    return _area(
        observed_height, observed_width, height_breakpoints, poly_fits,
        area_median_flow, fit_width_std**2, fit_height_std**2,
        cov_height_width, num_obs)

def _area(
    observed_height, observed_width, height_breakpoints, poly_fits,
    area_median_flow, fit_width_var, fit_height_var, cov_height_width,
    num_obs):
    """
    Computes cross-sectional area from fit, based on CalculatedAEIV.m at
    https://github.com/mikedurand/SWOTAprimeCalcs

    observed_height - swot observed height for this reach
    observed_width - swot observed width for this reach
    height_breakpoints - boundaries for fits in height
    poly_fits - polynominal coeffs for the fits
    area_median_flow - cross-sectional area at median flow
    fit_width_var - width error std**2
    fit_height_var - height error std**2
    cov_height_width - covariance matrix for width / height
    """
    poly_ints = []
    for poly_fit in poly_fits:
        poly_ints.append(np.polyint(poly_fit))
    poly_ints = np.array(poly_ints)

    height_fits_ll = height_breakpoints[0:-1]
    height_fits_ul = height_breakpoints[1:]

    ifit = np.argwhere(np.logical_and(
        observed_height >= height_fits_ll,
        observed_height < height_fits_ul))

    low_height_snr = (
        cov_height_width[1, 1] - fit_height_var)/fit_height_var < 2

    if ifit.size == 0:
        observed_height_hat = np.nan
        observed_width_hat = observed_width
        if observed_height > height_breakpoints.max():
            delta_area_hat = (
                np.polyval(poly_ints[-1], height_breakpoints[-1]) -
                np.polyval(poly_ints[-1], height_breakpoints[-2]) +
                area_median_flow)
            dAunc = np.sqrt(
                fit_height_var*observed_width**2 +
                2*fit_width_var*(observed_height-height_breakpoints[-1])**2)

        else:
            delta_area_hat = (
                - area_median_flow - ((height_breakpoints[0]-observed_height)
                * (observed_width + poly_fits[0][0]*height_breakpoints[0]
                + poly_fits[0][1])/2))
            dAunc = np.sqrt(
                fit_height_var*observed_width**2 +
                2*fit_width_var*(observed_height-height_breakpoints[0])**2)

    else:
        ifit = ifit[0][0]
        if low_height_snr:
            observed_height_hat = observed_height
        else:
            observed_height_hat = estimate_height(
                observed_width, observed_height, poly_fits[ifit],
                fit_width_var, fit_height_var)

        ifit_hat = np.argwhere(np.logical_and(
            observed_height_hat >= height_fits_ll,
            observed_height_hat < height_fits_ul))

        if ifit_hat.size > 0:
            ifit = ifit_hat[0][0]
            observed_height_hat = estimate_height(
                observed_width, observed_height, poly_fits[ifit],
                fit_width_var, fit_height_var)

        if low_height_snr:
            observed_width_hat = observed_width
        else:
            observed_width_hat = np.polyval(
                poly_fits[ifit], observed_height_hat)

        delta_area_hat = 0
        for poly_int, height_ll, height_ul in zip(
            poly_ints[:ifit+1], height_fits_ll[:ifit+1],
            height_fits_ul[:ifit+1]):

            delta_area_hat += (
                np.polyval(poly_int, np.min([observed_height_hat, height_ul]))
                - np.polyval(poly_int, height_ll))

        delta_area_hat -= area_median_flow

        if poly_fits[ifit][0] == 0:
            dAunc = poly_fits[ifit][1] * np.sqrt(fit_height_var)
        else:
            mu = (np.sqrt(
                poly_fits[ifit][0]/2) *
                (observed_height_hat - height_fits_ul[ifit]) + np.polyval(
                poly_fits[ifit], height_fits_ul[ifit]) / np.sqrt(
                2 * poly_fits[ifit][0]))
            sigma = np.sqrt(poly_fits[ifit][0]/2) * np.sqrt(fit_height_var)
            dAunc = np.sqrt(4*mu**2*sigma**2 + 2*sigma**4);

    return delta_area_hat, observed_width_hat, observed_height_hat, dAunc

def estimate_height(observed_width, observed_height, poly_fit, fit_width_var,
                    fit_height_var):
    """Estimates optimal height using error in variables approach"""
    #note this implements eqn. 1.3.17 in Fuller, assuming sigma_eu=0
    sigma_vv = fit_width_var + poly_fit[0]**2 * fit_height_var
    sigma_uv = -poly_fit[0] * fit_height_var
    v = observed_width - poly_fit[1] - poly_fit[0] * observed_height
    return observed_height - v * sigma_uv/sigma_vv
