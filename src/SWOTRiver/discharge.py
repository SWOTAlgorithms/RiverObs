"""
Module for computing discharge for river reaches
"""
import numpy as np
import warnings

from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9


def compute(reach, reach_height, reach_width, reach_slope):
    """Computes the discharge models"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        area_fit_outputs = area(
            reach_height, reach_width, reach.metadata['area_fits'])

    d_x_area = area_fit_outputs[0]
    if d_x_area < -10000000 or np.ma.is_masked(d_x_area):
        d_x_area = MISSING_VALUE_FLT

    d_x_area_u = area_fit_outputs[3]
    if d_x_area_u < 0 or np.ma.is_masked(d_x_area_u):
        d_x_area_u = MISSING_VALUE_FLT

    outputs = {'d_x_area': d_x_area, 'd_x_area_u': d_x_area_u}
    for key, models in reach.metadata['discharge_models'].items():

        metro_ninf = models['MetroMan']['ninf']
        metro_Abar = models['MetroMan']['Abar']
        metro_p = models['MetroMan']['p']

        if (reach_width > 0 and reach_slope > 0 and metro_Abar+d_x_area >= 0
                and metro_Abar > 0 and metro_ninf > 0):

            metro_n = metro_ninf * (
                (d_x_area+metro_Abar) / reach_width)**metro_p
            if reach.metadata['p_low_slp']:
                # Low slope flag is TRUE in PRD. Use different flow law.
                metro_q = metro_n * (reach_height - metro_Abar)**metro_p
            else:
                metro_q = (
                    (d_x_area+metro_Abar)**(5/3) * reach_width**(-2/3) *
                    reach_slope**(1/2)) / metro_n
        else:
            metro_q = MISSING_VALUE_FLT

        # 3: Compute BAM model
        bam_n = models['BAM']['n']
        bam_Abar = models['BAM']['Abar']

        if (reach_width > 0 and reach_slope > 0 and bam_Abar+d_x_area >= 0 and
            bam_Abar > 0 and bam_n > 0):

            bam_q = (
                (d_x_area+bam_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) / bam_n
        else:
            bam_q = MISSING_VALUE_FLT

        # 4: Compute HiVDI model
        hivdi_Abar = models['HiVDI']['Abar']
        hivdi_alpha = models['HiVDI']['alpha']
        hivdi_beta = models['HiVDI']['beta']

        if (reach_width > 0 and reach_slope > 0 and hivdi_Abar+d_x_area >= 0 and
            hivdi_Abar > 0 and hivdi_alpha > 0):
            hivdi_n_inv = hivdi_alpha * (
                (d_x_area+hivdi_Abar)/reach_width)**hivdi_beta
            hivdi_q = (
                (d_x_area+hivdi_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) * hivdi_n_inv
        else:
            hivdi_q = MISSING_VALUE_FLT

        # 5: Compute MOMMA model
        momma_B = models['MOMMA']['B']
        momma_H = models['MOMMA']['H']
        momma_Save = models['MOMMA']['Save']
        momma_r = 2

        momma_nb = 0.11 * momma_Save**0.18
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            log_factor = np.log10((momma_H-momma_B)/(reach_height-momma_B))

        if reach_height <= momma_H:
            momma_n = momma_nb*(1+log_factor)
            log_check = log_factor > -1
        else:
            momma_n = momma_nb*(1-log_factor)
            log_check = log_factor < 1

        if (reach_width > 0 and reach_slope > 0 and momma_n > 0 and
            momma_Save > 0 and momma_H > momma_B and momma_nb > 0
                and log_check):

            momma_q = (
                ((reach_height - momma_B)*(momma_r/(1+momma_r)))**(5/3) *
                reach_width * reach_slope**(1/2)) / momma_n
        else:
            momma_q = MISSING_VALUE_FLT

        # 6: Compute SADS model
        sads_Abar = models['SADS']['Abar']
        sads_n = models['SADS']['n']

        if (reach_width > 0 and reach_slope > 0 and sads_Abar+d_x_area >= 0 and
            sads_Abar > 0 and sads_n > 0):
            sads_q = (
                (d_x_area+sads_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) / sads_n
        else:
            sads_q = MISSING_VALUE_FLT

        # 7: Compute SIC4DVar model
        sic4dvar_n = models['SIC4DVar']['n']
        sic4dvar_Abar = models['SIC4DVar']['Abar']

        if (reach_width > 0 and reach_slope > 0 and sic4dvar_Abar+d_x_area >= 0
                and sic4dvar_Abar > 0 and sic4dvar_n > 0):

            sic4dvar_q = (
                (d_x_area+sic4dvar_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) / sic4dvar_n
        else:
            sic4dvar_q = MISSING_VALUE_FLT

        if key == 'constrained':
            outputs['metro_q_c'] = metro_q
            outputs['bam_q_c'] = bam_q
            outputs['hivdi_q_c'] = hivdi_q
            outputs['momma_q_c'] = momma_q
            outputs['sads_q_c'] = sads_q
            outputs['sic4dvar_q_c'] = sic4dvar_q
        elif key == 'unconstrained':
            outputs['metro_q_uc'] = metro_q
            outputs['bam_q_uc'] = bam_q
            outputs['hivdi_q_uc'] = hivdi_q
            outputs['momma_q_uc'] = momma_q
            outputs['sads_q_uc'] = sads_q
            outputs['sic4dvar_q_uc'] = sic4dvar_q

    # populate the constrained height and width outputs
    if np.isnan(area_fit_outputs[1]) or np.isnan(area_fit_outputs[2]):
        outputs['width_c'] = MISSING_VALUE_FLT
        outputs['width_c_u'] = MISSING_VALUE_FLT
        outputs['height_c'] = MISSING_VALUE_FLT
        outputs['height_c_u'] = MISSING_VALUE_FLT
    else:
        outputs['width_c'] = area_fit_outputs[1]
        outputs['width_c_u'] = MISSING_VALUE_FLT
        outputs['height_c'] = area_fit_outputs[2]
        outputs['height_c_u'] = MISSING_VALUE_FLT

    return outputs


def area(observed_height, observed_width, area_fits):
    """
    Provides a nicer interface for _area wrapping up the unpacking of prior
    db area_fits into CalculatedAEIV.m inputs.

    observed_height - swot observed height for this reach
    observed_width - swot observed width for this reach
    area_fits - dictionary of things extracted from prior DB
    """
    height_breakpoints = np.squeeze(area_fits['h_break'])
    poly_fits = [
        np.squeeze(area_fits['fit_coeffs'])[:, 0],
        np.squeeze(area_fits['fit_coeffs'])[:, 1],
        np.squeeze(area_fits['fit_coeffs'])[:, 2]]

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

    Inputs
    observed_height - swot observed height for this reach
    observed_width - swot observed width for this reach
    height_breakpoints - boundaries for fits in height
    poly_fits - polynominal coeffs for the fits
    area_median_flow - cross-sectional area at median flow
    fit_width_var - width error std**2
    fit_height_var - height error std**2
    cov_height_width - covariance matrix for width / height

    Outputs
    delta_area_hat - estimated cross-sectional area
    observed_width_hat - estimated width using height-width fit
    observed_height_hat - estimated height using height-width fit
    dAunc - uncertainty in the cross-sectional area

    TO-DO:
    add uncertainties observed_width_hat_u and observed_height_hat_u
    """
    poly_ints = np.array([np.polyint(item) for item in poly_fits])

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
    # TO-DO: add uncertainties for width_hat and height_hat
    return delta_area_hat, observed_width_hat, observed_height_hat, dAunc


def estimate_height(observed_width, observed_height, poly_fit, fit_width_var,
                    fit_height_var):
    """Estimates optimal height using error in variables approach"""
    #note this implements eqn. 1.3.17 in Fuller, assuming sigma_eu=0
    sigma_vv = fit_width_var + poly_fit[0]**2 * fit_height_var
    sigma_uv = -poly_fit[0] * fit_height_var
    v = observed_width - poly_fit[1] - poly_fit[0] * observed_height
    return observed_height - v * sigma_uv/sigma_vv
