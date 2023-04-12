"""
Module for computing discharge for river reaches
"""
import numpy as np
import warnings

from RiverObs.RiverObs import \
    MISSING_VALUE_FLT, MISSING_VALUE_INT4, MISSING_VALUE_INT9


def compute(reach, reach_height, reach_height_u, reach_width, reach_width_u,
            reach_slope, reach_slope_u):
    """Computes the discharge models"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d_x_area, d_x_area_u, width_c, width_c_u, height_c, height_c_u = area(
            reach_height, reach_height_u, reach_width, reach_width_u,
            reach.metadata['area_fits'])

    if d_x_area < -10000000 or np.ma.is_masked(d_x_area):
        d_x_area = MISSING_VALUE_FLT

    if d_x_area_u < 0 or np.ma.is_masked(d_x_area_u):
        d_x_area_u = MISSING_VALUE_FLT

    outputs = {'d_x_area': d_x_area, 'd_x_area_u': d_x_area_u}
    for key, models in reach.metadata['discharge_models'].items():

        metro_ninf = models['MetroMan']['ninf']
        metro_Abar = models['MetroMan']['Abar']
        metro_p = models['MetroMan']['p']
        metro_s_rel_u = models['MetroMan']['sbQ_rel'].item()

        if (reach_width > 0 and reach_slope > 0 and metro_Abar+d_x_area >= 0
                and metro_Abar > 0 and metro_ninf > 0):

            metro_n = metro_ninf * (
                (d_x_area+metro_Abar) / reach_width)**metro_p
            if reach.metadata['p_low_slp']:
                # Low slope flag is TRUE in PRD. Use different flow law.
                metro_q = metro_n * (reach_height - metro_Abar)**metro_p
                # to-do, random error uncertainty placeholder
                metro_r_u = MISSING_VALUE_FLT
            else:
                metro_q = (
                    (d_x_area+metro_Abar)**(5/3) * reach_width**(-2/3) *
                    reach_slope**(1/2)) / metro_n
                # metroman opt 4 uncertainties
                metro_width_u = (
                        ((metro_p - (2/3)) / reach_width) * reach_width_u)
                metro_d_x_area_u = ((5/3) - metro_p) / (
                        metro_Abar + d_x_area) * d_x_area_u
                metro_slp_u = reach_slope_u / (2 * reach_slope)
                metro_r_u = np.sqrt(
                    metro_width_u**2 + metro_d_x_area_u**2 + metro_slp_u**2)
            if 0 <= metro_s_rel_u < 1:
                metro_s_u, metro_u = discharge_uncertainty(metro_s_rel_u,
                                                           metro_r_u)
            else:
                metro_s_rel_u = MISSING_VALUE_FLT
                metro_s_u = MISSING_VALUE_FLT
                metro_u = MISSING_VALUE_FLT

        else:
            metro_q = MISSING_VALUE_FLT
            metro_s_rel_u = MISSING_VALUE_FLT
            metro_s_u = MISSING_VALUE_FLT
            metro_r_u = MISSING_VALUE_FLT
            metro_u = MISSING_VALUE_FLT

        # 3: Compute BAM model
        bam_n = models['BAM']['n']
        bam_Abar = models['BAM']['Abar']
        bam_s_rel_u = models['BAM']['sbQ_rel'].item()

        if (reach_width > 0 and reach_slope > 0 and bam_Abar+d_x_area >= 0 and
            bam_Abar > 0 and bam_n > 0):

            bam_q = (
                (d_x_area+bam_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) / bam_n
            bam_width_u = (2 * reach_width_u) / (3 * reach_width)
            bam_slp_u = reach_slope_u / (2 * reach_slope)
            bam_d_x_area_u = 5 * d_x_area_u / (3 * (bam_Abar + d_x_area))
            bam_r_u = np.sqrt(bam_width_u**2 + bam_slp_u**2 +
                              bam_d_x_area_u**2)
            if 0 <= bam_s_rel_u < 1:
                bam_s_u, bam_u = discharge_uncertainty(bam_s_rel_u, bam_r_u)
            else:
                bam_s_rel_u = MISSING_VALUE_FLT
                bam_s_u = MISSING_VALUE_FLT
                bam_u = MISSING_VALUE_FLT

        else:
            bam_q = MISSING_VALUE_FLT
            bam_s_rel_u = MISSING_VALUE_FLT
            bam_s_u = MISSING_VALUE_FLT
            bam_r_u = MISSING_VALUE_FLT
            bam_u = MISSING_VALUE_FLT

        # 4: Compute HiVDI model
        hivdi_Abar = models['HiVDI']['Abar']
        hivdi_alpha = models['HiVDI']['alpha']
        hivdi_beta = models['HiVDI']['beta']
        hivdi_s_rel_u = models['HiVDI']['sbQ_rel'].item()

        if (reach_width > 0 and reach_slope > 0 and hivdi_Abar+d_x_area >= 0 and
            hivdi_Abar > 0 and hivdi_alpha > 0):

            hivdi_n_inv = hivdi_alpha * (
                (d_x_area+hivdi_Abar)/reach_width)**hivdi_beta
            hivdi_q = (
                (d_x_area+hivdi_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) * hivdi_n_inv
            hivdi_width_u = reach_width_u * (2/3 + hivdi_beta) / reach_width
            hivdi_slp_u = reach_slope_u / (2 * reach_slope)
            hivdi_d_x_area_u = (
                    d_x_area_u * (5/3 + hivdi_beta) / (hivdi_Abar + d_x_area))
            hivdi_r_u = np.sqrt(hivdi_width_u**2 + hivdi_slp_u**2 +
                                hivdi_d_x_area_u**2)
            if 0 <= hivdi_s_rel_u < 1:
                hivdi_s_u, hivdi_u = discharge_uncertainty(hivdi_s_rel_u,
                                                           hivdi_r_u)
            else:
                hivdi_s_rel_u = MISSING_VALUE_FLT
                hivdi_s_u = MISSING_VALUE_FLT
                hivdi_u = MISSING_VALUE_FLT

        else:
            hivdi_q = MISSING_VALUE_FLT
            hivdi_s_rel_u = MISSING_VALUE_FLT
            hivdi_s_u = MISSING_VALUE_FLT
            hivdi_r_u = MISSING_VALUE_FLT
            hivdi_u = MISSING_VALUE_FLT

        # 5: Compute MOMMA model
        momma_B = models['MOMMA']['B']
        momma_H = models['MOMMA']['H']
        momma_Save = models['MOMMA']['Save']
        momma_r = 2
        momma_s_rel_u = models['MOMMA']['sbQ_rel'].item()

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
            momma_width_u = reach_width_u / reach_width
            momma_slp_u = reach_slope_u / (2 * reach_slope)
            momma_wse_u = 5 * reach_height_u / (3 * (reach_height - momma_B))
            momma_r_u = np.sqrt(momma_width_u**2 + momma_slp_u**2 +
                                momma_wse_u**2)
            if 0 <= momma_s_rel_u < 1:
                momma_s_u, momma_u = discharge_uncertainty(momma_s_rel_u,
                                                           momma_r_u)
            else:
                momma_s_rel_u = MISSING_VALUE_FLT
                momma_s_u = MISSING_VALUE_FLT
                momma_u = MISSING_VALUE_FLT

        else:
            momma_q = MISSING_VALUE_FLT
            momma_s_rel_u = MISSING_VALUE_FLT
            momma_s_u = MISSING_VALUE_FLT
            momma_r_u = MISSING_VALUE_FLT
            momma_u = MISSING_VALUE_FLT

        # 6: Compute SADS model
        sads_Abar = models['SADS']['Abar']
        sads_n = models['SADS']['n']
        sads_s_rel_u = models['SADS']['sbQ_rel'].item()

        if (reach_width > 0 and reach_slope > 0 and sads_Abar+d_x_area >= 0 and
            sads_Abar > 0 and sads_n > 0):

            sads_q = (
                (d_x_area+sads_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) / sads_n
            sads_width_u = (2 * reach_width_u) / (3 * reach_width)
            sads_slp_u = reach_slope_u / (2 * reach_slope)
            sads_d_x_area_u = 5 * d_x_area_u / (3 * (sads_Abar + d_x_area))
            sads_r_u = np.sqrt(sads_width_u**2 + sads_slp_u**2 +
                               sads_d_x_area_u**2)
            if 0 <= sads_s_rel_u < 1:
                sads_s_u, sads_u = discharge_uncertainty(sads_s_rel_u,
                                                         sads_r_u)
            else:
                sads_s_rel_u = MISSING_VALUE_FLT
                sads_s_u = MISSING_VALUE_FLT
                sads_u = MISSING_VALUE_FLT

        else:
            sads_q = MISSING_VALUE_FLT
            sads_s_rel_u = MISSING_VALUE_FLT
            sads_s_u = MISSING_VALUE_FLT
            sads_r_u = MISSING_VALUE_FLT
            sads_u = MISSING_VALUE_FLT

        # 7: Compute SIC4DVar model
        sic4dvar_n = models['SIC4DVar']['n']
        sic4dvar_Abar = models['SIC4DVar']['Abar']
        sic4dvar_s_rel_u = models['SIC4DVar']['sbQ_rel'].item()

        if (reach_width > 0 and reach_slope > 0 and sic4dvar_Abar+d_x_area >= 0
                and sic4dvar_Abar > 0 and sic4dvar_n > 0):

            sic4dvar_q = (
                (d_x_area+sic4dvar_Abar)**(5/3) * reach_width**(-2/3) *
                (reach_slope)**(1/2)) / sic4dvar_n
            sic4dvar_width_u = (2 * reach_width_u) / (3 * reach_width)
            sic4dvar_slp_u = reach_slope_u / (2 * reach_slope)
            sic4dvar_d_x_area_u = 5 * d_x_area_u / (
                        3 * (sic4dvar_Abar + d_x_area))
            sic4dvar_r_u = np.sqrt(
                sic4dvar_width_u**2 + sic4dvar_slp_u**2 +
                sic4dvar_d_x_area_u**2)
            if 0 <= sic4dvar_s_rel_u < 1:
                sic4dvar_s_u, sic4dvar_u = discharge_uncertainty(
                    sic4dvar_s_rel_u, sic4dvar_r_u)
            else:
                sic4dvar_s_rel_u = MISSING_VALUE_FLT
                sic4dvar_s_u = MISSING_VALUE_FLT
                sic4dvar_u = MISSING_VALUE_FLT

        else:
            sic4dvar_q = MISSING_VALUE_FLT
            sic4dvar_s_rel_u = MISSING_VALUE_FLT
            sic4dvar_s_u = MISSING_VALUE_FLT
            sic4dvar_r_u = MISSING_VALUE_FLT
            sic4dvar_u = MISSING_VALUE_FLT

        # 8: Compute consensus discharge and its uncertainties
        q_results = np.ma.masked_values([metro_q, bam_q, hivdi_q, momma_q,
                                         sads_q, sic4dvar_q],
                                        MISSING_VALUE_FLT)
        q_r_u = np.ma.masked_values([metro_r_u, bam_r_u, hivdi_r_u, momma_r_u,
                                     sads_r_u, sic4dvar_r_u],
                                    MISSING_VALUE_FLT)
        q_s_u = np.ma.masked_values([metro_s_u, bam_s_u, hivdi_s_u,
                                     momma_s_u, sads_s_u, sic4dvar_s_u],
                                    MISSING_VALUE_FLT)
        nalgo = np.sum(q_results.mask == False)
        if nalgo >= 1:
            consensus_q = np.median(q_results)
            consensus_s_u = np.sqrt(np.pi / 2 * np.mean(q_s_u)**2 / nalgo)
            consensus_u = np.sqrt(consensus_s_u**2 + np.median(q_r_u)**2)
            consensus_s_rel_u = consensus_s_u / consensus_u
        else:
            consensus_q = MISSING_VALUE_FLT
            consensus_s_rel_u = MISSING_VALUE_FLT
            consensus_u = MISSING_VALUE_FLT

        # Note we have no algorithm yet for _q (qual)
        # discharge vars, they are set to missing_value.
        if key == 'constrained':
            outputs['dschg_gm'] = metro_q
            outputs['dschg_gm_u'] = metro_u
            outputs['dschg_gmsf'] = metro_s_rel_u
            outputs['dschg_gm_q'] = MISSING_VALUE_INT4

            outputs['dschg_gb'] = bam_q
            outputs['dschg_gb_u'] = bam_u
            outputs['dschg_gbsf'] = bam_s_rel_u
            outputs['dschg_gb_q'] = MISSING_VALUE_INT4

            outputs['dschg_gh'] = hivdi_q
            outputs['dschg_gh_u'] = hivdi_u
            outputs['dschg_ghsf'] = hivdi_s_rel_u
            outputs['dschg_gh_q'] = MISSING_VALUE_INT4

            outputs['dschg_go'] = momma_q
            outputs['dschg_go_u'] = momma_u
            outputs['dschg_gosf'] = momma_s_rel_u
            outputs['dschg_go_q'] = MISSING_VALUE_INT4

            outputs['dschg_gs'] = sads_q
            outputs['dschg_gs_u'] = sads_u
            outputs['dschg_gssf'] = sads_s_rel_u
            outputs['dschg_gs_q'] = MISSING_VALUE_INT4

            outputs['dschg_gi'] = sic4dvar_q
            outputs['dschg_gi_u'] = sic4dvar_u
            outputs['dschg_gisf'] = sic4dvar_s_rel_u
            outputs['dschg_gi_q'] = MISSING_VALUE_INT4

            outputs['dschg_gc'] = consensus_q
            outputs['dschg_gc_u'] = consensus_u
            outputs['dschg_gcsf'] = consensus_s_rel_u
            outputs['dschg_gc_q'] = MISSING_VALUE_INT4

        elif key == 'unconstrained':
            outputs['dschg_m'] = metro_q
            outputs['dschg_m_u'] = metro_u
            outputs['dschg_msf'] = metro_s_rel_u
            outputs['dschg_m_q'] = MISSING_VALUE_INT4

            outputs['dschg_b'] = bam_q
            outputs['dschg_b_u'] = bam_u
            outputs['dschg_bsf'] = bam_s_rel_u
            outputs['dschg_b_q'] = MISSING_VALUE_INT4

            outputs['dschg_h'] = hivdi_q
            outputs['dschg_h_u'] = hivdi_u
            outputs['dschg_hsf'] = hivdi_s_rel_u
            outputs['dschg_h_q'] = MISSING_VALUE_INT4

            outputs['dschg_o'] = momma_q
            outputs['dschg_o_u'] = momma_u
            outputs['dschg_osf'] = momma_s_rel_u
            outputs['dschg_o_q'] = MISSING_VALUE_INT4

            outputs['dschg_s'] = sads_q
            outputs['dschg_s_u'] = sads_u
            outputs['dschg_ssf'] = sads_s_rel_u
            outputs['dschg_s_q'] = MISSING_VALUE_INT4

            outputs['dschg_i'] = sic4dvar_q
            outputs['dschg_i_u'] = sic4dvar_u
            outputs['dschg_isf'] = sic4dvar_s_rel_u
            outputs['dschg_i_q'] = MISSING_VALUE_INT4

            outputs['dschg_c'] = consensus_q
            outputs['dschg_c_u'] = consensus_u
            outputs['dschg_csf'] = consensus_s_rel_u
            outputs['dschg_c_q'] = MISSING_VALUE_INT4

    # populate the constrained height and width outputs
    if np.isnan(width_c) or np.isnan(height_c):
        outputs['width_c'] = MISSING_VALUE_FLT
        outputs['width_c_u'] = MISSING_VALUE_FLT
        outputs['height_c'] = MISSING_VALUE_FLT
        outputs['height_c_u'] = MISSING_VALUE_FLT
    else:
        outputs['width_c'] = width_c
        outputs['width_c_u'] = width_c_u
        outputs['height_c'] = height_c
        outputs['height_c_u'] = height_c_u

    return outputs


def area(observed_height, observed_height_u, observed_width, observed_width_u,
         area_fits):
    """
    Provides a nicer interface for _area wrapping up the unpacking of prior
    db area_fits into CalculatedAEIV.m inputs.

    observed_height - swot observed height for this reach
    observed_height_u - swot observed height uncertainty for this reach
    observed_width - swot observed width for this reach
    observed_width_u - swot observed width uncertainty for this reach
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
        observed_height, observed_height_u, observed_width, observed_width_u,
        height_breakpoints, poly_fits, area_median_flow, fit_width_std**2,
        fit_height_std**2, cov_height_width, num_obs)


def _area(
    observed_height, observed_height_u, observed_width, observed_width_u,
        height_breakpoints, poly_fits, area_median_flow, fit_width_var,
        fit_height_var, cov_height_width, num_obs):
    """
    Computes cross-sectional area from fit, based on CalculatedAEIV.m at
    https://github.com/mikedurand/SWOTAprimeCalcs

    Inputs
    observed_height - swot observed height for this reach
    observed_height_u - swot observed height uncertainty for this reach
    observed_width - swot observed width for this reach
    observed_width_u - swot observed width uncertainty for this reach
    height_breakpoints - boundaries for fits in height
    poly_fits - polynominal coeffs for the fits
    area_median_flow - cross-sectional area at median flow
    fit_width_var - width error std**2
    fit_height_var - height error std**2
    cov_height_width - covariance matrix for width / height

    Outputs
    delta_area_hat - estimated cross-sectional area
    dAunc - uncertainty in the cross-sectional area
    observed_width_hat - estimated width using height-width fit
    observed_width_hat_u - observed_width_hat uncertainty
    observed_height_hat - estimated height using height-width fit
    observed_height_hat_u - observed_height_hat uncertainty
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
        observed_height_hat_u = np.nan
        observed_width_hat = observed_width
        observed_width_hat_u = observed_width_u
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
            observed_height_hat_u = observed_height_u
        else:
            observed_height_hat, observed_height_hat_u = estimate_height(
                observed_width, observed_height, poly_fits[ifit],
                fit_width_var, fit_height_var)

        ifit_hat = np.argwhere(np.logical_and(
            observed_height_hat >= height_fits_ll,
            observed_height_hat < height_fits_ul))

        if ifit_hat.size > 0:
            ifit = ifit_hat[0][0]
            observed_height_hat, observed_height_hat_u = estimate_height(
                observed_width, observed_height, poly_fits[ifit],
                fit_width_var, fit_height_var)

        if low_height_snr:
            observed_width_hat = observed_width
            observed_width_hat_u = observed_width_u
        else:
            observed_width_hat = np.polyval(
                poly_fits[ifit], observed_height_hat)
            observed_width_hat_u = observed_height_hat_u * poly_fit[ifit][0]

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
            dAunc = np.sqrt(4*mu**2*sigma**2 + 2*sigma**4)

    return delta_area_hat, dAunc, observed_width_hat, observed_width_hat_u,\
           observed_height_hat, observed_height_hat_u


def estimate_height(observed_width, observed_height, poly_fit, fit_width_var,
                    fit_height_var):
    """Estimates optimal height using error in variables approach"""
    #note this implements eqn. 1.3.17 in Fuller, assuming sigma_eu=0
    sigma_vv = fit_width_var + poly_fit[0]**2 * fit_height_var
    sigma_uv = -poly_fit[0] * fit_height_var
    v = observed_width - poly_fit[1] - poly_fit[0] * observed_height
    observed_height_hat = observed_height - v * sigma_uv/sigma_vv
    observed_height_hat_u = np.sqrt(fit_height_var - sigma_uv**2/sigma_vv)
    return observed_height_hat, observed_height_hat_u


def discharge_uncertainty(sbq_rel_u, q_random_u):
    """
    Compute systematic and total discharge uncertainties.

    Inputs:
    sbq_rel_u - Fractional systematic uncertainty in discharge.
    q_random_u - Random uncertainty in discharge.

    Outputs:
    syst_u - Systematic discharge uncertainty.
    tot_u - Total discharge uncertainty.
    """
    syst_u = sbq_rel_u * q_random_u / np.sqrt(1 - sbq_rel_u**2)
    tot_u = np.sqrt(q_random_u**2 + syst_u**2)
    return syst_u, tot_u
