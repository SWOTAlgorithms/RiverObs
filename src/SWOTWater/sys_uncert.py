import numpy as np

import SWOTWater.aggregate as agg

A = 2.20158638e-7 # Uncertainty bell curve coeff A
B = 1.39409027e-7 # Uncertainty bell curve coeff B
U_RES = 5 # Residual uncertainty over the ocean (in cm^2)
S = 0.02561919 # Slope in cross-track direction
GROUND_SPEED = 6.392002633921294 # SWOT ground speed (in km/s)

HEIGHT_SYS_UNCERT = 0.0482  # m
DEFAULT_XOVER_HEIGHT_SYS_UNCERT = 0.0755 # m
SLOPE_SYS_UNCERT = 0.000003  # m/m

def uncert_bell(p):
    """
    Uncertainty bell curve: p represents the position within the segment
    (-1: start, 0: middle, 1: end)
    """
    w_prev = 1 - np.sin(np.pi * p / 2)
    w_next = 1 + np.sin(np.pi * p / 2)
    return A * (w_prev * (p + 1) ** 2 + w_next * (p - 1) ** 2) + B

def xover_uncert_model(p, d, x):
    """
    The XCal uncertainty model takes 3 parameters as input
    p: position within the segment between xovers (in km)
    d: length of the segment (distance between consecutive xovers) (in km)
    x: cross-track distance (in km)
    returns the uncertainty resulting from the residual uncertainty and the
    model (in cm)
    """
    return np.ma.sqrt(U_RES + uncert_bell(p) * (S * x * d) ** 2)

def xover_uncert(xtrk_dist, time_from_prev_xover, time_to_next_xover):
    """ Compute systematic xover uncertainty (in m) """
    dist_prev = np.abs(GROUND_SPEED * time_from_prev_xover)
    dist_next = np.abs(GROUND_SPEED * time_to_next_xover)
    p = (dist_prev - dist_next) / (dist_prev + dist_next)
    d = dist_prev + dist_next
    return xover_uncert_model(p, d, xtrk_dist / 1000) / 100

def height_systematic_uncert(
        xtrk_dist, good, time_from_prev_xover=None, time_to_next_xover=None):
    """
    Get systematic uncertainty for height as the simple mean of the
    RSS of constant and xover terms
    """
    if good.any():
        if time_from_prev_xover is None or time_to_next_xover is None:
            xover_uncert_term = DEFAULT_XOVER_HEIGHT_SYS_UNCERT
        else:
            xover_uncert_term = xover_uncert(
                xtrk_dist[good], time_from_prev_xover[good],
                time_to_next_xover[good])

            # Fill any missing values with the default
            if np.ma.isMaskedArray(xover_uncert_term):
                xover_uncert_term = xover_uncert_term.filled(
                    DEFAULT_XOVER_HEIGHT_SYS_UNCERT)

        px_sys_uncert = np.sqrt(HEIGHT_SYS_UNCERT**2 + xover_uncert_term**2)
        sys_uncert = agg.simple(px_sys_uncert, metric='mean')
    else:
        # no good pixels in aggregation, return NaN
        sys_uncert = np.nan

    return sys_uncert
