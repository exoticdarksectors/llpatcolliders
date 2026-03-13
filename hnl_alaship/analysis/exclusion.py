"""
analysis/exclusion.py

Extract (m_N, U²) exclusion contour from N_signal scans.

The exclusion region is an "island" in (m_N, U²) space:
- Too small U²: insufficient production rate
- Too large U²: HNL decays before reaching the detector
"""

import numpy as np

from analysis.constants import N_THRESHOLD


def find_exclusion_band(u2_grid, N_grid, N_threshold=N_THRESHOLD):
    """
    Find the U² range where N_signal >= N_threshold.

    Parameters
    ----------
    u2_grid : (M,) log-spaced U² values
    N_grid : (M,) corresponding N_signal values
    N_threshold : float — minimum signal events for exclusion

    Returns
    -------
    dict with:
        u2_min : float — lower edge of exclusion band (NaN if none)
        u2_max : float — upper edge of exclusion band (NaN if none)
        peak_N : float — maximum N_signal over the scan
        peak_u2 : float — U² at peak sensitivity
        has_sensitivity : bool
    """
    peak_idx = np.argmax(N_grid)
    peak_N = float(N_grid[peak_idx])
    peak_u2 = float(u2_grid[peak_idx])

    mask = N_grid >= N_threshold
    if not np.any(mask):
        return {
            "u2_min": np.nan, "u2_max": np.nan,
            "peak_N": peak_N, "peak_u2": peak_u2,
            "has_sensitivity": False,
        }

    idx = np.where(mask)[0]
    i_lo, i_hi = idx[0], idx[-1]

    u2_min = _interpolate_threshold(u2_grid, N_grid, i_lo, N_threshold, "lower")
    u2_max = _interpolate_threshold(u2_grid, N_grid, i_hi, N_threshold, "upper")

    return {
        "u2_min": u2_min, "u2_max": u2_max,
        "peak_N": peak_N, "peak_u2": peak_u2,
        "has_sensitivity": True,
    }


def _interpolate_threshold(u2, N, idx, N_thr, side):
    """Interpolate in log(U²) space to find where N crosses the threshold."""
    if side == "lower" and idx > 0:
        N_below, N_above = N[idx - 1], N[idx]
        if N_above > N_below:
            frac = (N_thr - N_below) / (N_above - N_below)
            log_lo = np.log10(u2[idx - 1])
            log_hi = np.log10(u2[idx])
            return 10.0 ** (log_lo + frac * (log_hi - log_lo))
    elif side == "upper" and idx < len(u2) - 1:
        N_above, N_below = N[idx], N[idx + 1]
        if N_above > N_below:
            frac = (N_above - N_thr) / (N_above - N_below)
            log_lo = np.log10(u2[idx])
            log_hi = np.log10(u2[idx + 1])
            return 10.0 ** (log_lo + frac * (log_hi - log_lo))
    return float(u2[idx])
