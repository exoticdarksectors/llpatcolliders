"""
analysis/sensitivity.py

Signal-yield scan using a cached FairShip decay-acceptance kernel.
"""

from __future__ import annotations

import numpy as np

from analysis.constants import LOG_U2_MAX, LOG_U2_MIN, N_U2_POINTS


def _event_decay_probabilities(position_centers, position_widths, acceptance, beta_gamma, ctau_u2_1, u2):
    decay_lengths = beta_gamma * (ctau_u2_1 / u2)
    safe = (decay_lengths > 0.0) & np.isfinite(decay_lengths)
    result = np.zeros(len(beta_gamma), dtype=np.float64)
    if not np.any(safe):
        return result
    dl = decay_lengths[safe, None]
    density = np.exp(-position_centers[safe] / dl) / dl
    result[safe] = np.sum(density * acceptance[safe] * position_widths[safe], axis=1)
    return result


def compute_n_signal(cache, ctau_u2_1, u2, L_int_pb):
    """
    Compute the expected geometry-only signal yield at a single ``U^2`` point.
    """
    if u2 <= 0.0 or ctau_u2_1 <= 0.0 or cache["n_sampled"] == 0:
        return 0.0

    probabilities = _event_decay_probabilities(
        position_centers=cache["position_centers"],
        position_widths=cache["position_widths"],
        acceptance=cache["acceptance"],
        beta_gamma=cache["beta_gamma"],
        ctau_u2_1=ctau_u2_1,
        u2=u2,
    )

    estimator = cache["estimator"]
    if estimator == "exact":
        weighted_sum = float(np.dot(cache["sample_weights"], probabilities))
    elif estimator == "weighted_resample":
        weighted_sum = float(cache["total_hit_weight"] * probabilities.mean())
    else:
        weighted_sum = 0.0

    return L_int_pb * u2 * weighted_sum


def scan_u2(
    cache,
    ctau_u2_1,
    L_int_pb,
    log_u2_min=LOG_U2_MIN,
    log_u2_max=LOG_U2_MAX,
    n_points=N_U2_POINTS,
):
    """
    Compute ``N_signal(U^2)`` over a log-spaced scan.
    """
    u2_grid = np.logspace(log_u2_min, log_u2_max, int(n_points))
    if cache["n_sampled"] == 0 or ctau_u2_1 <= 0.0:
        return u2_grid, np.zeros_like(u2_grid)

    n_grid = np.array([compute_n_signal(cache, ctau_u2_1, u2, L_int_pb) for u2 in u2_grid])
    return u2_grid, n_grid
