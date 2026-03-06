"""
analysis/sensitivity.py

Core HNL sensitivity calculation for GARGOYLE with full 2-body acceptance.

Replicates the cos θ* acceptance bounds from decayProbPerEvent_2body.py
(track separation and momentum cuts) and integrates over decay position
within the detector fiducial volume.

Physics:
    N_signal = L_int × U² × Σᵢ [ wᵢ × P_accepted_i ]

where P_accepted_i integrates the decay density × acceptance fraction
along the LLP path through the detector.
"""

import numpy as np
from scipy.integrate import quad

from analysis.constants import (
    SEP_MIN, SEP_MAX, P_CUT, M_ELECTRON,
    LOG_U2_MIN, LOG_U2_MAX, N_U2_POINTS,
)


# =========================================================================
# 2-body acceptance bounds (matching decayProbPerEvent_2body.py exactly)
# =========================================================================

def compute_c_upper(gamma, beta, mass, m_daughter=M_ELECTRON, p_cut=P_CUT):
    """
    Upper bound on |cos θ*| from momentum cut and forward requirement.

    The softer daughter in the lab frame must satisfy |p| > p_cut.
    Also capped at β (forward hemisphere in lab ↔ |cosθ*| < β).
    """
    E_min = np.sqrt(p_cut**2 + m_daughter**2)
    c_P = (1.0 - 2.0 * E_min / (gamma * mass)) / beta
    c_P = min(c_P, 1.0)
    return min(c_P, beta)


def compute_c_S(d_remaining, gamma, beta, sep_min=SEP_MIN):
    """
    Lower bound on |cos θ*| from minimum separation requirement.

    Tracks must be separated by at least sep_min at the detector wall.
    """
    if d_remaining <= 0:
        return 1.0
    theta_min = sep_min / d_remaining
    if theta_min <= 0:
        return 0.0
    cos_theta_min = np.cos(min(theta_min, np.pi))
    denom = gamma**2 * (1.0 - cos_theta_min)
    if denom <= 0:
        return 0.0
    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0
    return np.sqrt(val)


def compute_c_max_sep(d_remaining, gamma, beta, sep_max=SEP_MAX):
    """
    Upper bound on |cos θ*| from maximum separation requirement.

    Tracks must not be separated by more than sep_max at the detector wall.
    """
    if d_remaining <= 0:
        return 1.0
    theta_max = sep_max / d_remaining
    if theta_max >= np.pi:
        return 1.0
    cos_theta_max = np.cos(theta_max)
    denom = gamma**2 * (1.0 - cos_theta_max)
    if denom <= 0:
        return 1.0
    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0
    if val >= 1.0:
        return 1.0
    return np.sqrt(val)


# =========================================================================
# Vectorized 2-body acceptance bounds (for batch scan)
# =========================================================================

def compute_c_upper_vec(gamma, beta, mass, m_daughter=M_ELECTRON, p_cut=P_CUT):
    """Vectorized c_upper: (n_hits,) arrays → (n_hits,)."""
    E_min = np.sqrt(p_cut**2 + m_daughter**2)
    c_P = (1.0 - 2.0 * E_min / (gamma * mass)) / beta
    c_P = np.minimum(c_P, 1.0)
    return np.minimum(c_P, beta)


def compute_c_S_vec(d_remaining, gamma, beta, sep_min=SEP_MIN):
    """Vectorized c_S: arbitrary-shape arrays → same shape."""
    result = np.zeros_like(d_remaining)
    dead = d_remaining <= 0
    result[dead] = 1.0
    live = ~dead
    if not np.any(live):
        return result
    d_rem = d_remaining[live]
    g = gamma[live] if np.ndim(gamma) > 0 else gamma
    b = beta[live] if np.ndim(beta) > 0 else beta
    theta_min = np.minimum(sep_min / d_rem, np.pi)
    cos_theta_min = np.cos(theta_min)
    denom = g**2 * (1.0 - cos_theta_min)
    val = np.where(denom > 0, (1.0 - 2.0 / denom) / b**2, 0.0)
    result[live] = np.where(val > 0, np.sqrt(np.maximum(val, 0.0)), 0.0)
    return result


def compute_c_max_sep_vec(d_remaining, gamma, beta, sep_max=SEP_MAX):
    """Vectorized c_max_sep: arbitrary-shape arrays → same shape."""
    result = np.ones_like(d_remaining)
    dead = d_remaining <= 0
    # dead entries already 1.0
    live = ~dead
    if not np.any(live):
        return result
    d_rem = d_remaining[live]
    g = gamma[live] if np.ndim(gamma) > 0 else gamma
    b = beta[live] if np.ndim(beta) > 0 else beta
    theta_max = sep_max / d_rem
    large_theta = theta_max >= np.pi
    # For large_theta entries, result stays 1.0
    need_calc = ~large_theta
    if np.any(need_calc):
        t_max = theta_max[need_calc]
        g_c = g[need_calc] if np.ndim(g) > 0 else g
        b_c = b[need_calc] if np.ndim(b) > 0 else b
        cos_theta_max = np.cos(t_max)
        denom = g_c**2 * (1.0 - cos_theta_max)
        val = np.where(denom > 0, (1.0 - 2.0 / denom) / b_c**2, 0.0)
        calc_result = np.where(val <= 0, 0.0,
                               np.where(val >= 1.0, 1.0, np.sqrt(np.maximum(val, 0.0))))
        # Write back into the live positions that need calculation
        live_idx = np.flatnonzero(live)
        calc_idx = live_idx[need_calc]
        result.flat[calc_idx] = calc_result
    return result


# =========================================================================
# Per-event decay probability with acceptance
# =========================================================================

def acceptance_weighted_decay_prob(entry_d, exit_d, gamma, beta, mass,
                                  decay_length, m_daughter=M_ELECTRON,
                                  p_cut=P_CUT, sep_min=SEP_MIN,
                                  sep_max=SEP_MAX):
    """
    Probability that an LLP decays inside the detector AND passes
    2-body acceptance cuts.

    Integrates (1/λ) exp(-d/λ) × A(d) from entry_d to exit_d,
    where A(d) = max(0, c_upper - c_lower) is the accepted fraction
    of rest-frame cos θ* at decay position d.
    """
    if decay_length <= 0 or entry_d >= exit_d:
        return 0.0

    c_P_upper = compute_c_upper(gamma, beta, mass, m_daughter, p_cut)
    if c_P_upper <= 0:
        return 0.0

    def integrand(d):
        d_remaining = exit_d - d
        c_lower = compute_c_S(d_remaining, gamma, beta, sep_min)
        c_upper_sep = compute_c_max_sep(d_remaining, gamma, beta, sep_max)
        c_upper = min(c_P_upper, c_upper_sep)
        A = max(0.0, c_upper - c_lower)
        return (1.0 / decay_length) * np.exp(-d / decay_length) * A

    result, _ = quad(integrand, entry_d, exit_d, limit=100,
                     epsabs=1e-15, epsrel=1e-10)
    return result


def unweighted_decay_prob(entry_d, exit_d, decay_length):
    """
    Simple decay probability without acceptance cuts.

    P = exp(-entry_d/λ) × (1 - exp(-(exit_d - entry_d)/λ))
    """
    if decay_length <= 0 or entry_d >= exit_d:
        return 0.0
    return np.exp(-entry_d / decay_length) * (
        -np.expm1(-(exit_d - entry_d) / decay_length)
    )


# =========================================================================
# Signal computation
# =========================================================================

def compute_n_signal(weight, beta_gamma, gamma, beta,
                     hits, entry_d, exit_d,
                     ctau_u2_1, u2, m_N, L_int_pb,
                     m_daughter=M_ELECTRON,
                     use_acceptance=True):
    """
    Compute expected signal events at a given (m_N, U²) point.

    Parameters
    ----------
    weight : (N,) production weight in pb at U²=1
    beta_gamma : (N,) p / m_N
    gamma, beta : (N,) Lorentz factors
    hits : (N,) bool mask — events that hit the detector
    entry_d, exit_d : (N,) entry/exit distances in meters
    ctau_u2_1 : float — proper decay length at U²=1 in meters
    u2 : float — mixing angle squared
    m_N : float — HNL mass in GeV
    L_int_pb : float — integrated luminosity in pb⁻¹
    m_daughter : float — daughter particle mass for acceptance
    use_acceptance : bool — if False, use simple P_decay (no cuts)

    Returns
    -------
    float : N_signal
    """
    if u2 <= 0 or ctau_u2_1 <= 0:
        return 0.0

    ctau = ctau_u2_1 / u2
    mask = hits & np.isfinite(entry_d) & np.isfinite(exit_d)
    if not np.any(mask):
        return 0.0

    idx = np.where(mask)[0]
    total = 0.0
    for i in idx:
        lam = beta_gamma[i] * ctau  # lab-frame decay length
        if lam <= 0:
            continue

        if use_acceptance:
            p_decay = acceptance_weighted_decay_prob(
                entry_d[i], exit_d[i], gamma[i], beta[i], m_N,
                lam, m_daughter)
        else:
            p_decay = unweighted_decay_prob(entry_d[i], exit_d[i], lam)

        total += weight[i] * p_decay

    return L_int_pb * u2 * total


def scan_u2(weight, beta_gamma, gamma, beta,
            hits, entry_d, exit_d,
            ctau_u2_1, m_N, L_int_pb,
            m_daughter=M_ELECTRON,
            use_acceptance=True,
            log_u2_min=LOG_U2_MIN, log_u2_max=LOG_U2_MAX,
            n_points=N_U2_POINTS,
            br_vis=1.0):
    """
    Compute N_signal over a log-spaced grid of U² values.

    Parameters
    ----------
    use_acceptance : bool or str
        True  → vectorized trapezoid quadrature (fast, default)
        False → no acceptance cuts (pure exponential, fastest)
        "exact" → scipy.quad per event (slow, for validation)

    Returns
    -------
    u2_grid : (n_points,)
    N_grid : (n_points,)
    """
    u2_grid = np.logspace(log_u2_min, log_u2_max, n_points)

    # Pre-filter to only events that hit the detector
    mask = hits & np.isfinite(entry_d) & np.isfinite(exit_d)
    if not np.any(mask):
        return u2_grid, np.zeros(n_points)

    idx = np.where(mask)[0]
    w_hit = weight[idx]
    bg_hit = beta_gamma[idx]
    g_hit = gamma[idx]
    b_hit = beta[idx]
    entry_hit = entry_d[idx]
    exit_hit = exit_d[idx]

    if not use_acceptance:
        return _scan_u2_fast(
            w_hit, bg_hit, entry_hit, exit_hit,
            ctau_u2_1, L_int_pb, u2_grid, br_vis=br_vis)

    if use_acceptance == "exact":
        return _scan_u2_exact(
            w_hit, bg_hit, g_hit, b_hit, entry_hit, exit_hit,
            ctau_u2_1, m_N, L_int_pb, m_daughter, u2_grid, br_vis=br_vis)

    # Default: vectorized trapezoid quadrature
    return _scan_u2_vectorized(
        w_hit, bg_hit, g_hit, b_hit, entry_hit, exit_hit,
        ctau_u2_1, m_N, L_int_pb, m_daughter, u2_grid, br_vis=br_vis)


def _scan_u2_exact(w_hit, bg_hit, g_hit, b_hit, entry_hit, exit_hit,
                   ctau_u2_1, m_N, L_int_pb, m_daughter, u2_grid, br_vis=1.0):
    """Slow exact path using scipy.quad per event (for validation)."""
    n_points = len(u2_grid)
    N_grid = np.zeros(n_points)
    for iu, u2 in enumerate(u2_grid):
        ctau = ctau_u2_1 / u2
        total = 0.0
        for j in range(len(w_hit)):
            lam = bg_hit[j] * ctau
            if lam <= 0:
                continue
            p_acc = acceptance_weighted_decay_prob(
                entry_hit[j], exit_hit[j], g_hit[j], b_hit[j],
                m_N, lam, m_daughter)
            total += w_hit[j] * p_acc
        N_grid[iu] = L_int_pb * u2 * br_vis * total
    return u2_grid, N_grid


def _scan_u2_vectorized(w_hit, bg_hit, g_hit, b_hit, entry_hit, exit_hit,
                        ctau_u2_1, m_N, L_int_pb, m_daughter, u2_grid,
                        n_quad=30, br_vis=1.0):
    """
    Vectorized U² scan with 2-body acceptance using trapezoid quadrature.

    Replaces the double loop (U² × events × scipy.quad) with numpy arrays.
    """
    n_hits = len(w_hit)
    n_u2 = len(u2_grid)

    # Quadrature grid: (n_hits, n_quad) positions along each event's path
    t = np.linspace(0.0, 1.0, n_quad)[None, :]  # (1, n_quad)
    path_len = exit_hit - entry_hit               # (n_hits,)
    d_grid = entry_hit[:, None] + t * path_len[:, None]  # (n_hits, n_quad)
    d_remaining = exit_hit[:, None] - d_grid              # (n_hits, n_quad)

    # Acceptance fraction A(d) — independent of U²
    # c_P_upper: per-event (n_hits,) — depends only on gamma, beta, mass
    c_P = compute_c_upper_vec(g_hit, b_hit, m_N, m_daughter)  # (n_hits,)

    # c_S and c_max_sep: per-event per-position (n_hits, n_quad)
    g_2d = g_hit[:, None] * np.ones((1, n_quad))  # broadcast to (n_hits, n_quad)
    b_2d = b_hit[:, None] * np.ones((1, n_quad))
    c_lower = compute_c_S_vec(d_remaining, g_2d, b_2d)
    c_upper_sep = compute_c_max_sep_vec(d_remaining, g_2d, b_2d)
    c_upper = np.minimum(c_P[:, None], c_upper_sep)  # (n_hits, n_quad)
    A = np.maximum(0.0, c_upper - c_lower)           # (n_hits, n_quad)

    # Scan over U² grid
    N_grid = np.zeros(n_u2)
    for iu in range(n_u2):
        u2 = u2_grid[iu]
        lam = bg_hit * ctau_u2_1 / u2  # (n_hits,) lab-frame decay length
        # Decay density: (1/λ) exp(-d/λ)
        inv_lam = 1.0 / lam[:, None]                          # (n_hits, 1)
        density = inv_lam * np.exp(-d_grid * inv_lam)          # (n_hits, n_quad)
        integrand = density * A                                 # (n_hits, n_quad)
        P_acc = np.trapz(integrand, d_grid, axis=1)            # (n_hits,)
        N_grid[iu] = L_int_pb * u2 * br_vis * np.dot(w_hit, P_acc)

    return u2_grid, N_grid


def _scan_u2_fast(w_hit, bg_hit, entry_hit, exit_hit,
                  ctau_u2_1, L_int_pb, u2_grid, br_vis=1.0):
    """Vectorized U² scan without acceptance cuts (broadcasting)."""
    # lam[i, j] = bg[j] * ctau_u2_1 / u2[i]
    lam = bg_hit[None, :] * ctau_u2_1 / u2_grid[:, None]

    arg_entry = -entry_hit[None, :] / lam
    arg_path = -(exit_hit - entry_hit)[None, :] / lam
    P_decay = np.exp(arg_entry) * (-np.expm1(arg_path))

    N_grid = L_int_pb * u2_grid * br_vis * np.sum(w_hit[None, :] * P_decay, axis=1)
    return u2_grid, N_grid
