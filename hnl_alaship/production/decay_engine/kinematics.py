"""
production/decay_engine/kinematics.py

Vectorized phase-space decay functions for meson → HNL production.

Supports:
  - 2-body decay: M → ℓ N  (isotropic in parent rest frame)
  - 3-body decay: M → H' ℓ N  (flat Dalitz phase-space, no ME weighting)

All functions are vectorized over arrays of parent 4-vectors.
"""

import numpy as np


def _boost_to_lab(p4_rest, parent_E, parent_px, parent_py, parent_pz, parent_m):
    """
    Boost 4-vectors from parent rest frame to lab frame.

    Parameters
    ----------
    p4_rest : ndarray, shape (N, 4)
        4-vectors (E, px, py, pz) in parent rest frame.
    parent_E, parent_px, parent_py, parent_pz : ndarray, shape (N,)
        Parent 4-momentum in lab frame.
    parent_m : ndarray, shape (N,)
        Parent mass.

    Returns
    -------
    p4_lab : ndarray, shape (N, 4)
        Boosted 4-vectors in lab frame.
    """
    # Boost vector: β = p_parent / E_parent
    E_p = parent_E
    beta_x = parent_px / E_p
    beta_y = parent_py / E_p
    beta_z = parent_pz / E_p
    beta2 = beta_x**2 + beta_y**2 + beta_z**2
    beta2 = np.clip(beta2, 0, 1 - 1e-12)
    gamma = 1.0 / np.sqrt(1.0 - beta2)

    # Components in rest frame
    e_r = p4_rest[:, 0]
    px_r = p4_rest[:, 1]
    py_r = p4_rest[:, 2]
    pz_r = p4_rest[:, 3]

    # p · β
    bp = beta_x * px_r + beta_y * py_r + beta_z * pz_r

    # Handle beta2 ~ 0 (parent at rest)
    at_rest = beta2 < 1e-30

    e_lab = gamma * (e_r + bp)
    px_lab = px_r.copy()
    py_lab = py_r.copy()
    pz_lab = pz_r.copy()

    moving = ~at_rest
    if np.any(moving):
        k = (gamma[moving] - 1.0) * bp[moving] / beta2[moving]
        px_lab[moving] += beta_x[moving] * k + gamma[moving] * beta_x[moving] * e_r[moving]
        py_lab[moving] += beta_y[moving] * k + gamma[moving] * beta_y[moving] * e_r[moving]
        pz_lab[moving] += beta_z[moving] * k + gamma[moving] * beta_z[moving] * e_r[moving]

    return np.column_stack([e_lab, px_lab, py_lab, pz_lab])


def decay_2body(parent_E, parent_px, parent_py, parent_pz, m_parent, m1, m2, rng=None):
    """
    Isotropic 2-body decay: parent → daughter1 (mass m1) + daughter2 (mass m2).

    Parameters
    ----------
    parent_E, parent_px, parent_py, parent_pz : ndarray, shape (N,)
        Parent 4-momenta in lab frame.
    m_parent : float
        Parent mass in GeV.
    m1, m2 : float
        Daughter masses in GeV.
    rng : numpy.random.Generator, optional

    Returns
    -------
    d1 : ndarray, shape (N, 4) — (E, px, py, pz) of daughter 1 in lab frame
    d2 : ndarray, shape (N, 4) — (E, px, py, pz) of daughter 2 in lab frame
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)

    # Momentum magnitude in parent rest frame (2-body kinematics)
    M = m_parent
    p_star = np.sqrt(
        np.maximum((M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2), 0.0)
    ) / (2.0 * M)

    # Isotropic direction in rest frame
    cos_theta = rng.uniform(-1.0, 1.0, N)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    phi = rng.uniform(0.0, 2 * np.pi, N)

    px_star = p_star * sin_theta * np.cos(phi)
    py_star = p_star * sin_theta * np.sin(phi)
    pz_star = p_star * cos_theta

    E1_star_val = np.sqrt(p_star**2 + m1**2)
    E2_star_val = np.sqrt(p_star**2 + m2**2)
    E1_star = np.full(N, E1_star_val)
    E2_star = np.full(N, E2_star_val)

    d1_rest = np.column_stack([E1_star, px_star, py_star, pz_star])
    d2_rest = np.column_stack([E2_star, -px_star, -py_star, -pz_star])

    # Boost to lab
    m_arr = np.full(N, M)
    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz, m_arr)
    d2_lab = _boost_to_lab(d2_rest, parent_E, parent_px, parent_py, parent_pz, m_arr)

    return d1_lab, d2_lab


def decay_3body_flat(parent_E, parent_px, parent_py, parent_pz,
                     m_parent, m1, m2, m3, rng=None):
    """
    Flat 3-body phase-space decay (Dalitz): parent → d1(m1) + d2(m2) + d3(m3).

    No matrix-element weighting — uniform sampling in phase space using the
    sequential 2-body method:
      1. Sample m23 uniformly in [(m2+m3), (M-m1)]
      2. Decay parent → d1 + (23-system) as 2-body
      3. Decay (23-system) → d2 + d3 as 2-body

    This produces flat phase-space (up to the Jacobian |p*_1| × |p*_23|
    which is already the phase-space measure).

    Parameters
    ----------
    parent_E, parent_px, parent_py, parent_pz : ndarray, shape (N,)
    m_parent : float
    m1, m2, m3 : float
        Daughter masses in GeV.
    rng : numpy.random.Generator, optional

    Returns
    -------
    d1, d2, d3 : ndarray, shape (N, 4) — lab frame 4-vectors
    """
    if rng is None:
        rng = np.random.default_rng()

    N = len(parent_E)
    M = m_parent

    # Sample invariant mass of (23) system with correct phase-space weight.
    # The phase-space weight for flat sampling in m23 is |p*_1| × |p*_23|.
    # We use accept-reject on m23 to get the correct distribution.
    m23_min = m2 + m3
    m23_max = M - m1
    if m23_max <= m23_min:
        raise RuntimeError("No available 3-body phase space (at or below threshold).")

    # Use accept-reject: weight = p_star_1(m23) × p_star_23(m23)
    # Accept-reject envelope: max weight (found numerically per-event, but
    # since all events have same masses, one max suffices).

    # Precompute max weight by sampling densely
    m23_test = np.linspace(m23_min + 1e-6, m23_max - 1e-6, 1000)
    lam1 = (M**2 - m1**2 - m23_test**2)**2 - 4 * m1**2 * m23_test**2
    lam23 = (m23_test**2 - m2**2 - m3**2)**2 - 4 * m2**2 * m3**2
    p1_test = np.sqrt(np.maximum(lam1, 0.0)) / (2.0 * M)
    p23_test = np.sqrt(np.maximum(lam23, 0.0)) / (2.0 * np.maximum(m23_test, 1e-15))
    w_test = p1_test * p23_test
    w_max = w_test.max() * 1.01  # small safety margin
    if w_max <= 0:
        raise RuntimeError("No available 3-body phase space weight (at threshold).")

    # Accept-reject loop
    m23 = np.empty(N)
    accepted = np.zeros(N, dtype=bool)
    n_iter = 0
    max_iter = 10000
    while not accepted.all():
        n_iter += 1
        if n_iter > max_iter:
            raise RuntimeError("Failed to sample 3-body decay after max iterations.")
        n_need = (~accepted).sum()
        m23_cand = rng.uniform(m23_min, m23_max, n_need)
        lam1_c = (M**2 - m1**2 - m23_cand**2)**2 - 4 * m1**2 * m23_cand**2
        lam23_c = (m23_cand**2 - m2**2 - m3**2)**2 - 4 * m2**2 * m3**2
        p1_cand = np.sqrt(np.maximum(lam1_c, 0.0)) / (2.0 * M)
        p23_cand = np.sqrt(np.maximum(lam23_c, 0.0)) / (2.0 * np.maximum(m23_cand, 1e-15))
        w_cand = p1_cand * p23_cand
        u = rng.random(n_need) * w_max
        accept_mask = u < w_cand
        idx = np.where(~accepted)[0]
        m23[idx[accept_mask]] = m23_cand[accept_mask]
        accepted[idx[accept_mask]] = True

    # Step 1: decay parent → d1 + system23 in parent rest frame
    # Momentum of d1 in parent rest frame
    p1_star = np.sqrt(np.maximum(
        (M**2 - (m1 + m23)**2) * (M**2 - (m1 - m23)**2), 0.0
    )) / (2.0 * M)

    cos_th1 = rng.uniform(-1.0, 1.0, N)
    sin_th1 = np.sqrt(1.0 - cos_th1**2)
    phi1 = rng.uniform(0.0, 2 * np.pi, N)

    d1_px = p1_star * sin_th1 * np.cos(phi1)
    d1_py = p1_star * sin_th1 * np.sin(phi1)
    d1_pz = p1_star * cos_th1
    d1_E = np.sqrt(p1_star**2 + m1**2)

    sys23_E = np.sqrt(p1_star**2 + m23**2)
    sys23_px = -d1_px
    sys23_py = -d1_py
    sys23_pz = -d1_pz

    # Boost d1 and sys23 from parent rest frame to lab
    m_arr = np.full(N, M)
    d1_rest = np.column_stack([d1_E, d1_px, d1_py, d1_pz])
    d1_lab = _boost_to_lab(d1_rest, parent_E, parent_px, parent_py, parent_pz, m_arr)

    sys23_rest = np.column_stack([sys23_E, sys23_px, sys23_py, sys23_pz])
    sys23_lab = _boost_to_lab(sys23_rest, parent_E, parent_px, parent_py, parent_pz, m_arr)

    # Step 2: decay sys23 → d2 + d3 in sys23 rest frame
    p23_star = np.sqrt(np.maximum(
        (m23**2 - (m2 + m3)**2) * (m23**2 - (m2 - m3)**2), 0.0
    )) / (2.0 * m23)

    cos_th2 = rng.uniform(-1.0, 1.0, N)
    sin_th2 = np.sqrt(1.0 - cos_th2**2)
    phi2 = rng.uniform(0.0, 2 * np.pi, N)

    d2_px_r = p23_star * sin_th2 * np.cos(phi2)
    d2_py_r = p23_star * sin_th2 * np.sin(phi2)
    d2_pz_r = p23_star * cos_th2
    d2_E_r = np.sqrt(p23_star**2 + m2**2)
    d3_E_r = np.sqrt(p23_star**2 + m3**2)

    d2_in_sys23 = np.column_stack([d2_E_r, d2_px_r, d2_py_r, d2_pz_r])
    d3_in_sys23 = np.column_stack([d3_E_r, -d2_px_r, -d2_py_r, -d2_pz_r])

    # Boost d2, d3 from sys23 rest frame to lab
    d2_lab = _boost_to_lab(
        d2_in_sys23,
        sys23_lab[:, 0], sys23_lab[:, 1], sys23_lab[:, 2], sys23_lab[:, 3],
        m23
    )
    d3_lab = _boost_to_lab(
        d3_in_sys23,
        sys23_lab[:, 0], sys23_lab[:, 1], sys23_lab[:, 2], sys23_lab[:, 3],
        m23
    )

    return d1_lab, d2_lab, d3_lab
