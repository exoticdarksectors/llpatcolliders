"""
production/fonll/meson_sampler.py

Inverse-CDF sampler from 2D FONLL dσ/dpT/dy grids.

Generates meson 4-vectors (E, px, py, pz) drawn from the FONLL
differential cross-section, with species assigned according to
fragmentation fractions.
"""

import numpy as np
from .fonll_parser import parse_fonll_file, FONLL_FILES
from ..constants import QUARK_MESON_MAP, MESON_MASSES


def _build_cdf(pt_arr, y_arr, dsigma_2d):
    """
    Build a 1D CDF from the 2D dσ/dpT/dy grid for inverse-CDF sampling.

    Each bin probability is dσ/dpT/dy × ΔpT × Δy (bin area).

    Returns
    -------
    cdf : ndarray, shape (N_pt * N_y,)
        Cumulative distribution function (normalized).
    pt_edges : ndarray
        Bin edges in pT.
    y_edges : ndarray
        Bin edges in rapidity.
    flat_idx_to_ij : callable
        Maps flat index → (i_pt, i_y) tuple.
    """
    n_pt = len(pt_arr)
    n_y = len(y_arr)

    # Compute bin widths (non-uniform spacing)
    dpt = np.diff(pt_arr)
    dy = np.diff(y_arr)

    # CDF weights: bin area = dsigma * dpt * dy.
    # Edge bins get half-width toward the boundary (pT=0 boundary for first bin).
    pt_cdf_widths = np.zeros(n_pt)
    pt_cdf_widths[0] = dpt[0] / 2.0
    pt_cdf_widths[-1] = dpt[-1]
    pt_cdf_widths[1:-1] = (dpt[:-1] + dpt[1:]) / 2.0

    y_cdf_widths = np.zeros(n_y)
    y_cdf_widths[0] = dy[0]
    y_cdf_widths[-1] = dy[-1]
    y_cdf_widths[1:-1] = (dy[:-1] + dy[1:]) / 2.0

    # Smearing widths: for sub-bin uniform smearing. The first pT bin
    # needs the full gap width so that the smeared range [0 - dpt/2, 0 + dpt/2]
    # clipped at pT >= 0 covers the correct interval [0, dpt/2].
    pt_smear_widths = pt_cdf_widths.copy()
    pt_smear_widths[0] = dpt[0]

    # 2D bin probabilities (unnormalized) — use CDF widths for correct area
    prob_2d = dsigma_2d * pt_cdf_widths[:, None] * y_cdf_widths[None, :]
    prob_2d = np.maximum(prob_2d, 0.0)  # clip negatives

    # Flatten and build CDF
    prob_flat = prob_2d.ravel()
    cdf = np.cumsum(prob_flat)
    norm = cdf[-1]
    if norm > 0:
        cdf /= norm

    return cdf, pt_arr, y_arr, pt_smear_widths, y_cdf_widths


def sample_meson_4vectors(n_events, quark, rng=None):
    """
    Sample meson 4-vectors from the FONLL dσ/dpT/dy distribution.

    For each sample:
      1. Draw (pT, y) bin from inverse CDF
      2. Uniform sub-bin smearing within the drawn bin
      3. Uniform φ in [0, 2π)
      4. Assign meson species from fragmentation fractions
      5. Compute 4-vector (E, px, py, pz) from (pT, y, φ, m)

    Parameters
    ----------
    n_events : int
        Number of meson 4-vectors to generate.
    quark : str
        "bottom" or "charm"
    rng : numpy.random.Generator, optional
        Random number generator (for reproducibility).

    Returns
    -------
    dict with keys:
        'E', 'px', 'py', 'pz' : ndarray, shape (n_events,)
            Meson 4-vectors in the lab frame.
        'species_pdg' : ndarray of int, shape (n_events,)
            PDG ID of the assigned meson species.
        'pt', 'y', 'phi' : ndarray, shape (n_events,)
            Kinematic variables used to construct the 4-vector.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Parse FONLL grid and build CDF
    path = FONLL_FILES[quark]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)
    cdf, _, _, pt_widths, y_widths = _build_cdf(pt_arr, y_arr, dsigma_2d)

    n_pt = len(pt_arr)
    n_y = len(y_arr)

    # Draw random samples from CDF
    u = rng.random(n_events)
    flat_idx = np.searchsorted(cdf, u)
    flat_idx = np.clip(flat_idx, 0, n_pt * n_y - 1)

    i_pt = flat_idx // n_y
    i_y = flat_idx % n_y

    # Sub-bin smearing: uniform within the bin
    pt_sampled = pt_arr[i_pt] + rng.uniform(-0.5, 0.5, n_events) * pt_widths[i_pt]
    y_sampled = y_arr[i_y] + rng.uniform(-0.5, 0.5, n_events) * y_widths[i_y]
    pt_sampled = np.maximum(pt_sampled, 0.0)  # pT >= 0

    # Uniform azimuthal angle
    phi = rng.uniform(0, 2 * np.pi, n_events)

    # Assign meson species from fragmentation fractions
    species_list = QUARK_MESON_MAP[quark]
    pdg_ids = np.array([s[0] for s in species_list])
    fractions = np.array([s[1] for s in species_list])
    frac_cumsum = np.cumsum(fractions)
    frac_cumsum /= frac_cumsum[-1]  # normalize

    u_species = rng.random(n_events)
    species_idx = np.searchsorted(frac_cumsum, u_species)
    species_idx = np.clip(species_idx, 0, len(pdg_ids) - 1)
    species_pdg = pdg_ids[species_idx]

    # Meson mass for each event
    m = np.array([MESON_MASSES[pdg] for pdg in species_pdg])

    # 4-vector from (pT, y, φ, m):
    #   mT = sqrt(pT² + m²)
    #   pz = mT * sinh(y)
    #   E  = mT * cosh(y)
    #   px = pT * cos(φ)
    #   py = pT * sin(φ)
    mt = np.sqrt(pt_sampled**2 + m**2)
    pz = mt * np.sinh(y_sampled)
    E = mt * np.cosh(y_sampled)
    px = pt_sampled * np.cos(phi)
    py = pt_sampled * np.sin(phi)

    return {
        'E': E, 'px': px, 'py': py, 'pz': pz,
        'species_pdg': species_pdg,
        'pt': pt_sampled, 'y': y_sampled, 'phi': phi,
    }
