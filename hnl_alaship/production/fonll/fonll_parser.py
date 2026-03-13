"""
production/fonll/fonll_parser.py

Parse vendored FONLL meson-level differential cross-section tables.

Files: 100×100 grid in (pT, y), 3 columns: pT  y  dσ/dpT/dy [pb/GeV]
pT: 0–50 GeV (100 bins), y: −3 to +3 (100 bins).
pT varies slowly (outer loop), y varies fast (inner loop).
"""

import numpy as np
from pathlib import Path

# Vendored FONLL table paths (relative to hnl_alaship/ root)
_FONLL_GRIDS_DIR = Path(__file__).parent.parent.parent / "vendored" / "fonll_grids"

FONLL_FILES = {
    "bottom": _FONLL_GRIDS_DIR / "fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat",
    "charm":  _FONLL_GRIDS_DIR / "fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat",
}


def parse_fonll_file(path):
    """
    Parse a FONLL meson-level dσ/dpT/dy table.

    Parameters
    ----------
    path : str or Path
        Path to the .dat file.

    Returns
    -------
    pt_arr : ndarray, shape (N_pt,)
        Unique pT values in GeV.
    y_arr : ndarray, shape (N_y,)
        Unique rapidity values.
    dsigma_2d : ndarray, shape (N_pt, N_y)
        dσ/dpT/dy in pb/GeV.
    """
    data = np.loadtxt(path, comments="#")
    pt_all = data[:, 0]
    y_all = data[:, 1]
    dsigma_all = data[:, 2]

    pt_unique = np.unique(pt_all)
    y_unique = np.unique(y_all)
    n_pt = len(pt_unique)
    n_y = len(y_unique)

    # Reshape: pT varies slowly (outer), y varies fast (inner)
    dsigma_2d = dsigma_all.reshape(n_pt, n_y)

    return pt_unique, y_unique, dsigma_2d


def get_sigma_total(quark):
    """
    Integrate FONLL dσ/dpT/dy over the full (pT, y) grid → total σ in pb.

    Uses the trapezoidal rule on the 2D grid.

    Parameters
    ----------
    quark : str
        "bottom" or "charm"

    Returns
    -------
    float
        Total meson-level cross-section in pb.
    """
    path = FONLL_FILES[quark]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)

    # Trapezoidal integration: first over y, then over pT
    _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # np.trapz removed in NumPy 2.0
    integral_over_y = _trapz(dsigma_2d, y_arr, axis=1)  # shape (N_pt,)
    sigma_total = _trapz(integral_over_y, pt_arr)        # scalar (pb)

    return sigma_total
