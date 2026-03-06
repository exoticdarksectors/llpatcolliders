"""
analysis/reference_curves.py

Load reference exclusion curves from other PBC experiments
(MATHUSLA, ANUBIS, CODEX-b, etc.) for comparison plotting.

Expected file format: space-separated, 3 columns (mass_GeV, u2_min, u2_max).
Lines starting with '#' are comments.

Files stored in: vendored/reference_curves/{experiment}_{flavor}.dat
"""

import numpy as np
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "vendored" / "reference_curves"


def load_reference_curve(experiment, flavor):
    """
    Load a reference (m_N, U²_min, U²_max) exclusion contour.

    Parameters
    ----------
    experiment : str
        e.g. "MATHUSLA", "ANUBIS", "CODEX-b"
    flavor : str
        "Ue", "Umu", or "Utau"

    Returns
    -------
    dict with mass, u2_min, u2_max arrays, or None if file not found.
    """
    path = DATA_DIR / f"{experiment}_{flavor}.dat"
    if not path.exists():
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        return None
    return {
        "mass": data[:, 0],
        "u2_min": data[:, 1],
        "u2_max": data[:, 2],
    }


def load_all_references(flavors=None):
    """
    Discover and load all reference curves in DATA_DIR.

    Returns
    -------
    dict : {experiment: {flavor: {mass, u2_min, u2_max}}}
    """
    if flavors is None:
        flavors = ["Ue", "Umu", "Utau"]

    if not DATA_DIR.exists():
        return {}

    # Discover experiments from filenames
    experiments = set()
    for f in DATA_DIR.glob("*.dat"):
        name = f.stem  # e.g. "MATHUSLA_Ue"
        for flav in flavors:
            if name.endswith(f"_{flav}"):
                exp = name[: -(len(flav) + 1)]
                experiments.add(exp)

    result = {}
    for exp in sorted(experiments):
        curves = {}
        for flav in flavors:
            c = load_reference_curve(exp, flav)
            if c is not None:
                curves[flav] = c
        if curves:
            result[exp] = curves

    return result
