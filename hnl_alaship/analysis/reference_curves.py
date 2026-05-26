"""
analysis/reference_curves.py

Load reference exclusion curves from other PBC experiments
(MATHUSLA, ANUBIS, CODEX-b, PastExclusion) for comparison plotting.

File format: space-separated, 2 columns (mass_GeV, u2). Lines starting with
'#' are comments. A '# kind: contour' or '# kind: envelope' header marks
how the curve should be drawn:

- contour: rows trace a closed exclusion contour (lower edge L->R, upper
  edge R->L). Plot as a single connected line; fill the interior.
- envelope: rows are mass-sorted; the curve is the lower bound of an
  "everything above is excluded" region. Shade above the curve up to ymax.
"""

import numpy as np
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "vendored" / "reference_curves"


def _parse_kind(path):
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if "kind:" in line:
            return line.split("kind:", 1)[1].strip().lower()
    return "contour"


def load_reference_curve(experiment, flavor):
    """
    Load a reference exclusion curve.

    Returns dict {mass, u2, kind} or None if file not found.
    """
    path = DATA_DIR / f"{experiment}_{flavor}.dat"
    if not path.exists():
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        return None
    return {
        "mass": data[:, 0],
        "u2": data[:, 1],
        "kind": _parse_kind(path),
    }


def load_all_references(flavors=None):
    """
    Discover and load all reference curves in DATA_DIR.

    Returns dict: {experiment: {flavor: {mass, u2, kind}}}
    """
    if flavors is None:
        flavors = ["Ue", "Umu", "Utau"]

    if not DATA_DIR.exists():
        return {}

    experiments = set()
    for f in DATA_DIR.glob("*.dat"):
        name = f.stem
        for flav in flavors:
            if name.endswith(f"_{flav}"):
                experiments.add(name[: -(len(flav) + 1)])

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
