"""
analysis/constants.py

Analysis-level constants for GARGOYLE HNL sensitivity calculation.
"""

# HL-LHC integrated luminosity
L_INT_FB = 3000.0           # fb⁻¹
L_INT_PB = L_INT_FB * 1e3   # pb⁻¹ = 3×10⁶ pb⁻¹

# Exclusion threshold: N_signal >= N_THRESHOLD
# 95% CL Poisson upper limit assuming zero background
N_THRESHOLD = 3.0

# U² scan range (log10 space)
LOG_U2_MIN = -12.0
LOG_U2_MAX = -1.0
N_U2_POINTS = 200

# CMS IP5 origin for ray-casting
CMS_ORIGIN = (0.0, 0.0, 0.0)

# Geometry-only charged-track acceptance cuts (from geometry/detector_cuts.py)
import sys as _sys
from pathlib import Path as _Path
_geom_dir = str(_Path(__file__).resolve().parent.parent.parent / "geometry")
if _geom_dir not in _sys.path:
    _sys.path.insert(0, _geom_dir)
from detector_cuts import P_CUT, SEP_MIN, SEP_MAX  # noqa: E402

# Default decay-kernel sampling settings
DEFAULT_MOTHER_SAMPLES = 5_000
DEFAULT_POSITION_BINS = 12
DEFAULT_DECAY_TEMPLATES = 20_000
DEFAULT_DECAYS_PER_BIN = 10
DEFAULT_RANDOM_SEED = 12_345
MIN_CHARGED_DAUGHTERS = 2

# FONLL meson channels close around 5 GeV, but WZ extends to ~80 GeV.
# Scan the full mass grid to include WZ-only high-mass region.
FONLL_MASS_MAX = 10.0  # GeV (full grid; WZ contributes above ~5 GeV)

# Flavors
FLAVORS = ["Ue", "Umu", "Utau"]
DEFAULT_ANALYSIS_FLAVORS = ["Umu"]
