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

# 2-body acceptance cuts (must match decayProbPerEvent_2body.py)
SEP_MIN = 0.001    # m — minimum track separation (1 mm)
SEP_MAX = 1.0      # m — maximum track separation
P_CUT = 0.600      # GeV/c — minimum daughter momentum

# Electron mass for 2-body acceptance kinematics
M_ELECTRON = 0.000511  # GeV/c²

# FONLL mass range: channels effective for m_N below ~5 GeV
FONLL_MASS_MAX = 5.0  # GeV

# Flavors
FLAVORS = ["Ue", "Umu", "Utau"]
