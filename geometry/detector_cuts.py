"""
Charged-track acceptance cuts shared across all GARGOYLE portals.

These govern which daughter tracks are counted as reconstructable.
Kept in a standalone module with no heavy imports so any code path
can import without pulling in trimesh / ROOT / etc.
"""

P_CUT = 0.600      # GeV/c — minimum daughter momentum
SEP_MIN = 0.001    # m — minimum track separation (1 mm)
SEP_MAX = 1.0      # m — maximum track separation
