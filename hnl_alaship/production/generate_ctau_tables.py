#!/usr/bin/env python3
"""
production/generate_ctau_tables.py

Compute HNL proper decay length (cτ) vs mass for each flavor at U²=1.

Uses FairShip hnl.py (via the analysis layer) for decay widths → lifetime → cτ.
This replaces HNLCalc for the decay side (which has a double-counting bug above 1 GeV).

Output: output/ctau/ctau_{Ue,Umu,Utau}.dat
Format: mN_GeV  ctau_meters_at_Usq1
"""

import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID
from analysis.fairship_decay import ctau_u2eq1

OUTPUT_DIR = PROJECT_ROOT / "output" / "ctau"


def compute_ctau_table(flavor, masses):
    ctau_vals = [ctau_u2eq1(m, flavor=flavor) for m in masses]
    return np.column_stack([masses, ctau_vals])


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate HNL cτ tables (FairShip decay backend)")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    args = parser.parse_args()

    masses = args.masses if args.masses else MASS_GRID
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for flavor in args.flavor:
        print(f"Computing cτ for {flavor} ({len(masses)} mass points) [FairShip backend]...")
        table = compute_ctau_table(flavor, masses)
        out_path = OUTPUT_DIR / f"ctau_{flavor}.dat"
        np.savetxt(out_path, table, header="mN_GeV  ctau_meters_at_Usq1",
                   fmt=["%.4f", "%.6e"])
        print(f"  → {out_path}")
        for i in [0, len(masses)//4, len(masses)//2, -1]:
            print(f"    m_N = {table[i,0]:.2f} GeV → cτ = {table[i,1]:.3e} m")

    print("Done.")


if __name__ == "__main__":
    main()
