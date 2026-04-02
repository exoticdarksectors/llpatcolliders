#!/usr/bin/env python3
"""
generate_br_table.py  --  regenerate dp_brs_deliver.csv from a local DeLiVeR clone.

Usage (from dark_photon/generator/):
    conda run -n llpatcolliders python generate_br_table.py \\
        --deliver-path /path/to/DeLiVeR/clone \\
        --outfile br_tables/dp_brs_deliver.csv

DeLiVeR clone:  git clone https://github.com/preimitz/DeLiVeR
Ref:            arXiv:2201.01788

scipy >= 1.14 compat: scipy.integrate.quadrature was removed.
This script patches the DeLiVeR form-factor files in-place (only in /tmp copy)
before importing, so the original clone is NOT modified.
"""

import argparse, os, re, shutil, sys, tempfile
import numpy as np


def _patch_deliver(src_dir, dst_dir):
    """Copy DeLiVeR src/ to dst_dir and patch quadrature -> quad."""
    shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)
    ff_dir = os.path.join(dst_dir, "src", "form_factors")
    for fname in os.listdir(ff_dir):
        if not fname.endswith(".py"):
            continue
        path = os.path.join(ff_dir, fname)
        txt = open(path).read()
        if "scipy.integrate.quadrature" not in txt:
            continue
        # Step 1: rename quadrature -> quad, map tol= -> epsabs=, maxiter= -> limit=
        txt = txt.replace("scipy.integrate.quadrature", "scipy.integrate.quad")
        txt = re.sub(r",\s*tol=([\de\-\.]+)", r", epsabs=\1", txt)
        txt = re.sub(r",\s*maxiter=(\d+)", r", limit=\1", txt)
        # Step 2: wrap integrand so scalar x from quad is converted to [x]
        # Pattern A: args=scalar  (e.g., args=s)
        txt = re.sub(
            r"scipy\.integrate\.quad\((\w+),([\w]+),([\w]+),(args=\w+), epsabs=([^,]+), limit=(\d+)\)\[0\]",
            r"scipy.integrate.quad(lambda _x,*_a: \1([_x],*_a)[0],\2,\3,\4, epsabs=\5, limit=\6)[0]",
            txt,
        )
        # Pattern B: args=tuple  (e.g., args=(s,mode))
        txt = re.sub(
            r"scipy\.integrate\.quad\((\w+),([\w]+),([\w]+),(args=\([^)]+\)), epsabs=([^,]+), limit=(\d+)\)\[0\]",
            r"scipy.integrate.quad(lambda _x,*_a: \1([_x],*_a)[0],\2,\3,\4, epsabs=\5, limit=\6)[0]",
            txt,
        )
        open(path, "w").write(txt)


def run(deliver_path, outfile, mmax=1.75):
    import matplotlib
    matplotlib.use("Agg")

    # Work in a temp copy so original is untouched
    tmp = tempfile.mkdtemp(prefix="deliver_patched_")
    try:
        _patch_deliver(deliver_path, tmp)
        orig_dir = os.getcwd()
        os.chdir(tmp)
        sys.path.insert(0, tmp)

        import src.vecdecays as vd
        import src.pars as par

        modelDP = vd.Model("DP")
        modelDP.set_charges([
            -par.ge / 3, par.ge * 2.0 / 3,
            -par.ge / 3, par.ge * 2.0 / 3,
            -par.ge / 3, par.ge * 2.0 / 3,
            -par.ge, -par.ge, -par.ge,
            0.0, 0.0, 0.0,
        ])
        modelDP.gQ = 1.0
        modelDP.set_DMtype(DM="No")
        modelDP.set_folders()

        widths = vd.Widths(modelDP)
        widths.calc(mmax=mmax)
        for sub in ["KK_c", "KK_n", "4pi_c", "4pi_n"]:
            widths.calc_single(sub)
        widths.save()

        brs = vd.Branching(widths)
        brs.calc()
        for sub in ["KK_c", "KK_n", "4pi_c", "4pi_n"]:
            brs.calc_single(sub)
        brs.save()

        masses = np.array(widths.masses)

        os.chdir(orig_dir)
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)

        MINOR_HAD = [
            "EtaOmega", "EtaPhi", "PhiPi", "OmegaPion",
            "EtaPiPi", "EtaPrimePiPi", "6pi",
            "KKpi", "KKpipi", "PhiPiPi", "OmPiPi", "nnbar",
        ]

        header = (
            "# Dark photon (A') branching ratios from DeLiVeR (VMD + PDG R-ratio)\n"
            "# Ref: arXiv:2201.01788, https://github.com/preimitz/DeLiVeR\n"
            "# Model: kinetic mixing portal (EM charges x e per SM fermion), gQ=epsilon\n"
            f"# mmax={mmax} GeV; scipy>=1.14 compat patch applied\n"
            "# Columns: mass_GeV,ee,mumu,tau,pipi,pi3,pi4c,pi4n,KKc,KKn,"
            "PiGamma,EtaGamma,ppbar,had_other,BRqcd\n"
        )

        rows = []
        for i, m in enumerate(masses):
            ee  = brs.BRslep["elec"][i]
            mu  = brs.BRslep["muon"][i]
            tt  = brs.BRslep["tau"][i]
            p2  = brs.BRshad["2pi"][i]
            p3  = brs.BRshad["3pi"][i]
            p4c = brs.BRsinglehad.get("4pi_c", [0] * len(masses))[i]
            p4n = brs.BRsinglehad.get("4pi_n", [0] * len(masses))[i]
            KKc = brs.BRsinglehad.get("KK_c",  [0] * len(masses))[i]
            KKn = brs.BRsinglehad.get("KK_n",  [0] * len(masses))[i]
            PG  = brs.BRshad["PiGamma"][i]
            EG  = brs.BRshad["EtaGamma"][i]
            pp  = brs.BRshad["ppbar"][i]
            oth = sum(brs.BRshad[c][i] for c in MINOR_HAD)
            qcd = brs.BRqcd[i]
            rows.append([m, ee, mu, tt, p2, p3, p4c, p4n, KKc, KKn, PG, EG, pp, oth, qcd])

        with open(outfile, "w") as f:
            f.write(header)
            f.write("mass_GeV,ee,mumu,tau,pipi,pi3,pi4c,pi4n,KKc,KKn,"
                    "PiGamma,EtaGamma,ppbar,had_other,BRqcd\n")
            for row in rows:
                f.write(",".join(f"{v:.8f}" for v in row) + "\n")

        print(f"Saved {len(rows)}-row BR table to {outfile}")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--deliver-path", required=True,
                   help="Path to DeLiVeR git clone")
    p.add_argument("--outfile", default="br_tables/dp_brs_deliver.csv")
    p.add_argument("--mmax", type=float, default=1.75,
                   help="Max mass in GeV (default 1.75; DeLiVeR is VMD-based below ~1.7 GeV)")
    args = p.parse_args()
    run(args.deliver_path, args.outfile, args.mmax)


if __name__ == "__main__":
    main()
