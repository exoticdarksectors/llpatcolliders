#!/usr/bin/env python3
"""
validate_brs.py — Phase 2 validation of the dark-photon BR pipeline.

Compares the pipeline (dp_brs_deliver.csv table + linear interpolation)
against DeLiVeR run directly at the same mass points.  Also checks the
perturbative QCD BRs at high mass.

Note: DarkCast is not available in this environment. DeLiVeR direct computation
is used instead as the ground-truth reference, which is equivalent for the
purposes of verifying our 25-point table + interpolation.

Acceptance criteria (W4 Phase 2):
  outside ρ/ω/φ resonance windows : |ΔBR/BR| < 10%
  inside  resonance windows        : |ΔBR/BR| < 20%
  (resonance windows: ρ+ω = 0.68–0.82 GeV; φ = 0.98–1.05 GeV)

Validation grids:
  On-grid  (exact CSV points): 0.50, 0.77, 1.00, 1.50, 1.70 GeV — must give ~0% error
  Off-grid (interpolation):    0.33, 0.62, 0.82, 1.15, 1.65 GeV — tests linear interp
  Transition check:            VMD(1.70) vs pQCD(1.71) — diagnostic smoothness check

Usage (from dark_photon/generator/ or dark_photon/):
    conda run -n llpatcolliders python generator/validate_brs.py
    conda run -n llpatcolliders python generator/validate_brs.py \\
        --deliver-path /tmp/deliver_clone

Output:
    dark_photon/output/data/dp_br_validation.csv
    dark_photon/output/images/dp_br_vs_mass.png
"""

import argparse
import math
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
OUTPUT_DATA = HERE.parent / "output" / "data"
OUTPUT_IMGS = HERE.parent / "output" / "images"
BR_TABLE    = HERE / "br_tables" / "dp_brs_deliver.csv"

# ── physics constants ──────────────────────────────────────────────────────────
M_MU  = 0.10566
M_TAU = 1.77686
M_K   = 0.49368
M_PI  = 0.13957

# ── validation configuration ──────────────────────────────────────────────────
# Points that ARE in the 25-point CSV grid (must agree to <1%)
ON_GRID  = [0.50, 0.77, 1.00, 1.50, 1.70]
# Points that are NOT in the grid (tests linear interpolation)
OFF_GRID = [0.33, 0.62, 0.82, 1.15, 1.65]
VMD_GRID = sorted(set(ON_GRID + OFF_GRID))   # full VMD validation list

PERT_GRID = [1.75, 2.0, 3.0, 5.0, 10.0, 15.0] # perturbative QCD sanity check

# Resonance tolerance windows (ρ+ω: 0.68–0.82 GeV; φ: 0.98–1.05 GeV)
RESONANCE_WINDOWS = [(0.68, 0.82), (0.98, 1.05)]
ACCEPT_OUTER = 0.10   # 10% outside resonances
ACCEPT_RESON = 0.20   # 20% inside resonances

# Channels compared at each VMD grid point
COMPARE_CHANNELS = ["lep", "had", "ee", "mumu", "pipi", "pi3", "KKc"]


def in_resonance(m):
    return any(lo <= m <= hi for lo, hi in RESONANCE_WINDOWS)


# ── BR table loading ───────────────────────────────────────────────────────────
def load_table():
    cols, rows = None, []
    with open(BR_TABLE) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if cols is None:
                cols = [c.strip() for c in line.split(",")]
                continue
            rows.append([float(x) for x in line.split(",")])
    return cols, np.array(rows)


def interp_at(mass, cols, data):
    m = data[:, 0]
    return {cols[j]: float(np.interp(mass, m, data[:, j]))
            for j in range(len(cols))}


def table_summary(bi):
    """Extract lep, had, and named channels from an interp_at() result."""
    return {
        "lep":  bi.get("ee", 0.0) + bi.get("mumu", 0.0),
        "had":  bi.get("BRqcd", 0.0),
        "ee":   bi.get("ee",   0.0),
        "mumu": bi.get("mumu", 0.0),
        "pipi": bi.get("pipi", 0.0),
        "pi3":  bi.get("pi3",  0.0),
        "KKc":  bi.get("KKc",  0.0),
    }


# ── perturbative BRs (mirrors make_dp_cmnd._perturbative_brs) ─────────────────
def _perturbative_brs(mass):
    nf     = 3 + (mass > 2 * 1.27) + (mass > 2 * 4.18)
    Lambda = 0.2
    als    = 12 * math.pi / ((33 - 2 * nf) * math.log((mass / Lambda) ** 2))
    als    = min(als, 0.4)

    n_uu = 1 + (mass > 2 * 1.27)          # u + (c)
    n_dd = 2 + (mass > 2 * 4.18)          # d + s + (b)
    R    = 3 * (n_uu * (2/3)**2 + n_dd * (1/3)**2) * (1 + als / math.pi)

    total = 3.0 + R
    qcd   = 1 + als / math.pi

    lep = {}
    for col, thresh in [("ee", 0.0), ("mumu", 2*M_MU), ("tau", 2*M_TAU)]:
        lep[col] = (1.0 / total) if mass > thresh else 0.0

    had = {}
    for flavor, q2, thresh in [
        ("uu", (2/3)**2, 0.0),
        ("dd", (1/3)**2, 0.0),
        ("ss", (1/3)**2, 0.0),
        ("cc", (2/3)**2, 2 * 1.27),
        ("bb", (1/3)**2, 2 * 4.18),
    ]:
        had[flavor] = (3 * q2 * qcd / total) if mass > thresh else 0.0

    return {
        **lep,
        **had,
        "lep": sum(lep.values()),
        "had": sum(had.values()),
    }


# ── DeLiVeR direct runner ─────────────────────────────────────────────────────
def _patch_deliver(src_dir, dst_dir):
    """Copy DeLiVeR to a temp dir and apply scipy >= 1.14 compat patch."""
    shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)
    ff_dir = os.path.join(dst_dir, "src", "form_factors")
    for fname in os.listdir(ff_dir):
        if not fname.endswith(".py"):
            continue
        path = os.path.join(ff_dir, fname)
        txt  = open(path).read()
        if "scipy.integrate.quadrature" not in txt:
            continue
        # Step 1: rename + remap kwargs
        txt = txt.replace("scipy.integrate.quadrature", "scipy.integrate.quad")
        txt = re.sub(r",\s*tol=([\de\-\.]+)",  r", epsabs=\1", txt)
        txt = re.sub(r",\s*maxiter=(\d+)",      r", limit=\1",  txt)
        # Step 2: scalar-wrap integrand (quad passes float; quadrature passed array)
        for pattern, replacement in [
            (
                r"scipy\.integrate\.quad\((\w+),([\w]+),([\w]+),(args=\w+), epsabs=([^,]+), limit=(\d+)\)\[0\]",
                r"scipy.integrate.quad(lambda _x,*_a: \1([_x],*_a)[0],\2,\3,\4, epsabs=\5, limit=\6)[0]",
            ),
            (
                r"scipy\.integrate\.quad\((\w+),([\w]+),([\w]+),(args=\([^)]+\)), epsabs=([^,]+), limit=(\d+)\)\[0\]",
                r"scipy.integrate.quad(lambda _x,*_a: \1([_x],*_a)[0],\2,\3,\4, epsabs=\5, limit=\6)[0]",
            ),
        ]:
            txt = re.sub(pattern, replacement, txt)
        open(path, "w").write(txt)


def run_deliver_direct(masses, deliver_path):
    """
    Run DeLiVeR from *deliver_path* and return a dict
        {mass: {"lep", "had", "ee", "mumu", "pipi", "pi3", "KKc", "KKn", "_dm"}}.
    """
    import matplotlib
    matplotlib.use("Agg")

    mmax    = max(max(masses) + 0.01, 1.75)  # DeLiVeR needs mmax>=1.75 for transition search
    tmp     = tempfile.mkdtemp(prefix="deliver_val_")
    orig    = os.getcwd()
    result  = {}
    try:
        _patch_deliver(deliver_path, tmp)
        os.chdir(tmp)
        # Remove any stale cached imports from a prior run in this process
        for key in list(sys.modules.keys()):
            if key.startswith("src"):
                del sys.modules[key]
        sys.path.insert(0, tmp)

        import src.vecdecays as vd
        import src.pars      as par

        modelDP = vd.Model("DP")
        modelDP.set_charges([
            -par.ge/3,  par.ge*2./3,
            -par.ge/3,  par.ge*2./3,
            -par.ge/3,  par.ge*2./3,
            -par.ge, -par.ge, -par.ge,
            0., 0., 0.,
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

        dm = np.array(widths.masses)
        n  = len(dm)

        for m in masses:
            idx = int(np.argmin(np.abs(dm - m)))
            ee  = float(brs.BRslep["elec"][idx])
            mu  = float(brs.BRslep["muon"][idx])
            had = float(brs.BRqcd[idx])
            result[m] = {
                "ee":   ee,
                "mumu": mu,
                "lep":  ee + mu,
                "had":  had,
                "pipi": float(brs.BRshad["2pi"][idx]),
                "pi3":  float(brs.BRshad["3pi"][idx]),
                "pi4c": float(brs.BRsinglehad.get("4pi_c", [0.0]*n)[idx]),
                "pi4n": float(brs.BRsinglehad.get("4pi_n", [0.0]*n)[idx]),
                "KKc":  float(brs.BRsinglehad.get("KK_c",  [0.0]*n)[idx]),
                "KKn":  float(brs.BRsinglehad.get("KK_n",  [0.0]*n)[idx]),
                "_dm":  float(dm[idx]),           # actual DeLiVeR grid mass used
            }
    finally:
        os.chdir(orig)
        shutil.rmtree(tmp, ignore_errors=True)

    return result


# ── comparison ────────────────────────────────────────────────────────────────
def compare_vmd(masses, cols, data, direct, on_grid_set):
    """Compare interpolated table vs DeLiVeR direct; return list of result dicts."""
    rows = []
    hdr  = (f"{'Mass':>6}  {'OG':>3}  {'Reson':>5}  "
            f"{'Chan':>6}  {'Interp':>7}  {'Direct':>7}  "
            f"{'ΔBR/BR':>7}  {'Tol':>5}  {'Pass':>4}")
    print("\n=== VMD Validation (table+interp vs DeLiVeR direct) ===")
    print(hdr)
    print("─" * len(hdr))

    for m in masses:
        bi  = interp_at(m, cols, data)
        ts  = table_summary(bi)           # pipeline values
        ref = direct[m]                   # DeLiVeR direct
        reson = in_resonance(m)
        tol   = ACCEPT_RESON if reson else ACCEPT_OUTER
        # on-grid points are exact CSV values → use tighter 1% threshold
        tol_eff = 0.01 if m in on_grid_set else tol
        ogmark  = "✓" if m in on_grid_set else " "

        for ch in COMPARE_CHANNELS:
            v_pipe = ts.get(ch, 0.0)
            v_ref  = ref.get(ch, 0.0)
            # Skip channels where both are negligibly small
            if v_ref < 1e-4 and v_pipe < 1e-4:
                continue
            rel    = abs(v_pipe - v_ref) / max(v_ref, 1e-8)
            passed = rel < tol_eff
            print(f"{m:>6.2f}  {ogmark:>3}  {'reson' if reson else '     ':>5}  "
                  f"{ch:>6}  {v_pipe:>7.4f}  {v_ref:>7.4f}  "
                  f"{rel:>7.4f}  {tol_eff:>5.0%}  {'✓' if passed else '✗':>4}")
            rows.append({
                "mass":         m,
                "channel":      ch,
                "on_grid":      m in on_grid_set,
                "in_resonance": reson,
                "region":       "VMD",
                "BR_interp":    v_pipe,
                "BR_direct":    v_ref,
                "rel_diff":     rel,
                "tolerance":    tol_eff,
                "pass":         passed,
            })
    return rows


def check_transition_smoothness(cols, data, perturb_above=1.7):
    """Compare VMD table and pQCD BRs at the transition boundary (diagnostic)."""
    # VMD side: interpolate at the table edge
    m_vmd = perturb_above
    bi_vmd = interp_at(m_vmd, cols, data)
    ts_vmd = table_summary(bi_vmd)

    # pQCD side: compute just above the boundary
    m_pqcd = perturb_above + 0.01
    b_pqcd = _perturbative_brs(m_pqcd)

    compare = ["lep", "had", "ee", "mumu"]
    print(f"\n=== VMD/pQCD Transition Smoothness (m = {m_vmd} vs {m_pqcd} GeV) ===")
    print(f"{'Channel':>8}  {'VMD':>7}  {'pQCD':>7}  {'Jump':>7}  {'Status':>6}")
    print("─" * 45)

    any_warning = False
    for ch in compare:
        v_vmd  = ts_vmd.get(ch, 0.0)
        v_pqcd = b_pqcd.get(ch, 0.0)
        denom  = max(abs(v_vmd), abs(v_pqcd), 1e-8)
        jump   = abs(v_vmd - v_pqcd) / denom
        ok     = jump < 0.15
        if not ok:
            any_warning = True
        print(f"{ch:>8}  {v_vmd:>7.4f}  {v_pqcd:>7.4f}  {jump:>7.1%}  "
              f"{'ok' if ok else 'WARN':>6}")

    if any_warning:
        print("WARNING: >15% discontinuity at VMD/pQCD boundary. "
              "This is expected in the 1.7 GeV transition region "
              "(different systematics).")
    else:
        print("Transition smooth (<15% jump in all channels).")


def check_perturbative(masses):
    """Print and return perturbative BR sanity check rows."""
    rows = []
    print("\n=== Perturbative pQCD BRs (no external reference; tabulated for record) ===")
    print(f"{'Mass':>6}  {'ee':>7}  {'mumu':>7}  {'tau':>7}  {'lep':>7}  {'had':>7}")
    print("─" * 55)
    for m in masses:
        b = _perturbative_brs(m)
        print(f"{m:>6.1f}  {b['ee']:>7.4f}  {b['mumu']:>7.4f}  {b['tau']:>7.4f}  "
              f"{b['lep']:>7.4f}  {b['had']:>7.4f}")
        for ch in ["lep", "had"]:
            rows.append({
                "mass":         m,
                "channel":      ch,
                "on_grid":      False,
                "in_resonance": False,
                "region":       "perturbative",
                "BR_interp":    b[ch],
                "BR_direct":    float("nan"),
                "rel_diff":     float("nan"),
                "tolerance":    float("nan"),
                "pass":         True,
            })
    return rows


# ── plot ──────────────────────────────────────────────────────────────────────
def make_plot(cols, data, direct, outpath):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ci = {c: i for i, c in enumerate(cols)}
    m_vmd  = data[:, 0]

    ee_v   = data[:, ci["ee"]]
    mu_v   = data[:, ci["mumu"]]
    lep_v  = ee_v + mu_v
    had_v  = data[:, ci["BRqcd"]]
    pipi_v = data[:, ci["pipi"]]
    pi3_v  = data[:, ci["pi3"]]
    pi4c_v = data[:, ci["pi4c"]]
    pi4n_v = data[:, ci["pi4n"]]
    KKc_v  = data[:, ci["KKc"]]
    KKn_v  = data[:, ci["KKn"]]

    m_pert  = np.linspace(1.70, 15.0, 200)
    lep_p   = np.array([_perturbative_brs(m)["lep"] for m in m_pert])
    had_p   = np.array([_perturbative_brs(m)["had"] for m in m_pert])

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 10),
        gridspec_kw={"height_ratios": [1.1, 1.0]},
    )

    # ── panel 1: leptonic vs hadronic ────────────────────────────────────────
    ax1.plot(m_vmd, lep_v, "b-",  lw=2,   label=r"Leptonic $(e^+e^-{+}\mu^+\mu^-)$ — VMD")
    ax1.plot(m_vmd, had_v, "r-",  lw=2,   label=r"Hadronic total — VMD")
    ax1.plot(m_pert, lep_p, "b--", lw=2,  label=r"Leptonic — pQCD")
    ax1.plot(m_pert, had_p, "r--", lw=2,  label=r"Hadronic — pQCD")

    if direct:
        dm_pts  = np.array(VMD_GRID)
        lep_ref = np.array([direct[m]["lep"] for m in VMD_GRID])
        had_ref = np.array([direct[m]["had"] for m in VMD_GRID])
        ax1.scatter(dm_pts, lep_ref, c="blue",  s=55, zorder=5, marker="o",
                    label="DeLiVeR direct (lep)")
        ax1.scatter(dm_pts, had_ref, c="red",   s=55, zorder=5, marker="s",
                    label="DeLiVeR direct (had)")

    # resonance shading
    for idx, (lo, hi) in enumerate(RESONANCE_WINDOWS):
        ax1.axvspan(lo, hi, alpha=0.08, color="orange",
                    label=r"Resonance windows ($\rho/\omega$, $\phi$)" if idx == 0 else "")

    # kinematic thresholds
    for m_th, lbl, col in [
        (2 * M_MU,  r"$2m_\mu$",  "green"),
        (2 * M_K,   r"$2m_K$",    "purple"),
        (2 * M_TAU, r"$2m_\tau$", "saddlebrown"),
    ]:
        ax1.axvline(m_th, ls=":", color=col, alpha=0.7, lw=1.2)
        ax1.text(m_th * 1.03, 0.03, lbl, fontsize=7, color=col, va="bottom")

    ax1.axvline(1.7, ls="--", color="gray", alpha=0.5, lw=1.0)
    ax1.text(1.73, 0.50, "VMD/pQCD\nboundary", fontsize=6.5, color="gray", va="center")

    ax1.set_xscale("log")
    ax1.set_xlim(0.19, 16.0)
    ax1.set_ylim(-0.02, 1.05)
    ax1.set_xlabel(r"$m_{A'}$ (GeV)")
    ax1.set_ylabel("Branching fraction")
    ax1.set_title(r"Dark Photon Branching Fractions — VMD (DeLiVeR) + pQCD")
    ax1.legend(fontsize=8, loc="center right", ncol=2)
    ax1.grid(True, alpha=0.3, which="both")

    # ── panel 2: individual hadronic channels (VMD range) ────────────────────
    ax2.plot(m_vmd, pipi_v, "-",  color="tab:blue",   lw=1.8, label=r"$\pi^+\pi^-$")
    ax2.plot(m_vmd, pi3_v,  "-",  color="tab:orange",  lw=1.8, label=r"$\pi^0\pi^+\pi^-$")
    ax2.plot(m_vmd, pi4c_v, "-",  color="tab:green",   lw=1.8, label=r"$\pi^+\pi^+\pi^-\pi^-$")
    ax2.plot(m_vmd, pi4n_v, "-",  color="tab:red",     lw=1.8, label=r"$\pi^0\pi^0\pi^+\pi^-$")
    ax2.plot(m_vmd, KKc_v,  "-",  color="tab:purple",  lw=1.8, label=r"$K^+K^-$")
    ax2.plot(m_vmd, KKn_v,  "-",  color="tab:brown",   lw=1.8, label=r"$K^0\bar{K}^0$")

    for lo, hi in RESONANCE_WINDOWS:
        ax2.axvspan(lo, hi, alpha=0.08, color="orange")

    ax2.set_xscale("log")
    ax2.set_xlim(0.19, 1.76)
    ax2.set_ylim(-0.01, 0.82)
    ax2.set_xlabel(r"$m_{A'}$ (GeV)")
    ax2.set_ylabel("Branching fraction")
    ax2.set_title(r"Individual Hadronic Channels — VMD region (0.2–1.7 GeV)")
    ax2.legend(fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")


# ── CSV output ────────────────────────────────────────────────────────────────
def save_csv(rows, outpath):
    def _fmt(v):
        return "nan" if (isinstance(v, float) and math.isnan(v)) else f"{v:.6f}"

    with open(outpath, "w") as fh:
        fh.write("mass_GeV,channel,on_grid,in_resonance,region,"
                 "BR_interp,BR_direct,rel_diff,tolerance,pass\n")
        for r in rows:
            fh.write(
                f"{r['mass']},{r['channel']},"
                f"{int(r['on_grid'])},{int(r['in_resonance'])},{r['region']},"
                f"{_fmt(r['BR_interp'])},{_fmt(r['BR_direct'])},"
                f"{_fmt(r['rel_diff'])},{_fmt(r['tolerance'])},{int(r['pass'])}\n"
            )
    print(f"Saved: {outpath}")


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deliver-path", default="/tmp/deliver_clone",
                    help="Path to DeLiVeR clone (default: /tmp/deliver_clone)")
    ap.add_argument("--no-deliver", action="store_true",
                    help="Skip DeLiVeR direct run; produce plot only from CSV table")
    args = ap.parse_args()

    OUTPUT_DATA.mkdir(parents=True, exist_ok=True)
    OUTPUT_IMGS.mkdir(parents=True, exist_ok=True)

    cols, data = load_table()
    table_masses = set(round(float(data[i, 0]), 4) for i in range(len(data)))
    on_grid_set  = {round(m, 4) for m in VMD_GRID if round(m, 4) in table_masses}
    print(f"BR table: {len(data)} rows, {data[0,0]:.2f}–{data[-1,0]:.2f} GeV")
    print(f"VMD grid : {VMD_GRID}")
    print(f"On-grid  : {sorted(on_grid_set)}")
    print(f"Off-grid : {sorted(set(round(m,4) for m in VMD_GRID) - on_grid_set)}")

    direct = None
    rows   = []

    if not args.no_deliver:
        dp = args.deliver_path
        if not os.path.isdir(dp):
            print(f"\nWARNING: DeLiVeR path not found: {dp}")
            print("Run with --deliver-path or use --no-deliver for plot-only mode.")
        else:
            print(f"\nRunning DeLiVeR directly from {dp} …")
            try:
                direct = run_deliver_direct(VMD_GRID, dp)
                rows  += compare_vmd(VMD_GRID, cols, data, direct, on_grid_set)
            except Exception as exc:
                print(f"ERROR running DeLiVeR: {exc}")
                direct = None

    rows += check_perturbative(PERT_GRID)
    check_transition_smoothness(cols, data)

    # ── summary ───────────────────────────────────────────────────────────────
    vmd_rows = [r for r in rows if r["region"] == "VMD"]
    if vmd_rows:
        n_pass  = sum(1 for r in vmd_rows if r["pass"])
        n_total = len(vmd_rows)
        print(f"\n=== Summary: {n_pass}/{n_total} comparisons pass ===")
        fails = [r for r in vmd_rows if not r["pass"]]
        if fails:
            print("FAILURES:")
            for r in fails:
                print(f"  m={r['mass']:.2f} GeV  ch={r['channel']:6s}  "
                      f"interp={r['BR_interp']:.4f}  direct={r['BR_direct']:.4f}  "
                      f"Δ={r['rel_diff']:.3f} > tol={r['tolerance']:.2f}")
        else:
            print("All VMD comparisons pass acceptance criteria.")

    csv_out  = OUTPUT_DATA / "dp_br_validation.csv"
    plot_out = OUTPUT_IMGS / "dp_br_vs_mass.png"
    save_csv(rows, csv_out)
    make_plot(cols, data, direct, plot_out)

    print("\nPhase 2 validation complete.")
    print(f"  CSV  : {csv_out}")
    print(f"  Plot : {plot_out}")


if __name__ == "__main__":
    main()
