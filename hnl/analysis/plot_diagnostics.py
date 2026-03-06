#!/usr/bin/env python3
"""
analysis/plot_diagnostics.py

Diagnostic investigation plots for understanding exclusion contour features.

Produces:
  1. BR_vis vs m_N (all 3 flavors)
  2. cτ vs m_N (all 3 flavors)
  3. Event count & hit count vs m_N (from sensitivity CSV)
  4. Production weight sum vs m_N per channel
  5. N_signal vs U² at selected mass points

Usage:
    python analysis/plot_diagnostics.py
    python analysis/plot_diagnostics.py --flavor Umu
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from analysis.constants import FLAVORS, FONLL_MASS_MAX, N_THRESHOLD

CTAU_DIR = PROJECT_ROOT / "output" / "ctau"
CSV_DIR = PROJECT_ROOT / "output" / "llp_4vectors"
ANALYSIS_DIR = PROJECT_ROOT / "output" / "analysis"
CHANNELS = ["Bmeson", "Dmeson", "Bc", "tau", "WZ"]

_FLAVOR_COLOR = {"Ue": "#e41a1c", "Umu": "#377eb8", "Utau": "#4daf4a"}
_FLAVOR_LABEL = {
    "Ue":   r"$|U_e|^2$",
    "Umu":  r"$|U_\mu|^2$",
    "Utau": r"$|U_\tau|^2$",
}
_CHAN_STYLE = {
    "Bmeson": {"color": "#e41a1c", "marker": "o", "label": "B meson"},
    "Dmeson": {"color": "#377eb8", "marker": "s", "label": "D meson"},
    "Bc":     {"color": "#4daf4a", "marker": "^", "label": r"$B_c$"},
    "tau":    {"color": "#984ea3", "marker": "D", "label": r"$\tau$"},
    "WZ":     {"color": "#ff7f00", "marker": "v", "label": "W/Z"},
}


def _load_table(path):
    """Load a 2-column .dat file (mass, value), skipping comments."""
    masses, values = [], []
    for line in path.read_text().strip().split("\n"):
        if line.startswith("#"):
            continue
        parts = line.split()
        masses.append(float(parts[0]))
        values.append(float(parts[1]))
    return np.array(masses), np.array(values)


def plot_br_vis(flavors, output_dir):
    """Plot 1: BR_vis vs m_N for each flavor."""
    fig, ax = plt.subplots(figsize=(10, 6))
    for flavor in flavors:
        path = CTAU_DIR / f"br_vis_{flavor}.dat"
        if not path.exists():
            continue
        masses, br_vis = _load_table(path)
        mask = masses <= FONLL_MASS_MAX
        ax.plot(masses[mask], br_vis[mask], "-o", markersize=2,
                color=_FLAVOR_COLOR[flavor],
                label=f"{flavor}: {_FLAVOR_LABEL[flavor]}", linewidth=1.5)
    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel(r"$\mathrm{BR}_\mathrm{vis}$", fontsize=13)
    ax.set_title("Visible branching ratio (fraction with charged tracks)", fontsize=13)
    ax.set_xlim([0.15, 5.2])
    ax.set_ylim([0, 1.0])
    ax.axhline(0.5, color="gray", ls=":", alpha=0.5)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    out = output_dir / "diag_br_vis.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_ctau(flavors, output_dir):
    """Plot 2: cτ(U²=1) vs m_N for each flavor."""
    fig, ax = plt.subplots(figsize=(10, 6))
    for flavor in flavors:
        path = CTAU_DIR / f"ctau_{flavor}.dat"
        if not path.exists():
            continue
        masses, ctau = _load_table(path)
        mask = masses <= FONLL_MASS_MAX
        ax.plot(masses[mask], ctau[mask], "-", markersize=2,
                color=_FLAVOR_COLOR[flavor],
                label=f"{flavor}", linewidth=1.5)

    # Mark HNLCalc hadron/quark transition at ~1 GeV
    ax.axvline(1.0, color="gray", ls="--", alpha=0.5, label="hadron/quark transition")

    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel(r"$c\tau$ at $U^2=1$ [m]", fontsize=13)
    ax.set_title(r"HNL proper decay length at $U^2 = 1$", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.15, 5.2])
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.2)
    out = output_dir / "diag_ctau.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_event_counts(flavors, output_dir):
    """Plot 3: Event count and hit count vs m_N from sensitivity CSV."""
    csv_path = ANALYSIS_DIR / "gargoyle_hnl_sensitivity.csv"
    if not csv_path.exists():
        print(f"  Skipping event counts: {csv_path} not found")
        return
    df = pd.read_csv(csv_path)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: event counts
    ax = axes[0]
    for flavor in flavors:
        sel = df[df["flavor"] == flavor].sort_values("mass_GeV")
        ax.plot(sel["mass_GeV"], sel["n_events"], "-o", markersize=3,
                color=_FLAVOR_COLOR[flavor], label=flavor, linewidth=1)
    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel("Number of events", fontsize=13)
    ax.set_title("Combined 4-vector pool size", fontsize=13)
    ax.set_xscale("log")
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.2)

    # Right: hit counts (geometry acceptance)
    ax = axes[1]
    for flavor in flavors:
        sel = df[df["flavor"] == flavor].sort_values("mass_GeV")
        ax.plot(sel["mass_GeV"], sel["n_hits"], "-o", markersize=3,
                color=_FLAVOR_COLOR[flavor], label=flavor, linewidth=1)
    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel("Number of detector hits", fontsize=13)
    ax.set_title("Events hitting GARGOYLE", fontsize=13)
    ax.set_xscale("log")
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out = output_dir / "diag_event_counts.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_weight_sums(flavors, output_dir):
    """Plot 4: Production cross-section (weight sum) vs m_N per channel."""
    fonll_masses = [m for m in MASS_GRID if m <= FONLL_MASS_MAX]

    for flavor in flavors:
        fig, ax = plt.subplots(figsize=(10, 6))
        channel_data = {}

        for chan in CHANNELS:
            chan_dir = CSV_DIR / flavor / chan
            if not chan_dir.exists():
                continue
            masses_found = []
            wt_sums = []
            for mass in fonll_masses:
                mlabel = format_mass_for_filename(mass)
                csv_path = chan_dir / f"mN_{mlabel}.csv"
                if not csv_path.exists():
                    continue
                try:
                    data = np.loadtxt(csv_path, delimiter=",")
                    if data.ndim == 1:
                        data = data.reshape(1, -1)
                    wt_sum = data[:, 0].sum()
                    masses_found.append(mass)
                    wt_sums.append(wt_sum)
                except Exception:
                    continue
            if masses_found:
                channel_data[chan] = (np.array(masses_found), np.array(wt_sums))
                style = _CHAN_STYLE.get(chan, {})
                ax.plot(masses_found, wt_sums, "-",
                        marker=style.get("marker", "o"), markersize=3,
                        color=style.get("color", "gray"),
                        label=style.get("label", chan), linewidth=1.2)

        # Also plot combined
        comb_dir = CSV_DIR / flavor / "combined"
        if comb_dir.exists():
            masses_found = []
            wt_sums = []
            for mass in fonll_masses:
                mlabel = format_mass_for_filename(mass)
                csv_path = comb_dir / f"mN_{mlabel}.csv"
                if not csv_path.exists():
                    continue
                try:
                    data = np.loadtxt(csv_path, delimiter=",")
                    if data.ndim == 1:
                        data = data.reshape(1, -1)
                    wt_sum = data[:, 0].sum()
                    masses_found.append(mass)
                    wt_sums.append(wt_sum)
                except Exception:
                    continue
            if masses_found:
                ax.plot(masses_found, wt_sums, "k-", linewidth=2,
                        label="Combined", zorder=10)

        ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
        ax.set_ylabel(r"$\sum w_i$ [pb] (at $U^2=1$)", fontsize=13)
        ax.set_title(f"Production weight sum — {flavor}", fontsize=13)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim([0.15, 5.2])
        ax.legend(fontsize=10, loc="upper right")
        ax.grid(True, which="both", alpha=0.2)
        out = output_dir / f"diag_weight_sum_{flavor}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out}")


def plot_sensitivity_details(flavors, output_dir):
    """Plot 5: Exclusion band boundaries + peak N_signal vs m_N."""
    csv_path = ANALYSIS_DIR / "gargoyle_hnl_sensitivity.csv"
    if not csv_path.exists():
        print(f"  Skipping sensitivity details: {csv_path} not found")
        return
    df = pd.read_csv(csv_path)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: exclusion band (u2_min, u2_max) vs mass
    ax = axes[0]
    for flavor in flavors:
        sel = df[(df["flavor"] == flavor) & (df["has_sensitivity"] == True)].sort_values("mass_GeV")
        if len(sel) == 0:
            continue
        ax.fill_between(sel["mass_GeV"], sel["u2_min"], sel["u2_max"],
                        alpha=0.2, color=_FLAVOR_COLOR[flavor])
        ax.plot(sel["mass_GeV"], sel["u2_min"], "-",
                color=_FLAVOR_COLOR[flavor], linewidth=1, label=f"{flavor} lower")
        ax.plot(sel["mass_GeV"], sel["u2_max"], "--",
                color=_FLAVOR_COLOR[flavor], linewidth=1, label=f"{flavor} upper")
    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel(r"$U^2$", fontsize=13)
    ax.set_title("Exclusion band boundaries (raw, unsmoothed)", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.15, 5.2])
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, which="both", alpha=0.2)

    # Right: peak N_signal vs mass
    ax = axes[1]
    for flavor in flavors:
        sel = df[df["flavor"] == flavor].sort_values("mass_GeV")
        ax.plot(sel["mass_GeV"], sel["peak_N"], "-o", markersize=2,
                color=_FLAVOR_COLOR[flavor], label=flavor, linewidth=1)
    ax.axhline(N_THRESHOLD, color="gray", ls="--", alpha=0.7,
                label=rf"$N_{{\mathrm{{thr}}}}={N_THRESHOLD:.0f}$")
    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel(r"Peak $N_\mathrm{signal}$", fontsize=13)
    ax.set_title("Peak signal yield (at optimal U²)", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.15, 5.2])
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.2)

    fig.tight_layout()
    out = output_dir / "diag_sensitivity_details.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_br_vis_comparison(output_dir):
    """Bonus: Compare our BR_vis with MATHUSLA's BrNtoinvistoMATHUSLA for Umu."""
    mathusla_path = (PROJECT_ROOT / "vendored" /
                     "MATHUSLA_LLPfiles_RHN_Umu" / "RHNbrUmu_exclusive.dat")
    our_path = CTAU_DIR / "br_vis_Umu.dat"

    if not mathusla_path.exists() or not our_path.exists():
        print("  Skipping BR_vis comparison: missing MATHUSLA or our data")
        return

    # Load our BR_vis
    our_m, our_br = _load_table(our_path)

    # Load MATHUSLA data
    # Columns: mN/GeV  BrNtovee  BrNtovmumu  BrNtovemu  BrNto2tracks  BrNtoinvistoMATHUSLA
    math_data = []
    for line in mathusla_path.read_text().strip().split("\n"):
        if line.startswith("#") or line.startswith("mN"):
            continue
        parts = line.split()
        if len(parts) >= 6:
            math_data.append([float(parts[0]), float(parts[5])])
    if not math_data:
        return
    math_data = np.array(math_data)
    math_m = math_data[:, 0]
    math_br_invis = math_data[:, 1]
    math_br_vis = 1.0 - math_br_invis

    fig, ax = plt.subplots(figsize=(10, 6))
    mask = our_m <= 1.0
    ax.plot(our_m[mask], our_br[mask], "b-o", markersize=3, linewidth=1.5,
            label="Our BR_vis (Umu)")
    ax.plot(math_m, math_br_vis, "r--s", markersize=3, linewidth=1.5,
            label="MATHUSLA BR_vis (Umu)")
    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=13)
    ax.set_ylabel(r"$\mathrm{BR}_\mathrm{vis}$", fontsize=13)
    ax.set_title(r"BR$_\mathrm{vis}$ comparison: ours vs MATHUSLA ($U_\mu$ coupling)", fontsize=13)
    ax.set_xlim([0.1, 1.0])
    ax.set_ylim([0, 1.0])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    out = output_dir / "diag_br_vis_comparison.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    parser = argparse.ArgumentParser(description="Diagnostic investigation plots")
    parser.add_argument("--flavor", nargs="+", default=None, choices=FLAVORS)
    args = parser.parse_args()

    flavors = args.flavor or FLAVORS
    output_dir = ANALYSIS_DIR / "diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Diagnostic Investigation Plots ===\n")

    print("1. BR_vis vs m_N")
    plot_br_vis(flavors, output_dir)

    print("\n2. cτ vs m_N")
    plot_ctau(flavors, output_dir)

    print("\n3. Event counts vs m_N")
    plot_event_counts(flavors, output_dir)

    print("\n4. Weight sums vs m_N per channel")
    plot_weight_sums(flavors, output_dir)

    print("\n5. Sensitivity details vs m_N")
    plot_sensitivity_details(flavors, output_dir)

    print("\n6. BR_vis comparison (ours vs MATHUSLA)")
    plot_br_vis_comparison(output_dir)

    print(f"\nAll diagnostics saved to: {output_dir}")


if __name__ == "__main__":
    main()
