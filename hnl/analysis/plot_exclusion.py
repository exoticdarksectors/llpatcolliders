"""
analysis/plot_exclusion.py

Publication-quality (m_N, U²) exclusion plots for GARGOYLE HNL sensitivity.

Produces a 1×3 panel figure (Ue, Umu, Utau) with GARGOYLE exclusion island
and reference contours from MATHUSLA, ANUBIS, CODEX-b.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.signal import savgol_filter

from analysis.reference_curves import load_all_references
from analysis.constants import FLAVORS, N_THRESHOLD

# Styling for reference experiments
_REF_STYLE = {
    "MATHUSLA":  {"color": "#2166ac", "ls": "--",  "lw": 1.5},
    "ANUBIS":    {"color": "#4dac26", "ls": "-.",  "lw": 1.5},
    "CODEX-b":   {"color": "#e08214", "ls": ":",   "lw": 1.8},
    "SHiP":      {"color": "#7b3294", "ls": "--",  "lw": 1.2},
}

_FLAVOR_LABEL = {
    "Ue":   r"$|U_e|^2$",
    "Umu":  r"$|U_\mu|^2$",
    "Utau": r"$|U_\tau|^2$",
}


def plot_exclusion(results_csv, output_dir=None):
    """
    Create the exclusion money plot from sensitivity results.

    Parameters
    ----------
    results_csv : str or Path
        CSV file with columns: mass_GeV, flavor, u2_min, u2_max, peak_N,
                               peak_u2, has_sensitivity
    output_dir : str or Path, optional
        Where to save plots. Defaults to same directory as results_csv.
    """
    results_csv = Path(results_csv)
    df = pd.read_csv(results_csv)
    if output_dir is None:
        output_dir = results_csv.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_curves = load_all_references()

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    for idx, flavor in enumerate(FLAVORS):
        ax = axes[idx]
        _plot_single_panel(ax, df, flavor, ref_curves)

    fig.tight_layout(w_pad=2.0)
    out_path = output_dir / "gargoyle_hnl_exclusion.pdf"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    out_png = output_dir / "gargoyle_hnl_exclusion.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    print(f"Saved: {out_png}")


def _smooth_log(values, window=7):
    """Smooth values in log-space using Savitzky-Golay filter."""
    if len(values) < window:
        return values
    # Ensure window is odd
    if window % 2 == 0:
        window += 1
    log_vals = np.log10(values)
    smoothed = savgol_filter(log_vals, window, polyorder=3)
    return 10.0 ** smoothed


def _build_island_polygon(sel):
    """
    Build a closed polygon for the exclusion island from the results DataFrame.

    Interpolates closure points at mass boundaries where sensitivity
    transitions, so the island closes smoothly instead of with vertical walls.
    Applies Savitzky-Golay smoothing in log-space to remove jaggedness.

    Returns arrays (mass_poly, u2_poly) for a closed polygon, or None.
    """
    valid = sel[sel["has_sensitivity"] == True].copy()
    if len(valid) == 0:
        return None

    mass = valid["mass_GeV"].values
    u2_lo = valid["u2_min"].values
    u2_hi = valid["u2_max"].values

    # Smooth the boundaries in log-space
    u2_lo = _smooth_log(u2_lo)
    u2_hi = _smooth_log(u2_hi)

    # High-mass closure: find the first non-sensitive point after the last sensitive
    all_sorted = sel.sort_values("mass_GeV")
    last_sens_mass = mass[-1]
    after = all_sorted[
        (all_sorted["mass_GeV"] > last_sens_mass) &
        (all_sorted["has_sensitivity"] == False)
    ]

    # Build polygon: upper contour left→right, then lower contour right→left
    m_upper = list(mass)
    u2_upper = list(u2_hi)
    m_lower = list(mass)
    u2_lower = list(u2_lo)

    # Add high-mass closure tip
    if len(after) > 0:
        row = after.iloc[0]
        peak_N_last = valid.iloc[-1]["peak_N"]
        peak_N_next = row["peak_N"]
        if peak_N_next < peak_N_last and peak_N_last > 0:
            frac = (N_THRESHOLD - peak_N_last) / (peak_N_next - peak_N_last)
            frac = np.clip(frac, 0.01, 0.99)
            m_tip = last_sens_mass + frac * (row["mass_GeV"] - last_sens_mass)
        else:
            m_tip = last_sens_mass
        u2_tip = row["peak_u2"]
        m_upper.append(m_tip)
        u2_upper.append(u2_tip)
        m_lower.append(m_tip)
        u2_lower.append(u2_tip)
    else:
        m_upper.append(mass[-1])
        u2_upper.append(np.sqrt(u2_lo[-1] * u2_hi[-1]))
        m_lower.append(mass[-1])
        u2_lower.append(np.sqrt(u2_lo[-1] * u2_hi[-1]))

    # Low-mass closure tip
    first_sens_mass = mass[0]
    before = all_sorted[
        (all_sorted["mass_GeV"] < first_sens_mass) &
        (all_sorted["has_sensitivity"] == False)
    ]
    if len(before) > 0:
        row = before.iloc[-1]
        peak_N_first = valid.iloc[0]["peak_N"]
        peak_N_prev = row["peak_N"]
        if peak_N_prev < peak_N_first and peak_N_first > 0:
            frac = (N_THRESHOLD - peak_N_first) / (peak_N_prev - peak_N_first)
            frac = np.clip(frac, 0.01, 0.99)
            m_tip = first_sens_mass + frac * (row["mass_GeV"] - first_sens_mass)
        else:
            m_tip = first_sens_mass
        u2_tip = row["peak_u2"]
        m_upper.insert(0, m_tip)
        u2_upper.insert(0, u2_tip)
        m_lower.insert(0, m_tip)
        u2_lower.insert(0, u2_tip)

    # Construct closed polygon: upper L→R then lower R→L
    m_poly = np.array(m_upper + m_lower[::-1])
    u2_poly = np.array(u2_upper + u2_lower[::-1])
    return m_poly, u2_poly


def _plot_single_panel(ax, df, flavor, ref_curves):
    """Plot one flavor panel."""
    sel = df[df["flavor"] == flavor].sort_values("mass_GeV")
    valid = sel[sel["has_sensitivity"] == True]

    # GARGOYLE exclusion island as closed polygon
    if len(valid) > 0:
        poly = _build_island_polygon(sel)
        if poly is not None:
            m_poly, u2_poly = poly
            ax.fill(m_poly, u2_poly, alpha=0.25, color="red",
                    label="GARGOYLE", zorder=5)
            ax.plot(m_poly, u2_poly, "r-", linewidth=1.8, zorder=6)

    # Reference contours
    for exp, curves in ref_curves.items():
        if flavor not in curves:
            continue
        c = curves[flavor]
        style = _REF_STYLE.get(exp, {"color": "gray", "ls": "-", "lw": 1.0})
        ax.fill_between(c["mass"], c["u2_min"], c["u2_max"],
                        alpha=0.08, color=style["color"])
        ax.plot(c["mass"], c["u2_min"],
                color=style["color"], ls=style["ls"], lw=style["lw"],
                label=exp)
        ax.plot(c["mass"], c["u2_max"],
                color=style["color"], ls=style["ls"], lw=style["lw"])

    ax.set_xlabel(r"$m_N$ [GeV]", fontsize=14)
    if flavor == "Ue":
        ax.set_ylabel(r"$|U|^2$", fontsize=14)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([0.15, 6.0])
    ax.set_ylim([1e-12, 1e-1])
    ax.grid(True, which="both", alpha=0.2, linewidth=0.5)
    ax.legend(fontsize=10, loc="upper right")
    ax.set_title(f"HNL {_FLAVOR_LABEL[flavor]}", fontsize=14)


def plot_nsignal_vs_u2(u2_grid, N_grid, mass, flavor, output_dir):
    """
    Diagnostic plot: N_signal vs U² for a single mass point.
    Useful for verifying the island shape.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(u2_grid, N_grid, "b-", linewidth=1.5)
    ax.axhline(N_THRESHOLD, color="r", ls="--", alpha=0.7,
                label=rf"$N_{{\mathrm{{thr}}}} = {N_THRESHOLD:.0f}$")
    ax.set_xlabel(r"$U^2$", fontsize=13)
    ax.set_ylabel(r"$N_\mathrm{signal}$", fontsize=13)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"{flavor}, $m_N$ = {mass:.2f} GeV", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.2)

    out = output_dir / f"nsignal_{flavor}_mN_{mass:.2f}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
