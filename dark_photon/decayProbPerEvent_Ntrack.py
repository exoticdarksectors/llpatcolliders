"""
Dark photon displaced-vertex analysis: ≥N charged tracks.

Signal definition (matching MATHUSLA/ANUBIS/CODEX-b):
  LLP decays inside fiducial volume AND ≥ N_MIN_TRACKS charged
  daughters pass p > P_CUT.

Inputs:
  - LLP CSV:      event,id,pt,eta,phi,momentum,mass
  - Daughter CSV:  event,llp_idx,daughter_pid,pt,eta,phi,px,py,pz,energy,charge
  - Meta JSON:    n_generated, llp_pdg_id

Usage:
  python decayProbPerEvent_Ntrack.py output/dp_heavy_m15.csv \
      --xsec 60000 --lumi 3000 --outdir output/
"""

import os
import sys
import json
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm

from gargoyle_geometry import (
    SPEED_OF_LIGHT,
    calculate_decay_length, cache_geometry,
    mesh_fiducial,
)

# Analysis cuts
P_CUT = 0.600         # GeV/c — minimum daughter momentum
N_MIN_TRACKS = 2      # minimum charged tracks per DV


def count_passing_daughters(dau_df, p_cut=P_CUT):
    """
    For each (event, llp_idx), count how many charged daughters pass p > p_cut.

    Returns dict: (event, llp_idx) -> n_passing
    """
    dau_df = dau_df.copy()
    dau_df['p'] = np.sqrt(dau_df['px']**2 + dau_df['py']**2 + dau_df['pz']**2)
    passing = dau_df[dau_df['p'] >= p_cut]
    counts = passing.groupby(['event', 'llp_idx']).size().to_dict()
    return counts


def analyze_Ntrack(llp_csv, dau_csv, geo_cache, lifetime_range,
                   xsec_fb=60e3, lumi_fb=3000,
                   n_generated=None, p_cut=P_CUT, n_min=N_MIN_TRACKS):
    """
    Scan over lifetimes, computing expected signal with ≥n_min track cut.

    For each LLP that hits the fiducial volume:
      P_detect = P_decay_in_volume × track_accept
    where track_accept = 1 if n_charged_passing >= n_min, else 0.

    The track acceptance is fixed per LLP (depends on decay products, not
    on where the LLP decays). The decay probability is the exponential
    integral over the path through the fiducial volume.
    """
    llp_df = pd.read_csv(llp_csv)
    llp_df.columns = llp_df.columns.str.strip()
    dau_df = pd.read_csv(dau_csv)
    dau_df.columns = dau_df.columns.str.strip()

    n_llp_events = llp_df['event'].nunique()
    if n_generated is None:
        n_generated = n_llp_events
    if n_generated < n_llp_events:
        print(f"WARNING: n_generated ({n_generated}) < n_llp_events "
              f"({n_llp_events}); using n_llp_events.")
        n_generated = n_llp_events

    # Count passing daughters per LLP
    dau_counts = count_passing_daughters(dau_df, p_cut)

    # Build per-LLP track acceptance: 1 if enough tracks, 0 otherwise.
    # LLP rows in llp_df are indexed by (event, llp_idx_within_event).
    n_llp = len(llp_df)
    track_accept = np.zeros(n_llp)
    llp_idx_in_event = np.zeros(n_llp, dtype=int)

    # Assign llp_idx within each event (0, 1, ...)
    event_counts = {}
    for i, row in llp_df.iterrows():
        ev = int(row['event'])
        idx = event_counts.get(ev, 0)
        event_counts[ev] = idx + 1
        llp_idx_in_event[i] = idx
        key = (ev, idx)
        n_pass = dau_counts.get(key, 0)
        track_accept[i] = 1.0 if n_pass >= n_min else 0.0

    hits = geo_cache['hits']
    entry_d = geo_cache['entry_d']
    exit_d = geo_cache['exit_d']
    momentum = geo_cache['momentum']
    mass = geo_cache['mass']

    n_track_pass = int(track_accept.sum())
    n_hits = int(hits.sum())
    n_both = int((track_accept * hits).sum())
    print(f"\nTrack acceptance: {n_track_pass}/{n_llp} LLPs have "
          f"≥{n_min} charged tracks with p>{p_cut*1000:.0f} MeV")
    print(f"Geometric acceptance: {n_hits}/{n_llp} LLPs hit fiducial volume")
    print(f"Both: {n_both}/{n_llp}")

    results = {
        'lifetimes': lifetime_range,
        'mean_decay_prob_Ntrack': [],
        'mean_decay_prob_no_cuts': [],
        'mean_at_least_one': [],
        'mean_at_least_one_no_cuts': [],
        'exclusion': [],
        'exclusion_no_cuts': [],
        'n_generated': n_generated,
        'n_llp_events': n_llp_events,
    }

    for lifetime in tqdm(lifetime_range, desc="Scanning lifetimes"):
        # Per-LLP decay probability (no cuts)
        decay_prob_raw = np.zeros(n_llp)
        decay_prob_Ntrack = np.zeros(n_llp)

        for i in range(n_llp):
            if not hits[i]:
                continue
            dl = calculate_decay_length(momentum[i], mass[i], lifetime)
            path = exit_d[i] - entry_d[i]
            p_decay = np.exp(-entry_d[i] / dl) * (1.0 - np.exp(-path / dl))
            decay_prob_raw[i] = p_decay
            decay_prob_Ntrack[i] = p_decay * track_accept[i]

        # Per-event: P(≥1 LLP detected)
        p_event = np.zeros(n_llp_events)
        p_event_nc = np.zeros(n_llp_events)
        events_unique = llp_df['event'].unique()

        for j, ev in enumerate(events_unique):
            mask = llp_df['event'].values == ev
            probs = decay_prob_Ntrack[mask]
            probs_nc = decay_prob_raw[mask]
            # P(≥1) = 1 - prod(1 - p_i)
            p_event[j] = 1.0 - np.prod(1.0 - probs)
            p_event_nc[j] = 1.0 - np.prod(1.0 - probs_nc)

        mean_p1 = p_event.sum() / n_generated
        mean_p1_nc = p_event_nc.sum() / n_generated

        hit_mask = hits.astype(bool)
        mean_single = decay_prob_Ntrack[hit_mask].mean() if hit_mask.any() else 0
        mean_single_nc = decay_prob_raw[hit_mask].mean() if hit_mask.any() else 0

        denom = mean_p1 * lumi_fb * xsec_fb
        denom_nc = mean_p1_nc * lumi_fb * xsec_fb

        results['mean_decay_prob_Ntrack'].append(mean_single)
        results['mean_decay_prob_no_cuts'].append(mean_single_nc)
        results['mean_at_least_one'].append(mean_p1)
        results['mean_at_least_one_no_cuts'].append(mean_p1_nc)
        results['exclusion'].append(3 / denom if denom > 0 else np.inf)
        results['exclusion_no_cuts'].append(3 / denom_nc if denom_nc > 0 else np.inf)

    return results


# ==================================================================
# Main
# ==================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dark photon displaced vertex analysis (≥N charged tracks)")
    parser.add_argument("csv_file", help="LLP CSV (event,id,pt,eta,phi,momentum,mass)")
    parser.add_argument("--xsec", type=float, default=60e3,
                        help="production cross-section in fb (default: 60e3 for pp→h)")
    parser.add_argument("--lumi", type=float, default=3000,
                        help="integrated luminosity in fb⁻¹ (default: 3000)")
    parser.add_argument("--outdir", default="output",
                        help="output directory for plots and CSVs")
    parser.add_argument("--n-events", type=int, default=None,
                        help="total generated events (auto-read from _meta.json)")
    parser.add_argument("--n-min-tracks", type=int, default=N_MIN_TRACKS,
                        help=f"minimum charged tracks per DV (default: {N_MIN_TRACKS})")
    parser.add_argument("--p-cut", type=float, default=P_CUT,
                        help=f"minimum daughter momentum in GeV (default: {P_CUT})")
    args = parser.parse_args()

    xsec_arg_given = any(
        tok == "--xsec" or tok.startswith("--xsec=")
        for tok in sys.argv[1:]
    )
    if not xsec_arg_given:
        print("WARNING: --xsec not provided; using default 60e3 fb (pp→h).")
    os.makedirs(args.outdir, exist_ok=True)

    sample_csv = args.csv_file
    output_tag = os.path.basename(sample_csv).replace('.csv', '')

    # Daughter CSV: same name with _daughters suffix
    dau_csv = sample_csv.replace('.csv', '_daughters.csv')
    if not os.path.isfile(dau_csv):
        print(f"ERROR: daughter CSV not found: {dau_csv}")
        print("  Re-run production with updated generator.cc to get daughter output.")
        sys.exit(1)

    # Read metadata
    n_generated = args.n_events
    llp_pdg_id = None
    meta_path = sample_csv.replace('.csv', '_meta.json')
    if os.path.isfile(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
        if n_generated is None:
            n_generated = meta['n_generated']
        llp_pdg_id = meta.get('llp_pdg_id')
        print(f"Read n_generated={n_generated}, llp_pdg_id={llp_pdg_id} "
              f"from {meta_path}")
    else:
        print(f"WARNING: no _meta.json found ({meta_path}).")

    origin = [0, 0, 0]
    geo_cache = cache_geometry(sample_csv, mesh_fiducial, origin)

    # Lifetime scan
    print("\n" + "=" * 50)
    print(f"LIFETIME SCAN (≥{args.n_min_tracks} tracks, "
          f"p > {args.p_cut*1000:.0f} MeV)")
    print("=" * 50)
    print(f"  σ = {args.xsec:.3g} fb,  L = {args.lumi:.0f} fb⁻¹")
    print(f"  n_generated = {n_generated}")

    lifetimes = np.logspace(-9.5, -4.5, 20)
    scan = analyze_Ntrack(
        sample_csv, dau_csv, geo_cache, lifetimes,
        xsec_fb=args.xsec, lumi_fb=args.lumi,
        n_generated=n_generated,
        p_cut=args.p_cut, n_min=args.n_min_tracks)

    # === Plotting ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: P(≥1 detected) vs lifetime
    ax1 = axes[0]
    ax1.loglog(lifetimes * 1e9, scan['mean_at_least_one'],
               'r-', linewidth=2, label=f'≥{args.n_min_tracks} tracks')
    ax1.loglog(lifetimes * 1e9, scan['mean_at_least_one_no_cuts'],
               'r--', linewidth=2, alpha=0.5, label='Geometric only')
    ax1.set_xlabel('Lifetime (ns)')
    ax1.set_ylabel('Mean P(≥1 LLP detected)')
    ax1.set_title(f'Detection probability\n'
                  f'(≥{args.n_min_tracks} charged tracks, '
                  f'p > {args.p_cut*1000:.0f} MeV)')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend()

    # Plot 2: Exclusion curve (BR vs cτ)
    ax2 = axes[1]
    ax2.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion'],
               color='blue', linewidth=2,
               label=f"PX56 (≥{args.n_min_tracks} tracks)")
    ax2.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion_no_cuts'],
               color='blue', linewidth=2, linestyle='--', alpha=0.5,
               label="PX56 (geometric only)")
    ax2.set_xlabel(r'$c\tau$ (m)')
    ax2.set_ylabel(r"BR$(h \to A'A')_{\min}$")
    ax2.set_title('Exclusion sensitivity (3 signal events)')
    ax2.grid(True, which="both", ls="-", alpha=0.2)

    # External comparison curves
    ext_curves = [
        ("MATHUSLA",    "external/MATHUSLA.csv",          "green",   "-"),
        ("CODEX-b",     "external/CODEX.csv",             "cyan",    "-"),
        ("ANUBIS",      "external/ANUBIS.csv",            "purple",  "-"),
        ("ANUBIS Opt",  "external/ANUBISOpt.csv",         "purple",  "--"),
        ("ANUBIS Cons", "external/ANUBISUpdateCons.csv",  "magenta", "--"),
    ]
    for label, path, color, ls in ext_curves:
        try:
            data = np.loadtxt(path, delimiter=",")
            ax2.loglog(data[:, 0], data[:, 1],
                       color=color, linewidth=2, linestyle=ls, label=label)
        except (FileNotFoundError, OSError):
            pass  # external curves not yet in dark_photon/external/

    ax2.legend(fontsize=8, loc='upper right')
    plt.tight_layout()

    plot_name = f"exclusion_Ntrack_{output_tag}.png"
    plt.savefig(os.path.join(args.outdir, plot_name), dpi=150)

    print(f"\nPlot saved: {os.path.join(args.outdir, plot_name)}")
    print(f"Normalization: n_generated={scan['n_generated']}, "
          f"n_llp_events={scan['n_llp_events']} "
          f"(0-LLP fraction: "
          f"{1 - scan['n_llp_events']/scan['n_generated']:.1%})")
