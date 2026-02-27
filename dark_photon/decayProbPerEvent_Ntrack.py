"""
Dark photon displaced-vertex analysis with realistic multi-body daughters.

Signal definition:
  LLP decays inside fiducial volume AND at least one daughter pair
  (from >= N_MIN_TRACKS charged daughters with p > P_CUT) has opening
  that yields separation inside [SEP_MIN, SEP_MAX] at the detector.

Inputs:
  - LLP CSV:      event,id,pt,eta,phi,momentum,mass
  - Daughter CSV: event,llp_idx,daughter_pid,pt,eta,phi,px,py,pz,energy,charge
  - Meta JSON:    n_generated, llp_pdg_id

Usage:
  python decayProbPerEvent_Ntrack.py output/data/dp_heavy_m15.csv \
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
    calculate_decay_length,
    cache_geometry,
    mesh_fiducial,
)

# Resolve paths relative to this script so external curves are found
# regardless of the working directory the user runs from.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Analysis cuts
P_CUT = 0.600          # GeV/c — minimum daughter momentum
N_MIN_TRACKS = 2       # minimum charged daughters to form candidate pairs
SEP_MIN = 0.001        # m — minimum separation at detector (1 mm)
SEP_MAX = 1.0          # m — maximum separation at detector (1 m)
HIST_LIFETIME_NS = 1000.0  # ns — default lifetime for separation histogram


def assign_llp_indices(llp_df):
    """
    Assign per-event LLP index (0, 1, ...) matching generator output ordering.

    Returns list of keys (event, llp_idx) aligned with llp_df rows.
    """
    keys = []
    event_counts = {}
    for ev in llp_df['event'].astype(int).values:
        idx = event_counts.get(ev, 0)
        event_counts[ev] = idx + 1
        keys.append((ev, idx))
    return keys


def build_event_row_indices(llp_df):
    """Build list of row-index arrays grouped by event."""
    event_rows = {}
    for i, ev in enumerate(llp_df['event'].astype(int).values):
        event_rows.setdefault(ev, []).append(i)
    return [np.array(rows, dtype=int) for rows in event_rows.values()]


def build_pair_thetas(dau_df, p_cut=P_CUT, n_min=N_MIN_TRACKS):
    """
    Build candidate daughter-pair opening angles per (event, llp_idx).

    Only daughters with p >= p_cut are kept.
    Pairs are built only if n_passing >= max(2, n_min).

    Returns:
      pair_thetas:   dict (event, llp_idx) -> np.array(theta_ij)
      passing_count: dict (event, llp_idx) -> n passing daughters
    """
    dau = dau_df.copy()
    dau['p'] = np.sqrt(dau['px']**2 + dau['py']**2 + dau['pz']**2)
    passing = dau[dau['p'] >= p_cut]

    pair_thetas = {}
    passing_count = {}
    min_needed = max(2, n_min)

    for (ev, llp_idx), group in passing.groupby(['event', 'llp_idx']):
        key = (int(ev), int(llp_idx))

        mom = group[['px', 'py', 'pz']].to_numpy(dtype=float)
        norms = np.linalg.norm(mom, axis=1)
        valid = np.isfinite(norms) & (norms > 0)
        mom = mom[valid]
        norms = norms[valid]

        n_pass = mom.shape[0]
        passing_count[key] = n_pass
        if n_pass < min_needed:
            continue

        thetas = []
        for i in range(n_pass - 1):
            for j in range(i + 1, n_pass):
                # arctan2(|cross|, dot) is numerically stable at all angles,
                # including the very small angles relevant here (~µrad–mrad).
                # arccos(dot) loses precision when cos ≈ 1.
                u = mom[i] / norms[i]
                v = mom[j] / norms[j]
                theta = float(np.arctan2(np.linalg.norm(np.cross(u, v)),
                                         np.dot(u, v)))
                if np.isfinite(theta):
                    thetas.append(theta)

        if thetas:
            pair_thetas[key] = np.array(thetas, dtype=float)

    return pair_thetas, passing_count


def integrate_decay_prob(d_start, d_end, decay_length):
    """
    Exponential decay probability integrated over [d_start, d_end].

    P = ∫ (1/lambda) exp(-d/lambda) dd = exp(-d_start/lambda)-exp(-d_end/lambda)
    """
    if decay_length <= 0 or d_end <= d_start:
        return 0.0
    return float(np.exp(-d_start / decay_length) - np.exp(-d_end / decay_length))


def merge_intervals(intervals):
    """
    Merge a list of (lo, hi) float intervals into disjoint sorted intervals.
    Intervals with lo >= hi are dropped.
    """
    valid = sorted((lo, hi) for lo, hi in intervals if lo < hi)
    if not valid:
        return []
    merged = [list(valid[0])]
    for lo, hi in valid[1:]:
        if lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return merged


def union_pair_window_probability(entry_d, exit_d, thetas, decay_length,
                                   sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    P(decay at a position where AT LEAST ONE daughter pair has its
    separation in [sep_min, sep_max]).

    For each pair with opening angle theta_ij:
      sep(d) = theta_ij * (exit_d - d)
      in-window: exit_d - sep_max/theta <= d <= exit_d - sep_min/theta

    Per-pair d-intervals are merged (union), then the exponential decay
    probability is integrated over the merged set — exact, not the
    pessimistic max-pair approximation.

    Also returns the theta of the single pair with highest individual
    probability (used only for visualisation in the histogram).
    """
    if not len(thetas) or sep_max <= sep_min or decay_length <= 0:
        return 0.0, np.nan

    intervals = []
    best_p = -1.0
    best_t = np.nan

    for theta in thetas:
        if theta <= 0:
            continue
        d_s = max(entry_d, exit_d - sep_max / theta)
        d_e = min(exit_d,  exit_d - sep_min / theta)
        if d_s >= d_e:
            continue
        intervals.append((d_s, d_e))
        p = integrate_decay_prob(d_s, d_e, decay_length)
        if p > best_p:
            best_p = p
            best_t = theta

    if not intervals:
        return 0.0, np.nan

    p_union = sum(integrate_decay_prob(lo, hi, decay_length)
                  for lo, hi in merge_intervals(intervals))
    return p_union, best_t


def evaluate_llps_for_lifetime(geo_cache, llp_keys, pair_thetas_by_llp,
                               lifetime_seconds,
                               sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    Evaluate per-LLP probabilities at one lifetime.

    Returns:
      decay_prob_raw:        P(decay in fiducial), no pair/separation cuts
      decay_prob_selected:   P(decay in fiducial AND at least one pair in sep window)
      best_theta:            highest-P single-pair theta (for histogram only)
      best_conditional:      P(any pair in window | decay in fiducial)
      has_pair_candidates:   LLP has >=2 passing daughters (and >=n_min total)
    """
    hits = geo_cache['hits'].astype(bool)
    entry_d = geo_cache['entry_d']
    exit_d = geo_cache['exit_d']
    momentum = geo_cache['momentum']
    mass = geo_cache['mass']

    n_llp = len(llp_keys)
    decay_prob_raw = np.zeros(n_llp)
    decay_prob_selected = np.zeros(n_llp)
    best_theta = np.full(n_llp, np.nan)
    best_conditional = np.zeros(n_llp)
    has_pair_candidates = np.zeros(n_llp, dtype=bool)

    for i, key in enumerate(llp_keys):
        if not hits[i]:
            continue

        dl = calculate_decay_length(momentum[i], mass[i], lifetime_seconds)
        p_raw = integrate_decay_prob(entry_d[i], exit_d[i], dl)
        decay_prob_raw[i] = p_raw

        thetas = pair_thetas_by_llp.get(key)
        if thetas is None or thetas.size == 0:
            continue
        has_pair_candidates[i] = True

        p_union, best_t = union_pair_window_probability(
            entry_d[i], exit_d[i], thetas, dl,
            sep_min=sep_min, sep_max=sep_max)

        decay_prob_selected[i] = p_union
        best_theta[i] = best_t
        if p_raw > 0:
            best_conditional[i] = p_union / p_raw

    return decay_prob_raw, decay_prob_selected, best_theta, best_conditional, has_pair_candidates


def compute_event_probabilities(event_row_indices, llp_prob):
    """Compute per-event P(>=1 LLP detected) from per-LLP probabilities."""
    p_event = np.zeros(len(event_row_indices))
    for j, rows in enumerate(event_row_indices):
        p_event[j] = 1.0 - np.prod(1.0 - llp_prob[rows])
    return p_event


def sample_selected_pair_separations(geo_cache, selected_theta, lifetime_seconds,
                                     n_samples_per_particle=200, rng_seed=42):
    """
    Sample separation distribution using the selected best pair per LLP.

    For each LLP with finite selected theta that hits the fiducial volume:
      - sample decay position inside [entry, exit] from truncated exponential
      - compute separation = theta_selected * d_remaining
      - weight each sample by P(decay in fiducial)/n_samples
    """
    rng = np.random.default_rng(rng_seed)

    hits = geo_cache['hits'].astype(bool)
    entry_d = geo_cache['entry_d']
    exit_d = geo_cache['exit_d']
    momentum = geo_cache['momentum']
    mass = geo_cache['mass']

    all_seps = []
    all_weights = []
    all_momenta = []

    for i in np.where(hits)[0]:
        theta = selected_theta[i]
        if not np.isfinite(theta) or theta <= 0:
            continue

        entry = entry_d[i]
        exit_ = exit_d[i]
        decay_length = calculate_decay_length(momentum[i], mass[i], lifetime_seconds)
        if decay_length <= 0:
            continue

        exp_entry = np.exp(-entry / decay_length)
        exp_exit = np.exp(-exit_ / decay_length)
        denom = exp_entry - exp_exit
        if denom < 1e-300:
            continue

        u = rng.uniform(0.0, 1.0, n_samples_per_particle)
        d_samples = -decay_length * np.log(exp_entry - u * denom)
        d_remaining = exit_ - d_samples

        sep = theta * d_remaining
        p_decay = float(exp_entry * (1.0 - np.exp(-(exit_ - entry) / decay_length)))
        w = p_decay / n_samples_per_particle

        all_seps.append(sep)
        all_weights.append(np.full(n_samples_per_particle, w))
        all_momenta.append(np.full(n_samples_per_particle, momentum[i]))

    if not all_seps:
        return np.array([]), np.array([]), np.array([])

    return (np.concatenate(all_seps),
            np.concatenate(all_weights),
            np.concatenate(all_momenta))


def analyze_Ntrack(llp_csv, dau_csv, geo_cache, lifetime_range,
                   xsec_fb=60e3, lumi_fb=3000,
                   n_generated=None, p_cut=P_CUT, n_min=N_MIN_TRACKS,
                   sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    Scan over lifetimes, computing expected signal with pair+separation cuts.

    For each LLP that hits fiducial volume:
      P_detect = P(decay in fiducial AND selected-best pair in sep window)

    The selected pair is the pair with maximal in-window probability at each
    lifetime among daughters passing p > p_cut.
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

    llp_keys = assign_llp_indices(llp_df)
    event_row_indices = build_event_row_indices(llp_df)
    pair_thetas_by_llp, passing_count = build_pair_thetas(
        dau_df, p_cut=p_cut, n_min=n_min)

    hits = geo_cache['hits'].astype(bool)
    n_llp = len(llp_df)

    n_track_pass = int(sum(passing_count.get(key, 0) >= n_min for key in llp_keys))
    n_pair_candidates = int(sum(key in pair_thetas_by_llp for key in llp_keys))
    n_hits = int(hits.sum())
    n_pair_hits = int(sum(hits[i] and llp_keys[i] in pair_thetas_by_llp
                          for i in range(n_llp)))

    print(f"\nTrack preselection: {n_track_pass}/{n_llp} LLPs have "
          f"≥{n_min} charged daughters with p>{p_cut*1000:.0f} MeV")
    print(f"Pair candidates: {n_pair_candidates}/{n_llp} LLPs provide at least one pair")
    print(f"Geometric acceptance: {n_hits}/{n_llp} LLPs hit fiducial volume")
    print(f"Pair+geometry overlap: {n_pair_hits}/{n_llp}")
    print(f"Separation window: {sep_min*1000:.1f} mm <= sep <= {sep_max:.2f} m")

    results = {
        'lifetimes': lifetime_range,
        'mean_decay_prob_Ntrack': [],
        'mean_decay_prob_no_cuts': [],
        'mean_at_least_one': [],
        'mean_at_least_one_no_cuts': [],
        'exclusion': [],
        'exclusion_no_cuts': [],
        'mean_pair_acceptance': [],
        'n_generated': n_generated,
        'n_llp_events': n_llp_events,
    }

    for lifetime in tqdm(lifetime_range, desc="Scanning lifetimes"):
        decay_prob_raw, decay_prob_selected, _, best_conditional, has_pairs = (
            evaluate_llps_for_lifetime(
                geo_cache,
                llp_keys,
                pair_thetas_by_llp,
                lifetime_seconds=lifetime,
                sep_min=sep_min,
                sep_max=sep_max,
            )
        )

        p_event = compute_event_probabilities(event_row_indices, decay_prob_selected)
        p_event_nc = compute_event_probabilities(event_row_indices, decay_prob_raw)

        mean_p1 = p_event.sum() / n_generated
        mean_p1_nc = p_event_nc.sum() / n_generated

        hit_mask = hits
        mean_single = decay_prob_selected[hit_mask].mean() if hit_mask.any() else 0.0
        mean_single_nc = decay_prob_raw[hit_mask].mean() if hit_mask.any() else 0.0

        pair_hit_mask = hit_mask & has_pairs
        mean_pair_acc = (
            best_conditional[pair_hit_mask].mean() if pair_hit_mask.any() else 0.0
        )

        denom = mean_p1 * lumi_fb * xsec_fb
        denom_nc = mean_p1_nc * lumi_fb * xsec_fb

        results['mean_decay_prob_Ntrack'].append(mean_single)
        results['mean_decay_prob_no_cuts'].append(mean_single_nc)
        results['mean_at_least_one'].append(mean_p1)
        results['mean_at_least_one_no_cuts'].append(mean_p1_nc)
        results['exclusion'].append(3 / denom if denom > 0 else np.inf)
        results['exclusion_no_cuts'].append(3 / denom_nc if denom_nc > 0 else np.inf)
        results['mean_pair_acceptance'].append(mean_pair_acc)

    state = {
        'llp_df': llp_df,
        'llp_keys': llp_keys,
        'event_row_indices': event_row_indices,
        'pair_thetas_by_llp': pair_thetas_by_llp,
    }
    return results, state


# ==================================================================
# Main
# ==================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dark photon displaced-vertex analysis with pair separation window")
    parser.add_argument("csv_file", help="LLP CSV (event,id,pt,eta,phi,momentum,mass)")
    parser.add_argument("--xsec", type=float, default=60e3,
                        help="production cross-section in fb (default: 60e3 for pp→h)")
    parser.add_argument("--lumi", type=float, default=3000,
                        help="integrated luminosity in fb⁻¹ (default: 3000)")
    parser.add_argument("--outdir", default="output",
                        help="output root directory (writes images to <outdir>/images)")
    parser.add_argument("--n-events", type=int, default=None,
                        help="total generated events (auto-read from _meta.json)")
    parser.add_argument("--n-min-tracks", type=int, default=N_MIN_TRACKS,
                        help=f"minimum charged daughters per LLP (default: {N_MIN_TRACKS})")
    parser.add_argument("--p-cut", type=float, default=P_CUT,
                        help=f"minimum daughter momentum in GeV (default: {P_CUT})")
    parser.add_argument("--sep-min", type=float, default=SEP_MIN,
                        help=f"minimum separation in m (default: {SEP_MIN})")
    parser.add_argument("--sep-max", type=float, default=SEP_MAX,
                        help=f"maximum separation in m (default: {SEP_MAX})")
    parser.add_argument("--hist-lifetime-ns", type=float, default=HIST_LIFETIME_NS,
                        help=f"lifetime in ns for separation histogram (default: {HIST_LIFETIME_NS})")
    args = parser.parse_args()

    xsec_arg_given = any(
        tok == "--xsec" or tok.startswith("--xsec=")
        for tok in sys.argv[1:]
    )
    if not xsec_arg_given:
        print("WARNING: --xsec not provided; using default 60e3 fb (pp→h).")

    if args.sep_min <= 0 or args.sep_max <= args.sep_min:
        print("ERROR: require 0 < --sep-min < --sep-max")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    data_dir = os.path.join(args.outdir, "data")
    image_dir = os.path.join(args.outdir, "images")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(image_dir, exist_ok=True)

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
    print(f"LIFETIME SCAN (>= {args.n_min_tracks} daughters, p > {args.p_cut*1000:.0f} MeV, "
          f"{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)")
    print("=" * 50)
    print(f"  sigma = {args.xsec:.3g} fb,  L = {args.lumi:.0f} fb^-1")
    print(f"  n_generated = {n_generated}")

    lifetimes = np.logspace(-9.5, -4.5, 20)
    scan, state = analyze_Ntrack(
        sample_csv, dau_csv, geo_cache, lifetimes,
        xsec_fb=args.xsec, lumi_fb=args.lumi,
        n_generated=n_generated,
        p_cut=args.p_cut,
        n_min=args.n_min_tracks,
        sep_min=args.sep_min,
        sep_max=args.sep_max,
    )

    # === Exclusion plotting ===
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax1 = axes[0]
    ax1.loglog(lifetimes * 1e9, scan['mean_at_least_one'],
               'r-', linewidth=2, label='Pair+window selection')
    ax1.loglog(lifetimes * 1e9, scan['mean_at_least_one_no_cuts'],
               'r--', linewidth=2, alpha=0.5, label='Geometric only')
    ax1.set_xlabel('Lifetime (ns)')
    ax1.set_ylabel('Mean P(>=1 LLP detected)')
    ax1.set_title(
        'Detection probability\n'
        f'(>= {args.n_min_tracks} daughters, p > {args.p_cut*1000:.0f} MeV, '\
        f'{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)'
    )
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend()

    ax2 = axes[1]
    ax2.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion'],
               color='blue', linewidth=2,
               label='PX56 (pair+window)')
    ax2.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion_no_cuts'],
               color='blue', linewidth=2, linestyle='--', alpha=0.5,
               label='PX56 (geometric only)')
    ax2.set_xlabel(r'$c\tau$ (m)')
    ax2.set_ylabel(r"BR$(h \to A'A')_{\min}$")
    ax2.set_title('Exclusion sensitivity (3 signal events)')
    ax2.grid(True, which="both", ls="-", alpha=0.2)

    ext_curves = [
        ("MATHUSLA",    "external/MATHUSLA.csv",          "green",   "-"),
        ("CODEX-b",     "external/CODEX.csv",             "cyan",    "-"),
        ("ANUBIS",      "external/ANUBIS.csv",            "purple",  "-"),
        ("ANUBIS Opt",  "external/ANUBISOpt.csv",         "purple",  "--"),
        ("ANUBIS Cons", "external/ANUBISUpdateCons.csv",  "magenta", "--"),
    ]
    for label, path, color, ls in ext_curves:
        try:
            data = np.loadtxt(os.path.join(_SCRIPT_DIR, path), delimiter=",")
            ax2.loglog(data[:, 0], data[:, 1],
                       color=color, linewidth=2, linestyle=ls, label=label)
        except (FileNotFoundError, OSError):
            pass

    ax2.legend(fontsize=8, loc='upper right')
    plt.tight_layout()

    exclusion_plot_name = f"exclusion_Ntrack_{output_tag}.png"
    exclusion_plot_path = os.path.join(image_dir, exclusion_plot_name)
    plt.savefig(exclusion_plot_path, dpi=150)
    plt.close(fig)

    print(f"\nPlot saved: {exclusion_plot_path}")
    print(f"Normalization: n_generated={scan['n_generated']}, "
          f"n_llp_events={scan['n_llp_events']} "
          f"(0-LLP fraction: "
          f"{1 - scan['n_llp_events']/scan['n_generated']:.1%})")

    # === Separation histogram using best pair at chosen lifetime ===
    hist_lifetime_s = args.hist_lifetime_ns * 1e-9
    _, _, best_theta_hist, _, _ = evaluate_llps_for_lifetime(
        geo_cache,
        state['llp_keys'],
        state['pair_thetas_by_llp'],
        lifetime_seconds=hist_lifetime_s,
        sep_min=args.sep_min,
        sep_max=args.sep_max,
    )

    seps, weights, momenta = sample_selected_pair_separations(
        geo_cache,
        best_theta_hist,
        lifetime_seconds=hist_lifetime_s,
        n_samples_per_particle=200,
    )

    separation_plot_name = f"separation_histogram_{output_tag}.png"
    separation_plot_path = os.path.join(image_dir, separation_plot_name)

    if len(seps) > 0:
        fig_sep, axes_sep = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes_sep[0]
        lin_xmax = 1.05 * args.sep_max
        bins = np.linspace(0, lin_xmax, 100)
        ax.hist(seps, bins=bins, weights=weights, color='steelblue',
                edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.axvline(args.sep_min, color='red', linestyle='--', linewidth=2,
                   label=f'min sep = {args.sep_min*1000:.0f} mm')
        ax.axvline(args.sep_max, color='red', linestyle='-', linewidth=2,
                   label=f'max sep = {args.sep_max*100:.0f} cm')
        ax.axvspan(args.sep_min, args.sep_max, color='green', alpha=0.1,
                   label='Accepted window')
        ax.set_xlabel('Separation at detector (m)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Selected best-pair separation (tau = {args.hist_lifetime_ns:.0f} ns)')
        ax.legend(fontsize=9)
        ax.set_xlim(0, lin_xmax)

        ax2 = axes_sep[1]
        bins_log = np.logspace(-4, 1, 100)
        ax2.hist(seps, bins=bins_log, weights=weights, color='steelblue',
                 edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.axvline(args.sep_min, color='red', linestyle='--', linewidth=2,
                    label=f'min sep = {args.sep_min*1000:.0f} mm')
        ax2.axvline(args.sep_max, color='red', linestyle='-', linewidth=2,
                    label=f'max sep = {args.sep_max*100:.0f} cm')
        ax2.axvspan(args.sep_min, args.sep_max, color='green', alpha=0.1,
                    label='Accepted window')
        ax2.set_xscale('log')
        ax2.set_xlabel('Separation at detector (m)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('Log-scale separation')
        ax2.legend(fontsize=9)

        ax3 = axes_sep[2]
        mask_finite = np.isfinite(seps) & (seps > 0)
        h = ax3.hist2d(momenta[mask_finite], seps[mask_finite] * 100,
                       bins=[np.linspace(0, 500, 50), np.linspace(0, 500, 50)],
                       weights=weights[mask_finite],
                       cmap='viridis', cmin=1e-20)
        ax3.axhline(args.sep_min * 100, color='red', linestyle='--', linewidth=2,
                    label=f'min = {args.sep_min*1000:.0f} mm')
        ax3.axhline(args.sep_max * 100, color='red', linestyle='-', linewidth=2,
                    label=f'max = {args.sep_max*100:.0f} cm')
        ax3.set_xlabel('LLP momentum (GeV/c)')
        ax3.set_ylabel('Separation at detector (cm)')
        ax3.set_title('Selected best-pair separation vs LLP momentum')
        ax3.legend(fontsize=9, loc='upper right')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig(separation_plot_path, dpi=150)
        plt.close(fig_sep)

        if weights.sum() > 0:
            in_window = (seps >= args.sep_min) & (seps <= args.sep_max)
            frac_accepted = weights[in_window].sum() / weights.sum()
            print(f"Separation summary (selected best-pair @ {args.hist_lifetime_ns:.0f} ns):")
            print(f"  Fraction in window: {frac_accepted:.3f}")
            print(f"  Median separation (all): {np.median(seps)*100:.1f} cm")
            accepted_seps = seps[in_window]
            if len(accepted_seps) > 0:
                print(f"  Median separation (in window): {np.median(accepted_seps)*100:.1f} cm")
    else:
        fig_sep, ax = plt.subplots(1, 1, figsize=(8, 5))
        ax.text(0.5, 0.5,
                'No selected pair separations available\nfor the requested sample/lifetime.',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(separation_plot_path, dpi=150)
        plt.close(fig_sep)

    print(f"Separation plot saved: {separation_plot_path}")
