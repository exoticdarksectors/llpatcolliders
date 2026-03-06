"""
N-track displaced-vertex analysis for multi-body LLP decays.

Signal definition:
  LLP decays inside fiducial volume AND at least one daughter pair
  (from >= N_MIN_TRACKS charged daughters with p > P_CUT) has opening
  that yields separation inside [SEP_MIN, SEP_MAX] at the detector.

Inputs:
  - LLP CSV:      event,id,pt,eta,phi,momentum,mass
  - Daughter CSV:  event,llp_idx,daughter_pid,pt,eta,phi,px,py,pz,energy,charge
  - Meta JSON:     n_generated, llp_pdg_id

Usage:
  python decayProbPerEvent_Ntrack.py <csv_file> --xsec <fb> [options]
"""

import os
import sys
import json
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'geometry'))
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

from utils import infer_sample_mass, overlay_mass_matched_external_curves, exclusion_ylabel

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

# Analysis cuts
P_CUT = 0.600          # GeV/c — minimum daughter momentum
N_MIN_TRACKS = 2       # minimum charged daughters to form candidate pairs
SEP_MIN = 0.001        # m — minimum separation at detector (1 mm)
SEP_MAX = 1.0          # m — maximum separation at detector (1 m)
HIST_LIFETIME_NS = 1000.0  # ns — default lifetime for separation histogram
DEFAULT_LIFETIME_MIN_NS = 10**-2
DEFAULT_LIFETIME_MAX_NS = 10**5.5
DEFAULT_LIFETIME_POINTS = 80
DEFAULT_EPS2_MIN = 1e-15
DEFAULT_EPS2_MAX = 1e-12
DEFAULT_EPS2_POINTS = 60


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

    P = exp(-d_start/lambda) - exp(-d_end/lambda)

    Uses expm1 reparametrisation to avoid catastrophic cancellation
    when (d_end - d_start) << decay_length.
    """
    if decay_length <= 0 or d_end <= d_start:
        return 0.0
    lead = -d_start / decay_length
    delta = -(d_end - d_start) / decay_length
    return float(np.exp(lead) * (-np.expm1(delta)))


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


def analyze_meson_Ntrack(llp_csv, dau_csv, geo_cache, eps2_range,
                         dp_mass, parent_meson=None, br_func=None,
                         xsec_fb=80e9, lumi_fb=3000,
                         n_generated=None, p_cut=P_CUT,
                         n_min=N_MIN_TRACKS,
                         sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    ε² scan mode: scan over ε² (kinetic mixing squared).

    Used for meson production (eta/omega → A') and Drell-Yan (qq-bar → A').

    For each ε²:
      1. Compute cτ = ℏc / (ε² × Γ₀(m_A'))
      2. Compute BR_prod(ε²) via br_func
      3. Evaluate geometric acceptance at cτ
      4. N_signal = BR_prod × σ × L × <P(detect)>
      5. Store ε²_min for N_signal >= 3
    """
    sys.path.insert(0, os.path.join(_SCRIPT_DIR, '..', 'dark_photon'))
    from dp_meson_brs import (
        dp_ctau_m, dp_width_eps1,
        br_eta_to_dp_gamma, br_omega_to_dp_pi0,
    )

    if br_func is None:
        br_func = {
            "eta": br_eta_to_dp_gamma,
            "omega": br_omega_to_dp_pi0,
        }[parent_meson]

    llp_df = pd.read_csv(llp_csv)
    llp_df.columns = llp_df.columns.str.strip()

    dau_df = pd.read_csv(dau_csv)
    dau_df.columns = dau_df.columns.str.strip()

    n_llp_events = llp_df['event'].nunique()
    if n_generated is None:
        n_generated = n_llp_events
    if n_generated < n_llp_events:
        n_generated = n_llp_events

    llp_keys = assign_llp_indices(llp_df)
    event_row_indices = build_event_row_indices(llp_df)
    pair_thetas_by_llp, passing_count = build_pair_thetas(
        dau_df, p_cut=p_cut, n_min=n_min)

    hits = geo_cache['hits'].astype(bool)
    n_llp = len(llp_df)

    gamma0 = dp_width_eps1(dp_mass)
    mode_label = f"parent={parent_meson}" if parent_meson else "Drell-Yan"
    print(f"\nε² scan mode: {mode_label}, m_A'={dp_mass} GeV")
    print(f"  Γ₀(A') [ε=1] = {gamma0:.4e} GeV")
    print(f"  {n_llp} LLPs, {n_llp_events} events with LLP, "
          f"{n_generated} total events")

    results = {
        'eps2_range': eps2_range,
        'ctau_m': [],
        'br_production': [],
        'mean_at_least_one': [],
        'mean_at_least_one_no_cuts': [],
        'n_signal': [],
        'n_signal_no_cuts': [],
        'n_generated': n_generated,
        'n_llp_events': n_llp_events,
    }

    for eps2 in tqdm(eps2_range, desc="Scanning ε²"):
        ctau = dp_ctau_m(eps2, dp_mass)
        lifetime_s = ctau / SPEED_OF_LIGHT
        br_prod = br_func(eps2, dp_mass)

        decay_prob_raw, decay_prob_selected, _, _, _ = (
            evaluate_llps_for_lifetime(
                geo_cache,
                llp_keys,
                pair_thetas_by_llp,
                lifetime_seconds=lifetime_s,
                sep_min=sep_min,
                sep_max=sep_max,
            )
        )

        p_event = compute_event_probabilities(event_row_indices, decay_prob_selected)
        p_event_nc = compute_event_probabilities(event_row_indices, decay_prob_raw)

        mean_p1 = p_event.sum() / n_generated
        mean_p1_nc = p_event_nc.sum() / n_generated

        n_sig = br_prod * xsec_fb * lumi_fb * mean_p1
        n_sig_nc = br_prod * xsec_fb * lumi_fb * mean_p1_nc

        results['ctau_m'].append(ctau)
        results['br_production'].append(br_prod)
        results['mean_at_least_one'].append(mean_p1)
        results['mean_at_least_one_no_cuts'].append(mean_p1_nc)
        results['n_signal'].append(n_sig)
        results['n_signal_no_cuts'].append(n_sig_nc)

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
    parser.add_argument("--xsec", type=float, default=54700,
                        help="production cross-section in fb (default: 54700 for σ(gg→h) at 14 TeV)")
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
    parser.add_argument("--lifetime-min-ns", type=float, default=DEFAULT_LIFETIME_MIN_NS,
                        help="minimum lifetime in ns for Higgs-mode scan")
    parser.add_argument("--lifetime-max-ns", type=float, default=DEFAULT_LIFETIME_MAX_NS,
                        help="maximum lifetime in ns for Higgs-mode scan")
    parser.add_argument("--lifetime-points", type=int, default=DEFAULT_LIFETIME_POINTS,
                        help="number of lifetime points in Higgs-mode scan")
    parser.add_argument("--eps2-min", type=float, default=DEFAULT_EPS2_MIN,
                        help="minimum epsilon^2 for meson-mode scan")
    parser.add_argument("--eps2-max", type=float, default=DEFAULT_EPS2_MAX,
                        help="maximum epsilon^2 for meson-mode scan")
    parser.add_argument("--eps2-points", type=int, default=DEFAULT_EPS2_POINTS,
                        help="number of epsilon^2 points in meson-mode scan")
    parser.add_argument("--production", choices=["higgs", "meson", "drell_yan"],
                        default="higgs",
                        help="production mode: higgs (default, BR_min vs cτ), "
                             "meson (ε² scan for meson→A'), or "
                             "drell_yan (ε² scan for qq̄→A')")
    parser.add_argument("--parent-meson", choices=["eta", "omega"],
                        help="parent meson (required for --production meson)")
    parser.add_argument("--dp-mass", type=float,
                        help="A' mass in GeV (required for --production meson/drell_yan, "
                             "for ε²→cτ conversion)")
    parser.add_argument("--external-dir", default=None,
                        help="path to external/ comparison curves (default: none)")
    args = parser.parse_args()

    if args.production in ("meson", "drell_yan"):
        if args.dp_mass is None:
            print(f"ERROR: --dp-mass required for --production {args.production}")
            sys.exit(1)
    if args.production == "meson":
        if args.parent_meson is None:
            print("ERROR: --parent-meson required for --production meson")
            sys.exit(1)

    xsec_arg_given = any(
        tok == "--xsec" or tok.startswith("--xsec=")
        for tok in sys.argv[1:]
    )
    if not xsec_arg_given:
        if args.production == "meson":
            print("WARNING: --xsec not provided; using default 54700 fb (σ(gg→h) at 14 TeV).\n"
                  "  For meson production, use --xsec 80000000000 (σ_inel ≈ 80 mb).")
        elif args.production == "drell_yan":
            print("WARNING: --xsec not provided; using default 54700 fb (σ(gg→h) at 14 TeV).\n"
                  "  For Drell-Yan, use σ(pp→A') at ε²=1 from Pythia8 log.")
        else:
            print("WARNING: --xsec not provided; using default 54700 fb (σ(gg→h) at 14 TeV).")

    if args.sep_min <= 0 or args.sep_max <= args.sep_min:
        print("ERROR: require 0 < --sep-min < --sep-max")
        sys.exit(1)
    if args.lifetime_min_ns <= 0 or args.lifetime_max_ns <= args.lifetime_min_ns:
        print("ERROR: require 0 < --lifetime-min-ns < --lifetime-max-ns")
        sys.exit(1)
    if args.eps2_min <= 0 or args.eps2_max <= args.eps2_min:
        print("ERROR: require 0 < --eps2-min < --eps2-max")
        sys.exit(1)
    if args.lifetime_points < 2 or args.eps2_points < 2:
        print("ERROR: --lifetime-points and --eps2-points must both be >= 2")
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
    sample_mass = infer_sample_mass(geo_cache['mass'])
    if sample_mass is not None:
        print(f"Inferred LLP mass from sample: {sample_mass:.3g} GeV")

    if args.production == "meson":
        # ── Meson production mode: scan ε² ──
        print("\n" + "=" * 50)
        print(f"MESON ε² SCAN (>= {args.n_min_tracks} daughters, "
              f"p > {args.p_cut*1000:.0f} MeV, "
              f"{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)")
        print("=" * 50)
        print(f"  σ_inel = {args.xsec:.3g} fb,  L = {args.lumi:.0f} fb^-1")
        print(f"  parent meson = {args.parent_meson},  m_A' = {args.dp_mass} GeV")
        print(f"  n_generated = {n_generated}")

        eps2_range = np.logspace(np.log10(args.eps2_min),
                                 np.log10(args.eps2_max),
                                 args.eps2_points)
        scan, state = analyze_meson_Ntrack(
            sample_csv, dau_csv, geo_cache, eps2_range,
            dp_mass=args.dp_mass, parent_meson=args.parent_meson,
            xsec_fb=args.xsec, lumi_fb=args.lumi,
            n_generated=n_generated,
            p_cut=args.p_cut,
            n_min=args.n_min_tracks,
            sep_min=args.sep_min,
            sep_max=args.sep_max,
        )

        ctau_arr = np.array(scan['ctau_m'])
        n_sig_arr = np.array(scan['n_signal'])
        n_sig_nc_arr = np.array(scan['n_signal_no_cuts'])

        # === Meson exclusion plotting ===
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        ax1 = axes[0]
        ax1.loglog(eps2_range, scan['mean_at_least_one'],
                   'r-', linewidth=2, label='Pair+window selection')
        ax1.loglog(eps2_range, scan['mean_at_least_one_no_cuts'],
                   'r--', linewidth=2, alpha=0.5, label='Geometric only')
        ax1.set_xlabel(r'$\varepsilon^2$')
        ax1.set_ylabel('Mean P(>=1 LLP detected)')
        ax1.set_title(
            f'Detection probability (m_A\' = {args.dp_mass} GeV)\n'
            f'(>= {args.n_min_tracks} daughters, p > {args.p_cut*1000:.0f} MeV, '
            f'{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)'
        )
        ax1.grid(True, which="both", ls="-", alpha=0.2)
        ax1.legend()

        # Secondary x-axis: cτ
        ax1b = ax1.twiny()
        ctau_ticks = [1e-6, 1e-4, 1e-2, 1, 100]
        ax1b.set_xscale('log')
        ax1b.set_xlim(ctau_arr[0], ctau_arr[-1])
        ax1b.set_xlabel(r'$c\tau$ (m)')

        ax2 = axes[1]
        ax2.loglog(ctau_arr, n_sig_arr,
                   color='blue', linewidth=2,
                   label='PX56 (pair+window)')
        ax2.loglog(ctau_arr, n_sig_nc_arr,
                   color='blue', linewidth=2, linestyle='--', alpha=0.5,
                   label='PX56 (geometric only)')
        ax2.axhline(3, color='red', linestyle=':', linewidth=1.5,
                    label='3 signal events')
        ax2.set_xlabel(r'$c\tau$ (m)')
        meson_label = {"eta": r"$\eta$", "omega": r"$\omega$"}[args.parent_meson]
        ax2.set_ylabel(f'Expected signal events ({meson_label} channel)')
        ax2.set_title(
            f'Signal yield vs $c\\tau$ '
            f'($m_{{A\'}}$ = {args.dp_mass} GeV, {meson_label} production)')
        ax2.grid(True, which="both", ls="-", alpha=0.2)
        ax2.legend(fontsize=8, loc='upper right')

        # Data-driven y-limits: ~4 decades padding around signal peak
        _yvals = np.concatenate([n_sig_arr[n_sig_arr > 0],
                                 n_sig_nc_arr[n_sig_nc_arr > 0]])
        if len(_yvals) > 0:
            _ymax = _yvals.max()
            ax2.set_ylim(_ymax * 1e-4, max(_ymax * 1e4, 10))

        plt.tight_layout()

        exclusion_plot_name = f"exclusion_meson_{output_tag}.png"
        exclusion_plot_path = os.path.join(image_dir, exclusion_plot_name)
        plt.savefig(exclusion_plot_path, dpi=150)
        plt.close(fig)

        print(f"\nPlot saved: {exclusion_plot_path}")
        print(f"Normalization: n_generated={scan['n_generated']}, "
              f"n_llp_events={scan['n_llp_events']}")

        # Print peak sensitivity
        if np.any(n_sig_arr > 0):
            i_peak = np.argmax(n_sig_arr)
            print(f"Peak signal: N={n_sig_arr[i_peak]:.2e} at "
                  f"ε²={eps2_range[i_peak]:.2e}, "
                  f"cτ={ctau_arr[i_peak]:.2e} m")
        else:
            print("No signal detected at any ε² value.")

        # Save scan data
        scan_csv_name = f"meson_scan_{output_tag}.csv"
        scan_csv_path = os.path.join(data_dir, scan_csv_name)
        np.savetxt(scan_csv_path,
                   np.column_stack([eps2_range, ctau_arr,
                                    scan['br_production'],
                                    scan['mean_at_least_one'],
                                    n_sig_arr]),
                   delimiter=",",
                   header="eps2,ctau_m,br_production,mean_p1,n_signal",
                   comments="")
        print(f"Scan data saved: {scan_csv_path}")

        # Use peak-signal cτ for separation histogram
        if np.any(n_sig_arr > 0):
            hist_lifetime_s = ctau_arr[np.argmax(n_sig_arr)] / SPEED_OF_LIGHT
            print(f"Using peak-sensitivity lifetime for separation histogram: "
                  f"cτ={ctau_arr[np.argmax(n_sig_arr)]:.2e} m "
                  f"({hist_lifetime_s * 1e9:.2e} ns)")
        else:
            hist_lifetime_s = args.hist_lifetime_ns * 1e-9

    elif args.production == "drell_yan":
        # ── Drell-Yan production mode: scan ε² ──
        print("\n" + "=" * 50)
        print(f"DRELL-YAN ε² SCAN (>= {args.n_min_tracks} daughters, "
              f"p > {args.p_cut*1000:.0f} MeV, "
              f"{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)")
        print("=" * 50)
        print(f"  σ(pp→A')|_{{ε²=1}} = {args.xsec:.3g} fb,  "
              f"L = {args.lumi:.0f} fb^-1")
        print(f"  m_A' = {args.dp_mass} GeV")
        print(f"  n_generated = {n_generated}")

        eps2_range = np.logspace(np.log10(args.eps2_min),
                                 np.log10(args.eps2_max),
                                 args.eps2_points)

        # Drell-Yan: br_prod = ε² (σ already at ε²=1, so scaling is just ε²)
        dy_br_func = lambda eps2, m: eps2

        scan, state = analyze_meson_Ntrack(
            sample_csv, dau_csv, geo_cache, eps2_range,
            dp_mass=args.dp_mass, br_func=dy_br_func,
            xsec_fb=args.xsec, lumi_fb=args.lumi,
            n_generated=n_generated,
            p_cut=args.p_cut,
            n_min=args.n_min_tracks,
            sep_min=args.sep_min,
            sep_max=args.sep_max,
        )

        ctau_arr = np.array(scan['ctau_m'])
        n_sig_arr = np.array(scan['n_signal'])
        n_sig_nc_arr = np.array(scan['n_signal_no_cuts'])

        # === Drell-Yan exclusion plotting ===
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        ax1 = axes[0]
        ax1.loglog(eps2_range, scan['mean_at_least_one'],
                   'r-', linewidth=2, label='Pair+window selection')
        ax1.loglog(eps2_range, scan['mean_at_least_one_no_cuts'],
                   'r--', linewidth=2, alpha=0.5, label='Geometric only')
        ax1.set_xlabel(r'$\varepsilon^2$')
        ax1.set_ylabel('Mean P(>=1 LLP detected)')
        ax1.set_title(
            f'Detection probability (m_A\' = {args.dp_mass} GeV, Drell-Yan)\n'
            f'(>= {args.n_min_tracks} daughters, p > {args.p_cut*1000:.0f} MeV, '
            f'{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)'
        )
        ax1.grid(True, which="both", ls="-", alpha=0.2)
        ax1.legend()

        # Secondary x-axis: cτ
        ax1b = ax1.twiny()
        ax1b.set_xscale('log')
        ax1b.set_xlim(ctau_arr[0], ctau_arr[-1])
        ax1b.set_xlabel(r'$c\tau$ (m)')

        ax2 = axes[1]
        ax2.loglog(ctau_arr, n_sig_arr,
                   color='blue', linewidth=2,
                   label='PX56 (pair+window)')
        ax2.loglog(ctau_arr, n_sig_nc_arr,
                   color='blue', linewidth=2, linestyle='--', alpha=0.5,
                   label='PX56 (geometric only)')
        ax2.axhline(3, color='red', linestyle=':', linewidth=1.5,
                    label='3 signal events')
        ax2.set_xlabel(r'$c\tau$ (m)')
        ax2.set_ylabel(r'Expected signal events (Drell-Yan)')
        ax2.set_title(
            f'Signal yield vs $c\\tau$ '
            r'($m_{A\prime}$' + f' = {args.dp_mass} GeV, Drell-Yan)')
        ax2.grid(True, which="both", ls="-", alpha=0.2)
        ax2.legend(fontsize=8, loc='upper right')

        # Data-driven y-limits: ~4 decades padding around signal peak
        _yvals = np.concatenate([n_sig_arr[n_sig_arr > 0],
                                 n_sig_nc_arr[n_sig_nc_arr > 0]])
        if len(_yvals) > 0:
            _ymax = _yvals.max()
            ax2.set_ylim(_ymax * 1e-4, max(_ymax * 1e4, 10))

        plt.tight_layout()

        exclusion_plot_name = f"exclusion_dy_{output_tag}.png"
        exclusion_plot_path = os.path.join(image_dir, exclusion_plot_name)
        plt.savefig(exclusion_plot_path, dpi=150)
        plt.close(fig)

        print(f"\nPlot saved: {exclusion_plot_path}")
        print(f"Normalization: n_generated={scan['n_generated']}, "
              f"n_llp_events={scan['n_llp_events']}")

        # Print peak sensitivity
        if np.any(n_sig_arr > 0):
            i_peak = np.argmax(n_sig_arr)
            print(f"Peak signal: N={n_sig_arr[i_peak]:.2e} at "
                  f"ε²={eps2_range[i_peak]:.2e}, "
                  f"cτ={ctau_arr[i_peak]:.2e} m")
        else:
            print("No signal detected at any ε² value.")

        # Save scan data
        scan_csv_name = f"dy_scan_{output_tag}.csv"
        scan_csv_path = os.path.join(data_dir, scan_csv_name)
        np.savetxt(scan_csv_path,
                   np.column_stack([eps2_range, ctau_arr,
                                    scan['br_production'],
                                    scan['mean_at_least_one'],
                                    n_sig_arr]),
                   delimiter=",",
                   header="eps2,ctau_m,br_production,mean_p1,n_signal",
                   comments="")
        print(f"Scan data saved: {scan_csv_path}")

        # Use peak-signal cτ for separation histogram
        if np.any(n_sig_arr > 0):
            hist_lifetime_s = ctau_arr[np.argmax(n_sig_arr)] / SPEED_OF_LIGHT
            print(f"Using peak-sensitivity lifetime for separation histogram: "
                  f"cτ={ctau_arr[np.argmax(n_sig_arr)]:.2e} m "
                  f"({hist_lifetime_s * 1e9:.2e} ns)")
        else:
            hist_lifetime_s = args.hist_lifetime_ns * 1e-9

    else:
        # ── Higgs production mode: scan cτ (original behavior) ──
        print("\n" + "=" * 50)
        print(f"LIFETIME SCAN (>= {args.n_min_tracks} daughters, "
              f"p > {args.p_cut*1000:.0f} MeV, "
              f"{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)")
        print("=" * 50)
        print(f"  sigma = {args.xsec:.3g} fb,  L = {args.lumi:.0f} fb^-1")
        print(f"  n_generated = {n_generated}")

        lifetimes = np.logspace(np.log10(args.lifetime_min_ns),
                                np.log10(args.lifetime_max_ns),
                                args.lifetime_points) * 1e-9
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
        mean_p1 = np.array(scan['mean_at_least_one'])
        mean_p1_nc = np.array(scan['mean_at_least_one_no_cuts'])
        # Plot geometric-only first (behind), then pair+window on top
        ax1.loglog(lifetimes * 1e9, mean_p1_nc,
                   'b--', linewidth=2, alpha=0.5, label='Geometric only')
        ax1.loglog(lifetimes * 1e9, mean_p1,
                   'r-', linewidth=2, label='Pair+window selection')
        ax1.set_xlabel('Lifetime (ns)')
        ax1.set_ylabel('Mean P(>=1 LLP detected)')
        ax1.set_title(
            'Detection probability\n'
            f'(>= {args.n_min_tracks} daughters, p > {args.p_cut*1000:.0f} MeV, '
            f'{args.sep_min*1000:.1f} mm <= sep <= {args.sep_max:.2f} m)'
        )
        ax1.grid(True, which="both", ls="-", alpha=0.2)
        ax1.legend()
        # Annotate if curves overlap (pair+window ≈ geometric)
        _ratio1 = mean_p1 / np.clip(mean_p1_nc, 1e-30, None)
        _valid1 = _ratio1[mean_p1 > 0]
        if len(_valid1) > 0 and np.median(_valid1) > 0.95:
            ax1.text(0.05, 0.05,
                     'Curves overlap:\npair+window \u2248 geometric',
                     transform=ax1.transAxes, fontsize=8,
                     style='italic', alpha=0.6)

        ax2 = axes[1]
        excl = np.array(scan['exclusion'])
        excl_nc = np.array(scan['exclusion_no_cuts'])
        # Plot geometric-only first (behind), then pair+window on top
        ax2.loglog(lifetimes * SPEED_OF_LIGHT, excl_nc,
                   color='blue', linewidth=2, linestyle='--', alpha=0.5,
                   label='PX56 (geometric only)')
        ax2.loglog(lifetimes * SPEED_OF_LIGHT, excl,
                   color='blue', linewidth=2,
                   label='PX56 (pair+window)')
        ax2.set_xlabel(r'$c\tau$ (m)')
        ax2.set_ylabel(exclusion_ylabel(llp_pdg_id))
        ax2.set_title('Exclusion sensitivity (3 signal events)')
        ax2.set_ylim(1e-6, 1)
        ax2.grid(True, which="both", ls="-", alpha=0.2)
        # Annotate if curves overlap
        _ratio2 = excl / np.clip(excl_nc, 1e-30, None)
        _valid2 = _ratio2[(excl > 0) & (excl_nc > 0)]
        if len(_valid2) > 0 and np.median(_valid2) < 1.05:
            ax2.text(0.05, 0.05,
                     'Curves overlap:\npair+window \u2248 geometric',
                     transform=ax2.transAxes, fontsize=8,
                     style='italic', alpha=0.6)

        if args.external_dir:
            overlay_mass_matched_external_curves(ax2, sample_mass, args.external_dir,
                                                 llp_pdg_id)

        ax2.legend(fontsize=8, loc='upper right')

        mass_str = f'{sample_mass:.1f} GeV' if sample_mass is not None else '?'
        fig.suptitle(f'{output_tag}  |  $m = {mass_str}$  |  '
                     f'$\\sigma = {args.xsec:.3g}$ fb  |  '
                     f'$L = {args.lumi:.0f}$ fb$^{{-1}}$',
                     fontsize=12)
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        exclusion_plot_name = f"exclusion_Ntrack_{output_tag}.png"
        exclusion_plot_path = os.path.join(image_dir, exclusion_plot_name)
        plt.savefig(exclusion_plot_path, dpi=150)
        plt.close(fig)

        print(f"\nPlot saved: {exclusion_plot_path}")
        print(f"Normalization: n_generated={scan['n_generated']}, "
              f"n_llp_events={scan['n_llp_events']} "
              f"(0-LLP fraction: "
              f"{1 - scan['n_llp_events']/scan['n_generated']:.1%})")

        hist_lifetime_s = args.hist_lifetime_ns * 1e-9

    # === Separation histogram using best pair at chosen lifetime ===
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

    hist_lifetime_ns = hist_lifetime_s * 1e9

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
        ax.set_title(f'Selected best-pair separation (tau = {hist_lifetime_ns:.2e} ns)')
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
        mom_arr = momenta[mask_finite]
        sep_cm_arr = seps[mask_finite] * 100
        p99_mom = np.percentile(mom_arr, 99) if len(mom_arr) > 0 else 500
        p99_sep = np.percentile(sep_cm_arr, 99) if len(sep_cm_arr) > 0 else 500
        mom_bins = np.linspace(0, p99_mom * 1.1, 50)
        sep_bins = np.linspace(0, p99_sep * 1.1, 50)
        h = ax3.hist2d(mom_arr, sep_cm_arr,
                       bins=[mom_bins, sep_bins],
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
            print(f"Separation summary (selected best-pair @ {hist_lifetime_ns:.2e} ns):")
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
