"""
Two-body decay acceptance analysis for LLP sensitivity studies.

Computes acceptance-weighted decay probabilities for LLP -> 2-body decays
(e.g. a -> mu+ mu-) in the PX56 drainage tunnel detector, scanning over
proper lifetime to produce exclusion curves.

Physics model:
  - LLP (mass M, momentum p) decays isotropically in rest frame
  - Two-body kinematics with finite daughter mass
  - Acceptance window: daughter momentum cut, min/max separation at detector
  - Decay probability integrated over fiducial path with acceptance weighting

Usage:
  python decayProbPerEvent_2body.py <csv_file> --xsec <fb> [options]
"""

import os
import sys
import json
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'geometry'))
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR = os.path.join(_SCRIPT_DIR, '..')

from utils import infer_sample_mass, overlay_mass_matched_external_curves, exclusion_ylabel

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import quad
from gargoyle_geometry import (
    SPEED_OF_LIGHT,
    calculate_decay_length, cache_geometry,
    mesh_fiducial,
)

M_DAUGHTER = 0.10566  # GeV/c^2 (muon mass)

# Analysis cuts
P_CUT   = 0.600    # GeV/c -- minimum daughter momentum
SEP_MIN = 0.001    # m -- minimum separation at detector (1 mm)
SEP_MAX = 10.      # m -- maximum separation at detector

# Default lifetime scan range (extended for full U-curve)
DEFAULT_LIFETIME_MIN_NS = 10**-2
DEFAULT_LIFETIME_MAX_NS = 10**5.5
DEFAULT_LIFETIME_POINTS = 80


# ============================================================
# Two-body decay acceptance (analytical)
# ============================================================
#
# LLP (mass M, momentum p_LLP) -> daughter+ daughter-
#
# Rest frame: E* = M/2,  p* = sqrt(M^2/4 - m_daughter^2)
# Isotropic in cos(theta*), phi*
#
# Lab frame (boost along LLP direction):
#   E_{1,2}  = gamma M/2 (1 +/- beta cos(theta*))
#   |p|_{1,2} ~ E_{1,2}  (since E >> m_daughter)
#
# Opening angle:
#   cos(theta_12) = 1 - 2 / [gamma^2 (1 - beta^2 cos^2(theta*))]
#
# Accepted region of |cos(theta*)|:
#   c_S(d) < |cos(theta*)| < min(c_P, c_max_sep(d), beta)
#

def compute_c_upper(gamma, beta, mass, p_cut=P_CUT):
    """
    Upper limit on |cos(theta*)| from momentum cut and forward requirement.

    Momentum cut on softer daughter:
      |p_2| > p_cut  ->  E_2 > sqrt(p_cut^2 + m_daughter^2)
      gamma M/2 (1 - beta|cos(theta*)|) > sqrt(p_cut^2 + m_daughter^2)
      |cos(theta*)| < (1 - 2*sqrt(p_cut^2 + m_daughter^2) / (gamma*M)) / beta

    Forward requirement: both daughters forward -> |cos(theta*)| < beta
    """
    E_min = np.sqrt(p_cut**2 + M_DAUGHTER**2)
    c_P = (1.0 - 2.0 * E_min / (gamma * mass)) / beta
    c_P = min(c_P, 1.0)
    return min(c_P, beta)


def compute_c_S(d_remaining, gamma, beta, sep_min=SEP_MIN):
    """
    Minimum |cos(theta*)| from minimum separation requirement.

    Separation ~ theta_12 * d_remaining > sep_min
    Need theta_12 > theta_min -> need |cos(theta*)| > c_S

    From cos(theta_12) = 1 - 2/(gamma^2 (1 - beta^2 cos^2(theta*))):
      theta_12 > theta_min  <->  cos(theta_12) < cos(theta_min)
      beta^2 cos^2(theta*) > 1 - 2/(gamma^2 (1 - cos(theta_min)))
    """
    if d_remaining <= 0:
        return 1.0

    theta_min = sep_min / d_remaining
    if theta_min <= 0:
        return 0.0

    cos_theta_min = np.cos(min(theta_min, np.pi))
    denom = gamma**2 * (1.0 - cos_theta_min)
    if denom <= 0:
        return 0.0

    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0
    return np.sqrt(val)


def compute_c_max_sep(d_remaining, gamma, beta, sep_max=SEP_MAX):
    """
    Maximum |cos(theta*)| from maximum separation requirement.

    Separation ~ theta_12 * d_remaining < sep_max
    Need theta_12 < theta_max -> need |cos(theta*)| < c_max_sep

    From cos(theta_12) = 1 - 2/(gamma^2 (1 - beta^2 cos^2(theta*))):
      theta_12 < theta_max  <->  cos(theta_12) > cos(theta_max)
      beta^2 cos^2(theta*) < 1 - 2/(gamma^2 (1 - cos(theta_max)))
    """
    if d_remaining <= 0:
        return 1.0

    theta_max = sep_max / d_remaining
    if theta_max >= np.pi:
        return 1.0

    cos_theta_max = np.cos(theta_max)
    denom = gamma**2 * (1.0 - cos_theta_max)
    if denom <= 0:
        return 1.0

    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0
    if val >= 1.0:
        return 1.0
    return np.sqrt(val)


def acceptance_weighted_decay_prob(entry_d, exit_d, gamma, beta, mass,
                                    decay_length, p_cut=P_CUT,
                                    sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    Compute decay probability weighted by two-body decay acceptance.

    P = integral_{entry}^{exit} (1/lambda) e^{-d/lambda} * A(d) dd

    A(d) = max(0, c_upper(d) - c_lower(d))
    where:
      c_lower = c_S(d)                          [min separation]
      c_upper = min(c_P, c_max_sep(d), beta)    [momentum + max sep + forward]
    """
    if mass < 2 * M_DAUGHTER:
        return 0.0

    c_P_upper = compute_c_upper(gamma, beta, mass, p_cut)

    if c_P_upper <= 0:
        return 0.0

    def integrand(d):
        d_remaining = exit_d - d
        c_lower = compute_c_S(d_remaining, gamma, beta, sep_min)
        c_upper_sep = compute_c_max_sep(d_remaining, gamma, beta, sep_max)
        c_upper = min(c_P_upper, c_upper_sep)
        A = max(0.0, c_upper - c_lower)
        return (1.0 / decay_length) * np.exp(-d / decay_length) * A

    result, _ = quad(integrand, entry_d, exit_d, limit=100,
                     epsabs=1e-15, epsrel=1e-10)
    return result


def unweighted_decay_prob(entry_d, exit_d, decay_length):
    """Decay probability with no acceptance cuts."""
    path_length = exit_d - entry_d
    return np.exp(-entry_d / decay_length) * (1.0 - np.exp(-path_length / decay_length))


def process_with_acceptance(csv_file_or_df, lifetime_seconds, geo_cache,
                             p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX):
    if isinstance(csv_file_or_df, pd.DataFrame):
        df = csv_file_or_df.copy()
    else:
        df = pd.read_csv(csv_file_or_df)
        df.columns = df.columns.str.strip()
    hits = geo_cache['hits']

    df['hits_tube'] = hits
    df['decay_probability'] = 0.0
    df['decay_probability_no_cuts'] = 0.0
    df['acceptance'] = 0.0

    for idx in np.where(hits)[0]:
        p = geo_cache['momentum'][idx]
        m = geo_cache['mass'][idx]
        gamma = geo_cache['gamma'][idx]
        beta = geo_cache['beta'][idx]
        entry = geo_cache['entry_d'][idx]
        exit_ = geo_cache['exit_d'][idx]
        decay_length = calculate_decay_length(p, m, lifetime_seconds)

        p_with = acceptance_weighted_decay_prob(
            entry, exit_, gamma, beta, m, decay_length, p_cut, sep_min, sep_max)
        p_no = unweighted_decay_prob(entry, exit_, decay_length)

        df.iat[idx, df.columns.get_loc('decay_probability')] = p_with
        df.iat[idx, df.columns.get_loc('decay_probability_no_cuts')] = p_no
        if p_no > 0:
            df.iat[idx, df.columns.get_loc('acceptance')] = p_with / p_no

    event_stats = []
    for event_idx, group in df.groupby('event'):
        n_part = len(group)
        p1 = group.iloc[0]['decay_probability'] if n_part >= 1 else 0
        p2 = group.iloc[1]['decay_probability'] if n_part >= 2 else 0
        p1_nc = group.iloc[0]['decay_probability_no_cuts'] if n_part >= 1 else 0
        p2_nc = group.iloc[1]['decay_probability_no_cuts'] if n_part >= 2 else 0

        p_at_least_one = 1 - (1 - p1) * (1 - p2)
        p_at_least_one_nc = 1 - (1 - p1_nc) * (1 - p2_nc)

        event_stats.append({
            'event': event_idx,
            'n_particles_hitting_tube': group['hits_tube'].sum(),
            'prob_at_least_one_decays': p_at_least_one,
            'prob_at_least_one_no_cuts': p_at_least_one_nc,
            'prob_both_decay': p1 * p2,
            'particle1_decay_prob': p1,
            'particle2_decay_prob': p2,
        })

    event_df = pd.DataFrame(event_stats)
    return df, event_df


def analyze_decay_vs_lifetime(csv_file, geo_cache, lifetime_range,
                               p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                               xsec_fb=54700, lumi_fb=3000,
                               n_generated=None):
    df_base = pd.read_csv(csv_file)
    df_base.columns = df_base.columns.str.strip()
    n_llp_events = df_base['event'].nunique()

    # n_generated = total Pythia events (including those with 0 LLPs).
    # Events with 0 LLPs contribute P(>=1 decay)=0 to the average.
    if n_generated is None:
        n_generated = n_llp_events
    if n_generated < n_llp_events:
        print(f"WARNING: n_generated ({n_generated}) < n_llp_events "
              f"({n_llp_events}); using n_llp_events as denominator.")
        n_generated = n_llp_events

    results = {
        'lifetimes': lifetime_range,
        'mean_single_particle_decay_prob': [],
        'mean_single_no_cuts': [],
        'mean_at_least_one_decay_prob': [],
        'mean_at_least_one_no_cuts': [],
        'mean_both_decay_prob': [],
        'exclusion': [],
        'exclusion_no_cuts': [],
        'mean_acceptance': [],
        'n_generated': n_generated,
        'n_llp_events': n_llp_events,
    }

    for lifetime in tqdm(lifetime_range, desc="Scanning lifetimes"):
        df, event_df = process_with_acceptance(
            df_base, lifetime, geo_cache, p_cut, sep_min, sep_max)

        hits = df[df['hits_tube']]
        mean_single = hits['decay_probability'].mean() if len(hits) > 0 else 0
        mean_single_nc = hits['decay_probability_no_cuts'].mean() if len(hits) > 0 else 0
        mean_acc = hits['acceptance'].mean() if len(hits) > 0 else 0

        # Average over ALL generated events: events without LLPs
        # contribute P(>=1)=0, so sum(P1) / n_generated is correct.
        mean_p1 = event_df['prob_at_least_one_decays'].sum() / n_generated
        mean_p1_nc = event_df['prob_at_least_one_no_cuts'].sum() / n_generated
        mean_both = event_df['prob_both_decay'].sum() / n_generated

        results['mean_single_particle_decay_prob'].append(mean_single)
        results['mean_single_no_cuts'].append(mean_single_nc)
        results['mean_at_least_one_decay_prob'].append(mean_p1)
        results['mean_at_least_one_no_cuts'].append(mean_p1_nc)
        results['mean_both_decay_prob'].append(mean_both)
        denom = mean_p1 * lumi_fb * xsec_fb
        denom_nc = mean_p1_nc * lumi_fb * xsec_fb
        results['exclusion'].append(3 / denom if denom > 0 else np.inf)
        results['exclusion_no_cuts'].append(3 / denom_nc if denom_nc > 0 else np.inf)
        results['mean_acceptance'].append(mean_acc)

    return results


def sample_separations(geo_cache, lifetime_seconds, n_samples_per_particle=100,
                       rng_seed=42):
    """
    Monte Carlo sample decay positions and rest-frame angles to build
    a distribution of daughter-pair separations at the detector.

    For each particle that hits the fiducial volume:
      1. Sample decay position d from (1/lambda) exp(-d/lambda) within [entry, exit]
         using inverse CDF sampling
      2. Sample |cos(theta*)| uniformly in [0, 1] (isotropic decay)
      3. Compute opening angle theta_12 from
         cos(theta_12) = 1 - 2/(gamma^2 (1 - beta^2 cos^2(theta*)))
      4. Compute separation = theta_12 * d_remaining

    Each sample is weighted by the exponential decay probability so the
    histogram represents the physical separation distribution.

    Returns:
        separations: array of separation values (m)
        weights: array of per-sample weights (decay probability contribution)
        momenta: array of parent LLP momentum for each sample
    """
    rng = np.random.default_rng(rng_seed)

    hits = geo_cache['hits']
    hit_idx = np.where(hits)[0]

    all_seps = []
    all_weights = []
    all_momenta = []

    for idx in hit_idx:
        entry = geo_cache['entry_d'][idx]
        exit_ = geo_cache['exit_d'][idx]
        gamma = geo_cache['gamma'][idx]
        beta = geo_cache['beta'][idx]
        mass = geo_cache['mass'][idx]
        p_llp = geo_cache['momentum'][idx]

        decay_length = calculate_decay_length(p_llp, mass, lifetime_seconds)
        path_length = exit_ - entry

        # Inverse CDF sampling of decay position within [entry, exit]
        u = rng.uniform(0, 1, n_samples_per_particle)
        exp_entry = np.exp(-entry / decay_length)
        exp_exit = np.exp(-exit_ / decay_length)
        denom = exp_entry - exp_exit
        if denom < 1e-300:
            continue

        d_samples = -decay_length * np.log(exp_entry - u * denom)
        d_remaining = exit_ - d_samples

        # Weight: overall probability that the particle decays in the fiducial volume
        p_decay = exp_entry * (1 - np.exp(-path_length / decay_length))
        w = p_decay / n_samples_per_particle

        # Sample |cos(theta*)| uniformly in [0, 1]
        cos_theta_star = rng.uniform(0, 1, n_samples_per_particle)

        # Compute opening angle
        denom_angle = gamma**2 * (1 - beta**2 * cos_theta_star**2)
        cos_theta_12 = 1 - 2.0 / denom_angle
        cos_theta_12 = np.clip(cos_theta_12, -1, 1)
        theta_12 = np.arccos(cos_theta_12)

        # Separation at detector
        sep = theta_12 * d_remaining

        all_seps.append(sep)
        all_weights.append(np.full(n_samples_per_particle, w))
        all_momenta.append(np.full(n_samples_per_particle, p_llp))

    if not all_seps:
        return np.array([]), np.array([]), np.array([])

    return (np.concatenate(all_seps),
            np.concatenate(all_weights),
            np.concatenate(all_momenta))


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Two-body decay acceptance analysis for LLP sensitivity.")
    parser.add_argument("csv_file", nargs="?", default="LLP.csv")
    parser.add_argument("--xsec", type=float, default=None,
                        help="production cross-section in fb (REQUIRED). "
                             "Examples: heavy ALP 54700, light ALP 373e6, "
                             "DP Higgs 54700")
    parser.add_argument("--lumi", type=float, default=3000,
                        help="integrated luminosity in fb^-1 (default: 3000)")
    parser.add_argument("--outdir", default="output",
                        help="output directory for plots and CSVs (default: output)")
    parser.add_argument("--n-events", type=int, default=None,
                        help="total generated events (including 0-LLP events). "
                             "Auto-read from <csv>_meta.json if available.")
    parser.add_argument("--external-dir", default=None,
                        help="directory containing external comparison curves "
                             "(default: None, no overlay)")
    parser.add_argument("--lifetime-min-ns", type=float,
                        default=DEFAULT_LIFETIME_MIN_NS,
                        help="minimum lifetime in ns for scan "
                             f"(default: {DEFAULT_LIFETIME_MIN_NS:.0e})")
    parser.add_argument("--lifetime-max-ns", type=float,
                        default=DEFAULT_LIFETIME_MAX_NS,
                        help="maximum lifetime in ns for scan "
                             f"(default: {DEFAULT_LIFETIME_MAX_NS:.1e})")
    parser.add_argument("--lifetime-points", type=int,
                        default=DEFAULT_LIFETIME_POINTS,
                        help="number of lifetime scan points "
                             f"(default: {DEFAULT_LIFETIME_POINTS})")
    args = parser.parse_args()

    if args.xsec is None:
        print("ERROR: --xsec is required. Example values:")
        print("  Heavy ALP (gg->h, 14 TeV):      --xsec 54700")
        print("  Light ALP (pp->bb, 14 TeV):      --xsec 373e6")
        print("  Dark photon (gg->h, 14 TeV):     --xsec 54700")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    data_dir = os.path.join(args.outdir, "data")
    image_dir = os.path.join(args.outdir, "images")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(image_dir, exist_ok=True)

    sample_csv = args.csv_file
    output_tag = os.path.basename(sample_csv).replace('.csv', '')

    # Resolve n_generated and llp_pdg_id: CLI > meta.json > CSV fallback
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
        print(f"WARNING: no _meta.json found ({meta_path}). "
              "Using CSV event count (0-LLP events missing). "
              "Pass --n-events or re-run production to fix.")

    if llp_pdg_id is None:
        _tmp = pd.read_csv(sample_csv, nrows=1)
        _tmp.columns = _tmp.columns.str.strip()
        llp_pdg_id = int(abs(_tmp['id'].iloc[0]))
        print(f"Inferred llp_pdg_id={llp_pdg_id} from CSV")

    print(f"\nInput: {sample_csv}")
    print(f"  xsec = {args.xsec:.3g} fb, lumi = {args.lumi:.0f} fb^-1")
    print(f"  n_generated = {n_generated}, llp_pdg_id = {llp_pdg_id}")
    print(f"  output -> {args.outdir}")

    separation_plot_name = f"separation_histogram_{output_tag}.png"
    exclusion_plot_name = f"exclusion_2body_{output_tag}.png"
    origin = [0, 0, 0]

    geo_cache = cache_geometry(sample_csv, mesh_fiducial, origin)

    # --- Single lifetime for separation histogram ---
    lifetime = 100e-8
    df_results, event_df = process_with_acceptance(
        sample_csv, lifetime, geo_cache)

    # --- Separation histogram ---
    seps, weights, momenta = sample_separations(
        geo_cache, lifetime, n_samples_per_particle=200)

    if len(seps) > 0:
        fig_sep, axes_sep = plt.subplots(1, 3, figsize=(18, 5))

        # (a) All separations, decay-probability weighted
        ax = axes_sep[0]
        lin_xmax = 1.05 * SEP_MAX
        bins = np.linspace(0, lin_xmax, 100)
        ax.hist(seps, bins=bins, weights=weights, color='steelblue',
                edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.axvline(SEP_MIN, color='red', linestyle='--', linewidth=2,
                   label=f'min sep = {SEP_MIN*1000:.0f} mm')
        ax.axvline(SEP_MAX, color='red', linestyle='-', linewidth=2,
                   label=f'max sep = {SEP_MAX*100:.0f} cm')
        ax.axvspan(SEP_MIN, SEP_MAX, color='green', alpha=0.1, label='Accepted')
        ax.set_xlabel('Separation at detector (m)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Separation distribution (tau = {lifetime*1e9:.0f} ns)\n'
                     f'All decays, weighted by P(decay)')
        ax.legend(fontsize=9)
        ax.set_xlim(0, lin_xmax)

        # (b) Log-scale zoom
        ax2 = axes_sep[1]
        bins_log = np.logspace(-4, 1, 100)
        ax2.hist(seps, bins=bins_log, weights=weights, color='steelblue',
                 edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.axvline(SEP_MIN, color='red', linestyle='--', linewidth=2,
                    label=f'min sep = {SEP_MIN*1000:.0f} mm')
        ax2.axvline(SEP_MAX, color='red', linestyle='-', linewidth=2,
                    label=f'max sep = {SEP_MAX*100:.0f} cm')
        ax2.axvspan(SEP_MIN, SEP_MAX, color='green', alpha=0.1, label='Accepted')
        ax2.set_xscale('log')
        ax2.set_xlabel('Separation at detector (m)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('Log-scale separation\n(showing full range)')
        ax2.legend(fontsize=9)

        # (c) Separation vs LLP momentum (2D)
        ax3 = axes_sep[2]
        mask_finite = np.isfinite(seps) & (seps > 0)
        h = ax3.hist2d(momenta[mask_finite], seps[mask_finite] * 100,
                       bins=[np.linspace(0, 500, 50), np.linspace(0, 500, 50)],
                       weights=weights[mask_finite],
                       cmap='viridis', cmin=1e-20)
        ax3.axhline(SEP_MIN * 100, color='red', linestyle='--', linewidth=2,
                    label=f'min = {SEP_MIN*1000:.0f} mm')
        ax3.axhline(SEP_MAX * 100, color='red', linestyle='-', linewidth=2,
                    label=f'max = {SEP_MAX*100:.0f} cm')
        ax3.set_xlabel('LLP momentum (GeV/c)')
        ax3.set_ylabel('Separation at detector (cm)')
        ax3.set_title('Separation vs LLP momentum')
        ax3.legend(fontsize=9, loc='upper right')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig(os.path.join(image_dir, separation_plot_name), dpi=150)

    # --- Lifetime scan ---
    lifetimes = np.logspace(np.log10(args.lifetime_min_ns * 1e-9),
                            np.log10(args.lifetime_max_ns * 1e-9),
                            args.lifetime_points)
    scan = analyze_decay_vs_lifetime(sample_csv, geo_cache, lifetimes,
                                     xsec_fb=args.xsec, lumi_fb=args.lumi,
                                     n_generated=n_generated)

    # Print lifetime scan summary
    exclusion_arr = np.array(scan['exclusion'])
    finite_mask = np.isfinite(exclusion_arr)
    if finite_mask.any():
        best_idx = np.argmin(exclusion_arr[finite_mask])
        best_lt = lifetimes[finite_mask][best_idx]
        best_br = exclusion_arr[finite_mask][best_idx]
        print(f"\nLifetime scan: {args.lifetime_points} points, "
              f"[{args.lifetime_min_ns:.0e}, {args.lifetime_max_ns:.1e}] ns")
        print(f"  Peak signal at tau = {best_lt*1e9:.2g} ns "
              f"(ctau = {best_lt*SPEED_OF_LIGHT:.2g} m)")
        print(f"  BR_min = {best_br:.2e}")
        print(f"  n_generated = {scan['n_generated']}, "
              f"n_llp_events = {scan['n_llp_events']} "
              f"(0-LLP fraction: {1 - scan['n_llp_events']/scan['n_generated']:.1%})")

    # === Plotting ===
    sample_mass = infer_sample_mass(geo_cache['mass'])
    if sample_mass is not None and sample_mass < 2 * M_DAUGHTER:
        print(f"WARNING: sample mass {sample_mass:.4f} GeV < 2*M_DAUGHTER "
              f"({2*M_DAUGHTER:.4f} GeV). Decay to mu+mu- is kinematically "
              f"forbidden. All acceptance-weighted probabilities will be zero.")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Acceptance vs lifetime
    ax1 = axes[0, 0]
    ax1.semilogx(lifetimes * 1e9, scan['mean_acceptance'],
                 'b-', linewidth=2)
    ax1.set_xlabel('Lifetime (ns)')
    ax1.set_ylabel('Mean Acceptance')
    ax1.set_title(f'Two-body acceptance\n(p > {P_CUT*1000:.0f} MeV/c, '
                  f'{SEP_MIN*1000:.0f} mm < sep < {SEP_MAX*100:.0f} cm)')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.set_ylim(0, 1.05)

    # Plot 2: With vs without cuts
    ax2 = axes[0, 1]
    ax2.loglog(lifetimes * 1e9, scan['mean_at_least_one_decay_prob'],
               'r-', linewidth=2, label='With cuts')
    ax2.loglog(lifetimes * 1e9, scan['mean_at_least_one_no_cuts'],
               'r--', linewidth=2, alpha=0.5, label='Without cuts')
    ax2.set_xlabel('Lifetime (ns)')
    ax2.set_ylabel('Mean P(>=1 decay)')
    ax2.set_title('Effect of analysis cuts')
    ax2.grid(True, which="both", ls="-", alpha=0.2)
    ax2.legend()

    # Plot 3: Single particle
    ax3 = axes[1, 0]
    ax3.loglog(lifetimes * 1e9, scan['mean_single_particle_decay_prob'],
               'b-', linewidth=2, label='With cuts')
    ax3.loglog(lifetimes * 1e9, scan['mean_single_no_cuts'],
               'b--', linewidth=2, alpha=0.5, label='Without cuts')
    ax3.set_xlabel('Lifetime (ns)')
    ax3.set_ylabel('Mean Decay Probability')
    ax3.set_title('Single Particle: With vs Without Cuts')
    ax3.grid(True, which="both", ls="-", alpha=0.2)
    ax3.legend()

    # Plot 4: Exclusion curves
    ax4 = axes[1, 1]
    ax4.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion'],
               color='blue', linewidth=2, label="PX56 (with cuts)")
    ax4.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion_no_cuts'],
               color='blue', linewidth=2, linestyle='--', alpha=0.5,
               label="PX56 (no cuts)")
    ax4.set_xlabel(r'$c\tau$ (m)')
    ax4.grid(True, which="both", ls="-", alpha=0.2)

    ax4.set_ylabel(exclusion_ylabel(llp_pdg_id))

    # External comparison curves (mass-matched overlay)
    if args.external_dir:
        overlay_mass_matched_external_curves(
            ax4, sample_mass, args.external_dir, llp_pdg_id)
    ax4.set_ylim(1e-6, 1)

    # Decay mode annotation
    if np.isclose(M_DAUGHTER, 0.10566, atol=0.001):
        decay_mode = r'$a \to \mu^+\mu^-$'
    else:
        decay_mode = f'$m_{{d}}={M_DAUGHTER:.4g}$ GeV'
    ax4.text(0.05, 0.95, decay_mode, transform=ax4.transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax4.legend(fontsize=8, loc='upper right')

    # Suptitle with sample info
    mass_str = f'{sample_mass:.1f} GeV' if sample_mass is not None else '?'
    fig.suptitle(f'{output_tag}  |  $m = {mass_str}$  |  '
                 f'{decay_mode}  |  '
                 f'$\\sigma = {args.xsec:.3g}$ fb  |  '
                 f'$L = {args.lumi:.0f}$ fb$^{{-1}}$',
                 fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(image_dir, exclusion_plot_name), dpi=150)

    df_results.to_csv(os.path.join(data_dir,
                      "particle_decay_results_2body.csv"), index=False)
    event_df.to_csv(os.path.join(data_dir,
                    "event_decay_statistics_2body.csv"), index=False)
    print(f"\nPlots: {image_dir}/{exclusion_plot_name}, "
          f"{image_dir}/{separation_plot_name}")
    print(f"Data:  {data_dir}/particle_decay_results_2body.csv, "
          f"event_decay_statistics_2body.csv")
