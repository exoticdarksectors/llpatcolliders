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

M_DAUGHTER = 0.10566  # GeV/c² (muon mass, matching generator decay a → μ⁺μ⁻)

# Analysis cuts
P_CUT   = 0.600    # GeV/c — minimum daughter momentum
SEP_MIN = 0.001    # m — minimum separation at detector (1 mm)
SEP_MAX = 1.0     # m — maximum separation at detector

# ============================================================
# Two-body decay acceptance (analytical)
# ============================================================
#
# LLP (mass M, momentum p_LLP) → μ+ μ-
#
# Rest frame: E* = M/2,  p* = sqrt(M²/4 - m_μ²)
# Isotropic in cosθ*, φ*
#
# Lab frame (boost along LLP direction, m_μ → 0 limit):
#   E_{1,2}  = γ M/2 (1 ± β cosθ*)
#   |p|_{1,2} ≈ E_{1,2}  (since E >> m_μ)
#
# Opening angle:
#   cos(θ_12) = 1 - 2 / [γ²(1 - β² cos²θ*)]
#
# Accepted region of |cosθ*|:
#
#   c_S(d)      < |cosθ*| < min(c_P, c_max_sep(d), β)
#   ↑ min sep                 ↑ momentum cut  ↑ max sep  ↑ forward
#
# Accepted fraction: A(d) = max(0, c_upper(d) - c_S(d))
#
# Detectable decay probability:
#   P = ∫_{d_entry}^{d_exit} (1/λ) e^{-d/λ} · A(d) dd
#

def compute_c_upper(gamma, beta, mass, p_cut=P_CUT):
    """
    Upper limit on |cosθ*| from momentum cut and forward requirement.
    
    Momentum cut on softer muon:
      |p_2| > p_cut  →  E_2 > sqrt(p_cut² + m_μ²)
      γ M/2 (1 - β|cosθ*|) > sqrt(p_cut² + m_μ²)
      |cosθ*| < (1 - 2·sqrt(p_cut² + m_μ²) / (γM)) / β

    Forward requirement: both muons forward → |cosθ*| < β
    """
    E_min = np.sqrt(p_cut**2 + M_DAUGHTER**2)
    c_P = (1.0 - 2.0 * E_min / (gamma * mass)) / beta
    c_P = min(c_P, 1.0)
    return min(c_P, beta)


def compute_c_S(d_remaining, gamma, beta, sep_min=SEP_MIN):
    """
    Minimum |cosθ*| from minimum separation requirement.
    
    Separation ≈ θ_12 × d_remaining > sep_min
    Need θ_12 > θ_min → need |cosθ*| > c_S
    
    From cos(θ_12) = 1 - 2/(γ²(1 - β²cos²θ*)):
      θ_12 > θ_min  ↔  cos(θ_12) < cos(θ_min)
      β²cos²θ* > 1 - 2/(γ²(1 - cos(θ_min)))
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
    Maximum |cosθ*| from maximum separation requirement.
    
    Separation ≈ θ_12 × d_remaining < sep_max
    Need θ_12 < θ_max → need |cosθ*| < c_max_sep
    
    From cos(θ_12) = 1 - 2/(γ²(1 - β²cos²θ*)):
      θ_12 < θ_max  ↔  cos(θ_12) > cos(θ_max)
      β²cos²θ* < 1 - 2/(γ²(1 - cos(θ_max)))
    """
    if d_remaining <= 0:
        return 1.0  # No constraint at exit point
    
    theta_max = sep_max / d_remaining
    if theta_max >= np.pi:
        return 1.0  # No constraint
    
    cos_theta_max = np.cos(theta_max)
    denom = gamma**2 * (1.0 - cos_theta_max)
    if denom <= 0:
        return 1.0
    
    val = (1.0 - 2.0 / denom) / beta**2
    if val <= 0:
        return 0.0  # Even cosθ*=0 exceeds max separation
    if val >= 1.0:
        return 1.0  # No constraint
    return np.sqrt(val)


def acceptance_weighted_decay_prob(entry_d, exit_d, gamma, beta, mass,
                                    decay_length, p_cut=P_CUT,
                                    sep_min=SEP_MIN, sep_max=SEP_MAX):
    """
    Compute decay probability weighted by two-body decay acceptance.
    
    P = ∫_{entry}^{exit} (1/λ) e^{-d/λ} · A(d) dd
    
    A(d) = max(0, c_upper(d) - c_lower(d))
    where:
      c_lower = c_S(d)                          [min separation]
      c_upper = min(c_P, c_max_sep(d), β)       [momentum + max sep + forward]
    """
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
                               xsec_fb=60E3, lumi_fb=3000,
                               n_generated=None):
    df_base = pd.read_csv(csv_file)
    df_base.columns = df_base.columns.str.strip()
    n_llp_events = df_base['event'].nunique()

    # n_generated = total Pythia events (including those with 0 LLPs).
    # Events with 0 LLPs contribute P(≥1 decay)=0 to the average.
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
        # contribute P(≥1)=0, so sum(P1) / n_generated is correct.
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
    a distribution of muon-pair separations at the detector.
    
    For each particle that hits the fiducial volume:
      1. Sample decay position d from (1/λ) exp(-d/λ) within [entry, exit]
         using inverse CDF sampling
      2. Sample |cosθ*| uniformly in [0, 1] (isotropic decay)
      3. Compute opening angle θ_12 from cos(θ_12) = 1 - 2/(γ²(1-β²cos²θ*))
      4. Compute separation = θ_12 × d_remaining
    
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
        # CDF: F(d) = [exp(-entry/λ) - exp(-d/λ)] / [exp(-entry/λ) - exp(-exit/λ)]
        # Inverse: d = -λ ln(exp(-entry/λ) - u·[exp(-entry/λ) - exp(-exit/λ)])
        u = rng.uniform(0, 1, n_samples_per_particle)
        exp_entry = np.exp(-entry / decay_length)
        exp_exit = np.exp(-exit_ / decay_length)
        denom = exp_entry - exp_exit
        if denom < 1e-300:
            continue  # Negligible decay probability
        
        d_samples = -decay_length * np.log(exp_entry - u * denom)
        d_remaining = exit_ - d_samples
        
        # Weight: overall probability that the particle decays in the fiducial volume
        # Each sample gets equal weight = P_decay / n_samples
        p_decay = exp_entry * (1 - np.exp(-path_length / decay_length))
        w = p_decay / n_samples_per_particle
        
        # Sample |cosθ*| uniformly in [0, 1]
        cos_theta_star = rng.uniform(0, 1, n_samples_per_particle)
        
        # Compute opening angle
        # cos(θ_12) = 1 - 2/(γ²(1 - β²cos²θ*))
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
    import argparse
    import os
    import sys
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", nargs="?", default="LLP.csv")
    parser.add_argument("--xsec", type=float, default=60E3,
                        help="production cross-section in fb "
                             "(heavy ALP: ~60e3 for pp→h; "
                             "light ALP: ~373e6 for inclusive pp→bb̄)")
    parser.add_argument("--lumi", type=float, default=3000,
                        help="integrated luminosity in fb⁻¹ (default: 3000)")
    parser.add_argument("--outdir", default="output",
                        help="output directory for plots and CSVs (default: output)")
    parser.add_argument("--n-events", type=int, default=None,
                        help="total generated events (including 0-LLP events). "
                             "Auto-read from <csv>_meta.json if available.")
    args = parser.parse_args()
    xsec_arg_given = any(
        tok == "--xsec" or tok.startswith("--xsec=")
        for tok in sys.argv[1:]
    )
    if not xsec_arg_given:
        print("WARNING: --xsec not provided; using default 60e3 fb (pp→h).")
        print("         For light ALP (pp→bb̄): use --xsec 373e6.")
        print("         Always pass benchmark-specific --xsec explicitly.")
    os.makedirs(args.outdir, exist_ok=True)
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
        print(f"Read n_generated={n_generated}, llp_pdg_id={llp_pdg_id} from {meta_path}")
    else:
        print(f"WARNING: no _meta.json found ({meta_path}).")
        print("         Using CSV event count (0-LLP events missing!).")
        print("         Pass --n-events or re-run production to fix.")
    # Fallback: infer PDG ID from CSV data if meta not available
    if llp_pdg_id is None:
        _tmp = pd.read_csv(sample_csv, nrows=1)
        _tmp.columns = _tmp.columns.str.strip()
        llp_pdg_id = int(abs(_tmp['id'].iloc[0]))
        print(f"Inferred llp_pdg_id={llp_pdg_id} from CSV")
    separation_plot_name = f"separation_histogram_{output_tag}.png"
    exclusion_plot_name = f"exclusion_2body_{output_tag}.png"
    origin = [0, 0, 0]
    
    geo_cache = cache_geometry(sample_csv, mesh_fiducial, origin)
    
    # --- Diagnostics ---
    print("\n" + "="*50)
    print("ACCEPTANCE DIAGNOSTICS")
    print("="*50)
    
    test_mass = geo_cache['mass'][0]
    print(f"\nM = {test_mass:.2f} GeV → μ⁺μ⁻, p_cut = {P_CUT*1000:.0f} MeV/c, "
          f"sep: {SEP_MIN*1000:.0f} mm – {SEP_MAX*100:.0f} cm")
    print(f"\n{'p_LLP':>8} {'γ':>6} {'c_P':>7} {'θ_min':>10} "
          f"{'d(1mm)':>8} {'d(10cm)':>9} {'θ_12=10cm/1m':>14}")

    for p_test in [20, 50, 100, 200, 500]:
        E = np.sqrt(p_test**2 + test_mass**2)
        g = E / test_mass
        b = p_test / E
        c_P = compute_c_upper(g, b, test_mass)
        theta_min = 2.0 / g
        d_min_sep = SEP_MIN / theta_min   # distance where min sep = θ_min
        d_max_sep = SEP_MAX / theta_min   # distance where max sep = θ_min
        # At 1m remaining, what is max sep from θ_12(cosθ*=0)?
        sep_at_1m = theta_min * 1.0
        print(f"{p_test:>8.0f} {g:>6.1f} {c_P:>7.4f} {theta_min*1000:>8.1f} mrad "
              f"{d_min_sep*100:>7.2f}cm {d_max_sep:>7.1f}m {sep_at_1m*100:>12.1f}cm")

    # Show how max sep constrains high-γ particles
    g_lo = np.sqrt(20**2 + test_mass**2) / test_mass
    g_hi = np.sqrt(500**2 + test_mass**2) / test_mass
    print(f"\nMax sep cut effect: at d_remaining = 1m, θ_max = {SEP_MAX}m / 1m = 100 mrad")
    print("  → kills decays far from detector for low-γ (large opening angle)")
    print(f"  → for γ={g_lo:.1f} (p=20 GeV): θ_min={2/g_lo*1000:.0f} mrad "
          f"{'>> 100 mrad → heavily constrained' if 2/g_lo*1000 > 100 else '<= 100 mrad → mostly unconstrained'}")
    print(f"  → for γ={g_hi:.1f} (p=500 GeV): θ_min={2/g_hi*1000:.0f} mrad "
          f"{'>> 100 mrad → heavily constrained' if 2/g_hi*1000 > 100 else '<= 100 mrad → mostly unconstrained'}")
    
    # Single lifetime
    print("\n" + "="*50)
    print("SINGLE LIFETIME ANALYSIS")
    print("="*50)
    
    lifetime = 100e-8
    df_results, event_df = process_with_acceptance(
        sample_csv, lifetime, geo_cache)
    
    hits = df_results[df_results['hits_tube']]
    print(f"\nτ = {lifetime*1e9:.1f} ns:")
    print(f"  Particles hitting fiducial: {len(hits)}")
    if len(hits) > 0:
        print(f"  Mean acceptance: {hits['acceptance'].mean():.4f}")
        print(f"  Mean P_decay (with cuts):    {hits['decay_probability'].mean():.6f}")
        print(f"  Mean P_decay (without cuts): {hits['decay_probability_no_cuts'].mean():.6f}")
    n_llp_events = event_df.shape[0]
    n_denom = n_generated if n_generated is not None else n_llp_events
    mean_p1_all = event_df['prob_at_least_one_decays'].sum() / n_denom
    mean_p1_nc_all = event_df['prob_at_least_one_no_cuts'].sum() / n_denom
    print(f"  Events with LLPs: {n_llp_events} / {n_denom} generated")
    print(f"  Mean P(≥1) with cuts:    {mean_p1_all:.6e}")
    print(f"  Mean P(≥1) without cuts: {mean_p1_nc_all:.6e}")
    
    # --- Separation histogram ---
    print("\n" + "="*50)
    print("SEPARATION DISTRIBUTION")
    print("="*50)
    
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
        ax.set_title(f'Separation distribution (τ = {lifetime*1e9:.0f} ns)\n'
                     f'All decays, weighted by P(decay)')
        ax.legend(fontsize=9)
        ax.set_xlim(0, lin_xmax)
        
        # (b) Log-scale zoom to see the tails and cut region
        ax2 = axes_sep[1]
        bins_log = np.logspace(-4, 1, 100)  # 0.1 mm to 10 m
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
        plt.savefig(os.path.join(args.outdir, separation_plot_name), dpi=150)
        
        # Print summary
        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        frac_accepted = weights[in_window].sum() / weights.sum()
        print(f"  Fraction of decays in separation window: {frac_accepted:.3f}")
        print(f"  Median separation (all decays): {np.median(seps)*100:.1f} cm")
        accepted_seps = seps[in_window]
        if len(accepted_seps) > 0:
            print(f"  Median separation (accepted):  {np.median(accepted_seps)*100:.1f} cm")

    # Lifetime scan
    print("\n" + "="*50)
    print("LIFETIME SCAN")
    print("="*50)
    print(f"  σ = {args.xsec:.3g} fb,  L = {args.lumi:.0f} fb⁻¹")
    print(f"  n_generated = {n_generated if n_generated is not None else '(unknown — using CSV events only)'}")

    lifetimes = np.logspace(-9.5, -4.5, 20)
    scan = analyze_decay_vs_lifetime(sample_csv, geo_cache, lifetimes,
                                     xsec_fb=args.xsec, lumi_fb=args.lumi,
                                     n_generated=n_generated)
    
    # === Plotting ===
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
    ax2.set_ylabel('Mean P(≥1 decay)')
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
    
    # External comparison curves — model depends on production channel.
    # PDG 9000001 = light ALP (B→K(*)a): use CODEX B→KS curves.
    # PDG 6000113 = heavy ALP (h→aa):    use dark-Higgs H(125)→SS curves
    #                                     (valid comparison, same production).
    LIGHT_ALP_PDG = 9000001
    # BKS curves available at m = 0.5, 1.0, 2.0, 3.0 GeV
    _bks_masses = {0.5: "external/BKS/CODEX_BKS_m05.csv",
                   1.0: "external/BKS/CODEX_BKS_m1.csv",
                   2.0: "external/BKS/CODEX_BKS_m2.csv",
                   3.0: "external/BKS/CODEX_BKS_m3.csv"}

    if llp_pdg_id == LIGHT_ALP_PDG:
        # Pick the nearest available mass
        alp_mass = geo_cache['mass'][0]
        nearest_m = min(_bks_masses, key=lambda m: abs(m - alp_mass))
        bks_path = _bks_masses[nearest_m]
        print(f"  Light ALP (PDG {LIGHT_ALP_PDG}): overlaying CODEX B→KS "
              f"curve at m={nearest_m} GeV (ALP mass={alp_mass:.2f} GeV)")
        try:
            data = np.loadtxt(bks_path, delimiter=",")
            ax4.loglog(data[:, 0], data[:, 1],
                       color="cyan", linewidth=2, linestyle="-",
                       label=f"CODEX-b B→KS (m={nearest_m} GeV)")
        except (FileNotFoundError, OSError):
            print(f"  Note: BKS curve '{bks_path}' not found, skipping.")
        ax4.set_ylabel(r'BR$(B \to K^{(*)} a)_{\min}$')
    else:
        # Heavy ALP: dark-Higgs curves (H(125)→SS, m_S ~ 1 and 15 GeV).
        # TODO P5: confirm per-file mass and only overlay mass-matched subset.
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
                ax4.loglog(data[:, 0], data[:, 1],
                           color=color, linewidth=2, linestyle=ls, label=label)
            except (FileNotFoundError, OSError):
                print(f"  Note: external curve '{path}' not found, skipping.")
        ax4.set_ylabel(r'BR$(h \to aa)_{\min}$')

    ax4.legend(fontsize=8, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, exclusion_plot_name), dpi=150)

    df_results.to_csv(os.path.join(args.outdir, "particle_decay_results_2body.csv"), index=False)
    event_df.to_csv(os.path.join(args.outdir, "event_decay_statistics_2body.csv"), index=False)
    print("\nResults saved to", args.outdir)
    print(f"Plots: {exclusion_plot_name}, {separation_plot_name}")
    print(f"\nNormalization: n_generated={scan['n_generated']}, "
          f"n_llp_events={scan['n_llp_events']} "
          f"(0-LLP fraction: {1 - scan['n_llp_events']/scan['n_generated']:.1%})")
