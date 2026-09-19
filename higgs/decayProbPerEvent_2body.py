import os
import sys
import numpy as np
import pandas as pd
import matplotlib
# Batch mode: use the non-interactive Agg backend unless the user has set
# MPLBACKEND or passes --interactive. Default is batch so the script runs
# on headless machines without a DISPLAY.
INTERACTIVE = ('--interactive' in sys.argv)
if not INTERACTIVE and not os.environ.get('MPLBACKEND'):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import quad


def show_or_close():
    """plt.show() in interactive mode, plt.close('all') in batch mode."""
    if INTERACTIVE:
        plt.show()
    else:
        plt.close('all')

from grendel_geometry import (
    SPEED_OF_LIGHT, DETECTOR_THICKNESS,
    calculate_decay_length, cache_geometry, mesh_fiducial,
    points_on_tracker, classify_points_with_basis,
    path_3d_fiducial, tunnel_profile_points,
)
from reco_common import SIGMA_T_DEFAULT, CHI2_TIMING_MAX

M_ELECTRON = 0.000511  # GeV/c²

# Signal normalisation. Inclusive SM Higgs production at sqrt(s) = 14 TeV,
# m_H = 125.09 GeV, from the LHC Higgs XS WG (YR4): ggF 54.61 + VBF 4.275 +
# WH 1.510 + ZH 0.984 + ttH 0.613 + bbH 0.596 + tH 0.112 pb.
SIGMA_H_FB = 62.7e3    # fb — inclusive sigma(pp -> H) at 14 TeV
LUMI_FB    = 3000.0    # fb^-1 — HL-LHC integrated luminosity
N_SIG_EXCL = 3.0       # signal events for the quoted exclusion

# Analysis cuts
P_CUT   = 0.100    # GeV/c — minimum electron momentum (track must be energetic
                   # enough to reconstruct as a straight track despite multiple
                   # scattering in the scintillator; ~2 deg Highland angle at 100 MeV)
SEP_MIN = 0.01    # m — minimum separation at detector (1 cm)
SEP_MAX = 10.0      # m — maximum separation at detector
DCA_CUT = 0.1  # m (10 cm) — maximum DCA between reconstructed tracks
# Conditional max-sep_outer cut, motivated by the geometric statement that
# a real vertex satisfies sep_outer = open_angle × (D + L) with D bounded
# by the fiducial path length. For near-parallel tracks (open_angle at the
# angular-resolution noise floor), the largest sep_outer a real vertex can
# produce is ~θ_noise × max_path ≈ 30 mrad × 50 m ≈ 1.5 m. Above that,
# we are looking at two physically separated parallel particles
# (e.g., beam-halo or cosmic coincidences), not a vertex.
# When open_angle is clearly above the noise floor (diverging tracks),
# no upper bound on sep_outer is applied.
THETA_PARALLEL      = 0.100  # rad (100 mrad — covers parallel + moderate-opening regime)
SEP_OUT_MAX_PARALLEL = 0.30  # m — max sep_outer applied only when open_angle < THETA_PARALLEL

# IP-muon-transit veto. A single muon shooting outward from the CMS IP makes
# four hits that lie on a single straight line. The "collinearity" is the RMS
# perpendicular residual of all 4 hits to their best-fit line (SVD of the 4×3
# centered hit matrix; sqrt((S[1]^2 + S[2]^2) / 4)). For a real two-body decay
# the 4 hits trace a "V", so collinearity ≈ (sep_in + sep_out)/4 — large. For
# a muon transit it sits at the hit-noise floor (~σ_hit). To protect collimated
# signal (where the V is narrow and naturally looks like a line), the cut only
# fires when sep_outer > SEP_OUT_COLLIN_GATE, i.e. the topology is not in the
# collimated regime where signal and muon are geometrically indistinguishable.
COLLIN_MIN          = 0.030  # m (30 mm) — min collinearity for events above the gate
SEP_OUT_COLLIN_GATE = 0.30   # m — collinearity cut applies only when sep_outer > this

# L-scaled ("fractional") collinearity cut — the ACTIVE collinearity veto,
# shared with the cosmic-decay background (cosmic_decay_check imports COLLIN_FRAC
# from here). When sep_outer > SEP_OUT_GATE require collinearity > COLLIN_FRAC * L
# (L = DETECTOR_THICKNESS). Loosened to 0.40L: at the 0.5 ns timing default the
# cosmic decay-in-flight background is timing-handled, so the looser collinearity
# recovers high-mass signal (40 GeV) the tighter 0.48L would cost. COLLIN_MIN /
# SEP_OUT_COLLIN_GATE above are the older static version, kept for comparison.
COLLIN_FRAC = 0.40

# sep_outer gate for the collinearity veto: the cut fires when sep_outer >
# SEP_OUT_GATE. Lowered below L (= DETECTOR_THICKNESS = 24 cm) to 16 cm so the
# collinearity cut reaches the cosmic decay-in-flight "crack" (narrow events with
# sep_out just under L), at a small cost to collimated low-mass signal. Shared
# with the cosmic background (single source).
SEP_OUT_GATE = 0.16   # m (16 cm)

# Conditional tight pointing. For collimated decays (small sep_inner) the
# bisector points back to the IP very well, so a tight pointing cut kills
# background with little signal loss. For wider-sep_inner topologies the
# pointing resolution degrades and we leave it open.
POINT_TIGHT_SEP_IN  = 0.050  # rad (50 mrad) — max pointing when sep_in < gate
SEP_IN_POINT_GATE   = 0.10   # m (10 cm) — pointing cut applies only when sep_in < this

# Global pointing cut for the well-separated (sep_in >= 10 cm) topologies, on
# top of the tight 50 mrad applied to close (sep_in < 10 cm) decays above. This
# is the single-region selection: pointing < 50 mrad for close, < 800 mrad for
# well separated. 800 mrad captures the high-mass signal tail (15/40 GeV pointing
# reach ~660 mrad at 99%); the cosmic decay-in-flight background it admits is
# handled by the timing chi2 cut (default 0.5 ns). Shared with the cosmic
# background (single source).
POINT_GLOBAL = 0.8   # rad (800 mrad)

# Reconstructed vertex stand-off from the inner tracking layer. The
# reconstructed distance from the decay vertex to the inner tracking layer is
# d_implied = sep_inner / open_angle (the geometric identity
# sep_inner = open_angle × (vertex → inner-layer distance)). Requiring this to
# exceed VTX_INNER_MIN enforces a genuinely displaced vertex, standing off from
# the first instrumented layer. 0 disables the cut.
VTX_INNER_MIN = 0.0   # m (default off; 0.30 m = 30 cm stand-off)

outString = "15GeVFixANUBIS"
sample_csv = "LLPSmall.csv"

# Tracking resolution
HIT_RESOLUTION = 0.003  # m (3 mm per layer)
N_LAYERS       = 2       # number of tracking layers (stations)

# ============================================================
# Two-body decay acceptance (analytical)
# ============================================================
#
# LLP (mass M, momentum p_LLP) → e+ e-
#
# Rest frame: E* = M/2,  p* = sqrt(M²/4 - m_e²)
# Isotropic in cosθ*, φ*
#
# Lab frame (boost along LLP direction, m_e → 0 limit):
#   E_{1,2}  = γ M/2 (1 ± β cosθ*)
#   |p|_{1,2} ≈ E_{1,2}  (since E >> m_e)
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
    
    Momentum cut on softer electron:
      |p_2| > p_cut  →  E_2 > sqrt(p_cut² + m_e²)
      γ M/2 (1 - β|cosθ*|) > sqrt(p_cut² + m_e²)
      |cosθ*| < (1 - 2·sqrt(p_cut² + m_e²) / (γM)) / β
    
    Forward requirement: both electrons forward → |cosθ*| < β
    """
    E_min = np.sqrt(p_cut**2 + M_ELECTRON**2)
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


# ============================================================
# Processing functions
# ============================================================

def process_with_acceptance(csv_file, lifetime_seconds, geo_cache,
                             p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX):
    df = pd.read_csv(csv_file)
    n = len(df)
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
                               p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX):
    df_base = pd.read_csv(csv_file)
    n_events = df_base['event'].nunique()
    
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
        'total_events': n_events
    }
    
    for lifetime in tqdm(lifetime_range, desc="Scanning lifetimes"):
        df, event_df = process_with_acceptance(
            csv_file, lifetime, geo_cache, p_cut, sep_min, sep_max)
        
        hits = df[df['hits_tube']]
        mean_single = hits['decay_probability'].mean() if len(hits) > 0 else 0
        mean_single_nc = hits['decay_probability_no_cuts'].mean() if len(hits) > 0 else 0
        mean_acc = hits['acceptance'].mean() if len(hits) > 0 else 0
        
        mean_p1 = event_df['prob_at_least_one_decays'].mean()
        mean_p1_nc = event_df['prob_at_least_one_no_cuts'].mean()
        mean_both = event_df['prob_both_decay'].mean()
        
        results['mean_single_particle_decay_prob'].append(mean_single)
        results['mean_single_no_cuts'].append(mean_single_nc)
        results['mean_at_least_one_decay_prob'].append(mean_p1)
        results['mean_at_least_one_no_cuts'].append(mean_p1_nc)
        results['mean_both_decay_prob'].append(mean_both)
        results['exclusion'].append(
            N_SIG_EXCL / (mean_p1 * LUMI_FB * SIGMA_H_FB))
        results['exclusion_no_cuts'].append(
            N_SIG_EXCL / (mean_p1_nc * LUMI_FB * SIGMA_H_FB))
        results['mean_acceptance'].append(mean_acc)
    
    return results


def _first_forward_hit(mesh, origins, dirs):
    """
    Batch ray-cast: for each ray return the nearest forward intersection
    with *mesh*, or NaN if the ray misses. origins/dirs are (M, 3).
    """
    out = np.full((len(origins), 3), np.nan)
    locs, ray_idx, _ = mesh.ray.intersects_location(
        ray_origins=origins, ray_directions=dirs)
    if len(locs) == 0:
        return out
    signed = np.einsum('ij,ij->i', locs - origins[ray_idx], dirs[ray_idx])
    fwd = signed > 1e-6
    locs, ray_idx, signed = locs[fwd], ray_idx[fwd], signed[fwd]
    if len(locs) == 0:
        return out
    order = np.lexsort((signed, ray_idx))
    locs, ray_idx = locs[order], ray_idx[order]
    _, first = np.unique(ray_idx, return_index=True)
    out[ray_idx[first]] = locs[first]
    return out


def sample_separations(geo_cache, lifetime_seconds, n_samples_per_particle=100,
                       rng_seed=42, hit_resolution=HIT_RESOLUTION,
                       n_layers=N_LAYERS, use_3d_reco=True,
                       sigma_t=SIGMA_T_DEFAULT):
    """
    Monte Carlo sample decay positions and rest-frame angles to build
    distributions of electron-pair separations, pointing angles, DCA, and
    the landing surface of each daughter.

    For each particle that hits the fiducial volume:
      1. Sample decay position d *uniformly* within [entry, exit]
      2. Sample |cosθ*| uniformly in [0, 1] and the decay-plane azimuth φ
      3. Compute lab-frame electron angles from boost kinematics
      4. Place true hits at inner/outer tracking layers (separated by
         DETECTOR_THICKNESS), then smear by hit_resolution (local frame)
      5. Reconstruct tracks from smeared hits → separation, pointing, DCA
      6. Check that reconstructed vertex (PCA) is inside the fiducial volume
      7. Build world-frame daughter directions, ray-cast both to the wall,
         and flag whether BOTH land on a tracker surface (scintillator veto)

    Decay positions are sampled uniformly (not exponentially) so the same
    pass can be reweighted to any lifetime by w(d;τ) = (path/N)·(1/λ)e^{-d/λ}.
    The returned ``weights`` are evaluated at ``lifetime_seconds``.

    Set hit_resolution=0 to recover truth-level distributions (DCA=0).

    Returns a dict of aligned 1-D arrays (one entry per sampled decay):
        sep, sep_outer, momenta, pointing, p_soft, dca, vtx_in, open_angle,
        d_implied, collin, on_tracker, weights, event, pid, d, path_len,
        betagamma, plus the scalar n_per (= n_samples_per_particle).
    """
    rng = np.random.default_rng(rng_seed)

    hits = geo_cache['hits']
    hit_idx = np.where(hits)[0]

    all_seps, all_weights, all_momenta, all_pointing = [], [], [], []
    all_p_soft, all_dca, all_vtx_in, all_open = [], [], [], []
    all_sep_out, all_dimp, all_collin = [], [], []
    all_event, all_pid, all_d, all_path, all_bg = [], [], [], [], []
    all_fwd = []
    decay_list, dir1_list, dir2_list = [], [], []

    N = n_samples_per_particle

    for pid, idx in enumerate(hit_idx):
        entry = geo_cache['entry_d'][idx]
        exit_ = geo_cache['exit_d'][idx]
        gamma = geo_cache['gamma'][idx]
        beta = geo_cache['beta'][idx]
        mass = geo_cache['mass'][idx]
        p_llp = geo_cache['momentum'][idx]
        z_hat = geo_cache['direction'][idx]
        event = geo_cache['event'][idx]

        decay_length = calculate_decay_length(p_llp, mass, lifetime_seconds)
        path_length = exit_ - entry

        # Uniform decay-position sampling (lifetime-independent), with the
        # decay-probability density folded into the per-sample weight.
        d_samples = rng.uniform(entry, exit_, N)
        D = exit_ - d_samples  # distance from vertex to inner tracking layer
        w = (path_length / N) * (1.0 / decay_length) \
            * np.exp(-d_samples / decay_length)

        # Sample rest-frame decay angle and decay-plane azimuth
        cos_theta_star = rng.uniform(0, 1, N)
        sin_theta_star = np.sqrt(1 - cos_theta_star**2)
        phi = rng.uniform(0, 2 * np.pi, N)

        # Softer electron momentum (truth, not affected by resolution)
        p_soft = gamma * mass / 2 * (1 - beta * cos_theta_star)

        # Exact lab-frame momentum components (units of M/2; m_e -> 0):
        #   pz1 = gamma (cos0* + beta)        [harder, always forward]
        #   pz2 = gamma (beta - cos0*)        [softer; NEGATIVE if cos0* > beta]
        #   pt  = sin0*                       [same for both]
        # pz2 < 0 means the softer daughter goes BACKWARD in the lab and
        # can never make hits in the forward two-layer tracker. These
        # samples are real decays (their decay-probability weight is
        # kept) but are flagged unreconstructable via `fwd`, mirroring
        # the |cos0*| < beta cap in compute_c_upper for the analytic
        # acceptance.
        pz1 = gamma * (cos_theta_star + beta)
        pz2 = gamma * (beta - cos_theta_star)
        pt = sin_theta_star
        fwd = pz2 > 0

        # Lab-frame angles (small-angle regime), used only by the
        # idealized local-frame fallback below. tan_a2 is clamped so a
        # near-90-degree topology (cos0* -> beta) cannot overflow the
        # batched SVD; the 3D ray-cast path uses the exact dir1/dir2
        # built further down and is unaffected by the clamp.
        tan_a1 = pt / pz1
        with np.errstate(divide='ignore', invalid='ignore'):
            tan_a2 = np.where(np.abs(pz2) > 1e-12, pt / pz2,
                              np.sign(pz2 + 1e-300) * 1e12)
        tan_a2 = np.clip(tan_a2, -1e9, 1e9)

        L = DETECTOR_THICKNESS

        # ---- 3D reconstruction in local frame ----
        # Origin = true vertex, ẑ = LLP direction, x̂ = decay plane
        # True hit positions at inner (z=D) and outer (z=D+L) tracking layers
        x1_in_true  = D * tan_a1
        x2_in_true  = -D * tan_a2
        x1_out_true = (D + L) * tan_a1
        x2_out_true = -(D + L) * tan_a2
        # y = 0 for both (in the decay plane)

        if hit_resolution > 0:
            noise = rng.normal(0, hit_resolution, (8, N))
            x1_in  = x1_in_true  + noise[0]
            y1_in  = noise[1]
            x1_out = x1_out_true + noise[2]
            y1_out = noise[3]
            x2_in  = x2_in_true  + noise[4]
            y2_in  = noise[5]
            x2_out = x2_out_true + noise[6]
            y2_out = noise[7]
        else:
            x1_in  = x1_in_true;  y1_in  = np.zeros(N)
            x1_out = x1_out_true; y1_out = np.zeros(N)
            x2_in  = x2_in_true;  y2_in  = np.zeros(N)
            x2_out = x2_out_true; y2_out = np.zeros(N)

        # Reconstructed track directions (unnormalised)
        dx1 = x1_out - x1_in;  dy1 = y1_out - y1_in   # dz = L
        dx2 = x2_out - x2_in;  dy2 = y2_out - y2_in

        # --- Separation at inner and outer tracking layers ---
        sep     = np.sqrt((x1_in  - x2_in )**2 + (y1_in  - y2_in )**2)
        sep_out = np.sqrt((x1_out - x2_out)**2 + (y1_out - y2_out)**2)

        # --- Pointing angle (bisector vs ẑ = LLP direction) ---
        mag1 = np.sqrt(dx1**2 + dy1**2 + L**2)
        mag2 = np.sqrt(dx2**2 + dy2**2 + L**2)
        bx = dx1 / mag1 + dx2 / mag2
        by = dy1 / mag1 + dy2 / mag2
        bz = L / mag1 + L / mag2
        bmag = np.sqrt(bx**2 + by**2 + bz**2)
        pointing = np.arccos(np.clip(bz / bmag, -1, 1))

        # --- Opening angle between the two reconstructed tracks ---
        cos_open = (dx1 * dx2 + dy1 * dy2 + L * L) / (mag1 * mag2)
        open_angle = np.arccos(np.clip(cos_open, -1, 1))

        # --- Implied vertex distance: d_implied = sep_inner / open_angle ---
        # Geometric identity for two tracks emerging from a common point:
        # sep_inner = open_angle × (vertex → inner-layer distance).
        # Robust against the PCA-degeneracy that hits near-parallel tracks.
        d_implied = np.where(open_angle > 1e-9, sep / open_angle, np.inf)

        # --- Collinearity: RMS perpendicular residual of the 4 hits to their
        # best-fit line. ≈ σ_hit for a single-line topology (IP-muon transit),
        # ≈ (sep_inner + sep_outer)/4 for a true V-shaped decay vertex.
        pts = np.stack([
            np.stack([x1_in,  y1_in,  np.zeros(N)],   axis=-1),
            np.stack([x1_out, y1_out, np.full(N, L)], axis=-1),
            np.stack([x2_in,  y2_in,  np.zeros(N)],   axis=-1),
            np.stack([x2_out, y2_out, np.full(N, L)], axis=-1),
        ], axis=1)  # (N, 4, 3)
        ctr = pts.mean(axis=1, keepdims=True)
        _, S_svd, _ = np.linalg.svd(pts - ctr, full_matrices=False)
        collin = np.sqrt((S_svd[:, 1]**2 + S_svd[:, 2]**2) / 4)

        # --- DCA between reconstructed tracks ---
        # Lines: P_i + t_i * d_i, with P_i at inner-layer hits (z = D)
        # Cross product n = d1 × d2
        nx = dy1 * L - L * dy2
        ny = L * dx2 - dx1 * L
        nz = dx1 * dy2 - dy1 * dx2
        n_mag = np.sqrt(nx**2 + ny**2 + nz**2)

        wx = x1_in - x2_in
        wy = y1_in - y2_in
        # wz = 0 (both at z = D)

        dca = np.where(n_mag > 1e-30,
                       np.abs(wx * nx + wy * ny) / n_mag,
                       np.sqrt(wx**2 + wy**2))

        # --- PCA (point of closest approach) → reconstructed vertex ---
        a_coeff = dx1**2 + dy1**2 + L**2
        b_coeff = dx1 * dx2 + dy1 * dy2 + L**2
        c_coeff = dx2**2 + dy2**2 + L**2
        d_coeff = dx1 * wx + dy1 * wy
        e_coeff = dx2 * wx + dy2 * wy
        det = a_coeff * c_coeff - b_coeff**2
        det_safe = np.where(np.abs(det) < 1e-30, 1e-30, det)
        t1 = (b_coeff * e_coeff - c_coeff * d_coeff) / det_safe
        t2 = (a_coeff * e_coeff - b_coeff * d_coeff) / det_safe

        # PCA midpoint z in local frame (vertex at z=0)
        pca_z = D + 0.5 * (t1 + t2) * L

        # Reconstructed vertex distance from IP along LLP direction
        d_vtx = d_samples + pca_z
        vtx_in = (d_vtx > entry) & (d_vtx < exit_)

        # ---- World-frame daughter directions for surface ray-casting ----
        # Build a transverse basis (x̂,ŷ) ⟂ the LLP direction ẑ and orient
        # the decay plane at azimuth φ. The harder/softer lab angles and the
        # opposite-side convention match the local-frame reconstruction above.
        ref = np.array([0., 1., 0.]) if abs(z_hat[1]) < 0.9 \
            else np.array([1., 0., 0.])
        x_hat = np.cross(z_hat, ref); x_hat /= np.linalg.norm(x_hat)
        y_hat = np.cross(z_hat, x_hat)
        u = np.cos(phi)[:, None] * x_hat + np.sin(phi)[:, None] * y_hat

        # Exact unit vectors from the signed momentum components — keeps
        # the correct (negative) z-component when the softer daughter is
        # backward, and stays finite at cos0* = beta where tan_a2 blows up.
        n1 = np.sqrt(pz1**2 + pt**2)
        n2 = np.sqrt(pz2**2 + pt**2)
        dir1 = (pz1 / n1)[:, None] * z_hat + (pt / n1)[:, None] * u
        dir2 = (pz2 / n2)[:, None] * z_hat - (pt / n2)[:, None] * u
        decay_pos = d_samples[:, None] * z_hat  # IP at the origin

        decay_list.append(decay_pos)
        dir1_list.append(dir1)
        dir2_list.append(dir2)

        all_seps.append(sep)
        all_weights.append(w)
        all_momenta.append(np.full(N, p_llp))
        all_pointing.append(pointing)
        all_p_soft.append(p_soft)
        all_dca.append(dca)
        all_vtx_in.append(vtx_in)
        all_open.append(open_angle)
        all_sep_out.append(sep_out)
        all_dimp.append(d_implied)
        all_collin.append(collin)
        all_event.append(np.full(N, event))
        all_pid.append(np.full(N, pid))
        all_d.append(d_samples)
        all_path.append(np.full(N, path_length))
        all_bg.append(np.full(N, gamma * beta))
        all_fwd.append(fwd)

    if not all_seps:
        empty = np.array([])
        empty3 = np.empty((0, 3))
        return {k: empty for k in (
            'sep', 'sep_outer', 'momenta', 'pointing', 'p_soft', 'dca',
            'vtx_in', 'open_angle', 'd_implied', 'collin', 'on_tracker',
            'xy_disp_1', 'xy_disp_2',
            'weights', 'event', 'pid', 'd', 'path_len', 'betagamma')} \
            | {k: empty3 for k in ('decay_pos', 'dir1', 'dir2',
                                   'exit_pt_1', 'exit_pt_2')} \
            | {'n_per': N}

    # ---- Batch ray-cast both daughters and classify landing surface ----
    decay_all = np.concatenate(decay_list)
    dir1_all = np.concatenate(dir1_list)
    dir2_all = np.concatenate(dir2_list)
    n_tot = len(decay_all)

    origins = np.concatenate([decay_all, decay_all])
    dirs = np.concatenate([dir1_all, dir2_all])
    print(f"  Ray-casting {2 * n_tot} daughter rays for surface "
          f"classification...")
    exit_pts = _first_forward_hit(mesh_fiducial, origins, dirs)

    on_trk = np.zeros(2 * n_tot, dtype=bool)
    valid = ~np.isnan(exit_pts[:, 0])
    if valid.any():
        on_trk[valid] = points_on_tracker(exit_pts[valid])
    # Backward softer daughters (cos0* > beta) cannot make the two
    # forward-layer hits the reconstruction assumes: fold the forward
    # requirement into on_tracker so every selection path (cutflow,
    # selection_mask, MC exclusion) inherits the cap automatically.
    fwd_all = np.concatenate(all_fwd)
    on_tracker = on_trk[:n_tot] & on_trk[n_tot:] & fwd_all

    # ---- Per-daughter on-surface displacement between tracker layers ----
    # In the local cavern basis at the daughter's wall hit p1, the
    # layer1->layer2 displacement decomposes into a radial component
    # (= L, the layer spacing) and two tangential components: one along
    # the centreline direction ("down the tunnel", Delta_s) and one
    # around the cross-section profile ("along the arc", Delta_arc).
    # The experimentally observed lateral separation between the two hits
    # is sqrt(Delta_s^2 + Delta_arc^2). Daughters whose predicted layer-2
    # hit lands off the tracker are marked NaN.
    xy_disp = np.full(2 * n_tot, np.nan)
    if valid.any():
        theta_e, _, tangent_e, right_e, up_e = classify_points_with_basis(
            exit_pts[valid])
        cos_th = np.cos(theta_e)[:, None]
        sin_th = np.sin(theta_e)[:, None]
        n_hat = cos_th * right_e + sin_th * up_e
        t_arc = -sin_th * right_e + cos_th * up_e
        d_valid = dirs[valid]
        d_dot_n   = np.einsum('ij,ij->i', d_valid, n_hat)
        d_dot_tan = np.einsum('ij,ij->i', d_valid, tangent_e)
        d_dot_arc = np.einsum('ij,ij->i', d_valid, t_arc)

        safe_n = np.where(np.abs(d_dot_n) < 1e-6, np.nan, d_dot_n)
        lateral = DETECTOR_THICKNESS \
            * np.sqrt(d_dot_tan**2 + d_dot_arc**2) / np.abs(safe_n)

        # Filter: predicted layer-2 hit must also land on the tracker
        lam = DETECTOR_THICKNESS / safe_n
        p2_pred = exit_pts[valid] + lam[:, None] * d_valid
        finite_p2 = np.all(np.isfinite(p2_pred), axis=1)
        p2_on_trk = np.zeros(len(p2_pred), dtype=bool)
        if finite_p2.any():
            p2_on_trk[finite_p2] = points_on_tracker(p2_pred[finite_p2])
        lateral = np.where(p2_on_trk, lateral, np.nan)
        xy_disp[valid] = lateral
    xy_disp_1 = xy_disp[:n_tot]
    xy_disp_2 = xy_disp[n_tot:]

    # ---- Geometric observables ----
    # By default (use_3d_reco) recompute sep/open/DCA/collin/pointing/vertex from
    # the REAL 3D daughter hits via the shared reco_common routine — the same
    # code the cosmic-decay background uses — so the two are guaranteed
    # consistent. Set use_3d_reco=False to keep the old idealized local-frame
    # values (the previous behaviour, retained for comparison).
    sep_arr   = np.concatenate(all_seps)
    sepo_arr  = np.concatenate(all_sep_out)
    point_arr = np.concatenate(all_pointing)
    dca_arr   = np.concatenate(all_dca)
    vtx_arr   = np.concatenate(all_vtx_in)
    open_arr  = np.concatenate(all_open)
    collin_arr = np.concatenate(all_collin)

    timing_arr = np.full(n_tot, np.inf)   # non-reconstructable -> fail timing
    if use_3d_reco:
        import reco_common as _rc
        in1, out1 = _rc.wall_inner_outer(exit_pts[:n_tot], dir1_all,
                                         DETECTOR_THICKNESS)
        in2, out2 = _rc.wall_inner_outer(exit_pts[n_tot:], dir2_all,
                                         DETECTOR_THICKNESS)
        fin = (np.all(np.isfinite(in1), 1) & np.all(np.isfinite(out1), 1)
               & np.all(np.isfinite(in2), 1) & np.all(np.isfinite(out2), 1))
        if fin.any():
            g = _rc.reconstruct_3d(out1[fin], in1[fin], in2[fin], out2[fin],
                                   hit_resolution, rng)
            sep_arr[fin]    = g['sep']
            sepo_arr[fin]   = g['sep_outer']
            open_arr[fin]   = g['open_angle']
            dca_arr[fin]    = g['dca']
            collin_arr[fin] = g['collin']
            point_arr[fin]  = g['pointing']
            vtx_arr[fin]    = g['vtx_in']
            # Timing chi2, identical model to the cosmic background: both
            # daughters are OUTGOING from the vertex (sign +1) at c (beta ~ 1
            # for GeV electrons). True hits set the true times; the smeared
            # hits + reco vertex set the predicted times.
            # NOTE: this is a single frozen random draw per particle, reused for
            # every lifetime point in the exclusion scan (sample_separations is
            # called once and reweighted). Unbiased, but it adds correlated noise
            # along the lifetime curve near the chi2 cut boundary -- only relevant
            # if few-percent wiggles appear in the exclusion curve.
            nf = int(fin.sum())
            true_hits = np.stack([out1[fin], in1[fin], in2[fin], out2[fin]], 1)
            smeared = np.stack([g['H_out1'], g['H_in1'],
                                g['H_in2'], g['H_out2']], 1)
            tchi2, _, _, _ = _rc.timing_chi2_4hit(
                true_hits, decay_all[fin], np.ones((nf, 4)), np.ones((nf, 4)),
                smeared, g['V_reco'], sigma_t, rng)
            timing_arr[fin] = tchi2
        vtx_arr[~fin] = False     # no well-defined 3D hits -> not reconstructable
        # Fold the 4-hit reconstructability (both daughters give a valid bounded
        # tracker pair) into on_tracker, so grazing / would-exit tracks die at the
        # "both daughters on tracker" stage -- matching the cosmic cutflow, where
        # ok_pair is folded into on_tracker. The net selection is unchanged; this
        # only makes the stage-by-stage efficiencies directly comparable.
        on_tracker = on_tracker & fin
    else:
        # Idealized local-frame fallback (kept only for compare_reco.py): it
        # builds no 3D hits, so the timing chi2 cannot be simulated. Treat timing
        # as a no-op (pass) rather than leaving timing_arr = +inf, which
        # selection_mask (apply_timing=True by default) would otherwise turn into
        # zero acceptance for the whole sample.
        timing_arr[:] = 0.0

    # Backward softer daughter: no valid forward vertex in either the
    # 3D-reco or the idealized local-frame path.
    vtx_arr[~fwd_all] = False

    return {
        'sep': sep_arr,
        'sep_outer': sepo_arr,
        'momenta': np.concatenate(all_momenta),
        'pointing': point_arr,
        'p_soft': np.concatenate(all_p_soft),
        'dca': dca_arr,
        'vtx_in': vtx_arr,
        'open_angle': open_arr,
        'd_implied': np.concatenate(all_dimp),
        'collin': collin_arr,
        'timing_chi2': timing_arr,
        'xy_disp_1': xy_disp_1,
        'xy_disp_2': xy_disp_2,
        'on_tracker': on_tracker,
        'weights': np.concatenate(all_weights),
        'event': np.concatenate(all_event),
        'pid': np.concatenate(all_pid),
        'd': np.concatenate(all_d),
        'path_len': np.concatenate(all_path),
        'betagamma': np.concatenate(all_bg),
        'decay_pos': decay_all,
        'dir1': dir1_all,
        'dir2': dir2_all,
        'exit_pt_1': exit_pts[:n_tot],
        'exit_pt_2': exit_pts[n_tot:],
        'n_per': N,
    }




def build_cutflow(seps, pointing, weights, momenta, p_soft, dca, vtx_in,
                  open_angle, sep_outer, collinearity, on_tracker=None,
                  p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                  dca_cut=DCA_CUT, theta_parallel=THETA_PARALLEL,
                  sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                  collin_min=COLLIN_MIN,
                  sep_out_collin_gate=SEP_OUT_COLLIN_GATE,
                  collin_frac=COLLIN_FRAC,
                  sep_out_gate=SEP_OUT_GATE,
                  point_tight_sep_in=POINT_TIGHT_SEP_IN,
                  sep_in_point_gate=SEP_IN_POINT_GATE,
                  point_global=POINT_GLOBAL,
                  timing_chi2=None, chi2_timing_max=CHI2_TIMING_MAX,
                  pointing_cut=None):
    """
    Apply signal selection cuts sequentially and return a cutflow table.

    Parameters
    ----------
    on_tracker : array of bool or None
        Whether both daughters land on a tracker surface (Arch/Ceiling or
        Left Wall). When provided, a "both daughters on tracker" cut is
        applied right after the fiducial-decay requirement — pairs landing
        on scintillator (Floor / Right Wall) cannot be reconstructed.
    dca_cut : float
        Maximum DCA between reconstructed tracks (m).
    theta_parallel, sep_out_max_parallel : float
        Conditional max-sep_outer cut: when open_angle < theta_parallel,
        require sep_outer < sep_out_max_parallel. For larger open_angle
        no upper bound on sep_outer is applied.
    collin_min, sep_out_collin_gate : float
        Retired static IP-muon-transit veto. Accepted for signature
        compatibility but NOT applied — the active veto is the L-scaled
        version (collin_frac / sep_out_gate): when sep_outer >
        sep_out_gate, require collinearity > collin_frac * L.
    pointing_cut : float or None
        Maximum pointing angle (rad). If None, no pointing cut is applied.

    Returns
    -------
    cutflow : list of dicts with keys 'cut', 'efficiency',
              'marginal_efficiency', 'weighted_yield'
    """
    total_w = weights.sum()
    if total_w == 0:
        return []

    mask = np.ones(len(seps), dtype=bool)
    rows = []

    def add_row(name, m):
        w = weights[m].sum()
        prev_w = weights[prev_mask].sum()
        rows.append({
            'cut': name,
            'weighted_yield': w,
            'efficiency': w / total_w,
            'marginal_efficiency': w / prev_w if prev_w > 0 else 0,
        })

    prev_mask = mask.copy()
    add_row('Decay in fiducial', mask)

    if on_tracker is not None:
        prev_mask = mask.copy()
        mask = mask & on_tracker
        add_row('Both daughters on tracker', mask)

    prev_mask = mask.copy()
    mask = mask & (p_soft >= p_cut)
    add_row(f'p_soft > {p_cut*1000:.0f} MeV/c', mask)

    prev_mask = mask.copy()
    mask = mask & (seps >= sep_min) & (sep_outer >= sep_min)
    add_row(f'sep_in & sep_out > {sep_min*1000:.0f} mm', mask)

    prev_mask = mask.copy()
    mask = mask & (seps <= sep_max)
    add_row(f'sep_in < {sep_max*100:.0f} cm', mask)

    prev_mask = mask.copy()
    mask = mask & (dca <= dca_cut)
    add_row(f'DCA < {dca_cut*100:.1f} cm', mask)

    # Conditional max sep_outer: tight bound only when tracks are
    # nearly parallel (open_angle below the angular-resolution floor)
    prev_mask = mask.copy()
    is_parallel = open_angle < theta_parallel
    parallel_ok = (~is_parallel) | (sep_outer < sep_out_max_parallel)
    mask = mask & parallel_ok
    add_row(f'sep_out<{sep_out_max_parallel:.1f}m if θ<{theta_parallel*1000:.0f}mrad',
            mask)

    # Conditional collinearity cut (L-scaled): when sep_outer > SEP_OUT_GATE
    # require collinearity > COLLIN_FRAC * L. The gate stays closed for
    # collimated signal (small sep_outer); shared with the cosmic-decay background.
    prev_mask = mask.copy()
    _Lc = DETECTOR_THICKNESS
    gated = sep_outer > sep_out_gate
    collin_ok = (~gated) | (collinearity > collin_frac * _Lc)
    mask = mask & collin_ok
    add_row(f'collin>{collin_frac:.2f}L if sep_out>{sep_out_gate*100:.0f}cm '
            f'({collin_frac*_Lc*1000:.0f}mm)', mask)

    prev_mask = mask.copy()
    mask = mask & vtx_in
    add_row('Vertex in fiducial (PCA)', mask)

    # Conditional tight pointing: kills off-IP background for collimated
    # decays without touching wider-sep topologies whose pointing is noisy.
    prev_mask = mask.copy()
    gated_pt = seps < sep_in_point_gate
    point_ok = (~gated_pt) | (pointing < point_tight_sep_in)
    mask = mask & point_ok
    add_row(f'pointing<{point_tight_sep_in*1000:.0f}mrad if sep_in<{sep_in_point_gate*100:.0f}cm',
            mask)

    # Global pointing cut (all sep_in): clips the >1 rad tail.
    prev_mask = mask.copy()
    mask = mask & (pointing < point_global)
    add_row(f'global pointing < {point_global*1000:.0f} mrad', mask)

    # Timing chi2 (same model as the cosmic background): signal daughters are
    # outgoing from the vertex at c, so this is smearing-limited (= chi2(3) CDF,
    # ~0.89 at chi2<6, ~0.97 at chi2<9) and sigma_t-independent, flat in mass.
    if timing_chi2 is not None:
        prev_mask = mask.copy()
        mask = mask & (timing_chi2 < chi2_timing_max)
        add_row(f'timing chi2 < {chi2_timing_max:.0f}', mask)

    if pointing_cut is not None:
        prev_mask = mask.copy()
        mask = mask & (pointing <= pointing_cut)
        add_row(f'pointing < {pointing_cut*1000:.1f} mrad', mask)

    return rows


def selection_mask(mc, p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                   dca_cut=DCA_CUT, theta_parallel=THETA_PARALLEL,
                   sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                   collin_min=COLLIN_MIN,
                   sep_out_collin_gate=SEP_OUT_COLLIN_GATE,
                   collin_frac=COLLIN_FRAC,
                   sep_out_gate=SEP_OUT_GATE,
                   point_tight_sep_in=POINT_TIGHT_SEP_IN,
                   sep_in_point_gate=SEP_IN_POINT_GATE,
                   point_global=POINT_GLOBAL,
                   vtx_inner_min=VTX_INNER_MIN,
                   apply_timing=True, chi2_timing_max=CHI2_TIMING_MAX):
    """
    Boolean per-sample mask for the full signal selection (everything except
    the implicit 'decay in fiducial', which holds for all samples by
    construction). Includes the both-daughters-on-tracker requirement.
    """
    sep        = mc['sep']
    sep_outer  = mc['sep_outer']
    p_soft     = mc['p_soft']
    dca        = mc['dca']
    open_angle = mc['open_angle']
    collin     = mc['collin']
    vtx_in     = mc['vtx_in']
    on_tracker = mc['on_tracker']
    pointing   = mc['pointing']

    m = (p_soft >= p_cut)
    m &= (sep >= sep_min) & (sep_outer >= sep_min)
    m &= (sep <= sep_max)
    m &= (dca <= dca_cut)
    is_parallel = open_angle < theta_parallel
    m &= (~is_parallel) | (sep_outer < sep_out_max_parallel)
    gated = sep_outer > sep_out_gate
    m &= (~gated) | (collin > collin_frac * DETECTOR_THICKNESS)
    gated_pt = sep < sep_in_point_gate
    m &= (~gated_pt) | (pointing < point_tight_sep_in)
    m &= pointing < point_global
    if vtx_inner_min > 0:
        # Reconstructed vertex must stand off from the inner tracking layer.
        m &= (mc['d_implied'] >= vtx_inner_min)
    if apply_timing and 'timing_chi2' in mc:
        m &= mc['timing_chi2'] < chi2_timing_max
    m &= vtx_in
    m &= on_tracker
    return m


def mc_exclusion_vs_lifetime(mc, lifetimes, total_events, **cut_kwargs):
    """
    Build the exclusion curve from a single uniform-sampled MC pass.

    The geometry of each sampled decay (separations, surfaces, vertex, …) is
    lifetime-independent, so the full selection mask is evaluated once. For
    each lifetime the decay-probability weight of every sample is recomputed
    as w(d;τ) = (path/N)·(1/λ)·e^{-d/λ} with λ = βγcτ, and summed per particle
    to give that particle's decay-and-pass probability. Particles are then
    combined per event into P(≥1 decays & passes), averaged over ALL events
    (non-hitting events contribute 0), and converted to an excluded BR.

    Returns a dict mirroring analyze_decay_vs_lifetime: 'lifetimes',
    'mean_at_least_one_decay_prob', 'exclusion', 'mean_single_pass_prob'.
    """
    sel = selection_mask(mc, **cut_kwargs)
    pid = mc['pid']
    d = mc['d']
    path_len = mc['path_len']
    bg = mc['betagamma']
    N = mc['n_per']

    n_part = int(pid.max()) + 1 if len(pid) else 0
    # Event id per particle (constant within a particle's samples)
    pid_event = np.zeros(n_part, dtype=mc['event'].dtype)
    if n_part:
        first = np.unique(pid, return_index=True)[1]
        pid_event[pid[first]] = mc['event'][first]

    mean_p1, exclusion, mean_single = [], [], []

    for tau in lifetimes:
        decay_length = bg * SPEED_OF_LIGHT * tau
        w = (path_len / N) * (1.0 / decay_length) * np.exp(-d / decay_length)
        w_pass = np.where(sel, w, 0.0)

        # Per-particle decay-and-pass probability
        p_part = np.bincount(pid, weights=w_pass, minlength=n_part)

        # Combine particles into per-event P(≥1) via 1 - Π(1 - p_i)
        log_surv = np.bincount(pid_event,
                               weights=np.log1p(-np.clip(p_part, 0, 1 - 1e-15)),
                               minlength=0) if n_part else np.array([])
        p_event = 1.0 - np.exp(log_surv) if len(log_surv) else np.array([])

        mean_at_least_one = p_event.sum() / total_events
        mean_p1.append(mean_at_least_one)
        mean_single.append(p_part.mean() if n_part else 0.0)
        if mean_at_least_one > 0:
            exclusion.append(
                N_SIG_EXCL / (mean_at_least_one * LUMI_FB * SIGMA_H_FB))
        else:
            exclusion.append(np.inf)

    return {
        'lifetimes': lifetimes,
        'mean_at_least_one_decay_prob': np.array(mean_p1),
        'exclusion': np.array(exclusion),
        'mean_single_pass_prob': np.array(mean_single),
    }


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    import sys
    # Optional per-mass override:  python decayProbPerEvent_2body.py <csv> <outString>
    if len(sys.argv) >= 3:
        sample_csv, outString = sys.argv[1], sys.argv[2]
    origin = [0, 0, 0]

    # All outputs go under <outString>/, event displays in
    # <outString>/event_displays/.
    OUT_DIR = outString
    EV_DIR = os.path.join(OUT_DIR, 'event_displays')
    os.makedirs(EV_DIR, exist_ok=True)

    print(f"  Output directory: {OUT_DIR}/")
    print(f"  Detector thickness: {DETECTOR_THICKNESS*100:.0f} cm")
    print(f"  Cuts: p_e > {P_CUT*1000:.0f} MeV/c, "
          f"{SEP_MIN*1000:.0f} mm < separation < {SEP_MAX*100:.0f} cm")
    
    geo_cache = cache_geometry(sample_csv, mesh_fiducial, origin)
    
    # --- Diagnostics ---
    print("\n" + "="*50)
    print("ACCEPTANCE DIAGNOSTICS")
    print("="*50)
    
    test_mass = 15.0
    print(f"\nM = {test_mass} GeV → e+e-, p_cut = {P_CUT*1000:.0f} MeV/c, "
          f"sep: {SEP_MIN*1000:.0f} mm – {SEP_MAX*100:.0f} cm")
    print(f"\n{'p_LLP':>8} {'γ':>6} {'c_P':>7} {'θ_min':>10} "
          f"{'d(1mm)':>8} {'d(10cm)':>9} {'θ_12=10cm/1m':>14}")
    
    for p_test in [20, 50, 100, 200, 500]:
        E = np.sqrt(p_test**2 + test_mass**2)
        g = E / test_mass
        b = p_test / E
        c_P = compute_c_upper(g, b, test_mass)
        theta_min = 2.0 / g
        d_min_sep = SEP_MIN / theta_min
        d_max_sep = SEP_MAX / theta_min
        sep_at_1m = theta_min * 1.0
        print(f"{p_test:>8.0f} {g:>6.1f} {c_P:>7.4f} {theta_min*1000:>8.1f} mrad "
              f"{d_min_sep*100:>7.2f}cm {d_max_sep:>7.1f}m {sep_at_1m*100:>12.1f}cm")
    
    print(f"\nMax sep cut effect: at d_remaining = 1m, θ_max = {SEP_MAX}m / 1m = 100 mrad")
    print("  → kills decays far from detector for low-γ (large opening angle)")
    print(f"  → for γ=1.67 (p=20 GeV): θ_min=1200 mrad >> 100 mrad → heavily constrained")
    print(f"  → for γ=33.4 (p=500 GeV): θ_min=60 mrad < 100 mrad → mostly unconstrained")
    
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
    print(f"  Mean P(≥1) with cuts:    {event_df['prob_at_least_one_decays'].mean():.6e}")
    print(f"  Mean P(≥1) without cuts: {event_df['prob_at_least_one_no_cuts'].mean():.6e}")
    
    # --- Separation histogram ---
    print("\n" + "="*50)
    print("SEPARATION DISTRIBUTION")
    print("="*50)
    
    mc = sample_separations(geo_cache, lifetime, n_samples_per_particle=200)
    seps         = mc['sep']
    weights      = mc['weights']
    momenta      = mc['momenta']
    pointing     = mc['pointing']
    p_soft       = mc['p_soft']
    dca          = mc['dca']
    vtx_in       = mc['vtx_in']
    open_angle   = mc['open_angle']
    sep_outer    = mc['sep_outer']
    d_implied    = mc['d_implied']
    collinearity = mc['collin']
    timing_chi2  = mc['timing_chi2']
    xy_disp_1    = mc['xy_disp_1']
    xy_disp_2    = mc['xy_disp_2']
    on_tracker   = mc['on_tracker']

    print(f"  Hit resolution: {HIT_RESOLUTION*1000:.1f} mm, "
          f"N layers: {N_LAYERS}, DCA cut: {DCA_CUT*100:.1f} cm")
    print(f"  Conditional max sep_outer: < {SEP_OUT_MAX_PARALLEL:.1f} m "
          f"when open_angle < {THETA_PARALLEL*1000:.0f} mrad")
    print(f"  Collinearity veto: collinearity > {COLLIN_FRAC*DETECTOR_THICKNESS*1000:.0f} mm "
          f"({COLLIN_FRAC:.2f}L) when sep_outer > {SEP_OUT_GATE*100:.0f} cm")
    print(f"  Global pointing cut: < {POINT_GLOBAL*1000:.0f} mrad "
          f"(tight {POINT_TIGHT_SEP_IN*1000:.0f} mrad when sep_in < {SEP_IN_POINT_GATE*100:.0f} cm)")
    if len(seps) > 0:
        w_on = weights[on_tracker].sum() / weights.sum() if weights.sum() else 0
        print(f"  Scintillator veto: both daughters on a tracker surface "
              f"(Arch/Ceiling + Left Wall) = {w_on*100:.1f}% of decay weight")

    if len(seps) > 0:
        fig_sep, axes_sep = plt.subplots(1, 3, figsize=(18, 5))
        
        ax = axes_sep[0]
        bins = np.linspace(0, 0.5, 100)
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
        ax.set_xlim(0, 0.5)
        
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
        plt.savefig(os.path.join(OUT_DIR, 'separation_histogram'+outString+'.png'), dpi=150)
        show_or_close()
        
        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        frac_accepted = weights[in_window].sum() / weights.sum()
        print(f"  Fraction of decays in separation window: {frac_accepted:.3f}")
        print(f"  Median separation (all decays): {np.median(seps)*100:.1f} cm")
        accepted_seps = seps[in_window]
        if len(accepted_seps) > 0:
            print(f"  Median separation (accepted):  {np.median(accepted_seps)*100:.1f} cm")

    # --- Pointing angle distribution ---
    if len(pointing) > 0:
        print("\n" + "="*50)
        print("POINTING ANGLE DISTRIBUTION")
        print("="*50)

        # Apply separation acceptance to pointing angles
        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        pointing_acc = pointing[in_window]
        weights_acc = weights[in_window]
        momenta_acc = momenta[in_window]

        pointing_mrad = pointing_acc * 1000

        fig_pt, axes_pt = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes_pt[0]
        bins_pt = np.linspace(0, np.percentile(pointing_mrad, 99.5), 80)
        ax.hist(pointing_mrad, bins=bins_pt, weights=weights_acc,
                color='darkorange', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.set_xlabel('Pointing angle (mrad)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Pointing angle: bisector vs LLP direction\n'
                     f'(τ = {lifetime*1e9:.0f} ns, accepted decays)')
        median_pt = np.median(pointing_mrad)
        ax.axvline(median_pt, color='red', linestyle='--', linewidth=2,
                   label=f'Median = {median_pt:.2f} mrad')
        ax.axvline(POINT_GLOBAL*1000, color='purple', linestyle='-', linewidth=2,
                   label=f'global cut = {POINT_GLOBAL*1000:.0f} mrad')
        ax.axvline(POINT_TIGHT_SEP_IN*1000, color='blue', linestyle=':', linewidth=2,
                   label=f'tight cut (sep_in<{SEP_IN_POINT_GATE*100:.0f}cm) = {POINT_TIGHT_SEP_IN*1000:.0f} mrad')
        ax.legend(fontsize=8)

        ax2 = axes_pt[1]
        bins_log_pt = np.logspace(np.log10(max(pointing_mrad.min(), 1e-3)),
                                  np.log10(pointing_mrad.max()), 80)
        ax2.hist(pointing_mrad, bins=bins_log_pt, weights=weights_acc,
                 color='darkorange', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.axvline(POINT_GLOBAL*1000, color='purple', linestyle='-', linewidth=2,
                    label=f'global cut = {POINT_GLOBAL*1000:.0f} mrad')
        ax2.axvline(POINT_TIGHT_SEP_IN*1000, color='blue', linestyle=':', linewidth=2,
                    label=f'tight cut = {POINT_TIGHT_SEP_IN*1000:.0f} mrad')
        ax2.set_xscale('log')
        ax2.set_xlabel('Pointing angle (mrad)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('Log-scale pointing angle')
        ax2.legend(fontsize=8)

        ax3 = axes_pt[2]
        mask_fin = np.isfinite(pointing_mrad) & (pointing_mrad > 0)
        h = ax3.hist2d(momenta_acc[mask_fin], pointing_mrad[mask_fin],
                       bins=[np.linspace(0, 500, 50),
                             np.linspace(0, np.percentile(pointing_mrad, 99), 50)],
                       weights=weights_acc[mask_fin],
                       cmap='viridis', cmin=1e-20)
        ax3.set_xlabel('LLP momentum (GeV/c)')
        ax3.set_ylabel('Pointing angle (mrad)')
        ax3.set_title('Pointing angle vs LLP momentum')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'pointing_angle'+outString+'.png'), dpi=150)
        show_or_close()

        print(f"  Median pointing angle (accepted): {median_pt:.3f} mrad")
        p90 = np.percentile(pointing_mrad, 90)
        p99 = np.percentile(pointing_mrad, 99)
        print(f"  90th percentile: {p90:.3f} mrad")
        print(f"  99th percentile: {p99:.3f} mrad")

        # Convert to spatial resolution at the IP distance
        mean_entry = np.nanmean(geo_cache['entry_d'][geo_cache['hits']])
        print(f"  At mean entry distance ({mean_entry:.0f} m):")
        print(f"    Median miss distance: {median_pt * mean_entry:.1f} mm")
        print(f"    90th pct miss distance: {p90 * mean_entry:.1f} mm")

    # --- DCA distribution ---
    if len(dca) > 0:
        print("\n" + "="*50)
        print("DCA DISTRIBUTION")
        print("="*50)

        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        dca_acc = dca[in_window]
        weights_dca = weights[in_window]
        momenta_dca = momenta[in_window]

        dca_cm = dca_acc * 100

        fig_dca, axes_dca = plt.subplots(1, 3, figsize=(18, 5))

        ax = axes_dca[0]
        pct99 = np.percentile(dca_cm, 99.5) if len(dca_cm) > 0 else 10
        bins_dca = np.linspace(0, max(pct99, DCA_CUT * 100 * 2), 80)
        ax.hist(dca_cm, bins=bins_dca, weights=weights_dca,
                color='seagreen', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax.axvline(DCA_CUT * 100, color='red', linestyle='--', linewidth=2,
                   label=f'DCA cut = {DCA_CUT*100:.1f} cm')
        ax.set_xlabel('DCA (cm)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Track DCA (τ = {lifetime*1e9:.0f} ns, '
                     f'σ_hit = {HIT_RESOLUTION*1000:.1f} mm)')
        ax.legend(fontsize=9)

        ax2 = axes_dca[1]
        bins_log_dca = np.logspace(np.log10(max(dca_cm[dca_cm > 0].min(), 1e-4))
                                   if np.any(dca_cm > 0) else -4,
                                   np.log10(max(dca_cm.max(), 1)), 80)
        ax2.hist(dca_cm, bins=bins_log_dca, weights=weights_dca,
                 color='seagreen', edgecolor='black', linewidth=0.3, alpha=0.8)
        ax2.axvline(DCA_CUT * 100, color='red', linestyle='--', linewidth=2,
                    label=f'DCA cut = {DCA_CUT*100:.1f} cm')
        ax2.set_xscale('log')
        ax2.set_xlabel('DCA (cm)')
        ax2.set_ylabel('Weighted counts (decay prob.)')
        ax2.set_title('DCA (log scale)')
        ax2.legend(fontsize=9)

        ax3 = axes_dca[2]
        mask_fin = np.isfinite(dca_cm) & (dca_cm > 0)
        if np.any(mask_fin):
            h = ax3.hist2d(momenta_dca[mask_fin], dca_cm[mask_fin],
                           bins=[np.linspace(0, 500, 50),
                                 np.linspace(0, max(pct99, DCA_CUT * 100 * 2), 50)],
                           weights=weights_dca[mask_fin],
                           cmap='viridis', cmin=1e-20)
            ax3.axhline(DCA_CUT * 100, color='red', linestyle='--', linewidth=2,
                        label=f'DCA cut = {DCA_CUT*100:.1f} cm')
            ax3.set_xlabel('LLP momentum (GeV/c)')
            ax3.set_ylabel('DCA (cm)')
            ax3.set_title('DCA vs LLP momentum')
            ax3.legend(fontsize=9, loc='upper right')
            plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'dca_' + outString + '.png'), dpi=150)
        show_or_close()

        median_dca = np.median(dca_cm)
        w_total = weights_dca.sum()
        frac_dca = weights_dca[dca_acc <= DCA_CUT].sum() / w_total \
            if w_total > 0 else 0
        vtx_in_acc = vtx_in[in_window]
        frac_vtx = weights_dca[vtx_in_acc].sum() / w_total \
            if w_total > 0 else 0
        print(f"  Median DCA: {median_dca:.4f} cm")
        print(f"  Fraction passing DCA < {DCA_CUT*100:.1f} cm: {frac_dca:.4f}")
        print(f"  Fraction with vertex in fiducial: {frac_vtx:.4f}")
        print(f"  (fractions above computed after separation cuts)")

    # --- Per-track on-surface displacement between tracker layer 1 and 2 ---
    # sqrt(Delta_s^2 + Delta_arc^2): the experimentally observed lateral
    # separation between layer-1 and layer-2 hits on the tracker surface
    # (down-tunnel + around-the-arc components). Radial L is excluded.
    # Daughters whose predicted layer-2 hit falls off the tracker are
    # dropped (NaN).
    if len(xy_disp_1) > 0:
        print("\n" + "="*50)
        print("PER-TRACK LAYER-TO-LAYER XY DISPLACEMENT (on tracker surface)")
        print("="*50)

        sel_mask    = selection_mask(mc)
        xy_pool     = np.concatenate([xy_disp_1[sel_mask],
                                      xy_disp_2[sel_mask]]) * 100   # cm
        weight_pool = np.concatenate([weights[sel_mask],
                                      weights[sel_mask]])

        finite = np.isfinite(xy_pool)
        xy_pool, weight_pool = xy_pool[finite], weight_pool[finite]

        if len(xy_pool) > 0 and weight_pool.sum() > 0:
            fig_tg, axes_tg = plt.subplots(1, 2, figsize=(12, 5))

            pct99 = np.percentile(xy_pool, 99.5)
            pct_axis = np.percentile(xy_pool, 90.0)
            bins_lin = np.linspace(0, max(pct_axis, 1.0), 80)
            ax = axes_tg[0]
            ax.hist(xy_pool, bins=bins_lin, weights=weight_pool,
                    color='darkorange', edgecolor='black', linewidth=0.3,
                    alpha=0.85)
            ax.set_xlabel('xy displacement (cm)  [√(Δs² + Δarc²)]')
            ax.set_ylabel('Weighted counts (decay prob.)')
            ax.set_title(f'Per-track layer1→layer2 on-surface displacement '
                         f'(after full selection, linear axis to 90%)\n'
                         f'(L = {DETECTOR_THICKNESS*100:.0f} cm radial, '
                         f'τ = {lifetime*1e9:.0f} ns)')

            ax2 = axes_tg[1]
            t_pos = xy_pool[xy_pool > 0]
            if len(t_pos) > 0:
                bins_log = np.logspace(np.log10(max(t_pos.min(), 1e-4)),
                                       np.log10(max(t_pos.max(), 1.0)), 80)
                ax2.hist(xy_pool, bins=bins_log, weights=weight_pool,
                         color='darkorange', edgecolor='black', linewidth=0.3,
                         alpha=0.85)
                ax2.set_xscale('log')
            ax2.set_xlabel('xy displacement (cm)')
            ax2.set_ylabel('Weighted counts (decay prob.)')
            ax2.set_title('Same, log scale')

            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, 'xy_disp_' + outString + '.png'), dpi=150)
            show_or_close()

            w_tot = weight_pool.sum()
            mean_xy = np.average(xy_pool, weights=weight_pool)
            sort_idx = np.argsort(xy_pool)
            cw = np.cumsum(weight_pool[sort_idx]) / w_tot
            median_xy = xy_pool[sort_idx][np.searchsorted(cw, 0.5)]
            print(f"  Pool size: {len(xy_pool)} tracks "
                  f"(both daughters of pairs passing full selection)")
            print(f"  Weighted mean:    {mean_xy:.3f} cm")
            print(f"  Weighted median:  {median_xy:.3f} cm")
            print(f"  99.5th percentile: {pct99:.3f} cm")
            for thr in (5, 10, 20, 50):
                frac = weight_pool[xy_pool < thr].sum() / w_tot
                print(f"  Fraction with xy disp. < {thr} cm: {frac:.4f}")

    # --- Pointing angle split by inner separation (low vs high sep_inner) ---
    # Compare pointing for tight (sep_inner < SPLIT) and wide (sep_inner >=
    # SPLIT) decay topologies, with a common kinematic baseline (on_tracker
    # + momentum + outer-sep window). Surfaces whether pointing discriminates
    # the two topologies independently of the inner-separation cut.
    SEP_INNER_SPLIT = 0.10  # m  — split point for the comparison
    if len(pointing) > 0:
        print("\n" + "="*50)
        print(f"POINTING ANGLE: sep_inner < {SEP_INNER_SPLIT*100:.0f} cm "
              f"vs ≥ {SEP_INNER_SPLIT*100:.0f} cm")
        print("="*50)

        m_base = on_tracker & (p_soft >= P_CUT) \
            & (sep_outer >= SEP_MIN) & (seps <= SEP_MAX)
        m_low  = m_base & (seps <  SEP_INNER_SPLIT)
        m_high = m_base & (seps >= SEP_INNER_SPLIT)

        pt_low_mrad  = pointing[m_low]  * 1000
        w_low        = weights[m_low]
        pt_high_mrad = pointing[m_high] * 1000
        w_high       = weights[m_high]

        if (len(pt_low_mrad) > 0 and w_low.sum() > 0
                and len(pt_high_mrad) > 0 and w_high.sum() > 0):
            fig_pf, axes_pf = plt.subplots(1, 2, figsize=(12, 5))

            all_pt = np.concatenate([pt_low_mrad, pt_high_mrad])
            all_w  = np.concatenate([w_low, w_high])
            pct99  = np.percentile(all_pt, 99.5)
            bins_lin = np.linspace(0, max(pct99, 50.0), 80)

            label_low  = (f'sep_in < {SEP_INNER_SPLIT*100:.0f} cm  '
                          f'(yield {w_low.sum():.2e})')
            label_high = (f'sep_in ≥ {SEP_INNER_SPLIT*100:.0f} cm  '
                          f'(yield {w_high.sum():.2e})')

            ax = axes_pf[0]
            ax.hist(pt_low_mrad, bins=bins_lin, weights=w_low,
                    color='crimson', edgecolor='black', linewidth=0.3,
                    alpha=0.75, density=True, label=label_low)
            ax.hist(pt_high_mrad, bins=bins_lin, weights=w_high,
                    color='steelblue', edgecolor='black', linewidth=0.3,
                    alpha=0.5, density=True, label=label_high)
            ax.set_xlabel('Pointing angle (mrad)')
            ax.set_ylabel('Density')
            ax.set_title(f'Pointing angle by inner separation '
                         f'(τ = {lifetime*1e9:.0f} ns)\n'
                         f'(on_trk + p_soft + sep_out baseline)')
            ax.legend(fontsize=9)

            ax2 = axes_pf[1]
            all_pos = all_pt[all_pt > 0]
            if len(all_pos) > 0:
                lo = max(all_pos.min(), 1e-3)
                hi = max(all_pos.max(), 1.0)
                bins_log = np.logspace(np.log10(lo), np.log10(hi), 80)
                ax2.hist(pt_low_mrad, bins=bins_log, weights=w_low,
                         color='crimson', edgecolor='black', linewidth=0.3,
                         alpha=0.75, density=True,
                         label=f'sep_in < {SEP_INNER_SPLIT*100:.0f} cm')
                ax2.hist(pt_high_mrad, bins=bins_log, weights=w_high,
                         color='steelblue', edgecolor='black',
                         linewidth=0.3, alpha=0.5, density=True,
                         label=f'sep_in ≥ {SEP_INNER_SPLIT*100:.0f} cm')
                ax2.set_xscale('log')
            ax2.set_xlabel('Pointing angle (mrad)')
            ax2.set_ylabel('Density')
            ax2.set_title('Same, log scale')
            ax2.legend(fontsize=9)

            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, 'pointing_by_inner_sep_' + outString + '.png'),
                        dpi=150)
            show_or_close()

            def _wmed(x, w):
                wt = w.sum()
                if wt <= 0:
                    return float('nan')
                idx = np.argsort(x)
                return x[idx][np.searchsorted(np.cumsum(w[idx]) / wt, 0.5)]

            for label, pt, w in (("sep_in < split", pt_low_mrad, w_low),
                                  ("sep_in ≥ split", pt_high_mrad, w_high)):
                print(f"  {label}: samples={len(pt)}, "
                      f"yield={w.sum():.4e}, "
                      f"median={_wmed(pt, w):.3f} mrad, "
                      f"mean={np.average(pt, weights=w):.3f} mrad")

    # --- sep_outer distribution, split by open-angle band ---
    # Motivates the conditional sep_outer<X-when-parallel cut:
    # nearly-parallel signal stays at small sep_outer; clearly diverging
    # signal can occupy arbitrary sep_outer (legitimately).
    if len(sep_outer) > 0:
        print("\n" + "="*50)
        print("SEP_OUTER BY OPEN-ANGLE BAND")
        print("="*50)

        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        so   = sep_outer[in_window]
        op   = open_angle[in_window]
        ww   = weights[in_window]

        bands = [
            (f'open_angle < {THETA_PARALLEL*1000:.0f} mrad',
                 (None, THETA_PARALLEL), 'steelblue'),
            (f'{THETA_PARALLEL*1000:.0f} ≤ open_angle < 100 mrad',
                 (THETA_PARALLEL, 0.100), 'darkorange'),
            ('open_angle ≥ 100 mrad',
                 (0.100, None), 'crimson'),
        ]

        fig_so, axes_so = plt.subplots(1, 2, figsize=(14, 5))

        bins_lin = np.linspace(0, max(2.0, SEP_OUT_MAX_PARALLEL * 1.5), 80)
        positive = so[so > 0]
        lo = np.log10(positive.min()) if len(positive) else -3
        hi = np.log10(max(so.max(), 10.0))
        bins_log = np.logspace(lo, hi, 80)

        for ax, bins, scale in [(axes_so[0], bins_lin, 'linear'),
                                (axes_so[1], bins_log, 'log')]:
            for label, (lo_b, hi_b), color in bands:
                mask = np.ones(len(op), dtype=bool)
                if lo_b is not None: mask &= (op >= lo_b)
                if hi_b is not None: mask &= (op <  hi_b)
                if mask.sum() == 0:
                    continue
                w_sum = ww[mask].sum()
                ax.hist(so[mask], bins=bins, weights=ww[mask],
                        histtype='step', linewidth=2, color=color,
                        label=f'{label}  (w={w_sum:.2e})')
            ax.axvline(SEP_OUT_MAX_PARALLEL, color='gray', linestyle='--',
                       linewidth=1.5,
                       label=f'parallel max = {SEP_OUT_MAX_PARALLEL:.1f} m')
            ax.set_xlabel('sep_outer (m)')
            ax.set_ylabel('Weighted counts (decay prob.)')
            ax.set_xscale(scale)
            if scale == 'log':
                ax.set_yscale('log')
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(True, which='both', alpha=0.3)

        axes_so[0].set_title(f'sep_outer by open-angle band  (τ = {lifetime*1e9:.0f} ns)')
        axes_so[1].set_title('Log-log view')
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'sep_outer_bands_' + outString + '.png'), dpi=150)
        show_or_close()

        # Numerical summary
        def wpct(x, w, q):
            if len(x) == 0: return 0.0
            idx = np.argsort(x); xs, ws = x[idx], w[idx]
            c = np.cumsum(ws); return xs[np.searchsorted(c, q*c[-1])]
        print(f"  {'band':<32} {'med':>10} {'90%':>9} {'frac>parallel_max':>20}")
        for label, (lo_b, hi_b), _ in bands:
            mask = np.ones(len(op), dtype=bool)
            if lo_b is not None: mask &= (op >= lo_b)
            if hi_b is not None: mask &= (op <  hi_b)
            if mask.sum() == 0:
                continue
            x_, w_ = so[mask], ww[mask]
            med = wpct(x_, w_, 0.50); p90 = wpct(x_, w_, 0.90)
            frac = w_[x_ > SEP_OUT_MAX_PARALLEL].sum() / w_.sum()
            print(f"  {label:<32} {med*100:>8.2f}cm {p90*100:>7.2f}cm "
                  f"{frac:>20.4f}")

    # --- Collinearity distribution (IP-muon-transit veto diagnostic) ---
    if len(collinearity) > 0:
        print("\n" + "="*50)
        print("COLLINEARITY (IP-muon-transit veto)")
        print("="*50)

        in_window = (seps >= SEP_MIN) & (seps <= SEP_MAX)
        co = collinearity[in_window]
        so = sep_outer[in_window]
        ww = weights[in_window]

        # Active collinearity veto thresholds (single-sourced):
        gate_cm = SEP_OUT_GATE * 100
        collin_cut_mm = COLLIN_FRAC * DETECTOR_THICKNESS * 1000
        below = so <= SEP_OUT_GATE
        above = so >  SEP_OUT_GATE

        fig_co, axes_co = plt.subplots(1, 3, figsize=(18, 5))

        co_mm = co * 1000
        ax = axes_co[0]
        bins = np.linspace(0, max(np.percentile(co_mm, 99.5), collin_cut_mm*1.3), 80)
        for lbl, m, color in [
            (f'sep_out ≤ {gate_cm:.0f} cm (gate closed)',
                 below, 'steelblue'),
            (f'sep_out > {gate_cm:.0f} cm (gate open)',
                 above, 'crimson'),
        ]:
            if m.sum() == 0:
                continue
            ax.hist(co_mm[m], bins=bins, weights=ww[m],
                    histtype='step', linewidth=2, color=color,
                    label=f'{lbl}  (w={ww[m].sum():.2e})')
        ax.axvline(collin_cut_mm, color='gray', linestyle='--', linewidth=1.5,
                   label=f'cut = {COLLIN_FRAC:.2f}L = {collin_cut_mm:.0f} mm')
        ax.axvline(DETECTOR_THICKNESS/2*1000, color='green', linestyle=':',
                   linewidth=1.5, label=f'signal ceiling L/2 = {DETECTOR_THICKNESS/2*1000:.0f} mm')
        ax.set_xlabel('collinearity (mm)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Collinearity by sep_outer gate (τ = {lifetime*1e9:.0f} ns)')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)

        ax2 = axes_co[1]
        positive = co_mm[co_mm > 0]
        lo = np.log10(max(positive.min(), 0.1)) if len(positive) else -1
        hi = np.log10(max(co_mm.max(), 100.0))
        bins_log = np.logspace(lo, hi, 80)
        for lbl, m, color in [
            ('gate closed', below, 'steelblue'),
            ('gate open',   above, 'crimson'),
        ]:
            if m.sum() == 0:
                continue
            ax2.hist(co_mm[m], bins=bins_log, weights=ww[m],
                     histtype='step', linewidth=2, color=color, label=lbl)
        ax2.axvline(collin_cut_mm, color='gray', linestyle='--', linewidth=1.5,
                    label=f'cut = {collin_cut_mm:.0f} mm')
        ax2.axvline(DETECTOR_THICKNESS/2*1000, color='green', linestyle=':',
                    linewidth=1.5, label=f'L/2 = {DETECTOR_THICKNESS/2*1000:.0f} mm')
        ax2.set_xscale('log'); ax2.set_yscale('log')
        ax2.set_xlabel('collinearity (mm)')
        ax2.set_ylabel('Weighted counts')
        ax2.set_title('Log-log view')
        ax2.legend(fontsize=8, loc='upper left')
        ax2.grid(True, which='both', alpha=0.3)

        ax3 = axes_co[2]
        # 2D heat map: collinearity vs sep_outer with rejection region shaded
        mask_fin = np.isfinite(co) & np.isfinite(so) & (co > 0) & (so > 0)
        h = ax3.hist2d(so[mask_fin]*100, co_mm[mask_fin],
                       bins=[np.logspace(0, 3, 50), np.logspace(-1, 3, 50)],
                       weights=ww[mask_fin], cmap='viridis', cmin=1e-30)
        ax3.axvline(gate_cm, color='red', linestyle='--',
                    linewidth=1.5, label=f'gate = {gate_cm:.0f} cm')
        ax3.axhline(collin_cut_mm, color='red', linestyle='--', linewidth=1.5,
                    label=f'cut = {collin_cut_mm:.0f} mm')
        # Shade rejection region (sep_out > gate AND collin < cut)
        ax3.fill_between([gate_cm, 1e3], 1e-1, collin_cut_mm,
                         color='red', alpha=0.15, label='reject (muon-like)')
        ax3.set_xscale('log'); ax3.set_yscale('log')
        ax3.set_xlabel('sep_outer (cm)')
        ax3.set_ylabel('collinearity (mm)')
        ax3.set_title('Conditional cut region')
        ax3.legend(fontsize=8, loc='lower right')
        plt.colorbar(h[3], ax=ax3, label='Weighted counts')

        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'collinearity_' + outString + '.png'), dpi=150)
        show_or_close()

        # Numerical summary
        gated = above
        if gated.sum() > 0:
            x_, w_ = co[gated], ww[gated]
            def wpct(x, w, q):
                idx = np.argsort(x); xs, ws = x[idx], w[idx]
                c = np.cumsum(ws); return xs[np.searchsorted(c, q*c[-1])]
            med = wpct(x_, w_, 0.50); p1 = wpct(x_, w_, 0.01)
            _ccut = COLLIN_FRAC * DETECTOR_THICKNESS
            frac_fail = w_[x_ < _ccut].sum() / w_.sum()
            print(f"  Above gate (sep_out > {SEP_OUT_GATE*100:.0f} cm): "
                  f"median = {med*1000:.1f} mm, 1% = {p1*1000:.1f} mm")
            print(f"  Signal failing collinearity > {_ccut*1000:.0f} mm "
                  f"in this band: {frac_fail*100:.3f}%")
        else:
            print(f"  No events above the gate (sep_out > "
                  f"{SEP_OUT_GATE*100:.0f} cm) for this sample.")

    # --- Reconstructed vertex stand-off distribution (VTX_INNER_MIN cut) ---
    # Reconstructed along-LLP distance from the decay vertex to the inner
    # tracking layer, two estimators compared:
    #   small-angle:  d = sep_in / open_angle          (the current cut variable)
    #   exact:        d = sep_in / (2 tan(open/2))      (wide-angle-correct; the
    #                 bisector bisects open_angle, so each track sits at open/2)
    # 2 tan(x/2) -> x as x -> 0, so they agree for collimated decays and the
    # small-angle form overestimates the stand-off as the opening angle widens.
    # Shown for signal passing every OTHER cut, isolating the marginal effect.
    if len(d_implied) > 0:
        print("\n" + "="*50)
        print("RECO VERTEX STAND-OFF: small-angle vs exact estimator")
        print("="*50)

        pass_other = selection_mask(mc, vtx_inner_min=0.0)
        sep_p = seps[pass_other]
        open_p = open_angle[pass_other]
        wv = weights[pass_other]
        with np.errstate(divide='ignore', invalid='ignore'):
            d_small = np.where(open_p > 1e-9, sep_p / open_p, np.inf)
            d_exact = np.where(open_p > 1e-9,
                               sep_p / (2.0 * np.tan(open_p / 2.0)), np.inf)
        marks = [0.30, 0.60, 1.00]
        tot = wv.sum()

        ests = [('small-angle  $s/\\theta$', d_small, 'steelblue'),
                ('exact  $s/2\\tan(\\theta/2)$', d_exact, 'darkorange')]

        fig_vx, axes_vx = plt.subplots(1, 2, figsize=(12, 5))

        # Left: decay-weighted distributions of both estimators.
        ax = axes_vx[0]
        xmax = 2.0
        bins = np.linspace(0, xmax * 100, 60)
        for lbl, dv, color in ests:
            fin = np.isfinite(dv)
            ax.hist(np.clip(dv[fin] * 100, 0, xmax * 100), bins=bins,
                    weights=wv[fin], histtype='step', linewidth=2, color=color,
                    label=lbl)
        for vm in marks:
            ax.axvline(vm * 100, color='crimson', linestyle='--', linewidth=1.2)
        if VTX_INNER_MIN > 0:
            ax.axvline(VTX_INNER_MIN * 100, color='black', linewidth=2,
                       label=f'active cut = {VTX_INNER_MIN*100:.0f} cm')
        ax.set_xlabel('reco vertex stand-off along LLP (cm)')
        ax.set_ylabel('Weighted counts (decay prob.)')
        ax.set_title(f'Vertex stand-off (τ = {lifetime*1e9:.0f} ns), passing all other cuts')
        ax.set_xlim(0, xmax * 100)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        # Right: signal efficiency retained vs the stand-off requirement.
        ax2 = axes_vx[1]
        thr = np.linspace(0, xmax, 200)
        for lbl, dv, color in ests:
            eff = np.array([wv[dv >= t].sum() / tot for t in thr]) if tot > 0 \
                else np.zeros_like(thr)
            ax2.plot(thr * 100, eff, color=color, linewidth=2, label=lbl)
            for vm in marks:
                e = wv[dv >= vm].sum() / tot if tot > 0 else 0
                ax2.plot([vm * 100], [e], 'o', color=color, ms=5)
        for vm in marks:
            ax2.axvline(vm * 100, color='crimson', linestyle='--', linewidth=1)
        ax2.set_xlabel('reco vertex stand-off requirement (cm)')
        ax2.set_ylabel('signal efficiency retained')
        ax2.set_title('Marginal impact: small-angle vs exact')
        ax2.set_xlim(0, xmax * 100)
        ax2.set_ylim(0, 1.02)
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'vtx_standoff_' + outString + '.png'),
                    dpi=150)
        show_or_close()

        if tot > 0:
            for lbl, dv, _ in ests:
                lbl_txt = 'small-angle' if 'small' in lbl else 'exact      '
                print(f"  {lbl_txt} eff vs stand-off: "
                      + ", ".join(f"{v*100:.0f}cm={wv[dv>=v].sum()/tot:.3f}"
                                  for v in marks))

    # --- Cutflow table ---
    if len(seps) > 0:
        print("\n" + "="*50)
        print("CUTFLOW TABLE")
        print("="*50)

        cutflow = build_cutflow(seps, pointing, weights, momenta, p_soft,
                                dca, vtx_in, open_angle, sep_outer,
                                collinearity, on_tracker=on_tracker,
                                timing_chi2=timing_chi2)
        cutflow_df = pd.DataFrame(cutflow)
        print(f"\n{'Cut':<25} {'Eff (cumul.)':<15} {'Eff (margin.)':<15}")
        print("-" * 55)
        for row in cutflow:
            print(f"{row['cut']:<25} {row['efficiency']:<15.4f} "
                  f"{row['marginal_efficiency']:<15.4f}")

        cutflow_df.to_csv(os.path.join(OUT_DIR, 'cutflow_' + outString + '.csv'), index=False)
        print(f"\nCutflow saved to {OUT_DIR}/cutflow_{outString}.csv")

    # --- Decay position heatmap in tunnel cross-section ---
    # For events passing the full signal selection, project each decay
    # position onto the nearest tunnel-centreline segment and bin the
    # resulting (x_local, y_local) cross-section coordinates weighted by
    # the decay-probability weight. Shows where in the horseshoe the
    # accepted decays sit relative to the tunnel walls and fiducial inset.
    if len(mc.get('decay_pos', [])) > 0:
        print("\n" + "="*50)
        print("DECAY POSITION HEATMAP (TUNNEL CROSS-SECTION)")
        print("="*50)

        sel_mask = selection_mask(mc)
        if sel_mask.sum() > 0 and weights[sel_mask].sum() > 0:
            decay_sel = mc['decay_pos'][sel_mask]
            w_sel = weights[sel_mask]

            # Project each decay position onto the nearest centreline
            # segment to get the local cross-section (x, y) coordinates.
            m = len(decay_sel)
            best_d2 = np.full(m, np.inf)
            x_local = np.zeros(m)
            y_local = np.zeros(m)
            for i in range(len(path_3d_fiducial) - 1):
                seg = path_3d_fiducial[i + 1] - path_3d_fiducial[i]
                seg_len = np.linalg.norm(seg)
                if seg_len == 0:
                    continue
                seg_hat = seg / seg_len
                world_up = np.array([0., 1., 0.]) if abs(seg_hat[1]) < 0.9 \
                    else np.array([0., 0., 1.])
                right = np.cross(seg_hat, world_up)
                right /= np.linalg.norm(right)
                up = np.cross(right, seg_hat)
                up /= np.linalg.norm(up)
                rel = decay_sel - path_3d_fiducial[i]
                t = np.clip(rel @ seg_hat, 0, seg_len)
                closest = path_3d_fiducial[i] + np.outer(t, seg_hat)
                diff = decay_sel - closest
                d2 = np.einsum('ij,ij->i', diff, diff)
                upd = d2 < best_d2
                best_d2[upd] = d2[upd]
                x_local[upd] = diff[upd] @ right
                y_local[upd] = diff[upd] @ up

            outer = tunnel_profile_points(inset=0.0)
            inner = tunnel_profile_points(inset=DETECTOR_THICKNESS)
            outer_loop = np.vstack([outer, outer[:1]])
            inner_loop = np.vstack([inner, inner[:1]])

            x_lim = max(np.abs(outer[:, 0]).max(),
                        np.abs(x_local).max() if len(x_local) else 0) + 0.2
            y_lo = min(outer[:, 1].min(), y_local.min() if len(y_local) else 0) - 0.1
            y_hi = max(outer[:, 1].max(), y_local.max() if len(y_local) else 0) + 0.1

            fig_xs, ax_xs = plt.subplots(figsize=(8, 8))
            h = ax_xs.hist2d(
                x_local, y_local, bins=60,
                range=[[-x_lim, x_lim], [y_lo, y_hi]],
                weights=w_sel, cmap='magma', cmin=1e-30)
            ax_xs.plot(outer_loop[:, 0], outer_loop[:, 1],
                       color='white', linewidth=2, label='Tunnel wall')
            ax_xs.plot(inner_loop[:, 0], inner_loop[:, 1],
                       color='white', linewidth=1.2, linestyle='--',
                       label=f'Fiducial (inset {DETECTOR_THICKNESS*100:.0f} cm)')
            ax_xs.set_aspect('equal')
            ax_xs.set_xlabel('x_local (m)  [horizontal across profile]')
            ax_xs.set_ylabel('y_local (m)  [up]')
            ax_xs.set_title(
                f'Decay position in tunnel cross-section\n'
                f'(events passing full selection, τ = {lifetime*1e9:.0f} ns)')
            ax_xs.legend(loc='upper right', fontsize=9, framealpha=0.7)
            plt.colorbar(h[3], ax=ax_xs, label='Weighted yield (decay prob.)')
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, 'decay_xsec_heatmap_' + outString + '.png'), dpi=150)
            show_or_close()

            print(f"  Selected samples: {int(sel_mask.sum())}, "
                  f"weighted yield: {w_sel.sum():.4e}")
            print(f"  Saved: {OUT_DIR}/decay_xsec_heatmap_{outString}.png")
        else:
            print("  No samples pass the full selection — skipping plot.")

    # --- Save MC distributions for overlay plotting ---
    if len(seps) > 0:
        np.savez(os.path.join(OUT_DIR, 'mc_distributions_' + outString + '.npz'),
                 seps=seps, pointing=pointing, weights=weights,
                 momenta=momenta, p_soft=p_soft,
                 dca=dca, vtx_in=vtx_in, open_angle=open_angle,
                 sep_outer=sep_outer, d_implied=d_implied,
                 collinearity=collinearity, on_tracker=on_tracker,
                 lifetime=lifetime, label=outString,
                 p_cut=P_CUT, sep_min=SEP_MIN, sep_max=SEP_MAX,
                 dca_cut=DCA_CUT,
                 theta_parallel=THETA_PARALLEL,
                 sep_out_max_parallel=SEP_OUT_MAX_PARALLEL,
                 collin_min=COLLIN_MIN,
                 sep_out_collin_gate=SEP_OUT_COLLIN_GATE,
                 collin_frac=COLLIN_FRAC, sep_out_gate=SEP_OUT_GATE,
                 point_tight_sep_in=POINT_TIGHT_SEP_IN,
                 sep_in_point_gate=SEP_IN_POINT_GATE,
                 point_global=POINT_GLOBAL,
                 hit_resolution=HIT_RESOLUTION, n_layers=N_LAYERS)
        print(f"MC distributions saved to {OUT_DIR}/mc_distributions_{outString}.npz")

    # --- Event displays (3D + top-down + side, per-selection) ---
    # Edit ``display_selections`` to control which event populations get
    # rendered. Each entry is a per-sample boolean mask; the first
    # ``n_per`` events with at least one passing sample are drawn.
    if len(mc.get('decay_pos', [])) > 0:
        import event_display
        sel_full = selection_mask(mc)
        base_kin = (mc['p_soft'] >= P_CUT) \
            & (mc['sep_outer'] >= SEP_MIN) & (mc['sep'] <= SEP_MAX)
        display_selections = {
            'accepted':      sel_full,
            'off_tracker':   (~mc['on_tracker']) & base_kin,
            'fail_sep_in':   mc['on_tracker'] & base_kin & (mc['sep'] < SEP_MIN),
            'fail_dca':      mc['on_tracker'] & base_kin
                & (mc['sep'] >= SEP_MIN) & (mc['dca'] > DCA_CUT),
            'fail_vtx':      mc['on_tracker'] & base_kin
                & (mc['sep'] >= SEP_MIN) & ~mc['vtx_in'],
        }
        event_display.make_event_displays(
            mc, display_selections, n_per=10,
            out_prefix=f'event_display_{outString}',
            out_dir=EV_DIR,
        )

    # Lifetime scan
    print("\n" + "="*50)
    print("LIFETIME SCAN")
    print("="*50)
    
    lifetimes = np.logspace(-10.5, -3.5, 20)
    scan = analyze_decay_vs_lifetime(sample_csv, geo_cache, lifetimes)

    # MC-based exclusion: reweight the single uniform-sampled MC pass to each
    # lifetime, with the full reconstruction-level selection (DCA, PCA vertex,
    # collinearity, outer separation, smearing) and the scintillator veto
    # applied. This is the headline GRENDEL curve. The analytic 'scan' applies
    # only the acceptance-level cuts (momentum + min/max separation + forward
    # cap) under full-coverage tracking; it is the dashed overlay. The gap
    # between the two is dominated by the extra selection cuts, not the veto.
    mc_scan = mc_exclusion_vs_lifetime(mc, lifetimes, scan['total_events'])

    best_mc = np.nanmin(mc_scan['exclusion'])
    best_an = np.nanmin(scan['exclusion'])
    print(f"  Best excluded BR  — MC (full selection): {best_mc:.3e}")
    print(f"  Best excluded BR  — analytic (acceptance only): {best_an:.3e}")
    if best_mc > 0 and np.isfinite(best_mc):
        print(f"  Full selection vs acceptance-only at best point: "
              f"x{best_mc / best_an:.2f} weaker")

    # === Plotting ===
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    ax1 = axes[0, 0]
    ax1.semilogx(lifetimes * 1e9, scan['mean_acceptance'],
                 'b-', linewidth=2)
    ax1.set_xlabel('Lifetime (ns)')
    ax1.set_ylabel('Mean Acceptance')
    ax1.set_title(f'Two-body acceptance\n(p > {P_CUT*1000:.0f} MeV/c, '
                  f'{SEP_MIN*1000:.0f} mm < sep < {SEP_MAX*100:.0f} cm)')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.set_ylim(0, 1.05)
    
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
    
    ax4 = axes[1, 1]
    ax4.loglog(lifetimes * SPEED_OF_LIGHT, mc_scan['exclusion'],
               color='blue', linewidth=2, label="GRENDEL (full selection)")
    ax4.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion'],
               color='blue', linewidth=2, linestyle=':', alpha=0.5,
               label="GRENDEL (acceptance only)")
    ax4.set_xlabel(r'$c\tau$ (m)')
    ax4.set_ylabel('BR')
    ax4.grid(True, which="both", ls="-", alpha=0.2)
    
    ext = {}
    if df_results['mass'].iloc[0] == 15:
        ext["MATHUSLA"] = np.loadtxt("external/MATHUSLA.csv", delimiter=",")
        ext["CODEX"] = np.loadtxt("external/CODEX.csv", delimiter=",")
        ext["CMS"] = np.loadtxt("external/CMS.csv", delimiter=",")
        ext["ANUBISOpt"] = np.loadtxt("external/ANUBISPBC.csv", delimiter=",")

        ax4.loglog(ext["MATHUSLA"][:, 0], ext["MATHUSLA"][:, 1],
                   color="green", linewidth=2, linestyle="--", label="MATHUSLA")
        ax4.loglog(ext["CODEX"][:, 0], ext["CODEX"][:, 1],
                   color="cyan", linewidth=2, linestyle="--",label="CODEX-b")
        # ax4.loglog(ext["CMS"][:, 0], ext["CMS"][:, 1],
        #            color="purple", linewidth=2, label="CMS")
        ax4.loglog(ext["ANUBISOpt"][:, 0], ext["ANUBISOpt"][:, 1],
                    color="magenta", linewidth=2, linestyle="--", label="ANUBIS (PBC)")
        ax4.legend(fontsize=8, loc='lower right')
        ax4.set_ylim([1E-5,1])
    elif df_results['mass'].iloc[0] == 0.5:
        ext["CODEX"] = np.loadtxt("external/CODEX0p5_2.csv", delimiter=",")
    
        ax4.loglog(ext["CODEX"][:, 0], ext["CODEX"][:, 1],
                   color="cyan", linewidth=2,linestyle="--", label="CODEX-b")
        ax4.legend(fontsize=8, loc='lower right')
        
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'exclusion_2body'+outString+'.png'), dpi=150)
    show_or_close()
    
    df_results.to_csv(os.path.join(OUT_DIR, "particle_decay_results_2body.csv"), index=False)
    event_df.to_csv(os.path.join(OUT_DIR, "event_decay_statistics_2body.csv"), index=False)
    print("\nResults saved.")
    print(f"Plots in {OUT_DIR}/: exclusion_2body{outString}.png, "
          f"separation_histogram{outString}.png, ...")
    print(f"Event displays in {EV_DIR}/")
