import numpy as np
import pandas as pd
import trimesh
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import quad

# Speed of light in m/s
SPEED_OF_LIGHT = 299792458.0  # m/s
M_ELECTRON = 0.000511  # GeV/c²

# ============================================================
# Tunnel cross-section definition (Position C measurements)
# ============================================================
TUNNEL_ALPHA = 2.90  # floor width (m)
TUNNEL_BETA  = 3.15  # total height (m)
TUNNEL_GAMMA = 2.90  # arch width at springline (m)
TUNNEL_DELTA = 1.90  # arch height (m)
TUNNEL_WALL_HEIGHT = TUNNEL_BETA - TUNNEL_DELTA  # 1.25m

DETECTOR_THICKNESS = 0.24  # 24 cm

# Analysis cuts
P_CUT   = 0.600    # GeV/c — minimum electron momentum
SEP_MIN = 0.001    # m — minimum separation at detector (1 mm)
SEP_MAX = 1.0     # m — maximum separation at detector (10 cm)

outString = "1GeV"

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
# Geometry functions
# ============================================================

def tunnel_profile_points(n_arch=32, n_wall=4, inset=0.0, inset_floor=False):
    half_w = TUNNEL_GAMMA / 2 - inset
    half_floor = TUNNEL_ALPHA / 2 - inset
    wall_h = TUNNEL_WALL_HEIGHT
    a = TUNNEL_GAMMA / 2 - inset
    b = TUNNEL_DELTA - inset
    floor_y = inset if inset_floor else 0.0
    points = []
    points.append([-half_floor, floor_y])
    points.append([half_floor, floor_y])
    for i in range(1, n_wall + 1):
        frac = i / n_wall
        y = floor_y + (wall_h - floor_y) * frac if inset_floor else frac * wall_h
        x = half_floor + (half_w - half_floor) * frac
        points.append([x, y])
    for i in range(1, n_arch):
        angle = np.pi * i / n_arch
        x = a * np.cos(angle)
        y = wall_h + b * np.sin(angle)
        points.append([x, y])
    for i in range(n_wall, 0, -1):
        frac = i / n_wall
        y = floor_y + (wall_h - floor_y) * frac if inset_floor else frac * wall_h
        x = half_floor + (half_w - half_floor) * frac
        points.append([-x, y])
    points = np.array(points)
    rect_area = TUNNEL_ALPHA * TUNNEL_WALL_HEIGHT
    rect_cy = TUNNEL_WALL_HEIGHT / 2
    a0 = TUNNEL_GAMMA / 2
    b0 = TUNNEL_DELTA
    ellipse_area = np.pi * a0 * b0 / 2
    ellipse_cy = TUNNEL_WALL_HEIGHT + 4 * b0 / (3 * np.pi)
    total_area = rect_area + ellipse_area
    centroid_y = (rect_area * rect_cy + ellipse_area * ellipse_cy) / total_area
    points[:, 1] -= centroid_y
    return points


def create_profile_mesh(path_points, profile_2d):
    n_profile = len(profile_2d)
    vertices = []
    faces = []
    for i in range(len(path_points)):
        if i == 0:
            tangent = path_points[1] - path_points[0]
        elif i == len(path_points) - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        tangent = tangent / np.linalg.norm(tangent)
        if abs(tangent[2]) < 0.9:
            world_up = np.array([0, 0, 1])
        else:
            world_up = np.array([1, 0, 0])
        right = np.cross(tangent, world_up)
        right = right / np.linalg.norm(right)
        up = np.cross(right, tangent)
        up = up / np.linalg.norm(up)
        for j in range(n_profile):
            offset = profile_2d[j, 0] * right + profile_2d[j, 1] * up
            vertices.append(path_points[i] + offset)
        if i > 0:
            for j in range(n_profile):
                v1 = (i-1) * n_profile + j
                v2 = (i-1) * n_profile + (j + 1) % n_profile
                v3 = i * n_profile + (j + 1) % n_profile
                v4 = i * n_profile + j
                faces.append([v1, v4, v3])
                faces.append([v1, v3, v2])
    center_start = len(vertices)
    vertices.append(path_points[0].copy())
    for j in range(n_profile):
        faces.append([center_start, (j + 1) % n_profile, j])
    center_end = len(vertices)
    vertices.append(path_points[-1].copy())
    last = (len(path_points) - 1) * n_profile
    for j in range(n_profile):
        faces.append([center_end, last + j, last + (j + 1) % n_profile])
    return np.array(vertices), np.array(faces)


def eta_phi_to_direction(eta, phi):
    theta = 2 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    return np.array([dx, dy, dz])


def calculate_decay_length(momentum, mass, lifetime):
    energy = np.sqrt(momentum**2 + mass**2)
    beta = momentum / energy
    gamma = energy / mass
    return gamma * beta * SPEED_OF_LIGHT * lifetime


# ============================================================
# Geometry caching and processing
# ============================================================

def cache_geometry(csv_file, mesh, origin):
    df = pd.read_csv(csv_file)
    n = len(df)
    origin = np.array(origin)
    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)
    momentum = df['momentum'].values
    mass = df['mass'].values
    energy = np.sqrt(momentum**2 + mass**2)
    gamma = energy / mass
    beta = momentum / energy
    
    print(f"Caching fiducial volume geometry for {n} particles...")
    for idx, row in tqdm(df.iterrows(), total=n, desc="Ray-casting"):
        direction = eta_phi_to_direction(row['eta'], row['phi'])
        locations, _, _ = mesh.ray.intersects_location(
            ray_origins=[origin], ray_directions=[direction])
        if len(locations) >= 2:
            hits[idx] = True
            distances = sorted([np.linalg.norm(loc - origin) for loc in locations])
            entry_d[idx] = distances[0]
            exit_d[idx] = distances[1]
    
    n_hits = hits.sum()
    print(f"  {n_hits} / {n} particles hit fiducial volume ({n_hits/n*100:.1f}%)")
    if n_hits > 0:
        print(f"  Mean path length: {(exit_d[hits] - entry_d[hits]).mean():.2f} m")
    
    return {
        'hits': hits, 'entry_d': entry_d, 'exit_d': exit_d,
        'gamma': gamma, 'beta': beta, 'momentum': momentum, 'mass': mass
    }


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
        results['exclusion'].append(3 / (mean_p1 * 3000 * 52E3))
        results['exclusion_no_cuts'].append(3 / (mean_p1_nc * 3000 * 52E3))
        results['mean_acceptance'].append(mean_acc)
    
    return results


def sample_separations(geo_cache, lifetime_seconds, n_samples_per_particle=100,
                       rng_seed=42):
    """
    Monte Carlo sample decay positions and rest-frame angles to build
    a distribution of electron-pair separations at the detector.
    
    For each particle that hits the fiducial volume:
      1. Sample decay position d from (1/λ) exp(-d/λ) within [entry, exit]
         using inverse CDF sampling
      2. Sample |cosθ*| uniformly in [0, β] (forward requirement)
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
# Tunnel centerline
# ============================================================

correctedVert = [
(-86.57954338701529, 0.1882163986665546  ),
(-1731.590867740335, 3.764327973349282   ),
(-3549.761278867689, 7.716872345365118   ),
(-5887.408950317142, 12.798715109387558  ),
(-8053.403266181902, -504.23173203003535 ),
(-10046.991360867298, -1282.5065405198511),
(-11783.350377373874, -2930.9057600491833),
(-12913.652590171332, -4580.622494369192 ),
(-13095.344153684957, -7536.749251839814 ),
(-13099.610392054752, -9015.000846973791 ),
(-13278.792403586143, -11101.567842600896),
(-13372.39869252341, -13536.146959364076 ),
(-13292.093029091975, -15710.234580371536),
(-12779.140603923677, -17972.21925955668   ), 
(-11659.12755425337, -19887.69754879509    ),
(-10105.714877251532, -21630.204967658145  ),
(-7512.845769209047, -23201.0590309365     ),
(-5262.530506741277, -23466.820585854904   ),
(-2751.72374851779, -23472.278861416264    ),
(-241.41890069074725, -23651.64908934632   ),
(1749.6596420124115, -23742.93404270002    ),
(3827.568683300815, -23747.45123626804     ),
(6078.6368113632525, -23752.344862633392   ),
(8502.613071001502, -23844.570897980426    ),
(11446.568501358292, -23764.01427935077    ),
(13438.399909656131, -23594.431304151418   ),
(15777.051401898476, -23251.689242178036   ),
(18289.614846509525, -22648.455684448927   ),
(20889.761655300477, -21697.58643838109    ),
(23143.841245741598, -20659.00835053422    ),
(25486.006110759066, -19098.88262197991    ),
(27742.09334278597, -17364.656724658227    ),
(28871.391734790544, -16062.763895075637   ),
(30781.662703665817, -14153.873179790575   ),
(32518.021720172394, -12505.473960261239   ),
(34513.49197884447, -11075.029330388788    ),
(36636.57295581305, -10427.47081077351     ),
(38759.40297758341, -9866.868267342572     ),
(41357.416667189485, -9655.12481884172     ),
(43694.93886103982, -9703.684649697909     ),
(46379.03018363646, -9666.041369964427     ),
(49409.43967978114, -9629.150955825604     ),
(51660.88424064092, -9503.610617914434     ),
(54258.0195870532, -9596.213086058811      ),
(57028.564975437745, -9602.236010816167    ),
(59539.87364405768, -9433.782334008818     ),
(62050.42944708294, -9526.196585754526     )]

correctedVertWithShift = []
for x, y in correctedVert:
    correctedVertWithShift.append(
        ((x - 11908.8279764855) / 1000, (y + 13591.106147774964) / 1000))

Z_POSITION = 22
path_3d = np.array([[x, Z_POSITION,y] for x, y in correctedVertWithShift])

# Build fiducial volume mesh
print("Building fiducial volume mesh...")
profile_fiducial = tunnel_profile_points(inset=DETECTOR_THICKNESS, inset_floor=False)
verts, faces = create_profile_mesh(path_3d, profile_fiducial)
mesh_fiducial = trimesh.Trimesh(vertices=verts, faces=faces)
if mesh_fiducial.volume < 0:
    mesh_fiducial.invert()

print(f"  Fiducial volume: {mesh_fiducial.volume:.1f} m³")
print(f"  Detector thickness: {DETECTOR_THICKNESS*100:.0f} cm")
print(f"  Cuts: p_e > {P_CUT*1000:.0f} MeV/c, "
      f"{SEP_MIN*1000:.0f} mm < separation < {SEP_MAX*100:.0f} cm")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    sample_csv = "LLP.csv"
    origin = [0, 0, 0]
    
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
        d_min_sep = SEP_MIN / theta_min   # distance where min sep = θ_min
        d_max_sep = SEP_MAX / theta_min   # distance where max sep = θ_min
        # At 1m remaining, what is max sep from θ_12(cosθ*=0)?
        sep_at_1m = theta_min * 1.0
        print(f"{p_test:>8.0f} {g:>6.1f} {c_P:>7.4f} {theta_min*1000:>8.1f} mrad "
              f"{d_min_sep*100:>7.2f}cm {d_max_sep:>7.1f}m {sep_at_1m*100:>12.1f}cm")
    
    # Show how max sep constrains high-γ particles
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
    
    seps, weights, momenta = sample_separations(
        geo_cache, lifetime, n_samples_per_particle=200)
    
    if len(seps) > 0:
        fig_sep, axes_sep = plt.subplots(1, 3, figsize=(18, 5))
        
        # (a) All separations, decay-probability weighted
        ax = axes_sep[0]
        bins = np.linspace(0, 0.5, 100)  # 0 to 50 cm
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
        plt.savefig('separation_histogram'+outString+'.png', dpi=150)
        plt.show()
        
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
    
    lifetimes = np.logspace(-9.5, -4.5, 20)
    scan = analyze_decay_vs_lifetime(sample_csv, geo_cache, lifetimes)
    
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
               color='blue', linewidth=2, label="milliQan (with cuts)")
    ax4.loglog(lifetimes * SPEED_OF_LIGHT, scan['exclusion_no_cuts'],
               color='blue', linewidth=2, linestyle='--', alpha=0.5,
               label="milliQan (no cuts)")
    ax4.set_xlabel(r'$c\tau$ (m)')
    ax4.set_ylabel('BR')
    ax4.grid(True, which="both", ls="-", alpha=0.2)
    
    ext = {}
    ext["MATHUSLA"] = np.loadtxt("external/MATHUSLA.csv", delimiter=",")
    ext["CODEX"] = np.loadtxt("external/CODEX.csv", delimiter=",")
    ext["ANUBIS"] = np.loadtxt("external/ANUBIS.csv", delimiter=",")
    ext["ANUBISOpt"] = np.loadtxt("external/ANUBISOpt.csv", delimiter=",")
    ext["ANUBISCons"] = np.loadtxt("external/ANUBISUpdateCons.csv", delimiter=",")
    
    ax4.loglog(ext["MATHUSLA"][:, 0], ext["MATHUSLA"][:, 1],
               color="green", linewidth=2, label="MATHUSLA")
    ax4.loglog(ext["CODEX"][:, 0], ext["CODEX"][:, 1],
               color="cyan", linewidth=2, label="CODEX-b")
    ax4.loglog(ext["ANUBIS"][:, 0], ext["ANUBIS"][:, 1],
               color="purple", linewidth=2, label="ANUBIS")
    ax4.loglog(ext["ANUBISOpt"][:, 0], ext["ANUBISOpt"][:, 1],
               color="purple", linewidth=2, linestyle="--", label="ANUBIS Opt")
    ax4.loglog(ext["ANUBISCons"][:, 0], ext["ANUBISCons"][:, 1],
               color="magenta", linewidth=2, linestyle="--", label="ANUBIS Cons")
    ax4.legend(fontsize=8, loc='upper right')
    
    plt.tight_layout()
    plt.savefig('exclusion_2body'+outString+'.png', dpi=150)
    plt.show()
    
    df_results.to_csv("particle_decay_results_2body.csv", index=False)
    event_df.to_csv("event_decay_statistics_2body.csv", index=False)
    print("\nResults saved.")
    print("Plots: exclusion_2body"+outString+".png, separation_histogram"+outString+".png")
