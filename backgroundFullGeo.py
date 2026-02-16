"""
Beam Muon Trident Background MC — Full Tunnel Geometry
=======================================================

Uses the identical tunnel mesh and ray-casting as the signal code.
Muons originate from the IP (origin), traverse through rock + CMS,
and enter the tunnel fiducial volume at various (eta, phi).

The muon flux is assumed to depend only on distance from IP,
anchored to the milliQan measurement:
  Phi = 0.3e-3 pb^-1 cm^-2 = 0.3 fb^-1 cm^-2  at r = 33 m

For a muon at distance r:  Phi(r) = Phi_mQ × (33/r)^2

Approach:
  1. Build the tunnel fiducial mesh (identical to signal code)
  2. Sample muon (eta, phi) isotropically in the solid angle 
     subtended by the tunnel
  3. Ray-cast from IP to get entry/exit points and path length
  4. Weight each muon by the flux at its entry distance
  5. Compute trident probability along the air path
  6. For each trident, sample pair kinematics and apply cuts
  7. Sum weights → background rate per fb^-1

Compatible with: decayProbPerEvent_2body.py
"""

import numpy as np
import trimesh
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import quad
from scipy.interpolate import interp1d
from tqdm import tqdm
from dataclasses import dataclass

# =============================================================================
# Constants (same as signal code + background code)
# =============================================================================
SPEED_OF_LIGHT = 299792458.0
ALPHA_EM       = 1.0 / 137.036
R_E            = 2.8179e-15    # m
M_ELECTRON     = 0.000511      # GeV
M_MUON         = 0.10566       # GeV
N_A            = 6.022e23

# milliQan flux
MILLIQAN_FLUX_PER_FB_CM2 = 0.3   # muons per fb^-1 per cm^2 at r=33m
MILLIQAN_DISTANCE = 33.0          # m
E_THRESHOLD = 15.0                # GeV production threshold

# Analysis cuts
E_CUT   = 0.600   # GeV
SEP_MIN = 0.001   # m (1 cm)

# Tunnel geometry constants (from signal code)
TUNNEL_ALPHA = 2.90
TUNNEL_BETA  = 3.15
TUNNEL_GAMMA = 2.90
TUNNEL_DELTA = 1.90
TUNNEL_WALL_HEIGHT = TUNNEL_BETA - TUNNEL_DELTA
DETECTOR_THICKNESS = 0.24

FIDUCIAL_WIDTH = TUNNEL_GAMMA - 2 * DETECTOR_THICKNESS
SEP_MAX = FIDUCIAL_WIDTH  # ~2.42 m

# Air properties
AIR_Z   = 7.3
AIR_A   = 14.5
AIR_RHO = 1.205e-3  # g/cm^3
AIR_N   = N_A * AIR_RHO / AIR_A  # atoms/cm^3


# =============================================================================
# Tunnel geometry (copied from signal code)
# =============================================================================

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
            tangent = path_points[i] - path_points[i - 1]
        else:
            tangent = path_points[i + 1] - path_points[i - 1]
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
                v1 = (i - 1) * n_profile + j
                v2 = (i - 1) * n_profile + (j + 1) % n_profile
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


def build_tunnel_mesh():
    """Build the fiducial volume mesh (identical to signal code)."""
    correctedVert = [
        (-86.57954338701529, 0.1882163986665546),
        (-1731.590867740335, 3.764327973349282),
        (-3549.761278867689, 7.716872345365118),
        (-5887.408950317142, 12.798715109387558),
        (-8053.403266181902, -504.23173203003535),
        (-10046.991360867298, -1282.5065405198511),
        (-11783.350377373874, -2930.9057600491833),
        (-12913.652590171332, -4580.622494369192),
        (-13095.344153684957, -7536.749251839814),
        (-13099.610392054752, -9015.000846973791),
        (-13278.792403586143, -11101.567842600896),
        (-13372.39869252341, -13536.146959364076),
        (-13292.093029091975, -15710.234580371536),
        (-12779.140603923677, -17972.21925955668),
        (-11659.12755425337, -19887.69754879509),
        (-10105.714877251532, -21630.204967658145),
        (-7512.845769209047, -23201.0590309365),
        (-5262.530506741277, -23466.820585854904),
        (-2751.72374851779, -23472.278861416264),
        (-241.41890069074725, -23651.64908934632),
        (1749.6596420124115, -23742.93404270002),
        (3827.568683300815, -23747.45123626804),
        (6078.6368113632525, -23752.344862633392),
        (8502.613071001502, -23844.570897980426),
        (11446.568501358292, -23764.01427935077),
        (13438.399909656131, -23594.431304151418),
        (15777.051401898476, -23251.689242178036),
        (18289.614846509525, -22648.455684448927),
        (20889.761655300477, -21697.58643838109),
        (23143.841245741598, -20659.00835053422),
        (25486.006110759066, -19098.88262197991),
        (27742.09334278597, -17364.656724658227),
        (28871.391734790544, -16062.763895075637),
        (30781.662703665817, -14153.873179790575),
        (32518.021720172394, -12505.473960261239),
        (34513.49197884447, -11075.029330388788),
        (36636.57295581305, -10427.47081077351),
        (38759.40297758341, -9866.868267342572),
        (41357.416667189485, -9655.12481884172),
        (43694.93886103982, -9703.684649697909),
        (46379.03018363646, -9666.041369964427),
        (49409.43967978114, -9629.150955825604),
        (51660.88424064092, -9503.610617914434),
        (54258.0195870532, -9596.213086058811),
        (57028.564975437745, -9602.236010816167),
        (59539.87364405768, -9433.782334008818),
        (62050.42944708294, -9526.196585754526),
    ]

    correctedVertWithShift = []
    for x, y in correctedVert:
        correctedVertWithShift.append(
            ((x - 11908.8279764855) / 1000, (y + 13591.106147774964) / 1000))

    Z_POSITION = 22
    path_3d = np.array([[x, y, Z_POSITION] for x, y in correctedVertWithShift])

    profile_fiducial = tunnel_profile_points(inset=DETECTOR_THICKNESS, inset_floor=False)
    verts, faces = create_profile_mesh(path_3d, profile_fiducial)
    mesh = trimesh.Trimesh(vertices=verts, faces=faces)
    if mesh.volume < 0:
        mesh.invert()

    return mesh, path_3d


# =============================================================================
# Muon flux model (distance-dependent, anchored to milliQan)
# =============================================================================

class MuonFlux:
    """
    Muon flux depending only on distance from IP.
    
    Phi(r) = Phi_mQ × (r_mQ / r)^2
    
    The energy spectrum shape is a power law with exponential cutoff,
    normalized to the total flux.
    """
    def __init__(self):
        self._alpha = 2.7
        self._E0 = 200.0
        self._norm = self._compute_norm()

    def _spectrum_shape(self, E):
        if E <= E_THRESHOLD:
            return 0.0
        return E**(-self._alpha) * np.exp(-E / self._E0)

    def _compute_norm(self):
        integral, _ = quad(self._spectrum_shape, E_THRESHOLD, 2000, limit=200)
        return MILLIQAN_FLUX_PER_FB_CM2 / integral if integral > 0 else 0.0

    def total_flux_at_r(self, r):
        """Total muon flux (per fb^-1 per cm^2) at distance r from IP."""
        return MILLIQAN_FLUX_PER_FB_CM2 * (MILLIQAN_DISTANCE / r)**2

    def differential_flux(self, E, r):
        """dPhi/dE at distance r.  Units: muons / GeV / cm^2 / fb^-1."""
        scalar = np.ndim(E) == 0
        E = np.atleast_1d(np.array(E, dtype=float))
        flux = np.zeros_like(E)
        mask = E > E_THRESHOLD
        geo = (MILLIQAN_DISTANCE / r)**2
        flux[mask] = self._norm * geo * E[mask]**(-self._alpha) * np.exp(-E[mask] / self._E0)
        return float(flux[0]) if scalar else flux

    def mean_energy(self):
        num, _ = quad(lambda E: E * self._spectrum_shape(E), E_THRESHOLD, 2000, limit=200)
        den, _ = quad(self._spectrum_shape, E_THRESHOLD, 2000, limit=200)
        return num / den if den > 0 else 0

    def sample_energy(self, n, rng):
        """Inverse-CDF sampling of muon energies."""
        E_grid = np.linspace(E_THRESHOLD, 500, 2000)
        pdf = np.array([self._spectrum_shape(e) for e in E_grid])
        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        inv = interp1d(cdf, E_grid, bounds_error=False, fill_value=(E_THRESHOLD, 500))
        return inv(rng.uniform(0, 1, n))


# =============================================================================
# Trident physics
# =============================================================================

def trident_cross_section_air(E_muon):
    """Trident cross section per air atom (cm^2)."""
    if E_muon <= 2 * M_ELECTRON + M_MUON:
        return 0.0
    Z = AIR_Z
    L = np.log(183.0 / Z**(1.0/3))
    log_ratio = np.log(E_muon / M_MUON)
    return (ALPHA_EM**3) * Z**2 * (R_E * 100)**2 * (28.0/9.0) * L * log_ratio**2 / np.pi


def trident_probability(E_muon, path_length_m):
    """P(trident) for a muon traversing path_length_m of air."""
    sigma = trident_cross_section_air(E_muon)
    return AIR_N * sigma * (path_length_m * 100)  # path in cm


# =============================================================================
# MC Engine
# =============================================================================

class TunnelBackgroundMC:
    """
    Full-geometry MC for beam muon trident backgrounds.
    
    Strategy:
      1. Sample muon directions isotropically over a bounding solid angle
      2. Ray-cast from IP → fiducial mesh to find entry/exit
      3. For muons that hit: weight by flux at entry distance
      4. Sample trident production along the air path
      5. Simulate pair kinematics, apply cuts
      6. Sum flux-weighted events → rate per fb^-1
    """
    
    def __init__(self, mesh, path_3d, origin=np.array([0, 0, 0]), seed=42):
        self.mesh = mesh
        self.path_3d = path_3d
        self.origin = np.array(origin, dtype=float)
        self.rng = np.random.default_rng(seed)
        self.flux = MuonFlux()
        
        # Determine the bounding solid angle of the tunnel seen from the IP.
        # Use the mesh bounding box to define an eta-phi window.
        self._setup_angular_bounds()
    
    def _setup_angular_bounds(self):
        """
        Compute the eta-phi bounding box of the tunnel as seen from IP.
        This lets us sample muon directions efficiently.
        """
        # Use tunnel centerline path_3d to get angular range
        vecs = self.path_3d - self.origin
        r = np.linalg.norm(vecs, axis=1)
        
        # Convert to eta, phi
        theta = np.arccos(np.clip(vecs[:, 2] / r, -1, 1))
        eta = -np.log(np.tan(theta / 2 + 1e-10))
        phi = np.arctan2(vecs[:, 1], vecs[:, 0])
        
        # Add margin for the tunnel cross-section (~1.5 m radius at ~30 m)
        margin_angle = np.arctan(2.0 / 30.0)  # ~4 degrees
        margin_eta = margin_angle / np.sin(np.mean(theta))  # approximate
        
        self.eta_min = eta.min() - abs(margin_eta)
        self.eta_max = eta.max() + abs(margin_eta)
        self.phi_min = phi.min() - margin_angle
        self.phi_max = phi.max() + margin_angle
        
        # Solid angle of the bounding box (approximate)
        d_eta = self.eta_max - self.eta_min
        d_phi = self.phi_max - self.phi_min
        # dΩ = dη dφ (in eta-phi space, this is exact for the rapidity measure)
        self.solid_angle_eta_phi = d_eta * d_phi
        
        print(f"  Tunnel angular range from IP:")
        print(f"    eta: [{self.eta_min:.3f}, {self.eta_max:.3f}]")
        print(f"    phi: [{self.phi_min:.4f}, {self.phi_max:.4f}] rad "
              f"([{np.degrees(self.phi_min):.2f}, {np.degrees(self.phi_max):.2f}] deg)")
        print(f"    Bounding deta×dphi = {d_eta:.3f} × {d_phi:.4f} = {self.solid_angle_eta_phi:.4f}")
    
    def _eta_phi_to_direction(self, eta, phi):
        """Convert (eta, phi) to unit direction vector."""
        theta = 2 * np.arctan(np.exp(-eta))
        dx = np.sin(theta) * np.cos(phi)
        dy = np.sin(theta) * np.sin(phi)
        dz = np.cos(theta)
        return np.array([dx, dy, dz])
    
    def _sample_pair_forced(self, E_muon, path_length):
        """
        Sample trident e+e- pair kinematics (forced production).
        
        Every call produces a pair. The caller is responsible for
        weighting the event by P(trident).
        
        Args:
            E_muon: muon energy (GeV)
            path_length: total air path through fiducial volume (m)
        
        Returns: dict of pair properties
        """
        # Vertex position: uniform along the path
        frac = self.rng.uniform()
        d_remaining = path_length * (1.0 - frac)
        
        # Sample pair energy fraction v = E_pair / E_muon
        v = self._sample_v(E_muon)
        E_pair = v * E_muon
        
        # Sample asymmetry rho: dσ/drho ~ 1/(1-rho^2) roughly
        rho = self.rng.uniform(-0.85, 0.85)
        E_plus = E_pair * (1 + rho) / 2
        E_minus = E_pair * (1 - rho) / 2
        
        # Opening angle between e+ and e-
        theta_open = 2 * M_ELECTRON / E_pair / max(np.sqrt(1 - rho**2), 0.01)
        
        # Each electron's angle from muon direction
        theta_eplus = M_ELECTRON / max(E_plus, M_ELECTRON)
        theta_eminus = M_ELECTRON / max(E_minus, M_ELECTRON)
        
        # Separation at the wall (using remaining distance)
        separation = theta_open * d_remaining
        
        # Muon-electron separations
        sep_mu_eplus = theta_eplus * d_remaining
        sep_mu_eminus = theta_eminus * d_remaining
        
        # Apply cuts
        passes_E = (E_plus > E_CUT) and (E_minus > E_CUT)
        passes_sep_min = separation > SEP_MIN
        passes_sep_max = separation < SEP_MAX
        passes_muon_sep = sep_mu_eplus < SEP_MIN or sep_mu_eminus < SEP_MIN
        passes_all = passes_E and passes_sep_min and passes_sep_max and passes_muon_sep
        
        return {
            'E_muon': E_muon,
            'E_pair': E_pair, 'E_plus': E_plus, 'E_minus': E_minus,
            'v': v, 'rho': rho,
            'theta_opening_mrad': theta_open * 1000,
            'path_length': path_length,
            'd_remaining': d_remaining,
            'separation_mm': separation * 1000,
            'sep_mu_eplus_mm': sep_mu_eplus * 1000,
            'sep_mu_eminus_mm': sep_mu_eminus * 1000,
            'passes_energy': passes_E,
            'passes_sep_min': passes_sep_min,
            'passes_sep_max': passes_sep_max,
            'passes_muon_sep': passes_muon_sep,
            'passes_all': passes_all,
        }
    
    def _sample_v(self, E_muon):
        """Rejection-sample v from dσ/dv ~ (1/v)(1 - 4v/3 + v^2)."""
        v_min = 2 * M_ELECTRON / E_muon
        v_max = 0.95
        while True:
            v = v_min * (v_max / v_min)**self.rng.uniform()
            if self.rng.uniform() < (1.0 - 4.0*v/3.0 + v**2):
                return v
    
    def run(self, n_rays):
        """
        Run the full MC with forced trident production.
        
        Every muon that hits the tunnel produces exactly one trident,
        weighted by P(trident). This is far more efficient than 
        Poisson sampling when P(trident) ~ 10^-6.
        
        Args:
            n_rays: number of muon directions to sample
        
        Returns: dict with results
        """
        print(f"\nRunning full-geometry background MC (forced trident)...")
        print(f"  Sampling {n_rays:,} muon directions in bounding box")
        
        # Sample eta, phi uniformly in the bounding box
        eta_samples = self.rng.uniform(self.eta_min, self.eta_max, n_rays)
        phi_samples = self.rng.uniform(self.phi_min, self.phi_max, n_rays)
        
        # Sample muon energies
        E_samples = self.flux.sample_energy(n_rays, self.rng)
        
        # Ray-cast all muons
        n_hit = 0
        n_pass = 0
        
        flux_times_r2 = MILLIQAN_FLUX_PER_FB_CM2 * 1e4 * MILLIQAN_DISTANCE**2
        
        d_eta = self.eta_max - self.eta_min
        d_phi = self.phi_max - self.phi_min
        solid_angle_per_ray_base = d_eta * d_phi / n_rays
        
        hit_data = []
        trident_data = []
        
        weighted_total_muons = 0.0
        weighted_trident = 0.0
        weighted_pass = 0.0
        
        for i in tqdm(range(n_rays), desc="Ray-casting"):
            eta = eta_samples[i]
            phi = phi_samples[i]
            E_mu = E_samples[i]
            
            direction = self._eta_phi_to_direction(eta, phi)
            
            locations, _, _ = self.mesh.ray.intersects_location(
                ray_origins=[self.origin], ray_directions=[direction])
            
            if len(locations) < 2:
                continue
            
            distances = sorted([np.linalg.norm(loc - self.origin) for loc in locations])
            entry_d = distances[0]
            exit_d = distances[1]
            path_length = exit_d - entry_d
            
            if path_length < 0.01:
                continue
            
            n_hit += 1
            
            theta = 2 * np.arctan(np.exp(-eta))
            sin2_theta = np.sin(theta)**2
            ray_weight = flux_times_r2 * solid_angle_per_ray_base * sin2_theta
            
            weighted_total_muons += ray_weight
            
            hit_data.append({
                'eta': eta, 'phi': phi, 'E_mu': E_mu,
                'entry_d': entry_d, 'exit_d': exit_d,
                'path_length': path_length,
                'weight': ray_weight,
            })
            
            # Forced trident: always produce a pair, weight by P(trident)
            P_tri = trident_probability(E_mu, path_length)
            pair = self._sample_pair_forced(E_mu, path_length)
            
            pair['eta'] = eta
            pair['phi'] = phi
            pair['entry_d'] = entry_d
            pair['exit_d'] = exit_d
            pair['weight'] = ray_weight * P_tri
            
            trident_data.append(pair)
            weighted_trident += pair['weight']
            
            if pair['passes_all']:
                n_pass += 1
                weighted_pass += pair['weight']
        
        print(f"\n  Rays sampled:    {n_rays:,}")
        print(f"  Hit tunnel:      {n_hit:,} ({n_hit/n_rays*100:.2f}%)")
        print(f"  Trident events:  {n_hit} (forced, 1 per muon)")
        print(f"  Pass all cuts:   {n_pass}")
        
        print(f"\n  Weighted muons through tunnel: {weighted_total_muons:.2e} per fb^-1")
        print(f"  Weighted trident (all):        {weighted_trident:.2e} per fb^-1")
        print(f"  Weighted trident (pass cuts):  {weighted_pass:.2e} per fb^-1")
        
        print(f"\n  Cross-check:")
        print(f"    Mean entry distance:  {np.mean([h['entry_d'] for h in hit_data]):.1f} m")
        print(f"    Mean path length:     {np.mean([h['path_length'] for h in hit_data]):.2f} m")
        
        return {
            'n_rays': n_rays,
            'n_hit': n_hit,
            'n_trident': n_hit,
            'n_pass': n_pass,
            'weighted_muons_per_fb': weighted_total_muons,
            'weighted_trident_per_fb': weighted_trident,
            'weighted_pass_per_fb': weighted_pass,
            'hit_data': hit_data,
            'trident_data': trident_data,
        }


# =============================================================================
# Visualization
# =============================================================================

def plot_results(results, save_path=None):
    """Comprehensive results plot."""
    
    hit_data = results['hit_data']
    trident_data = results['trident_data']
    passing = [t for t in trident_data if t['passes_all']]
    
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    fig.suptitle(
        f'Full-Geometry Beam Muon Trident MC\n'
        f'{results["n_rays"]:,} rays → {results["n_hit"]:,} hit → '
        f'{results["n_trident"]} tridents → {results["n_pass"]} pass cuts\n'
        f'Weighted: {results["weighted_pass_per_fb"]:.2e} per fb$^{{-1}}$ (before muon veto)',
        fontsize=13, fontweight='bold')
    
    # --- 1: eta-phi map of tunnel hits ---
    ax1 = fig.add_subplot(gs[0, 0])
    if hit_data:
        eta_h = [h['eta'] for h in hit_data]
        phi_h = [np.degrees(h['phi']) for h in hit_data]
        ax1.scatter(eta_h, phi_h, s=1, alpha=0.3, color='blue', label='Hits')
        if trident_data:
            eta_t = [t['eta'] for t in trident_data]
            phi_t = [np.degrees(t['phi']) for t in trident_data]
            ax1.scatter(eta_t, phi_t, s=20, color='orange', zorder=5, label='Trident')
        if passing:
            eta_p = [t['eta'] for t in passing]
            phi_p = [np.degrees(t['phi']) for t in passing]
            ax1.scatter(eta_p, phi_p, s=60, marker='*', color='red', 
                       zorder=10, label='Pass cuts')
    ax1.set_xlabel(r'$\eta$'); ax1.set_ylabel(r'$\phi$ (deg)')
    ax1.set_title(r'Tunnel in $\eta$–$\phi$ from IP')
    ax1.legend(fontsize=8, markerscale=2)
    
    # --- 2: Path length distribution ---
    ax2 = fig.add_subplot(gs[0, 1])
    if hit_data:
        paths = [h['path_length'] for h in hit_data]
        ax2.hist(paths, bins=50, alpha=0.7, color='steelblue')
        ax2.axvline(np.mean(paths), color='red', ls='--',
                    label=f'Mean = {np.mean(paths):.2f} m')
        ax2.axvline(FIDUCIAL_WIDTH, color='orange', ls=':',
                    label=f'Fiducial width = {FIDUCIAL_WIDTH:.2f} m')
    ax2.set_xlabel('Path Length Through Air (m)')
    ax2.set_ylabel('Muon Rays')
    ax2.set_title('Air Path Length Distribution')
    ax2.legend(fontsize=8)
    
    # --- 3: Entry distance distribution ---
    ax3 = fig.add_subplot(gs[0, 2])
    if hit_data:
        entry_ds = [h['entry_d'] for h in hit_data]
        ax3.hist(entry_ds, bins=50, alpha=0.7, color='steelblue')
        ax3.axvline(MILLIQAN_DISTANCE, color='red', ls='--',
                    label=f'milliQan = {MILLIQAN_DISTANCE} m')
    ax3.set_xlabel('Distance from IP to Entry (m)')
    ax3.set_ylabel('Muon Rays')
    ax3.set_title('Entry Distance')
    ax3.legend(fontsize=8)
    
    # --- 4: Pair energy ---
    ax4 = fig.add_subplot(gs[1, 0])
    if trident_data:
        Ep = [t['E_pair'] for t in trident_data]
        ax4.hist(Ep, bins=30, alpha=0.6, color='blue', label='All tridents')
        if passing:
            Epp = [t['E_pair'] for t in passing]
            ax4.hist(Epp, bins=30, alpha=0.7, color='red', label='Pass cuts')
        ax4.axvline(2*E_CUT, color='k', ls='--', label=f'2×E_cut')
    ax4.set_xlabel('Pair Energy (GeV)'); ax4.set_ylabel('Events')
    ax4.set_title('e⁺e⁻ Pair Energy'); ax4.legend(fontsize=8)
    ax4.set_yscale('log')
    
    # --- 5: Separation ---
    ax5 = fig.add_subplot(gs[1, 1])
    if trident_data:
        sep = [t['separation_mm'] / 1000 for t in trident_data]  # m
        smuon = min(SEP_MAX * 1.5, max(sep) * 1.1) if sep else SEP_MAX * 1.5
        ax5.hist(sep, bins=40, range=(0, smuon), alpha=0.6, color='blue',
                 label='All tridents')
        if passing:
            sep_p = [t['separation_mm'] / 1000 for t in passing]
            ax5.hist(sep_p, bins=40, range=(0, smuon), alpha=0.7, color='red',
                     label='Pass cuts')
        ax5.axvline(SEP_MIN, color='k', ls='--', lw=1.5, label=f'sep_min')
        ax5.axvline(SEP_MAX, color='darkred', ls='--', lw=1.5, label=f'sep_max')
    ax5.set_xlabel('Separation (m)'); ax5.set_ylabel('Events')
    ax5.set_title('e⁺e⁻ Separation at Wall'); ax5.legend(fontsize=8)
    
    # --- 6: Muon-electron separation ---
    ax6 = fig.add_subplot(gs[1, 2])
    if trident_data:
        mu_close = [min(t['sep_mu_eplus_mm'], t['sep_mu_eminus_mm']) for t in trident_data]
        ax6.hist(mu_close, bins=30, alpha=0.6, color='blue', label='All tridents')
        if passing:
            mu_close_p = [min(t['sep_mu_eplus_mm'], t['sep_mu_eminus_mm']) for t in passing]
            ax6.hist(mu_close_p, bins=30, alpha=0.7, color='red', label='Pass cuts')
    ax6.set_xlabel('μ to Closer Electron (mm)'); ax6.set_ylabel('Events')
    ax6.set_title('Parent Muon – Electron Separation'); ax6.legend(fontsize=8)
    
    # --- 7: Cut flow ---
    ax7 = fig.add_subplot(gs[2, 0])
    if trident_data:
        n_t = len(trident_data)
        n_E = sum(1 for t in trident_data if t['passes_energy'])
        n_svert = sum(1 for t in trident_data 
                     if t['passes_energy'] and t['passes_sep_min'] and t['passes_sep_max'])
        n_smuon = sum(1 for t in trident_data
                     if t['passes_energy'] and t['passes_sep_min'] and t['passes_sep_max'] and t['passes_muon_sep'])
        labels = ['All\ntri.', f'E>{E_CUT*1e3:.0f}\nMeV',
                  f'sep>\n{SEP_MIN*1e3:.0f}mm', f'sep<\n{SEP_MAX:.2f}m']
        counts = [n_t, n_E, n_svert, n_smuon]
        cols = ['#2196F3', '#FF9800', '#4CAF50', '#F44336']
        bars = ax7.bar(labels, counts, color=cols, alpha=0.8)
        for b, c in zip(bars, counts):
            ax7.text(b.get_x() + b.get_width()/2, c + 0.3, str(c),
                     ha='center', fontsize=10, fontweight='bold')
        ax7.set_ylabel('Events (unweighted)')
        ax7.set_title('Cut Flow')
    
    # --- 8: Weighted cut flow ---
    ax8 = fig.add_subplot(gs[2, 1])
    if trident_data:
        w_all = sum(t['weight'] for t in trident_data)
        w_E = sum(t['weight'] for t in trident_data if t['passes_energy'])
        w_svert = sum(t['weight'] for t in trident_data
                     if t['passes_energy'] and t['passes_sep_min'] and t['passes_sep_max'])
        w_smuon = sum(t['weight'] for t in trident_data
                     if t['passes_energy'] and t['passes_sep_min'] and t['passes_sep_max'] and t['passes_muon_sep'])
        labels = ['All\ntri.', f'E>{E_CUT*1e3:.0f}\nMeV',
                  f'sep>\n{SEP_MIN*1e3:.0f}mm', f'sep<\n{SEP_MAX:.2f}m']
        weights = [w_all, w_E, w_svert, w_smuon]
        bars = ax8.bar(labels, weights, color=cols, alpha=0.8)
        for b, w in zip(bars, weights):
            ax8.text(b.get_x() + b.get_width()/2, w * 1.05, f'{w:.1e}',
                     ha='center', fontsize=8)
        ax8.set_ylabel('Weighted rate (per fb⁻¹)')
        ax8.set_title('Weighted Cut Flow')
        ax8.set_yscale('log')
    
    # --- 9: Summary text ---
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')
    
    mean_path = np.mean([h['path_length'] for h in hit_data]) if hit_data else 0
    mean_entry = np.mean([h['entry_d'] for h in hit_data]) if hit_data else 0
    mean_E = np.mean([h['E_mu'] for h in hit_data]) if hit_data else 0
    
    text = (
        f"RESULTS SUMMARY\n"
        f"{'─'*35}\n"
        f"Rays sampled:      {results['n_rays']:>10,}\n"
        f"Hit tunnel:        {results['n_hit']:>10,}\n"
        f"Mean path length:  {mean_path:>10.2f} m\n"
        f"Mean entry dist:   {mean_entry:>10.1f} m\n"
        f"Mean muon energy:  {mean_E:>10.1f} GeV\n"
        f"{'─'*35}\n"
        f"Trident events:    {results['n_trident']:>10}\n"
        f"Pass all cuts:     {results['n_pass']:>10}\n"
        f"{'─'*35}\n"
        f"RATES PER fb⁻¹:\n"
        f"  Muons:       {results['weighted_muons_per_fb']:>12.2e}\n"
        f"  Trident:     {results['weighted_trident_per_fb']:>12.2e}\n"
        f"  Pass cuts:   {results['weighted_pass_per_fb']:>12.2e}\n"
        f"{'─'*35}\n"
        f"PER 3000 fb⁻¹ (HLLHC):\n"
        f"  Before veto: {results['weighted_pass_per_fb']*3000:>12.2e}\n"
        f"  After 99.9%: {results['weighted_pass_per_fb']*3000*1e-3:>12.2e}\n"
    )
    ax9.text(0.05, 0.95, text, transform=ax9.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.savefig(save_path or 'background_fullgeo_mc.png', dpi=150, bbox_inches='tight')
    print(f"  Plot saved to {save_path or 'background_fullgeo_mc.png'}")
    return fig


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    
    print("Building tunnel fiducial mesh...")
    mesh, path_3d = build_tunnel_mesh()
    print(f"  Fiducial volume: {mesh.volume:.1f} m^3")
    
    origin = np.array([0, 0, 0])
    
    mc = TunnelBackgroundMC(mesh, path_3d, origin, seed=123)
    
    # Run with oversampling to get more trident statistics.
    # oversample_trident=100 means we attempt 100 trident trials per muon
    # that hits the tunnel, each with weight/100.
    # This doesn't change the expected rate, just reduces variance.
    results = mc.run(n_rays=500_000)
    
    print(f"\n{'='*60}")
    print(f"  TRIDENT-IN-AIR BACKGROUND (full geometry)")
    print(f"{'='*60}")
    print(f"  Before cuts:  {results['weighted_trident_per_fb']:.2e} per fb^-1")
    print(f"  After cuts:   {results['weighted_pass_per_fb']:.2e} per fb^-1")
    print(f"  Per 3000 fb^-1 (Run 3, before muon veto): "
          f"{results['weighted_pass_per_fb'] * 3000:.2e}")
    print(f"  Per 3000 fb^-1 (after 99.9% muon veto):   "
          f"{results['weighted_pass_per_fb'] * 3000 * 1e-3:.2e}")
    
    print("\nGenerating plots...")
    fig = plot_results(results, save_path='background_fullgeo_mc.png')
    
    print("\nDone!")
