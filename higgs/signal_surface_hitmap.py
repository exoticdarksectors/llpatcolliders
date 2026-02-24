"""
Signal Surface Hit Map Analysis (v2)
=====================================

Determines which tunnel surfaces receive LLP decay products,
allowing optimization of tracking coverage for a target efficiency.

Changes from v1:
  - Correct surface labels derived from actual profile geometry
  - Cumulative efficiency uses contiguous angular regions grown from
    the signal peak, not individual sorted bins
  - Coarser binning option for low-statistics samples

Uses the SAME geometry framework and CSV input format as
decayProbPerEvent_2body.py.

Usage:
  python signal_surface_hitmap_v2.py LLP.csv              # path-length weighting
  python signal_surface_hitmap_v2.py LLP.csv 1e-7         # decay prob at τ = 100 ns
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm

from gargoyle_geometry import (
    TUNNEL_ALPHA, TUNNEL_GAMMA, TUNNEL_DELTA,
    TUNNEL_WALL_HEIGHT, DETECTOR_THICKNESS,
    tunnel_profile_points, eta_phi_to_direction,
    mesh_fiducial, path_3d_fiducial,
)

E_CUT = 0.600       # GeV
SEP_MIN = 0.001     # m


# =============================================================================
# Acceptance (same as decayProbPerEvent_2body.py)
# =============================================================================
def compute_acceptance_limits(gamma, beta, mass, E_cut=E_CUT):
    c_E = (1.0 - 2.0 * E_cut / (gamma * mass)) / beta
    c_E = min(c_E, 1.0)
    return min(c_E, beta)

def compute_c_S(d_remaining, gamma, beta, sep_min=SEP_MIN):
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

def acceptance_at_exit(gamma, beta, mass, d_remaining=0.01,
                       E_cut=E_CUT, sep_min=SEP_MIN):
    c_max = compute_acceptance_limits(gamma, beta, mass, E_cut)
    if c_max <= 0:
        return 0.0
    c_S = compute_c_S(d_remaining, gamma, beta, sep_min)
    return max(0.0, c_max - c_S)


# =============================================================================
# Tunnel path derived quantities
# =============================================================================
path_3d = path_3d_fiducial

seg_lengths = np.array([np.linalg.norm(path_3d[i+1] - path_3d[i])
                         for i in range(len(path_3d)-1)])
cumulative_length = np.concatenate([[0], np.cumsum(seg_lengths)])
total_length = cumulative_length[-1]


# =============================================================================
# Profile geometry for surface labels
# =============================================================================
profile_pts = tunnel_profile_points(inset=0.0)
n_prof = len(profile_pts)

# Compute the angle θ in local frame for each profile vertex
profile_angles = np.array([np.arctan2(p[1], p[0]) for p in profile_pts])
profile_angles = np.where(profile_angles < 0,
                           profile_angles + 2*np.pi, profile_angles)

# Identify key profile vertices and their angles
# Profile order: [bottom-left, bottom-right, right wall up..., arch..., left wall down...]
n_wall = 4
n_arch_pts = 31  # n_arch-1 interior points for n_arch=32

idx_floor_left = 0
idx_floor_right = 1
idx_rwall_top = 1 + n_wall
idx_arch_top = 1 + n_wall + (n_arch_pts // 2)
idx_lwall_top = 1 + n_wall + n_arch_pts
idx_lwall_bot = n_prof - 1

theta_floor_right  = profile_angles[idx_floor_right]
theta_rwall_top    = profile_angles[idx_rwall_top]
theta_arch_top     = profile_angles[idx_arch_top]
theta_lwall_top    = profile_angles[idx_lwall_top]
theta_floor_left   = profile_angles[idx_floor_left]

print(f"\nProfile boundary angles (degrees):")
print(f"  Floor right corner:    {np.degrees(theta_floor_right):.1f}°")
print(f"  Right wall top (springline): {np.degrees(theta_rwall_top):.1f}°")
print(f"  Arch apex:             {np.degrees(theta_arch_top):.1f}°")
print(f"  Left wall top (springline):  {np.degrees(theta_lwall_top):.1f}°")
print(f"  Floor left corner:     {np.degrees(theta_floor_left):.1f}°")

# Surface boundaries
surface_boundaries = [
    (theta_floor_right, theta_rwall_top, 'Right Wall'),
    (theta_rwall_top,   theta_lwall_top, 'Arch/Ceiling'),
    (theta_lwall_top,   theta_floor_left, 'Left Wall'),
]
theta_floor_center = (theta_floor_left + theta_floor_right) / 2
surface_boundaries.append(
    (theta_floor_left, theta_floor_right, 'Floor'))

# Perimeter
perimeter = sum(np.linalg.norm(profile_pts[(i+1) % n_prof] - profile_pts[i])
                for i in range(n_prof))
total_surface_area = perimeter * total_length
print(f"  Tunnel perimeter: {perimeter:.2f} m")
print(f"  Tunnel length: {total_length:.0f} m")
print(f"  Total surface area: {total_surface_area:.0f} m²")


# =============================================================================
# Surface classification using geometry-derived boundaries
# =============================================================================
def classify_by_theta(theta):
    """Classify a profile angle to a surface name using geometry boundaries."""
    t = theta % (2 * np.pi)

    # Right wall: from floor_right up to springline
    if theta_floor_right <= theta_rwall_top:
        if theta_floor_right <= t < theta_rwall_top:
            return 'Right Wall'
    else:
        if t >= theta_floor_right or t < theta_rwall_top:
            return 'Right Wall'

    # Arch/ceiling
    if theta_rwall_top <= theta_lwall_top:
        if theta_rwall_top <= t < theta_lwall_top:
            return 'Arch/Ceiling'
    else:
        if t >= theta_rwall_top or t < theta_lwall_top:
            return 'Arch/Ceiling'

    # Left wall
    if theta_lwall_top <= theta_floor_left:
        if theta_lwall_top <= t < theta_floor_left:
            return 'Left Wall'
    else:
        if t >= theta_lwall_top or t < theta_floor_left:
            return 'Left Wall'

    return 'Floor'


def classify_exit_point(point, path_3d, cumulative_length):
    """
    For a 3D point on the tunnel wall, compute (s, theta, x_local, y_local).
    """
    best_s = 0.0
    best_dist_sq = np.inf
    best_x = 0.0
    best_y = 0.0

    for i in range(len(path_3d) - 1):
        seg = path_3d[i+1] - path_3d[i]
        seg_len = np.linalg.norm(seg)
        if seg_len == 0:
            continue
        seg_hat = seg / seg_len
        t = np.clip(np.dot(point - path_3d[i], seg_hat), 0, seg_len)
        closest = path_3d[i] + t * seg_hat
        diff = point - closest
        dist_sq = np.dot(diff, diff)

        if dist_sq < best_dist_sq:
            best_dist_sq = dist_sq
            best_s = cumulative_length[i] + t
            tangent = seg_hat
            # CMS convention: Y is up
            if abs(tangent[1]) < 0.9:
                world_up = np.array([0., 1., 0.])
            else:
                world_up = np.array([0., 0., 1.])
            right = np.cross(tangent, world_up)
            right /= np.linalg.norm(right)
            up = np.cross(right, tangent)
            up /= np.linalg.norm(up)
            best_x = np.dot(diff, right)
            best_y = np.dot(diff, up)

    theta = np.arctan2(best_y, best_x)
    if theta < 0:
        theta += 2 * np.pi

    return best_s, theta, best_x, best_y


# =============================================================================
# Contiguous angular region efficiency
# =============================================================================
def contiguous_efficiency(theta_values, weights, n_sectors=72):
    """
    Compute signal efficiency as a function of contiguous angular coverage.

    Divides the profile into n_sectors angular bins. For each possible
    contiguous arc (starting point × width), computes the enclosed signal
    fraction. Returns the optimal arc for each target width.

    Returns
    -------
    arc_widths : array of angular widths (radians), length n_sectors
    best_eff   : best efficiency achievable for each width
    best_start : optimal starting angle for each width
    sector_signal : signal weight per sector
    sector_edges  : angular bin edges
    """
    sector_edges = np.linspace(0, 2*np.pi, n_sectors + 1)
    sector_signal, _ = np.histogram(theta_values, bins=sector_edges, weights=weights)

    best_eff = np.zeros(n_sectors)
    best_start_idx = np.zeros(n_sectors, dtype=int)

    doubled = np.concatenate([sector_signal, sector_signal])

    for width in range(1, n_sectors + 1):
        cumsum = np.cumsum(doubled)
        cumsum = np.insert(cumsum, 0, 0)
        window_sums = cumsum[width:width+n_sectors] - cumsum[:n_sectors]
        best_idx = np.argmax(window_sums)
        best_eff[width - 1] = window_sums[best_idx]
        best_start_idx[width - 1] = best_idx

    arc_widths = np.arange(1, n_sectors + 1) * (2*np.pi / n_sectors)
    arc_fractions = arc_widths / (2*np.pi)
    arc_areas = arc_fractions * perimeter * total_length

    return arc_widths, arc_areas, best_eff, best_start_idx, sector_signal, sector_edges


# =============================================================================
# Main analysis
# =============================================================================
def run_analysis(csv_file, lifetime_seconds=None):
    origin = np.array([0.0, 0.0, 0.0])

    # --- Load CSV ---
    df = pd.read_csv(csv_file)
    n = len(df)
    print(f"\nLoaded {n} particles from {csv_file}")
    print(f"  Events: {df['event'].nunique()}")
    print(f"  Mass: {df['mass'].iloc[0]:.4f} GeV")
    print(f"  Momentum range: {df['momentum'].min():.1f} – "
          f"{df['momentum'].max():.1f} GeV")

    # --- Ray-cast ---
    print("\nRay-casting LLPs against fiducial volume...")
    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)
    entry_pts = np.full((n, 3), np.nan)
    exit_pts = np.full((n, 3), np.nan)

    for idx, row in tqdm(df.iterrows(), total=n, desc="Ray-casting"):
        direction = eta_phi_to_direction(row['eta'], row['phi'])
        locations, _, _ = mesh_fiducial.ray.intersects_location(
            ray_origins=[origin], ray_directions=[direction])
        if len(locations) >= 2:
            dists = sorted([(np.linalg.norm(loc - origin), loc) for loc in locations])
            hits[idx] = True
            entry_d[idx] = dists[0][0]
            exit_d[idx] = dists[-1][0]
            entry_pts[idx] = dists[0][1]
            exit_pts[idx] = dists[-1][1]

    n_hits = hits.sum()
    print(f"  Hits: {n_hits} / {n} ({n_hits/n*100:.1f}%)")
    if n_hits == 0:
        print("No particles hit fiducial volume. Exiting.")
        return

    path_l = exit_d[hits] - entry_d[hits]
    print(f"  Mean path length: {path_l.mean():.2f} m")

    # --- Weights ---
    momentum = df['momentum'].values
    mass = df['mass'].values
    energy = np.sqrt(momentum**2 + mass**2)
    gamma = energy / mass
    beta = momentum / energy

    if lifetime_seconds is not None:
        from scipy.integrate import quad
        print(f"\nComputing decay probabilities (τ = {lifetime_seconds:.2e} s)...")
        decay_prob = np.zeros(n)
        for idx in tqdm(np.where(hits)[0], desc="Decay prob"):
            g = gamma[idx]; b = beta[idx]; m = mass[idx]
            lam = g * b * 299792458.0 * lifetime_seconds
            d_entry = entry_d[idx]; d_exit = exit_d[idx]
            def integrand(d):
                return (1.0/lam) * np.exp(-d/lam) * \
                       acceptance_at_exit(g, b, m, d_exit - d)
            result, _ = quad(integrand, d_entry, d_exit, limit=100)
            decay_prob[idx] = result
        weights = decay_prob[hits]
        weight_label = f'Decay prob (τ={lifetime_seconds:.1e}s)'
    else:
        weights = path_l.copy()
        weight_label = 'Path length weight'

    if weights.sum() == 0:
        print("All weights are zero. Exiting.")
        return
    weights_norm = weights / weights.sum()

    # --- Classify exit & entry points ---
    print("\nClassifying surface points...")
    hit_indices = np.where(hits)[0]

    exit_s     = np.zeros(n_hits)
    exit_theta = np.zeros(n_hits)
    exit_xl    = np.zeros(n_hits)
    exit_yl    = np.zeros(n_hits)

    entry_s     = np.zeros(n_hits)
    entry_theta = np.zeros(n_hits)
    entry_xl    = np.zeros(n_hits)
    entry_yl    = np.zeros(n_hits)

    for i, idx in enumerate(tqdm(hit_indices, desc="Exit points")):
        s, th, xl, yl = classify_exit_point(exit_pts[idx], path_3d, cumulative_length)
        exit_s[i] = s; exit_theta[i] = th; exit_xl[i] = xl; exit_yl[i] = yl

    for i, idx in enumerate(tqdm(hit_indices, desc="Entry points")):
        s, th, xl, yl = classify_exit_point(entry_pts[idx], path_3d, cumulative_length)
        entry_s[i] = s; entry_theta[i] = th; entry_xl[i] = xl; entry_yl[i] = yl

    exit_surface  = np.array([classify_by_theta(t) for t in exit_theta])
    entry_surface = np.array([classify_by_theta(t) for t in entry_theta])

    # --- Surface fractions ---
    surfaces = ['Floor', 'Right Wall', 'Arch/Ceiling', 'Left Wall']
    print(f"\n{'='*60}")
    print("SIGNAL SURFACE DISTRIBUTION")
    print(f"{'='*60}")
    print(f"\n{'Surface':<20} {'Exit count':>10} {'Exit wtd%':>10} {'Entry wtd%':>10}")
    print("-" * 55)
    for surf in surfaces:
        ex_mask = exit_surface == surf
        en_mask = entry_surface == surf
        print(f"{surf:<20} {ex_mask.sum():>10} "
              f"{weights_norm[ex_mask].sum()*100:>9.1f}% "
              f"{weights_norm[en_mask].sum()*100:>9.1f}%")

    # --- Contiguous angular efficiency ---
    print(f"\n{'='*60}")
    print("CONTIGUOUS ANGULAR REGION EFFICIENCY")
    print(f"{'='*60}")

    n_sectors = 72  # 5° per sector
    arc_widths, arc_areas, best_eff, best_start, sector_sig, sector_edges = \
        contiguous_efficiency(exit_theta, weights_norm, n_sectors)

    targets = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]
    area_for_target = {}
    angle_for_target = {}

    print(f"\n{'Eff':>6} {'Arc (°)':>8} {'Area (m²)':>10} {'% total':>8} "
          f"{'4-layer':>10} {'~Cost':>8} {'Arc center':>12}")
    print("-" * 70)
    for tgt in targets:
        i = np.searchsorted(best_eff, tgt)
        if i < len(arc_areas):
            a = arc_areas[i]
            w_deg = np.degrees(arc_widths[i])
            start_deg = np.degrees(sector_edges[best_start[i]])
            center_deg = start_deg + w_deg / 2
            cost = a / 500 * 7
            area_for_target[tgt] = a
            angle_for_target[tgt] = (sector_edges[best_start[i]],
                                      arc_widths[i])
            print(f"{tgt:>5.0%} {w_deg:>7.0f}° {a:>10.0f} {a/total_surface_area:>7.0%} "
                  f"{4*a:>10.0f} {'${:.1f}M'.format(cost):>8} "
                  f"{center_deg:>8.0f}°")

    # Along-tunnel efficiency
    print(f"\n--- Along-tunnel coverage (full perimeter tracking) ---")
    n_z_bins = 50
    z_edges = np.linspace(0, total_length, n_z_bins + 1)
    z_signal, _ = np.histogram(exit_s, bins=z_edges, weights=weights_norm)
    z_cum = np.cumsum(np.sort(z_signal)[::-1])
    doubled_z = np.concatenate([z_signal, z_signal])
    z_best_eff = np.zeros(n_z_bins)
    for width in range(1, n_z_bins + 1):
        cumsum = np.cumsum(np.insert(doubled_z, 0, 0))
        window = cumsum[width:width+n_z_bins] - cumsum[:n_z_bins]
        z_best_eff[width-1] = window.max()
    z_lengths = np.arange(1, n_z_bins+1) * (total_length / n_z_bins)

    print(f"\n{'Eff':>6} {'Length (m)':>10} {'% of tunnel':>12}")
    print("-" * 30)
    for tgt in [0.50, 0.80, 0.90, 0.95]:
        i = np.searchsorted(z_best_eff, tgt)
        if i < len(z_lengths):
            print(f"{tgt:>5.0%} {z_lengths[i]:>10.0f} {z_lengths[i]/total_length:>11.0%}")

    # =================================================================
    # 2D histogram
    # =================================================================
    n_s_bins = 40
    n_theta_bins = 36
    s_edges = np.linspace(0, total_length, n_s_bins + 1)
    theta_edges_2d = np.linspace(0, 2*np.pi, n_theta_bins + 1)

    exit_hist, _, _ = np.histogram2d(
        exit_s, exit_theta, bins=[s_edges, theta_edges_2d], weights=weights_norm)
    entry_hist, _, _ = np.histogram2d(
        entry_s, entry_theta, bins=[s_edges, theta_edges_2d], weights=weights_norm)

    # =================================================================
    # Plotting
    # =================================================================
    print("\nGenerating plots...")

    s_centers = (s_edges[:-1] + s_edges[1:]) / 2
    th_centers = np.degrees((theta_edges_2d[:-1] + theta_edges_2d[1:]) / 2)

    def circular_midpoint_deg(start_rad, end_rad):
        """Midpoint of arc going counterclockwise from start to end."""
        s = np.degrees(start_rad) % 360
        e = np.degrees(end_rad) % 360
        if e >= s:
            return (s + e) / 2
        else:
            span = (360 - s) + e
            mid = s + span / 2
            return mid % 360

    label_positions = {
        'Right Wall':  circular_midpoint_deg(theta_floor_right, theta_rwall_top),
        'Arch/Ceiling':circular_midpoint_deg(theta_rwall_top,   theta_lwall_top),
        'Left Wall':   circular_midpoint_deg(theta_lwall_top,   theta_floor_left),
        'Floor':       circular_midpoint_deg(theta_floor_left,  theta_floor_right),
    }

    boundary_angles_deg = sorted([
        np.degrees(theta_floor_right) % 360,
        np.degrees(theta_rwall_top) % 360,
        np.degrees(theta_lwall_top) % 360,
        np.degrees(theta_floor_left) % 360,
    ])

    print(f"\n  Label positions (degrees):")
    for name, pos in label_positions.items():
        print(f"    {name:20s}: {pos:.1f}°")
    print(f"  Boundary lines at: {', '.join(f'{a:.0f}°' for a in boundary_angles_deg)}")

    fig = plt.figure(figsize=(18, 24))
    gs = GridSpec(5, 2, figure=fig, hspace=0.35, wspace=0.30)

    # ---- Panel 1: Signal exit heat map ----
    ax = fig.add_subplot(gs[0, :])
    im = ax.pcolormesh(s_centers, th_centers, exit_hist.T,
                        cmap='hot', shading='auto')
    plt.colorbar(im, ax=ax, label='Signal weight')
    ax.set_xlabel('Distance along tunnel (m)')
    ax.set_ylabel('Angle around profile (deg)')
    ax.set_title('Signal Hit Distribution (exit points — where e⁺e⁻ hit wall)')
    for angle_deg in boundary_angles_deg:
        ax.axhline(y=angle_deg, color='cyan', alpha=0.6, ls='--', lw=1)
    for name, pos in label_positions.items():
        ax.text(total_length * 0.98, pos, name, color='cyan',
                fontsize=9, fontweight='bold', ha='right', va='center',
                bbox=dict(facecolor='black', alpha=0.6, pad=2))

    # ---- Panel 2: Entry heat map ----
    ax = fig.add_subplot(gs[1, :])
    im2 = ax.pcolormesh(s_centers, th_centers, entry_hist.T,
                         cmap='Blues', shading='auto')
    plt.colorbar(im2, ax=ax, label='Entry weight')
    ax.set_xlabel('Distance along tunnel (m)')
    ax.set_ylabel('Angle around profile (deg)')
    ax.set_title('LLP Entry Distribution (where LLPs enter fiducial volume)')
    for angle_deg in boundary_angles_deg:
        ax.axhline(y=angle_deg, color='gray', alpha=0.5, ls='--', lw=1)
    for name, pos in label_positions.items():
        ax.text(total_length * 0.98, pos, name, color='black',
                fontsize=9, fontweight='bold', ha='right', va='center',
                bbox=dict(facecolor='white', alpha=0.6, pad=2))

    # ---- Panel 3: 80% contiguous region overlay ----
    ax = fig.add_subplot(gs[2, :])
    ax.pcolormesh(s_centers, th_centers, exit_hist.T,
                   cmap='hot', shading='auto', alpha=0.6)
    if 0.80 in angle_for_target:
        start_rad, width_rad = angle_for_target[0.80]
        end_rad = start_rad + width_rad
        start_deg = np.degrees(start_rad)
        end_deg = np.degrees(end_rad)
        a80 = area_for_target[0.80]
        if end_deg <= 360:
            ax.axhspan(start_deg, end_deg, color='green', alpha=0.25)
            ax.axhline(y=start_deg, color='green', lw=2)
            ax.axhline(y=end_deg, color='green', lw=2)
        else:
            ax.axhspan(start_deg, 360, color='green', alpha=0.25)
            ax.axhspan(0, end_deg - 360, color='green', alpha=0.25)
            ax.axhline(y=start_deg, color='green', lw=2)
            ax.axhline(y=end_deg - 360, color='green', lw=2)
        ax.set_title(f'80% Efficiency: contiguous {np.degrees(width_rad):.0f}° arc, '
                     f'{a80:.0f} m² of {total_surface_area:.0f} m²')
    else:
        ax.set_title('80% Efficiency Region')
    ax.set_xlabel('Distance along tunnel (m)')
    ax.set_ylabel('Angle around profile (deg)')
    for angle_deg in boundary_angles_deg:
        ax.axhline(y=angle_deg, color='cyan', alpha=0.4, ls='--', lw=1)
    for name, pos in label_positions.items():
        ax.text(total_length * 0.98, pos, name, color='cyan', fontsize=8,
                fontweight='bold', ha='right', va='center',
                bbox=dict(facecolor='black', alpha=0.5, pad=2))

    # ---- Panel 4: Contiguous efficiency curve ----
    ax = fig.add_subplot(gs[3, 0])
    ax.plot(arc_areas, best_eff * 100, 'b-', lw=2, label='Angular arc')
    for tgt, col in [(0.80, 'r'), (0.90, 'orange'), (0.95, 'green')]:
        ax.axhline(y=tgt*100, color=col, ls='--', alpha=0.5, label=f'{tgt:.0%}')
        if tgt in area_for_target:
            a_t = area_for_target[tgt]
            cost = a_t / 500 * 7
            ax.annotate(f'{a_t:.0f} m²\n~${cost:.1f}M',
                        xy=(a_t, tgt*100), fontsize=8,
                        xytext=(a_t + total_surface_area*0.08, tgt*100 - 4),
                        arrowprops=dict(arrowstyle='->', color=col), color=col)
    ax.set_xlabel('Instrumented area (m²)')
    ax.set_ylabel('Signal efficiency (%)')
    ax.set_title('Efficiency vs Contiguous Tracking Area')
    ax.legend(fontsize=8)
    ax.set_xlim(0, total_surface_area * 1.05)
    ax.set_ylim(0, 102)
    ax.grid(True, alpha=0.3)

    # ---- Panel 5: Signal fraction by surface ----
    ax = fig.add_subplot(gs[3, 1])
    fracs_exit = [weights_norm[exit_surface == s].sum() * 100 for s in surfaces]
    fracs_entry = [weights_norm[entry_surface == s].sum() * 100 for s in surfaces]
    x_pos = np.arange(len(surfaces))
    w = 0.35
    colors = ['#8B4513', '#228B22', '#FFD700', '#4169E1']
    bars1 = ax.bar(x_pos - w/2, fracs_exit, w, label='Signal exit',
                    color=colors, edgecolor='black')
    ax.bar(x_pos + w/2, fracs_entry, w, label='LLP entry',
           color='none', edgecolor='red', linewidth=2, linestyle='--')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(surfaces, fontsize=9)
    ax.set_ylabel('Fraction (%)')
    ax.set_title('Signal vs Entry by Surface')
    ax.legend()
    for bar, frac in zip(bars1, fracs_exit):
        if frac > 0.5:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f'{frac:.1f}%', ha='center', fontsize=9, fontweight='bold')

    # ---- Panel 6: Cross-section density profile ----
    ax = fig.add_subplot(gs[4, 0])

    # Bin densities around the profile angle
    n_band = 120
    band_edges = np.linspace(0, 2*np.pi, n_band + 1)
    band_centers = (band_edges[:-1] + band_edges[1:]) / 2

    exit_density, _ = np.histogram(exit_theta, bins=band_edges,
                                    weights=weights_norm)
    entry_density, _ = np.histogram(entry_theta, bins=band_edges,
                                     weights=weights_norm)

    # Map each angular bin to a point on the fiducial wall + outward normal
    # Use the outer profile (inset=0) as reference wall
    prof_outer = tunnel_profile_points(inset=0.0)
    prof_outer_c = np.vstack([prof_outer, prof_outer[0]])

    # For each band center, find the wall point by interpolating the profile
    # Profile points are ordered CCW; compute their angles
    prof_ang = np.array([np.arctan2(p[1], p[0]) for p in prof_outer])
    prof_ang = np.where(prof_ang < 0, prof_ang + 2*np.pi, prof_ang)

    # Sort by angle for interpolation
    sort_idx = np.argsort(prof_ang)
    prof_ang_sorted = prof_ang[sort_idx]
    prof_xy_sorted = prof_outer[sort_idx]

    # Wrap for interpolation continuity
    prof_ang_wrap = np.concatenate([prof_ang_sorted - 2*np.pi,
                                    prof_ang_sorted,
                                    prof_ang_sorted + 2*np.pi])
    prof_x_wrap = np.tile(prof_xy_sorted[:, 0], 3)
    prof_y_wrap = np.tile(prof_xy_sorted[:, 1], 3)

    wall_x = np.interp(band_centers, prof_ang_wrap, prof_x_wrap)
    wall_y = np.interp(band_centers, prof_ang_wrap, prof_y_wrap)

    # Outward normal at each point (from centroid at ~origin)
    r = np.sqrt(wall_x**2 + wall_y**2)
    nx, ny = wall_x / r, wall_y / r

    # Draw density as colored thick lines offset from wall
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    offset_out = 0.18   # exit: outside the wall
    offset_in  = -0.12  # entry: inside the wall
    lw = 4.5

    # Exit density (outside, red colormap)
    ex_x = wall_x + offset_out * nx
    ex_y = wall_y + offset_out * ny
    segments_ex = [[(ex_x[i], ex_y[i]), (ex_x[(i+1) % n_band], ex_y[(i+1) % n_band])]
                   for i in range(n_band)]
    norm_ex = Normalize(vmin=0, vmax=max(exit_density.max(), 1e-10))
    lc_ex = LineCollection(segments_ex, cmap='Reds', norm=norm_ex,
                           linewidths=lw, zorder=4)
    lc_ex.set_array(exit_density)
    ax.add_collection(lc_ex)
    plt.colorbar(lc_ex, ax=ax, label='Signal exit', shrink=0.7,
                 pad=0.02, aspect=30)

    # Entry density (inside, blue colormap)
    en_x = wall_x + offset_in * nx
    en_y = wall_y + offset_in * ny
    segments_en = [[(en_x[i], en_y[i]), (en_x[(i+1) % n_band], en_y[(i+1) % n_band])]
                   for i in range(n_band)]
    norm_en = Normalize(vmin=0, vmax=max(entry_density.max(), 1e-10))
    lc_en = LineCollection(segments_en, cmap='Blues', norm=norm_en,
                           linewidths=lw, zorder=3)
    lc_en.set_array(entry_density)
    ax.add_collection(lc_en)
    plt.colorbar(lc_en, ax=ax, label='LLP entry', shrink=0.7,
                 pad=0.02, aspect=30)

    # Tunnel outline
    ax.plot(prof_outer_c[:, 0], prof_outer_c[:, 1], 'k-', lw=1.5)
    prof_in = tunnel_profile_points(inset=DETECTOR_THICKNESS)
    prof_in_c = np.vstack([prof_in, prof_in[0]])
    ax.plot(prof_in_c[:, 0], prof_in_c[:, 1], 'k--', lw=0.8, alpha=0.5)

    ax.set_xlabel('Local x (m)')
    ax.set_ylabel('Local y (m)')
    ax.set_title('Cross-section Density:\nSignal exit (outer/red) — LLP entry (inner/blue)')
    ax.set_aspect('equal')
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.8, 2.0)
    ax.grid(True, alpha=0.2)

    # ---- Panel 7: Angular distribution (1D) with surface boundaries ----
    ax = fig.add_subplot(gs[4, 1])
    sector_centers_deg = np.degrees(
        (sector_edges[:-1] + sector_edges[1:]) / 2)
    ax.bar(sector_centers_deg, sector_sig * 100,
           width=360/n_sectors * 0.9, color='red', alpha=0.6,
           label='Signal exit')
    entry_sector, _ = np.histogram(entry_theta, bins=sector_edges,
                                    weights=weights_norm)
    ax.bar(sector_centers_deg, entry_sector * 100,
           width=360/n_sectors * 0.9, color='blue', alpha=0.2,
           label='LLP entry')
    for angle_deg in boundary_angles_deg:
        ax.axvline(x=angle_deg, color='gray', alpha=0.5, ls='--')
    for name, pos in label_positions.items():
        ax.annotate(name, xy=(pos, 1.0), xycoords=('data', 'axes fraction'),
                    fontsize=7, ha='center', va='top', color='gray',
                    rotation=90)
    ax.set_xlabel('Profile angle (deg)')
    ax.set_ylabel('Signal fraction (%)')
    ax.set_title('Angular Distribution with Surface Boundaries')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 360)

    fig.suptitle(f'Signal Surface Coverage — {csv_file}\n'
                 f'M = {df["mass"].iloc[0]:.3f} GeV, '
                 f'{n_hits} hits, {weight_label}',
                 fontsize=13, fontweight='bold', y=1.01)

    outname = csv_file.replace('.csv', '') + '_surface_hitmap_v2.png'
    plt.savefig(outname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {outname}")
    plt.close('all')

    # =================================================================
    # Summary
    # =================================================================
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Input: {csv_file}")
    print(f"LLP mass: {df['mass'].iloc[0]:.4f} GeV")
    print(f"Particles hitting fiducial: {n_hits}/{n} ({n_hits/n*100:.1f}%)")
    print(f"Mean path length: {path_l.mean():.2f} m")
    print(f"\nSignal exit by surface (weighted):")
    for surf in surfaces:
        mask = exit_surface == surf
        print(f"  {surf:20s}: {weights_norm[mask].sum():.1%}")

    for tgt in [0.80, 0.90, 0.95]:
        if tgt in area_for_target:
            a = area_for_target[tgt]
            w_deg = np.degrees(angle_for_target[tgt][1])
            cost = a / 500 * 7
            remaining = total_surface_area - a
            print(f"\n{tgt:.0%} efficiency: {a:.0f} m² tracking "
                  f"(contiguous {w_deg:.0f}° arc, "
                  f"{a/total_surface_area:.0%} of total, ~${cost:.1f}M)")
            print(f"  Remaining {remaining:.0f} m²: veto-only (~$0.5-1M)")

    print("\nDone!")


# =============================================================================
if __name__ == "__main__":
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "LLPSmall.csv"
    lifetime = float(sys.argv[2]) if len(sys.argv) > 2 else None
    run_analysis(csv_file, lifetime_seconds=lifetime)
