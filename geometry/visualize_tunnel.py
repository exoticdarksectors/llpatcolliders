"""
Visualize the GARGOYLE tunnel geometry.

Produces a 4-panel static figure and an interactive 3D matplotlib window.
CMS convention: X = horizontal transverse, Y = vertical (up), Z = beam axis.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from gargoyle_geometry import (
    mesh_fiducial, path_3d_fiducial,
    tunnel_profile_points, create_profile_mesh,
    Y_POSITION, DETECTOR_THICKNESS,
    TUNNEL_ALPHA, TUNNEL_BETA, TUNNEL_GAMMA, TUNNEL_DELTA, TUNNEL_WALL_HEIGHT,
)

# ── helpers ──────────────────────────────────────────────────────────

def _profile_ring_3d(path_3d, idx, profile_2d):
    """Project a 2D profile into 3D at a given centreline index."""
    if idx == 0:
        tangent = path_3d[1] - path_3d[0]
    elif idx == len(path_3d) - 1:
        tangent = path_3d[idx] - path_3d[idx - 1]
    else:
        tangent = path_3d[idx + 1] - path_3d[idx - 1]
    tangent = tangent / np.linalg.norm(tangent)

    if abs(tangent[1]) < 0.9:
        world_up = np.array([0.0, 1.0, 0.0])
    else:
        world_up = np.array([0.0, 0.0, 1.0])

    right = np.cross(tangent, world_up)
    right /= np.linalg.norm(right)
    up = np.cross(right, tangent)
    up /= np.linalg.norm(up)

    pts = np.array([
        path_3d[idx] + p[0] * right + p[1] * up for p in profile_2d
    ])
    # close the loop
    return np.vstack([pts, pts[0]])


# ── build mesh vertices for the interactive plot ─────────────────────

profile_outer = tunnel_profile_points(inset=0)
profile_fiducial = tunnel_profile_points(inset=DETECTOR_THICKNESS)
n_profile = len(profile_fiducial)

verts_fid, _ = create_profile_mesh(path_3d_fiducial, profile_fiducial)

origin = np.array([0.0, 0.0, 0.0])
path = path_3d_fiducial

# ── Panel figure ─────────────────────────────────────────────────────

fig = plt.figure(figsize=(16, 12))

# ── Panel 1: 3D mesh view ───────────────────────────────────────────

ax1 = fig.add_subplot(221, projection='3d')

# Subsample mesh faces and draw wireframe
# NOTE: swap Y↔Z for display so matplotlib's vertical axis = CMS Y (up)
n_sample = min(1000, len(mesh_fiducial.faces))
face_idx = np.random.choice(len(mesh_fiducial.faces), n_sample, replace=False)
for fi in face_idx:
    tri = mesh_fiducial.vertices[mesh_fiducial.faces[fi]]
    tri = np.vstack([tri, tri[0]])
    ax1.plot(tri[:, 0], tri[:, 2], tri[:, 1], 'b-', alpha=0.1, linewidth=0.5)

# Centreline and origin
ax1.plot(path[:, 0], path[:, 2], path[:, 1], 'r-', linewidth=3, label='Centreline')
ax1.scatter(origin[0], origin[2], origin[1], color='green', s=200, marker='o', label='Origin (IP)')

# Coordinate axes
arrow = 10
ax1.quiver(origin[0], origin[2], origin[1], arrow, 0, 0, color='red',   arrow_length_ratio=0.1)
ax1.quiver(origin[0], origin[2], origin[1], 0, 0, arrow, color='green', arrow_length_ratio=0.1)
ax1.quiver(origin[0], origin[2], origin[1], 0, arrow, 0, color='blue',  arrow_length_ratio=0.1)
ax1.text(origin[0]+arrow, origin[2], origin[1], 'X', fontsize=10)
ax1.text(origin[0], origin[2], origin[1]+arrow, 'Y (up)', fontsize=10)
ax1.text(origin[0], origin[2]+arrow, origin[1], 'Z (beam)', fontsize=10)

ax1.set_xlabel('X (m)')
ax1.set_ylabel('Z (beam, m)')
ax1.set_zlabel('Y (up, m)')
ax1.set_title('3D Tunnel Mesh')
ax1.legend(fontsize=8)

# ── Panel 2: Top view (X-Z plane, bird's eye) ───────────────────────

ax2 = fig.add_subplot(222)

# Draw projected profile outlines at intervals along centreline
for i in np.linspace(0, len(path)-1, 15, dtype=int):
    ring = _profile_ring_3d(path, i, profile_outer)
    ax2.plot(ring[:, 0], ring[:, 2], 'b-', alpha=0.3, linewidth=0.5)

# Centreline
ax2.plot(path[:, 0], path[:, 2], 'r-', linewidth=2, label='Centreline')
ax2.scatter(origin[0], origin[2], color='green', s=200, marker='o', label='Origin')

# Nearest / farthest annotations
dists_xz = np.sqrt((path[:, 0] - origin[0])**2 + (path[:, 2] - origin[2])**2)
near_i = np.argmin(dists_xz)
far_i  = np.argmax(dists_xz)
ax2.plot([origin[0], path[near_i, 0]], [origin[2], path[near_i, 2]],
         'g--', alpha=0.5, label=f'Nearest: {dists_xz[near_i]:.1f} m')
ax2.plot([origin[0], path[far_i, 0]], [origin[2], path[far_i, 2]],
         'r--', alpha=0.5, label=f'Farthest: {dists_xz[far_i]:.1f} m')

ax2.set_xlabel('X (m)')
ax2.set_ylabel('Z (beam, m)')
ax2.set_title('Top View (X-Z Plane)')
ax2.grid(True, alpha=0.3)
ax2.axis('equal')
ax2.legend(fontsize=8)

# ── Panel 3: Side view (X-Y plane, shows height) ────────────────────

ax3 = fig.add_subplot(223)

# Envelope: centreline Y ± half-height of outer profile
profile_y_vals = profile_outer[:, 1]
y_lo = profile_y_vals.min()
y_hi = profile_y_vals.max()

ax3.fill_between(path[:, 0],
                 path[:, 1] + y_lo,
                 path[:, 1] + y_hi,
                 alpha=0.3, color='blue', label='Tunnel envelope')
ax3.plot(path[:, 0], path[:, 1], 'r-', linewidth=2, label='Centreline')
ax3.scatter(origin[0], origin[1], color='green', s=200, marker='o', label='Origin')
ax3.axhline(y=origin[1], color='green', linestyle=':', alpha=0.5)
ax3.axhline(y=Y_POSITION, color='red', linestyle=':', alpha=0.5,
            label=f'Y = {Y_POSITION} m')
ax3.annotate(f'{Y_POSITION} m', xy=(path[0, 0], Y_POSITION),
             xytext=(path[0, 0] - 5, Y_POSITION + 1), fontsize=9)

ax3.set_xlabel('X (m)')
ax3.set_ylabel('Y (up, m)')
ax3.set_title('Side View (X-Y Plane)')
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=8)
ax3.set_ylim(-5, Y_POSITION + 10)

# ── Panel 4: Cross-section view ─────────────────────────────────────

ax4 = fig.add_subplot(224)

# Outer tunnel wall and fiducial inset
outer = tunnel_profile_points(inset=0)
inner = tunnel_profile_points(inset=DETECTOR_THICKNESS)
outer_closed = np.vstack([outer, outer[0]])
inner_closed = np.vstack([inner, inner[0]])

ax4.plot(outer_closed[:, 0], outer_closed[:, 1], 'b-', linewidth=2,
         label='Tunnel wall')
ax4.plot(inner_closed[:, 0], inner_closed[:, 1], 'r--', linewidth=2,
         label=f'Fiducial (inset {DETECTOR_THICKNESS*100:.0f} cm)')
ax4.fill(outer_closed[:, 0], outer_closed[:, 1], alpha=0.08, color='blue')
ax4.fill(inner_closed[:, 0], inner_closed[:, 1], alpha=0.12, color='red')

# Surface labels
ax4.text(0, outer[:, 1].min() - 0.25, 'Floor', ha='center', fontsize=9)
ax4.text(outer[:, 0].max() + 0.15, 0, 'Right\nWall', ha='left',
         va='center', fontsize=9)
ax4.text(outer[:, 0].min() - 0.15, 0, 'Left\nWall', ha='right',
         va='center', fontsize=9)
ax4.text(0, outer[:, 1].max() + 0.15, 'Arch / Ceiling', ha='center',
         fontsize=9)

# Dimension annotations
half_w = TUNNEL_ALPHA / 2
bot_y = outer[:, 1].min()
top_y = outer[:, 1].max()
ax4.annotate('', xy=(half_w, bot_y - 0.35), xytext=(-half_w, bot_y - 0.35),
             arrowprops=dict(arrowstyle='<->', color='gray'))
ax4.text(0, bot_y - 0.55, f'{TUNNEL_ALPHA:.2f} m', ha='center', fontsize=8,
         color='gray')
ax4.annotate('', xy=(half_w + 0.35, top_y), xytext=(half_w + 0.35, bot_y),
             arrowprops=dict(arrowstyle='<->', color='gray'))
ax4.text(half_w + 0.55, (top_y + bot_y) / 2, f'{TUNNEL_BETA:.2f} m',
         ha='left', fontsize=8, color='gray', rotation=90, va='center')

ax4.set_xlabel('Transverse X (m)')
ax4.set_ylabel('Vertical Y (m)')
ax4.set_title('Cross-Section (Arch-on-Wall)')
ax4.grid(True, alpha=0.3)
ax4.axis('equal')
ax4.legend(fontsize=8)

fig.suptitle('GARGOYLE Tunnel Geometry', fontsize=14, y=1.01)
plt.tight_layout()

# ── Print statistics ─────────────────────────────────────────────────

centreline_len = np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))
print("\nTunnel statistics:")
print(f"  Centreline length: {centreline_len:.1f} m")
print(f"  Fiducial volume:   {mesh_fiducial.volume:.1f} m^3")
print(f"  Surface area:      {mesh_fiducial.area:.1f} m^2")
print(f"  Y offset (height): {Y_POSITION} m")
print(f"  Cross-section:     {TUNNEL_ALPHA:.2f} x {TUNNEL_BETA:.2f} m "
      f"(wall {TUNNEL_WALL_HEIGHT:.2f} m + arch {TUNNEL_DELTA:.2f} m)")

# ── Interactive 3D figure ────────────────────────────────────────────

fig2 = plt.figure(figsize=(12, 10))
ax = fig2.add_subplot(111, projection='3d')

# Draw profile rings at regular intervals
# NOTE: swap Y↔Z for display so matplotlib's vertical axis = CMS Y (up)
n_rings = 30
ring_indices = np.linspace(0, len(path)-1, n_rings, dtype=int)
for i in ring_indices:
    ring = _profile_ring_3d(path, i, profile_fiducial)
    ax.plot(ring[:, 0], ring[:, 2], ring[:, 1], 'b-', alpha=0.5, linewidth=1)

# Longitudinal lines connecting corresponding profile vertices
for j in range(0, n_profile, max(1, n_profile // 8)):
    line = []
    for i in range(len(path)):
        idx = i * n_profile + j
        if idx < len(verts_fid):
            line.append(verts_fid[idx])
    if line:
        line = np.array(line)
        ax.plot(line[:, 0], line[:, 2], line[:, 1], 'b-', alpha=0.3,
                linewidth=0.5)

# Centreline
ax.plot(path[:, 0], path[:, 2], path[:, 1], 'r-', linewidth=3,
        label='Centreline')

# Origin and axes
ax.scatter(origin[0], origin[2], origin[1], color='green', s=300, marker='o',
           edgecolors='black', linewidth=2, label='Origin (IP)')

axis_len = 15
ax.quiver(origin[0], origin[2], origin[1], axis_len, 0, 0, color='red',   arrow_length_ratio=0.1, linewidth=2)
ax.quiver(origin[0], origin[2], origin[1], 0, 0, axis_len, color='green', arrow_length_ratio=0.1, linewidth=2)
ax.quiver(origin[0], origin[2], origin[1], 0, axis_len, 0, color='blue',  arrow_length_ratio=0.1, linewidth=2)
ax.text(origin[0]+axis_len, origin[2], origin[1], 'X', fontsize=12)
ax.text(origin[0], origin[2], origin[1]+axis_len, 'Y (up)', fontsize=12)
ax.text(origin[0], origin[2]+axis_len, origin[1], 'Z (beam)', fontsize=12)

ax.set_xlabel('X (m)', fontsize=10)
ax.set_ylabel('Z (beam, m)', fontsize=10)
ax.set_zlabel('Y (up, m)', fontsize=10)
ax.set_title('GARGOYLE Tunnel — Interactive 3D (rotate with mouse)', fontsize=14)
ax.view_init(elev=20, azim=45)
ax.set_box_aspect([1, 1, 0.5])
ax.legend(fontsize=10)

plt.show()
