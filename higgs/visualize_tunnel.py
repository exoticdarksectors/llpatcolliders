"""
Visualize the GRENDEL tunnel geometry.

Produces a 4-panel static figure and an interactive 3D matplotlib window.
CMS convention: X = horizontal transverse, Y = vertical (up), Z = beam axis.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from grendel_geometry import (
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


def set_equal_aspect_3d(ax, expand=(1.0, 1.0, 1.0)):
    """Make 1 m the same length on all three 3D axes (no stretching).

    Sets the box aspect proportional to the autoscaled data extents, so the
    tunnel keeps its true proportions instead of being forced into a cube.
    `expand` grows each axis range about its centre before fixing the aspect
    (e.g. 1.2 on the CMS-Z/beam axis = +20% range). Axis order is the matplotlib
    (x, y, z) display order: here x=X, y=Z (beam), z=Y (up)."""
    setters = (ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d)
    lims = (ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d())
    ranges = []
    for (lo, hi), f, setter in zip(lims, expand, setters):
        c, h = 0.5 * (lo + hi), 0.5 * (hi - lo) * f
        setter(c - h, c + h)
        ranges.append(2 * h)
    ax.set_box_aspect(ranges)


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
ax1.set_title('3D Active Decay Volume (inner tracker boundary)')
ax1.legend(fontsize=8)
set_equal_aspect_3d(ax1, expand=(1.0, 1.2, 1.0))  # +20% on CMS-Z (beam)

# ── Panel 2: Top view (X-Z plane, bird's eye) ───────────────────────

ax2 = fig.add_subplot(222)

# Tunnel width as a solid blue band: centreline offset by ±half-width
# perpendicular to the path in the X-Z plane.
half_width = profile_outer[:, 0].max()
xz = path[:, [0, 2]]
tang = np.gradient(xz, axis=0)
tang /= np.linalg.norm(tang, axis=1, keepdims=True)
perp = np.column_stack([-tang[:, 1], tang[:, 0]])  # in-plane normal
left = xz + half_width * perp
right = xz - half_width * perp
band = np.vstack([left, right[::-1]])
ax2.fill(band[:, 0], band[:, 1], color='blue', alpha=0.3, label='Tunnel width')

# Centreline
ax2.plot(path[:, 0], path[:, 2], 'r-', linewidth=2, label='Centreline')
ax2.scatter(origin[0], origin[2], color='green', s=200, marker='o', label='Origin')

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
#
# Geometry: two tracker layers (separated by DETECTOR_THICKNESS) on the
# walls + ceiling, veto scintillator on the floor. The ACTIVE DECAY
# VOLUME is the air interior bounded by the inner tracker layer (top
# and sides) and the floor scintillator (bottom).
#
# Visualization layers (from outside in):
#   blue solid line     — outer tracker layer (at tunnel wall)
#   blue dashed line    — inner tracker layer (24 cm inside outer)
#   blue hatched band   — tracker shell between the two layers
#   red filled region   — active decay volume (interior)
#   green strip on y=0  — veto scintillator on the floor
#
# Note: mesh_fiducial is built from the inner-tracker profile and the
# tunnel floor, so its volume IS the active decay volume — NOT the
# tracker shell.

ax4 = fig.add_subplot(224)

outer = tunnel_profile_points(inset=0)
inner = tunnel_profile_points(inset=DETECTOR_THICKNESS)
outer_closed = np.vstack([outer, outer[0]])
inner_closed = np.vstack([inner, inner[0]])

# Active decay volume (interior of inner tracker layer + floor)
ax4.fill(inner_closed[:, 0], inner_closed[:, 1], alpha=0.35, color='red',
         label='Active decay volume')

# Tracker shell (between outer and inner layers) drawn as a hatched ring
shell_x = np.concatenate([outer_closed[:, 0], inner_closed[::-1, 0]])
shell_y = np.concatenate([outer_closed[:, 1], inner_closed[::-1, 1]])
ax4.fill(shell_x, shell_y, facecolor='none', edgecolor='steelblue',
         hatch='///', linewidth=0, alpha=0.7,
         label=f'Tracker shell ({DETECTOR_THICKNESS*100:.0f} cm, 2 layers)')

# The two tracker-layer outlines
ax4.plot(outer_closed[:, 0], outer_closed[:, 1], '-', color='navy',
         linewidth=2, label='Outer tracker layer (at wall)')
ax4.plot(inner_closed[:, 0], inner_closed[:, 1], '--', color='navy',
         linewidth=1.5, label='Inner tracker layer')

# Veto scintillator: a thin band sitting just below the floor of the
# active volume (full tunnel width at y = floor)
floor_y = outer[:, 1].min()
scint_thick = 0.08
ax4.fill_between([-TUNNEL_ALPHA / 2, TUNNEL_ALPHA / 2],
                 floor_y - scint_thick, floor_y,
                 color='green', alpha=0.6, label='Veto scintillator (floor)')

# Dimension annotations (kept; surface labels removed to avoid clutter
# now that everything is named in the legend)
half_w = TUNNEL_ALPHA / 2
bot_y = outer[:, 1].min()
top_y = outer[:, 1].max()

# Horizontal: floor width (placed below the scintillator)
ax4.annotate('', xy=(half_w, bot_y - 0.45), xytext=(-half_w, bot_y - 0.45),
             arrowprops=dict(arrowstyle='<->', color='gray'))
ax4.text(0, bot_y - 0.70, f'{TUNNEL_ALPHA:.2f} m', ha='center', fontsize=8,
         color='gray')

# Vertical: total height (placed well outside the right wall labels)
dim_x = half_w + 0.55
ax4.annotate('', xy=(dim_x, top_y), xytext=(dim_x, bot_y),
             arrowprops=dict(arrowstyle='<->', color='gray'))
ax4.text(dim_x + 0.18, (top_y + bot_y) / 2, f'{TUNNEL_BETA:.2f} m',
         ha='left', fontsize=8, color='gray', rotation=90, va='center')

# Leave extra room above so the legend doesn't sit on the arch
ax4.set_xlim(-half_w - 1.1, half_w + 1.4)
ax4.set_ylim(bot_y - 0.95, top_y + 0.4)

ax4.set_xlabel('Transverse X (m)')
ax4.set_ylabel('Vertical Y (m)')
ax4.set_title('Cross-Section — tracker shell, active volume, floor veto')
ax4.grid(True, alpha=0.3)
ax4.set_aspect('equal', adjustable='box')
ax4.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.18),
           ncol=2, frameon=True, framealpha=0.95)

fig.suptitle('GRENDEL Tunnel Geometry — active decay volume bounded by '
             'inner tracker layer (walls/arch) and floor veto',
             fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig('tunnel_4panel.png', dpi=120, bbox_inches='tight')

# ── Print statistics ─────────────────────────────────────────────────

centreline_len = np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))

# Tunnel envelope (full cross-section) for shell-volume bookkeeping
import trimesh  # noqa: E402  (local — only used here)
verts_outer_mesh, faces_outer_mesh = create_profile_mesh(
    path_3d_fiducial, tunnel_profile_points(inset=0))
mesh_tunnel = trimesh.Trimesh(vertices=verts_outer_mesh,
                              faces=faces_outer_mesh)
mesh_tunnel.fix_normals()  # see grendel_geometry.build_fiducial_mesh
shell_volume = mesh_tunnel.volume - mesh_fiducial.volume

print("\nTunnel statistics:")
print(f"  Centreline length:           {centreline_len:.1f} m")
print(f"  Active decay volume:         {mesh_fiducial.volume:.1f} m^3 "
      f"(between inner tracker and floor)")
print(f"  Tracker shell volume:        {shell_volume:.1f} m^3 "
      f"(between the two tracker layers, {DETECTOR_THICKNESS*100:.0f} cm)")
print(f"  Tunnel envelope volume:      {mesh_tunnel.volume:.1f} m^3")
print(f"  Active-volume surface area:  {mesh_fiducial.area:.1f} m^2")
print(f"  Y offset (height):           {Y_POSITION} m")
print(f"  Cross-section:               {TUNNEL_ALPHA:.2f} x {TUNNEL_BETA:.2f} m "
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
ax.set_title('GRENDEL Tunnel — Active Decay Volume (inner tracker layer extrusion)',
             fontsize=13)
ax.view_init(elev=20, azim=45)
set_equal_aspect_3d(ax, expand=(1.0, 1.2, 1.0))  # +20% on CMS-Z (beam)
ax.legend(fontsize=10)

plt.show()
