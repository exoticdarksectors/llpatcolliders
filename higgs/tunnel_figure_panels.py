"""
Render the two tunnel schematic panels (3D view and top view) as SEPARATE
image files, for use as standalone figures:
    tunnel_3d.png       -- 3D view of the tunnel + CMS IP + axes
    tunnel_topview.png  -- top view (X-Z plane) with tunnel width band

Geometry is imported from grendel_geometry (CMS convention: X horizontal
transverse, Y up, Z beam). Physics-neutral schematic; rock/concrete not shown.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d proj)

from grendel_geometry import (
    mesh_fiducial, path_3d_fiducial, tunnel_profile_points, DETECTOR_THICKNESS,
)

profile_outer = tunnel_profile_points(inset=0)
path = path_3d_fiducial
origin = np.array([0.0, 0.0, 0.0])


def set_equal_aspect_3d(ax, expand=(1.0, 1.0, 1.0)):
    setters = (ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d)
    lims = (ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d())
    ranges = []
    for (lo, hi), f, setter in zip(lims, expand, setters):
        c, h = 0.5 * (lo + hi), 0.5 * (hi - lo) * f
        setter(c - h, c + h)
        ranges.append(2 * h)
    ax.set_box_aspect(ranges)


def make_3d(outfile="tunnel_3d.png"):
    fig = plt.figure(figsize=(8, 7))
    ax1 = fig.add_subplot(111, projection="3d")

    # Subsampled wireframe of the fiducial mesh (Y<->Z swapped for display so
    # matplotlib's vertical axis shows CMS Y (up)).
    rng = np.random.default_rng(0)
    n_sample = min(1000, len(mesh_fiducial.faces))
    face_idx = rng.choice(len(mesh_fiducial.faces), n_sample, replace=False)
    for fi in face_idx:
        tri = mesh_fiducial.vertices[mesh_fiducial.faces[fi]]
        tri = np.vstack([tri, tri[0]])
        ax1.plot(tri[:, 0], tri[:, 2], tri[:, 1], "b-", alpha=0.1, linewidth=0.5)

    ax1.plot(path[:, 0], path[:, 2], path[:, 1], "r-", linewidth=3,
             label="Centreline")
    ax1.scatter(origin[0], origin[2], origin[1], color="green", s=200,
                marker="o", label="Origin (IP)")

    arrow = 10
    ax1.quiver(0, 0, 0, arrow, 0, 0, color="red", arrow_length_ratio=0.1)
    ax1.quiver(0, 0, 0, 0, 0, arrow, color="green", arrow_length_ratio=0.1)
    ax1.quiver(0, 0, 0, 0, arrow, 0, color="blue", arrow_length_ratio=0.1)
    ax1.text(arrow, 0, 0, "X", fontsize=11)
    ax1.text(0, 0, arrow, "Y (up)", fontsize=11)
    ax1.text(0, arrow, 0, "Z (beam)", fontsize=11)

    ax1.set_xlabel("X (m)")
    ax1.set_ylabel("Z (beam, m)")
    ax1.set_zlabel("Y (up, m)")
    set_equal_aspect_3d(ax1, expand=(1.0, 1.2, 1.0))
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print("Saved", outfile)


def make_topview(outfile="tunnel_topview.png"):
    fig, ax2 = plt.subplots(figsize=(8, 7))

    # Tunnel width as a solid band: centreline +/- half-width perpendicular to
    # the path in the X-Z plane.
    half_width = profile_outer[:, 0].max()
    xz = path[:, [0, 2]]
    tang = np.gradient(xz, axis=0)
    tang /= np.linalg.norm(tang, axis=1, keepdims=True)
    perp = np.column_stack([-tang[:, 1], tang[:, 0]])
    left = xz + half_width * perp
    right = xz - half_width * perp
    band = np.vstack([left, right[::-1]])
    ax2.fill(band[:, 0], band[:, 1], color="blue", alpha=0.3,
             label="Tunnel width")

    ax2.plot(path[:, 0], path[:, 2], "r-", linewidth=2, label="Centreline")
    ax2.scatter(origin[0], origin[2], color="green", s=200, marker="o",
                label="Origin (IP)")

    ax2.set_xlabel("X (m)")
    ax2.set_ylabel("Z (beam, m)")
    ax2.set_title("Top View (X-Z Plane)")
    ax2.grid(True, alpha=0.3)
    ax2.axis("equal")
    ax2.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print("Saved", outfile)


if __name__ == "__main__":
    make_3d()
    make_topview()
