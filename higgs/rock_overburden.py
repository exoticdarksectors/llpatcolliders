"""
Rock overburden between the CMS interaction point and the GRENDEL tunnel.

For each point along the tunnel centreline, the straight line of sight from the
CMS IP passes first through the air of the experimental cavern (UXC55) and then
through rock until it reaches the tunnel. This script computes that rock
thickness as a function of position along the tunnel and produces
`rock_overburden.png`.

Cavern model
------------
UXC55 is approximated as a cylinder whose long axis runs along the beam (Z),
following the CMS convention of `grendel_geometry`:
    - radius   R_CAV  = 13.25 m   (26.5 m diameter)
    - half-len HALF_L = 26.5  m   (53 m long)
The IP does NOT sit at the cavern's transverse centre: the beam line sits low in
the cavern (~8.5 m above the floor). The vertical offset IP_BELOW is anchored to
the milliQan Letter of Intent (arXiv:1607.04669), whose PX56-gallery site is
33 m from the IP with 17 m of rock (=> 16 m of cavern air) at 43.1 deg
elevation. IP_BELOW = 4.7 m reproduces that (16.0 m air / 17.0 m rock).

Only the *transverse* (vertical) IP offset matters: the tunnel is directly
overhead, so every line of sight leaves the cavern through its curved side, never
the end-caps -- the 53 m length and the longitudinal IP position are irrelevant
to the result (verified: 0/47 tunnel nodes exit via an end-cap).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from grendel_geometry import path_3d_fiducial, cumulative_length

# ---------------------------------------------------------------------------
# Cavern parameters (current, milliQan-anchored assumptions)
# ---------------------------------------------------------------------------
R_CAV = 26.5 / 2.0    # cavern radius (m)
HALF_L = 53.0 / 2.0   # cavern half-length along beam Z (m)
IP_BELOW = 4.7        # IP sits this far below the cavern vertical centre (m)
MIN_ROCK = 7.0        # floor on the overburden (m). The cavern is not
                      # perfectly cylindrical, so near the thinnest (overhead)
                      # point the true rock (~7 m) exceeds the bare cylinder
                      # (~4.5 m). A floor corrects that region while leaving the
                      # milliQan-anchored points (rock >> 7 m) untouched -- a
                      # flat additive term cannot do this without an unphysical
                      # IP position (it would preserve the too-large cylinder
                      # variation between the minimum and the milliQan point).


def cavern_exit_distance(u):
    """Distance from the IP to the cavern boundary along unit vector u.

    Cavern = cylinder of radius R_CAV about the beam axis (Z), with its circular
    cross-section centred at Y = +IP_BELOW (i.e. the IP lies IP_BELOW below the
    axis). Returns the smaller of the curved-side and end-cap intersection
    distances.
    """
    ux, uy, uz = u
    a = ux * ux + uy * uy                      # transverse-to-Z component^2
    if a < 1e-12:
        t_side = np.inf
    else:
        b = -2.0 * uy * IP_BELOW
        c = IP_BELOW ** 2 - R_CAV ** 2
        t_side = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
    t_cap = np.inf if abs(uz) < 1e-12 else HALF_L / abs(uz)
    return min(t_side, t_cap)


def compute_overburden(path):
    """Return (rock, total_los, cavern_air) arrays for each point on `path`."""
    rock = np.empty(len(path))
    total = np.empty(len(path))
    air = np.empty(len(path))
    for i, P in enumerate(path):
        d = np.linalg.norm(P)
        u = P / d
        de = cavern_exit_distance(u)
        total[i] = d
        air[i] = de
        rock[i] = max(d - de, MIN_ROCK)
    return rock, total, air


def main():
    path = path_3d_fiducial
    s = cumulative_length
    rock, total, air = compute_overburden(path)

    print("Rock overburden (IP %.1f m below cavern centre, milliQan-anchored):"
          % IP_BELOW)
    print("  range %.1f - %.1f m, mean %.1f m" % (rock.min(), rock.max(),
                                                  rock.mean()))
    print("  thinnest at s=%.1f m (Z=%.1f), thickest at s=%.1f m (X=%.1f)"
          % (s[np.argmin(rock)], path[np.argmin(rock), 2],
             s[np.argmax(rock)], path[np.argmax(rock), 0]))

    fig, ax = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left: rock / air / total vs arc-length ---------------------------
    ax[0].plot(s, rock, "b-", lw=2.5, label="Rock")
    ax[0].plot(s, total, "k--", lw=1, alpha=.6,
               label="Total IP-tunnel line of sight")
    ax[0].plot(s, air, "g:", lw=1.5, label="Cavern air")
    ax[0].axhline(rock.mean(), color="b", ls=":", alpha=.4)
    ax[0].text(2, rock.mean() + 0.6, "mean %.1f m" % rock.mean(),
               color="b", fontsize=9)
    ax[0].set_xlabel("Arc-length along tunnel centreline s (m)")
    ax[0].set_ylabel("Distance (m)")
    ax[0].set_title("Rock overburden vs tunnel position\n"
                    "(IP %.1f m below cavern centre, milliQan-anchored)"
                    % IP_BELOW)
    ax[0].grid(alpha=.3)
    ax[0].legend()

    # --- Right: top view (X-Z) coloured by rock thickness -----------------
    sc = ax[1].scatter(path[:, 0], path[:, 2], c=rock, cmap="viridis", s=45)
    ax[1].scatter(0, 0, color="red", s=140, marker="*", label="CMS IP")
    # Beam-parallel (Z-axis) cylinder projects to a rectangle in the top view.
    ax[1].add_patch(Rectangle((-R_CAV, -HALF_L), 2 * R_CAV, 2 * HALF_L,
                              fill=False, edgecolor="red", alpha=.4, lw=1.5,
                              label="Cavern footprint (%.1f x %.0f m)"
                              % (2 * R_CAV, 2 * HALF_L)))
    ax[1].set_xlabel("X (m)")
    ax[1].set_ylabel("Z beam (m)")
    ax[1].set_aspect("equal")
    ax[1].set_title("Top view (X-Z), coloured by rock thickness")
    ax[1].legend(fontsize=8)
    plt.colorbar(sc, ax=ax[1], label="Rock (m)")

    plt.tight_layout()
    out = "rock_overburden.png"
    fig.savefig(out, dpi=120)
    print("Saved %s" % out)


if __name__ == "__main__":
    main()
