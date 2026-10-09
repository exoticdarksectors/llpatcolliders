"""Per-face hit fraction on the septagon tracker vs LLP mass (0.5 / 15 / 40 GeV).

For each LLP that decays in the fiducial volume, both e+e- daughters are
propagated to the septagon wall and classified onto one of the 7 faces
(floor + right wall = veto; left wall + 4 roof facets = tracker). Shows how the
instrument-able face set depends on mass: low-mass pairs are collimated and pile
onto the central roof, high-mass pairs spread across the arch and leak onto the
veto walls. Writes facefrac_<mass>.npy and face_hitmap3.png.
"""
import numpy as np, trimesh
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.cm as cm, matplotlib.colors as mcolors
from grendel_geometry import (path_3d_fiducial, create_profile_mesh,
    tunnel_profile_points, cache_geometry, local_transverse_xy,
    DETECTOR_THICKNESS, TUNNEL_ALPHA, TUNNEL_WALL_HEIGHT)

P_CUT = 0.100          # GeV/c, min softer-electron momentum
ORIGIN = [0, 0, 0]
N = 40                 # decay samples per LLP
rng = np.random.default_rng(1)

# --- septagon (inner/fiducial) cross-section: floor + 2 walls + 4 roof facets --
fid = tunnel_profile_points(inset=DETECTOR_THICKNESS)
hw = TUNNEL_ALPHA / 2 - DETECTOR_THICKNESS
yf = fid[:, 1].min(); yt = fid[:, 1].max(); ysp = yf + TUNNEL_WALL_HEIGHT
a = hw; b = yt - ysp
roof = np.array([[a * np.cos(p), ysp + b * np.sin(p)]
                 for p in np.linspace(0, np.pi, 5)])   # 4 roof facets
V = np.vstack([[-hw, yf], [hw, yf], roof])              # 7 vertices
FACE_NAMES = ["floor (veto)", "right wall (veto)", "roof-1", "roof-2",
              "roof-3", "roof-4", "left wall"]
TRACKER = {2, 3, 4, 5, 6}                               # roof facets + left wall
pv, pf = create_profile_mesh(path_3d_fiducial, V)
mesh = trimesh.Trimesh(vertices=pv, faces=pf); mesh.fix_normals()

MASSES = {"0.5 GeV": "LLP0p5GeVSmall.csv",
          "15 GeV":  "LLPSmall.csv",
          "40 GeV":  "LLP40GeVSmall.csv"}


def classify(pts):
    """Nearest-edge face index for each wall hit point (local cross-section)."""
    x, y = local_transverse_xy(pts); P = np.column_stack([x, y]); nV = len(V)
    best = np.full(len(P), 1e18); face = np.zeros(len(P), int)
    for i in range(nV):
        aa = V[i]; bb = V[(i + 1) % nV]; ab = bb - aa; L2 = ab @ ab
        t = np.clip(((P - aa) @ ab) / L2, 0, 1); proj = aa + t[:, None] * ab
        d2 = np.sum((P - proj) ** 2, axis=1)
        upd = d2 < best; best[upd] = d2[upd]; face[upd] = i
    return face


def first_hit(orig, dr):
    loc, ri, _ = mesh.ray.intersects_location(ray_origins=orig, ray_directions=dr)
    out = np.full((len(orig), 3), np.nan)
    if len(loc) == 0:
        return out
    d = np.linalg.norm(loc - orig[ri], axis=1)
    for k in range(len(orig)):
        m = ri == k
        if m.any():
            out[k] = loc[m][np.argmin(d[m])]
    return out


def run(csv):
    """Path-length-weighted per-face fraction of forward daughter hits."""
    gc = cache_geometry(csv, mesh, ORIGIN)
    O, D, W = [], [], []
    for idx in np.where(gc['hits'])[0]:
        en, ex = gc['entry_d'][idx], gc['exit_d'][idx]
        g, be, m = gc['gamma'][idx], gc['beta'][idx], gc['mass'][idx]
        z = gc['direction'][idx]; path = ex - en
        d = rng.uniform(en, ex, N); cs = rng.uniform(0, 1, N); sn = np.sqrt(1 - cs * cs)
        phi = rng.uniform(0, 2 * np.pi, N)
        p_soft = g * m / 2 * (1 - be * cs)
        pz1 = g * (cs + be); pz2 = g * (be - cs); pt = sn
        ref = np.array([0., 1., 0.]) if abs(z[1]) < 0.9 else np.array([1., 0., 0.])
        xh = np.cross(z, ref); xh /= np.linalg.norm(xh); yh = np.cross(z, xh)
        u = np.cos(phi)[:, None] * xh + np.sin(phi)[:, None] * yh
        n1 = np.sqrt(pz1 ** 2 + pt ** 2); n2 = np.sqrt(pz2 ** 2 + pt ** 2)
        dir1 = (pz1 / n1)[:, None] * z + (pt / n1)[:, None] * u
        dir2 = (pz2 / n2)[:, None] * z - (pt / n2)[:, None] * u
        dpos = d[:, None] * z
        fwd = (pz2 > 0) & (p_soft > P_CUT)   # softer daughter forward + energetic
        for j in range(N):
            if not fwd[j]:
                continue
            O.append(dpos[j]); D.append(dir1[j]); W.append(path)
            O.append(dpos[j]); D.append(dir2[j]); W.append(path)
    O = np.array(O); D = np.array(D); W = np.array(W)
    hp = first_hit(O, D); ok = ~np.isnan(hp[:, 0])
    faces = np.full(len(O), -1); faces[ok] = classify(hp[ok])
    frac = np.array([W[faces == i].sum() for i in range(len(V))])
    return frac / frac.sum()


def main():
    results = {}
    for label, csv in MASSES.items():
        print("\n== %s ==" % label)
        fr = run(csv); results[label] = fr
        for i, nm in enumerate(FACE_NAMES):
            print("  %-18s [%s]: %5.1f%%"
                  % (nm, "TRK" if i in TRACKER else "veto", 100 * fr[i]))
        print("  -> tracker %.1f%% | veto-lost %.1f%%"
              % (100 * sum(fr[i] for i in TRACKER),
                 100 * sum(fr[i] for i in range(len(V)) if i not in TRACKER)))
        np.save("facefrac_%s.npy" % label.replace(" ", ""), fr)

    # --- 3-panel figure: cross-section faces coloured by hit fraction ---
    norm = mcolors.Normalize(0, 45); cmap = cm.get_cmap('inferno')
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))
    for ax, label in zip(axes, MASSES):
        f = results[label]
        ax.plot(*np.vstack([fid, fid[:1]]).T, color='0.6', lw=1)
        for i in range(len(V)):
            A = V[i]; B = V[(i + 1) % len(V)]
            ax.add_collection(LineCollection([[A, B]],
                              colors=[cmap(norm(100 * f[i]))], linewidths=9))
            mid = (A + B) / 2; off = mid / (np.linalg.norm(mid) + 1e-9) * 0.33
            ax.text(mid[0] + off[0], mid[1] + off[1], '%.0f%%' % (100 * f[i]),
                    ha='center', va='center', fontsize=8, weight='bold')
        ax.set_title('m=%s   tracker %.0f%% / veto-lost %.0f%%'
                     % (label, 100 * f[2:].sum(), 100 * (f[0] + f[1])))
        ax.set_aspect('equal'); ax.set_xlabel('X (m)'); ax.grid(alpha=.25)
    axes[0].set_ylabel('Y (m)')
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    fig.colorbar(sm, ax=axes, label='% of forward daughter hits', fraction=.025)
    fig.suptitle('Septagon per-face hit fraction vs mass', fontsize=13)
    fig.savefig("face_hitmap3.png", dpi=130, bbox_inches='tight')
    print("\nsaved face_hitmap3.png")


if __name__ == "__main__":
    main()
