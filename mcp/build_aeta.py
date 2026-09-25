#!/usr/bin/env python3
"""
Build A(eta) live from the GRENDEL geometry in ../higgs/grendel_geometry.py.

A(eta) is the azimuthal fraction, at each pseudorapidity, of directions from
the CMS IP that cross the detector.  It is the only geometry input to this
directory; analyze.py and per_event.py fold the chi eta spectrum through it.

    python3 build_aeta.py                      # DEFAULT: whole fiducial
                                               # volume -> external/A_eta.npy
    python3 build_aeta.py --surface tracker    # arch/ceiling + right wall only
    python3 build_aeta.py --inset 0            # un-inset tunnel wall

Output is the same (2, 48) layout as external/A_eta.npy: row 0 bin centres,
row 1 the fraction.

Surface definitions come from grendel_geometry itself:
    tracker  points_on_tracker(), i.e. TRACKER_SURFACES = arch/ceiling +
             right wall (the surfaces a particle from the IP exits through);
             the complement, floor + left wall, is the IP-facing veto.
    full     any crossing of the fiducial volume, ignoring which face.

The default is `full`: the whole fiducial volume, ignoring which face is
crossed.  The superseded hand-made array is kept as
external/A_eta_2026-09-02_preflip.npy -- it predates the beam-Z flip (commit
31f8512) so its eta axis is reversed.  Its closure (0.02384) matches the
corrected `tracker` option to 0.2%.  See the "Geometry provenance" section of
README.md.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "higgs"))


def build(surface="tracker", inset=0.24, n_per_bin=20000, nbins=48,
          eta_max=1.2, seed=12345):
    import grendel_geometry as G

    mesh = (G.mesh_fiducial if inset == G.DETECTOR_THICKNESS
            else G.build_fiducial_mesh(detector_thickness=inset)[0])
    rng = np.random.default_rng(seed)

    edges = np.linspace(-eta_max, eta_max, nbins + 1)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    d = edges[1] - edges[0]
    frac = np.zeros(nbins)

    for i, e0 in enumerate(ctr):
        eta = e0 + (rng.random(n_per_bin) - 0.5) * d
        phi = rng.random(n_per_bin) * 2 * np.pi
        dirs = np.array([G.eta_phi_to_direction(a, b) for a, b in zip(eta, phi)])
        org = np.zeros_like(dirs)
        if surface == "full":
            frac[i] = mesh.ray.intersects_any(org, dirs).mean()
            continue
        loc, iray, _ = mesh.ray.intersects_location(org, dirs, multiple_hits=True)
        if len(loc) == 0:
            continue
        frac[i] = len(np.unique(iray[G.points_on_tracker(loc)])) / n_per_bin

    return np.vstack([ctr, frac])


def closure(a):
    """Isotropic per-particle acceptance: int A(eta) * 0.5 sech^2(eta) d_eta."""
    ctr, f = a
    return float((f * 0.5 / np.cosh(ctr) ** 2 * (ctr[1] - ctr[0])).sum())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--surface", choices=("tracker", "full"), default="full")
    p.add_argument("--inset", type=float, default=0.24,
                   help="fiducial inset from the tunnel wall (m); 0 = wall itself")
    p.add_argument("-n", "--n-per-bin", type=int, default=20000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--out", default=None)
    a = p.parse_args()

    out = a.out or os.path.join(
        HERE, "external",
        "A_eta.npy" if a.surface == "full" else "A_eta_live_%s.npy" % a.surface)
    arr = build(a.surface, a.inset, a.n_per_bin, seed=a.seed)
    np.save(out, arr)
    c = closure(arr)
    print("wrote %s" % out)
    print("  surface=%s  inset=%.2f m  %d rays/bin" %
          (a.surface, a.inset, a.n_per_bin))
    print("  isotropic closure = %.5f  (%.4f sr, %.2f%% of 4pi)" %
          (c, c * 4 * np.pi, 100 * c))
    ref = os.path.join(HERE, "external", "A_eta.npy")
    if os.path.exists(ref) and os.path.abspath(out) != os.path.abspath(ref):
        print("  frozen external/A_eta.npy closure = %.5f  (ratio %.3f)"
              % (closure(np.load(ref)), c / closure(np.load(ref))))


if __name__ == "__main__":
    main()
