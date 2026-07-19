"""B4 -- Bc production-normalization nuisance on the lower edge.

The Bc channel owns a large share of the accepted lower-edge yield near the dome
(C9: ~28-44% over 1.8-3.6 GeV), and its normalization is the frozen
`SIGMA_BC_PB`, anchored to the LHCb f(Bc)/f(B) ratio but carrying a real
uncertainty (`SIGMA_BC_REL_UNCERT`, default 0.40). Scaling the Bc channel by
`(1 +/- delta_Bc)` moves the lower edge `u2_min` wherever Bc contributes.

Because the signal is additive in channels, this is a linear reweight of one
channel -- no full re-run. For each (flavor, mass) we scan the *combined* and the
*Bc-only* events over the U^2 grid and read off
  u2_min_central = boundary(N_tot),
  u2_min_hi      = boundary(N_tot + delta*N_Bc)   (Bc up   -> stronger limit),
  u2_min_lo      = boundary(N_tot - delta*N_Bc)   (Bc down -> weaker  limit).

  HNL_TMP_DIR=tmp HNL_RUN_TAG=<central_run> \
    python -m analysis.bc_nuisance --out tmp/runs/bc_nuisance_band.csv

SCOPE: this bands the DIRECT Bc channel (Bc -> l N) only. SIGMA_BC_PB also
normalizes the Bc -> tau nu source inside `induced_tau`, but that source is just
0.062% of the tau pool (sigma_Bc << sigma_charm/bottom, even with the large
BR(Bc->tau nu)), so co-varying it shifts induced_tau by <0.03% -- negligible
versus the ~0.2 dex direct-Bc band, even where induced_tau dominates (Utau low
mass). The band is therefore direct-Bc only and is not understated in practice.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HNL_ROOT = Path(__file__).resolve().parent.parent
if str(HNL_ROOT) not in sys.path:
    sys.path.insert(0, str(HNL_ROOT))

from config_mass_grid import format_mass_for_filename                  # noqa: E402
from production.constants import SIGMA_BC_REL_UNCERT                   # noqa: E402
from production.paths import LLP_VECTORS_DIR                          # noqa: E402
from analysis.constants import (                                      # noqa: E402
    L_INT_PB, LOG_U2_MAX, LOG_U2_MIN, N_THRESHOLD, N_U2_POINTS)
from analysis.exclusion import find_exclusion_band                   # noqa: E402
from analysis.format_bridge import load_combined_csv                 # noqa: E402
from analysis.decay_reco_acceptance import build_event_mc, scan_u2   # noqa: E402
from analysis._engine import (                               # noqa: E402
    _get_mesh, _select_hit_sample, _seed_for, _eta_phi_to_directions_batch,
    compute_geometry, load_decay_templates)

DEFAULT_MASSES = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6]


def _scan_csv(csv, flavor, mass, tag, mesh, templates, u2_grid,
              max_hit_events, decay_samples, seed_salt):
    """N(U^2) for one channel/combined CSV (zeros if empty)."""
    if not csv.exists() or csv.stat().st_size == 0:
        return np.zeros(len(u2_grid))
    data = load_combined_csv(csv, mass)
    if len(data["weight"]) == 0:
        return np.zeros(len(u2_grid))
    hits, entry_d, exit_d = compute_geometry(data["eta"], data["phi"], mesh)
    idx_all = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
    if len(idx_all) == 0:
        return np.zeros(len(u2_grid))
    rng = np.random.default_rng(_seed_for(
        flavor, f"{format_mass_for_filename(mass)}_{tag}", seed_salt))
    idx, w, _ = _select_hit_sample(idx_all, data["weight"][idx_all], max_hit_events, rng)
    direction = _eta_phi_to_directions_batch(data["eta"][idx], data["phi"][idx])
    p4 = np.column_stack([data["gamma"][idx] * mass,
                          (data["beta_gamma"][idx] * mass)[:, None] * direction])
    d, passed, _ = build_event_mc(p4, direction, entry_d[idx], exit_d[idx],
                                  templates, decay_samples, rng)
    _, N = scan_u2(d, passed, exit_d[idx] - entry_d[idx], w,
                   data["beta_gamma"][idx], float(templates["ctau_m_u2eq1"]),
                   L_INT_PB, u2_grid)
    return N


def run(flavors, masses, delta_bc, max_hit_events, decay_samples,
        seed_salt=""):
    mesh = _get_mesh()
    u2_grid = np.logspace(LOG_U2_MIN, LOG_U2_MAX, N_U2_POINTS)
    rows = []
    for flavor in flavors:
        for mass in masses:
            ml = format_mass_for_filename(mass)
            templates = load_decay_templates(flavor, ml)
            if templates is None:
                continue
            base = LLP_VECTORS_DIR / flavor
            N_tot = _scan_csv(base / "combined" / f"mN_{ml}.csv", flavor, mass, "comb",
                              mesh, templates, u2_grid, max_hit_events,
                              decay_samples, seed_salt)
            N_bc = _scan_csv(base / "Bc" / f"mN_{ml}.csv", flavor, mass, "Bc",
                             mesh, templates, u2_grid, max_hit_events,
                             decay_samples, seed_salt)

            def u2min(N):
                r = find_exclusion_band(u2_grid, N, N_THRESHOLD)
                return r["u2_min"] if r["has_sensitivity"] else np.nan

            c = u2min(N_tot)
            hi = u2min(N_tot + delta_bc * N_bc)   # Bc up   -> lower u2_min (stronger)
            lo = u2min(N_tot - delta_bc * N_bc)   # Bc down -> higher u2_min (weaker)
            # Bc share of the accepted yield at the central lower edge
            f_bc = np.nan
            if np.isfinite(c) and c > 0:
                iu = int(np.argmin(np.abs(u2_grid - c)))
                f_bc = N_bc[iu] / N_tot[iu] if N_tot[iu] > 0 else 0.0
            rows.append({"flavor": flavor, "mass_GeV": mass, "delta_bc": delta_bc,
                         "u2_min": c, "u2_min_bc_lo": lo, "u2_min_bc_hi": hi,
                         "bc_frac_at_u2min": f_bc,
                         "decay_samples": decay_samples,
                         "max_hit_events": max_hit_events,
                         "seed_salt": seed_salt})
            shift = (np.log10(lo / hi) if np.isfinite(lo) and np.isfinite(hi) and hi > 0
                     else np.nan)
            print(f"  {flavor:5s} m={mass:<4} Bc_frac={f_bc*100:4.0f}%  "
                  f"u2_min={c:.3e}  band=[{hi:.3e},{lo:.3e}]  width={shift:.3f} dex",
                  flush=True)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", nargs="+", type=float, default=None)
    ap.add_argument("--delta-bc", type=float, default=SIGMA_BC_REL_UNCERT)
    ap.add_argument("--max-hit-events", type=int, default=4000)
    ap.add_argument("--decay-samples", type=int, default=50)
    ap.add_argument(
        "--seed-salt", default="",
        help="Optional deterministic salt for independent numerical-control repeats")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    df = run(args.flavor, args.mass or DEFAULT_MASSES, args.delta_bc,
             args.max_hit_events, args.decay_samples, args.seed_salt)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} rows, delta_bc={args.delta_bc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
