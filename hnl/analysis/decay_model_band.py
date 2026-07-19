"""Decay-model uncertainty band on the HNL exclusion contour.

The decay-model band is ONE coherent nuisance -- the seam-derived hadronic-width
uncertainty ``delta(m)`` (``analysis/width_band.py``) -- driven through both:

* **lifetime leg (axis a):** ``ctau = hbar/Gamma_tot`` -> ``P_decay``; re-scan
  the central run's already-built acceptance MC with ``ctau`` rescaled by
  ``1/(1+/-delta)``. Moves the upper edge ``u2_max``.
* **composition leg (axis b):** ``Gamma_had`` shifts the visible/invisible mix
  -> ``vis_frac``; reweight the decay templates by mode from the SAME varied
  widths. (Added on top of the lifetime leg; shares ``Gamma_tot`` and partly
  self-cancels, so propagated coherently with the lifetime leg -- never summed
  in quadrature as an independent axis.)

Scope: production (4-vectors) and geometry are reused unchanged; reconstruction
selection / detector response / background / ``N>=3`` are the partner layer and
are held fixed (idealized). The band therefore covers production + decay physics
only -- see the provenance note on the figure.

  HNL_TMP_DIR=tmp HNL_RUN_TAG=<central_run> \
    python -m analysis.decay_model_band --flavor Ue Umu Utau --out band.csv
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

from analysis import width_band                                    # noqa: E402
from analysis._engine import process_mass_point, _get_mesh  # noqa: E402

# densified at 3.0-4.4 GeV where the dome closes (band changes fast there).
# All entries MUST be in config_mass_grid.MASS_GRID or the central run has no
# 4-vectors for them and the point is dropped (see the warning in run()).
DEFAULT_MASSES = [0.305, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8,
                  3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4]


def _f(v):
    return float(v) if v is not None else float("nan")


def run(flavors, masses, decay_samples, max_hit_events, event_chunk,
        seed_salt) -> pd.DataFrame:
    mesh = _get_mesh()
    rows = []
    for flavor in flavors:
        for mass in masses:
            d, hf = width_band.delta_and_hadfrac(flavor, mass)
            r = process_mass_point(
                flavor, mass, mesh,
                decay_samples=decay_samples,
                max_hit_events=max_hit_events,
                event_chunk=event_chunk,
                seed_salt=seed_salt,
                width_delta=d,
                had_frac=hf,
            )
            if r is None:
                print(f"  {flavor:5s} m={mass}: SKIPPED -- no 4-vectors "
                      f"(mass not in MASS_GRID / not produced in the central run)",
                      flush=True)
                continue
            if not r.get("has_sensitivity", False):
                continue
            rows.append(r)
            print(f"  {flavor:5s} m={mass:<4} delta={d:.3f} had_frac={hf:.2f}  "
                  f"u2_max={_f(r.get('u2_max')):.3e}  "
                  f"decay-band=[{_f(r.get('u2_max_dm_lo')):.3e}, "
                  f"{_f(r.get('u2_max_dm_hi')):.3e}]", flush=True)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", nargs="+", type=float, default=None)
    ap.add_argument("--decay-samples", type=int, default=50)
    ap.add_argument("--max-hit-events", type=int, default=None,
                    help="Optional weighted-resampling cap; omit for exact hits")
    ap.add_argument("--event-chunk", type=int, default=1000,
                    help="Exact-hit memory chunk (default: 1000)")
    ap.add_argument("--seed-salt", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    if args.decay_samples < 1:
        ap.error("--decay-samples must be >= 1")
    if args.max_hit_events is not None and args.max_hit_events < 1:
        ap.error("--max-hit-events must be >= 1 when provided")
    if args.event_chunk is not None and args.event_chunk < 1:
        ap.error("--event-chunk must be >= 1 when provided")

    df = run(args.flavor, args.mass or DEFAULT_MASSES,
             args.decay_samples, args.max_hit_events, args.event_chunk,
             args.seed_salt)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out} ({len(df)} rows with sensitivity)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
