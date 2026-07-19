#!/usr/bin/env python3
"""
hnl/analysis/run_sensitivity.py

Main driver: compute HNL sensitivity from the combined 4-vector CSVs that
``hnl/production/`` writes, and produce the (m_N, U^2) money plot.

Default flavor is Umu only (HNL pattern 010 -- the muon-mixing scenario).

Usage examples:

    python -m analysis.run_sensitivity                       # Umu, full grid
    python -m analysis.run_sensitivity --flavor Ue Umu       # multi-flavor
    python -m analysis.run_sensitivity --mass 1.0 2.0        # custom masses
    python -m analysis.run_sensitivity --workers 3           # parallel workers
    python -m analysis.run_sensitivity --plot-only           # re-plot only

Inputs (read via ``production.paths``):

    {LLP_VECTORS_DIR}/{flavor}/combined/mN_{mass}.csv    -- 4-vectors per point
    {TMP_DIR}/decay_templates/{flavor}/templates_{mass}.npz  -- FairShip decays + ctau
    ../higgs/grendel_geometry.py  (sibling package, PR #13)  -- shared mesh

Outputs:

    {ANALYSIS_DIR}/hnl_sensitivity.csv
    {ANALYSIS_DIR}/hnl_exclusion.{pdf,png}
    {ANALYSIS_DIR}/geometry_cache/{flavor}/geom_{mass}.npz
    {ANALYSIS_DIR}/diagnostics/                          -- with --diagnostics
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd

HNL_ROOT = Path(__file__).resolve().parent.parent
if str(HNL_ROOT) not in sys.path:
    sys.path.insert(0, str(HNL_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.paths import ANALYSIS_DIR, LLP_VECTORS_DIR
from analysis.constants import (
    DEFAULT_FLAVORS,
    FLAVORS,
    ANALYSIS_MASS_MAX,
    L_INT_PB,
    N_THRESHOLD,
)
from analysis.plot_exclusion import plot_exclusion
from analysis.decay_reco_acceptance import P_CUT
from analysis._engine import (
    process_mass_point,
    _get_mesh,
    GEOM_CACHE_DIR,
    DECAY_SAMPLES,
    TEMPLATE_DIR,
)


# =========================================================================
# Parallel worker plumbing
# =========================================================================

_WORKER_MESH = None


def _worker_init():
    global _WORKER_MESH
    _WORKER_MESH = _get_mesh()


def _worker_process_point(args):
    (flavor, mass, force_geom, save_diag, decay_samples, max_hit_events,
     event_chunk, seed_salt) = args
    t0 = time.time()
    result = process_mass_point(
        flavor, mass, _WORKER_MESH,
        force_geometry=force_geom, save_diagnostic=save_diag,
        decay_samples=decay_samples, max_hit_events=max_hit_events,
        event_chunk=event_chunk, seed_salt=seed_salt)
    elapsed = time.time() - t0
    mass_label = format_mass_for_filename(mass)
    tag = f"{flavor}/mN_{mass_label}"
    if result and result.get("has_sensitivity"):
        print(f"  {tag}: sensitivity found ({elapsed:.1f}s)", flush=True)
    elif result:
        print(f"  {tag}: done, no sensitivity ({elapsed:.1f}s)", flush=True)
    else:
        print(f"  {tag}: skipped ({elapsed:.1f}s)", flush=True)
    return result


# =========================================================================
# Driver
# =========================================================================

def run(flavors, masses, n_workers=1, force_geometry=False,
        save_diagnostics=False, decay_samples=DECAY_SAMPLES,
        max_hit_events=None, event_chunk=None, seed_salt="",
        checkpoint_every=0):
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    work_items = []
    for flavor in flavors:
        for mass in masses:
            work_items.append((
                flavor, mass, force_geometry, save_diagnostics,
                decay_samples, max_hit_events, event_chunk, seed_salt))

    n_total = len(work_items)
    print(f"\nProcessing {n_total} mass points "
          f"with {n_workers} worker(s)...\n", flush=True)

    t_start = time.time()
    results = []
    status_path = ANALYSIS_DIR / "scan_status.json"
    partial_csv = ANALYSIS_DIR / "hnl_sensitivity.partial.csv"

    def _write_partial(done):
        if not results:
            return
        rows = sorted(results, key=lambda r: (r["flavor"], r["mass_GeV"]))
        pd.DataFrame(rows).to_csv(partial_csv, index=False)
        if checkpoint_every and done % checkpoint_every == 0:
            plot_exclusion(partial_csv, ANALYSIS_DIR,
                           basename="hnl_exclusion_partial")

    def _write_status(done, n_sens):
        elapsed = time.time() - t_start
        eta = (elapsed / max(done, 1)) * (n_total - done) if done > 0 else 0
        status = {
            "ts": datetime.now().isoformat(),
            "done": done, "total": n_total,
            "n_sensitive": n_sens,
            "elapsed_s": round(elapsed, 1),
            "eta_s": round(eta, 1),
        }
        status_path.write_text(json.dumps(status, indent=2))

    if n_workers <= 1:
        mesh = _get_mesh()
        for i, item in enumerate(work_items):
            flavor, mass, fg, sd, ds, mhe, echunk, salt = item
            t0 = time.time()
            r = process_mass_point(flavor, mass, mesh,
                                   force_geometry=fg, save_diagnostic=sd,
                                   decay_samples=ds, max_hit_events=mhe,
                                   event_chunk=echunk, seed_salt=salt)
            elapsed = time.time() - t0
            mass_label = format_mass_for_filename(mass)
            if r is not None:
                results.append(r)
                _write_partial(i + 1)
            n_sens = sum(1 for r in results if r.get("has_sensitivity"))
            _write_status(i + 1, n_sens)
            print(f"  [{i+1}/{n_total}] {flavor}/mN_{mass_label} "
                  f"({elapsed:.1f}s)", flush=True)
    else:
        with ProcessPoolExecutor(
            max_workers=n_workers, initializer=_worker_init
        ) as pool:
            futures = {
                pool.submit(_worker_process_point, item): i
                for i, item in enumerate(work_items)
            }
            done_count = 0
            for future in as_completed(futures):
                done_count += 1
                r = future.result()
                if r is not None:
                    results.append(r)
                    _write_partial(done_count)
                n_sens = sum(1 for r in results if r.get("has_sensitivity"))
                _write_status(done_count, n_sens)
                if done_count % 10 == 0 or done_count == n_total:
                    elapsed = time.time() - t_start
                    print(f"  [{done_count}/{n_total}] "
                          f"{n_sens} sensitive, "
                          f"{elapsed:.0f}s elapsed", flush=True)

    total_time = time.time() - t_start

    # Honest coverage accounting: every requested (flavor, mass) that produced
    # no result is recorded with the reason it was dropped, so a missing point
    # never silently disappears from the money plot.
    processed = {(r["flavor"], r["mass_GeV"]) for r in results}
    skipped = []
    for flavor in flavors:
        for mass in masses:
            if (flavor, mass) in processed:
                continue
            label = format_mass_for_filename(mass)
            csv = LLP_VECTORS_DIR / flavor / "combined" / f"mN_{label}.csv"
            tmpl = TEMPLATE_DIR / flavor / f"templates_{label}.npz"
            if not csv.exists() or csv.stat().st_size == 0:
                reason = "missing_or_empty_combined_csv"
            elif not tmpl.exists():
                reason = "missing_decay_templates"
            else:
                reason = "no_acceptance_or_nonpositive_ctau"
            skipped.append({"flavor": flavor, "mass_GeV": mass, "reason": reason})

    n_requested = len(flavors) * len(masses)
    n_processed = len(results)
    meta = {
        "timestamp": datetime.now().isoformat(),
        "flavors": flavors,
        "masses_GeV": [float(mass) for mass in masses],
        "n_masses": len(masses),
        "n_points_requested": n_requested,
        "n_points_processed": n_processed,
        "n_points_skipped": len(skipped),
        "skipped": skipped,
        "n_workers": n_workers,
        "decay_samples": decay_samples,
        "max_hit_events": max_hit_events,
        "event_chunk": event_chunk,
        "seed_salt": seed_salt,
        "track_momentum_cut_GeV": P_CUT,
        "checkpoint_every": checkpoint_every,
        "n_results": n_processed,
        "n_sensitive": 0,
        "mass_range": [float(min(masses)), float(max(masses))],
        "total_time_s": round(total_time, 1),
        "llp_vectors_dir": str(LLP_VECTORS_DIR),
    }

    if skipped:
        by_reason = {}
        for s in skipped:
            by_reason[s["reason"]] = by_reason.get(s["reason"], 0) + 1
        print(f"\nWARNING: {len(skipped)}/{n_requested} requested points "
              f"skipped (recorded in run_metadata.json):")
        for reason, n in sorted(by_reason.items()):
            print(f"    {n:4d}  {reason}")

    if not results:
        print("\nNo results produced -- every requested point was skipped.")
        (ANALYSIS_DIR / "run_metadata.json").write_text(json.dumps(meta, indent=2))
        return 0

    results.sort(key=lambda r: (r["flavor"], r["mass_GeV"]))
    df = pd.DataFrame(results)
    out_csv = ANALYSIS_DIR / "hnl_sensitivity.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved: {out_csv}")

    n_sens = int(df["has_sensitivity"].sum())
    meta["n_sensitive"] = n_sens
    print(f"Sensitivity found at {n_sens}/{len(df)} mass points")
    print(f"Total time: {total_time:.0f}s ({total_time/60:.1f} min)")

    (ANALYSIS_DIR / "run_metadata.json").write_text(json.dumps(meta, indent=2))

    plot_exclusion(out_csv, ANALYSIS_DIR)
    return n_processed


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="HNL sensitivity + (m_N, U^2) money plot")
    parser.add_argument(
        "--flavor", nargs="+", default=None, choices=FLAVORS,
        help=f"Flavors to process (default: {DEFAULT_FLAVORS})")
    parser.add_argument(
        "--mass", nargs="+", type=float, default=None,
        help=f"Specific masses in GeV (default: full grid <= {ANALYSIS_MASS_MAX} GeV)")
    parser.add_argument(
        "--plot-only", action="store_true",
        help="Re-plot from existing results CSV")
    parser.add_argument(
        "--force-geometry", action="store_true",
        help="Recompute geometry cache (ignore existing NPZ)")
    parser.add_argument(
        "--diagnostics", action="store_true",
        help="Save N_signal vs U^2 plots for each mass point")
    parser.add_argument(
        "--workers", type=int, default=3,
        help="Number of parallel workers (default: 3)")
    parser.add_argument(
        "--decay-samples", type=int,
        default=int(os.environ.get("HNL_DECAY_SAMPLES", DECAY_SAMPLES)),
        help=f"Decay positions sampled per hit event (default: {DECAY_SAMPLES}; env HNL_DECAY_SAMPLES)")
    parser.add_argument(
        "--max-hit-events", type=int,
        default=(int(os.environ["HNL_MAX_HIT_EVENTS"]) if os.environ.get("HNL_MAX_HIT_EVENTS") else None),
        help="Approximate mode: weighted-resample at most this many hit events per mass point (env HNL_MAX_HIT_EVENTS)")
    parser.add_argument(
        "--event-chunk", type=int,
        default=(int(os.environ["HNL_EVENT_CHUNK"]) if os.environ.get("HNL_EVENT_CHUNK") else None),
        help="Bound peak memory by processing this many selected hit events at a time without changing their statistical sample (env HNL_EVENT_CHUNK)")
    parser.add_argument(
        "--seed-salt", default=os.environ.get("HNL_ANALYSIS_SEED_SALT", ""),
        help="Optional deterministic salt for independent numerical-control repeats (env HNL_ANALYSIS_SEED_SALT)")
    parser.add_argument(
        "--mass-stride", type=int, default=int(os.environ.get("HNL_MASS_STRIDE", "1")),
        help="Approximate mode: keep every Nth mass-grid point (default: 1; env HNL_MASS_STRIDE)")
    parser.add_argument(
        "--mass-offset", type=int, default=0,
        help="Approximate mode: offset used with --mass-stride (default: 0)")
    parser.add_argument(
        "--checkpoint-every", type=int, default=0,
        help="Write partial CSV always; additionally plot partial result every N completed points")
    args = parser.parse_args(argv)

    if args.decay_samples < 1:
        parser.error("--decay-samples must be >= 1")
    if args.max_hit_events is not None and args.max_hit_events < 1:
        parser.error("--max-hit-events must be >= 1 when provided")
    if args.event_chunk is not None and args.event_chunk < 1:
        parser.error("--event-chunk must be >= 1 when provided")
    if args.mass_stride < 1:
        parser.error("--mass-stride must be >= 1")
    if args.mass_offset < 0 or args.mass_offset >= args.mass_stride:
        parser.error("--mass-offset must satisfy 0 <= offset < stride")
    if args.checkpoint_every < 0:
        parser.error("--checkpoint-every must be >= 0")

    if args.plot_only:
        csv_path = ANALYSIS_DIR / "hnl_sensitivity.csv"
        if not csv_path.exists():
            print(f"No results file found: {csv_path}")
            return 1
        plot_exclusion(csv_path, ANALYSIS_DIR)
        return 0

    flavors = args.flavor or DEFAULT_FLAVORS
    masses = args.mass or [m for m in MASS_GRID if m <= ANALYSIS_MASS_MAX]
    if args.mass_stride > 1:
        masses = masses[args.mass_offset::args.mass_stride]
    if not masses:
        parser.error("selected mass grid is empty")

    print("HNL Sensitivity Analysis")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points "
          f"({min(masses):.2f} - {max(masses):.2f} GeV)")
    print(f"  Workers: {args.workers}")
    print(f"  Decay samples: {args.decay_samples}")
    if args.max_hit_events:
        print(f"  Hit-event cap: {args.max_hit_events} (weighted resampling)")
    else:
        print("  Hit-event sample: exact (all detector hits)")
    if args.event_chunk:
        print(f"  Event chunk: {args.event_chunk} (memory bound only)")
    if args.seed_salt:
        print(f"  Analysis seed salt: {args.seed_salt}")
    if args.mass_stride > 1:
        print(f"  Mass stride: every {args.mass_stride} point(s), offset {args.mass_offset}")
    print(f"  Luminosity: {L_INT_PB:.0f} pb^-1 ({L_INT_PB/1e3:.0f} fb^-1)")
    print(f"  Threshold: N_signal >= {N_THRESHOLD}")
    print(f"  Inputs:    {LLP_VECTORS_DIR}")
    print(f"  Outputs:   {ANALYSIS_DIR}")

    n_processed = run(flavors, masses,
                      n_workers=args.workers,
                      force_geometry=args.force_geometry,
                      save_diagnostics=args.diagnostics,
                      decay_samples=args.decay_samples,
                      max_hit_events=args.max_hit_events,
                      event_chunk=args.event_chunk,
                      seed_salt=args.seed_salt,
                      checkpoint_every=args.checkpoint_every)
    if not n_processed:
        print("ERROR: analysis produced no results; see run_metadata.json "
              "for the per-point skip reasons.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
