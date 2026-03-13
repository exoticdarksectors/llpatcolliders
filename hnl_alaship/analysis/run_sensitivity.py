#!/usr/bin/env python3
"""
analysis/run_sensitivity.py

Main FairShip-backed GARGOYLE HNL sensitivity driver.

Supports all three flavors (Ue, Umu, Utau) and uses:
- Shared geometry from REPO_ROOT/geometry/gargoyle_geometry.py
- FairShip hnl.py for decay widths / cτ (vendored)
- FairShip TPythia8/ROOT for rest-frame decay templates
- Full 3D Lorentz boost along actual parent 3-momentum direction
- Daughter ray-cast to GARGOYLE mesh for acceptance
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent  # llpatcolliders root
sys.path.insert(0, str(PROJECT_ROOT))

from analysis.constants import (
    CMS_ORIGIN,
    DEFAULT_ANALYSIS_FLAVORS,
    DEFAULT_DECAYS_PER_BIN,
    DEFAULT_DECAY_TEMPLATES,
    DEFAULT_MOTHER_SAMPLES,
    DEFAULT_POSITION_BINS,
    DEFAULT_RANDOM_SEED,
    FLAVORS,
    FONLL_MASS_MAX,
    L_INT_PB,
    LOG_U2_MAX,
    LOG_U2_MIN,
    N_THRESHOLD,
    N_U2_POINTS,
)
from analysis.decay_acceptance import load_or_build_decay_cache
from analysis.exclusion import find_exclusion_band
from analysis.fairship_decay import ctau_u2eq1
from analysis.format_bridge import load_combined_csv
from analysis.plot_exclusion import plot_exclusion, plot_nsignal_vs_u2
from analysis.sensitivity import scan_u2
from config_mass_grid import MASS_GRID, format_mass_for_filename

CSV_DIR = PROJECT_ROOT / "output" / "llp_4vectors"
CTAU_DIR = PROJECT_ROOT / "output" / "ctau"
OUTPUT_DIR = PROJECT_ROOT / "output" / "analysis"
GEOM_CACHE_DIR = OUTPUT_DIR / "geometry_cache"

_WORKER_MESH = None


def _geometry_helpers():
    geom_dir = str(REPO_ROOT / "geometry")
    if geom_dir not in sys.path:
        sys.path.insert(0, geom_dir)
    from gargoyle_geometry import mesh_fiducial, momenta_to_directions, ray_intersection_distances

    return mesh_fiducial, momenta_to_directions, ray_intersection_distances


def _get_mesh():
    mesh_fiducial, _, _ = _geometry_helpers()
    return mesh_fiducial


def compute_geometry(px, py, pz, mesh, origin=CMS_ORIGIN, batch_label=""):
    """Batch ray-cast HNL mothers from the CMS IP to the GARGOYLE fiducial mesh."""
    _, momenta_to_directions, ray_intersection_distances = _geometry_helpers()

    momenta = np.column_stack([px, py, pz])
    directions, valid = momenta_to_directions(momenta)
    origin_arr = np.asarray(origin, dtype=np.float64)

    hits = np.zeros(len(momenta), dtype=bool)
    entry_d = np.full(len(momenta), np.nan)
    exit_d = np.full(len(momenta), np.nan)

    candidates = np.flatnonzero(valid & (directions[:, 1] > 0.01))
    if len(candidates) == 0:
        return hits, entry_d, exit_d

    print(f"  Batch ray-casting {len(candidates)}/{len(momenta)} candidates {batch_label}...", flush=True)
    origins = np.repeat(origin_arr[None, :], len(candidates), axis=0)
    hit_local, first_d, second_d = ray_intersection_distances(mesh, origins, directions[candidates])
    good = hit_local & np.isfinite(second_d)

    hits[candidates[good]] = True
    entry_d[candidates[good]] = first_d[good]
    exit_d[candidates[good]] = second_d[good]

    n_hits = int(hits.sum())
    print(f"  {n_hits}/{len(momenta)} events hit detector ({100.0 * n_hits / max(len(momenta), 1):.2f}%)", flush=True)
    if n_hits > 0:
        print(f"  Mean path length: {(exit_d[hits] - entry_d[hits]).mean():.2f} m", flush=True)
    return hits, entry_d, exit_d


def _input_hash(px, py, pz):
    import hashlib

    digest = hashlib.sha256()
    for arr in (px, py, pz):
        digest.update(np.ascontiguousarray(arr).tobytes())
    return digest.hexdigest()[:16]


def load_or_compute_geometry(flavor, mass_label, px, py, pz, mesh, force=False):
    cache_dir = GEOM_CACHE_DIR / flavor
    cache_path = cache_dir / f"geom_{mass_label}.npz"
    input_sig = _input_hash(px, py, pz)

    if cache_path.exists() and not force:
        data = np.load(cache_path)
        cached_sig = str(data["input_sig"]) if "input_sig" in data else ""
        if cached_sig == input_sig and len(data["hits"]) == len(px):
            return data["hits"].astype(bool), data["entry_d"], data["exit_d"]

    hits, entry_d, exit_d = compute_geometry(px, py, pz, mesh, origin=CMS_ORIGIN, batch_label=f"[{flavor}/{mass_label}]")
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_path, hits=hits, entry_d=entry_d, exit_d=exit_d, input_sig=np.array(input_sig))
    return hits, entry_d, exit_d


def write_fairship_diagnostics(masses, flavor="Umu"):
    """Dump FairShip-backed cτ table for diagnostics."""
    CTAU_DIR.mkdir(parents=True, exist_ok=True)
    ctau_table = np.column_stack([masses, [ctau_u2eq1(mass, flavor=flavor) for mass in masses]])

    np.savetxt(
        CTAU_DIR / f"ctau_{flavor}.dat",
        ctau_table,
        header="mN_GeV  ctau_meters_at_Usq1",
        fmt=["%.4f", "%.6e"],
    )


def process_mass_point(
    flavor,
    mass,
    mesh,
    force_geometry=False,
    force_decay_cache=False,
    save_diagnostic=False,
    mother_samples=DEFAULT_MOTHER_SAMPLES,
    position_bins=DEFAULT_POSITION_BINS,
    decay_templates=DEFAULT_DECAY_TEMPLATES,
    decays_per_bin=DEFAULT_DECAYS_PER_BIN,
    random_seed=DEFAULT_RANDOM_SEED,
):
    """Process a single mass point."""
    mass_label = format_mass_for_filename(mass)
    csv_path = CSV_DIR / flavor / "combined" / f"mN_{mass_label}.csv"
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return None

    data = load_combined_csv(csv_path, mass)
    n_events = len(data["weight"])
    if n_events == 0:
        return None

    hits, entry_d, exit_d = load_or_compute_geometry(
        flavor=flavor,
        mass_label=mass_label,
        px=data["px"],
        py=data["py"],
        pz=data["pz"],
        mesh=mesh,
        force=force_geometry,
    )
    n_hits = int(hits.sum())
    if n_hits == 0:
        return {
            "mass_GeV": mass,
            "flavor": flavor,
            "u2_min": np.nan,
            "u2_max": np.nan,
            "peak_N": 0.0,
            "peak_u2": np.nan,
            "has_sensitivity": False,
            "n_events": n_events,
            "n_hits": 0,
            "n_sampled": 0,
            "decay_estimator": "empty",
        }

    decay_cache = load_or_build_decay_cache(
        flavor=flavor,
        mass_GeV=mass,
        mass_label=mass_label,
        data=data,
        hits=hits,
        entry_d=entry_d,
        exit_d=exit_d,
        mesh=mesh,
        mother_samples=mother_samples,
        position_bins=position_bins,
        decay_templates=decay_templates,
        decays_per_bin=decays_per_bin,
        random_seed=random_seed,
        force=force_decay_cache,
    )
    ctau_m = ctau_u2eq1(mass, flavor=flavor)
    u2_grid, n_grid = scan_u2(
        cache=decay_cache,
        ctau_u2_1=ctau_m,
        L_int_pb=L_INT_PB,
        log_u2_min=LOG_U2_MIN,
        log_u2_max=LOG_U2_MAX,
        n_points=N_U2_POINTS,
    )

    result = find_exclusion_band(u2_grid, n_grid, N_THRESHOLD)
    result.update(
        {
            "mass_GeV": mass,
            "flavor": flavor,
            "n_events": n_events,
            "n_hits": n_hits,
            "n_sampled": decay_cache["n_sampled"],
            "decay_estimator": decay_cache["estimator"],
        }
    )

    if save_diagnostic:
        diag_dir = OUTPUT_DIR / "diagnostics"
        plot_nsignal_vs_u2(u2_grid, n_grid, mass, flavor, diag_dir)

    return result


def _worker_init():
    global _WORKER_MESH
    _WORKER_MESH = _get_mesh()


def _worker_process_point(args):
    (
        flavor,
        mass,
        force_geometry,
        force_decay_cache,
        save_diagnostic,
        mother_samples,
        position_bins,
        decay_templates,
        decays_per_bin,
        random_seed,
    ) = args
    t0 = time.time()
    result = process_mass_point(
        flavor=flavor,
        mass=mass,
        mesh=_WORKER_MESH,
        force_geometry=force_geometry,
        force_decay_cache=force_decay_cache,
        save_diagnostic=save_diagnostic,
        mother_samples=mother_samples,
        position_bins=position_bins,
        decay_templates=decay_templates,
        decays_per_bin=decays_per_bin,
        random_seed=random_seed,
    )
    elapsed = time.time() - t0
    tag = f"{flavor}/mN_{format_mass_for_filename(mass)}"
    if result and result.get("has_sensitivity"):
        print(f"  {tag}: sensitivity found ({elapsed:.1f}s)", flush=True)
    elif result:
        print(f"  {tag}: done, no sensitivity ({elapsed:.1f}s)", flush=True)
    else:
        print(f"  {tag}: skipped ({elapsed:.1f}s)", flush=True)
    return result


def run(
    flavors,
    masses,
    n_workers=1,
    force_geometry=False,
    force_decay_cache=False,
    save_diagnostics=False,
    mother_samples=DEFAULT_MOTHER_SAMPLES,
    position_bins=DEFAULT_POSITION_BINS,
    decay_templates=DEFAULT_DECAY_TEMPLATES,
    decays_per_bin=DEFAULT_DECAYS_PER_BIN,
    random_seed=DEFAULT_RANDOM_SEED,
):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    for flavor in flavors:
        write_fairship_diagnostics(masses, flavor=flavor)

    work_items = [
        (
            flavor,
            mass,
            force_geometry,
            force_decay_cache,
            save_diagnostics,
            mother_samples,
            position_bins,
            decay_templates,
            decays_per_bin,
            random_seed,
        )
        for flavor in flavors
        for mass in masses
    ]

    n_total = len(work_items)
    print(f"\nProcessing {n_total} mass points with {n_workers} worker(s)...\n", flush=True)
    t_start = time.time()
    status_path = OUTPUT_DIR / "scan_status.json"
    results = []

    def _write_status(done, n_sens):
        elapsed = time.time() - t_start
        eta = (elapsed / max(done, 1)) * (n_total - done) if done > 0 else 0.0
        status_path.write_text(
            json.dumps(
                {
                    "ts": datetime.now().isoformat(),
                    "done": done,
                    "total": n_total,
                    "n_sensitive": n_sens,
                    "elapsed_s": round(elapsed, 1),
                    "eta_s": round(eta, 1),
                },
                indent=2,
            )
        )

    if n_workers <= 1:
        mesh = _get_mesh()
        for i, item in enumerate(work_items, start=1):
            t0 = time.time()
            result = process_mass_point(
                flavor=item[0],
                mass=item[1],
                mesh=mesh,
                force_geometry=item[2],
                force_decay_cache=item[3],
                save_diagnostic=item[4],
                mother_samples=item[5],
                position_bins=item[6],
                decay_templates=item[7],
                decays_per_bin=item[8],
                random_seed=item[9],
            )
            elapsed = time.time() - t0
            if result is not None:
                results.append(result)
            n_sensitive = sum(1 for row in results if row.get("has_sensitivity"))
            _write_status(i, n_sensitive)
            print(f"  [{i}/{n_total}] {item[0]}/mN_{format_mass_for_filename(item[1])} ({elapsed:.1f}s)", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_worker_init) as pool:
            futures = {pool.submit(_worker_process_point, item): item for item in work_items}
            done_count = 0
            for future in as_completed(futures):
                done_count += 1
                result = future.result()
                if result is not None:
                    results.append(result)
                n_sensitive = sum(1 for row in results if row.get("has_sensitivity"))
                _write_status(done_count, n_sensitive)
                if done_count % 10 == 0 or done_count == n_total:
                    print(f"  [{done_count}/{n_total}] {n_sensitive} sensitive", flush=True)

    total_time = time.time() - t_start
    if not results:
        print("No results produced.")
        return

    results.sort(key=lambda row: (row["flavor"], row["mass_GeV"]))
    df = pd.DataFrame(results)
    out_csv = OUTPUT_DIR / "gargoyle_hnl_sensitivity.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved: {out_csv}")
    print(f"Sensitivity found at {int(df['has_sensitivity'].sum())}/{len(df)} mass points")
    print(f"Total time: {total_time:.0f}s ({total_time / 60.0:.1f} min)")

    meta_path = OUTPUT_DIR / "run_metadata.json"
    meta_path.write_text(
        json.dumps(
            {
                "timestamp": datetime.now().isoformat(),
                "flavors": flavors,
                "n_masses": len(masses),
                "n_workers": n_workers,
                "n_results": len(results),
                "n_sensitive": int(df["has_sensitivity"].sum()),
                "mass_range": [float(min(masses)), float(max(masses))],
                "total_time_s": round(total_time, 1),
                "mother_samples": mother_samples,
                "position_bins": position_bins,
                "decay_templates": decay_templates,
                "decays_per_bin": decays_per_bin,
                "random_seed": random_seed,
            },
            indent=2,
        )
    )

    plot_exclusion(out_csv, OUTPUT_DIR)


def main():
    parser = argparse.ArgumentParser(description="GARGOYLE HNL sensitivity calculation (FairShip decays, all flavors)")
    parser.add_argument("--flavor", nargs="+", default=None, choices=FLAVORS, help="Flavors to process (default: Umu)")
    parser.add_argument("--mass", nargs="+", type=float, default=None, help="Specific masses in GeV (default: full FONLL grid)")
    parser.add_argument("--plot-only", action="store_true", help="Re-plot from an existing sensitivity CSV")
    parser.add_argument("--force-geometry", action="store_true", help="Recompute the geometry cache")
    parser.add_argument("--force-decay-cache", action="store_true", help="Recompute the FairShip decay-acceptance cache")
    parser.add_argument("--diagnostics", action="store_true", help="Save N_signal vs U^2 diagnostic plots")
    parser.add_argument("--workers", type=int, default=1, help="Number of worker processes")
    parser.add_argument("--mother-samples", type=int, default=DEFAULT_MOTHER_SAMPLES)
    parser.add_argument("--position-bins", type=int, default=DEFAULT_POSITION_BINS)
    parser.add_argument("--decay-templates", type=int, default=DEFAULT_DECAY_TEMPLATES)
    parser.add_argument("--decays-per-bin", type=int, default=DEFAULT_DECAYS_PER_BIN)
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED)
    args = parser.parse_args()

    if args.plot_only:
        csv_path = OUTPUT_DIR / "gargoyle_hnl_sensitivity.csv"
        if not csv_path.exists():
            print(f"No results file found: {csv_path}")
            sys.exit(1)
        plot_exclusion(csv_path, OUTPUT_DIR)
        return

    flavors = args.flavor or DEFAULT_ANALYSIS_FLAVORS
    masses = args.mass or [mass for mass in MASS_GRID if mass <= FONLL_MASS_MAX]

    print("GARGOYLE HNL Sensitivity Analysis (FairShip pipeline)")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points ({min(masses):.2f}-{max(masses):.2f} GeV)")
    print(f"  Workers: {args.workers}")
    print(f"  Luminosity: {L_INT_PB:.0f} pb^-1 ({L_INT_PB / 1e3:.0f} fb^-1)")
    print(f"  Threshold: N_signal >= {N_THRESHOLD}")
    print("  Decay backend: FairShip (hnl.py + TPythia8)")
    print("  Boost: full 3D Lorentz (along parent 3-momentum)")

    run(
        flavors=flavors,
        masses=masses,
        n_workers=args.workers,
        force_geometry=args.force_geometry,
        force_decay_cache=args.force_decay_cache,
        save_diagnostics=args.diagnostics,
        mother_samples=args.mother_samples,
        position_bins=args.position_bins,
        decay_templates=args.decay_templates,
        decays_per_bin=args.decays_per_bin,
        random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()
