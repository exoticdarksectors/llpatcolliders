#!/usr/bin/env python3
"""
analysis/run_sensitivity.py

Main driver: compute GARGOYLE HNL sensitivity and produce exclusion plots.

Usage:
    python analysis/run_sensitivity.py                          # Ue + Umu
    python analysis/run_sensitivity.py --flavor Ue              # single flavor
    python analysis/run_sensitivity.py --flavor Ue --mass 1.0   # single point
    python analysis/run_sensitivity.py --plot-only              # re-plot only
    python analysis/run_sensitivity.py --workers 8              # 8 parallel workers
"""

import sys
import argparse
import json
import time
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # hnl/
REPO_ROOT = PROJECT_ROOT.parent                        # llpatcolliders/
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from analysis.constants import (
    L_INT_PB, N_THRESHOLD, CMS_ORIGIN, FLAVORS, FONLL_MASS_MAX,
    LOG_U2_MIN, LOG_U2_MAX, N_U2_POINTS, M_ELECTRON,
)
from analysis.format_bridge import load_combined_csv
from analysis.sensitivity import scan_u2
from analysis.exclusion import find_exclusion_band
from analysis.plot_exclusion import plot_exclusion, plot_nsignal_vs_u2

# Directories
CSV_DIR = PROJECT_ROOT / "output" / "llp_4vectors"
CTAU_DIR = PROJECT_ROOT / "output" / "ctau"
OUTPUT_DIR = PROJECT_ROOT / "output" / "analysis"
GEOM_CACHE_DIR = OUTPUT_DIR / "geometry_cache"

# Default flavors for FONLL scan (tau on hold)
DEFAULT_FLAVORS = ["Ue", "Umu"]


# =========================================================================
# Geometry: ray-cast through GARGOYLE mesh (lazy import)
# =========================================================================

def _get_mesh():
    """Lazy-import GARGOYLE mesh from shared geometry/ directory."""
    geom_path = REPO_ROOT / "geometry"
    if not geom_path.exists():
        raise FileNotFoundError(
            f"geometry/ directory not found at {geom_path}.\n"
            "Copy gargoyle_geometry.py from the dark photons repo into geometry/."
        )
    sys.path.insert(0, str(geom_path))
    try:
        from gargoyle_geometry import mesh_fiducial
    except ImportError as e:
        raise ImportError(
            f"Failed to import gargoyle_geometry from {geom_path}: {e}\n"
            "Ensure gargoyle_geometry.py is present in geometry/."
        ) from e
    return mesh_fiducial


def _eta_phi_to_directions_batch(eta, phi):
    """Convert arrays of (eta, phi) to (N, 3) unit direction vectors."""
    theta = 2.0 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    return np.column_stack([dx, dy, dz])


def compute_geometry(eta, phi, mesh, origin=CMS_ORIGIN, batch_label=""):
    """
    Batch ray-cast all (eta, phi) directions against the GARGOYLE mesh.

    Returns hits (bool), entry_d (meters), exit_d (meters).
    """
    n = len(eta)
    origin_arr = np.array(origin, dtype=np.float64)
    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)

    # Compute all directions at once
    directions = _eta_phi_to_directions_batch(eta, phi)

    # Pre-filter: only ray-cast events with Y component > 0.01
    candidates = np.where(directions[:, 1] > 0.01)[0]
    n_cand = len(candidates)
    if n_cand == 0:
        return hits, entry_d, exit_d

    print(f"  Batch ray-casting {n_cand}/{n} candidates {batch_label}...",
          flush=True)

    # Batch ray-cast: all rays in one call
    cand_dirs = directions[candidates]
    origins = np.tile(origin_arr, (n_cand, 1))
    locations, ray_ids, _ = mesh.ray.intersects_location(
        ray_origins=origins, ray_directions=cand_dirs)

    if len(locations) == 0:
        print(f"  0/{n} events hit detector", flush=True)
        return hits, entry_d, exit_d

    # Compute distances from origin for each intersection
    dists = np.linalg.norm(locations - origin_arr, axis=1)

    # Group intersections by ray_id and extract entry/exit
    # Use sorting for efficient grouping
    order = np.argsort(ray_ids)
    sorted_ray_ids = ray_ids[order]
    sorted_dists = dists[order]

    # Find unique ray IDs and their boundaries
    unique_rays, start_idx, counts = np.unique(
        sorted_ray_ids, return_index=True, return_counts=True)

    # Only rays with >= 2 intersections have valid entry/exit
    valid = counts >= 2
    valid_rays = unique_rays[valid]
    valid_starts = start_idx[valid]
    valid_counts = counts[valid]

    for i in range(len(valid_rays)):
        ray_local = valid_rays[i]
        orig_idx = candidates[ray_local]
        s = valid_starts[i]
        e = s + valid_counts[i]
        ray_dists = sorted_dists[s:e]
        ray_dists.sort()
        hits[orig_idx] = True
        entry_d[orig_idx] = ray_dists[0]
        exit_d[orig_idx] = ray_dists[1]

    n_hits = int(hits.sum())
    print(f"  {n_hits}/{n} events hit detector ({n_hits / n * 100:.2f}%)",
          flush=True)
    if n_hits > 0:
        path_lens = exit_d[hits] - entry_d[hits]
        print(f"  Mean path length: {path_lens.mean():.2f} m", flush=True)
    return hits, entry_d, exit_d


def _input_hash(eta, phi):
    """Fast hash of input arrays for cache validation."""
    import hashlib
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(eta).tobytes())
    h.update(np.ascontiguousarray(phi).tobytes())
    return h.hexdigest()[:16]


def load_or_compute_geometry(flavor, mass_label, eta, phi, mesh,
                             force=False):
    """Load cached geometry or compute + cache to NPZ.

    Cache is invalidated when the input (eta, phi) arrays change.
    """
    cache_dir = GEOM_CACHE_DIR / flavor
    cache_path = cache_dir / f"geom_{mass_label}.npz"
    input_sig = _input_hash(eta, phi)

    if cache_path.exists() and not force:
        data = np.load(cache_path)
        cached_sig = str(data["input_sig"]) if "input_sig" in data else ""
        if cached_sig == input_sig and len(data["hits"]) == len(eta):
            return data["hits"].astype(bool), data["entry_d"], data["exit_d"]

    hits, entry_d, exit_d = compute_geometry(
        eta, phi, mesh, CMS_ORIGIN, batch_label=f"[{flavor}/{mass_label}]")

    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache_path, hits=hits, entry_d=entry_d, exit_d=exit_d,
                        input_sig=np.array(input_sig))
    return hits, entry_d, exit_d


# =========================================================================
# cτ tables
# =========================================================================

def load_ctau_table(flavor):
    """Load ctau(m_N) at U²=1. Returns {mass_GeV: ctau_meters}."""
    path = CTAU_DIR / f"ctau_{flavor}.dat"
    lines = path.read_text().strip().split("\n")
    table = {}
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split()
        table[float(parts[0])] = float(parts[1])
    return table


def lookup_ctau(ctau_table, mass):
    """Lookup ctau for a mass, with nearest-neighbor fallback."""
    if mass in ctau_table:
        return ctau_table[mass]
    masses = np.array(list(ctau_table.keys()))
    idx = np.argmin(np.abs(masses - mass))
    nearest = masses[idx]
    print(f"  WARNING: ctau not found for m_N={mass:.4f} GeV, "
          f"using nearest m_N={nearest:.4f} GeV")
    return ctau_table[nearest]


# =========================================================================
# BR_vis tables
# =========================================================================

def load_br_vis_table(flavor):
    """Load BR_vis(m_N) table. Returns {mass_GeV: br_vis}."""
    path = CTAU_DIR / f"br_vis_{flavor}.dat"
    if not path.exists():
        print(f"  WARNING: BR_vis table not found at {path}, using BR_vis=1.0")
        return None
    lines = path.read_text().strip().split("\n")
    table = {}
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.split()
        table[float(parts[0])] = float(parts[1])
    return table


def lookup_br_vis(br_vis_table, mass):
    """Lookup BR_vis for a mass. Returns 1.0 if table is None."""
    if br_vis_table is None:
        return 1.0
    if mass in br_vis_table:
        return br_vis_table[mass]
    masses = np.array(list(br_vis_table.keys()))
    idx = np.argmin(np.abs(masses - mass))
    nearest = masses[idx]
    return br_vis_table[nearest]


# =========================================================================
# Main pipeline
# =========================================================================

def process_mass_point(flavor, mass, ctau_table, mesh,
                       br_vis_table=None,
                       force_geometry=False, save_diagnostic=False):
    """Process a single (flavor, mass) point → exclusion band."""
    mass_label = format_mass_for_filename(mass)
    csv_path = CSV_DIR / flavor / "combined" / f"mN_{mass_label}.csv"

    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return None

    # Load 4-vectors
    data = load_combined_csv(csv_path, mass)
    n_events = len(data["weight"])
    if n_events == 0:
        return None

    # Geometry
    hits, entry_d, exit_d = load_or_compute_geometry(
        flavor, mass_label, data["eta"], data["phi"], mesh,
        force=force_geometry)

    n_hits = int(hits.sum())
    if n_hits == 0:
        return {
            "mass_GeV": mass, "flavor": flavor,
            "u2_min": np.nan, "u2_max": np.nan,
            "peak_N": 0.0, "peak_u2": np.nan,
            "has_sensitivity": False, "n_events": n_events, "n_hits": 0,
        }

    # cτ at U²=1
    ctau_u2_1 = lookup_ctau(ctau_table, mass)
    if ctau_u2_1 <= 0:
        return None

    # BR_vis
    br_vis = lookup_br_vis(br_vis_table, mass)

    # U² scan with vectorized 2-body acceptance
    u2_grid, N_grid = scan_u2(
        data["weight"], data["beta_gamma"], data["gamma"], data["beta"],
        hits, entry_d, exit_d,
        ctau_u2_1, mass, L_INT_PB,
        m_daughter=M_ELECTRON, use_acceptance=True,
        log_u2_min=LOG_U2_MIN, log_u2_max=LOG_U2_MAX,
        n_points=N_U2_POINTS,
        br_vis=br_vis)

    # Extract exclusion band
    result = find_exclusion_band(u2_grid, N_grid, N_THRESHOLD)
    result["mass_GeV"] = mass
    result["flavor"] = flavor
    result["n_events"] = n_events
    result["n_hits"] = n_hits

    # Optional diagnostic plot
    if save_diagnostic:
        diag_dir = OUTPUT_DIR / "diagnostics"
        plot_nsignal_vs_u2(u2_grid, N_grid, mass, flavor, diag_dir)

    return result


# =========================================================================
# Parallel worker functions
# =========================================================================

_WORKER_MESH = None


def _worker_init():
    """Each worker process loads its own mesh copy."""
    global _WORKER_MESH
    _WORKER_MESH = _get_mesh()


def _worker_process_point(args):
    """Worker function for parallel mass point processing."""
    flavor, mass, ctau_table, br_vis_table, force_geom, save_diag = args
    t0 = time.time()
    result = process_mass_point(
        flavor, mass, ctau_table, _WORKER_MESH,
        br_vis_table=br_vis_table,
        force_geometry=force_geom, save_diagnostic=save_diag)
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
# Run driver
# =========================================================================

def run(flavors, masses, n_workers=1, force_geometry=False,
        save_diagnostics=False):
    """Run the full sensitivity pipeline."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Build work items
    work_items = []
    for flavor in flavors:
        ctau = load_ctau_table(flavor)
        br_vis_tab = load_br_vis_table(flavor)
        for mass in masses:
            work_items.append(
                (flavor, mass, ctau, br_vis_tab, force_geometry, save_diagnostics))

    n_total = len(work_items)
    print(f"\nProcessing {n_total} mass points "
          f"with {n_workers} worker(s)...\n", flush=True)

    t_start = time.time()
    results = []

    # Status file for monitoring
    status_path = OUTPUT_DIR / "scan_status.json"

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
        # Sequential mode
        mesh = _get_mesh()
        for i, item in enumerate(work_items):
            flavor, mass, ctau, br_vis_tab, fg, sd = item
            t0 = time.time()
            r = process_mass_point(flavor, mass, ctau, mesh,
                                   br_vis_table=br_vis_tab,
                                   force_geometry=fg, save_diagnostic=sd)
            elapsed = time.time() - t0
            mass_label = format_mass_for_filename(mass)
            if r is not None:
                results.append(r)
            n_sens = sum(1 for r in results if r.get("has_sensitivity"))
            _write_status(i + 1, n_sens)
            print(f"  [{i+1}/{n_total}] {flavor}/mN_{mass_label} "
                  f"({elapsed:.1f}s)", flush=True)
    else:
        # Parallel mode
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
                n_sens = sum(1 for r in results if r.get("has_sensitivity"))
                _write_status(done_count, n_sens)
                if done_count % 10 == 0 or done_count == n_total:
                    elapsed = time.time() - t_start
                    print(f"  [{done_count}/{n_total}] "
                          f"{n_sens} sensitive, "
                          f"{elapsed:.0f}s elapsed", flush=True)

    total_time = time.time() - t_start

    if not results:
        print("No results produced.")
        return

    # Sort results by flavor + mass for consistent output
    results.sort(key=lambda r: (r["flavor"], r["mass_GeV"]))

    df = pd.DataFrame(results)
    out_csv = OUTPUT_DIR / "gargoyle_hnl_sensitivity.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nResults saved: {out_csv}")

    # Summary
    n_sens = df["has_sensitivity"].sum()
    print(f"Sensitivity found at {n_sens}/{len(df)} mass points")
    print(f"Total time: {total_time:.0f}s ({total_time/60:.1f} min)")

    # Write run metadata
    meta = {
        "timestamp": datetime.now().isoformat(),
        "flavors": flavors,
        "n_masses": len(masses),
        "n_workers": n_workers,
        "n_results": len(results),
        "n_sensitive": int(n_sens),
        "mass_range": [float(min(masses)), float(max(masses))],
        "total_time_s": round(total_time, 1),
    }
    meta_path = OUTPUT_DIR / "run_metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2))

    # Plot
    plot_exclusion(out_csv, OUTPUT_DIR)


def main():
    parser = argparse.ArgumentParser(
        description="GARGOYLE HNL sensitivity calculation")
    parser.add_argument(
        "--flavor", nargs="+", default=None,
        choices=FLAVORS,
        help="Flavors to process (default: Ue Umu)")
    parser.add_argument(
        "--mass", nargs="+", type=float, default=None,
        help="Specific masses in GeV (default: full FONLL grid)")
    parser.add_argument(
        "--plot-only", action="store_true",
        help="Re-plot from existing results CSV")
    parser.add_argument(
        "--force-geometry", action="store_true",
        help="Recompute geometry cache (ignore existing NPZ)")
    parser.add_argument(
        "--diagnostics", action="store_true",
        help="Save N_signal vs U² plots for each mass point")
    parser.add_argument(
        "--workers", type=int, default=6,
        help="Number of parallel workers (default: 6)")

    args = parser.parse_args()

    if args.plot_only:
        csv_path = OUTPUT_DIR / "gargoyle_hnl_sensitivity.csv"
        if not csv_path.exists():
            print(f"No results file found: {csv_path}")
            sys.exit(1)
        plot_exclusion(csv_path, OUTPUT_DIR)
        return

    flavors = args.flavor or DEFAULT_FLAVORS
    masses = args.mass or [m for m in MASS_GRID if m <= FONLL_MASS_MAX]

    print(f"GARGOYLE HNL Sensitivity Analysis")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points "
          f"({min(masses):.2f}–{max(masses):.2f} GeV)")
    print(f"  Workers: {args.workers}")
    print(f"  Luminosity: {L_INT_PB:.0f} pb⁻¹ ({L_INT_PB/1e3:.0f} fb⁻¹)")
    print(f"  Threshold: N_signal >= {N_THRESHOLD}")

    run(flavors, masses,
        n_workers=args.workers,
        force_geometry=args.force_geometry,
        save_diagnostics=args.diagnostics)


if __name__ == "__main__":
    main()
