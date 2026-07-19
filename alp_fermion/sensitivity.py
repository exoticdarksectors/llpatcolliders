#!/usr/bin/env python3
"""BC10 sensitivity: model-complete (m_a, 1/f) exclusion island for GRENDEL.

Unlike the model-agnostic (BR, ctau) scan that higgs/ uses, BC10 is
coupling-controlled: the single inverse decay constant 1/f fixes the production
yield, the lifetime ctau, and the visible BRs simultaneously -- exactly like the
HNL mixing |U|^2.  We exploit that by mapping the coupling to the HNL scan
variable

    u2 = ((1/f) / (1/f_ref))^2         (1/f_ref = model.INV_F_REF)

so that, off the single reference point at which production.py weights the a's
and templates.py stores ctau,

    production yield  ~ u2          (BR(B->K a) ~ (1/f)^2)
    ctau              ~ 1/u2        (Gamma_tot ~ (1/f)^2)
    visible BR ratios  independent of 1/f.

This is precisely the structure of the imported
``analysis.decay_reco_acceptance.scan_u2``, so the acceptance MC
(``build_event_mc`` -> best-two-track -> higgs/reco_common -> PR#13
``selection_mask``) and lifetime reweighting are reused verbatim. BC10-specific
production, decay templates, pole handling, and disconnected-component
extraction provide the inputs around that shared core. N_signal >= 3,
background-free, 3000 fb^-1. Each supported component has a lower edge (too
little production) and an upper edge (the ALP decays before reaching PX56);
unsupported eta/eta-prime pole rows are not bridged.

Usage:
    python -m alp_fermion.sensitivity                 # full grid
    python -m alp_fermion.sensitivity --mass 1.0 2.0  # subset
    python -m alp_fermion.sensitivity --plot-only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_ALP_ROOT = Path(__file__).resolve().parent
if str(_ALP_ROOT) not in sys.path:
    sys.path.insert(0, str(_ALP_ROOT))

import model  # noqa: E402
from paths import (  # noqa: E402
    LLP_VECTORS_DIR, TEMPLATE_DIR, ANALYSIS_DIR, GEOM_CACHE_DIR, add_hnl_to_path,
)
from alp_production import ALP_MASS_GRID, alp_csv_path  # noqa: E402

add_hnl_to_path()
from config_mass_grid import format_mass_for_filename  # noqa: E402
from analysis.constants import L_INT_PB, N_THRESHOLD, CMS_ORIGIN  # noqa: E402
from analysis.format_bridge import load_combined_csv  # noqa: E402
from analysis._engine import (  # noqa: E402
    compute_geometry, _eta_phi_to_directions_batch, _get_mesh,
)
from analysis.decay_reco_acceptance import (  # noqa: E402
    build_event_mc,
    scan_u2,
    signal_contribution_diagnostics,
)
from analysis.exclusion import find_exclusion_band_refined  # noqa: E402


# u2 = (1/f / 1/f_ref)^2 scan.  Wide enough to bracket the closed island for the
# whole mass range (1/f_ref = 1e-3 GeV^-1, so u2 in [1e-16,1] maps 1/f in
# [1e-11, 1e-3] GeV^-1).
LOG_U2_MIN, LOG_U2_MAX, N_U2 = -16.0, 0.0, 300
DECAY_SAMPLES = 60


def _write_checkpoint(rows, path):
    """Atomically persist completed mass rows for interruption-safe resume."""
    frame = pd.DataFrame(rows).sort_values("mass_GeV")
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def _invf_from_u2(u2):
    """Map the scan variable back to the physical coupling 1/f [GeV^-1]."""
    return model.INV_F_REF * np.sqrt(u2)


def _geometry(m_a, eta, phi, mesh, force=False, source_mtime=None):
    label = format_mass_for_filename(m_a)
    cache = GEOM_CACHE_DIR / f"geom_{label}.npz"
    fresh = cache.exists() and not force
    if fresh and source_mtime is not None:
        # invalidate when the input CSV is newer than the cache (a regenerated
        # pool must force a fresh ray-cast, not serve stale hits)
        fresh = cache.stat().st_mtime >= source_mtime
    if fresh:
        d = np.load(cache)
        return d["hits"].astype(bool), d["entry_d"], d["exit_d"]
    hits, entry_d, exit_d = compute_geometry(eta, phi, mesh, CMS_ORIGIN,
                                             batch_label=f"[mA_{label}]")
    GEOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, hits=hits, entry_d=entry_d, exit_d=exit_d)
    return hits, entry_d, exit_d


def process_mass(
    m_a,
    mesh,
    decay_samples=DECAY_SAMPLES,
    force_geom=False,
    reco_seed_offset=0,
):
    resonance = model.excluded_light_meson_resonance(m_a)
    if resonance is not None:
        return {
            "mass_GeV": m_a,
            "n_events": 0,
            "n_hits": 0,
            "has_sensitivity": False,
            "peak_N": np.nan,
            "invf_min": np.nan,
            "invf_max": np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
            "peak_invf": np.nan,
            "exclusion_reason": f"unsupported {resonance} resonance",
        }
    csv = alp_csv_path(m_a)
    if not csv.exists() or csv.stat().st_size == 0:
        return None
    data = load_combined_csv(csv, m_a)
    n_events = len(data["weight"])
    if n_events == 0:
        return None

    hits, entry_d, exit_d = _geometry(m_a, data["eta"], data["phi"], mesh,
                                      force_geom, source_mtime=csv.stat().st_mtime)
    n_hits = int(hits.sum())
    base = {"mass_GeV": m_a, "n_events": n_events, "n_hits": n_hits}
    if n_hits == 0:
        return {**base, "has_sensitivity": False, "peak_N": 0.0,
                "invf_min": np.nan, "invf_max": np.nan,
                "invf_min_open": False, "invf_max_open": False,
                "peak_invf": np.nan}

    tmpl_path = TEMPLATE_DIR / f"templates_{format_mass_for_filename(m_a)}.npz"
    if not tmpl_path.exists():
        print(
            f"  m_a={m_a:.3f}: missing templates; run "
            "generate_decay_templates_pythia.py first"
        )
        return None
    templates = np.load(tmpl_path)
    ctau_ref = float(templates["ctau_m_u2eq1"])
    if ctau_ref <= 0:
        return None

    idx = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
    direction = _eta_phi_to_directions_batch(data["eta"][idx], data["phi"][idx])
    p_mag = data["beta_gamma"][idx] * m_a
    energy = data["gamma"][idx] * m_a
    p4 = np.column_stack([energy, p_mag[:, None] * direction])
    # Memory-bounded acceptance MC: reco_common intermediates scale with
    # n_events * n_samples. Chunking keeps that allocation bounded and assigns
    # each chunk a deterministic sub-seed, making the published path exactly
    # reproducible for a fixed chunk size and configuration.
    EVENT_CHUNK = 256
    d_parts, passed_parts, template_index_parts = [], [], []
    base_seed = int(round(m_a * 1000)) * 100_000 + int(reco_seed_offset)
    for ci, cs in enumerate(range(0, len(idx), EVENT_CHUNK)):
        ce = min(cs + EVENT_CHUNK, len(idx))
        rng_c = np.random.default_rng(base_seed + ci)
        d_c, passed_c, template_index_c = build_event_mc(
            p4[cs:ce], direction[cs:ce],
            entry_d[idx[cs:ce]], exit_d[idx[cs:ce]],
            templates, decay_samples, rng_c, origin=CMS_ORIGIN)
        d_parts.append(d_c)
        passed_parts.append(passed_c)
        template_index_parts.append(template_index_c)
    d = np.concatenate(d_parts, axis=0)
    passed = np.concatenate(passed_parts, axis=0)
    template_indices = np.concatenate(template_index_parts, axis=0)

    # New Pythia templates sample the full exclusive BR mixture, including
    # neutral modes, so reconstruction itself supplies the visible fraction.
    # Retain compatibility only for explicitly requested legacy-proxy caches.
    full_branching = (
        "includes_full_branching" in templates.files
        and bool(templates["includes_full_branching"])
    )
    branching_factor = 1.0 if full_branching else model.visible_fraction(m_a)
    weights = data["weight"][idx] * branching_factor
    sample_weights = None
    if "matrix_element_weight" in templates.files:
        sample_weights = np.asarray(templates["matrix_element_weight"])[
            template_indices
        ]

    u2_grid = np.logspace(LOG_U2_MIN, LOG_U2_MAX, N_U2)
    u2_grid, N_grid = scan_u2(
        d, passed, exit_d[idx] - entry_d[idx], weights,
        data["beta_gamma"][idx], ctau_ref, L_INT_PB, u2_grid,
        sample_w=sample_weights)

    def evaluate(u2):
        _, signal = scan_u2(
            d, passed, exit_d[idx] - entry_d[idx], weights,
            data["beta_gamma"][idx], ctau_ref, L_INT_PB,
            np.asarray([u2]), sample_w=sample_weights,
        )
        return signal[0]

    band = find_exclusion_band_refined(
        u2_grid, N_grid, evaluate, N_THRESHOLD,
    )
    # Map the u2 island edges -> physical coupling 1/f.  Lower-u2 edge = too
    # little production -> SMALLER 1/f; higher-u2 edge = decays too early ->
    # LARGER 1/f.  So invf_min comes from u2_min, invf_max from u2_max.
    res = {
        **base,
        "has_sensitivity": band["has_sensitivity"],
        "peak_N": band["peak_N"],
        "peak_invf": _invf_from_u2(band["peak_u2"]),
        "invf_min": _invf_from_u2(band["u2_min"]) if np.isfinite(band["u2_min"]) else np.nan,
        "invf_max": _invf_from_u2(band["u2_max"]) if np.isfinite(band["u2_max"]) else np.nan,
        "invf_min_open": bool(band["u2_min_open"]),
        "invf_max_open": bool(band["u2_max_open"]),
    }
    diagnostic_points = {
        "invf_min": band["u2_min"],
        "peak": band["peak_u2"],
        "invf_max": band["u2_max"],
    }
    for label, u2 in diagnostic_points.items():
        if not np.isfinite(u2):
            for field in ("sample_ess", "event_ess", "max_event_fraction"):
                res[f"{label}_{field}"] = np.nan
            continue
        diagnostics = signal_contribution_diagnostics(
            d,
            passed,
            exit_d[idx] - entry_d[idx],
            weights,
            data["beta_gamma"][idx],
            ctau_ref,
            u2,
            sample_w=sample_weights,
        )
        for field in ("sample_ess", "event_ess", "max_event_fraction"):
            res[f"{label}_{field}"] = diagnostics[field]
    return res


def run(
    masses,
    force_geom=False,
    output=None,
    resume=False,
    reco_seed_offset=0,
    decay_samples=DECAY_SAMPLES,
):
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    out = Path(output) if output is not None else ANALYSIS_DIR / "bc10_sensitivity.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    mesh = _get_mesh()
    rows = []
    completed = set()
    if resume and out.exists():
        previous = pd.read_csv(out)
        rows = previous.to_dict("records")
        completed = set(previous["mass_GeV"].astype(float))
    for m_a in masses:
        if float(m_a) in completed:
            print(f"  m_a={m_a:.3f}: already checkpointed", flush=True)
            continue
        r = process_mass(
            m_a,
            mesh,
            decay_samples=decay_samples,
            force_geom=force_geom,
            reco_seed_offset=reco_seed_offset,
        )
        if r is None:
            print(f"  m_a={m_a:.3f}: skipped (no input)", flush=True)
            continue
        rows.append(r)
        # A mass point is expensive. Persist every completed row so an external
        # ROOT/Pythia or Python failure never discards the rest of the campaign.
        _write_checkpoint(rows, out)
        if r.get("exclusion_reason"):
            print(
                f"  m_a={m_a:.3f}: excluded ({r['exclusion_reason']})",
                flush=True,
            )
            continue
        if r["has_sensitivity"]:
            print(f"  m_a={m_a:.3f}: ISLAND peak_N={r['peak_N']:.1f}  "
                  f"1/f in [{r['invf_min']:.2e}, {r['invf_max']:.2e}] GeV^-1"
                  + ("  (min open)" if r["invf_min_open"] else "")
                  + ("  (max open)" if r["invf_max_open"] else ""), flush=True)
        else:
            print(f"  m_a={m_a:.3f}: no sensitivity (peak_N={r['peak_N']:.2f})",
                  flush=True)
    if not rows:
        print("No results.")
        return None
    df = pd.DataFrame(rows).sort_values("mass_GeV")
    _write_checkpoint(rows, out)
    print(f"\nSaved {out}")
    n_sens = int(df["has_sensitivity"].sum())
    print(f"Sensitivity at {n_sens}/{len(df)} masses; "
          f"max peak_N = {df['peak_N'].max():.1f}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="BC10 (m_a,1/f) island")
    ap.add_argument("--mass", type=float, nargs="+", default=None)
    ap.add_argument("--force-geometry", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--output", type=Path, default=None,
                    help="checkpoint/output CSV (default: tmp analysis CSV)")
    ap.add_argument("--resume", action="store_true",
                    help="keep rows already present in --output and skip them")
    ap.add_argument(
        "--reco-seed-offset", type=int, default=0,
        help="add this offset to every deterministic per-mass reconstruction seed",
    )
    ap.add_argument(
        "--decay-samples", type=int, default=DECAY_SAMPLES,
        help="decay/reconstruction samples per detector-entering ALP",
    )
    args = ap.parse_args(argv)
    if args.decay_samples <= 0:
        ap.error("--decay-samples must be positive")

    out = args.output or ANALYSIS_DIR / "bc10_sensitivity.csv"
    if not args.plot_only:
        masses = args.mass if args.mass else ALP_MASS_GRID
        print("BC10 fermiophilic-ALP sensitivity")
        print(f"  masses: {len(masses)}  L = {L_INT_PB/1e3:.0f} fb^-1  "
              f"N_thr = {N_THRESHOLD}")
        out = run(masses, force_geom=args.force_geometry,
                  output=out, resume=args.resume,
                  reco_seed_offset=args.reco_seed_offset,
                  decay_samples=args.decay_samples)
        if out is None:
            return 1
    if out and Path(out).exists():
        from plot import plot_island
        plot_island(out, ANALYSIS_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
