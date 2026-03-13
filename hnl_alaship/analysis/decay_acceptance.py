"""
analysis/decay_acceptance.py

Cache a geometry-only HNL decay-acceptance kernel using FairShip-sampled decays.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent  # llpatcolliders root (for shared geometry)

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(REPO_ROOT / "geometry") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "geometry"))

from analysis.constants import (
    CMS_ORIGIN,
    DEFAULT_DECAYS_PER_BIN,
    DEFAULT_DECAY_TEMPLATES,
    DEFAULT_MOTHER_SAMPLES,
    DEFAULT_POSITION_BINS,
    DEFAULT_RANDOM_SEED,
    MIN_CHARGED_DAUGHTERS,
    P_CUT,
    SEP_MAX,
    SEP_MIN,
)
from analysis.fairship_decay import DEFAULT_DECAY_CONFIG, boost_decay_to_lab, sample_rest_frame_decays
from gargoyle_geometry import momenta_to_directions, ray_intersection_distances

OUTPUT_DIR = PROJECT_ROOT / "output" / "analysis"
DECAY_CACHE_DIR = OUTPUT_DIR / "decay_cache"


def _cache_signature(
    data,
    hits,
    entry_d,
    exit_d,
    mother_samples,
    position_bins,
    decay_templates,
    decays_per_bin,
    random_seed,
):
    digest = hashlib.sha256()
    payloads = (
        data["weight"],
        data["E"],
        data["px"],
        data["py"],
        data["pz"],
        hits.astype(np.int8),
        np.nan_to_num(entry_d, nan=-1.0),
        np.nan_to_num(exit_d, nan=-1.0),
    )
    for payload in payloads:
        digest.update(np.ascontiguousarray(payload).tobytes())
    digest.update(str(Path(DEFAULT_DECAY_CONFIG).resolve()).encode())
    if Path(DEFAULT_DECAY_CONFIG).exists():
        digest.update(Path(DEFAULT_DECAY_CONFIG).read_bytes())
    digest.update(
        (
            f"{mother_samples}:{position_bins}:{decay_templates}:{decays_per_bin}:"
            f"{random_seed}:{P_CUT}:{SEP_MIN}:{SEP_MAX}"
        ).encode()
    )
    return digest.hexdigest()[:20]


def _pairwise_acceptance(points):
    n_points = len(points)
    if n_points < MIN_CHARGED_DAUGHTERS:
        return False

    i_idx, j_idx = np.triu_indices(n_points, k=1)
    diffs = points[i_idx] - points[j_idx]
    dists = np.linalg.norm(diffs, axis=1)
    return bool(np.any((dists >= SEP_MIN) & (dists <= SEP_MAX)))


def _accept_decay_at_position(mesh, decay_distance, exit_distance, mother_direction, parent_p4, template):
    origin = np.asarray(CMS_ORIGIN, dtype=np.float64)
    decay_position = origin + mother_direction * decay_distance
    daughters = boost_decay_to_lab(parent_p4, template)
    momenta = np.column_stack([daughters["px"], daughters["py"], daughters["pz"]])
    p_abs = np.linalg.norm(momenta, axis=1)

    charged_mask = daughters["stable"] & (np.abs(daughters["charge"]) > 0.5) & (p_abs > P_CUT)
    if int(charged_mask.sum()) < MIN_CHARGED_DAUGHTERS:
        return False

    charged_momenta = momenta[charged_mask]
    charged_p = p_abs[charged_mask]
    directions = charged_momenta / charged_p[:, None]
    origins = np.repeat(decay_position[None, :], len(directions), axis=0)

    track_hits, _, _ = ray_intersection_distances(mesh, origins, directions, eps=1e-6)
    if int(track_hits.sum()) < MIN_CHARGED_DAUGHTERS:
        return False

    survivor_directions = directions[track_hits]
    survivor_directions = survivor_directions[(survivor_directions @ mother_direction) > 1e-9]
    if len(survivor_directions) < MIN_CHARGED_DAUGHTERS:
        return False

    plane_distance = exit_distance - decay_distance
    plane_t = plane_distance / (survivor_directions @ mother_direction)
    forward = plane_t > 0.0
    if int(forward.sum()) < MIN_CHARGED_DAUGHTERS:
        return False

    plane_points = decay_position[None, :] + plane_t[forward, None] * survivor_directions[forward]
    return _pairwise_acceptance(plane_points)


def _select_sampled_mothers(data, hit_mask, max_samples, rng):
    hit_indices = np.flatnonzero(hit_mask)
    if len(hit_indices) == 0:
        return {
            "estimator": "empty",
            "indices": np.empty(0, dtype=np.int64),
            "weights": np.empty(0, dtype=np.float64),
            "total_weight": 0.0,
        }

    hit_weights = np.asarray(data["weight"][hit_indices], dtype=np.float64)
    total_weight = float(hit_weights.sum())
    if total_weight <= 0.0:
        raise ValueError("All hit-mother weights are non-positive")

    if len(hit_indices) <= max_samples:
        return {
            "estimator": "exact",
            "indices": hit_indices,
            "weights": hit_weights,
            "total_weight": total_weight,
        }

    probabilities = hit_weights / total_weight
    sampled_indices = rng.choice(hit_indices, size=max_samples, replace=True, p=probabilities)
    return {
        "estimator": "weighted_resample",
        "indices": sampled_indices.astype(np.int64),
        "weights": np.empty(0, dtype=np.float64),
        "total_weight": total_weight,
    }


def _save_cache(cache_path, cache):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        input_sig=np.array(cache["input_sig"]),
        estimator=np.array(cache["estimator"]),
        total_hit_weight=np.array(cache["total_hit_weight"]),
        sample_weights=np.asarray(cache["sample_weights"], dtype=np.float64),
        beta_gamma=np.asarray(cache["beta_gamma"], dtype=np.float64),
        entry_d=np.asarray(cache["entry_d"], dtype=np.float64),
        exit_d=np.asarray(cache["exit_d"], dtype=np.float64),
        position_centers=np.asarray(cache["position_centers"], dtype=np.float64),
        position_widths=np.asarray(cache["position_widths"], dtype=np.float64),
        acceptance=np.asarray(cache["acceptance"], dtype=np.float64),
        sampled_indices=np.asarray(cache["sampled_indices"], dtype=np.int64),
        n_hit=np.array(cache["n_hit"]),
        n_sampled=np.array(cache["n_sampled"]),
    )


def _load_cache(cache_path):
    data = np.load(cache_path, allow_pickle=False)
    return {
        "input_sig": str(data["input_sig"]),
        "estimator": str(data["estimator"]),
        "total_hit_weight": float(data["total_hit_weight"]),
        "sample_weights": data["sample_weights"],
        "beta_gamma": data["beta_gamma"],
        "entry_d": data["entry_d"],
        "exit_d": data["exit_d"],
        "position_centers": data["position_centers"],
        "position_widths": data["position_widths"],
        "acceptance": data["acceptance"],
        "sampled_indices": data["sampled_indices"],
        "n_hit": int(data["n_hit"]),
        "n_sampled": int(data["n_sampled"]),
    }


def load_or_build_decay_cache(
    flavor,
    mass_GeV,
    mass_label,
    data,
    hits,
    entry_d,
    exit_d,
    mesh,
    mother_samples=DEFAULT_MOTHER_SAMPLES,
    position_bins=DEFAULT_POSITION_BINS,
    decay_templates=DEFAULT_DECAY_TEMPLATES,
    decays_per_bin=DEFAULT_DECAYS_PER_BIN,
    random_seed=DEFAULT_RANDOM_SEED,
    force=False,
):
    """
    Load or build the per-mass decay-acceptance cache.

    The active-path estimator is exact for small hit samples and a weighted
    resampling estimator once the hit-mother pool becomes too large.
    """
    hit_mask = hits & np.isfinite(entry_d) & np.isfinite(exit_d) & (exit_d > entry_d)
    cache_dir = DECAY_CACHE_DIR / flavor
    cache_path = cache_dir / f"decay_acceptance_{mass_label}.npz"
    input_sig = _cache_signature(
        data=data,
        hits=hit_mask,
        entry_d=entry_d,
        exit_d=exit_d,
        mother_samples=mother_samples,
        position_bins=position_bins,
        decay_templates=decay_templates,
        decays_per_bin=decays_per_bin,
        random_seed=random_seed,
    )

    if cache_path.exists() and not force:
        cached = _load_cache(cache_path)
        if cached["input_sig"] == input_sig:
            return cached

    rng = np.random.default_rng(int(random_seed) + int(round(mass_GeV * 1_000)))
    sampled = _select_sampled_mothers(data, hit_mask, mother_samples, rng)
    sampled_indices = sampled["indices"]
    n_sampled = len(sampled_indices)
    n_hit = int(hit_mask.sum())

    if sampled["estimator"] == "empty":
        cache = {
            "input_sig": input_sig,
            "estimator": sampled["estimator"],
            "total_hit_weight": 0.0,
            "sample_weights": np.empty(0, dtype=np.float64),
            "beta_gamma": np.empty(0, dtype=np.float64),
            "entry_d": np.empty(0, dtype=np.float64),
            "exit_d": np.empty(0, dtype=np.float64),
            "position_centers": np.empty((0, position_bins), dtype=np.float64),
            "position_widths": np.empty((0, position_bins), dtype=np.float64),
            "acceptance": np.empty((0, position_bins), dtype=np.float64),
            "sampled_indices": np.empty(0, dtype=np.int64),
            "n_hit": 0,
            "n_sampled": 0,
        }
        _save_cache(cache_path, cache)
        return cache

    sample_weights = sampled["weights"]
    sample_entry = entry_d[sampled_indices]
    sample_exit = exit_d[sampled_indices]
    sample_beta_gamma = data["beta_gamma"][sampled_indices]
    sample_parent_p4 = np.column_stack(
        [
            data["E"][sampled_indices],
            data["px"][sampled_indices],
            data["py"][sampled_indices],
            data["pz"][sampled_indices],
        ]
    )
    sample_momenta = np.column_stack(
        [data["px"][sampled_indices], data["py"][sampled_indices], data["pz"][sampled_indices]]
    )
    sample_directions, valid_mothers = momenta_to_directions(sample_momenta)

    edge_grid = np.linspace(0.0, 1.0, position_bins + 1, dtype=np.float64)[None, :]
    sample_paths = sample_exit - sample_entry
    position_edges = sample_entry[:, None] + sample_paths[:, None] * edge_grid
    position_centers = 0.5 * (position_edges[:, :-1] + position_edges[:, 1:])
    position_widths = position_edges[:, 1:] - position_edges[:, :-1]

    acceptance = np.zeros((n_sampled, position_bins), dtype=np.float64)
    n_decay_draws = max(int(decays_per_bin), 1)
    templates = sample_rest_frame_decays(
        mass_GeV=mass_GeV,
        n_templates=decay_templates,
        flavor=flavor,
        random_seed=int(random_seed) + int(round(mass_GeV * 10_000)),
    )
    if not templates:
        raise RuntimeError(f"FairShip produced no HNL decay templates at m_N={mass_GeV:.3f} GeV")

    for sample_idx in range(n_sampled):
        if not valid_mothers[sample_idx]:
            continue
        mother_direction = sample_directions[sample_idx]
        parent_p4 = sample_parent_p4[sample_idx]
        for bin_idx in range(position_bins):
            accepted = 0
            for template_idx in rng.integers(len(templates), size=n_decay_draws):
                accepted += int(
                    _accept_decay_at_position(
                        mesh=mesh,
                        decay_distance=position_centers[sample_idx, bin_idx],
                        exit_distance=sample_exit[sample_idx],
                        mother_direction=mother_direction,
                        parent_p4=parent_p4,
                        template=templates[int(template_idx)],
                    )
                )
            acceptance[sample_idx, bin_idx] = accepted / float(n_decay_draws)

    cache = {
        "input_sig": input_sig,
        "estimator": sampled["estimator"],
        "total_hit_weight": sampled["total_weight"],
        "sample_weights": sample_weights,
        "beta_gamma": sample_beta_gamma,
        "entry_d": sample_entry,
        "exit_d": sample_exit,
        "position_centers": position_centers,
        "position_widths": position_widths,
        "acceptance": acceptance,
        "sampled_indices": sampled_indices,
        "n_hit": n_hit,
        "n_sampled": n_sampled,
    }
    _save_cache(cache_path, cache)
    return cache
