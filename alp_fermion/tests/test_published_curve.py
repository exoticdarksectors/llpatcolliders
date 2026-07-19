import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from alp_fermion.mass_grid import ALP_MASS_GRID


PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"
BUNDLE = PUBLISHED / "bundle"


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sensitive_runs(frame):
    runs = []
    start = None
    masses = frame["mass_GeV"].to_numpy(float)
    flags = frame["has_sensitivity"].astype(bool).to_numpy()
    for index, sensitive in enumerate(flags):
        if sensitive and start is None:
            start = index
        if start is not None and (not sensitive or index == len(flags) - 1):
            end = index if sensitive and index == len(flags) - 1 else index - 1
            runs.append([masses[start], masses[end]])
            start = None
    return runs


def test_published_curve_and_manifest_are_self_consistent():
    curve_path = PUBLISHED / "bc10_sensitivity.csv"
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    frame = pd.read_csv(curve_path)

    assert _sha256(curve_path) == manifest["csv_sha256"]
    assert not frame["mass_GeV"].duplicated().any()
    assert frame["mass_GeV"].is_monotonic_increasing
    assert np.array_equal(
        frame["mass_GeV"].to_numpy(float), np.asarray(ALP_MASS_GRID, dtype=float)
    )

    headline = manifest["headline_reach"]
    sensitive = frame[frame["has_sensitivity"].astype(bool)]
    best = sensitive.loc[sensitive["invf_min"].idxmin()]
    assert len(frame) == headline["n_grid_points"]
    assert len(sensitive) == headline["n_sensitive_points"]
    assert best["mass_GeV"] == pytest.approx(
        headline["best_invf_min_at_mass_GeV"]
    )
    assert best["invf_min"] == pytest.approx(
        headline["best_invf_min_GeV_inv"], rel=1e-12
    )

    last = sensitive.index[-1]
    following = last + 1
    fraction = (
        np.log(3.0) - np.log(frame.loc[last, "peak_N"])
    ) / (
        np.log(frame.loc[following, "peak_N"])
        - np.log(frame.loc[last, "peak_N"])
    )
    closure = frame.loc[last, "mass_GeV"] + fraction * (
        frame.loc[following, "mass_GeV"] - frame.loc[last, "mass_GeV"]
    )
    assert closure == pytest.approx(headline["mass_hi_closure_GeV"], rel=1e-12)
    assert _sensitive_runs(frame) == manifest["topology"]["sensitive_grid_runs_GeV"]


def test_uncertainty_bundle_is_complete_and_hash_linked():
    manifest = json.loads((BUNDLE / "UNCERTAINTY_MANIFEST.json").read_text())
    canonical_manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())

    for name, expected_sha256 in manifest["outputs"].items():
        path = BUNDLE / name
        assert path.is_file()
        assert _sha256(path) == expected_sha256

    assert manifest["inputs"]["canonical_central_curve"]["sha256"] == (
        canonical_manifest["csv_sha256"]
    )
    assert manifest["variation_counts"] == {
        "central": 1,
        "scale": 6,
        "pdf": 100,
        "mb": 2,
        "decay_gg": 3,
        "cbs": 2,
        "decay_structure": 1,
        "numerical_control": 2,
        "pointwise_halo_total": 114,
        "physics_plus_structural_total": 115,
        "total_with_controls": 117,
    }

    raw = pd.read_csv(BUNDLE / "bc10_uncertainty_variations.csv")
    band = pd.read_csv(BUNDLE / "bc10_single_source_variation_envelope.csv")
    dense = pd.read_csv(
        BUNDLE / "bc10_decay_2310_structural_alternative_dense.csv"
    )
    dense_curve = BUNDLE / "bc10_decay_2310_structural_curve_dense.csv"
    expected_axes = {
        "central": 1,
        "scale": 6,
        "pdf": 100,
        "mb": 2,
        "decay_gg": 3,
        "cbs": 2,
        "decay_structure": 1,
        "numerical_control": 2,
    }
    assert raw.groupby("axis")["variation"].nunique().to_dict() == expected_axes
    assert raw.groupby("variation").size().eq(99).all()
    assert len(raw) == 117 * 99
    assert len(band) == 99
    assert set(band["envelope_definition"]) == {
        "single_source_variation_envelope"
    }
    insensitive = ~band["has_sensitivity"].astype(bool)
    assert band.loc[
        insensitive,
        [
            "invf_min_envelope_lo",
            "invf_min_envelope_hi",
            "invf_max_envelope_lo",
            "invf_max_envelope_hi",
        ],
    ].isna().all().all()

    assert np.array_equal(
        dense["mass_GeV"].to_numpy(float), np.asarray(ALP_MASS_GRID, dtype=float)
    )
    assert _sha256(dense_curve) == manifest["inputs"][
        "dense_structural_curve"
    ]["sha256"]
    assert dense.loc[dense["restores_sensitivity"], "mass_GeV"].tolist() == [
        *(value / 100 for value in range(126, 140)),
        *(value / 100 for value in range(142, 147)),
    ]
    assert not dense["removes_sensitivity"].any()

    headline = manifest["headline"]
    assert headline["numerical_control_topology"] == {
        "n_differences": 1,
        "difference_masses_GeV": [3.3],
        "difference_variations_by_mass": {"3.3": "central_repeat_2"},
    }
    assert headline["dense_decay_structure"] == {
        "n_mass_points": 145,
        "n_sensitive": 122,
        "n_restored": 19,
        "n_removed": 0,
        "n_topology_differences": 19,
    }
