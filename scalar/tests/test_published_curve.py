import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scalar.production import MASS_GRID


PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"
BUNDLE = PUBLISHED / "bundle"


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_published_curve_and_manifest_are_self_consistent():
    curve_path = PUBLISHED / "bc4_island.csv"
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    frame = pd.read_csv(curve_path)

    assert _sha256(curve_path) == manifest["csv_sha256"]
    assert not frame["mass_GeV"].duplicated().any()
    assert frame["mass_GeV"].is_monotonic_increasing
    assert np.array_equal(
        frame["mass_GeV"].to_numpy(float), np.asarray(MASS_GRID, dtype=float)
    )

    headline = manifest["headline_reach"]
    sensitive = frame[frame["has_sensitivity"].astype(bool)]
    best = sensitive.loc[sensitive["u2_min"].idxmin()]
    assert len(frame) == headline["n_grid_points"]
    assert len(sensitive) == headline["n_sensitive_points"]
    assert best["mass_GeV"] == pytest.approx(
        headline["best_sin2theta_min_at_mass_GeV"]
    )
    assert best["u2_min"] == pytest.approx(
        headline["best_sin2theta_min"], rel=1e-12
    )
    assert sensitive["mass_GeV"].min() == pytest.approx(
        headline["mass_lo_grid_GeV"]
    )
    assert sensitive["mass_GeV"].max() == pytest.approx(
        headline["mass_hi_grid_GeV"]
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

    upper_open = frame[frame["u2_max_open"].astype(bool)]
    assert [upper_open["mass_GeV"].min(), upper_open["mass_GeV"].max()] == (
        manifest["topology"]["upper_open_grid_run_GeV"]
    )


def test_uncertainty_bundle_is_complete_and_hash_linked():
    manifest = json.loads((BUNDLE / "UNCERTAINTY_MANIFEST.json").read_text())
    canonical_manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())

    for name, metadata in manifest["outputs"].items():
        path = PUBLISHED / name
        if not path.is_file():
            path = BUNDLE / name
        assert path.is_file()
        assert _sha256(path) == metadata["sha256"]
        if "rows" in metadata:
            assert len(pd.read_csv(path)) == metadata["rows"]

    assert manifest["outputs"]["bc4_island.csv"]["sha256"] == (
        canonical_manifest["csv_sha256"]
    )
    assert manifest["collector"]["uncertainty_band_sha256"] == _sha256(
        PUBLISHED.parents[1] / "uncertainty_band.py"
    )

    raw = pd.read_csv(BUNDLE / "bc4_uncertainty_variations.csv")
    band = pd.read_csv(BUNDLE / "bc4_single_source_variation_envelope.csv")
    expected_axes = {
        "central": 1,
        "scale": 6,
        "pdf": 100,
        "mass": 2,
        "decay_model": 1,
        "numerical_control": 2,
    }
    assert raw.groupby("axis")["variation"].nunique().to_dict() == expected_axes
    assert raw.groupby("variation").size().eq(len(MASS_GRID)).all()
    assert len(raw) == sum(expected_axes.values()) * len(MASS_GRID)
    assert np.array_equal(
        band["mass_GeV"].to_numpy(float), np.asarray(MASS_GRID, dtype=float)
    )
    assert set(band["envelope_definition"]) == {
        "single_source_variation_envelope"
    }
    insensitive = ~band["has_sensitivity"].astype(bool)
    assert band.loc[
        insensitive,
        [
            "u2_min_envelope_lo",
            "u2_min_envelope_hi",
            "u2_max_envelope_lo",
            "u2_max_envelope_hi",
        ],
    ].isna().all().all()

    numerical = manifest["numerical_control_summary"]
    assert numerical["topology_vs_campaign"]["difference_masses_GeV"] == []
    assert numerical["topology_vs_canonical"]["difference_masses_GeV"] == []
    assert numerical["u2_min"]["masses_not_subdominant_GeV"] == []
    # Lower edge (the headline reach) is fully physics-dominated. On the noisier
    # short-lifetime upper edge, the 100k-parent-per-species Lambda_b bundle
    # leaves the numerical control not-subdominant at two masses (2.0, 2.3 GeV).
    # This is the accepted statistics level of the published bundle (not a runtime
    # or disk limit): there the same-physics repeat spread is comparable to the
    # physical variation, so those two upper-edge points are diagnostic-only and
    # must not be read as converged. The lower edge and deepest reach are
    # unaffected; a higher-statistics bundle would resolve them if the upper edge
    # ever entered the headline.
    assert numerical["u2_max"]["masses_not_subdominant_GeV"] == [2.0, 2.3]
    physical = manifest["physical_variation_topology_summary"]
    assert physical["topology_vs_campaign"]["difference_masses_GeV"] == [3.8]
    assert physical["topology_vs_canonical"]["difference_masses_GeV"] == [3.8]
