import json
from types import SimpleNamespace

import pandas as pd

from alp_fermion import record_dense_structural as record
from alp_fermion.uncertainty_campaign import sha256_file


def test_dense_structural_manifest_is_portable_and_pins_inputs(tmp_path, monkeypatch):
    mass_grid = tmp_path / "mass_grid.csv"
    central = tmp_path / "central.csv"
    central_manifest = tmp_path / "MANIFEST.json"
    structural = tmp_path / "structural.csv"
    vector_dir = tmp_path / "vectors"
    geometry_dir = tmp_path / "geometry"
    template_dir = tmp_path / "templates"
    for directory in (vector_dir, geometry_dir, template_dir):
        directory.mkdir()
    frame = pd.DataFrame({
        "mass_GeV": [1.30],
        "has_sensitivity": [True],
        "peak_N": [30.0],
        "peak_invf": [2e-7],
        "invf_min": [2e-8],
        "invf_max": [2e-6],
        "invf_min_open": [False],
        "invf_max_open": [False],
    })
    frame[["mass_GeV"]].to_csv(mass_grid, index=False)
    frame.to_csv(central, index=False)
    frame.to_csv(structural, index=False)
    central_manifest.write_text(json.dumps({
        "csv_sha256": sha256_file(central),
        "producer_git_sha": "b" * 40,
        "source_run": "1200000-event central production pool",
        "grid": {"n_masses": 1},
    }))

    monkeypatch.setattr(
        record,
        "_validate_vectors",
        lambda *args, **kwargs: {
            "n_vector_files": 1,
            "total_vector_bytes": 100,
            "vector_tree_sha256": "vectors",
        },
    )
    monkeypatch.setattr(
        record,
        "_validate_geometry",
        lambda *args, **kwargs: {
            "n_geometry_files": 1,
            "total_geometry_bytes": 50,
            "geometry_tree_sha256": "geometry",
        },
    )
    monkeypatch.setattr(
        record,
        "_template_provenance",
        lambda *args, **kwargs: {
            "path": str(template_dir),
            "n_files": 1,
            "tree_sha256": "templates",
            "total_bytes": 75,
            "decay_model": "2310_structural",
            "n_templates_per_mass": 20_000,
        },
    )
    monkeypatch.setattr(
        record,
        "_validate_sensitivity_csv",
        lambda *args, **kwargs: {
            "sensitivity_csv_sha256": sha256_file(structural),
            "n_sensitivity_rows": 1,
            "n_sensitive_rows": 1,
        },
    )
    monkeypatch.setattr(record, "_require_commit", lambda value: "a" * 40)
    monkeypatch.setattr(
        record,
        "_git_state",
        lambda: {
            "commit": "c" * 40,
            "tracked_diff_sha256": "0" * 64,
            "tracked_tree_clean": True,
        },
    )
    args = SimpleNamespace(
        mass_grid=mass_grid,
        central_curve=central,
        central_manifest=central_manifest,
        structural_curve=structural,
        vector_dir=vector_dir,
        geometry_dir=geometry_dir,
        template_dir=template_dir,
        producer_git_sha="a" * 40,
        command="python -m alp_fermion.sensitivity --decay-samples 60",
        n_pool=1_200_000,
        decay_samples=60,
        production_seed=42,
        template_seed=5234,
    )

    manifest = record.build_manifest(args)

    assert manifest["inputs"]["central_curve"]["sha256"] == sha256_file(central)
    assert manifest["inputs"]["central_publication_manifest"][
        "producer_git_sha"
    ] == "b" * 40
    assert manifest["inputs"]["production_vectors"]["vector_tree_sha256"] == "vectors"
    assert manifest["output"]["sha256"] == sha256_file(structural)
    assert manifest["producer_code"]["commit"] == "a" * 40
    assert manifest["recorder_code"]["commit"] == "c" * 40
    assert str(tmp_path) not in str(manifest)
