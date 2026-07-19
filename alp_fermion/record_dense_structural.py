#!/usr/bin/env python3
"""Record portable provenance for a completed refined-grid 2310 comparison."""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from alp_fermion.run_uncertainty_campaign import (
    REPO_ROOT,
    _git_state,
    _load_mass_grid,
    _template_provenance,
    _validate_geometry,
    _validate_sensitivity_csv,
    _validate_vectors,
)
from alp_fermion.uncertainty_campaign import (
    atomic_json,
    sha256_file,
    structural_alternative_from_curves,
)


def _require_commit(commit: str) -> str:
    resolved = subprocess.check_output(
        ["git", "rev-parse", f"{commit}^{{commit}}"],
        cwd=REPO_ROOT,
        text=True,
    ).strip()
    return resolved


def build_manifest(args) -> dict:
    mass_grid = _load_mass_grid(args.mass_grid)
    masses = mass_grid["masses_GeV"]
    central = pd.read_csv(args.central_curve)
    structural = pd.read_csv(args.structural_curve)
    central_manifest = json.loads(args.central_manifest.read_text())
    central_sha = sha256_file(args.central_curve)
    if central_manifest.get("csv_sha256") != central_sha:
        raise ValueError("canonical central curve/manifest checksum mismatch")
    if int(central_manifest.get("grid", {}).get("n_masses", -1)) != len(masses):
        raise ValueError("canonical central manifest uses a different mass grid")
    source_run = str(central_manifest.get("source_run", ""))
    if str(int(args.n_pool)) not in source_run.replace(",", ""):
        raise ValueError("canonical central manifest does not confirm --n-pool")
    comparison = structural_alternative_from_curves(central, structural)
    if list(comparison["mass_GeV"].astype(float)) != masses:
        raise ValueError("comparison curves do not exactly match --mass-grid order")

    vectors = _validate_vectors(args.vector_dir, masses, hash_outputs=True)
    geometry = _validate_geometry(args.geometry_dir, masses)
    templates = _template_provenance(
        args.template_dir,
        expected_decay_model="2310_structural",
        masses=masses,
    )
    sensitivity = _validate_sensitivity_csv(args.structural_curve, masses)
    producer_commit = _require_commit(args.producer_git_sha)

    central_sensitive = central["has_sensitivity"].fillna(False).astype(bool)
    structural_sensitive = structural["has_sensitivity"].fillna(False).astype(bool)
    return {
        "artifact": "GRENDEL BC10 refined-grid exact 2310 structural comparison",
        "schema_version": 1,
        "role": "one_sided_structural_model_comparison_outside_pointwise_halo",
        "producer_code": {
            "commit": producer_commit,
            "repository": "llpatcolliders (branch bc10-alp)",
        },
        "recorder_code": _git_state(),
        "command": args.command,
        "analysis": {
            "n_parent_pool": int(args.n_pool),
            "n_decay_samples_per_detector_entering_alp": int(args.decay_samples),
            "production_seed": int(args.production_seed),
            "template_seed": int(args.template_seed),
            "production_reuse": "exact canonical high-statistics ALP four-vectors",
            "geometry_reuse": "exact canonical ALP entry/exit cache",
            "decay_model": "exact SensCalc arXiv:2310.03524 structural alternative",
        },
        "inputs": {
            "central_curve": {
                "file": args.central_curve.name,
                "sha256": central_sha,
                "rows": len(central),
            },
            "central_publication_manifest": {
                "file": args.central_manifest.name,
                "sha256": sha256_file(args.central_manifest),
                "producer_git_sha": central_manifest.get("producer_git_sha"),
            },
            "mass_grid": {
                "file": args.mass_grid.name,
                "sha256": sha256_file(args.mass_grid),
                "canonical_sha256": mass_grid["canonical_sha256"],
                "n_masses": len(masses),
            },
            "production_vectors": vectors,
            "geometry_cache": geometry,
            "decay_templates": {
                key: value for key, value in templates.items() if key != "path"
            },
        },
        "output": {
            "file": args.structural_curve.name,
            "sha256": sensitivity["sensitivity_csv_sha256"],
            "rows": sensitivity["n_sensitivity_rows"],
            "n_sensitive": sensitivity["n_sensitive_rows"],
        },
        "topology_comparison": {
            "n_central_sensitive": int(central_sensitive.sum()),
            "n_structural_sensitive": int(structural_sensitive.sum()),
            "n_restored": int(comparison["restores_sensitivity"].sum()),
            "n_removed": int(comparison["removes_sensitivity"].sum()),
            "n_topology_differences": int(comparison["topology_differs"].sum()),
        },
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "platform": platform.platform(),
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--central-curve", type=Path, required=True)
    parser.add_argument("--central-manifest", type=Path, required=True)
    parser.add_argument("--structural-curve", type=Path, required=True)
    parser.add_argument("--mass-grid", type=Path, required=True)
    parser.add_argument("--vector-dir", type=Path, required=True)
    parser.add_argument("--geometry-dir", type=Path, required=True)
    parser.add_argument("--template-dir", type=Path, required=True)
    parser.add_argument("--producer-git-sha", required=True)
    parser.add_argument("--command", required=True)
    parser.add_argument("--n-pool", type=int, default=1_200_000)
    parser.add_argument("--decay-samples", type=int, default=60)
    parser.add_argument("--production-seed", type=int, default=42)
    parser.add_argument("--template-seed", type=int, default=5234)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = build_manifest(args)
    atomic_json(args.output, manifest)
    print(json.dumps(manifest["topology_comparison"], sort_keys=True))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
