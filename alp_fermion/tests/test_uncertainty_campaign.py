import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import uncertainty_campaign as campaign_definitions  # noqa: E402
import combine_uncertainty_band as campaign_collector  # noqa: E402
from uncertainty_campaign import (  # noqa: E402
    all_variations,
    combine_band,
    decay_structure_variation,
    numerical_control_variations,
    stable_seed,
    structural_alternative_from_curves,
    structural_alternative_table,
)
import run_uncertainty_campaign as campaign_runner  # noqa: E402
from run_uncertainty_campaign import _compact_completed_run  # noqa: E402


def test_stable_seed_is_reproducible_and_member_specific():
    assert stable_seed("pdf_0001") == stable_seed("pdf_0001")
    assert stable_seed("pdf_0001") != stable_seed("pdf_0002")
    assert 0 < stable_seed("scale_muR2_muF2") < 2**32


def test_numerical_controls_have_independent_production_and_reco_seeds(tmp_path):
    grid = tmp_path / "central.dat"
    grid.write_text("grid")
    controls = numerical_control_variations(grid)
    assert [item["name"] for item in controls] == [
        "central_repeat_1", "central_repeat_2"
    ]
    assert len({item["production_seed"] for item in controls}) == 2
    assert len({item["reco_seed_offset"] for item in controls}) == 2
    assert all(
        item["production_seed"] != item["reco_seed_offset"]
        for item in controls
    )


def test_decay_structure_uses_exact_central_vectors_and_named_templates(tmp_path):
    grid = tmp_path / "central.dat"
    grid.write_text("grid")
    variation = decay_structure_variation(grid)
    assert variation["name"] == "decay_2310_structural"
    assert variation["axis"] == "decay_structure"
    assert variation["decay_model"] == "2310_structural"
    assert variation["template_variant"] == "decay_2310_structural"
    assert variation["production_mode"] == "central_vectors_exact_reuse"
    assert variation["combination_role"].endswith("outside_halo")


def test_campaign_counts_keep_structural_model_outside_halo(tmp_path, monkeypatch):
    grid = tmp_path / "central.dat"
    grid.write_text("grid")
    fonll = [{
        "name": "central",
        "axis": "central",
        "campaign_axis": "fonll",
        "grid_path": str(grid),
    }]
    for axis, count in (("scale", 6), ("pdf", 100), ("mb", 2)):
        fonll.extend({
            "name": f"{axis}_{index}",
            "axis": axis,
            "campaign_axis": "fonll",
            "grid_path": str(grid),
        } for index in range(count))
    monkeypatch.setattr(
        campaign_definitions, "discover_fonll_variations", lambda _: fonll
    )
    variations = all_variations(tmp_path)
    counts = {
        axis: sum(item["axis"] == axis for item in variations)
        for axis in {
            "central", "scale", "pdf", "mb", "decay_gg", "cbs",
            "decay_structure", "numerical_control",
        }
    }
    assert counts == {
        "central": 1,
        "scale": 6,
        "pdf": 100,
        "mb": 2,
        "decay_gg": 3,
        "cbs": 2,
        "decay_structure": 1,
        "numerical_control": 2,
    }
    structural = next(v for v in variations if v["axis"] == "decay_structure")
    assert structural["combination_role"].endswith("outside_halo")


def test_structural_template_provenance_requires_20k_and_model_tag(
    tmp_path, monkeypatch
):
    paths = [tmp_path / "templates_1.npz", tmp_path / "templates_2.npz"]
    for path in paths:
        np.savez_compressed(
            path,
            decay_model=np.array("2310_structural"),
            n_templates=np.int32(20_000),
        )
    monkeypatch.setattr(
        campaign_runner,
        "_expected_template_paths",
        lambda directory, masses: paths,
    )
    result = campaign_runner._template_provenance(
        tmp_path, "2310_structural"
    )
    assert result["n_files"] == 2
    assert result["n_templates_per_mass"] == 20_000
    assert result["decay_model"] == "2310_structural"

    np.savez_compressed(
        paths[0],
        decay_model=np.array("2501"),
        n_templates=np.int32(20_000),
    )
    with pytest.raises(RuntimeError, match="do not match"):
        campaign_runner._template_provenance(tmp_path, "2310_structural")


def test_structural_reference_hashes_exact_central_vectors(tmp_path, monkeypatch):
    mass_grid = {
        "canonical_sha256": "mass-hash",
        "masses_GeV": [1.0],
    }
    marker_path = tmp_path / "production.complete.json"
    marker_path.write_text(json.dumps({
        "variation": {
            "name": "central",
            "grid_sha256": "grid-hash",
            "cbs_amplitude_scale": 1.0,
        },
        "n_pool": 600_000,
        "mass_grid": mass_grid,
        "vector_tree_sha256": "vector-hash",
    }))
    calls = []

    def validate(directory, masses, hash_outputs=True):
        assert masses == [1.0]
        calls.append(hash_outputs)
        return {
            "n_vector_files": 97,
            "total_vector_bytes": 123,
            "vector_tree_sha256": "vector-hash" if hash_outputs else None,
        }

    monkeypatch.setattr(campaign_runner, "_validate_vectors", validate)
    _, info = campaign_runner._validate_central_production_reference(
        marker_path,
        tmp_path / "vectors",
        {"grid_sha256": "grid-hash"},
        600_000,
        hash_vectors=True,
        mass_grid=mass_grid,
    )
    assert calls == [True]
    assert info["vector_tree_sha256"] == "vector-hash"

    with pytest.raises(RuntimeError, match="different FONLL grid"):
        campaign_runner._validate_central_production_reference(
            marker_path,
            tmp_path / "vectors",
            {"grid_sha256": "other-grid"},
            600_000,
            hash_vectors=True,
            mass_grid=mass_grid,
        )


def test_structural_reference_strictly_validates_legacy_central_grid(
    tmp_path, monkeypatch
):
    central_dir = tmp_path / "central"
    marker_path = central_dir / "production.complete.json"
    marker_path.parent.mkdir(parents=True)
    marker_path.write_text(json.dumps({
        "variation": {
            "name": "central",
            "grid_sha256": "grid-hash",
            "cbs_amplitude_scale": 1.0,
        },
        "n_pool": 600_000,
        "n_vector_files": 1,
        "vector_tree_sha256": "vector-hash",
    }))
    mass_grid = {
        "canonical_sha256": "mass-hash",
        "masses_GeV": [1.0],
    }
    monkeypatch.setattr(
        campaign_runner,
        "_validate_vectors",
        lambda *args, **kwargs: {
            "n_vector_files": 1,
            "total_vector_bytes": 123,
            "vector_tree_sha256": "vector-hash",
        },
    )
    calls = []

    def validate_curve(path, masses):
        calls.append((path, masses))
        return {
            "sensitivity_csv_sha256": "curve-hash",
            "n_sensitivity_rows": 1,
            "n_sensitive_rows": 1,
        }

    monkeypatch.setattr(
        campaign_runner, "_validate_sensitivity_csv", validate_curve
    )
    _, info = campaign_runner._validate_central_production_reference(
        marker_path,
        central_dir / "llp_4vectors",
        {"grid_sha256": "grid-hash"},
        600_000,
        hash_vectors=True,
        mass_grid=mass_grid,
    )

    assert calls == [
        (central_dir / "analysis" / "bc10_sensitivity.csv", [1.0])
    ]
    assert info["mass_grid_validation"] == {
        "mode": "legacy_vector_tree_and_sensitivity_csv",
        "central_sensitivity_csv": str(
            central_dir / "analysis" / "bc10_sensitivity.csv"
        ),
        "sensitivity_csv_sha256": "curve-hash",
        "n_sensitivity_rows": 1,
        "n_sensitive_rows": 1,
    }

    marker = json.loads(marker_path.read_text())
    marker["n_vector_files"] = 2
    marker_path.write_text(json.dumps(marker))
    with pytest.raises(RuntimeError, match="legacy central.*different mass grid"):
        campaign_runner._validate_central_production_reference(
            marker_path,
            central_dir / "llp_4vectors",
            {"grid_sha256": "grid-hash"},
            600_000,
            hash_vectors=True,
            mass_grid=mass_grid,
        )


def test_mass_grid_file_is_hashed_and_requires_strict_order(tmp_path):
    path = tmp_path / "dense_grid.csv"
    pd.DataFrame({"mass_GeV": [1.18, 1.19, 1.20]}).to_csv(path, index=False)

    result = campaign_runner._load_mass_grid(path)

    assert result["n_masses"] == 3
    assert result["masses_GeV"] == [1.18, 1.19, 1.20]
    assert result["source_path"] == str(path.resolve())
    assert result["source_sha256"] == campaign_definitions.sha256_file(path)
    assert len(result["canonical_sha256"]) == 64

    pd.DataFrame({"mass_GeV": [1.18, 1.20, 1.19]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match="strictly increasing"):
        campaign_runner._load_mass_grid(path)


def test_collector_records_separate_consistent_structural_code_state():
    halo_code = {"commit": "halo", "tracked_tree_clean": True}
    structural_code = {"commit": "structural", "tracked_tree_clean": True}
    registry = [
        {"variation": {"axis": "central"}, "code": halo_code},
        {"variation": {"axis": "scale"}, "code": halo_code},
        {
            "variation": {"axis": "decay_structure"},
            "code": structural_code,
        },
    ]
    observed_halo, observed_structural = (
        campaign_collector._campaign_code_states(registry)
    )
    assert observed_halo == halo_code
    assert observed_structural == structural_code

    registry[1]["code"] = {"commit": "mixed"}
    with pytest.raises(ValueError, match="pointwise-halo"):
        campaign_collector._campaign_code_states(registry)


def test_published_registry_paths_are_portable(tmp_path):
    scratch = tmp_path / "scratch"
    artifact = scratch / "runs" / "central" / "variation.complete.json"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("{}")
    variation = {
        "axis": "central",
        "grid_path": str(tmp_path / "fonll" / "central.dat"),
        "template_variant": "central",
    }
    template = {
        "path": str(tmp_path / "templates" / "central"),
        "tree_sha256": "template-hash",
    }

    payload = {
        "artifact": campaign_collector._scratch_relative(artifact, scratch),
        "variation": campaign_collector._portable_variation(variation),
        "template": campaign_collector._portable_template(template, variation),
    }

    assert payload["artifact"] == "runs/central/variation.complete.json"
    assert payload["variation"]["grid_file"] == "central.dat"
    assert payload["template"]["path_role"] == "central"
    assert str(tmp_path) not in json.dumps(payload)
    with pytest.raises(ValueError, match="outside scratch root"):
        campaign_collector._scratch_relative(tmp_path / "elsewhere", scratch)


def test_exact_campaign_builds_single_source_variation_envelope():
    xc = -8.0
    variations = [("central", "central", 0.0)]
    variations += [
        (f"scale_{index}", "scale", shift)
        for index, shift in enumerate((-0.10, -0.05, -0.01, 0.01, 0.05, 0.10))
    ]
    pdf_shifts = np.r_[np.linspace(-0.02, 0.02, 99), 1.0]
    variations += [
        (f"pdf_{index:04d}", "pdf", shift)
        for index, shift in enumerate(pdf_shifts, start=1)
    ]
    variations += [("mb_dn", "mb", -0.03), ("mb_up", "mb", 0.03)]
    variations += [
        ("gg_u", "decay_gg", -0.04),
        ("gg_d", "decay_gg", 0.0),
        ("gg_s", "decay_gg", 0.04),
        ("cbs_down", "cbs", -0.05),
        ("cbs_up", "cbs", 0.05),
        ("decay_2310_structural", "decay_structure", -0.50),
        ("central_repeat_1", "numerical_control", 0.12),
        ("central_repeat_2", "numerical_control", -0.007),
    ]
    rows = []
    for name, axis, shift in variations:
        rows.append({
            "mass_GeV": 1.0,
            "variation": name,
            "axis": axis,
            "has_sensitivity": True,
            "invf_min": 10.0 ** (xc + shift),
            "invf_max": 10.0 ** (xc + 1.0 + shift),
            "invf_min_open": False,
            "invf_max_open": False,
        })
    band = combine_band(pd.DataFrame(rows)).iloc[0]
    p16, p84 = np.quantile(pdf_shifts, [0.16, 0.84])
    assert band["envelope_definition"] == "single_source_variation_envelope"
    assert np.isclose(band["invf_min_pdf_p16"], 10.0 ** (xc + p16))
    assert np.isclose(band["invf_min_pdf_p84"], 10.0 ** (xc + p84))
    assert np.isclose(
        band["invf_min_pdf_log10_std"], np.std(pdf_shifts, ddof=1)
    )
    assert np.isclose(band["invf_min_envelope_lo"], 10.0 ** (xc - 0.10))
    assert np.isclose(band["invf_min_envelope_hi"], 10.0 ** (xc + 0.10))
    assert band["invf_min_envelope_lo_source"] == "scale"
    assert band["invf_min_envelope_hi_source"] == "scale"
    assert np.isclose(band["invf_min_decay_structure"], 10.0 ** (xc - 0.50))
    assert not bool(band["decay_structure_included_in_halo"])
    assert band["invf_min_pdf_n_finite"] == 100
    assert np.isclose(band["invf_min_repeat_max_abs_dex"], 0.12)
    assert bool(band["invf_min_repeat_not_subdominant"])
    assert not bool(band["invf_min_variation_missing"])

    reference = pd.DataFrame([{
        "mass_GeV": 1.0,
        "has_sensitivity": True,
        "invf_min": 2.0 * 10.0**xc,
        "invf_max": 2.0 * 10.0**(xc + 1.0),
        "invf_min_open": False,
        "invf_max_open": False,
    }])
    rebased = combine_band(pd.DataFrame(rows), reference).iloc[0]
    assert rebased["invf_min_campaign_central"] == pytest.approx(10.0**xc)
    assert rebased["invf_min_central"] == pytest.approx(2.0 * 10.0**xc)
    assert rebased["invf_min_envelope_lo"] == pytest.approx(
        2.0 * 10.0**(xc - 0.10)
    )
    assert rebased["invf_min_envelope_hi"] == pytest.approx(
        2.0 * 10.0**(xc + 0.10)
    )
    assert rebased["invf_min_decay_structure"] == pytest.approx(10.0**(xc - 0.50))
    assert rebased["envelope_reference"] == "canonical_high_statistics_central"
    assert bool(rebased["canonical_rebase_topology_compatible"])


def test_missing_variation_boundary_is_explicit():
    rows = []
    definitions = [("central", "central")]
    definitions += [(f"scale_{i}", "scale") for i in range(6)]
    definitions += [(f"pdf_{i}", "pdf") for i in range(100)]
    definitions += [("mb_dn", "mb"), ("mb_up", "mb")]
    definitions += [(f"gg_{q}", "decay_gg") for q in "uds"]
    definitions += [("cbs_down", "cbs"), ("cbs_up", "cbs")]
    definitions += [("decay_2310_structural", "decay_structure")]
    definitions += [
        ("central_repeat_1", "numerical_control"),
        ("central_repeat_2", "numerical_control"),
    ]
    for name, axis in definitions:
        sensitive = name != "scale_0"
        rows.append({
            "mass_GeV": 1.0,
            "variation": name,
            "axis": axis,
            "has_sensitivity": sensitive,
            "invf_min": 1e-8 if sensitive else np.nan,
            "invf_max": 1e-7 if sensitive else np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
        })
    band = combine_band(pd.DataFrame(rows)).iloc[0]
    assert bool(band["invf_min_variation_missing"])
    assert bool(band["invf_max_variation_missing"])
    assert np.isfinite(band["invf_min_envelope_lo"])


def test_central_insensitive_row_remains_an_explicit_gap():
    rows = []
    definitions = [("central", "central")]
    definitions += [(f"scale_{i}", "scale") for i in range(6)]
    definitions += [(f"pdf_{i}", "pdf") for i in range(100)]
    definitions += [("mb_dn", "mb"), ("mb_up", "mb")]
    definitions += [(f"gg_{q}", "decay_gg") for q in "uds"]
    definitions += [("cbs_down", "cbs"), ("cbs_up", "cbs")]
    definitions += [("decay_2310_structural", "decay_structure")]
    definitions += [
        ("central_repeat_1", "numerical_control"),
        ("central_repeat_2", "numerical_control"),
    ]
    for name, axis in definitions:
        sensitive = name == "decay_2310_structural"
        rows.append({
            "mass_GeV": 0.54,
            "variation": name,
            "axis": axis,
            "has_sensitivity": sensitive,
            "invf_min": 1e-8 if sensitive else np.nan,
            "invf_max": 1e-7 if sensitive else np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
        })
    band = combine_band(pd.DataFrame(rows)).iloc[0]
    assert not bool(band["has_sensitivity"])
    assert not bool(band["any_variation_sensitive"])
    assert not bool(band["any_halo_variation_sensitive"])
    assert bool(band["decay_structure_has_sensitivity"])
    assert bool(band["decay_structure_restores_sensitivity"])
    assert bool(band["decay_structure_topology_differs"])
    assert np.isnan(band["invf_min_envelope_lo"])
    assert np.isnan(band["invf_max_envelope_hi"])


def test_numerical_repeat_topology_difference_is_preserved_at_central_gap():
    definitions = [("central", "central")]
    definitions += [(f"scale_{i}", "scale") for i in range(6)]
    definitions += [(f"pdf_{i}", "pdf") for i in range(100)]
    definitions += [("mb_dn", "mb"), ("mb_up", "mb")]
    definitions += [(f"gg_{q}", "decay_gg") for q in "uds"]
    definitions += [("cbs_down", "cbs"), ("cbs_up", "cbs")]
    definitions += [("decay_2310_structural", "decay_structure")]
    definitions += [
        ("central_repeat_1", "numerical_control"),
        ("central_repeat_2", "numerical_control"),
    ]
    rows = []
    for name, axis in definitions:
        sensitive = name == "central_repeat_2"
        rows.append({
            "mass_GeV": 3.30,
            "variation": name,
            "axis": axis,
            "has_sensitivity": sensitive,
            "invf_min": 1e-8 if sensitive else np.nan,
            "invf_max": 1e-7 if sensitive else np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
        })

    band_frame = combine_band(pd.DataFrame(rows))
    band = band_frame.iloc[0]
    assert not bool(band["has_sensitivity"])
    assert not bool(band["any_halo_variation_sensitive"])
    assert bool(band["numerical_control_any_sensitive"])
    assert not bool(band["numerical_control_all_sensitive"])
    assert band["numerical_control_n_sensitive"] == 1
    assert band["numerical_control_sensitive_variations"] == "central_repeat_2"
    assert bool(band["numerical_control_topology_differs"])
    assert (
        band["numerical_control_topology_difference_variations"]
        == "central_repeat_2"
    )
    assert np.isnan(band["invf_min_central_repeat_1"])
    assert band["invf_min_central_repeat_2"] == pytest.approx(1e-8)
    assert not bool(band["numerical_control_included_in_halo"])

    topology = campaign_collector._headline(band_frame)[
        "numerical_control_topology"
    ]
    assert topology == {
        "n_differences": 1,
        "difference_masses_GeV": [3.3],
        "difference_variations_by_mass": {"3.3": "central_repeat_2"},
    }


def test_physical_variation_topology_difference_is_named_at_central_gap():
    definitions = [("central", "central")]
    definitions += [(f"scale_{i}", "scale") for i in range(6)]
    definitions += [(f"pdf_{i}", "pdf") for i in range(100)]
    definitions += [("mb_dn", "mb"), ("mb_up", "mb")]
    definitions += [(f"gg_{q}", "decay_gg") for q in "uds"]
    definitions += [("cbs_down", "cbs"), ("cbs_up", "cbs")]
    definitions += [("decay_2310_structural", "decay_structure")]
    definitions += [
        ("central_repeat_1", "numerical_control"),
        ("central_repeat_2", "numerical_control"),
    ]
    rows = []
    for name, axis in definitions:
        sensitive = name == "gg_u"
        rows.append({
            "mass_GeV": 3.30,
            "variation": name,
            "axis": axis,
            "has_sensitivity": sensitive,
            "invf_min": 1e-8 if sensitive else np.nan,
            "invf_max": 1e-7 if sensitive else np.nan,
            "invf_min_open": False,
            "invf_max_open": False,
        })

    band_frame = combine_band(pd.DataFrame(rows))
    band = band_frame.iloc[0]
    assert not bool(band["has_sensitivity"])
    assert bool(band["any_halo_variation_sensitive"])
    assert bool(band["halo_restores_sensitivity"])
    assert band["halo_restores_sensitivity_variations"] == "gg_u"
    assert not bool(band["halo_removes_sensitivity"])
    assert bool(band["halo_topology_differs"])
    assert band["halo_topology_difference_variations"] == "gg_u"
    assert np.isnan(band["invf_min_envelope_lo"])

    topology = campaign_collector._headline(band_frame)[
        "physical_variation_topology"
    ]
    assert topology == {
        "n_differences": 1,
        "difference_masses_GeV": [3.3],
        "difference_variations_by_mass": {"3.3": "gg_u"},
        "restored_masses_GeV": [3.3],
        "removed_masses_GeV": [],
    }


def test_structural_table_publishes_restored_heavy_pseudoscalar_topology():
    rows = []
    for mass in (1.30, 1.35, 1.45):
        rows.extend([
            {
                "variation": "central",
                "axis": "central",
                "mass_GeV": mass,
                "has_sensitivity": False,
                "peak_N": 1.0,
                "peak_invf": 1e-6,
                "invf_min": np.nan,
                "invf_max": np.nan,
                "invf_min_open": False,
                "invf_max_open": False,
            },
            {
                "variation": "decay_2310_structural",
                "axis": "decay_structure",
                "mass_GeV": mass,
                "has_sensitivity": True,
                "peak_N": 5.0,
                "peak_invf": 2e-6,
                "invf_min": 1e-7,
                "invf_max": 1e-5,
                "invf_min_open": False,
                "invf_max_open": False,
            },
        ])
    result = structural_alternative_table(pd.DataFrame(rows))
    assert list(result["mass_GeV"]) == [1.30, 1.35, 1.45]
    assert result["restores_sensitivity"].all()
    assert result["topology_differs"].all()
    assert not result["included_in_pointwise_halo"].any()
    assert (result["invf_min"] == 1e-7).all()
    assert (result["invf_max"] == 1e-5).all()


def test_dense_structural_table_requires_and_preserves_exact_mass_grid():
    central = pd.DataFrame({
        "mass_GeV": [1.30, 1.31],
        "has_sensitivity": [False, True],
        "peak_N": [0.5, 4.0],
        "peak_invf": [1e-7, 2e-7],
        "invf_min": [np.nan, 1e-8],
        "invf_max": [np.nan, 1e-6],
        "invf_min_open": [False, False],
        "invf_max_open": [False, False],
    })
    structural = central.copy()
    structural["has_sensitivity"] = True
    structural["peak_N"] = [30.0, 35.0]
    structural["invf_min"] = [2e-8, 2.1e-8]
    structural["invf_max"] = [2e-6, 2.1e-6]

    result = structural_alternative_from_curves(central, structural)

    assert list(result["mass_GeV"]) == [1.30, 1.31]
    assert list(result["restores_sensitivity"]) == [True, False]
    assert list(result["topology_differs"]) == [True, False]
    assert (result["variation"] == "decay_2310_structural").all()

    with pytest.raises(ValueError, match="mass grids differ"):
        structural_alternative_from_curves(
            central,
            structural[structural["mass_GeV"] != 1.31],
        )


def test_dense_structural_collector_pins_curve_and_central_hashes(tmp_path):
    columns = {
        "mass_GeV": [1.30],
        "has_sensitivity": [True],
        "peak_N": [30.0],
        "peak_invf": [2e-7],
        "invf_min": [2e-8],
        "invf_max": [2e-6],
        "invf_min_open": [False],
        "invf_max_open": [False],
    }
    central_path = tmp_path / "central.csv"
    structural_path = tmp_path / "structural.csv"
    provenance_path = tmp_path / "provenance.json"
    pd.DataFrame(columns).to_csv(central_path, index=False)
    pd.DataFrame(columns).to_csv(structural_path, index=False)
    provenance = {
        "inputs": {
            "central_curve": {
                "sha256": campaign_definitions.sha256_file(central_path),
            },
        },
        "output": {
            "sha256": campaign_definitions.sha256_file(structural_path),
            "rows": 1,
        },
    }
    provenance_path.write_text(json.dumps(provenance))

    dense, observed = campaign_collector._load_dense_structural(
        central_path,
        structural_path,
        provenance_path,
    )

    assert len(dense) == 1
    assert observed == provenance

    published_curve = tmp_path / "published" / "dense_curve.csv"
    campaign_collector._atomic_copy(structural_path, published_curve)
    assert published_curve.read_bytes() == structural_path.read_bytes()
    assert campaign_definitions.sha256_file(published_curve) == (
        provenance["output"]["sha256"]
    )

    provenance["output"]["sha256"] = "wrong"
    provenance_path.write_text(json.dumps(provenance))
    with pytest.raises(ValueError, match="checksum mismatch"):
        campaign_collector._load_dense_structural(
            central_path,
            structural_path,
            provenance_path,
        )


def test_post_validation_compaction_is_atomic_and_idempotent(tmp_path):
    run_dir = tmp_path / "runs" / "pdf_0001"
    vectors = run_dir / "llp_4vectors"
    geometry = run_dir / "analysis" / "geometry_cache"
    vectors.mkdir(parents=True)
    geometry.mkdir(parents=True)
    (vectors / "mA_1.csv").write_bytes(b"vectors")
    (geometry / "geom_1.npz").write_bytes(b"geometry")
    marker = {
        "variation": {
            "name": "pdf_0001",
            "production_mode": "fresh_600k",
        },
        "sensitivity_csv": str(run_dir / "analysis" / "bc10_sensitivity.csv"),
        "production_marker": str(run_dir / "production.complete.json"),
        "storage_state": "full",
    }
    completion = run_dir / "variation.complete.json"
    completion.write_text(json.dumps(marker))
    result = _compact_completed_run(completion, marker, run_dir, False)
    assert result["storage_state"] == "compacted"
    assert result["compaction"]["reclaimed_bytes"] == len(b"vectorsgeometry")
    assert not vectors.exists()
    assert not geometry.exists()
    assert json.loads(completion.read_text())["storage_state"] == "compacted"
    assert _compact_completed_run(completion, result, run_dir, False) == result
