import hashlib
import json

import numpy as np
import pytest

from alp_fermion import decay_matrix_elements as matrix_elements
from alp_fermion import exclusive_decays
from alp_fermion import model
from alp_fermion.decay_models import decay_data_dir, validate_decay_model
from alp_fermion.generate_decay_templates_pythia import (
    _masses_from_csv,
    _validate_resumed_template,
)


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_2310_export_manifest_pins_sources_and_every_output_hash():
    data_dir = decay_data_dir("2310_structural")
    manifest = json.loads((data_dir / "EXPORT_MANIFEST.json").read_text())
    assert manifest["senscalc_commit"] == (
        "0bca050633aae16e148d47f21840fa07ff4b8724"
    )
    assert manifest["decay_model"] == "2310_structural"
    assert manifest["phenomenology"] == "arXiv:2310.03524"
    assert manifest["source_sha256"]["widths"] == (
        "40b51e8d35a8edbb3027cfcc3a2bb369faff8fc009a958863d6264f227410e93"
    )
    for name, expected in manifest["output_sha256"].items():
        assert _sha256(data_dir / name) == expected


def test_2310_channel_taxonomy_matches_2501_positionally():
    keys = {
        "id", "canonical_name", "process_name", "decay_products",
        "process_input_form", "decay_products_input_form",
    }
    central = json.loads(
        (decay_data_dir("2501") / "decay_channels.json").read_text()
    )
    structural = json.loads(
        (decay_data_dir("2310_structural") / "decay_channels.json").read_text()
    )
    assert len(central) == len(structural) == 32
    assert [
        {key: row[key] for key in keys} for row in structural
    ] == [
        {key: row[key] for key in keys} for row in central
    ]


def test_2310_width_and_branching_anchors_and_perturbative_join():
    inv_f = model.INV_F_REF
    width = model.alp_total_width(
        1.30, inv_f, decay_model="2310_structural"
    )
    assert width / inv_f**2 == pytest.approx(
        0.012317291370015708, rel=1e-13
    )
    branching = model.alp_branchings(
        1.30, inv_f, decay_model="2310_structural"
    )
    assert branching["mumu"] == pytest.approx(
        0.05096158623008926, rel=1e-13
    )
    central = model.alp_total_width(1.30, inv_f)
    assert central / width == pytest.approx(88.789, rel=2e-5)
    assert model.alp_total_width(
        2.17, inv_f, decay_model="2310_structural"
    ) == pytest.approx(model.alp_total_width(2.17, inv_f), rel=1e-14)


def test_2310_exclusive_weights_are_independent_and_normalized():
    central = exclusive_decays.exclusive_branching_weights(1.30, "2501")
    structural = exclusive_decays.exclusive_branching_weights(
        1.30, "2310_structural"
    )
    assert sum(structural.values()) == pytest.approx(1.0, abs=1e-12)
    assert structural != central
    assert not exclusive_decays.DUPLICATE_CHANNEL_IDS.intersection(structural)


def test_2310_exact_matrix_element_matches_wolfram_anchor():
    value = matrix_elements.matrix_element_squared(
        "matrix_element_004",
        1.2,
        np.array([0.4]),
        np.array([0.3]),
        "2310_structural",
    )
    assert value[0] == pytest.approx(2.774159983255917, rel=1e-13)


def test_unknown_decay_model_is_rejected():
    with pytest.raises(ValueError, match="unknown decay model"):
        validate_decay_model("approximate")


def test_template_resume_is_model_count_and_surrogate_specific(tmp_path):
    path = tmp_path / "template.npz"
    np.savez_compressed(
        path,
        decay_model=np.array("2310_structural"),
        n_templates=np.int32(20_000),
        gluon_surrogate=np.array("uds"),
    )
    _validate_resumed_template(path, 20_000, "uds", "2310_structural")
    with pytest.raises(RuntimeError, match="model/count/surrogate"):
        _validate_resumed_template(path, 20_000, "uds", "2501")


def test_template_generator_accepts_exact_mass_grid_csv(tmp_path):
    path = tmp_path / "mass_grid.csv"
    path.write_text("mass_GeV\n1.18\n1.19\n1.20\n")
    assert _masses_from_csv(path) == [1.18, 1.19, 1.20]

    path.write_text("mass_GeV\n1.18\n1.20\n1.19\n")
    with pytest.raises(ValueError, match="strictly increasing"):
        _masses_from_csv(path)
