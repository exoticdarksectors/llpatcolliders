"""Named, pinned fermionic-ALP decay models used by BC10."""

from __future__ import annotations

from pathlib import Path


DEFAULT_DECAY_MODEL = "2501"
STRUCTURAL_DECAY_MODEL = "2310_structural"
DECAY_MODEL_SPECS = {
    DEFAULT_DECAY_MODEL: {
        "data_subdir": "senscalc_2501",
        "phenomenology": "arXiv:2501.04525",
        "role": "central",
    },
    STRUCTURAL_DECAY_MODEL: {
        "data_subdir": "senscalc_2310",
        "phenomenology": "arXiv:2310.03524",
        "role": "one_sided_heavy_pseudoscalar_structural_alternative",
    },
}


def validate_decay_model(decay_model: str) -> str:
    if decay_model not in DECAY_MODEL_SPECS:
        choices = ", ".join(DECAY_MODEL_SPECS)
        raise ValueError(f"unknown decay model {decay_model!r}; choose {choices}")
    return decay_model


def decay_data_dir(decay_model: str = DEFAULT_DECAY_MODEL) -> Path:
    model = validate_decay_model(decay_model)
    return (
        Path(__file__).resolve().parent
        / "data"
        / DECAY_MODEL_SPECS[model]["data_subdir"]
    )
