import pandas as pd

from scalar import production
from scalar.run_sensitivity import _reconstruction_seed, _write_checkpoint


def test_reconstruction_seed_is_stable_for_canonical_subsets():
    for index in (0, len(production.MASS_GRID) // 2, len(production.MASS_GRID) - 1):
        mass = production.MASS_GRID[index]
        assert _reconstruction_seed(mass) == 1000 + index


def test_checkpoint_is_sorted_and_replaced_atomically(tmp_path):
    output = tmp_path / "analysis" / "bc4_sensitivity.csv"
    _write_checkpoint(
        [{"mass_GeV": 1.0, "peak_N": 2.0}, {"mass_GeV": 0.5, "peak_N": 1.0}],
        output,
    )

    frame = pd.read_csv(output)
    assert frame["mass_GeV"].tolist() == [0.5, 1.0]
    assert not output.with_suffix(output.suffix + ".tmp").exists()
