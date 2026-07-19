"""Tests for the SensCalc LHC parent-spectrum boundary."""

import math

import numpy as np
import pytest

from alp_fermion.production_spectra import (
    LIGHT_MESON_PARENTS,
    TWO_BODY_PRODUCTION_CHANNELS,
    TabulatedSpectrum,
    pseudorapidity,
    sample_two_body_alps,
    transverse_theta_interval,
)


def _write_grid(path, theta, energy, density):
    rows = []
    for i, theta_value in enumerate(theta):
        for j, energy_value in enumerate(energy):
            rows.append((theta_value, energy_value, density[i, j]))
    np.savetxt(path, rows)


def test_load_integrate_and_select_polar_interval(tmp_path):
    path = tmp_path / "grid.txt"
    theta = np.array([0.2, 1.0, 2.0, 2.8])
    energy = np.array([1.0, 2.0, 4.0])
    density = np.ones((len(theta), len(energy)))
    _write_grid(path, theta, energy, density)

    spectrum = TabulatedSpectrum.load(path)

    assert spectrum.integral() == pytest.approx((2.8 - 0.2) * (4.0 - 1.0))
    assert spectrum.polar_fraction(1.0, 2.0) == pytest.approx(1.0 / 2.6)


def test_sample_returns_on_shell_four_vectors(tmp_path):
    path = tmp_path / "grid.txt"
    theta = np.array([0.4, 1.2, 2.0, 2.7])
    energy = np.array([0.4, 2.0, 4.0])
    density = np.ones((len(theta), len(energy)))
    _write_grid(path, theta, energy, density)
    spectrum = TabulatedSpectrum.load(path)

    mass = 0.5
    sample = spectrum.sample(20_000, mass, rng=np.random.default_rng(1234))
    momentum_squared = sample["px"] ** 2 + sample["py"] ** 2 + sample["pz"] ** 2

    assert len(sample["E"]) == 20_000
    assert np.all(sample["E"] >= mass)
    assert np.allclose(sample["E"] ** 2 - momentum_squared, mass**2)
    assert np.all((sample["theta"] >= theta.min()) & (sample["theta"] <= theta.max()))


def test_load_rejects_incomplete_cartesian_grid(tmp_path):
    path = tmp_path / "bad.txt"
    np.savetxt(path, [(0.1, 1.0, 1.0), (0.2, 2.0, 1.0)])
    with pytest.raises(ValueError, match="not a Cartesian grid"):
        TabulatedSpectrum.load(path)


def test_load_fragmentation_mass_slice(tmp_path):
    path = tmp_path / "fragmentation.csv"
    path.write_text(
        "mass_GeV,theta_rad,energy_GeV,density\n"
        "0.5,0.2,1.0,1.0\n"
        "0.5,0.2,2.0,2.0\n"
        "0.5,1.0,1.0,3.0\n"
        "0.5,1.0,2.0,4.0\n"
        "0.6,0.2,1.0,5.0\n"
        "0.6,0.2,2.0,6.0\n"
        "0.6,1.0,1.0,7.0\n"
        "0.6,1.0,2.0,8.0\n"
    )

    spectrum = TabulatedSpectrum.load_mass_slice(path, 0.6)

    assert np.array_equal(spectrum.theta, [0.2, 1.0])
    assert np.array_equal(spectrum.energy, [1.0, 2.0])
    assert np.array_equal(spectrum.density, [[5.0, 6.0], [7.0, 8.0]])
    with pytest.raises(ValueError, match="is not tabulated"):
        TabulatedSpectrum.load_mass_slice(path, 0.7)


def test_transverse_theta_interval_matches_pseudorapidity_definition():
    theta_min, theta_max = transverse_theta_interval(0.5)
    assert theta_min == pytest.approx(2.0 * math.atan(math.exp(-0.5)))
    assert theta_max == pytest.approx(math.pi - theta_min)


def test_two_body_alp_sample_is_on_shell_at_rest():
    parent = LIGHT_MESON_PARENTS["KS"]
    channel = TWO_BODY_PRODUCTION_CHANNELS["KS-to-ALP-Pi0"]
    n_events = 20_000
    parent_vectors = {
        "E": np.full(n_events, parent.mass_gev),
        "px": np.zeros(n_events),
        "py": np.zeros(n_events),
        "pz": np.zeros(n_events),
    }
    alp_mass = 0.22

    alp = sample_two_body_alps(
        parent_vectors,
        parent,
        channel,
        alp_mass,
        np.random.default_rng(1234),
    )

    momentum_squared = alp["px"] ** 2 + alp["py"] ** 2 + alp["pz"] ** 2
    expected_energy = (
        parent.mass_gev**2 + alp_mass**2 - channel.recoil_mass_gev**2
    ) / (2.0 * parent.mass_gev)
    assert np.allclose(alp["E"], expected_energy)
    assert np.allclose(alp["E"] ** 2 - momentum_squared, alp_mass**2)


def test_two_body_alp_sample_rejects_closed_channel():
    parent = LIGHT_MESON_PARENTS["KS"]
    channel = TWO_BODY_PRODUCTION_CHANNELS["KS-to-ALP-Pi0"]
    parent_vectors = {
        "E": np.array([parent.mass_gev]),
        "px": np.zeros(1),
        "py": np.zeros(1),
        "pz": np.zeros(1),
    }
    with pytest.raises(ValueError, match="is closed"):
        sample_two_body_alps(
            parent_vectors,
            parent,
            channel,
            parent.mass_gev - channel.recoil_mass_gev,
            np.random.default_rng(1234),
        )


def test_pseudorapidity_from_cartesian_momentum():
    expected = np.array([-0.5, 0.0, 0.5])
    vectors = {
        "px": np.ones(3),
        "py": np.zeros(3),
        "pz": np.sinh(expected),
    }
    assert np.allclose(pseudorapidity(vectors), expected)
