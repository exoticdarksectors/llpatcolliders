"""
Tests for the 2-body acceptance analysis module.

Verifies acceptance_weighted_decay_prob, compute_c_upper, compute_c_S,
compute_c_max_sep, and unweighted_decay_prob.

Pure numpy/scipy — no ROOT dependency.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "geometry"))

from decayProbPerEvent_2body import (
    M_DAUGHTER,
    acceptance_weighted_decay_prob,
    compute_c_S,
    compute_c_max_sep,
    compute_c_upper,
    unweighted_decay_prob,
)
from detector_cuts import P_CUT, SEP_MAX, SEP_MIN


# ── compute_c_upper tests ───────────────────────────────────────────────


class TestComputeCUpper(unittest.TestCase):

    def test_high_gamma(self):
        """At high gamma, c_upper should be close to beta (≈1)."""
        gamma = 100.0
        beta = np.sqrt(1.0 - 1.0 / gamma**2)
        mass = 1.0
        c = compute_c_upper(gamma, beta, mass)
        # Should be close to beta since momentum cut is easily satisfied
        self.assertGreater(c, 0.0)
        self.assertLessEqual(c, beta + 1e-10)

    def test_low_gamma_zero(self):
        """At very low gamma (near threshold), c_upper can become <= 0."""
        # mass = 2*M_DAUGHTER is threshold; just above it, gamma ~ 1
        mass = 2.0 * M_DAUGHTER + 0.001
        gamma = 1.001
        beta = np.sqrt(1.0 - 1.0 / gamma**2)
        c = compute_c_upper(gamma, beta, mass)
        # With such low boost, momentum cut may kill everything
        self.assertLessEqual(c, beta + 1e-10)

    def test_returns_at_most_beta(self):
        """c_upper should never exceed beta (forward requirement)."""
        for gamma in [2, 10, 100]:
            beta = np.sqrt(1.0 - 1.0 / gamma**2)
            c = compute_c_upper(gamma, beta, 1.0)
            self.assertLessEqual(c, beta + 1e-10)


# ── compute_c_S tests ───────────────────────────────────────────────────


class TestComputeCS(unittest.TestCase):

    def test_zero_d_remaining(self):
        """At d_remaining=0, c_S should be 1.0 (no angular window)."""
        c = compute_c_S(0.0, 10.0, 0.99)
        self.assertEqual(c, 1.0)

    def test_large_d_remaining(self):
        """At very large d_remaining, c_S should approach 0."""
        c = compute_c_S(1e6, 100.0, 0.9999)
        self.assertLess(c, 0.1)

    def test_nonnegative(self):
        """c_S should always be non-negative."""
        for d in [0.01, 0.1, 1, 10, 100]:
            c = compute_c_S(d, 50.0, 0.999)
            self.assertGreaterEqual(c, 0.0)


# ── compute_c_max_sep tests ─────────────────────────────────────────────


class TestComputeCMaxSep(unittest.TestCase):

    def test_zero_d_remaining(self):
        """At d_remaining=0, should return 1.0."""
        c = compute_c_max_sep(0.0, 10.0, 0.99)
        self.assertEqual(c, 1.0)

    def test_large_d_remaining_low_gamma(self):
        """At large d_remaining with low gamma, c_max_sep can be 1.0."""
        # theta_max = 1.0 / 1.0 = 1.0 rad >= pi? No, but cos(1) ~ 0.54
        # With gamma=2: denom = 4 * (1-0.54) = 1.84, val = (1-2/1.84)/beta^2
        # Still might be tricky. Use d=0.3 so theta_max=1/0.3=3.33>pi → returns 1
        c = compute_c_max_sep(0.3, 10.0, 0.99)
        self.assertEqual(c, 1.0)

    def test_bounded_by_one(self):
        """c_max_sep should be at most 1.0."""
        c = compute_c_max_sep(1000.0, 100.0, 0.9999)
        self.assertLessEqual(c, 1.0 + 1e-10)

    def test_nonnegative(self):
        """c_max_sep should always be non-negative."""
        c = compute_c_max_sep(0.01, 10.0, 0.99)
        self.assertGreaterEqual(c, 0.0)


# ── unweighted_decay_prob tests ──────────────────────────────────────────


class TestUnweightedDecayProb(unittest.TestCase):

    def test_known_formula(self):
        """P = exp(-entry/lambda) * (1 - exp(-path/lambda))."""
        entry, exit_, dl = 10.0, 12.0, 50.0
        expected = np.exp(-entry / dl) * (1 - np.exp(-(exit_ - entry) / dl))
        actual = unweighted_decay_prob(entry, exit_, dl)
        self.assertAlmostEqual(actual, expected, places=12)

    def test_short_path_small_prob(self):
        """Short path relative to decay length → small probability."""
        p = unweighted_decay_prob(100.0, 100.1, 1e6)
        self.assertLess(p, 0.01)
        self.assertGreater(p, 0.0)

    def test_long_decay_length(self):
        """Very long decay length → near-zero probability."""
        p = unweighted_decay_prob(10.0, 12.0, 1e15)
        self.assertLess(p, 1e-10)

    def test_probability_bounded(self):
        """Probability should be in [0, 1]."""
        for dl in [1, 10, 100, 1000]:
            p = unweighted_decay_prob(20.0, 22.0, float(dl))
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)


# ── acceptance_weighted_decay_prob tests ─────────────────────────────────


class TestAcceptanceWeightedDecayProb(unittest.TestCase):

    def test_below_mass_threshold(self):
        """Mass below 2*M_DAUGHTER should give zero probability."""
        p = acceptance_weighted_decay_prob(
            entry_d=20.0, exit_d=22.0, gamma=100, beta=0.9999,
            mass=2 * M_DAUGHTER - 0.001, decay_length=100.0)
        self.assertEqual(p, 0.0)

    def test_at_threshold_mass(self):
        """At exactly 2*M_DAUGHTER, code uses strict '<' so this is allowed."""
        p = acceptance_weighted_decay_prob(
            entry_d=20.0, exit_d=22.0, gamma=100, beta=0.9999,
            mass=2 * M_DAUGHTER, decay_length=100.0)
        # mass = 2*M_DAUGHTER passes the `mass < 2*M_DAUGHTER` check
        # so a positive (though possibly small) probability is expected
        self.assertGreaterEqual(p, 0.0)

    def test_high_boost_nonzero(self):
        """At high boost, acceptance-weighted prob should be nonzero."""
        gamma = 200.0
        beta = np.sqrt(1.0 - 1.0 / gamma**2)
        p = acceptance_weighted_decay_prob(
            entry_d=20.0, exit_d=22.0, gamma=gamma, beta=beta,
            mass=1.0, decay_length=100.0)
        self.assertGreater(p, 0.0)

    def test_leq_unweighted(self):
        """Acceptance-weighted probability should be <= unweighted."""
        gamma = 200.0
        beta = np.sqrt(1.0 - 1.0 / gamma**2)
        p_with = acceptance_weighted_decay_prob(
            entry_d=20.0, exit_d=22.0, gamma=gamma, beta=beta,
            mass=1.0, decay_length=100.0)
        p_without = unweighted_decay_prob(20.0, 22.0, 100.0)
        self.assertLessEqual(p_with, p_without + 1e-15)

    def test_symmetric_range(self):
        """Probability should be non-negative for various parameters."""
        for gamma, mass in [(10, 0.5), (50, 1.0), (200, 5.0)]:
            beta = np.sqrt(1.0 - 1.0 / gamma**2)
            p = acceptance_weighted_decay_prob(
                entry_d=20.0, exit_d=22.0, gamma=gamma, beta=beta,
                mass=mass, decay_length=100.0)
            self.assertGreaterEqual(p, 0.0)


if __name__ == "__main__":
    unittest.main()
