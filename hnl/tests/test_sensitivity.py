"""Tests for analysis.sensitivity — physics engine regression tests."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import unittest
import numpy as np

from analysis.sensitivity import (
    compute_c_upper, compute_c_S, compute_c_max_sep,
    acceptance_weighted_decay_prob, unweighted_decay_prob,
    compute_n_signal,
)


class TestAcceptanceBounds(unittest.TestCase):

    def test_c_upper_physical(self):
        """c_upper should be in [0, beta] for physical parameters."""
        gamma, beta, mass = 100.0, 0.99995, 1.0
        c = compute_c_upper(gamma, beta, mass)
        self.assertGreaterEqual(c, 0.0)
        self.assertLessEqual(c, beta)

    def test_c_upper_low_boost(self):
        """At very low boost, c_upper may go negative (no acceptance)."""
        gamma, beta, mass = 1.01, 0.14, 1.0
        c = compute_c_upper(gamma, beta, mass, p_cut=0.6)
        # At low boost, p_cut may not be satisfiable
        # c_upper can be negative → no accepted decays
        self.assertLessEqual(c, beta)

    def test_c_S_large_distance(self):
        """At large distance, minimum separation is easy → c_S near 0."""
        c = compute_c_S(d_remaining=100.0, gamma=100.0, beta=0.9999)
        self.assertAlmostEqual(c, 0.0, places=3)

    def test_c_S_small_distance(self):
        """At very small remaining distance, c_S should be large."""
        c = compute_c_S(d_remaining=0.001, gamma=10.0, beta=0.995)
        self.assertGreater(c, 0.5)

    def test_c_max_sep_small_distance(self):
        """At small remaining distance, max sep reached easily → c near 1."""
        c = compute_c_max_sep(d_remaining=0.5, gamma=10.0, beta=0.995)
        # With d_remaining=0.5m and sep_max=1.0m, theta_max >= pi → c = 1.0
        self.assertGreaterEqual(c, 0.99)


class TestDecayProbability(unittest.TestCase):

    def test_unweighted_basic(self):
        """Unweighted P_decay should be positive and < 1."""
        p = unweighted_decay_prob(entry_d=25.0, exit_d=27.0, decay_length=100.0)
        self.assertGreater(p, 0.0)
        self.assertLess(p, 1.0)

    def test_unweighted_zero_path(self):
        """Zero path length → zero probability."""
        p = unweighted_decay_prob(entry_d=25.0, exit_d=25.0, decay_length=100.0)
        self.assertEqual(p, 0.0)

    def test_acceptance_leq_unweighted(self):
        """Acceptance-weighted probability should be <= unweighted."""
        entry, exit_, gamma, beta, mass = 25.0, 27.0, 50.0, 0.9998, 1.0
        lam = 100.0

        p_acc = acceptance_weighted_decay_prob(
            entry, exit_, gamma, beta, mass, lam)
        p_raw = unweighted_decay_prob(entry, exit_, lam)

        self.assertLessEqual(p_acc, p_raw + 1e-15)

    def test_acceptance_high_boost(self):
        """At very high boost, acceptance should be non-zero for moderate path."""
        entry, exit_ = 22.0, 24.5
        gamma, beta, mass = 500.0, 0.999998, 1.0
        lam = 1000.0

        p_acc = acceptance_weighted_decay_prob(
            entry, exit_, gamma, beta, mass, lam)
        self.assertGreater(p_acc, 0.0)


class TestNSignal(unittest.TestCase):

    def test_zero_u2(self):
        """N_signal should be 0 at U²=0."""
        N = compute_n_signal(
            weight=np.array([1e-3]),
            beta_gamma=np.array([100.0]),
            gamma=np.array([100.005]),
            beta=np.array([0.99995]),
            hits=np.array([True]),
            entry_d=np.array([25.0]),
            exit_d=np.array([27.0]),
            ctau_u2_1=1.0,
            u2=0.0,
            m_N=1.0,
            L_int_pb=3e6,
        )
        self.assertEqual(N, 0.0)

    def test_n_signal_positive(self):
        """N_signal should be positive for reasonable parameters."""
        N = compute_n_signal(
            weight=np.array([1e-3]),
            beta_gamma=np.array([100.0]),
            gamma=np.array([100.005]),
            beta=np.array([0.99995]),
            hits=np.array([True]),
            entry_d=np.array([25.0]),
            exit_d=np.array([27.0]),
            ctau_u2_1=1.0,
            u2=1e-6,
            m_N=1.0,
            L_int_pb=3e6,
        )
        self.assertGreater(N, 0.0)

    def test_no_hits_gives_zero(self):
        """If no events hit the detector, N_signal = 0."""
        N = compute_n_signal(
            weight=np.array([1e-3, 2e-3]),
            beta_gamma=np.array([100.0, 200.0]),
            gamma=np.array([100.0, 200.0]),
            beta=np.array([0.9999, 0.9999]),
            hits=np.array([False, False]),
            entry_d=np.array([np.nan, np.nan]),
            exit_d=np.array([np.nan, np.nan]),
            ctau_u2_1=1.0,
            u2=1e-6,
            m_N=1.0,
            L_int_pb=3e6,
        )
        self.assertEqual(N, 0.0)

    def test_n_signal_scales_with_luminosity(self):
        """Doubling luminosity should double N_signal."""
        kwargs = dict(
            weight=np.array([1e-3]),
            beta_gamma=np.array([100.0]),
            gamma=np.array([100.005]),
            beta=np.array([0.99995]),
            hits=np.array([True]),
            entry_d=np.array([25.0]),
            exit_d=np.array([27.0]),
            ctau_u2_1=1.0,
            u2=1e-6,
            m_N=1.0,
            use_acceptance=False,
        )
        N1 = compute_n_signal(L_int_pb=3e6, **kwargs)
        N2 = compute_n_signal(L_int_pb=6e6, **kwargs)
        np.testing.assert_allclose(N2, 2.0 * N1, rtol=1e-10)


class TestVectorizedAcceptance(unittest.TestCase):

    def test_c_max_sep_vec_matches_scalar_2d(self):
        """Vectorized c_max_sep must match scalar for 2D input (quadrature grid)."""
        from analysis.sensitivity import compute_c_max_sep, compute_c_max_sep_vec
        rng = np.random.default_rng(42)
        n_hits, n_quad = 5, 30
        d_remaining = rng.uniform(0.5, 50.0, (n_hits, n_quad))
        gamma = rng.uniform(10, 500, (n_hits, n_quad))
        beta = 1.0 - 0.5 / gamma**2

        vec = compute_c_max_sep_vec(d_remaining, gamma, beta)
        scalar = np.empty_like(d_remaining)
        for i in range(n_hits):
            for j in range(n_quad):
                scalar[i, j] = compute_c_max_sep(
                    d_remaining[i, j], gamma[i, j], beta[i, j])

        np.testing.assert_allclose(vec, scalar, atol=1e-12)

    def test_c_max_sep_vec_matches_scalar_1d(self):
        """Vectorized c_max_sep must match scalar for 1D input."""
        from analysis.sensitivity import compute_c_max_sep, compute_c_max_sep_vec
        rng = np.random.default_rng(7)
        d_remaining = rng.uniform(0.5, 50.0, 20)
        gamma = rng.uniform(10, 500, 20)
        beta = 1.0 - 0.5 / gamma**2

        vec = compute_c_max_sep_vec(d_remaining, gamma, beta)
        scalar = np.array([
            compute_c_max_sep(d_remaining[i], gamma[i], beta[i])
            for i in range(20)])

        np.testing.assert_allclose(vec, scalar, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
