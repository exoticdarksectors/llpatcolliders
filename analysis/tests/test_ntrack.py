"""
Tests for the N-track analysis module.

Verifies integrate_decay_prob, merge_intervals, and
union_pair_window_probability.

Pure numpy — no ROOT dependency.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "geometry"))

from decayProbPerEvent_Ntrack import (
    integrate_decay_prob,
    merge_intervals,
    union_pair_window_probability,
)
from detector_cuts import SEP_MAX, SEP_MIN


# ── integrate_decay_prob tests ───────────────────────────────────────────


class TestIntegrateDecayProb(unittest.TestCase):

    def test_known_value(self):
        """exp(-10/50) - exp(-12/50) = known value."""
        d_start, d_end, dl = 10.0, 12.0, 50.0
        expected = np.exp(-d_start / dl) - np.exp(-d_end / dl)
        actual = integrate_decay_prob(d_start, d_end, dl)
        self.assertAlmostEqual(actual, expected, places=12)

    def test_zero_decay_length(self):
        """Zero decay length should return 0."""
        self.assertEqual(integrate_decay_prob(10.0, 12.0, 0.0), 0.0)

    def test_negative_decay_length(self):
        """Negative decay length should return 0."""
        self.assertEqual(integrate_decay_prob(10.0, 12.0, -1.0), 0.0)

    def test_reversed_interval(self):
        """d_end <= d_start should return 0."""
        self.assertEqual(integrate_decay_prob(12.0, 10.0, 50.0), 0.0)
        self.assertEqual(integrate_decay_prob(10.0, 10.0, 50.0), 0.0)

    def test_probability_bounded(self):
        """Result should be in [0, 1]."""
        for d_start, d_end, dl in [(0, 1, 10), (10, 12, 50), (100, 200, 1000)]:
            p = integrate_decay_prob(float(d_start), float(d_end), float(dl))
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)

    def test_expm1_accuracy_small_interval(self):
        """For tiny intervals << decay_length, expm1 should avoid cancellation."""
        # d_end - d_start = 1e-10, decay_length = 1e6
        d_start = 10.0
        d_end = 10.0 + 1e-10
        dl = 1e6
        p = integrate_decay_prob(d_start, d_end, dl)
        # Naive: exp(-10/1e6) - exp(-10.0000000001/1e6) would lose precision
        # With expm1: should be ~ exp(-10/1e6) * (1e-10 / 1e6) = ~1e-16
        expected = np.exp(-d_start / dl) * (1e-10 / dl)
        self.assertAlmostEqual(p, expected, delta=expected * 0.01)

    def test_full_survival_near_zero(self):
        """Starting from d=0 with short path, probability is well-defined."""
        p = integrate_decay_prob(0.0, 0.1, 100.0)
        expected = 1.0 - np.exp(-0.1 / 100.0)
        self.assertAlmostEqual(p, expected, places=10)


# ── merge_intervals tests ───────────────────────────────────────────────


class TestMergeIntervals(unittest.TestCase):

    def test_no_intervals(self):
        """Empty input returns empty."""
        self.assertEqual(merge_intervals([]), [])

    def test_single_interval(self):
        """Single valid interval is returned as-is."""
        result = merge_intervals([(1.0, 3.0)])
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0][0], 1.0)
        self.assertAlmostEqual(result[0][1], 3.0)

    def test_disjoint(self):
        """Disjoint intervals stay separate."""
        result = merge_intervals([(1.0, 2.0), (5.0, 6.0)])
        self.assertEqual(len(result), 2)

    def test_overlapping(self):
        """Overlapping intervals are merged."""
        result = merge_intervals([(1.0, 3.0), (2.0, 5.0)])
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0][0], 1.0)
        self.assertAlmostEqual(result[0][1], 5.0)

    def test_touching(self):
        """Touching intervals are merged."""
        result = merge_intervals([(1.0, 3.0), (3.0, 5.0)])
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0][0], 1.0)
        self.assertAlmostEqual(result[0][1], 5.0)

    def test_subset(self):
        """One interval contained in another."""
        result = merge_intervals([(1.0, 10.0), (3.0, 5.0)])
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0][0], 1.0)
        self.assertAlmostEqual(result[0][1], 10.0)

    def test_invalid_intervals_dropped(self):
        """Intervals with lo >= hi are dropped."""
        result = merge_intervals([(5.0, 3.0), (1.0, 2.0)])
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0][0], 1.0)
        self.assertAlmostEqual(result[0][1], 2.0)

    def test_unsorted_input(self):
        """Unsorted input should still be handled correctly."""
        result = merge_intervals([(5.0, 6.0), (1.0, 3.0), (2.0, 4.0)])
        self.assertEqual(len(result), 2)
        self.assertAlmostEqual(result[0][0], 1.0)
        self.assertAlmostEqual(result[0][1], 4.0)
        self.assertAlmostEqual(result[1][0], 5.0)
        self.assertAlmostEqual(result[1][1], 6.0)

    def test_all_invalid(self):
        """All invalid intervals returns empty."""
        result = merge_intervals([(5, 3), (4, 2), (10, 10)])
        self.assertEqual(len(result), 0)


# ── union_pair_window_probability tests ──────────────────────────────────


class TestUnionPairWindowProbability(unittest.TestCase):

    def test_empty_thetas(self):
        """No thetas should give zero probability."""
        p, best = union_pair_window_probability(20.0, 22.0, [], 100.0)
        self.assertEqual(p, 0.0)
        self.assertTrue(np.isnan(best))

    def test_zero_decay_length(self):
        """Zero decay length should give zero probability."""
        p, best = union_pair_window_probability(
            20.0, 22.0, np.array([0.01]), 0.0)
        self.assertEqual(p, 0.0)

    def test_single_theta_positive(self):
        """Single reasonable theta should give positive probability."""
        # theta = 0.01 rad, entry=20, exit=22, sep at exit ~ theta*(exit-d)
        # For sep to be in [0.001, 1.0], d must be in range
        # sep = 0.01 * (22 - d) → sep_max/theta = 100, sep_min/theta = 0.1
        # d_window: max(20, 22-100)=20 to min(22, 22-0.1)=21.9
        thetas = np.array([0.01])
        p, best = union_pair_window_probability(20.0, 22.0, thetas, 100.0)
        self.assertGreater(p, 0.0)
        self.assertAlmostEqual(best, 0.01, places=10)

    def test_multiple_thetas_union(self):
        """Two thetas should give >= max(individual probs)."""
        thetas1 = np.array([0.01])
        thetas2 = np.array([0.01, 0.005])
        p1, _ = union_pair_window_probability(20.0, 22.0, thetas1, 100.0)
        p2, _ = union_pair_window_probability(20.0, 22.0, thetas2, 100.0)
        self.assertGreaterEqual(p2, p1 - 1e-15)

    def test_zero_theta_ignored(self):
        """Theta = 0 should be skipped."""
        thetas = np.array([0.0, 0.01])
        p, best = union_pair_window_probability(20.0, 22.0, thetas, 100.0)
        self.assertGreater(p, 0.0)

    def test_no_valid_window(self):
        """If all windows are empty, probability is zero."""
        # theta very large → window starts after exit
        thetas = np.array([100.0])  # sep_min/theta ~ 0, exit - sep_max/theta ~ 22-0.01
        # Actually this gives a valid window. Let's use tiny entry/exit range
        # entry=20, exit=20.001, theta=0.0001
        # d_window: max(20, 20.001 - 10000) = 20, min(20.001, 20.001 - 10) = invalid
        thetas = np.array([0.0001])
        p, _ = union_pair_window_probability(20.0, 20.001, thetas, 100.0,
                                              sep_min=SEP_MIN, sep_max=SEP_MAX)
        # sep_min/theta = 10, exit - 10 = 10.001 < entry=20 → no window
        self.assertEqual(p, 0.0)

    def test_probability_bounded(self):
        """Probability should be in [0, 1]."""
        thetas = np.array([0.005, 0.01, 0.02])
        p, _ = union_pair_window_probability(20.0, 22.0, thetas, 50.0)
        self.assertGreaterEqual(p, 0.0)
        self.assertLessEqual(p, 1.0)


if __name__ == "__main__":
    unittest.main()
