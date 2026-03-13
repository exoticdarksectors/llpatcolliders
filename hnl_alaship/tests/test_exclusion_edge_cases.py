"""
Tests for exclusion band and island polygon edge cases.

Pure numpy — no ROOT dependency.
"""

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from analysis.exclusion import find_exclusion_band
from analysis.constants import N_THRESHOLD


# ── find_exclusion_band edge cases ───────────────────────────────────────


class TestFindExclusionBand(unittest.TestCase):

    def test_no_sensitivity(self):
        """All N_signal < threshold should return has_sensitivity=False."""
        u2 = np.logspace(-10, -3, 50)
        N = np.ones(50) * 0.5  # all below N_THRESHOLD=3
        result = find_exclusion_band(u2, N)
        self.assertFalse(result["has_sensitivity"])
        self.assertTrue(np.isnan(result["u2_min"]))
        self.assertTrue(np.isnan(result["u2_max"]))
        self.assertAlmostEqual(result["peak_N"], 0.5)

    def test_all_above_threshold(self):
        """All N_signal > threshold should give band spanning full grid."""
        u2 = np.logspace(-10, -3, 50)
        N = np.ones(50) * 100.0
        result = find_exclusion_band(u2, N)
        self.assertTrue(result["has_sensitivity"])
        # u2_min should be at edge (falls back to grid point)
        self.assertAlmostEqual(result["u2_min"], u2[0])
        self.assertAlmostEqual(result["u2_max"], u2[-1])

    def test_single_point_above(self):
        """Only one point above threshold should still give a band."""
        u2 = np.logspace(-10, -3, 50)
        N = np.ones(50) * 0.1
        N[25] = 10.0  # single peak above threshold
        result = find_exclusion_band(u2, N)
        self.assertTrue(result["has_sensitivity"])
        self.assertLessEqual(result["u2_min"], result["u2_max"])

    def test_peaked_island_shape(self):
        """Standard island: N rises, peaks, falls → u2_min < u2_max."""
        u2 = np.logspace(-10, -3, 200)
        # Physical island: N ∝ U² * exp(-ctau/(U² * beta_gamma * d))
        ctau = 1e-10
        N = 1e6 * u2 * np.exp(-ctau / (u2 * 50.0 * 22.0))
        result = find_exclusion_band(u2, N)
        if result["has_sensitivity"]:
            self.assertLess(result["u2_min"], result["u2_max"])
            self.assertGreater(result["peak_N"], N_THRESHOLD)

    def test_flat_at_threshold(self):
        """N_signal exactly at threshold everywhere should count as sensitive."""
        u2 = np.logspace(-10, -3, 50)
        N = np.ones(50) * N_THRESHOLD
        result = find_exclusion_band(u2, N)
        self.assertTrue(result["has_sensitivity"])

    def test_monotonically_decreasing(self):
        """Monotonically decreasing N should give lower edge only."""
        u2 = np.logspace(-10, -3, 50)
        N = np.logspace(2, -2, 50)  # decreasing from 100 to 0.01
        result = find_exclusion_band(u2, N)
        self.assertTrue(result["has_sensitivity"])
        # u2_min should be near start, u2_max somewhere in middle
        self.assertLess(result["u2_min"], result["u2_max"])

    def test_interpolation_at_edges(self):
        """Interpolation should give sub-grid-spacing precision."""
        u2 = np.logspace(-10, -3, 100)
        # Sharp step: N=10 for indices 30-70, N=0 otherwise
        N = np.zeros(100)
        N[30:71] = 10.0
        result = find_exclusion_band(u2, N)
        self.assertTrue(result["has_sensitivity"])
        # Interpolated u2_min should be between u2[29] and u2[30]
        self.assertGreater(result["u2_min"], u2[29])
        self.assertLessEqual(result["u2_min"], u2[30])


# ── _build_island_polygon edge cases ────────────────────────────────────


class TestBuildIslandPolygon(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from analysis.plot_exclusion import _build_island_polygon
        cls._build_island_polygon = staticmethod(_build_island_polygon)

    def test_no_sensitive_points(self):
        """Empty sensitivity → returns None."""
        df = pd.DataFrame({
            "mass_GeV": [0.5, 1.0, 1.5],
            "u2_min": [np.nan, np.nan, np.nan],
            "u2_max": [np.nan, np.nan, np.nan],
            "peak_N": [0.1, 0.5, 0.2],
            "peak_u2": [1e-6, 1e-5, 1e-6],
            "has_sensitivity": [False, False, False],
        })
        result = self._build_island_polygon(df)
        self.assertIsNone(result)

    def test_single_sensitive_point(self):
        """Single point with sensitivity should return a valid polygon or None."""
        df = pd.DataFrame({
            "mass_GeV": [0.5, 1.0, 1.5],
            "u2_min": [np.nan, 1e-8, np.nan],
            "u2_max": [np.nan, 1e-5, np.nan],
            "peak_N": [0.1, 10.0, 0.2],
            "peak_u2": [1e-6, 1e-6, 1e-6],
            "has_sensitivity": [False, True, False],
            "flavor": ["Umu", "Umu", "Umu"],
        })
        result = self._build_island_polygon(df)
        # With only 1 point, savgol_filter may fail (window > length)
        # The function should handle this gracefully
        if result is not None:
            m_poly, u2_poly = result
            self.assertEqual(len(m_poly), len(u2_poly))
            self.assertGreater(len(m_poly), 0)

    def test_normal_island(self):
        """Standard multi-point island produces a closed polygon."""
        masses = np.linspace(0.5, 3.0, 20)
        df = pd.DataFrame({
            "mass_GeV": masses,
            "u2_min": np.full(20, 1e-9),
            "u2_max": np.full(20, 1e-5),
            "peak_N": np.full(20, 50.0),
            "peak_u2": np.full(20, 1e-7),
            "has_sensitivity": np.ones(20, dtype=bool),
            "flavor": ["Umu"] * 20,
        })
        result = self._build_island_polygon(df)
        self.assertIsNotNone(result)
        m_poly, u2_poly = result
        # Polygon should be closed (first = last not required by code, but
        # it should have entries for both upper and lower contours)
        self.assertGreater(len(m_poly), 20)  # upper + lower + tips

    def test_gap_in_island(self):
        """Island with a gap (non-sensitive mass in the middle)."""
        masses = np.linspace(0.5, 3.0, 20)
        has_sens = np.ones(20, dtype=bool)
        has_sens[9:12] = False  # gap in the middle
        u2_min = np.where(has_sens, 1e-9, np.nan)
        u2_max = np.where(has_sens, 1e-5, np.nan)
        df = pd.DataFrame({
            "mass_GeV": masses,
            "u2_min": u2_min,
            "u2_max": u2_max,
            "peak_N": np.where(has_sens, 50.0, 1.0),
            "peak_u2": np.full(20, 1e-7),
            "has_sensitivity": has_sens,
            "flavor": ["Umu"] * 20,
        })
        # _build_island_polygon uses only the valid subset, so it
        # will connect across the gap (limitation, but shouldn't crash)
        result = self._build_island_polygon(df)
        # Should produce some polygon (possibly covering the gap)
        if result is not None:
            m_poly, u2_poly = result
            self.assertGreater(len(m_poly), 0)
            self.assertTrue(np.all(np.isfinite(m_poly)))
            self.assertTrue(np.all(np.isfinite(u2_poly)))


if __name__ == "__main__":
    unittest.main()
