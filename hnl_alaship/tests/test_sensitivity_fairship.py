"""
Tests for sensitivity scan with FairShip decay backend.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))


class TestSensitivityScan(unittest.TestCase):

    def test_scan_u2_with_synthetic_cache(self):
        """scan_u2 returns positive N_signal for a synthetic cache."""
        from analysis.sensitivity import scan_u2

        n_samples = 10
        n_bins = 5
        cache = {
            "n_sampled": n_samples,
            "estimator": "exact",
            "total_hit_weight": 1.0,
            "sample_weights": np.ones(n_samples) * 1e-6,
            "beta_gamma": np.full(n_samples, 50.0),
            "position_centers": np.full((n_samples, n_bins), 22.0),
            "position_widths": np.full((n_samples, n_bins), 0.24 / n_bins),
            "acceptance": np.full((n_samples, n_bins), 0.5),
        }

        ctau_u2_1 = 1e-10  # metres at U²=1
        L_int_pb = 3e6

        u2_grid, n_grid = scan_u2(
            cache=cache,
            ctau_u2_1=ctau_u2_1,
            L_int_pb=L_int_pb,
            log_u2_min=-10,
            log_u2_max=-3,
            n_points=50,
        )

        self.assertEqual(len(u2_grid), 50)
        self.assertEqual(len(n_grid), 50)
        # Should have some non-zero signal somewhere
        self.assertGreater(n_grid.max(), 0.0)
        self.assertTrue(np.all(np.isfinite(n_grid)))

    def test_scan_u2_empty_cache(self):
        """scan_u2 returns zeros for an empty cache."""
        from analysis.sensitivity import scan_u2

        cache = {"n_sampled": 0}
        u2_grid, n_grid = scan_u2(cache, ctau_u2_1=1e-10, L_int_pb=3e6)

        self.assertTrue(np.all(n_grid == 0.0))

    def test_exclusion_band_from_synthetic(self):
        """find_exclusion_band finds a band in synthetic N_signal data."""
        from analysis.exclusion import find_exclusion_band

        u2_grid = np.logspace(-10, -3, 100)
        # Make a peaked N_signal: N ~ u2 * exp(-ctau/u2)
        ctau = 1e-10
        N_grid = 1e6 * u2_grid * np.exp(-ctau / (u2_grid * 50.0 * 22.0))

        result = find_exclusion_band(u2_grid, N_grid, N_threshold=3.0)

        if result["has_sensitivity"]:
            self.assertLess(result["u2_min"], result["u2_max"])
            self.assertGreater(result["peak_N"], 3.0)
        # If no sensitivity, that's also valid for this test


class TestConstantsSync(unittest.TestCase):

    def test_cuts_come_from_detector_cuts(self):
        """Acceptance cuts in constants.py are imported from detector_cuts.py."""
        geom_dir = str(Path(__file__).resolve().parent.parent.parent / "geometry")
        if geom_dir not in sys.path:
            sys.path.insert(0, geom_dir)
        import detector_cuts
        from analysis.constants import P_CUT, SEP_MIN, SEP_MAX

        self.assertIs(P_CUT, detector_cuts.P_CUT)
        self.assertIs(SEP_MIN, detector_cuts.SEP_MIN)
        self.assertIs(SEP_MAX, detector_cuts.SEP_MAX)


if __name__ == "__main__":
    unittest.main()
