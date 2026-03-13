"""Tests for analysis.format_bridge — 4-vector to geometry-ready conversion."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import unittest
import tempfile
import numpy as np

from analysis.format_bridge import load_combined_csv


class TestFormatBridge(unittest.TestCase):

    def test_known_4vector_along_z(self):
        """Particle along +Z axis: eta -> large positive, phi ~ 0."""
        m_N = 1.0
        E, px, py, pz = 10.0, 0.0, 0.0, np.sqrt(10.0**2 - 1.0**2)
        weight = 1e-5

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(f"{weight},{E},{px},{py},{pz}\n")
            f.flush()
            data = load_combined_csv(f.name, m_N)

        self.assertEqual(len(data["weight"]), 1)
        self.assertAlmostEqual(data["weight"][0], weight, places=10)

        # eta should be large positive (forward)
        self.assertGreater(data["eta"][0], 2.0)

        # p should match
        expected_p = np.sqrt(px**2 + py**2 + pz**2)
        np.testing.assert_allclose(data["p"][0], expected_p, rtol=1e-10)

        # beta_gamma = p / m_N
        np.testing.assert_allclose(
            data["beta_gamma"][0], expected_p / m_N, rtol=1e-10)

    def test_known_4vector_transverse(self):
        """Particle in +Y direction: eta ~ 0, phi ~ pi/2."""
        m_N = 1.0
        E, px, py, pz = 5.0, 0.0, np.sqrt(5.0**2 - 1.0**2), 0.0
        weight = 2e-6

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(f"{weight},{E},{px},{py},{pz}\n")
            f.flush()
            data = load_combined_csv(f.name, m_N)

        # eta should be near zero (transverse)
        self.assertAlmostEqual(data["eta"][0], 0.0, places=1)

        # phi should be near pi/2
        self.assertAlmostEqual(data["phi"][0], np.pi / 2.0, places=5)

    def test_multi_event(self):
        """Multiple events load correctly."""
        m_N = 0.5
        lines = [
            "1e-5,3.0,0.1,0.2,2.5\n",
            "2e-5,4.0,-0.3,0.1,-3.8\n",
            "3e-5,2.0,0.5,-0.1,1.5\n",
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.writelines(lines)
            f.flush()
            data = load_combined_csv(f.name, m_N)

        self.assertEqual(len(data["weight"]), 3)
        np.testing.assert_allclose(
            data["weight"], [1e-5, 2e-5, 3e-5], rtol=1e-10)

        # All finite
        for key in ["eta", "phi", "p", "gamma", "beta", "beta_gamma"]:
            self.assertTrue(np.all(np.isfinite(data[key])),
                            f"{key} has non-finite values")

    def test_empty_file(self):
        """Empty or missing file returns empty arrays."""
        data = load_combined_csv("/nonexistent/path.csv", 1.0)
        self.assertEqual(len(data["weight"]), 0)
        self.assertEqual(data["eta"].shape, (0,))

    def test_gamma_beta_consistency(self):
        """gamma and beta satisfy E = gamma * m, p = beta * E."""
        m_N = 1.5
        E, px, py, pz = 10.0, 1.0, 2.0, 9.0
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(f"1e-6,{E},{px},{py},{pz}\n")
            f.flush()
            data = load_combined_csv(f.name, m_N)

        # gamma = E / m_N
        np.testing.assert_allclose(data["gamma"][0], E / m_N, rtol=1e-10)
        # beta = p / E
        p = np.sqrt(px**2 + py**2 + pz**2)
        np.testing.assert_allclose(data["beta"][0], p / E, rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
