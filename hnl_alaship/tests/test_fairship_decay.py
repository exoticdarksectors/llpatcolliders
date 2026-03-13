"""
Tests for the FairShip decay backend.

Requires ROOT with Pythia8 support. Tests are skipped if ROOT is not available.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import ROOT  # noqa: F401
    HAS_ROOT = True
except ImportError:
    HAS_ROOT = False


@unittest.skipUnless(HAS_ROOT, "ROOT not available")
class TestFairShipDecay(unittest.TestCase):

    def test_ctau_u2eq1_reasonable(self):
        """ctau_u2eq1 returns a positive finite value at m_N=1.0 GeV."""
        from analysis.fairship_decay import ctau_u2eq1
        ctau = ctau_u2eq1(1.0, flavor="Umu")
        self.assertGreater(ctau, 0.0)
        self.assertTrue(np.isfinite(ctau))
        # At 1 GeV with U²=1, FairShip gives ctau ~ O(1e-4) m
        self.assertLess(ctau, 1e-1)

    def test_ctau_flavor_dependence(self):
        """Different flavors give different ctau values."""
        from analysis.fairship_decay import ctau_u2eq1
        ctau_umu = ctau_u2eq1(1.0, flavor="Umu")
        ctau_ue = ctau_u2eq1(1.0, flavor="Ue")
        # These should differ (different coupling structure)
        self.assertNotAlmostEqual(ctau_umu, ctau_ue, places=15)

    def test_boost_decay_to_lab_preserves_4momentum(self):
        """Boosted daughters should conserve 4-momentum."""
        from analysis.fairship_decay import boost_decay_to_lab, DecayTemplate

        # Create a simple mock rest-frame decay: parent -> 2 equal daughters
        m_N = 1.0
        m_d = 0.1  # daughter mass
        p_star = np.sqrt((m_N / 2.0) ** 2 - m_d ** 2)
        template = DecayTemplate(
            pdg=np.array([13, -13], dtype=np.int32),
            px=np.array([p_star, -p_star]),
            py=np.array([0.0, 0.0]),
            pz=np.array([0.0, 0.0]),
            energy=np.array([m_N / 2.0, m_N / 2.0]),
            mass=np.array([m_d, m_d]),
            charge=np.array([-1.0, 1.0]),
            stable=np.array([True, True]),
            status=np.array([1, 1], dtype=np.int32),
        )

        # Boosted parent with invariant mass = m_N = 1.0
        # E=10, px=1, py=2, pz=sqrt(E²-px²-py²-m²) = sqrt(94)
        parent_p4 = np.array([10.0, 1.0, 2.0, np.sqrt(94.0)])
        result = boost_decay_to_lab(parent_p4, template)

        # Sum of daughter momenta should equal parent momentum
        total_E = result["energy"].sum()
        total_px = result["px"].sum()
        total_py = result["py"].sum()
        total_pz = result["pz"].sum()

        # The parent invariant mass should be m_N
        parent_mass2 = parent_p4[0]**2 - np.sum(parent_p4[1:]**2)
        daughter_mass2 = total_E**2 - (total_px**2 + total_py**2 + total_pz**2)

        np.testing.assert_allclose(daughter_mass2, parent_mass2, rtol=1e-8)

    def test_boost_at_rest(self):
        """Boost with parent at rest should leave daughters unchanged."""
        from analysis.fairship_decay import boost_decay_to_lab, DecayTemplate

        m_N = 1.0
        template = DecayTemplate(
            pdg=np.array([211, -211], dtype=np.int32),
            px=np.array([0.3, -0.3]),
            py=np.array([0.1, -0.1]),
            pz=np.array([0.2, -0.2]),
            energy=np.array([0.5, 0.5]),
            mass=np.array([0.1396, 0.1396]),
            charge=np.array([1.0, -1.0]),
            stable=np.array([True, True]),
            status=np.array([1, 1], dtype=np.int32),
        )

        parent_p4 = np.array([m_N, 0.0, 0.0, 0.0])
        result = boost_decay_to_lab(parent_p4, template)

        np.testing.assert_allclose(result["px"], template.px, atol=1e-12)
        np.testing.assert_allclose(result["py"], template.py, atol=1e-12)
        np.testing.assert_allclose(result["pz"], template.pz, atol=1e-12)


@unittest.skipUnless(HAS_ROOT, "ROOT not available")
class TestFairShipDecayBackend(unittest.TestCase):

    def test_sample_rest_frame_decay_produces_daughters(self):
        """FairShipDecayBackend produces non-empty decay templates."""
        from analysis.fairship_decay import FairShipDecayBackend
        backend = FairShipDecayBackend(mass_GeV=1.0, flavor="Umu", random_seed=42)
        template = backend.sample_rest_frame_decay()
        self.assertGreater(len(template.pdg), 0)
        self.assertTrue(np.isfinite(template.energy).all())

    def test_sample_multiple_decays(self):
        """Sampling multiple decays returns the right count."""
        from analysis.fairship_decay import sample_rest_frame_decays
        templates = sample_rest_frame_decays(1.0, n_templates=5, flavor="Umu", random_seed=99)
        self.assertEqual(len(templates), 5)
        for t in templates:
            self.assertGreater(len(t.pdg), 0)


if __name__ == "__main__":
    unittest.main()
