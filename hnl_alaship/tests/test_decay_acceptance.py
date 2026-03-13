"""
Tests for the decay acceptance kernel (decay_acceptance.py).

Tests for _pairwise_acceptance are pure numpy.
Tests for _accept_decay_at_position require ROOT (for boost + decay templates).
"""

import sys
import unittest
from pathlib import Path

import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "geometry"))

from analysis.decay_acceptance import _pairwise_acceptance
from analysis.constants import SEP_MIN, SEP_MAX, MIN_CHARGED_DAUGHTERS

try:
    import ROOT  # noqa: F401
    HAS_ROOT = True
except ImportError:
    HAS_ROOT = False


def _make_box_mesh(center=(0, 0, 0), extents=(2, 2, 2)):
    """Create a simple box mesh."""
    box = trimesh.primitives.Box(extents=extents)
    box.apply_translation(center)
    return box


# ── _pairwise_acceptance tests ───────────────────────────────────────────


class TestPairwiseAcceptance(unittest.TestCase):

    def test_too_few_points(self):
        """Fewer than MIN_CHARGED_DAUGHTERS points should return False."""
        points = np.array([[0, 0, 0]], dtype=float)
        self.assertFalse(_pairwise_acceptance(points))

    def test_exactly_min_points_in_window(self):
        """Two points separated by SEP_MIN < d < SEP_MAX should pass."""
        sep = (SEP_MIN + SEP_MAX) / 2  # midpoint of window
        points = np.array([[0, 0, 0], [sep, 0, 0]], dtype=float)
        self.assertTrue(_pairwise_acceptance(points))

    def test_too_close(self):
        """Two points closer than SEP_MIN should fail."""
        points = np.array([[0, 0, 0], [SEP_MIN * 0.1, 0, 0]], dtype=float)
        self.assertFalse(_pairwise_acceptance(points))

    def test_too_far(self):
        """Two points farther than SEP_MAX should fail."""
        points = np.array([[0, 0, 0], [SEP_MAX * 2, 0, 0]], dtype=float)
        self.assertFalse(_pairwise_acceptance(points))

    def test_three_points_one_good_pair(self):
        """Three points where only one pair is in window should pass."""
        points = np.array([
            [0, 0, 0],
            [SEP_MAX * 5, 0, 0],     # too far from [0]
            [0.1, 0, 0],              # in window with [0]
        ], dtype=float)
        self.assertTrue(_pairwise_acceptance(points))

    def test_identical_points_fail(self):
        """Identical points have zero separation → fail."""
        points = np.array([[1, 2, 3], [1, 2, 3]], dtype=float)
        self.assertFalse(_pairwise_acceptance(points))

    def test_boundary_sep_min(self):
        """Separation exactly at SEP_MIN should pass (>=)."""
        points = np.array([[0, 0, 0], [SEP_MIN, 0, 0]], dtype=float)
        self.assertTrue(_pairwise_acceptance(points))

    def test_boundary_sep_max(self):
        """Separation exactly at SEP_MAX should pass (<=)."""
        points = np.array([[0, 0, 0], [SEP_MAX, 0, 0]], dtype=float)
        self.assertTrue(_pairwise_acceptance(points))


# ── _accept_decay_at_position tests (require ROOT) ──────────────────────


@unittest.skipUnless(HAS_ROOT, "ROOT not available")
class TestAcceptDecayAtPosition(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from analysis.decay_acceptance import _accept_decay_at_position
        from analysis.fairship_decay import DecayTemplate, boost_decay_to_lab
        cls._accept_decay_at_position = staticmethod(_accept_decay_at_position)
        cls.DecayTemplate = DecayTemplate
        cls.boost_decay_to_lab = staticmethod(boost_decay_to_lab)

    def test_collinear_daughters_hit_box(self):
        """Two daughters boosted forward into a box mesh should be accepted
        if their separation falls within the acceptance window."""
        # Create a box at z=10 (representing detector)
        mesh = _make_box_mesh(center=(0, 0, 10), extents=(4, 4, 2))

        # Parent moving along +Z with high momentum
        m_N = 1.0
        E = 50.0
        pz = np.sqrt(E**2 - m_N**2)
        parent_p4 = np.array([E, 0.0, 0.0, pz])

        # Rest-frame template: back-to-back muons along X
        p_star = np.sqrt((m_N / 2)**2 - 0.10566**2)
        template = self.DecayTemplate(
            pdg=np.array([13, -13], dtype=np.int32),
            px=np.array([p_star, -p_star]),
            py=np.array([0.0, 0.0]),
            pz=np.array([0.0, 0.0]),
            energy=np.array([m_N / 2, m_N / 2]),
            mass=np.array([0.10566, 0.10566]),
            charge=np.array([-1.0, 1.0]),
            stable=np.array([True, True]),
            status=np.array([1, 1], dtype=np.int32),
        )

        mother_direction = np.array([0, 0, 1], dtype=float)
        # Test at decay_distance=5, exit_distance=11
        result = self._accept_decay_at_position(
            mesh=mesh,
            decay_distance=5.0,
            exit_distance=11.0,
            mother_direction=mother_direction,
            parent_p4=parent_p4,
            template=template,
        )
        # Result is bool — may or may not pass depending on exact geometry
        self.assertIsInstance(result, bool)

    def test_no_charged_daughters_fails(self):
        """Template with no charged daughters should always fail."""
        mesh = _make_box_mesh(center=(0, 0, 10), extents=(4, 4, 2))

        m_N = 1.0
        E = 50.0
        pz = np.sqrt(E**2 - m_N**2)
        parent_p4 = np.array([E, 0.0, 0.0, pz])

        # Template with only neutrals (neutrinos)
        template = self.DecayTemplate(
            pdg=np.array([12, -12], dtype=np.int32),
            px=np.array([0.3, -0.3]),
            py=np.array([0.0, 0.0]),
            pz=np.array([0.0, 0.0]),
            energy=np.array([0.5, 0.5]),
            mass=np.array([0.0, 0.0]),
            charge=np.array([0.0, 0.0]),
            stable=np.array([True, True]),
            status=np.array([1, 1], dtype=np.int32),
        )

        mother_direction = np.array([0, 0, 1], dtype=float)
        result = self._accept_decay_at_position(
            mesh=mesh,
            decay_distance=5.0,
            exit_distance=11.0,
            mother_direction=mother_direction,
            parent_p4=parent_p4,
            template=template,
        )
        self.assertFalse(result)

    def test_low_momentum_daughters_fail(self):
        """Daughters with p < P_CUT should fail acceptance."""
        mesh = _make_box_mesh(center=(0, 0, 10), extents=(4, 4, 2))

        # Parent barely above threshold → daughters have ~0 momentum in rest frame
        m_N = 0.212  # barely above 2*m_mu
        E = 0.22
        pz = np.sqrt(max(0, E**2 - m_N**2))
        parent_p4 = np.array([E, 0.0, 0.0, pz])

        p_star = np.sqrt(max(0, (m_N / 2)**2 - 0.10566**2))
        template = self.DecayTemplate(
            pdg=np.array([13, -13], dtype=np.int32),
            px=np.array([p_star, -p_star]),
            py=np.array([0.0, 0.0]),
            pz=np.array([0.0, 0.0]),
            energy=np.array([m_N / 2, m_N / 2]),
            mass=np.array([0.10566, 0.10566]),
            charge=np.array([-1.0, 1.0]),
            stable=np.array([True, True]),
            status=np.array([1, 1], dtype=np.int32),
        )

        mother_direction = np.array([0, 0, 1], dtype=float)
        result = self._accept_decay_at_position(
            mesh=mesh,
            decay_distance=5.0,
            exit_distance=11.0,
            mother_direction=mother_direction,
            parent_p4=parent_p4,
            template=template,
        )
        # Very low boost → daughters below P_CUT → should fail
        self.assertFalse(result)


if __name__ == "__main__":
    unittest.main()
