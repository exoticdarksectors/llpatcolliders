import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import unittest
import warnings
import numpy as np

from production.decay_engine.kinematics import decay_2body


class TestKinematicsRegressions(unittest.TestCase):
    def test_decay_2body_multievent_vectorized(self):
        rng = np.random.default_rng(42)
        n = 3
        parent_E = np.array([5.0, 5.2, 5.5])
        parent_px = np.array([0.8, -0.4, 0.2])
        parent_py = np.array([0.1, 0.3, -0.2])
        parent_pz = np.array([0.5, -0.7, 0.9])

        d1, d2 = decay_2body(
            parent_E, parent_px, parent_py, parent_pz,
            m_parent=5.0, m1=1.0, m2=1.0, rng=rng,
        )

        self.assertEqual(d1.shape, (n, 4))
        self.assertEqual(d2.shape, (n, 4))
        self.assertTrue(np.isfinite(d1).all())
        self.assertTrue(np.isfinite(d2).all())

    def test_decay_2body_rest_frame_no_runtime_warning(self):
        rng = np.random.default_rng(123)
        parent_E = np.array([5.0, 5.0])
        parent_px = np.array([0.0, 0.0])
        parent_py = np.array([0.0, 0.0])
        parent_pz = np.array([0.0, 0.0])

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            d1, d2 = decay_2body(
                parent_E, parent_px, parent_py, parent_pz,
                m_parent=5.0, m1=1.0, m2=1.0, rng=rng,
            )

        runtime_warnings = [w for w in caught if issubclass(w.category, RuntimeWarning)]
        self.assertEqual(len(runtime_warnings), 0)
        self.assertTrue(np.isfinite(d1).all())
        self.assertTrue(np.isfinite(d2).all())

    def test_decay_2body_single_event_still_works(self):
        rng = np.random.default_rng(7)
        parent_E = np.array([5.0])
        parent_px = np.array([0.3])
        parent_py = np.array([0.1])
        parent_pz = np.array([-0.4])

        d1, d2 = decay_2body(
            parent_E, parent_px, parent_py, parent_pz,
            m_parent=5.0, m1=1.0, m2=1.0, rng=rng,
        )

        self.assertEqual(d1.shape, (1, 4))
        self.assertEqual(d2.shape, (1, 4))
        self.assertTrue(np.isfinite(d1).all())
        self.assertTrue(np.isfinite(d2).all())

    def test_decay_2body_empty_batch(self):
        rng = np.random.default_rng(99)
        parent_E = np.array([], dtype=float)
        parent_px = np.array([], dtype=float)
        parent_py = np.array([], dtype=float)
        parent_pz = np.array([], dtype=float)

        d1, d2 = decay_2body(
            parent_E, parent_px, parent_py, parent_pz,
            m_parent=5.0, m1=1.0, m2=1.0, rng=rng,
        )

        self.assertEqual(d1.shape, (0, 4))
        self.assertEqual(d2.shape, (0, 4))


if __name__ == "__main__":
    unittest.main()
