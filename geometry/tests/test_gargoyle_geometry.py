"""
Tests for shared geometry module: ray-casting, direction conversion, decay length.

Uses simple box meshes for ray-casting tests — no ROOT dependency.
"""

import unittest

import numpy as np
import trimesh

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gargoyle_geometry import (
    SPEED_OF_LIGHT,
    calculate_decay_length,
    eta_phi_to_direction,
    momenta_to_directions,
    ray_intersection_distances,
    tunnel_profile_points,
)


def _make_box_mesh(center=(0, 0, 0), extents=(2, 2, 2)):
    """Create a simple axis-aligned box mesh for testing."""
    box = trimesh.primitives.Box(extents=extents)
    box.apply_translation(center)
    return box


# ── momenta_to_directions tests ──────────────────────────────────────────


class TestMomentaToDirections(unittest.TestCase):

    def test_unit_vectors(self):
        """Unit vectors should be returned as-is."""
        momenta = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        dirs, valid = momenta_to_directions(momenta)
        np.testing.assert_allclose(dirs, momenta, atol=1e-15)
        self.assertTrue(np.all(valid))

    def test_zero_vector(self):
        """Zero-momentum row should be flagged as invalid."""
        momenta = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
        dirs, valid = momenta_to_directions(momenta)
        self.assertFalse(valid[0])
        self.assertTrue(valid[1])
        np.testing.assert_allclose(dirs[0], [0, 0, 0])

    def test_normalization(self):
        """Non-unit vectors should be normalized."""
        momenta = np.array([[3, 4, 0], [0, 0, 5]], dtype=float)
        dirs, valid = momenta_to_directions(momenta)
        norms = np.linalg.norm(dirs[valid], axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-14)

    def test_batch_shape(self):
        """Output shapes should match input."""
        n = 10
        momenta = np.random.default_rng(42).normal(size=(n, 3))
        dirs, valid = momenta_to_directions(momenta)
        self.assertEqual(dirs.shape, (n, 3))
        self.assertEqual(valid.shape, (n,))

    def test_wrong_shape_raises(self):
        """Non-(N,3) input should raise ValueError."""
        with self.assertRaises(ValueError):
            momenta_to_directions(np.array([1, 2, 3]))
        with self.assertRaises(ValueError):
            momenta_to_directions(np.array([[1, 2]]))


# ── ray_intersection_distances tests ─────────────────────────────────────


class TestRayIntersectionDistances(unittest.TestCase):

    def setUp(self):
        """Create a 2x2x2 box centered at origin."""
        self.mesh = _make_box_mesh(center=(0, 0, 0), extents=(2, 2, 2))

    def test_ray_through_box_center(self):
        """Ray along +X from outside should hit both faces."""
        origins = np.array([[-5, 0, 0]], dtype=float)
        directions = np.array([[1, 0, 0]], dtype=float)
        hits, first_d, second_d = ray_intersection_distances(
            self.mesh, origins, directions)
        self.assertTrue(hits[0])
        # Entry at x=-1, exit at x=+1 → distances 4 and 6 from origin at x=-5
        np.testing.assert_allclose(first_d[0], 4.0, atol=0.05)
        np.testing.assert_allclose(second_d[0], 6.0, atol=0.05)

    def test_ray_misses_box(self):
        """Ray that misses the box should return no hit."""
        origins = np.array([[-5, 5, 0]], dtype=float)
        directions = np.array([[1, 0, 0]], dtype=float)
        hits, first_d, second_d = ray_intersection_distances(
            self.mesh, origins, directions)
        self.assertFalse(hits[0])
        self.assertTrue(np.isnan(first_d[0]))
        self.assertTrue(np.isnan(second_d[0]))

    def test_multiple_rays(self):
        """Batch of rays: some hit, some miss."""
        origins = np.array([
            [-5, 0, 0],   # hits
            [0, -5, 0],   # hits (along +Y)
            [-5, 5, 0],   # misses
        ], dtype=float)
        directions = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [1, 0, 0],
        ], dtype=float)
        hits, first_d, second_d = ray_intersection_distances(
            self.mesh, origins, directions)
        self.assertTrue(hits[0])
        self.assertTrue(hits[1])
        self.assertFalse(hits[2])

    def test_empty_rays(self):
        """Empty input should return empty output."""
        origins = np.empty((0, 3), dtype=float)
        directions = np.empty((0, 3), dtype=float)
        hits, first_d, second_d = ray_intersection_distances(
            self.mesh, origins, directions)
        self.assertEqual(len(hits), 0)

    def test_entry_before_exit(self):
        """First intersection distance should be less than second."""
        origins = np.array([[-5, 0, 0]], dtype=float)
        directions = np.array([[1, 0, 0]], dtype=float)
        hits, first_d, second_d = ray_intersection_distances(
            self.mesh, origins, directions)
        self.assertTrue(hits[0])
        self.assertLess(first_d[0], second_d[0])

    def test_diagonal_ray(self):
        """Ray along diagonal should also register two hits."""
        origins = np.array([[-5, -5, -5]], dtype=float)
        d = np.array([1, 1, 1], dtype=float)
        d /= np.linalg.norm(d)
        directions = d[None, :]
        hits, first_d, second_d = ray_intersection_distances(
            self.mesh, origins, directions)
        self.assertTrue(hits[0])
        self.assertTrue(np.isfinite(first_d[0]))
        self.assertTrue(np.isfinite(second_d[0]))

    def test_shape_mismatch_raises(self):
        """Mismatched origin/direction shapes should raise ValueError."""
        with self.assertRaises(ValueError):
            ray_intersection_distances(
                self.mesh,
                np.array([[0, 0, 0]]),
                np.array([[1, 0, 0], [0, 1, 0]]))

    def test_offset_box(self):
        """Ray-cast works with off-origin box."""
        mesh = _make_box_mesh(center=(10, 20, 0), extents=(2, 2, 2))
        origins = np.array([[10, 20, -5]], dtype=float)
        directions = np.array([[0, 0, 1]], dtype=float)
        hits, first_d, second_d = ray_intersection_distances(
            mesh, origins, directions)
        self.assertTrue(hits[0])
        np.testing.assert_allclose(first_d[0], 4.0, atol=0.05)
        np.testing.assert_allclose(second_d[0], 6.0, atol=0.05)


# ── calculate_decay_length tests ─────────────────────────────────────────


class TestCalculateDecayLength(unittest.TestCase):

    def test_known_value(self):
        """Verify βγcτ formula with known inputs."""
        momentum = 10.0  # GeV
        mass = 1.0       # GeV
        lifetime = 1e-9   # seconds
        energy = np.sqrt(momentum**2 + mass**2)
        beta = momentum / energy
        gamma = energy / mass
        expected = gamma * beta * SPEED_OF_LIGHT * lifetime
        actual = calculate_decay_length(momentum, mass, lifetime)
        self.assertAlmostEqual(actual, expected, places=5)

    def test_zero_momentum(self):
        """At rest, decay length should be zero (beta=0)."""
        dl = calculate_decay_length(0.0, 1.0, 1e-9)
        self.assertAlmostEqual(dl, 0.0, places=15)

    def test_proportional_to_lifetime(self):
        """Decay length should scale linearly with lifetime."""
        dl1 = calculate_decay_length(10.0, 1.0, 1e-9)
        dl2 = calculate_decay_length(10.0, 1.0, 2e-9)
        self.assertAlmostEqual(dl2 / dl1, 2.0, places=10)

    def test_ultrarelativistic(self):
        """For p >> m, decay length ≈ p/m * c * τ."""
        p, m, tau = 1000.0, 0.1, 1e-12
        dl = calculate_decay_length(p, m, tau)
        approx = (p / m) * SPEED_OF_LIGHT * tau
        self.assertAlmostEqual(dl / approx, 1.0, delta=0.001)


# ── eta_phi_to_direction tests ───────────────────────────────────────────


class TestEtaPhiToDirection(unittest.TestCase):

    def test_eta_zero_phi_zero(self):
        """eta=0, phi=0 should point along +X (transverse, theta=pi/2)."""
        d = eta_phi_to_direction(0.0, 0.0)
        np.testing.assert_allclose(d, [1, 0, 0], atol=1e-10)

    def test_eta_zero_phi_half_pi(self):
        """eta=0, phi=pi/2 should point along +Y."""
        d = eta_phi_to_direction(0.0, np.pi / 2)
        np.testing.assert_allclose(d, [0, 1, 0], atol=1e-10)

    def test_large_positive_eta(self):
        """Large positive eta should point near +Z (beam direction)."""
        d = eta_phi_to_direction(5.0, 0.0)
        self.assertGreater(d[2], 0.99)

    def test_large_negative_eta(self):
        """Large negative eta should point near -Z."""
        d = eta_phi_to_direction(-5.0, 0.0)
        self.assertLess(d[2], -0.99)

    def test_unit_norm(self):
        """Direction should always be a unit vector."""
        for eta, phi in [(0, 0), (1.5, 0.7), (-2, 3.14), (4, 0)]:
            d = eta_phi_to_direction(eta, phi)
            self.assertAlmostEqual(np.linalg.norm(d), 1.0, places=10)


# ── tunnel_profile_points tests ──────────────────────────────────────────


class TestTunnelProfile(unittest.TestCase):

    def test_profile_is_2d(self):
        """Profile points should have shape (N, 2)."""
        pts = tunnel_profile_points()
        self.assertEqual(pts.ndim, 2)
        self.assertEqual(pts.shape[1], 2)

    def test_inset_shrinks_profile(self):
        """Inset should produce a smaller cross-section."""
        outer = tunnel_profile_points(inset=0.0)
        inner = tunnel_profile_points(inset=0.24)
        # The bounding box of inner should be smaller
        outer_w = outer[:, 0].max() - outer[:, 0].min()
        inner_w = inner[:, 0].max() - inner[:, 0].min()
        self.assertLess(inner_w, outer_w)


if __name__ == "__main__":
    unittest.main()
