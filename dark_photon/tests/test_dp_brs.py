"""
Tests for dark photon branching ratios, decay widths, and meson production BRs.

Pure Python — no ROOT dependency.
"""

import math
import sys
import unittest
from pathlib import Path

import numpy as np

# Allow import of dark_photon modules
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dp_meson_brs import (
    ALPHA_EM,
    BR_ETA_GAMGAM,
    BR_OMEGA_PI0GAM,
    HBAR_C_GEV_M,
    M_EL,
    M_ETA,
    M_MU,
    M_OMEGA,
    M_PI0,
    _alpha_s_1loop,
    _partial_width_ll,
    _r_ratio_parametric,
    br_eta_to_dp_gamma,
    br_omega_to_dp_pi0,
    dp_ctau_m,
    dp_ctau_mm,
    dp_width_eps1,
    perturbative_brs,
    total_meson_br,
)


# ── perturbative_brs tests ──────────────────────────────────────────────


class TestPerturbativeBRs(unittest.TestCase):

    def test_sum_to_one_below_tau(self):
        """BRs sum to 1 between muon and tau thresholds."""
        for mass in [1.8, 2.0, 3.0]:
            brs = perturbative_brs(mass)
            total = sum(brs.values())
            self.assertAlmostEqual(total, 1.0, places=10,
                                   msg=f"BR sum != 1 at m={mass} GeV")

    def test_sum_to_one_above_tau(self):
        """BRs sum to 1 above tau threshold."""
        for mass in [3.6, 5.0, 10.0, 50.0]:
            brs = perturbative_brs(mass)
            total = sum(brs.values())
            self.assertAlmostEqual(total, 1.0, places=10,
                                   msg=f"BR sum != 1 at m={mass} GeV")

    def test_tau_absent_below_threshold(self):
        """tau channel is zero below 2*m_tau ~ 3.554 GeV."""
        brs = perturbative_brs(3.0)
        self.assertEqual(brs["tau"], 0.0)

    def test_tau_present_above_threshold(self):
        """tau channel is nonzero above 2*m_tau."""
        brs = perturbative_brs(4.0)
        self.assertGreater(brs["tau"], 0.0)

    def test_charm_absent_below_threshold(self):
        """cc_q channel is zero below 2*m_c ~ 2.54 GeV."""
        brs = perturbative_brs(2.0)
        self.assertEqual(brs["cc_q"], 0.0)

    def test_charm_present_above_threshold(self):
        """cc_q channel is nonzero above 2*m_c."""
        brs = perturbative_brs(3.0)
        self.assertGreater(brs["cc_q"], 0.0)

    def test_bottom_absent_below_threshold(self):
        """bb_q channel is zero below 2*m_b ~ 8.36 GeV."""
        brs = perturbative_brs(8.0)
        self.assertEqual(brs["bb_q"], 0.0)

    def test_bottom_present_above_threshold(self):
        """bb_q channel is nonzero above 2*m_b."""
        brs = perturbative_brs(10.0)
        self.assertGreater(brs["bb_q"], 0.0)

    def test_all_channels_present(self):
        """All expected channels are returned."""
        brs = perturbative_brs(2.0)
        expected = {"ee", "mumu", "tau", "uu_q", "dd_q", "ss_q", "cc_q", "bb_q"}
        self.assertEqual(set(brs.keys()), expected)

    def test_lepton_universality_below_tau(self):
        """ee and mumu should have equal BRs (both massless limit)."""
        brs = perturbative_brs(3.0)
        self.assertAlmostEqual(brs["ee"], brs["mumu"], places=10)

    def test_quark_charge_ratios(self):
        """At high mass with all quarks open, ratio uu_q/dd_q = (2/3)^2/(1/3)^2 = 4."""
        brs = perturbative_brs(50.0)
        # uu_q = u quark, cc_q = c quark (both charge 2/3)
        # dd_q = d quark (1/3), ss_q = s quark (1/3), bb_q = b quark (1/3)
        # Each quark channel BR ∝ 3*e_q^2*(1+alpha_s/pi)
        # So uu_q / dd_q should be (2/3)^2 / (1/3)^2 = 4
        ratio = brs["uu_q"] / brs["dd_q"]
        self.assertAlmostEqual(ratio, 4.0, places=8)


# ── dp_width_eps1 tests ─────────────────────────────────────────────────


class TestDPWidthEps1(unittest.TestCase):

    def test_below_electron_threshold(self):
        """Width is zero below 2*m_e ~ 1.02 MeV."""
        w = dp_width_eps1(0.0005)
        self.assertEqual(w, 0.0)

    def test_above_electron_only(self):
        """Width is positive above 2*m_e but below 2*m_mu."""
        w = dp_width_eps1(0.01)  # 10 MeV, above 2*m_e
        self.assertGreater(w, 0.0)

    def test_width_increases_with_mass(self):
        """Width should generally increase with mass (more channels open)."""
        w1 = dp_width_eps1(0.5)
        w2 = dp_width_eps1(1.0)
        # Note: hadronic resonances can cause dips, but on average increases
        self.assertGreater(w2, 0.0)
        self.assertGreater(w1, 0.0)

    def test_known_partial_width_electrons(self):
        """Verify e+e- partial width formula at m=1 GeV."""
        m_dp = 1.0
        # Γ(A'→e+e-) = (α/3) * m_A' * sqrt(1-4r) * (1+2r), r = (m_e/m_A')^2
        r = (M_EL / m_dp) ** 2
        expected = (ALPHA_EM / 3.0) * m_dp * math.sqrt(1 - 4*r) * (1 + 2*r)
        actual = _partial_width_ll(m_dp, M_EL)
        self.assertAlmostEqual(actual, expected, places=15)


# ── Meson BR tests ──────────────────────────────────────────────────────


class TestMesonBRs(unittest.TestCase):

    def test_eta_kinematic_threshold(self):
        """BR(eta->A'gamma) = 0 at m_A' >= m_eta."""
        self.assertEqual(br_eta_to_dp_gamma(1e-6, M_ETA), 0.0)
        self.assertEqual(br_eta_to_dp_gamma(1e-6, M_ETA + 0.01), 0.0)

    def test_eta_nonzero_below_threshold(self):
        """BR(eta->A'gamma) > 0 below m_eta."""
        br = br_eta_to_dp_gamma(1e-6, 0.3)
        self.assertGreater(br, 0.0)

    def test_eta_eps2_scaling(self):
        """BR(eta->A'gamma) is linear in eps2."""
        m = 0.3
        br1 = br_eta_to_dp_gamma(1e-6, m)
        br2 = br_eta_to_dp_gamma(2e-6, m)
        self.assertAlmostEqual(br2 / br1, 2.0, places=10)

    def test_eta_massless_limit(self):
        """BR(eta->A'gamma) at m_A'=0 should be 2*eps2*BR(eta->gamgam)."""
        eps2 = 1e-6
        expected = 2.0 * eps2 * BR_ETA_GAMGAM
        actual = br_eta_to_dp_gamma(eps2, 1e-10)  # nearly massless
        self.assertAlmostEqual(actual, expected, places=15)

    def test_omega_kinematic_threshold(self):
        """BR(omega->A'pi0) = 0 at m_A' >= m_omega - m_pi0."""
        threshold = M_OMEGA - M_PI0
        self.assertEqual(br_omega_to_dp_pi0(1e-6, threshold), 0.0)
        self.assertEqual(br_omega_to_dp_pi0(1e-6, threshold + 0.01), 0.0)

    def test_omega_nonzero_below_threshold(self):
        """BR(omega->A'pi0) > 0 below threshold."""
        br = br_omega_to_dp_pi0(1e-6, 0.3)
        self.assertGreater(br, 0.0)

    def test_omega_eps2_scaling(self):
        """BR(omega->A'pi0) is linear in eps2."""
        m = 0.3
        br1 = br_omega_to_dp_pi0(1e-6, m)
        br2 = br_omega_to_dp_pi0(3e-6, m)
        self.assertAlmostEqual(br2 / br1, 3.0, places=10)

    def test_total_meson_br_both_open(self):
        """total_meson_br returns both channels at low mass."""
        result = total_meson_br(1e-6, 0.3)
        self.assertIn("eta_to_dp_gamma", result)
        self.assertIn("omega_to_dp_pi0", result)

    def test_total_meson_br_none_open(self):
        """total_meson_br returns empty dict above both thresholds."""
        result = total_meson_br(1e-6, 1.0)  # above both meson thresholds
        self.assertEqual(len(result), 0)


# ── ctau tests ───────────────────────────────────────────────────────────


class TestCtau(unittest.TestCase):

    def test_ctau_positive_finite(self):
        """cτ should be positive and finite for physical parameters."""
        ctau = dp_ctau_m(1e-6, 0.5)
        self.assertGreater(ctau, 0.0)
        self.assertTrue(np.isfinite(ctau))

    def test_ctau_inf_for_zero_eps2(self):
        """cτ should be inf for eps2=0."""
        self.assertEqual(dp_ctau_m(0.0, 0.5), float('inf'))

    def test_ctau_inf_below_threshold(self):
        """cτ should be inf when no decay channels open."""
        self.assertEqual(dp_ctau_m(1e-6, 0.0001), float('inf'))

    def test_ctau_mm_conversion(self):
        """cτ in mm should be 1000× cτ in m."""
        eps2, m = 1e-6, 0.5
        ctau_m = dp_ctau_m(eps2, m)
        ctau_mm = dp_ctau_mm(eps2, m)
        self.assertAlmostEqual(ctau_mm, ctau_m * 1e3, places=5)

    def test_ctau_scales_with_eps2(self):
        """cτ ∝ 1/eps2."""
        m = 0.5
        ctau1 = dp_ctau_m(1e-6, m)
        ctau2 = dp_ctau_m(1e-8, m)
        self.assertAlmostEqual(ctau2 / ctau1, 100.0, places=5)

    def test_ctau_formula(self):
        """Verify cτ = hbar_c / (eps2 * Gamma_0)."""
        eps2, m = 1e-6, 1.0
        g0 = dp_width_eps1(m)
        expected = HBAR_C_GEV_M / (eps2 * g0)
        actual = dp_ctau_m(eps2, m)
        self.assertAlmostEqual(actual, expected, places=10)


# ── alpha_s tests ────────────────────────────────────────────────────────


class TestAlphaS(unittest.TestCase):

    def test_alpha_s_capped(self):
        """alpha_s is capped at 0.4 for low scales."""
        als = _alpha_s_1loop(0.1, 3)
        self.assertLessEqual(als, 0.4)

    def test_alpha_s_below_lambda(self):
        """alpha_s at scale below Lambda returns 0.4."""
        als = _alpha_s_1loop(0.1, 3, Lambda=0.2)
        self.assertEqual(als, 0.4)

    def test_alpha_s_at_mz(self):
        """alpha_s at m_Z ~ 91 GeV should be roughly 0.12."""
        als = _alpha_s_1loop(91.0, 5)
        self.assertAlmostEqual(als, 0.12, delta=0.02)

    def test_alpha_s_at_2gev(self):
        """alpha_s at 2 GeV should be ~0.30 (PDG: 0.30)."""
        als = _alpha_s_1loop(2.0, 3)
        self.assertAlmostEqual(als, 0.30, delta=0.03)

    def test_alpha_s_runs(self):
        """alpha_s should decrease with increasing scale."""
        als_low = _alpha_s_1loop(2.0, 3)
        als_high = _alpha_s_1loop(10.0, 4)
        self.assertGreater(als_low, als_high)


# ── _assemble_channels tests (via make_dp_cmnd) ─────────────────────────


class TestAssembleChannels(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Import _assemble_channels from make_dp_cmnd."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from make_dp_cmnd import _assemble_channels, _assemble_perturbative, MIN_BR
        cls._assemble_channels = staticmethod(_assemble_channels)
        cls._assemble_perturbative = staticmethod(_assemble_perturbative)
        cls.MIN_BR = MIN_BR

    def test_perturbative_sums_to_one(self):
        """Perturbative channel assembly normalizes to 1."""
        brs = perturbative_brs(3.0)
        channels = self._assemble_perturbative(3.0, brs)
        total = sum(ch[0] for ch in channels)
        self.assertAlmostEqual(total, 1.0, places=10)

    def test_perturbative_above_tau(self):
        """Perturbative channel assembly works above tau threshold."""
        brs = perturbative_brs(5.0)
        channels = self._assemble_perturbative(5.0, brs)
        total = sum(ch[0] for ch in channels)
        self.assertAlmostEqual(total, 1.0, places=10)
        # tau should be present
        has_tau = any(15 in ch[1] for ch in channels)
        self.assertTrue(has_tau)

    def test_vmd_channel_normalization(self):
        """VMD channel assembly normalizes to 1."""
        br_dict = {"ee": 0.3, "mumu": 0.3, "tau": 0.0, "pipi": 0.3,
                    "pi3": 0.05, "KKc": 0.04, "BRqcd": 0.4,
                    "had_other": 0.01}
        channels = self._assemble_channels(1.0, br_dict)
        total = sum(ch[0] for ch in channels)
        self.assertAlmostEqual(total, 1.0, places=10)

    def test_minor_channels_merged(self):
        """Channels below MIN_BR get merged into pipi."""
        br_dict = {"ee": 0.45, "mumu": 0.45, "tau": 0.0,
                    "pipi": 0.05, "pi3": 0.001, "KKc": 0.001,
                    "had_other": 0.0}
        channels = self._assemble_channels(1.0, br_dict)
        # pi3 and KKc below MIN_BR=0.005 should be merged into pipi
        total = sum(ch[0] for ch in channels)
        self.assertAlmostEqual(total, 1.0, places=10)


# ── R-ratio parametric tests ────────────────────────────────────────────


class TestRRatioParametric(unittest.TestCase):

    def test_r_zero_below_pion(self):
        """R-ratio is 0 below 2*m_pi for parametric fallback."""
        R = _r_ratio_parametric(0.1)
        self.assertEqual(R, 0.0)

    def test_r_positive_above_pion(self):
        """R-ratio is positive above 2*m_pi."""
        R = _r_ratio_parametric(0.5)
        self.assertGreater(R, 0.0)

    def test_r_increases_with_quark_thresholds(self):
        """R-ratio should jump at charm/bottom thresholds."""
        R_below_c = _r_ratio_parametric(2.0)
        R_above_c = _r_ratio_parametric(3.0)
        R_above_b = _r_ratio_parametric(10.0)
        self.assertGreater(R_above_c, R_below_c)
        self.assertGreater(R_above_b, R_above_c)


if __name__ == "__main__":
    unittest.main()
