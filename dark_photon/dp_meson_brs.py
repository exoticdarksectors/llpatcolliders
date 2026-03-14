#!/usr/bin/env python3
"""
dp_meson_brs.py  --  Dark photon production BRs from meson decay
                      and A' decay width / lifetime formulas.

Physics:
  Dark photon A' with kinetic mixing ε to the SM photon.
  Sub-GeV production at pp colliders via light-meson decays:
    - η  → A' γ    (pseudoscalar, m_A' < m_η ≈ 548 MeV)
    - ω  → A' π⁰   (vector,       m_A' < m_ω - m_π ≈ 648 MeV)

  A' decay width:
    Γ(A'→ℓ⁺ℓ⁻) = (α/3) ε² m_A' √(1 - 4m_ℓ²/m_A'²) (1 + 2m_ℓ²/m_A'²)
    Γ(A'→had)   = Γ(A'→μ⁺μ⁻) × R(m_A')   [R = σ(e⁺e⁻→had)/σ(e⁺e⁻→μ⁺μ⁻)]

References:
  - arXiv:2005.01515 (The Dark Photon, Fabbrichesi et al.)
  - arXiv:2403.04181 (η decay to dark photon in UPC)
  - arXiv:2011.05115 (SHiP dark photon sensitivity)
"""

import math
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.036      # fine structure constant
HBAR_C_GEV_M = 0.197326980e-15  # ℏc in GeV·m
HBAR_C_GEV_MM = HBAR_C_GEV_M * 1e3  # ℏc in GeV·mm

# ---------------------------------------------------------------------------
# Particle masses (GeV) — PDG 2024
# ---------------------------------------------------------------------------
from constants import M_EL, M_MU, M_PI, M_PI0, M_ETA, M_OMEGA

# ---------------------------------------------------------------------------
# PDG branching ratios
# ---------------------------------------------------------------------------
BR_ETA_GAMGAM = 0.3941       # BR(η → γγ)
BR_OMEGA_PI0GAM = 0.0828     # BR(ω → π⁰γ)


# ---------------------------------------------------------------------------
# Meson → A' production BRs
# ---------------------------------------------------------------------------

def br_eta_to_dp_gamma(eps2: float, m_dp: float) -> float:
    """
    BR(η → A'γ) for a dark photon with kinetic mixing ε.

    Formula (arXiv:2403.04181 Eq.10, arXiv:2005.01515 Sec.3):
        BR(η→A'γ) = 2 ε² BR(η→γγ) (1 - m_A'²/m_η²)³ |F(m_A'²)|²

    The transition form factor F ≈ 1 for m_A' << m_η (good to ~5% at 0.5 GeV).

    Parameters
    ----------
    eps2 : float
        Kinetic mixing squared (ε²).
    m_dp : float
        Dark photon mass in GeV.

    Returns
    -------
    float
        BR(η → A'γ).  Returns 0 if kinematically forbidden.
    """
    if m_dp >= M_ETA:
        return 0.0
    x = (m_dp / M_ETA) ** 2
    return 2.0 * eps2 * BR_ETA_GAMGAM * (1.0 - x) ** 3


def br_omega_to_dp_pi0(eps2: float, m_dp: float) -> float:
    """
    BR(ω → A'π⁰) for a dark photon with kinetic mixing ε.

    Formula (arXiv:2005.01515 Sec.3):
        BR(ω→A'π⁰) = ε² BR(ω→π⁰γ) (p_A'/p_γ)³

    where p_A' and p_γ are the daughter momenta in the ω rest frame:
        p_γ  = (m_ω² - m_π⁰²) / (2 m_ω)           [massless photon limit]
        p_A' = λ^{1/2}(m_ω², m_π⁰², m_A'²) / (2 m_ω)

    Transition form factor |F|² ≈ 1 (valid for m_A' << m_ω).

    Parameters
    ----------
    eps2 : float
        Kinetic mixing squared (ε²).
    m_dp : float
        Dark photon mass in GeV.

    Returns
    -------
    float
        BR(ω → A'π⁰).  Returns 0 if kinematically forbidden.
    """
    if m_dp >= M_OMEGA - M_PI0:
        return 0.0
    # Daughter momentum for ω → π⁰ γ (photon massless)
    p_gamma = (M_OMEGA ** 2 - M_PI0 ** 2) / (2.0 * M_OMEGA)
    # Daughter momentum for ω → π⁰ A' (Källén function)
    s = M_OMEGA ** 2
    s1 = M_PI0 ** 2
    s2 = m_dp ** 2
    lam = s ** 2 + s1 ** 2 + s2 ** 2 - 2 * s * s1 - 2 * s * s2 - 2 * s1 * s2
    if lam <= 0:
        return 0.0
    p_dp = math.sqrt(lam) / (2.0 * M_OMEGA)
    return eps2 * BR_OMEGA_PI0GAM * (p_dp / p_gamma) ** 3


# ---------------------------------------------------------------------------
# A' decay width
# ---------------------------------------------------------------------------

def _partial_width_ll(m_dp: float, m_l: float) -> float:
    """
    Γ(A' → ℓ⁺ℓ⁻) at ε=1, in GeV.

    Formula:
        Γ = (α/3) m_A' √(1 - 4m_ℓ²/m_A'²) (1 + 2m_ℓ²/m_A'²)
    """
    if m_dp <= 2.0 * m_l:
        return 0.0
    r = (m_l / m_dp) ** 2
    return (ALPHA_EM / 3.0) * m_dp * math.sqrt(1.0 - 4.0 * r) * (1.0 + 2.0 * r)


def _r_ratio_from_table(m_dp: float) -> float:
    """
    R(s) = σ(e⁺e⁻→hadrons)/σ(e⁺e⁻→μ⁺μ⁻) at √s = m_A'.

    Uses the DeLiVeR BR table if available, otherwise falls back to a
    simple parametrization based on known resonance structure.

    For m < 2m_π: R = 0 (no hadronic channels)
    For 2m_π < m < 1.7 GeV: read from DeLiVeR table (BRqcd/BRmumu)
    Above 1.7 GeV: perturbative R = Nc × Σ e_q² × (1 + α_s/π)
    """
    if m_dp < 2.0 * M_PI:
        return 0.0

    # Try to load DeLiVeR table
    table_path = Path(__file__).parent / "br_tables" / "dp_brs_deliver.csv"
    if table_path.exists():
        return _r_ratio_from_deliver(m_dp, table_path)

    # Fallback: rough parametrization
    return _r_ratio_parametric(m_dp)


def _r_ratio_from_deliver(m_dp: float, table_path: Path) -> float:
    """Compute R from DeLiVeR BR table: R = (1 - BR_ee - BR_mumu - BR_tau) / BR_mumu."""
    cols: list[str] = []
    rows = []
    with open(table_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if not cols:
                cols = [c.strip() for c in line.split(",")]
                continue
            rows.append([float(x) for x in line.split(",")])
    data = np.array(rows)
    ci = {c: i for i, c in enumerate(cols)}

    masses = data[:, ci["mass_GeV"]]
    m_max = masses[-1]

    if m_dp > m_max:
        return _r_ratio_parametric(m_dp)

    br_ee = float(np.interp(m_dp, masses, data[:, ci["ee"]]))
    br_mumu = float(np.interp(m_dp, masses, data[:, ci["mumu"]]))
    br_tau = float(np.interp(m_dp, masses, data[:, ci.get("tau", ci.get("tautau", 0))]))
    if "tau" not in ci and "tautau" not in ci:
        br_tau = 0.0

    br_had = 1.0 - br_ee - br_mumu - br_tau
    if br_had < 0:
        br_had = 0.0
    if br_mumu <= 0:
        return 0.0
    return br_had / br_mumu


def _alpha_s_1loop(mu: float, nf: int, Lambda: float = 0.2) -> float:
    """1-loop running alpha_s, capped at 0.4 for low scales."""
    if mu <= Lambda:
        return 0.4
    als = 12 * math.pi / ((33 - 2 * nf) * math.log((mu / Lambda) ** 2))
    return min(als, 0.4)


def _r_ratio_parametric(m_dp: float) -> float:
    """Perturbative R-ratio for m >> 1 GeV."""
    nf = 3 + (m_dp > 2 * 1.27) + (m_dp > 2 * 4.18)
    if m_dp <= 0.2:
        return 0.0
    als = _alpha_s_1loop(m_dp, nf)

    eq2_sum = (2.0 / 3.0) ** 2  # u
    eq2_sum += (1.0 / 3.0) ** 2  # d
    eq2_sum += (1.0 / 3.0) ** 2  # s
    if m_dp > 2 * 1.27:
        eq2_sum += (2.0 / 3.0) ** 2  # c
    if m_dp > 2 * 4.18:
        eq2_sum += (1.0 / 3.0) ** 2  # b

    return 3.0 * eq2_sum * (1.0 + als / math.pi)


def perturbative_brs(mass: float) -> dict:
    """
    Compute dark-photon BRs via perturbative R-ratio for mass >> 1 GeV.

    Includes alpha_s(mass) NLO correction.  Valid for m > ~2 GeV.
    Returns dict mapping channel label -> BR (normalised to 1).

    Canonical implementation — imported by make_dp_cmnd.py and validate_brs.py.

    The denominator counts only kinematically open lepton channels:
      total = n_open_leptons + R
    This ensures BRs sum to 1 at all masses, including 1.7 < m < 3.55 GeV
    where tau is kinematically forbidden.
    """
    nf = 3 + (mass > 2 * 1.27) + (mass > 2 * 4.18)
    als = _alpha_s_1loop(mass, nf)

    # R = Nc * sum(e_q^2) for active quarks
    active_quarks = {"uu": 0, "dd": 0}
    active_quarks["uu"] += 1             # u quark
    active_quarks["dd"] += 2             # d, s quarks
    if mass > 2 * 1.27:
        active_quarks["uu"] += 1         # c quark
    if mass > 2 * 4.18:
        active_quarks["dd"] += 1         # b quark

    R = 3 * ((active_quarks["uu"] * (2/3)**2 + active_quarks["dd"] * (1/3)**2)
             * (1 + als / math.pi))

    # Count only kinematically open lepton flavors in the denominator
    lepton_thresholds = [("ee", 0.0), ("mumu", 2 * M_MU), ("tau", 2 * 1.77686)]
    n_open_lep = sum(1 for _, thresh in lepton_thresholds if mass > thresh)
    total = n_open_lep + R

    brs: dict = {}
    for lep_col, thresh in lepton_thresholds:
        brs[lep_col] = (1.0 / total) if mass > thresh else 0.0

    qcd_corr = 1 + als / math.pi
    for flavor, q2 in [
        ("uu_q", (2/3)**2),
        ("dd_q", (1/3)**2),
        ("ss_q", (1/3)**2),
        ("cc_q", (2/3)**2),
        ("bb_q", (1/3)**2),
    ]:
        thresh_map = {
            "uu_q": 0.0, "dd_q": 0.0, "ss_q": 0.0,
            "cc_q": 2 * 1.27, "bb_q": 2 * 4.18,
        }
        if mass > thresh_map[flavor]:
            brs[flavor] = 3 * q2 * qcd_corr / total
        else:
            brs[flavor] = 0.0

    return brs


def dp_width_eps1(m_dp: float) -> float:
    """
    Total A' decay width at ε=1, in GeV.

    Γ₀(m_A') = Γ(A'→e⁺e⁻) + Γ(A'→μ⁺μ⁻) + Γ(A'→τ⁺τ⁻) + Γ(A'→hadrons)
    where Γ(A'→hadrons) = Γ(A'→μ⁺μ⁻) × R(m_A').

    All partial widths at ε=1; multiply by ε² for physical width.
    """
    g_ee = _partial_width_ll(m_dp, M_EL)
    g_mumu = _partial_width_ll(m_dp, M_MU)
    g_tautau = _partial_width_ll(m_dp, 1.77686)
    R = _r_ratio_from_table(m_dp)
    g_had = g_mumu * R
    return g_ee + g_mumu + g_tautau + g_had


def dp_ctau_m(eps2: float, m_dp: float) -> float:
    """
    Dark photon cτ in meters.

    cτ = ℏc / (ε² × Γ₀(m_A'))

    Parameters
    ----------
    eps2 : float
        Kinetic mixing squared (ε²).
    m_dp : float
        Dark photon mass in GeV.

    Returns
    -------
    float
        cτ in meters. Returns inf if width is zero.
    """
    g0 = dp_width_eps1(m_dp)
    if g0 <= 0 or eps2 <= 0:
        return float('inf')
    gamma_phys = eps2 * g0  # physical width in GeV
    return HBAR_C_GEV_M / gamma_phys


def dp_ctau_mm(eps2: float, m_dp: float) -> float:
    """Dark photon cτ in mm (for Pythia8 tau0 parameter)."""
    return dp_ctau_m(eps2, m_dp) * 1e3


# ---------------------------------------------------------------------------
# Convenience: total BR(meson→A'X) summed over contributing mesons
# ---------------------------------------------------------------------------

def total_meson_br(eps2: float, m_dp: float) -> dict:
    """
    Compute BR(meson→A'+X) for all kinematically open meson channels.

    Returns dict {channel_name: br_value}.
    """
    channels = {}
    br_eta = br_eta_to_dp_gamma(eps2, m_dp)
    if br_eta > 0:
        channels["eta_to_dp_gamma"] = br_eta
    br_omega = br_omega_to_dp_pi0(eps2, m_dp)
    if br_omega > 0:
        channels["omega_to_dp_pi0"] = br_omega
    return channels


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Dark photon meson production BRs and lifetime")
    print("=" * 60)

    for m in [0.3, 0.5, 0.6, 1.0]:
        print(f"\nm_A' = {m} GeV:")
        eps2 = 1e-6
        print(f"  ε² = {eps2:.0e}")

        br_eta = br_eta_to_dp_gamma(eps2, m)
        br_omega = br_omega_to_dp_pi0(eps2, m)
        print(f"  BR(η→A'γ)   = {br_eta:.4e}")
        print(f"  BR(ω→A'π⁰)  = {br_omega:.4e}")

        g0 = dp_width_eps1(m)
        ctau = dp_ctau_m(eps2, m)
        print(f"  Γ₀(A') [ε=1] = {g0:.4e} GeV")
        print(f"  cτ(ε²={eps2:.0e}) = {ctau:.2e} m")

        R = _r_ratio_from_table(m)
        print(f"  R(m_A')      = {R:.3f}")
