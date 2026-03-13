"""
production/constants.py

Physical constants for HNL production at LHC 14 TeV.

Masses via the `particle` package (PDG 2024).
Cross-sections from FONLL NLO+NLL.
Fragmentation fractions from PDG/HFLAV/ALICE.

References:
  - FONLL: Cacciari, Greco, Nason (NLO+NLL heavy-quark production)
  - Fragmentation fractions: PDG 2024, ALICE D-meson measurements, HFLAV
  - Bc: LHCb σ(Bc)/σ(B+) ≈ 2.6e-3 (CMS/LHCb)
  - K-factor: NLO/LO correction for W/Z → ℓN production
"""

from particle import Particle

# ==========================================================================
# Particle masses (GeV) from PDG via `particle` package
# ==========================================================================

M_B0 = Particle.from_pdgid(511).mass * 1e-3      # B0
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3    # B+
M_BS = Particle.from_pdgid(531).mass * 1e-3       # Bs
M_BC = Particle.from_pdgid(541).mass * 1e-3       # Bc+
M_D0 = Particle.from_pdgid(421).mass * 1e-3       # D0
M_DPLUS = Particle.from_pdgid(411).mass * 1e-3    # D+
M_DS = Particle.from_pdgid(431).mass * 1e-3       # Ds+
M_TAU = Particle.from_pdgid(15).mass * 1e-3       # tau
M_ELECTRON = Particle.from_pdgid(11).mass * 1e-3  # electron
M_MUON = Particle.from_pdgid(13).mass * 1e-3      # muon
M_PION = Particle.from_pdgid(211).mass * 1e-3     # pi+
M_KAON = Particle.from_pdgid(321).mass * 1e-3     # K+

# ==========================================================================
# FONLL inclusive cross-sections at 14 TeV (pb)
# ==========================================================================

# σ(pp → bb̄) ~ 495 μb at 14 TeV (FONLL NLO+NLL, LHCb-validated)
SIGMA_BBBAR_PB = 495.0e6

# σ(pp → cc̄) ~ 23.6 mb at 14 TeV (FONLL NLO+NLL)
SIGMA_CCBAR_PB = 23.6e9

# σ(pp → Bc) ~ 0.9 μb at 14 TeV (BCVEGPY/FONLL, CMS/LHCb)
SIGMA_BC_PB = 0.9e6

# ==========================================================================
# Fragmentation fractions
# ==========================================================================

# Beauty fragmentation (PDG/HFLAV)
# P(b → B-hadron species), particle+antiparticle combined
FRAG_B = {
    521: 0.408,   # B+/B-
    511: 0.408,   # B0/B0bar
    531: 0.100,   # Bs
}
# Bc fraction: σ(Bc)/σ(B+) ≈ 2.6e-3 from LHCb
F_BC = 2.6e-3

# Charm fragmentation (ALICE/PDG, relative fractions)
# These are fractions of c-quarks hadronizing into each species.
# The FONLL meson table already includes fragmentation → these fractions
# tell us how to split the FONLL charm meson cross-section among species.
FRAG_C = {
    421: 0.542,   # D0/D0bar (dominant)
    411: 0.225,   # D+/D-
    431: 0.080,   # Ds+/Ds-
}

# ==========================================================================
# EW K-factor
# ==========================================================================

K_FACTOR_EW = 1.3  # NLO/LO correction for W/Z → ℓN production

# ==========================================================================
# Meson species lookup tables
# ==========================================================================

# PDG ID → mass (GeV)
MESON_MASSES = {
    521: M_BPLUS, 511: M_B0, 531: M_BS, 541: M_BC,
    421: M_D0, 411: M_DPLUS, 431: M_DS,
}

# Lepton PDG ID → mass (GeV), keyed by flavor label
LEPTON_MASSES = {
    'Ue': M_ELECTRON,
    'Umu': M_MUON,
    'Utau': M_TAU,
}

# Flavor label → lepton PDG ID
FLAVOR_TO_LEPTON_PDG = {
    'Ue': 11,
    'Umu': 13,
    'Utau': 15,
}

# Flavor label → MadGraph flavor name
FLAVOR_TO_MG5 = {
    'Ue': 'electron',
    'Umu': 'muon',
    'Utau': 'tau',
}

# Quark type → list of (pdg_id, frag_fraction) for meson species
QUARK_MESON_MAP = {
    'bottom': [(521, FRAG_B[521]), (511, FRAG_B[511]), (531, FRAG_B[531])],
    'charm':  [(421, FRAG_C[421]), (411, FRAG_C[411]), (431, FRAG_C[431])],
}
