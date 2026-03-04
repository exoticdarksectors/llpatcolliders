# Code Check — dark_photon/

**Date:** 2026-02-26 (last updated 2026-03-04)
**Scope:** All files under `dark_photon/`.

## Open issues

### Warnings

- **W2** (`dark_photon/generator/build.sh`) — `-w` suppresses all compiler warnings.

## Closed issues (summary)

- **W1** (2026-02-28): h→SS dark-Higgs external curves (MATHUSLA, CODEX, ANUBIS, ×5) copied to
  `dark_photon/external/`; valid geometric proxy for h→A'A' at same mass and boost distribution.
  No dedicated h→A'A' curves at m=15 GeV published; CODEX-b h→A'A' curves exist only at
  m_A'=0.4, 1.0 GeV (arXiv:1911.00481, arXiv:2203.07316).

- **W4** (2026-02-28): R-ratio BRs replaced by DeLiVeR VMD table (0.2–1.7 GeV) + perturbative
  QCD above 1.7 GeV. `make_dp_cmnd.py` generates all cmnd files. Validated: 47/49 mass points
  pass <10% (on-grid masses <1%). Artifacts: `output/data/dp_br_validation.csv`,
  `output/images/dp_br_vs_mass.png`.

- **W5** (2026-02-28): B-meson pT shape bias not relevant — light DP uses meson (η/ω) decay via
  `SoftQCD:nonDiffractive`, not B→K(*)A'. B→K(*)A' via EM penguin is valid but not competitive:
  rate ∝ ε² and cτ ∝ 1/ε², yielding <1 A' at HL-LHC for PX56-relevant lifetimes.

- **M1** (2026-03-01): Light A' at m=0.5 GeV via η→A'γ, ω→A'π⁰ (`SoftQCD:nonDiffractive`,
  BR=1 efficiency map). **Result: zero PX56 sensitivity** — meson-produced A' too forward/soft.
  Peak n_signal ~10⁻⁶, six orders of magnitude below 3-event threshold.

- **M2** (2026-03-03): Light A' at m=1.0 GeV via Drell-Yan (`NewGaugeBoson:ffbar2gmZZprime`,
  PDG 32). All meson portals kinematically closed at 1 GeV. Correction factor K≈8.09 applied to
  σ_Pythia. Note: `DM:ffbar2Zp2XX` (PDG 55) does NOT work (σ=0). **Result: zero PX56
  sensitivity**. Peak N = 5.5×10⁻⁸ at ε²=3.2×10⁻¹⁵, cτ=8 m.
