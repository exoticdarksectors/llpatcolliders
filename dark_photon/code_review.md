# Code Check — dark_photon/

**Date:** 2026-02-26
**Scope:** All files under `dark_photon/`.
**Method:** Manual review of all dark_photon/ files.

## Open issues

### Warnings

- **W1** (`dark_photon/external/`) — No dedicated dark-photon comparison curves yet.

  **Status of external curves for h→A'A':**
  Published dark-Higgs h→SS exclusion contours (MATHUSLA, CODEX-b, ANUBIS) at m_S ~ 1 and
  15 GeV (BR(h→SS) vs cτ) are a valid geometric-acceptance comparison for h→A'A':
  same production mechanism (exotic Higgs decay), same mass range, same boost distribution.
  Scalar-vs-vector spin is a second-order effect on geometric acceptance.
  **Action:** Copy the relevant h→SS CSV files into `dark_photon/external/` for overlay.

  **Dedicated h→A'A' curves from the literature:**
  - CODEX-b: Figure 7 of [arXiv:1911.00481](https://arxiv.org/abs/1911.00481) and
    Figure 3 (top) of [arXiv:2203.07316](https://arxiv.org/abs/2203.07316) show
    h→A'A'→4 charged tracks, BR(h→A'A') vs cτ, for m_A' = 0.4 and 1.0 GeV.
    No digitised data available; would need manual digitisation.
  - MATHUSLA: [arXiv:2504.01999](https://arxiv.org/abs/2504.01999) (ESPP 2025) shows
    only hadronic LLP from h→XX, no explicit A' curves. The earlier physics case
    [arXiv:1806.07396](https://arxiv.org/abs/1806.07396) may contain A' projections
    in later sections but was not fully accessible.
  - ANUBIS: [arXiv:2510.26932](https://arxiv.org/abs/2510.26932) and
    [arXiv:1909.13022](https://arxiv.org/abs/1909.13022) show only h→ss (dark scalar);
    the HAHM dark photon Z_d is explicitly decoupled (mass set to ∞).

  **Conclusion:** The dark-Higgs h→SS curves provide a meaningful comparison and will be
  overlaid automatically once copied into `dark_photon/external/`.
  Dedicated h→A'A' digitisation is low priority unless the paper requires exact A' contours.
  The only published h→A'A' curves found are from CODEX-b at m_A' = 0.4, 1.0 GeV
  (too light for the default m=15 GeV benchmark).

- **W2** (`dark_photon/generator/build.sh`) — `-w` suppresses all compiler warnings.

- **W4** — R-ratio BRs in cmnd files are approximate (hand-computed from perturbative N_c×e_q²). At low masses (0.5–1 GeV), hadronic BRs are distorted by ρ/ω/φ resonances; the perturbative approximation is crude.

  **However, this is not an issue for the current `heavy_dp.cmnd` setup.** The h→A'A' production mode targets multi-GeV masses; the sub-GeV benchmarks (m=0.5, 1 GeV) with Higgs-portal production are placeholders (see M1). The ρ/ω/φ resonance region only matters once light DP production (B→K(\*)A') is implemented in `light_dp.cmnd`.

  **When `light_dp.cmnd` is added:** use tabulated BRs from [DeLiVeR](https://github.com/preimitz/deliver) ([arXiv:2201.01788](https://arxiv.org/abs/2201.01788)) or [DarkCast](https://gitlab.com/darkcast/releases), both of which implement Vector Meson Dominance to handle the resonance region correctly. The HAHM [BrTableData.txt](https://github.com/davidrcurtin/HAHM) (PDG experimental R-ratio data below 12 GeV) is the simplest drop-in table.

### Minor

- **M1** — Light A' production (B→K(*)A') not yet implemented. Sub-GeV mass points (m=0.5, 1 GeV) are physically better motivated by B-meson FCNC production than by exotic Higgs decay, but h→A'A' is used as placeholder until the B-physics cmnd file is ready.

## Verified correct

- `heavy_dp_m{05,1,15}.cmnd`: PDG 6000115, spinType=3 (vector), `isResonance=off`, `LLP:pdgId` override, clean `addChannel` syntax, R-ratio decay BRs. Three mass benchmarks: 0.5, 1, 15 GeV.
- `build.sh`: self-contained, compiles local `generator.cc`, requires `$PYTHIA8_DIR`.
- `produce.sh`, `clean.sh`: paths self-resolve to `dark_photon/output/` with `output/data` (CSV/meta) and `output/images` (plots); `llp_pdg_id` propagated, daughter CSV aggregated.
- `generator.cc`: reads `LLP:pdgId` from cmnd, writes LLP CSV + `_daughters.csv` (charged final-state particles from each LLP decay) + `_meta.json`.
- `decayProbPerEvent_Ntrack.py`: N-track displaced vertex analysis (>=2 charged daughters, p > 600 MeV, union of all valid daughter-pair separation windows [1 mm, 1 m]; opening angles via `arctan2` for numerical stability at µrad–mrad scale). Reads LLP CSV + daughters CSV. External curves resolved relative to script file.
- `decayProbPerEvent_2body.py`: legacy 2-body muon analysis (for cross-checks only).
- `signal_surface_hitmap_v2.py`, `gargoyle_geometry.py`, `visualize_tunnel.py`: self-contained copies; run from `dark_photon/`.
