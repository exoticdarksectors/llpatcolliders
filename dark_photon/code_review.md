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
  The analysis code (`decayProbPerEvent_2body.py` line 646 `else` branch) already overlays
  these curves for any PDG ID ≠ 9000001, so PDG 6000115 will get them automatically.
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

- ~~**W3**~~ FIXED: Signal definition generalized to ≥2 charged tracks (displaced vertex), matching MATHUSLA/ANUBIS/CODEX-b.
  - cmnd files now use R-ratio BRs (e⁺e⁻, μ⁺μ⁻, τ⁺τ⁻, qq̄ with Pythia8 hadronization) instead of forcing A'→μ⁺μ⁻.
  - `generator.cc` writes `_daughters.csv` with all final-state charged particles from each LLP decay.
  - `produce.sh` aggregates daughter CSVs across batches.
  - New analysis: `decayProbPerEvent_Ntrack.py` (geometric acceptance × decay probability × ≥N tracks with p > p_cut).
  - Old `decayProbPerEvent_2body.py` kept for backward compatibility/cross-checks.

- **W4** — R-ratio BRs in cmnd files are approximate (hand-computed from perturbative N_c×e_q²). At low masses (0.5–1 GeV), hadronic BRs are distorted by ρ/ω/φ resonances; the perturbative approximation is crude. For precision, use tabulated e⁺e⁻→hadrons R-ratio data or a dark photon BR calculator (e.g. from [arXiv:2005.01515](https://arxiv.org/abs/2005.01515)).

### Minor

- **M1** — Light A' production (B→K(*)A') not yet implemented. Sub-GeV mass points (m=0.5, 1 GeV) are physically better motivated by B-meson FCNC production than by exotic Higgs decay, but h→A'A' is used as placeholder until the B-physics cmnd file is ready.
- **M2** — Analysis scripts load external curves relative to CWD. Run from `dark_photon/` for `external/*.csv` to resolve.

## Fixed during smoke test (2026-02-26)

- ~~**B1**~~ (`produce.sh`) — Aggregated `_meta.json` was missing `llp_pdg_id`. Now captured from batch sidecars and written to final meta.
- ~~**B2**~~ (`heavy_dp.cmnd`) — `isResonance = on` was wrong for LLP treatment. Pythia8 `isResonance = on` means decay as part of the hard process (resonance width), not as a displaced particle. Fixed to `off` so the LLP decays after full event generation.
- ~~**B3**~~ (`gargoyle_geometry.py`) — `rtree` crash at high particle counts (~200k) due to degenerate ray directions (NaN/inf from extreme-eta particles). Added `np.isfinite` guard and try/except around `mesh.ray.intersects_location`.

## Verified correct

- `heavy_dp_m{05,1,15}.cmnd`: PDG 6000115, spinType=3 (vector), `isResonance=off`, `LLP:pdgId` override, clean `addChannel` syntax, R-ratio decay BRs. Three mass benchmarks: 0.5, 1, 15 GeV.
- `build.sh`: self-contained, compiles local `generator.cc`, requires `$PYTHIA8_DIR`.
- `produce.sh`, `clean.sh`: paths self-resolve to `dark_photon/output/`, `llp_pdg_id` propagated, daughter CSV aggregated.
- `generator.cc`: reads `LLP:pdgId` from cmnd, writes LLP CSV + `_daughters.csv` (charged final-state particles from each LLP decay) + `_meta.json`.
- `decayProbPerEvent_Ntrack.py`: N-track displaced vertex analysis (≥2 charged tracks, p > 600 MeV). Reads LLP CSV + daughters CSV.
- `decayProbPerEvent_2body.py`: legacy 2-body muon analysis (for cross-checks only).
- `signal_surface_hitmap_v2.py`, `gargoyle_geometry.py`, `visualize_tunnel.py`: self-contained copies; run from `dark_photon/`.
