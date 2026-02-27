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

- **W3** (`dark_photon/generator/heavy_dp_m{05,1,15}.cmnd`, `decayProbPerEvent_2body.py`) — Current signal definition is dimuon-only by construction.

  The generator forces `A' -> mu+ mu-` with `BR=1` in all three heavy mass-point cards. The analysis also assumes 2-body muons (`M_DAUGHTER = 0.10566` GeV), so current acceptance/exclusion outputs are a muon-channel efficiency map.

  For comparison to displaced-track searches (MATHUSLA/ANUBIS/CODEX-b) that accept generic visible charged decays, reinterpret yields as:
  `N_all-visible ~= N_mu-map x BR(A'->visible) / BR(A'->mu+mu-)`.

  In minimal kinetic-mixing scenarios without invisible dark-sector channels, `BR(A'->visible) ~= 1`, so the leading rescaling is approximately `1 / BR(A'->mu+mu-)` and is mass-dependent (hadronic/leptonic split set by the R-ratio).

  Caveat: this is first-order only. Hadronic decays change daughter multiplicity/topology relative to the current 2-body muon acceptance model.

  **Action:** for precision comparison, compute `BR(m_A')` from dedicated dark-photon BR tables/tools and/or generalize acceptance to `>=2` charged tracks.

### Minor

- **M1** — Light A' production (B→K(*)A') not yet implemented. Sub-GeV mass points (m=0.5, 1 GeV) are physically better motivated by B-meson FCNC production than by exotic Higgs decay, but h→A'A' is used as placeholder until the B-physics cmnd file is ready.
- **M2** — Analysis scripts load external curves relative to CWD. Run from `dark_photon/` for `external/*.csv` to resolve.

## Fixed during smoke test (2026-02-26)

- ~~**B1**~~ (`produce.sh`) — Aggregated `_meta.json` was missing `llp_pdg_id`. Now captured from batch sidecars and written to final meta.
- ~~**B2**~~ (`heavy_dp.cmnd`) — `isResonance = on` was wrong for LLP treatment. Pythia8 `isResonance = on` means decay as part of the hard process (resonance width), not as a displaced particle. Fixed to `off` so the LLP decays after full event generation.
- ~~**B3**~~ (`gargoyle_geometry.py`) — `rtree` crash at high particle counts (~200k) due to degenerate ray directions (NaN/inf from extreme-eta particles). Added `np.isfinite` guard and try/except around `mesh.ray.intersects_location`.

## Verified correct

- `heavy_dp_m{05,1,15}.cmnd`: PDG 6000115, spinType=3 (vector), `isResonance=off`, `LLP:pdgId` override, clean `addChannel` syntax (no index hack), `oneChannel` decay to μ⁺μ⁻. Three mass benchmarks: 0.5, 1, 15 GeV.
- `build.sh`: self-contained, compiles local `generator.cc`, requires `$PYTHIA8_DIR` (Pythia8 not a sibling of repo).
- `produce.sh`, `clean.sh`: paths self-resolve to `dark_photon/output/`, `llp_pdg_id` now propagated.
- `generator.cc`: reads `LLP:pdgId` from cmnd (overrides default 6000113), writes correct `llp_pdg_id` to per-batch `_meta.json`.
- Analysis scripts (`decayProbPerEvent_2body.py`, `signal_surface_hitmap_v2.py`, `gargoyle_geometry.py`, `visualize_tunnel.py`): self-contained copies in `dark_photon/`; run from `dark_photon/` for external curves to resolve.
