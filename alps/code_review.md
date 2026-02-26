# Code Check — Current Files

**Date:** 2026-02-25 (updated: adaptive 2D hitmap + ROOT audit + diff review)
**Scope:** All modified files in the repository.
**Method:** Manual review of current file state + full git diff audit.

## 1. Introduced by us

### 1.1 Physics

- **P5** ★ CRITICAL (`external/*.csv` + `decayProbPerEvent_2body.py`) — **External comparison curves are for dark Higgs only; unusable for alps.**
  The CSV files in `external/` (MATHUSLA.csv, CODEX.csv, ANUBIS.csv, ANUBISOpt.csv, ANUBISUpdateCons.csv) were digitized from **dark Higgs** (H(125)→SS) exclusion contours at m_S = 1 and 15 GeV. They plot BR(h→SS) vs cτ. All 5 curves have peak sensitivity at cτ ~ 4–30 m and are overlaid unconditionally on every exclusion plot regardless of the sample being analyzed.
  - **Heavy ALP (h→aa, m=15 GeV):** The curves are a **valid comparison** — same production mechanism (h→XX), same mass, similar boost distribution. Scalar-vs-pseudoscalar is a second-order effect for geometric acceptance. However, confirm with Matt which specific curves are at 15 GeV vs 1 GeV, and only overlay the mass-matched subset.
  - **Light ALP (B→Ka, m=1 GeV):** The curves are **invalid** even at 1 GeV. Production is h→SS (Higgs portal, isotropic in Higgs rest frame) vs B→Ka (FCNC, forward-boosted B mesons). Completely different boost/angular distributions and production cross-sections. The overlay is an apples-to-oranges comparison and **must not appear in the paper**.
  - For the alps branch, new external comparison curves must be digitized from published **light-scalar-via-B-decay** sensitivity studies:
    - MATHUSLA: [arXiv:1806.07396](https://arxiv.org/abs/1806.07396) §9.1 (SM+S scalar via B decay); data files at [MATHUSLA_LLPfiles_SMS](https://github.com/davidrcurtin/MATHUSLA_LLPfiles_SMS)
    - CODEX-b: [arXiv:1911.00481](https://arxiv.org/abs/1911.00481), [arXiv:2203.07316](https://arxiv.org/abs/2203.07316)
    - ANUBIS: [arXiv:1909.13022](https://arxiv.org/abs/1909.13022), [arXiv:2401.11604](https://arxiv.org/abs/2401.11604)
  - Note: published light-scalar curves use **sin²θ vs mS** axes, not BR vs cτ. Either convert our plot to that parameter space or map their contours to BR(B→KS) vs cτ using BR ∝ sin²θ and cτ ∝ 1/sin²θ.
  - **Action:** (1) Ask Matt which CSVs are at 1 GeV and which at 15 GeV. (2) For the higgs branch, only overlay the mass-matched curves. (3) For the alps branch, digitize B→KS exclusion contours from the above references, save as new CSVs (e.g. `external/MATHUSLA_BKS.csv`), and load those instead. (4) Do not use the existing `external/` CSVs for alps plots.

### 1.2 Portability

*(none remaining)*

## 2. Not introduced by us

### 2.1 Bugs / correctness

*(none remaining — see §3 for fixes applied this session)*

### 2.2 Warnings

- **W7** (`pythiaStuff/higgsLL.cmnd:22-23`) — Hardcoded channel index `25:76:*` is Pythia8-version-fragile. If Higgs decay table changes between versions, this could silently point at the wrong channel.

- **W10** (`pythiaStuff/make.sh:7`) — `-w` suppresses all compiler warnings.

### 2.3 Minor / style

- **M3** (`decayProbPerEvent_2body.py:162-163`) — `epsabs=1e-15` is overly strict for this physics application.
- **M6** (`pythiaStuff/higgsLL.cmnd:17`) — `mMin = m0 = 125` narrows the sampled Higgs lineshape to a delta function. *Intentional for fixed-mass setup.*
- **M7** (`pythiaStuff/make.sh`) — No `set -euo pipefail`; compilation errors may be silently ignored.

## 3. Fixed this session (2026-02-25)

- ~~**C3**~~ — CSV header now uses `,\t` separators matching data rows (main144.cc + parallel_produce.sh).
- ~~**E1**~~ — Axis mismatch fixed in signal_surface_hitmap_v2.py: `[x, y, Z]` → `[x, Z, y]`.
- ~~**E2**~~ — Docstring corrected: `[0, β]` → `[0, 1]`.
- ~~**N1**~~ — Stale example in parallel_produce.sh updated.
- ~~**N2**~~ — Stale "yield per B is ~1e-5" comment in alp_meson.cmnd updated.
- ~~**P6**~~ — "electron" → "muon"/"daughter" in comments.
- ~~**W2**~~ — Wrong comment "10 cm" removed (SEP_MAX = 1.0 m).
- ~~**W8**~~ — Debug `cout << "here"` removed from main144.cc.
- ~~**M1**~~ — Double semicolon removed from main144.cc.
- ~~**M2**~~ — `outString` removed; filenames derive from CSV basename.
- ~~**E3**~~ — `world_up` vector corrected from `[0,0,1]` (Z-up) to `[0,1,0]` (CMS Y-up) in `create_profile_mesh`. Was rotating the tunnel cross-section 90°. Fix is structural: all three alps scripts now import from `gargoyle_geometry.py` which uses the correct convention.
- ~~**E4**~~ — Axis-order bug in `backgroundFullGeo.py`: path was `[[x, y, Z_POSITION]]` (vertical at index 2) instead of `[[x, Y_POSITION, z]]` (vertical at index 1). Fixed by importing `path_3d_fiducial` from `gargoyle_geometry`.
- ~~**D1**~~ — Code duplication eliminated: geometry code (tunnel constants, profile functions, mesh builder, ray-caster, centerline vertices) was copy-pasted across three files. Replaced with a single shared module `gargoyle_geometry.py` (sourced from upstream higgs branch). `visualize_tunnel.py` also added for geometry validation.
- ~~**D2**~~ — `signal_surface_hitmap_v2.py` now aligned with `origin/main:higgs/signal_surface_hitmap.py`: added adaptive 2D region-growing (`adaptive_2d_efficiency`), dual efficiency curves (arc vs adaptive), density cross-section plot, and comparison summary. Alps-specific improvements kept (`--outdir`, `argparse`, `df.columns.str.strip()`).
- ~~**C3b**~~ — CSV separators changed from `,\t` to plain `,` everywhere (main144.cc header+data, parallel_produce.sh header+awk). Eliminates tab-in-column-name issue. All Python scripts now also strip column names on load (`df.columns.str.strip()`).
- ~~**N3**~~ — `decayProbPerEvent_2body.py` now warns when `--xsec` is not explicitly passed. README and howto.md updated to require explicit `--xsec` and `--lumi`.
- ~~**C2**~~ — `delete file, tree, evt;` → separate `delete` statements in main144.cc ROOT path.
- ~~**W1**~~ — Divide-by-zero guard: exclusion formula now returns `np.inf` when denominator is zero instead of relying on numpy's runtime handling.
- ~~**W5**~~ — `process_with_acceptance` now accepts a DataFrame directly; `analyze_decay_vs_lifetime` loads CSV once and passes the DataFrame for all ~20 lifetime steps.
- ~~**W6**~~ — External comparison curves (`MATHUSLA.csv` etc.) wrapped in try/except; missing files print a note and are skipped instead of crashing.
- ~~**W11**~~ — `((failed++))` → `((failed++)) || true` in parallel_produce.sh to prevent `set -e` abort on first job failure.
- ~~**W12**~~ — Added `MAX_ROUNDS=200` and consecutive-zero-progress guard (3 rounds) to parallel_produce.sh. Script now aborts with an error message instead of looping forever.
- ~~**W13**~~ — clean_production.sh now also cleans `*.root` fragments.
- ~~**E5**~~ — `make.sh` no longer hardcodes a user-specific `PYTHIA=` path. Now resolves via `$PYTHIA8_DIR` > sibling `../pythia8315` > `/usr/local` fallback. `howto.md` uses generic placeholders.
- ~~**P4**~~ — `pTHatMin = 20` removed from `alp_meson.cmnd`; generation is now inclusive (pTHatMin=0, standard practice). Existing samples must be regenerated.
- ~~**P4b**~~ — `--xsec` default and docs corrected: light ALP now uses σ(pp→bb̄, inclusive) = 3.73×10⁸ fb (from Pythia at 13.6 TeV), not the old placeholder 52000 fb which was ~7200× too small. Heavy ALP default unchanged at 60000 fb. Diagnostic table mass (`test_mass`) now reads from CSV instead of hardcoded 15 GeV.
- ~~**P7**~~ — **Light-ALP BR(B→Ka)≠1 due to `addChannel`**: `alp_meson.cmnd` used `addChannel` which *adds* ALP channels on top of SM decays, giving effective BR ≈ 0.72 instead of 1.0. Fixed by using `oneChannel` (wipes SM) + `addChannel` per B species. Now BR(B±/B⁰ → K(*)a) = 1.0 exactly. **Requires re-generation of light-ALP samples.**
- ~~**P8**~~ — **0-LLP event bias in normalization**: `decayProbPerEvent_2body.py` averaged `P(≥1 decay)` only over events present in CSV (i.e., events with ALPs). Events where both b-hadrons go to Bs/Λb produce no ALP and are missing from CSV (~18% of events). This inflated sensitivity by ~23%. Fixed: `main144.cc` now writes `_meta.json` sidecar with `n_generated`; `parallel_produce.sh` aggregates metadata; analysis reads it and averages over all generated events. Accepts `--n-events` CLI fallback.
- ~~**P9**~~ — **Inconsistent acceptance model**: `signal_surface_hitmap_v2.py` used `E_CUT` (energy cut) while `decayProbPerEvent_2body.py` used `P_CUT` (momentum cut with `E_min = sqrt(p²+m²)`). Hitmap also lacked `SEP_MAX`. Harmonized hitmap to use `P_CUT` with same formula. Hitmap intentionally omits `SEP_MAX` (not needed for relative spatial weighting).
- ~~**P10**~~ — **"e⁺e⁻" label on μ⁺μ⁻ decay**: surface hitmap title said "where e⁺e⁻ hit wall" but generator decay is a→μ⁺μ⁻. Fixed. Also fixed "electron" references in `decayProbPerEvent_2body.py` docstrings.

## 4. ROOT audit (2026-02-25)

**Current branch (alps):**
- Build is non-ROOT by default (`make.sh` has no `-DPY8ROOT`).
- `main144.cc` exits with error if `writeRoot=on` is set without ROOT build (line 215).
- CSV output is always written independently of ROOT (line 293, 325-332).
- ROOT output is optional extra I/O of the same filtered LLPs (requires `-DPY8ROOT` at compile).
- Verified: non-ROOT vs ROOT(writeRoot=off) vs ROOT(writeRoot=on) produce identical CSVs (hash-verified, same seed/events).

**origin/main:higgs behavior (historical):**
- Their build recipe linked ROOT and dictionary (`main144Dct.so`) with `-DPY8ROOT`.
- In that older `main144.cc`, CSV filling was inside `#ifdef PY8ROOT` and `if(writeRoot)`, so ROOT mode was effectively required for non-empty CSV.
- Old ROOT path had a filter bug (`if (!abs(prt.id()) == 6000113) continue;`) capturing far more particles than intended.

**For the paper:** use the current non-ROOT CSV path as baseline; treat ROOT output as auxiliary diagnostics only.
