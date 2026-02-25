# Code Check — Current Files

**Date:** 2026-02-25 (updated: adaptive 2D hitmap + ROOT audit + diff review)
**Scope:** All modified files in the repository.
**Method:** Manual review of current file state + full git diff audit.

## 1. Introduced by us

### 1.1 Physics

- **P4** (`pythiaStuff/alp_meson.cmnd:27`) — `PhaseSpace:pTHatMin = 20` introduces a kinematic bias without analysis-side event reweighting. *Note: biases toward high-pT B mesons. Acceptable if detector acceptance is dominated by high-pT ALPs, but worth being aware of.*

### 1.2 Portability

- **E5** (`pythiaStuff/make.sh:5`) — Hardcoded local `PYTHIA=` path reduces portability. *Note: acceptable for a single-user research repo.*

## 2. Not introduced by us

### 2.1 Bugs / correctness

*(none remaining)*

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
