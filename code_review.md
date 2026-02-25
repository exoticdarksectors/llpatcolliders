# Code Check — Current Files

**Date:** 2026-02-25 (updated after gargoyle_geometry refactor + world_up fix)
**Scope:** All modified files in the repository.
**Method:** Manual review of current file state.

## 1. Introduced by us

### 1.1 Physics

- **P4** (`pythiaStuff/alp_meson.cmnd:27`) — `PhaseSpace:pTHatMin = 20` introduces a kinematic bias without analysis-side event reweighting. *Note: biases toward high-pT B mesons. Acceptable if detector acceptance is dominated by high-pT ALPs, but worth being aware of.*

### 1.2 Portability

- **E5** (`pythiaStuff/make.sh:5`) — Hardcoded local `PYTHIA=` path reduces portability. *Note: acceptable for a single-user research repo.*

## 2. Not introduced by us

### 2.1 Bugs / correctness

- **C2** (`pythiaStuff/main144.cc:356`) — `delete file, tree, evt;` only deletes `file` (comma-operator). *Note: ROOT-only path (`#ifdef PY8ROOT`), doesn't affect CSV workflow.*

### 2.2 Warnings

- **W1** (`decayProbPerEvent_2body.py:403-404`) — Divide-by-zero if `mean_p1 == 0` in exclusion formula. *Note: produces `inf` in the plot, not a crash. Numpy handles this.*

- **W5** (`decayProbPerEvent_2body.py:369`) — CSV reloaded from disk in every call to `process_with_acceptance` during lifetime scan (20 times). Could pass DataFrame instead of path. *Minor perf issue.*

- **W6** (`decayProbPerEvent_2body.py:770-774`) — External CSV loads (`MATHUSLA.csv` etc.) have no error handling. Will crash if `external/` folder is missing.

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
