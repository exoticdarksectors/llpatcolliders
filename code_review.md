# Code Check — Current Files

**Date:** 2026-02-24 (updated after parallel-production + output/ refactor)
**Scope:** All modified files in the repository.
**Method:** Manual review of current file state.

## 1. Introduced by us

### 1.1 Bugs / correctness

- **C3** (`pythiaStuff/main144.cc:296,331`) — CSV header uses plain commas, data rows use `,\t` separators. Pandas handles this fine (treats tab as whitespace), but the `awk` in `parallel_produce.sh:56-57` depends on `,\t` being the exact separator. If main144 changes format, the renumbering breaks.

- **N1** (`pythiaStuff/parallel_produce.sh:11`) — Example comment still shows `alp_meson.cmnd 10000 alp_light_m1 8 50000` with explicit n_jobs and large batch. With BR=1.0 both ALPs fill fast; the example is stale.

- **N2** (`pythiaStuff/alp_meson.cmnd:86`) — Comment says "ALP yield per B is ~1e-5" but BR is now 1.0.

### 1.2 Physics

- **P4** (`pythiaStuff/alp_meson.cmnd:27`) — `PhaseSpace:pTHatMin = 20` introduces a kinematic bias without analysis-side event reweighting. *Note: this biases toward high-pT B mesons. Acceptable if the detector acceptance is dominated by high-pT ALPs, but worth being aware of.*

- **P6** (`pythiaStuff/higgsLL.cmnd:38`, `pythiaStuff/alp_meson.cmnd:41`, `decayProbPerEvent_2body.py:12`) — Generator decays ALP → μ⁺μ⁻ while some comments/variable names say "electron". `M_DAUGHTER` is correctly set to muon mass; the issue is cosmetic (comments at lines 26, 595 say "electron momentum").

### 1.3 Portability

- **E5** (`pythiaStuff/make.sh:5`) — Hardcoded local `PYTHIA=` path reduces portability. *Note: acceptable for a single-user research repo.*

## 2. Not introduced by us

### 2.1 Bugs / correctness

- **C2** (`pythiaStuff/main144.cc:357`) — `delete file, tree, evt;` only deletes `file` (comma-operator). *Note: ROOT-only path (`#ifdef PY8ROOT`), doesn't affect CSV workflow.*

- ~~**E1**~~ — **Fixed.** `signal_surface_hitmap_v2.py:225` had `[x, y, Z]` but the correct CMS convention is `[x, Z, y]` (survey x→CMS x, Z_POSITION=22m→CMS y=vertical, survey y→CMS z=beam). Now matches `decayProbPerEvent_2body.py:551`.

- **E2** (`decayProbPerEvent_2body.py:419` vs `469-470`) — Docstring says sample `|cosθ*|` in `[0, β]`, code samples `[0, 1]`. Code at line 470 is `cos_theta_star = rng.uniform(0, 1, ...)` which is correct for isotropic decay; the docstring at line 419 is wrong.

### 2.2 Warnings

- **W1** (`decayProbPerEvent_2body.py:403-404`) — Divide-by-zero if `mean_p1 == 0` in exclusion formula. *Note: would produce `inf` in the plot, not a crash. Numpy handles this.*

- **W2** (`decayProbPerEvent_2body.py:28`) — Comment says "10 cm" but `SEP_MAX = 1.0` m. Comment is wrong.

- **W4** (`decayProbPerEvent_2body.py:291`) — Per-row `iterrows()` ray-casting is slow. *Note: not a bug, just slow. Multiprocessing was rejected as too complex; this is fine for ~10k particles.*

- **W5** (`decayProbPerEvent_2body.py:369`) — CSV reloaded from disk in every call to `process_with_acceptance` during lifetime scan. *Note: lifetime scan calls this 20 times. Could pass DataFrame instead of path. Minor perf issue.*

- **W6** (`decayProbPerEvent_2body.py:770-774`) — External CSV loads (`MATHUSLA.csv` etc.) have no error handling. Will crash if `external/` folder is missing.

- **W7** (`pythiaStuff/higgsLL.cmnd:22-23`) — Hardcoded channel index `25:76:*` is Pythia8-version-fragile. If Higgs decay table changes between versions, this could silently point at the wrong channel.

- **W8** (`pythiaStuff/main144.cc:251,291`) — Debug `std::cout << "here"` remains.

- **W9** (`pythiaStuff/main144.cc:252-254`) — ROOT pointers (`file`, `tree`, `evt`) declared but uninitialized when `writeRoot` is false. *Note: only used inside `if (writeRoot)` blocks, so safe in practice.*

- **W10** (`pythiaStuff/make.sh:7`) — `-w` suppresses all compiler warnings.

### 2.3 Minor / style

- **M1** (`pythiaStuff/main144.cc:198`) — Double semicolon.
- ~~**M2**~~ — Removed: `outString` no longer exists; filenames now derive from CSV basename.
- **M3** (`decayProbPerEvent_2body.py:162-163`) — `epsabs=1e-15` is overly strict for this physics application.
- **M6** (`pythiaStuff/higgsLL.cmnd:17`) — `mMin = m0 = 125` narrows the sampled Higgs lineshape to a delta function.
- **M7** (`pythiaStuff/make.sh`) — No `set -euo pipefail`; compilation errors may be silently ignored.

### 2.4 Removed (fixed or no longer applicable)

- ~~**C4** (coordinate mismatch)~~ — Verified both scripts use same `(x, Z_POSITION, y)` convention.
- ~~**P3** (backgroundFullGeo.py muon-separation)~~ — Optional script, not in active workflow.
- ~~**W3** (backgroundFullGeo.py comment)~~ — Optional script, not in active workflow.
- ~~**P5** (SEP_MAX mismatch)~~ — backgroundFullGeo.py is a separate analysis with intentionally different cuts.
- ~~**M4** (mutable default argument)~~ — backgroundFullGeo.py, not in active workflow.
