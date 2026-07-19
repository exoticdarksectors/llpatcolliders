# BC4 — light dark scalar mixing with the Higgs (GRENDEL)

Purely-additive package implementing the PBC benchmark **BC4**: a light scalar
`S` that mixes with the Higgs through one angle `theta` (`sin^2 theta`),
produced inclusively in `b -> X_s S` (with two-body `B -> K S` recoil as a
kinematic proxy) and decaying `S -> SM` via the mixing. This is the
low-mass partner of the already-done BC5 (`higgs/`, `h -> SS`).

Everything geometric and reconstruction-related is **imported, not copied**, from
the shared single sources:

* `higgs/grendel_geometry.py` — GRENDEL PX56 mesh, ray-casting, `points_on_tracker`.
* `higgs/reco_common.py` — bounded 4-hit two-track vertex reco + timing chi^2.
* `hnl/analysis/decay_reco_acceptance.py` — `reconstruct_decays`, `selection_mask`
  (PR #13 selection), and `scan_u2` (the coupling reweighting).
* `hnl/production/fonll/` + `hnl/production/decay_engine/kinematics.py` — the FONLL
  B-meson sampler and 2-body decay used for the parent kinematics.

## Why the limit is set differently from `higgs/` (BC5)

BC5 is set model-agnostically by scanning `(BR, c*tau)` independently. **BC4 is
coupling-controlled**: the single `sin^2 theta`, at fixed `m_S`, fixes the
production rate, the lifetime `c*tau`, *and* the visible BR simultaneously. So,
like the HNL, the limit is set **model-completely**: for each `(m_S, sin^2 theta)`
we evaluate the actual `c*tau`, production yield and visible BR, require
`N_signal >= 3` (background-free, 3000 fb^-1), and report a **closed island** in
`(m_S, sin^2 theta)` with a lower edge (too little production) and an upper edge
(decays before reaching the detector).

## Layout

| file | role |
|------|------|
| `model.py` | **model layer** (the core deliverable): `c*tau(m_S, theta)`, `BR(B->K S)` / `B->X_s S` / `K->pi S`, and the visible BRs (`mu mu, ee, tau tau, pi pi, K K, s s, c c, g g`) from Winkler arXiv:1809.01876. The BC4 analogue of HNLCalc. |
| `production.py` | FONLL bottom pool -> inclusive `b -> X_s S` normalization over `B+`, `B0`, and `Bs`, with two-body `B -> K S` recoil -> weighted `S` four-vector CSVs. |
| `acceptance.py` | scalar decay engine → best-two-track → shared reco → PR#13 `selection_mask`; folds the visible BR in via the decay outcome (neutral sub-modes fail). |
| `run_sensitivity.py` | driver: produce → coupling scan (`N_signal >= 3`) → island CSV → plot. |
| `uncertainty_band.py` | independent full-statistics propagation of 109 coherent FONLL bottom grids, a fresh Winkler-vs-LO-ChPT/spectator decay-model run, and two excluded fresh-seed central numerical controls. |
| `plot_exclusion.py` | `(m_S, sin^2 theta)` island + competitor overlays. |
| `tests/test_model.py` | model-layer validation against published numbers. |

## Run

```
# Canonical publication configuration: first build the 1.2M-event/mass
# importance pools, then reuse them in the coupling scan.
PY=/Volumes/sandbox/conda/envs/llpatcolliders_FONLL/bin/python
$PY -m scalar.production --n-pool 400000 \
  --high-pt-tilt-scale 5 --nominal-mixture-fraction 0.5
$PY -m scalar.run_sensitivity --n-pool 400000 --resume

# Development-size run, a few masses, or re-plot only
python -m scalar.run_sensitivity
python -m scalar.run_sensitivity --masses 0.5 1.0 2.0 --n-pool 100000
python -m scalar.run_sensitivity --plot-only

# full theory-variation campaign (resumable per mass; requires embreex 4.4.0)
BC4_UNCERTAINTY_DIR=/Volumes/GRENDEL/extra_space/bc4_uncertainty \
  python -m scalar.uncertainty_band run --workers 2 \
    --grid-dir /path/to/fonll-local/output
python -m scalar.uncertainty_band status
python -m scalar.uncertainty_band collect

python -m pytest scalar/tests
```

Outputs land in `scalar/tmp/` (4-vector CSVs, `bc4_island.csv`,
`bc4_exclusion.{png,pdf}`). The canonical published curve (the one the paper
quotes) is committed at `scalar/data/published/` (CSV + MANIFEST; see its
README for the re-publish procedure).

The uncertainty campaign atomically checkpoints four-vectors, geometry, and
one result per variation/mass under `BC4_UNCERTAINTY_DIR` (by default the
external `/Volumes/GRENDEL/extra_space/bc4_uncertainty` scratch tree). Once a
complete variation has been checksummed, its raw vectors and geometry are
reclaimed while the compact results, stage tree hashes, seeds, provenance, and
log remain. The pinned Embree backend was cross-checked ray by ray against the
triangle intersector. Stable publication products are in
`scalar/data/published/bundle/`; its README records the non-probabilistic
single-source envelope prescription and limitations.

## References & caveats

* Winkler, Phys. Rev. D 99 (2019) 015018, arXiv:1809.01876 — scalar widths
  (eq. 12, 15, 21, 30-33), production (eqs. A2-A9).
* Unified FIP calculation arXiv:2311.00507 (used by the 2025 PBC report
  arXiv:2505.00947) — BC4 conventions and competitor curves.

The central decay model uses the supplied Winkler dispersive table, and
production uses the inclusive `B -> X_s S` rate. `EXTERNAL_INPUTS_NEEDED.md`
records the resolved input history and the centralized comparison-plot
policy. The analytic LO-ChPT/spectator calculation is retained as
an alternate-model uncertainty, not used for the central curve. The `K -> pi S`
production channel is implemented in the model layer but kaon-flux sampling at
the LHC IP is deferred.
