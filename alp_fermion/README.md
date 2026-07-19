# `alp_fermion/` — PBC benchmark BC10 (fermiophilic ALP) for GRENDEL

A pseudoscalar axion-like particle `a` with **universal coupling to SM
fermions** (single parameter `1/f`, Wilson coefficients `c_f = 1`), produced in
**B → K⁽ⁱ⁾ a** (the flavour-violating `b → s a` penguin, summed over the full
kaon tower as in GKOZ arXiv:2310.03524) and decaying to leptons and hadrons
using the arXiv:2501.04525 SensCalc tables. This is the genuine
charged-track "ALP" benchmark.

BC10 is **coupling-controlled**: the single `1/f` fixes the production yield, the
lifetime cτ, and the visible branching ratios *simultaneously*. The limit is
therefore set **model-completely** (like the HNL `|U|²`, not like the
model-agnostic `(BR, cτ)` scan of `higgs/`): for each `(m_a, 1/f)` the actual
cτ, production, and visible BR are evaluated and `N_signal ≥ 3` (background-free,
3000 fb⁻¹) is required, producing a **closed island** in the `(m_a, 1/f)` plane.

## What is reused (imported, never copied)
- `higgs/grendel_geometry.py`, `higgs/reco_common.py` — the shared GRENDEL mesh
  and the 4-hit two-track vertex reconstruction (`wall_inner_outer`,
  `reconstruct_3d`, `timing_chi2_4hit`) + PR#13 selection.
- `hnl/analysis/decay_reco_acceptance.py` — `build_event_mc`, `best_two_directions`,
  `selection_mask`, `scan_u2` (the acceptance MC + lifetime reweighting).
- `hnl/analysis/_engine.py`, `format_bridge.py`, `exclusion.py` — geometry
  ray-cast, CSV→kinematics bridge, island extraction.
- `hnl/production/fonll` + `decay_engine/kinematics.py` — FONLL bottom sampler
  (14 TeV) and the two-body decay sampler for B → K a.

The key trick that makes the reuse exact: map the coupling to the HNL scan
variable `u2 ≡ (1/f / 1/f_ref)²`. Off the single reference point
(`model.INV_F_REF`), production ∝ `u2` (since BR(B→K a) ∝ (1/f)²) and
cτ ∝ `1/u2` (since Γ_tot ∝ (1/f)²) — exactly the structure `scan_u2` already
implements. This exact scaling is why the shared acceptance and lifetime scan
can be reused. BC10-specific production, decay tables, full-branching
templates, matrix-element weights, pole handling, and component extraction
supply the model inputs around that shared core.

## Layout
```
alp_fermion/
  model.py          model layer: a→ff/γγ/hadrons widths (SensCalc 2501 tables), Γ_tot, cτ, BRs,
                    BR(B→K⁽ⁱ⁾ a) over the kaon tower (b→s a penguin, one-loop RG coefficient)
  alp_production.py B (FONLL) → B→K⁽ⁱ⁾ a → a four-vector CSVs (weight,E,px,py,pz)
  production_spectra.py  source-pinned SensCalc LHC light-parent sampler/audit
  generate_decay_templates_pythia.py  full-branching stable-particle templates
  templates.py      legacy two-track proxy, retained only for reproducibility
  exclusive_decays.py  decoded central/structural exclusive-channel to PDG mapping
  decay_models.py   named 2501 central and 2310 structural model registry
  mass_grid.py      shared dense-scan mass grid
  sensitivity.py    (m_a, 1/f) closed-island scan; N_signal ≥ 3, 3000 fb⁻¹
  plot.py           BC10 island in the (m_a, 1/f) plane (+ optional overlay curves)
  paths.py          output-dir policy (ALP_TMP_DIR override) + hnl/ import shim
  data/alpinist/    historical digitized GKOZ widths plus production-RG provenance
  data/senscalc_2501/ provenance and production audit for the pinned arXiv:2501.04525 upgrade
  data/senscalc_2310/ exact arXiv:2310.03524 structural widths, BRs, and MEs
  data/published/   canonical published sensitivity curve (CSV + MANIFEST; see its README)
  tools/compute_cbs_alpinist.py  regenerates the b→s a RG coefficient (BSD-3 ALPINIST port)
  tools/export_senscalc_2501.{py,wls}  verifies and exports the three SensCalc decay inputs
  tools/export_senscalc_2501_production.{py,wls}  verifies/decodes the LHC production inputs
  tests/test_model.py
  EXTERNAL_INPUTS_NEEDED.md   (input provenance, residual systematics, overlay policy)
```

## Run
```
PY=/Volumes/sandbox/conda/envs/llpatcolliders_FONLL/bin/python
PYROOT=/Volumes/sandbox/projects/aaaPHYSICSaaa/.venvs/fairship/bin/python
$PY -m alp_fermion.alp_production --n-pool 1200000 \
  --high-pt-tilt-scale 5 --nominal-mixture-fraction 0.5
$PYROOT -m alp_fermion.generate_decay_templates_pythia --n-templates 20000
$PY -m alp_fermion.sensitivity
python -m pytest alp_fermion/tests -q
```
Outputs go to `alp_fermion/tmp/` (gitignored): `analysis/bc10_sensitivity.csv`
and `analysis/bc10_island.{png,pdf}`.

The publication diagnostic campaign is driven by
`python -m alp_fermion.run_uncertainty_campaign`. It performs an independent
600k production, geometry, reconstruction, and sensitivity run for every
coherent FONLL scale/PDF/bottom-mass member, plus the explicit `a -> gg`
surrogate and `C_bs` scheme variations. It also runs the exact SensCalc 2310
decay model as a separate one-sided structural contour, outside the ordinary
pointwise envelope. These products are retained for audit and discussion; the
primary proposed-experiment comparison plots remain central-only. See
[`UNCERTAINTY.md`](UNCERTAINTY.md) for the prescription and reproduction
commands.

Sensitivity rows are checkpointed after every mass. For a long or partitioned
campaign, select a checkpoint with `--output path.csv`; rerun the same command
with `--resume` to skip masses already present in that file.

## Model layer references
- Coupling convention + a→f f̄ widths: Bauer–Neubert–Thamm, JHEP 12 (2017) 044
  (arXiv:1708.00443). Our 1/f axis is the BNT `g_aff = c_f m_f/f` convention
  (= 2/f in the GKOZ normalisation).
- Production phenomenology: GKOZ, "ALPs with universal fermion
  couplings — revisited" (arXiv:2310.03524), via the ALPINIST implementation
  (arXiv:2105.10806, BSD-3, pinned SHA in `data/alpinist/COMMIT_SHA.txt`):
  one-loop RG b→s a coefficient (`tools/compute_cbs_alpinist.py`), kaon-tower
  form factors (Boiarska et al., arXiv:1904.10447).
- Decay widths and exclusive branching ratios: arXiv:2501.04525 through the
  pinned SensCalc v1.3.3 binary inputs and reviewable exports in
  `data/senscalc_2501/`. Stable daughter decays and parton hadronization use
  Pythia 8.317. Three-body primary decays are generated in flat phase space
  and reweighted with the exact exported SensCalc squared matrix elements,
  normalized within each exclusive channel so its branching fraction is
  unchanged. The `a -> gg` mode uses an equal light-quark jet surrogate and
  remains an explicit decay-acceptance systematic.
- Heavy-pseudoscalar structural comparison: exact arXiv:2310.03524 widths,
  exclusive branching functions, and squared matrix elements from the same
  SensCalc pin, exported under `data/senscalc_2310/`. This alternative is
  propagated through independent full-decay templates and reconstruction and
  published separately rather than interpreted as a confidence interval.
- Benchmark definition: 2025 PBC report (arXiv:2505.00947); unified FIP
  calculation arXiv:2311.00507; state-of-the-art hadronic treatment
  arXiv:2501.04525.
See `EXTERNAL_INPUTS_NEEDED.md` for residual systematics (production
normalisation <~20% in 1/f; >3 GeV width from the matched perturbative
continuation) and the curated overlay policy.
