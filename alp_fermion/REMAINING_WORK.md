# Remaining work — BC10 fermiophilic ALP

This is the authoritative inventory of work still needed before the generated
GRENDEL BC10 exclusion island can be presented as publication-grade physics
results. It is derived from the code paths that produce the curve, not from
comments or earlier planning notes.

**Scope.** The `alp_fermion/` package only: the PBC BC10 pseudoscalar ALP with
universal SM-fermion couplings (`c_f = 1`), limit set on the single inverse
decay constant `1/f` in the BNT (arXiv:1708.00443) convention. Production is
`B -> K^(i) a` via the `b -> s a` penguin; decays use the pinned SensCalc
v1.3.3 arXiv:2501.04525 tables. The analysis core (geometry, reconstruction,
acceptance MC, island extraction) is *imported* from `higgs/` and
`hnl/analysis/`, never copied, so the P0 items below are shared with the other
GRENDEL benchmarks but have been re-verified here against the BC10 call chain.

The inventory covers three kinds of work:

- **P0 — result definition:** missing detector, background, or statistical
  inputs that can change the *meaning* of the exclusion.
- **P1 — uncertainty & completeness:** missing channels and uncertainty
  propagation needed for defensible central curves and bands.
- **P2 — validation & reproducibility:** convergence, provenance, and workflow
  checks needed to make the result repeatable and auditable.

An item is not complete merely because a central constant exists in the code.
It is complete when its source, uncertainty, propagation, and regression test
are present.

## Package status and what is already bounded

The 2026-07-15 publication campaign (`data/published/MANIFEST.json`) is a
145-point mass grid with a 1,200,000-event high-pT importance-sampled FONLL
bottom pool, 20,000 full-branching Pythia 8.317 templates per supported mass,
exact three-body matrix-element reweighting, and 60 decay/reconstruction
samples per detector-entering ALP. 103 of 145 grid points are sensitive.

The uncertainty campaign (`uncertainty_campaign.py`, `combine_uncertainty_band.py`)
propagates each nuisance through a **complete independent signal run** — no
event-level reweighting. Its axes and counts are asserted at
`uncertainty_campaign.py:95` and `combine_uncertainty_band.py:220`:

| axis | count | content |
|---|---|---|
| `central` | 1 | reference contour |
| `scale` | 6 | non-central coherent 7-point `(mu_R, mu_F)` bottom grids |
| `pdf` | 100 | NNPDF4.0 NLO replicas (member 0 excluded, `uncertainty_campaign.py:65`) |
| `mb` | 2 | `m_b = 4.5` / `5.0 GeV` around the `4.75 GeV` central |
| `decay_gg` | 3 | pure-`u` / pure-`d` / pure-`s` gluon surrogates |
| `cbs` | 2 | `+/-20%` `C_bs` amplitude (`CBS_SCHEME_RELATIVE_AMPLITUDE = 0.20`, `uncertainty_campaign.py:17`) |
| `decay_structure` | 1 | exact arXiv:2310.03524 model, one-sided, outside the halo |
| `numerical_control` | 2 | same-physics central repeats, excluded from the envelope |

The combination is explicitly a **one-source-at-a-time variation envelope, not
a confidence interval** (`UNCERTAINTY.md`, and the `method.combination` field
written at `combine_uncertainty_band.py:405`). Nothing is added in quadrature.
That is an honest and correctly labelled construction; the items below do not
ask for it to be relabelled.

### Resolved / not applicable — do not re-import from the HNL inventory

These are recorded so the HNL caveats are not copied onto BC10 by mistake.

1. **The charm `mu_F` PDF-extrapolation artifact does not apply.** BC10
   production is bottom-only: `alp_production.py:116,121` sample the FONLL
   `"bottom"` pool and take `get_sigma_total("bottom")`. Charm is never
   sampled, and `uncertainty_campaign.py:61` drops every non-bottom grid from
   the variation ensemble. `NNPDF40_nlo_as_01180` has `QMin = 1.65 GeV`, and
   charm's `mu_F = 0.5*sqrt(m_c^2 + pT^2)` dips below it for `pT < 2.94 GeV`,
   which is what distorts the HNL charm scale envelope's lower edge. Bottom's
   `mu_F = 0.5*sqrt(m_b^2 + pT^2) >= 2.375 GeV` always clears `QMin`, so BC10 is
   immune **by construction**, not by luck. The bottom scale envelope
   (`~+/-42%` at grid level) is a genuine perturbative uncertainty.

2. **There is no tau production chain, so tau spin correlations do not apply.**
   `alp_fermion/` contains zero references to `induced_tau` or `tau_decay`
   (verified by grep over the package). The tau machinery under `hnl/` is
   vendored shared infrastructure that the ALP package never calls. The HNL's
   "tau spin correlations" item is not in scope here. *(The ALP may still decay
   to taus — `channel_003 = (-15, 15)` in `exclusive_decays.py:27`. That is a
   separate question, treated on its own merits in item 12.)*

3. **The `a -> gg` surrogate is bounded and transparent — not an open gap.**
   Pythia's external-decay interface cannot fragment an isolated colour-singlet
   gluon pair, so `generate_decay_templates_pythia.py:167-180` replaces the
   imported `a -> gg` branching fraction (`channel_018`) with an equal `u/d/s`
   light-quark mixture. This is **already propagated**: `--gluon-surrogate`
   generates the pure-`u`, pure-`d`, and pure-`s` alternatives, each of which
   gets its own full geometry + reconstruction + sensitivity run
   (`uncertainty_campaign.py:105-117`), and they are retained as three **named
   full-simulation variations** rather than collapsed into a confidence
   interval. The surrogate changes decay *acceptance* only; the imported total
   width and branching ratios are untouched. Treat this as adequate and
   transparent. A better fix would require either a colour-connected ALP
   resonance inside Pythia or a physical exclusive hadronic model, and would
   mean regenerating decay templates (not production pools) — worth doing only
   if the `gg_u`/`gg_d`/`gg_s` spread ever becomes a leading term.

## P0: define the detector-level result

These three items are shared with every GRENDEL benchmark through the common
`higgs/` geometry + reconstruction and the shared exclusion code. **Verified
against the BC10 code path**, not assumed: `sensitivity.py:56-66` imports
`N_THRESHOLD`/`L_INT_PB`/`CMS_ORIGIN` from `hnl/analysis/constants.py`,
`build_event_mc`/`scan_u2` from `hnl/analysis/decay_reco_acceptance.py`, and
`find_exclusion_band_refined` from `hnl/analysis/exclusion.py`; the mesh comes
from `hnl/analysis/_engine.py::_get_mesh` over `higgs/grendel_geometry.py`.

### 1. Add detector response to the decay/reconstruction acceptance

**Current code:** `sensitivity.py:174` runs `build_event_mc` per detector-
entering ALP: it samples `DECAY_SAMPLES = 60` decay vertices along the in-volume
flight path (`sensitivity.py:73`), draws a full-branching Pythia template, boosts
to the lab, keeps the two highest-momentum charged stable daughters, and runs
the shared `higgs/reco_common` bounded 4-hit reconstruction plus the PR #13
selection (gate/pointing/collinearity/timing). Because the templates sample the
full exclusive mixture including neutral modes, reconstruction itself supplies
the visible fraction — `sensitivity.py:192` sets `branching_factor = 1.0` when
`includes_full_branching` is set, rather than applying an inclusive `BR_vis`.

Detector response is **geometric wall hits with Gaussian position/time smearing
only**: `higgs/reco_common.py:110-115` smears hit positions by `sigma_hit` in
the two transverse directions and `reco_common.py:54` smears times by
`SIGMA_T_DEFAULT = 0.5e-9 s` (`reco_common.py:30`). There is no material
interaction, tracking or vertexing inefficiency, particle threshold,
trigger/readout, pileup, dead region, or occupancy anywhere in the chain
(verified by grep over `higgs/reco_common.py`, `higgs/grendel_geometry.py`, and
`hnl/analysis/decay_reco_acceptance.py`). The Pythia templates are generator-
level stable final states with no detector simulation.

**Required work:**

- implement detector response beyond geometric hits + smearing: tracking,
  vertexing, particle thresholds, material interactions, trigger/readout, and
  event-selection efficiencies;
- include pileup, timing, dead regions, and occupancy if they affect the
  proposed detector;
- assign per-channel efficiency uncertainties. BC10 needs this per *decay
  channel*: the island spans `a -> ee/mumu/tautau` and 28 exclusive hadronic
  modes whose track multiplicities differ sharply, so a single flat efficiency
  is not a defensible summary.

**Completion test:** an efficiency map or detector simulation with versioned
inputs reproduces benchmark samples, and `sensitivity.py` consumes per-channel
efficiencies with uncertainty variations.

### 2. Validate the geometry and detector configuration

**Current code:** `higgs/grendel_geometry.py` builds a closed tunnel mesh from
one survey polyline with a fixed `Y_POSITION = 22` m vertical position
(`grendel_geometry.py:307`) and a fixed `DETECTOR_THICKNESS = 0.24` m wall inset
(`grendel_geometry.py:40`), classifying tracker vs scintillator surfaces. Rays
originate at `CMS_ORIGIN = (0, 0, 0)` (`hnl/analysis/constants.py:29`), used by
`sensitivity.py:100`. The mesh contains no supports, services, inactive regions,
material, alignment uncertainty, or configurable layout.

**Required work:**

- confirm the survey coordinates, CMS IP transform, tunnel cross-section, and
  proposed active volume with the detector/design owners;
- version the geometry and expose layout parameters rather than relying on
  import-time constants;
- model inactive and inaccessible regions and evaluate alignment/survey
  variations;
- independently compare ray intersections against a reference implementation.

**Completion test:** signed-off geometry inputs, reference intersection tests,
and a geometry-variation envelope propagated to the island.

### 3. Supply backgrounds and a statistical model

**Current code:** exclusion is `N_signal >= N_THRESHOLD` with
`N_THRESHOLD = 3.0` (`hnl/analysis/constants.py:21`), a zero-background 95% CL
approximation applied at `hnl/analysis/exclusion.py:46,122,130` and consumed by
`sensitivity.py:214`. Luminosity is a fixed `L_INT_PB = 3e6` (3000 fb^-1,
`constants.py:12`). There is **no background count in the statistical model**,
no control region, no nuisance parameter, no systematic uncertainty, and no
coverage calculation. Signal passes the PR #13 selection, but no background rate
after that selection enters the limit.

Note the BC10-specific consequence: the island's **upper** `1/f` edge is set by
the ALP decaying before reaching PX56, and the **lower** edge by too little
production (`sensitivity.py:218-228`). A background term would move the lower
edge inward but leaves the upper edge essentially untouched, so the two
boundaries do not degrade symmetrically once backgrounds exist. Do not assume a
uniform contour shift.

**Required work:**

- estimate beam-, collision-, cosmic-, neutrino-, and detector-induced
  backgrounds after the final selection;
- define control samples and background uncertainties;
- include luminosity, reconstruction, trigger, geometry, signal-model, and
  finite-MC nuisance parameters;
- choose and implement the statistical prescription (CLs or a documented
  Poisson construction) and validate its expected coverage;
- define how open island boundaries are reported (`invf_min_open` /
  `invf_max_open` already carry the open-edge state).

**Completion test:** the final likelihood/configuration is versioned, has toy or
asymptotic validation, and produces expected limits with uncertainty bands.

## P1: production model and its uncertainties

### 4. Bound the FONLL pT and rapidity truncation

**Current code:** every bottom grid stops at `pT = 50 GeV` and `|y| = 3` — the
bound is baked into the filename stem
`..._dsdpTdy_pt0-50_y-3to3` (`uncertainty_campaign.py:13-16`), which
`discover_fonll_variations` uses to locate each grid
(`uncertainty_campaign.py:69`). Nothing in the package bounds the contribution
from beyond those limits: a grep for `truncat`/`extrapolat`/`tail` across
`alp_fermion/` returns only the unrelated high-pT proposal note in
`data/published/README.md:62`. The GRENDEL acceptance is a tunnel *above* IP5,
so the accepted sample is not forward-dominated and the `|y| < 3` cut is
plausibly benign — but "plausibly benign" is not a bound, and this is stated
here as unquantified rather than negligible.

**Required work:**

- extend the `pT` and rapidity grids, or place a quantitative bound on the
  omitted contribution *after* GRENDEL acceptance (the production-level
  fraction alone is not the answer — acceptance must be folded in);
- state the bound in the published manifest next to the grid bounds.

**Completion test:** a documented accepted-yield bound on the `pT > 50 GeV` and
`|y| > 3` regions, or extended grids covering them.

### 5. Add a tracked alpha_s companion axis — negligible, for prescription completeness only

**Current code:** the campaign has **no alpha_s axis**. `GRID_STEM`
(`uncertainty_campaign.py:13-16`) hardcodes `as_01180`, and
`discover_fonll_variations` rebuilds every grid path from that stem
(`uncertainty_campaign.py:69`), so the `as_01170`/`as_01190` companion grids are
structurally unreachable even though they exist in the FONLL workspace. The
axis list at `uncertainty_campaign.py:92-97` and the assertion at
`combine_uncertainty_band.py:220` admit only `central/scale/pdf/mass`. This is
the exact analogue of the HNL's `run_variation_band.py` stem hardcoding.

**Measured size — do not over-prioritize.** The PDF4LHC alpha_s half-spread on
the **bottom** trapezoid-integrated grid is **1.49%**, against a bottom scale
envelope of `~+/-42%` and a bottom PDF 1sigma of `3.9%`. Folded in quadrature
with the terms it would join, it moves the band by `sqrt(42^2 + 1.49^2) - 42
= 0.03` percentage points — below the campaign's own numerical-control spread
and invisible on any plot. Run it for PDF4LHC prescription completeness and
provenance, **not** for physics. Unlike the HNL, BC10 has no charm `mu_F`
artifact competing for the same two-chain cost (see "Resolved" item 1), so
there is no larger item this one should queue behind — it is simply low value.

**Required work:**

- parametrize the grid stem per manifest entry instead of hardcoding
  `as_01180`, and read each grid's `path` field from the manifest;
- synthesize two coherent variation dicts for the companions (alpha_s is
  coherent across quarks, so it is structurally a scale point);
- run the two companion production+analysis chains and fold them in as a
  PDF4LHC half-difference.

**Completion test:** the two alpha_s companion contours are produced by tracked
code from the manifest and recorded in `UNCERTAINTY_MANIFEST.json`.

### 6. Do not implement `B_s -> phi a` without first sourcing a form factor

**Current code:** `PRODUCTION_PARENTS` is `(("B+", 521), ("B0", 511))`
(`model.py:548-551`) crossed with the nine-state `KAON_TOWER`
(`model.py:402-413`). `B_s -> phi a` is **absent**, and `model.py:59` records
this as deliberate: *"(Bs -> phi a is not included, matching ALPINIST.)"* The
channel is likewise absent from the pinned GKOZ (arXiv:2310.03524), ALPINIST,
and SensCalc v1.3.3 channel sets, so **adding it is not "finishing the
arXiv:2501.04525 implementation"** — it is new physics input this package would
have to source and own.

`EXTERNAL_INPUTS_NEEDED.md` §1 states the blocker precisely: Boiarska et al.
(arXiv:1904.10447) supplies the kaon-tower form factors for `B+`/`B0` but **no
corresponding `B_s` tower**. An earlier "about 10% yield" estimate is **not a
validated calculation**. The code supports the scepticism: `FRAG_B`
(`hnl/production/constants.py:66-70`) gives `f_s = 0.0883` against
`f_u + f_d = 0.724`, i.e. `f_s/(f_u+f_d) = 0.122`; with a ground-state
`B_s -> phi` rate comparable to `B -> K*(892)`, that points to a few-percent
yield correction, and a few-percent yield is a few-percent contour displacement
on a log-scale reach boundary — invisible.

**Required work (all of it, before writing any code):**

- establish a **citable** `B_s -> phi` form-factor source and compute the
  mass-dependent rate;
- normalize it consistently in the BNT convention at `Lambda_UV = 1 TeV`
  (`CBS_EFF = 3.518383e-4`, `model.py:118`);
- apply the `B_s` fragmentation fraction and two-body kinematics;
- rerun production and sensitivity.

**Completion test:** either a sourced, normalized, mass-dependent
`B_s -> phi a` rate in the tower with a production+sensitivity rerun, or a
committed note quantifying the omission against the accepted `B+`/`B0` yield.
Do not implement on the strength of the 10% estimate.

### 7. Close or bound the omitted bottom-hadron production fractions

**Current code:** production weights use `frag = FRAG_B[pdg]`
(`alp_production.py:153`) for `521` and `511` only. From
`hnl/production/constants.py:66-80`, that is `0.362 + 0.362 = 0.724` of the
bottom fragmentation. The remainder carries **no ALP production at all**:
`f_s = 0.0883` (item 6) and `OMITTED_FRAG_B["b_baryons"] = 0.1875` — together
**~27.6%** of bottom fragmentation. The b-baryon share is the larger of the two
and is not mentioned anywhere in `alp_fermion/`.

`Lambda_b -> Lambda a` is the same `b -> s a` penguin driving the kaon tower, so
it is not obviously negligible on physics grounds. **I could not determine from
the code or the committed provenance whether GKOZ/ALPINIST/SensCalc include a
b-baryon ALP channel**; `model.py:59` and `EXTERNAL_INPUTS_NEEDED.md` §1 address
only `B_s -> phi a` and are silent on baryons. This needs a literature check
before it can be either implemented or dismissed — it is stated here as an
open, unquantified omission, not as a known defect.

**Required work:**

- determine whether the pinned references provide a `Lambda_b -> Lambda^(i) a`
  form factor set; if so, implement it with the fragmentation fraction above;
- if not, bound the omission after GRENDEL acceptance and record the bound;
- state the closure explicitly in the manifest so 27.6% of bottom fragmentation
  is not silently absent.

**Completion test:** the bottom fragmentation accounting closes, and every
omitted component has either a generated channel or a documented bound on the
island.

### 8. Propagate form-factor and fragmentation uncertainties

**Current code:** two distinct gaps, neither currently varied.

*Form factors.* The kaon-tower amplitudes `_M_BP`/`_M_BS`/`_M_BV`/`_M_BA`/`_M_BT`
(`model.py:431-484`) hardcode the Boiarska et al. parameters as bare literals —
e.g. `f_0, m_fit = 0.33, 6.12` (`model.py:432`), the `K1` mixing angle
`th = -0.593` (`model.py:464`). No uncertainty, covariance, or source/version
table accompanies any of them. The campaign varies the *overall* `C_bs`
normalization by `+/-20%` (`uncertainty_campaign.py:118-133`), which covers the
penguin coefficient — `EXTERNAL_INPUTS_NEEDED.md` §1 bounds the `CBS_EFF`
scheme residual at `<~20%` in `1/f`, so that axis is adequately sized for what
it targets — but it does **not** cover per-channel form-factor shape
uncertainty, and the two are not the same nuisance.

*Fragmentation.* `FRAG_B[521]` and `FRAG_B[511]` are byte-identical constants
(`0.36205648081100655`) with no covariance and no `pT`/rapidity dependence.
Both parents also share **one** bottom FONLL shape: `alp_production.py:155`
calls `meson_4vec_from_kinematics(..., m_B)`, which re-derives `E` and `pz` at
the species mass from a single shared pool (documented as the "shared-shape
approximation" at `hnl/production/fonll/meson_sampler.py:48-51`). For `B+` vs
`B0` this is a mild approximation — their masses differ by ~0.3% — so it is far
less consequential here than in the HNL, where the same pool serves `Bs` and
`Lambda_b`. Flagged for completeness, not as a leading term.

**Required work:**

- attach sources, uncertainties, and correlations to the form-factor
  parameters, and propagate at least a per-channel normalization variation;
- propagate the experimental `FRAG_B` covariance;
- confirm or bound the `B+`/`B0` shared-shape approximation.

**Completion test:** form-factor and fragmentation variations appear as named
axes in `UNCERTAINTY_MANIFEST.json` with a documented source table.

### 9. External production-channel closure — audited, keep it that way

**Current code:** already in good shape; recorded so it is not redone.
`EXTERNAL_INPUTS_NEEDED.md` §3 and `data/senscalc_2501/PRODUCTION_AUDIT.md`
audit the broader arXiv:2501.04525 production policy (Drell-Yan, gluon fusion,
proton bremsstrahlung, light-meson decays, quark fragmentation) against the
GRENDEL acceptance, with `production_spectra.py` as the source-pinned sampler.
Findings: Drell-Yan negligible across the island; decoded generalized
fragmentation `<= 0.15%` of the B tower outside the pole windows; light-meson
rates over five orders below accepted B production at low mass. These are
quantified, not assumed.

**Required work:** re-run the audit if the mass grid extends or the island
moves; keep `PRODUCTION_AUDIT.md` in sync with the pinned SensCalc version.

**Completion test:** the audit's quantitative bounds still hold for the
published grid.

## P1: decay model

### 10. Unclassified neutral remainder — MEASURED 2026-07-17, negligible, closed

**Current code:** this is the decay-side finding most likely to be overlooked.
`exclusive_decays.py:100` computes
`weights[UNCLASSIFIED_CHANNEL_ID] = max(1.0 - known, 0.0)` — whatever exclusive
branching fraction the 32-channel table does not account for at a given mass is
swept into a single synthetic channel. At template generation
(`generate_decay_templates_pythia.py:184-186`) that remainder is materialized as
`(22, 22)` — a **diphoton**, i.e. zero charged tracks, therefore invisible to
the reconstruction.

This is deliberately conservative (`exclusive_decays.py:86` calls it "a
conservative neutral missing-width mode"): the remainder is assumed to
contribute no signal. But its **size is never reported and never varied**. If
the remainder is large at some masses, the island is being suppressed there by
an amount nobody has measured, and the conservatism is unquantified rather than
bounded. No axis in the campaign touches it.

**Resolution (2026-07-17): tabulated, negligible, item closed.** The
tabulation called for below was run over the full `ALP_MASS_GRID` (145 masses,
0.22--4.75 GeV) via `exclusive_branching_weights`:

- **maximum remainder 2.1%**, at `m_a = 2.15 GeV`;
- only **3 of 145** masses carry a remainder above 1%.

That maximum sits below the 4.1% median production-MC noise measured in the HNL
`audits/curve_impact_20260610/`, so the invisible-diphoton conservatism cannot
move the island above statistical noise anywhere on the grid. The conservatism
is now **bounded rather than unquantified**, which was the whole concern. No
variation axis is warranted.

Reproduce with:

```python
from alp_fermion.exclusive_decays import exclusive_branching_weights as w, UNCLASSIFIED_CHANNEL_ID as U
from alp_fermion.mass_grid import ALP_MASS_GRID as G
print(max(((w(m).get(U, 0.0), m) for m in G)))   # -> (0.021, 2.15)
```

**Residual work (optional, low priority):** persist the per-mass remainder
table as a published artifact so the bound is auditable from the bundle rather
than recomputed, and add a regression test asserting the maximum stays below a
declared threshold (say 5%) if the channel table is ever edited.

**Completion test:** met for the current channel table — the maximum is
recorded above. Re-run the one-liner if `exclusive_decays.py` changes.

### 11. Decide whether the 2310 structural comparison is the decay uncertainty

**Current code:** the only decay-model alternative is the exact arXiv:2310.03524
model (`decay_structure_variation`, `uncertainty_campaign.py:137-153`), and it
is deliberately excluded from the pointwise halo and published as a **one-sided
named structural comparison** — `combination_role` is
`"one_sided_structural_comparison_outside_halo"` (`uncertainty_campaign.py:152`).
`UNCERTAINTY.md` is explicit that it is "not a calibrated uncertainty". It is
run twice: once on the campaign grid and once directly on the refined
high-statistics central grid (`record_dense_structural.py`, the
`dense_decay_structure` block at `combine_uncertainty_band.py:360-367`), with
topology changes recorded rather than converted into displacements.

This is well-engineered and honestly labelled. What is **missing** is any
calibrated decay-side uncertainty: the widths, exclusive BRs, and matrix
elements from the pinned SensCalc 2501 tables enter with **no uncertainty at
all**. Since `1/f` controls production, lifetime, and visible BRs
simultaneously, a width error propagates to both island edges coherently, and
the lifetime leg is the more sensitive one — this is not a normalization-only
nuisance.

**Required work:**

- decide and document whether the 2310 contour *is* the decay uncertainty
  statement for the paper, or whether a calibrated band is required;
- if calibrated: source width uncertainties for the 2501 tables and propagate
  through both the lifetime leg (`ctau_at_reference_coupling`,
  `exclusive_decays.py:104`) and the composition leg;
- either way, state the choice in the published manifest.

**Completion test:** the manifest names the decay-uncertainty prescription, and
the paper text matches the data product's own labelling.

### 12. Validate the tau and charm decay legs

**Current code:** `a -> tau tau` (`channel_003`, `exclusive_decays.py:27`) and
`a -> c cbar` (`channel_019`, `exclusive_decays.py:43`) are passed to Pythia as
primary products; Pythia then decays the taus and hadronizes the charm pair
(`generate_decay_templates_pythia.py:181-192`). The ALP is declared as a scalar
resonance with `isResonance = false` and decays added as flat channels
(`generate_decay_templates_pythia.py:163-192`), so **Pythia has no information
about the parent's pseudoscalar nature when it decays the taus** — the tau
spin correlations implied by `a -> tau tau` are not communicated.

I could not determine from the code how much this matters. The effect is on the
charged-track multiplicity and momentum spectrum of tau decays, which feeds the
best-two-track selection. It is bounded above by the `tautau` branching
fraction, which is only open above `2 m_tau = 3.55 GeV` — **above the published
island's high-mass closure at `3.29 GeV`** (`MANIFEST.json`,
`headline_reach.mass_hi_closure_GeV`). On the current grid this is therefore
almost certainly irrelevant; it would only matter if the island extended above
`3.55 GeV`.

**Required work:**

- if the mass range is extended above `~3.55 GeV`, validate the tau leg against
  a spin-aware generator (TAUOLA or Pythia with correlations enabled);
- otherwise record that `tautau` is closed below the island's high-mass edge
  and close this item.

**Completion test:** either a documented statement that the `tautau` threshold
lies above the sensitive range, or a spin-correlation validation.

## P2: validation and reproducibility

### 13. Pin Pythia machine-readably in the template bundles

**Current code:** provenance is strong but has one narrow hole. Template bundles
record `n_templates`, `seed`, `decay_model`, `gluon_surrogate`, `flavor`,
`ctau_m_u2eq1`, and a `decay_backend` string
(`generate_decay_templates_pythia.py:320-339`), and `_validate_resumed_template`
(`generate_decay_templates_pythia.py:93-114`) refuses to reuse a bundle whose
model/count/surrogate differ. The campaign captures git state — commit, tracked
diff hash, and a clean flag — at `run_uncertainty_campaign.py:86-98`, and
FONLL grids are SHA-256 verified against the manifest
(`uncertainty_campaign.py:72-74`), as are sensitivity CSVs
(`combine_uncertainty_band.py:77-79`).

The hole: **the Pythia version is not stored in the template bundles**. It
appears only as prose in `data/published/MANIFEST.json`
(`physics_inputs.decays`: "Pythia 8.317"). A grep for a Pythia version constant
across `alp_fermion/*.py` returns nothing. So a bundle regenerated against a
different Pythia would pass `_validate_resumed_template` and be silently
consumed.

**Required work:**

- record the Pythia version (and ROOT/`TPythia8` build) in each template `.npz`
  and check it on resume/load;
- commit an environment lock for the ROOT/Pythia and FONLL Python environments;
- add a one-command clean-room reproduction of a representative mass point.

**Completion test:** regenerating a template under a different Pythia version
fails loudly rather than being reused.

### 14. Quantify finite-template statistics

**Current code:** templates are 20,000 events per mass
(`generate_decay_templates_pythia.py:360`), matched per mass label with no
interpolation between masses (`sensitivity.py:147`). The two `numerical_control`
repeats vary the **production** pool and the **reconstruction** RNG offset
(`uncertainty_campaign.py:156-175`) but keep `template_variant: "central"` —
i.e. **they reuse the same central template bundles**, so the published
numerical-control spread does not contain any template-statistics component.
Nothing else varies the template seed. `_load_dense_structural` and the ESS
diagnostics (`sensitivity.py:230-251`, minimum sensitive-edge `event_ess`
~1528 in `MANIFEST.json`) bound importance-weight degeneracy, not template
count.

**Required work:**

- regenerate templates at a different seed for a subset of masses and measure
  the contour shift;
- confirm 20,000 is converged, or raise it;
- confirm the mass-grid spacing needs no template interpolation.

**Completion test:** a declared template-statistics tolerance, measured on real
mass points, recorded next to the numerical controls.

### 15. Key the geometry cache by content, not mtime

**Current code:** `sensitivity.py:89-104` caches ray-cast results as
`geom_{mass}.npz`, keyed **only by mass label**. There is an mtime guard —
`fresh = cache.stat().st_mtime >= source_mtime` (`sensitivity.py:96`) — so a
regenerated production CSV forces a fresh ray-cast. But the cache embeds no
checksum of the input vectors, the geometry version, the origin, or the cuts,
so a *geometry-code* edit does not invalidate it; `--force-geometry` must be
passed by hand (`sensitivity.py:320`). Note the HNL's negative result applies
here too: geometry caches must not be reused across production variations,
since row `i` of a varied pool is not the same trajectory as row `i` of the
central pool. The campaign already avoids this by giving each variation its own
run tree and hashing vector/geometry trees into the completion marker.

**Required work:**

- key the cache by input checksum plus geometry/analysis configuration;
- verify no stale cache from another run tag can be consumed.

**Completion test:** editing the geometry code invalidates the cache
automatically.

### 16. Converge the scan and acceptance settings

**Current code:** the `u2` scan is 300 log-spaced points over `[1e-16, 1]`
(`sensitivity.py:72`), followed by `find_exclusion_band_refined`, which does
bracketed refinement against a callable `evaluate` (`sensitivity.py:206-216`) —
so the boundaries are already root-found rather than read off the grid, which
is better than a fixed-grid readout. `DECAY_SAMPLES = 60` (`sensitivity.py:73`)
and `EVENT_CHUNK = 256` (`sensitivity.py:168`) with deterministic per-chunk
sub-seeds (`sensitivity.py:170-176`) make the published path reproducible for a
fixed chunk size. The two `numerical_control` repeats measure the
production+reconstruction spread and are reported separately with a
`*_repeat_not_subdominant` flag (`combine_uncertainty_band.py:196-215`) that
fires when numerical noise rivals the physical envelope.

What is **not** demonstrated: that `DECAY_SAMPLES = 60` is converged against an
exact quadrature of the decay-in-volume integral, and that the mass grid is
fine enough at the pole windows and threshold tips (`MANIFEST.json` closes
threshold tips by log-yield/log-coupling interpolation).

**Required work:**

- vary `--decay-samples` and seeds; compare against exact quadrature;
- demonstrate stability under mass-grid refinement near the `eta`/`eta'` pole
  windows (`model.py:128-132`) and the high-mass closure;
- do not claim sub-grid precision.

**Completion test:** island boundaries stable under refinement within a declared
tolerance.

### 17. Extend regression coverage to the production path

**Current code:** 1,816 test lines across 12 files. Coverage is strong on the
campaign and data products — `test_uncertainty_campaign.py` (721 lines),
`test_published_curve.py` (158, checks manifest/curve self-consistency and
SHA-256), `test_export_senscalc_2501.py` (188), `test_decay_matrix_elements.py`
(cross-checks the exported `CForm` against Wolfram Engine anchors),
`test_production_spectra.py` (147). `test_model.py` and
`test_production_importance.py` are the only files referencing the production
model (`br_B_to_Ka`/`KAON_TOWER`/`alp_production`).

Gap: there is no regression test for the `alp_production.generate` weight
normalization end-to-end — the `2 * sigma_b * frag * br / n_each` construction
(`alp_production.py:160-164`) and the disjoint channel partitioning
(`_partition_parent_indices`, `alp_production.py:91-98`) are exercised only
indirectly.

**Required work:**

- add a channel-closure test asserting summed reference weights against an
  independent `br_B_tower_total` x `sigma_b` x `frag` calculation;
- assert the parent-index partition is disjoint and exhaustive.

**Completion test:** a factor-2/normalization bug in `generate` fails a test.

### 18. Decide the band's role in the published figure

**Current code:** `plot.py::plot_island` draws the **central island only** —
`fill_between(m, lo_fill, hi_fill, ...)` (`plot.py:118`) fills between
`invf_min` and `invf_max`, i.e. the island itself, **not** an uncertainty band.
Open edges fill to the axis limit and are not drawn as lines
(`plot.py:116-122`), which is honest open-boundary handling. The variation
envelope exists only as CSV data products under `data/published/bundle/`. This
is consistent with `UNCERTAINTY.md`, which states the envelope "is retained as
an audit diagnostic and is not drawn on the primary comparison plots" — a
deliberate choice, not an oversight. External comparison curves live in the
sibling `curves_PBC` repository, not here; `plot.py:128-133` accepts only an
optional long-format overlay CSV for focused validation.

**Required work:**

- confirm with the paper draft whether the published BC10 figure shows a band;
  if yes, wire `bc10_single_source_variation_envelope.csv` into the figure and
  caption it as a one-source-at-a-time envelope, **not** a confidence interval;
- if the 2310 structural contour is shown, keep it visually distinct
  (the campaign already publishes it as a separate dashed data product);
- add detector/background/numerical bands once P0 items 1-3 exist.

**Completion test:** every visible band traces to a stored variation output and
regenerates from the manifest; the caption matches the data product's labelling.

## Recommended execution order

1. **Define the detector-level result** (items 1-3). Until backgrounds and a
   statistical model exist, every other number is a diagnostic of a simplified
   analysis, not an uncertainty on a defined exclusion.
2. ~~**Tabulate the unclassified neutral remainder** (item 10).~~ **Done
   2026-07-17: max 2.1%, below MC noise, item closed.** No decay-acceptance
   axis is missing on this account.
3. **Settle the decay-uncertainty prescription** (item 11) and the tau/charm
   threshold statement (item 12), since these decide what the paper claims.
4. **Bound the omissions**: FONLL truncation (4), b-baryon/`B_s` fragmentation
   closure (7, 6). These are bounding exercises before they are implementation
   work — do the bound first, and only implement if the bound fails.
5. **Propagate the remaining production uncertainties** (item 8).
6. **Close reproducibility**: Pythia pinning (13), template statistics (14),
   cache keying (15), convergence (16), regression coverage (17).
7. **Figure/band decision** (item 18), once the above settle what there is to
   draw.
8. **alpha_s companion axis** (item 5) last, and only for PDF4LHC prescription
   completeness. At 1.49% against a 42% scale envelope it cannot move the
   contour.

Items 6 and 5 are the two most likely to be over-prioritized by someone reading
the HNL inventory across: `B_s -> phi a` looks like unfinished implementation
but is unsourced new physics with a few-percent expected effect, and alpha_s
looks like a missing PDF4LHC term but is numerically negligible for a
bottom-only benchmark.

## Note on the two package copies

At the time of writing, `llpatcolliders_BC10_PR18/alp_fermion/` and
`grendel_alp/alp_fermion/` are **identical** (`diff -rq` reports no difference
in any tracked file; the only difference is an untracked `tmp/` scratch
directory present in the former). Every claim above therefore holds equally in
both copies. If they diverge later, re-verify the `file.py:line` references,
which were read from `llpatcolliders_BC10_PR18/alp_fermion/`.
