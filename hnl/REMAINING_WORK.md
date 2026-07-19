# Remaining work

This is the authoritative inventory of work still needed before the generated
GRENDEL HNL exclusion curves can be presented as publication-grade physics
results. It is derived from the code paths that produce the curves, not from
comments or earlier planning notes.

The inventory covers four different kinds of work:

- **P0 -- result definition:** missing detector, background, statistical, or
  physics-model inputs that can change the meaning of the exclusion.
- **P1 -- uncertainty/completeness:** missing channels and uncertainty
  propagation needed for defensible central curves and bands.
- **P2 -- validation/reproducibility:** numerical convergence, provenance, and
  workflow checks needed to make the result repeatable and auditable.

An item is not complete merely because a central constant exists in the code.
It is complete when its source, uncertainty, propagation, and regression test
are present. The induced-tau central branching fractions in
`production/constants.py` have been manually reviewed against PDG, lattice,
and LHCb inputs, but that review is not yet reproducible from a machine-readable
source/version table. Their provenance, uncertainties, and correlations remain
outstanding below.

A companion diagnostic estimate of the curve impact is stored in
`audits/curve_impact_20260610/README.md`, with detailed results at every
configured mass point and in 0.3 GeV bins. Those ranges are explicitly
classified as calculations, literature-informed ranges, proxies, or unknowable
scenarios; they are not a statistical uncertainty band.

## Curve-visibility triage (production x analysis, by PBC scenario)

The items below are ordered by physics priority, not by their effect on the
*current* curves. As of the 2026-07-17 reconciliation, every PBC scenario --
100 (Ue), 010 (Umu), 001 (Utau) -- has complete production and decay templates
on the published 123-point grid; the live selection in
`analysis/decay_reco_acceptance.py` matches the upstream `higgs/` GRENDEL
reconstruction cut-for-cut, and the exclusion threshold is the standard
zero-background 95% CL value (`N >= 3`). No further large central Monte Carlo
campaign is currently indicated by the convergence checks. The remaining items
still matter: some refine production or decay theory, while detector response,
backgrounds, and statistical modelling can change the *meaning* of the result
and may ultimately move the plotted contour.

| layer \ scenario | 100 (Ue) | 010 (Umu) | 001 (Utau) |
|---|---|---|---|
| **production** | mesons + W/Z EW; Pythia charged-kaon flux + transport proxy | mesons + W/Z EW; Pythia charged-kaon flux + transport proxy | mesons + tau-parent + W/Z EW; tau Kallen & W-tau polarization already fixed |
| **analysis** | cuts/geometry/threshold single-sourced with `higgs/` -- aligned | aligned | aligned |

Robustness fixes that *prevent* a future silent curve corruption (not a current
defect) were applied on 2026-06-16: the geometry cache now invalidates when its
source CSV is newer (mtime guard, item 19); `run_sensitivity.py` records every
skipped `(flavor, mass)` with its reason in `run_metadata.json` and exits
non-zero when nothing is processed; and `run_full_all.sh` uses
`HNL_TEMPLATE_PYTHON` for the ROOT/Pythia stage when supplied, or consumes
pre-generated templates. A clean checkout without templates therefore fails
rather than emitting an empty plot. The
dead, stale `SEP_MIN/SEP_MAX/P_CUT` block was removed from
`analysis/constants.py` (the live cuts were always sourced from
`decay_reco_acceptance.py`).

## P0: define the detector-level result

### 1. Add detector response to the decay/reconstruction acceptance

**Current code:** the old analytic two-body acceptance (`analysis/sensitivity.py`)
is gone. `analysis/decay_reco_acceptance.py` now runs a per-event Monte Carlo:
for each HNL four-vector that crosses the fiducial volume it samples decay
vertices along the flight path, draws flavor-dependent FairShip rest-frame decay
templates, boosts them to the lab, keeps the two highest-momentum charged stable
daughters (best-two-track), and runs the shared `../higgs/reco_common` bounded
4-hit reconstruction plus the PR #13 selection (gate/pointing/collinearity/
timing). The visible branching fraction is implicit (the fraction of templates
with two reconstructable charged tracks), not a single inclusive `BR_vis`
factor. What is still idealized: detector response is geometric wall hits with
Gaussian position/time smearing only -- no material interactions, tracking or
vertexing inefficiency, trigger/readout, pileup, dead regions, or occupancy --
and the FairShip templates are Pythia-level final states with no detector
simulation.

**Required work:**

- implement detector response beyond geometric hits + smearing: tracking,
  vertexing, particle thresholds, material interactions, trigger/readout, and
  event-selection efficiencies;
- include pileup, timing, dead regions, and occupancy if they affect the
  proposed detector;
- assign per-channel / per-flavor efficiency uncertainties and validate the
  FairShip visible-final-state modelling against an independent generator.

**Completion test:** an efficiency map or detector simulation with versioned
inputs reproduces benchmark samples, and the sensitivity code consumes
per-channel efficiencies with uncertainty variations.

### 2. Validate the geometry and detector configuration

**Current code:** the shared `../higgs/grendel_geometry.py` (the PR #13 single
source) builds a closed tunnel mesh from one survey polyline, a fixed 22 m
vertical position, and a fixed 24 cm wall inset, and classifies tracker vs
scintillator surfaces. Rays originate at `(0, 0, 0)`. The mesh contains no
supports, services, inactive regions, material, alignment uncertainty, or
configurable detector layout.

**Required work:**

- confirm the survey coordinates, CMS IP transform, tunnel cross-section, and
  proposed active volume with the detector/design owners;
- version the geometry and expose layout parameters rather than relying on
  import-time constants;
- model inactive and inaccessible regions and evaluate alignment/survey
  variations;
- document the coordinate convention and independently compare ray
  intersections against a reference geometry implementation.

**Completion test:** signed-off geometry inputs, reference intersection tests,
and a geometry-variation envelope propagated to the curves.

### 3. Supply backgrounds and a statistical model

**Current code:** `analysis/constants.py` defines exclusion as
`N_signal >= 3`, a zero-background 95% CL approximation. Signal now passes the
PR #13 selection (gate/pointing/collinearity/timing), and a cosmic
decay-in-flight background is modelled by the shared `../higgs/` machinery, but
that background is not folded into the limit: there is still no background count
in the statistical model, no control region, nuisance parameter, systematic
uncertainty, or coverage calculation.

**Required work:**

- estimate beam-, collision-, cosmic-, neutrino-, and detector-induced
  backgrounds after the final selection;
- define control samples and background uncertainties;
- include luminosity, reconstruction, trigger, geometry, signal-model, and
  finite-MC nuisance parameters;
- choose and implement the statistical prescription (for example CLs or a
  documented Poisson construction) and validate its expected coverage;
- define how open upper/lower sensitivity boundaries are reported.

**Completion test:** the final likelihood/configuration is versioned, has toy
or asymptotic validation, and produces expected limits with uncertainty bands.

### 4. Fix the physics-hypothesis definition

**Current code:** production and lifetime tables are evaluated for the three
single-flavor patterns `Ue`, `Umu`, and `Utau`. The canonical
`data/published/MANIFEST.json` and diagnostic `run_metadata.json` now state the
single-flavour Majorana hypothesis. The convention itself is no longer open:
both
the HeavyN UFO (`N1` self-conjugate) and HNLCalc (charge-conjugate modes summed
per channel) are Majorana, so the result is self-consistently Majorana, matching
PBC BC6/7/8 and ANUBIS (see the 2026-06-29 resolution under the appended note).
The outstanding work is to complete the independent factor-of-two audit and,
if needed, support non-single-flavour hypotheses -- not to choose or document
the convention again.

**Required work:**

- retain the machine-readable Majorana declaration and independently validate
  the convention used by HNLCalc and the MadGraph UFO;
- audit all factors of two and charge-conjugate channels under that convention;
- support arbitrary mixing ratios if the intended result is not restricted to
  the three single-flavor hypotheses;
- regenerate lifetime, visible branching fraction, production, and decay
  inputs consistently for every advertised hypothesis.

**Completion test:** each output curve carries machine-readable model metadata
and agrees with independent benchmark rates for the same convention.

## P1: production model and its uncertainties

### 5. Produce FONLL scale, PDF, heavy-quark-mass, and grid-coverage variations

**Current code:** the committed bottom and charm grids in
`data/production/fonll/central/` carry one central scale choice and the central
`NNPDF40_nlo_as_01180` member, stopping at `pT = 50 GeV` and `|y| = 3`. The full
variation set has now been generated in the external NNPDF40 workspace -- a
218-grid SHA-256 `variation_manifest.json` (7-point scale, 100 NNPDF4.0 NLO
replicas, `m_b`/`m_c`) plus the `as_01170`/`as_01190` alpha_s companions. The
published per-mass band propagates the 111 coherent central/scale/PDF/heavy-mass
curves; the alpha_s production+analysis companion chains are supported by
post-processing but are not included in the current bundle. The `pT`/rapidity
truncation is still unbounded.

**Progress (2026-07-16): post-beta-fix exact FONLL campaign published.** The 111
coherent variations (6-point scale + 100 NNPDF4.0 replicas + 4 `m_b`/`m_c`)
were rerun end to end after the timing `beta=p/E` correction, using all detector
hits, 50 decay samples, `event_chunk=1000`, and two analysis workers.  All 111
curves retain all 54 anchors, and all 110 non-central members agree with central
on finite/open topology.  `analysis/combine_band.py` now records contributing
member counts and refuses to turn a topology change into a numeric ribbon.

The median combined half-widths are `-0.101/+0.131` dex on the lower edge
(scale dominated) and `-0.0066/+0.0081` dex on finite upper edges.  The largest
upper shift, 0.114 dex, survives independent exact-200 controls.  The previous
order-one upper structure was caused by the 4,000-hit weighted-resampling cap
and is retired.  The 5,994-row raw table, grid/campaign manifests, hashes,
combined band, and numerical controls are tracked in `data/published/bundle/`.
The diagnostic is separate from the central-only paper comparison.

**Required work:**

- correct the charm `mu_F = 0.5` scale prescription, which currently requests
  the NNPDF4.0 PDF below its `QMin = 1.65 GeV` for all `pT < 2.94 GeV` and so
  sets the charm scale envelope's `-76%` lower edge from an LHAPDF
  extrapolation rather than from perturbative uncertainty. This is the largest
  single term in the charm production budget and the highest-value item in this
  section; it is expected to *narrow* the band. See the workspace item 4 below
  for the measured ratios, the two-chain cost, and the accepted-yield check
  that should precede it;
- extend the `pT` and rapidity grids, or place a quantitative bound on the
  omitted contribution after GRENDEL acceptance (still outstanding);
- add a tracked driver for the alpha_s companion production+analysis chains
  (currently the `--alphas-lo/--alphas-hi` curves are run by hand). Measured
  effect is negligible (0.44%/1.49% half-spread); do this for PDF4LHC
  prescription completeness, not for physics -- see workspace item 5 below.

Use the same three-column rectangular table format as the committed grids:

```text
# pT_GeV  rapidity  d2sigma_dpTdy_pb_per_GeV
```

The tracked `data/published/bundle/FONLL_GRID_MANIFEST.json` and
`FONLL_MANIFEST.json` now contain the FONLL revision, PDF member, masses, scales,
beam energy, grid bounds, campaign configuration, and checksums.

**Completion test:** central and varied grids can be regenerated from the
manifest, and the full production plus analysis chain yields an uncertainty
band at every mass point.

### 6. Replace shared heavy-flavor shapes and static fragmentation fractions

**Current code:** `B+`, `B0`, and `Bs` share one bottom shape; `D0`, `D+`, and
`Ds` share one charm shape. Species masses are rebuilt after sampling.
`FRAG_B` and `FRAG_C` are constant fractions measured in restricted and
different acceptances, extrapolated over the full grid without covariance.
Fragmentation-function parameters and feeddown uncertainties are absent.

**Required work:**

- provide species-specific `d2sigma/dpTdy` grids or differential
  fragmentation fractions for each weakly decaying species;
- propagate the experimental fraction covariance and its `pT`/rapidity
  dependence;
- vary the fragmentation-function model and parameters used in the FONLL
  calculation;
- document strong-decay feeddown conventions so no parent is double counted.

**Completion test:** species fractions close consistently at each kinematic
point, reproduce their source measurements in the quoted acceptance, and
produce correlated curve variations.

### 7. Add omitted charmed-baryon HNL sources

**Current code:** `OMITTED_FRAG_C` records approximately 36.7% of the charm
fragmentation allocation as `Lambda_c+`, `Xi_c0`, `Xi_c+`, and `J/psi`, but
none contributes to direct HNL production or to the induced-tau pool.

**Measured impact of the dominant term (Lambda_c, 2026-07-17): negligible,
below MC noise.** Using HNLCalc's `get_3body_dbr_baryon` (the `Lambda_c -> Lambda`
form factors, `dq2dm122` integrator), the charm-sector yield boost from adding
`Lambda_c -> Lambda l N` is `frag_Lc * BR(Lc) / sum_D frag_D * BR(D)`. Folded
with the charm channel's accepted-yield fraction (audit baseline), the lower-edge
shift is **<= 1.2% (<= 0.005 dex)** across the whole charm-relevant range, for
both `Ue` and `Umu`:

| m_N (GeV) | Ue edge shift | Umu edge shift |
|---|---|---|
| 0.4-0.5 | -1.1% | -1.2% (max) |
| 0.7 | -0.6% | -0.5% |
| >= 1.0 | <= -0.03% | ~0 |

This is far below the 4.1% median production-MC noise, so **`Lambda_c` is not a
curve mover** and does not need to be generated. The reason it is so much smaller
than the analogous BC4 `Lambda_b` fix (`-10` to `-14%`): a baryon has **no
2-body leptonic mode**. `D+`/`Ds -> l N` is helicity-enhanced and two-body;
`Lambda_c -> l N` is forbidden (baryon number), leaving only the phase-space-
suppressed semileptonic `Lambda_c -> Lambda l N` (`BR ~ 1%`, confined to
`m_N < m_Lc - m_Lambda - m_l ~ 1.17 GeV`). This **supersedes** the audit's
`missing_charm_baryon_channels` proxy (median 0.1%, max 16.9%), which
overestimated by ~14x by scaling on fragmentation without the semileptonic
suppression. Caveat: the estimate is a production-BR ratio; the softer `Lambda_c`
HNLs would have somewhat lower acceptance, so `<= 1.2%` is an upper bound.

**Required work (residual):**

- `Lambda_c`: **bounded above (<= 1.2%, below noise)** per the measurement above;
  generation is optional. If ever generated, `get_3body_dbr_baryon('4122',
  '3122', lepton)` with the `dq2dm122` integrator supplies the rate, exactly as
  the `Bbaryon` channel already uses `Lambda_b -> Lambda_c l N`;
- `Xi_c0`/`Xi_c+`: even more phase-space suppressed (heavier `Xi` recoil) and a
  smaller fragmentation share than `Lambda_c`, so bounded below the `Lambda_c`
  number by the same argument -- a short explicit check would close them;
- evaluate whether charmonium (`J/psi`) HNL modes are relevant over the
  0.2--10 GeV mass grid and either implement or quantitatively dismiss them;
- audit any additional weakly decaying charm species omitted from the current
  closure.

**Completion test:** the charm fragmentation accounting is explicit and every
omitted component has either a generated channel or a documented negligible
bound on the final curves. `Lambda_c` now has a documented bound; `Xi_c`/`J/psi`
remain to be closed.

### 8. Replace the inclusive bottom-baryon closure approximation

**Current code:** the entire bottom-fragmentation remainder is treated as
`Lambda_b`, with the bottom FONLL shape rebuilt at the `Lambda_b` mass and only
`Lambda_b -> Lambda_c l N` generated. `Xi_b` and `Omega_b` are not separated.

**Required work:**

- provide differential production fractions and spectra for `Lambda_b`,
  `Xi_b`, and `Omega_b`;
- implement their relevant HNL and tau-producing decays with current form
  factors;
- propagate the large baryon-fraction kinematic dependence and covariance.

**Completion test:** bottom species close without an unidentified remainder,
and alternative baryon models are included in the production uncertainty.

### 9. Replace the borrowed Bc model

**Current code:** `Bc` uses the bottom FONLL shape at the `Bc` mass and an
unvalidated inclusive normalization `SIGMA_BC_PB = 0.9e6 pb`. That same model
also normalizes `Bc -> tau nu`.

**Required work:**

- generate a dedicated 14 TeV `Bc+ + Bc-` spectrum using BcVegPy, NRQCD, or a
  documented equivalent;
- provide normalization, scale, PDF, heavy-quark-mass, and NRQCD/model
  uncertainties;
- keep the charge convention explicit;
- propagate the correlated normalization to direct-Bc and induced-tau
  channels.

**Completion test:** a versioned grid and uncertainty set replace both the
shape reuse and `SIGMA_BC_PB`, with independent cross-section benchmarks.

### 10. Replace the kaon production and transport model -- DONE (2026-07-18)

**Status:** the default `generate_kaon_csvs.py` now uses the Pythia 8.315 SoftQCD
spectrum (`SIGMA_KAON_PB = 6.535e11 pb`) plus a charged-kaon transport survival
weight (`d_esc`), and BC6/BC7 low-mass was rerun and republished with it (see the
transport bullet above and `data/published/`). The paragraphs below record the
original (pre-2026-07-18) state and the decision gate that led here.

**Original code (superseded):** `generate_kaon_csvs.py` sampled only charged kaons
using a Tsallis `pT` model and a Gaussian rapidity with `SIGMA_KAON_PB = 3.0e11 pb`
(now the legacy `--spectrum tsallis` fallback). Neutral `K_S/K_L -> pi l N` modes
available in HNLCalc are omitted. Every kaon was decayed immediately, only the HNL
momentum stored, and analysis rayed every HNL from IP5. Parent lifetime, magnetic
bending, material survival, interaction losses, and the displaced HNL production
vertex were therefore absent. This model can dominate the lowest-mass result.

**Decision gate (measured 2026-07-17): transport dominates; do NOT ship a
spectrum-only fix.** A Pythia 8.315 `SoftQCD:inelastic` run at 14 TeV
(`sigma_inel = 78.9 mb`, `<n_K+-> = 8.28` per inelastic event) was compared
against the Tsallis stub through the actual HNL geometry (accepted-yield proxy
`mean(hit * path_len / beta_gamma)` at the long-lifetime edge, `Ue`/`Umu`,
`m_N = 0.2-0.4 GeV`):

- **Spectrum + normalization arm (both prompt-at-IP): the two errors nearly
  cancel.** The stub's normalization is 2.18x too low
  (`SIGMA_KAON_PB` should be ~6.54e11, not 3.0e11), but its Tsallis rapidity
  (`sigma = 2.5`) is too central -- the true Pythia spectrum is more forward, so
  its accepted fraction is only 0.50-0.57 of the stub's. Net yield ratio is
  1.08-1.25, i.e. a lower-edge shift of only **-4% to -10%** (a slight
  strengthening). At that decision point, the published BC6/BC7 was therefore
  accidentally close to the correct *prompt-IP* spectrum result.
- **Transport arm dominates and flips the sign.** Charged kaons are long-lived
  (`ctau = 3.7 m`) and the stub decays them promptly at IP5. Applying a realistic
  survival fraction `S = 0.1-0.7` weakens the lower edge by **+7% to +204%** --
  an order of magnitude larger than the spectrum arm and in the opposite
  direction.

**Consequence for scope.** Shipping the Pythia spectrum alone would move BC6/BC7
slightly *stronger* while silently omitting the larger transport loss that moves
it *weaker* -- a biased-optimistic curve. So the spectrum fix must NOT be shipped
on its own. This is `Ue`/`Umu` low-mass only (`m_N <~ 0.5 GeV`, the K -> l N
window); `Utau` has no kaon channel (`m_tau > m_K`). Two acceptable resolutions:
(1) **paper-final** -- implement charged-kaon transport (magnetic bending +
material survival + displaced HNL origin) together with the Pythia spectrum, or
(2) **interim** -- keep the then-current parametric flux and label BC6/BC7 below
~0.5 GeV explicitly provisional, transport-dominated. Neutral `K_S/K_L` remain a
smaller separate omission (audit additive proxy `<= 19%`). The Pythia SoftQCD driver is committed at
`production/decay_engine/kaon_softqcd.cc` (build with `pythia8-config`), for
whichever path is chosen.

**IMPLEMENTED 2026-07-17 -- Pythia spectrum + transport survival weight.** This is
the paper-final *approach* of option (1), with two deliberate, physics-justified
simplifications relative to that option's full wording: the material survival is a
fixed-`d_esc` *proxy* (not a propagated material map) and the kaon's magnetic
bending is not modelled (both second-order given the survival-weight equivalence;
see Remaining). `generate_kaon_csvs.py` now defaults to:
- **Pythia SoftQCD spectrum + normalization.** The committed
  `production/data/kaon_softqcd_spectrum.npz` (a 200x160 `(pT, y)` histogram from
  `kaon_softqcd.cc`, reproducible from tracked sources via
  `production/decay_engine/make_kaon_spectrum.py` at the pinned seed 42 with the
  vendored Pythia 8.315) replaces the Tsallis stub, and `SIGMA_KAON_PB = 6.535e11`
  (`sigma_inel = 78.93 mb` x `<n_K+-> = 8.28`) replaces the `3.0e11` stub. The
  legacy path stays behind `--spectrum tsallis`.
- **Charged-kaon transport as a survival weight.** Each kaon carries
  `w = 1 - exp(-d_esc / (beta*gamma * ctau_K))`, the probability it decays before
  being absorbed in dense material. The HNL is then cast from the IP rather than
  from the displaced kaon-decay point (median ~0.7 m, `<= d_esc ~ 1.5 m`): a good
  approximation because that shift is small and nearly collinear with the
  forward-boosted HNL compared with the ~22 m flight to the fiducial volume.
  `production/decay_engine/transport_control.py` casts each HNL from both origins
  against the real fiducial mesh (`production/data/transport_control.json`) and
  finds a displaced/IP accepted-yield ratio of `1.00` within `~3%` on the
  long-lifetime plateau that sets the sensitivity (`ctau_N >= 10 m`: 1.03, 1.00,
  0.99, 0.99), rising to `~1.1-1.2` only in the negligible short-lifetime tail --
  so no per-origin acceptance change is applied (small where it matters, not zero).
  `KAON_D_ESC = 1.5 m` is a proxy for the CMS material budget (calorimeter front
  ~1.3 m); `KAON_D_ESC_RANGE = [1, 3] m` is propagated in the published
  `data/published/bundle/kaon_desc_band.csv`. `--no-transport` restores the
  prompt-at-IP behaviour.
  `run_all.py` uses the new model by default.

**Pre-rerun impact estimate for BC6/BC7 (Ue/Umu, m_N < 0.5 GeV):** the transport
survival suppresses the accepted kaon yield to `~0.3-0.55` of the stub
(`d_esc = 1.5-2 m`), and the Pythia spectrum modestly raises it (net `x1.1-1.25`);
the product was expected to **weaken the low-mass lower edge by roughly +30% to
+80%**. The subsequent published rerun measured the larger x1.6-1.9 weakening
reported below. `Utau` is unaffected (`m_tau > m_K`).

**Remaining:** (1) DONE 2026-07-18 -- BC6/BC7 low-mass republished with the
Pythia 8.315 + transport kaon model (`d_esc=1.5 m`): the low-mass Ue/Umu edge
weakened by x1.6-1.9 and the Umu peak relocated out of the kaon window
(`8.40e-9 @0.365 -> 1.30e-8 @1.5 GeV`). See `data/published/MANIFEST.json`. The
`bundle/` FONLL band + channel/decay diagnostics were re-derived at the one
affected band mass (`0.305 GeV`, Ue/Umu) with the new kaon model (111-variation
FONLL campaign), so the bundle is consistent with the republished central.
(2) pin `d_esc` against a real CMS material map -- the dominant kaon-sector
uncertainty; the `KAON_D_ESC_RANGE = [1, 3] m` band is now published as
`data/published/bundle/kaon_desc_band.csv` (~+-20-25% on the low-mass edge), but
d_esc itself remains a proxy pending the material map;
(3) neutral `K_S/K_L -> pi l N` remain omitted (audit additive proxy `<= 19%`);
(4) magnetic bending of the kaon trajectory is not modelled (second-order given
the survival-weight equivalence, but a soft-kaon check is worthwhile).

**Required work:**

- produce charged- and neutral-kaon spectra from measured data and/or tuned
  Pythia 8 over the full kinematic region that can reach GRENDEL;
- state whether the input is a per-event yield or a cross section and provide
  the corresponding inelastic normalization;
- vary generator tune, soft-QCD model, multiplicity, charge ratio, and
  extrapolation outside measured coverage;
- propagate each kaon through the magnetic field and detector material,
  sample its decay or interaction point, and originate the HNL there;
- include `K_S` and `K_L` lifetimes, regeneration/interactions, and branching
  channels where relevant;
- include secondary kaon production/absorption only if it belongs in the
  signal definition, without double counting generator particles.

**Completion test:** data/generator comparisons, transport validation, and an
uncertainty envelope replace the prompt charged-kaon approximation in
kinematics, origin, and normalization.

### 11. Separate and upgrade electroweak production

**Current code:** LO MadGraph samples use one NNPDF central member and a flat
`K = 1.3`. Direct `W -> l N` and `Z -> nu N` rows are merged without a parent
tag, so all rows receive the W factor. Prompt-tau W and Drell-Yan origins are
separable but use the same central factor. Because the W and Drell-Yan factors
are currently identical (`K = 1.3`), the missing per-row W/Z split has **no
numerical effect on the present curves**; it is a logic refinement that only
matters once the factors differ.

**Required work:**

- carry the LHE mother/process identity into direct-HNL CSV rows;
- calculate process-, mass-, and preferably kinematics-dependent NLO QCD and
  electroweak corrections for W, Z, and Drell-Yan tau production;
- propagate MG scale/PDF variations and matching/shower uncertainties if a
  showered sample is used;
- check process definitions, interference, off-shell contributions, and
  overlap with heavy-flavor or top-origin W samples;
- implement or quantitatively dismiss Higgs-mediated `h -> nu N` production
  for every advertised model hypothesis.

**Completion test:** W, Z, and Drell-Yan yields can be varied independently
and reproduce reference inclusive and differential cross sections.

### 12. Close the heavy-hadron production-channel inventory

**Current code:** direct meson production is a hand-written list of exclusive
two- and three-body modes. It is not checked against the larger HNLCalc channel
inventory; for example, HNLCalc contains `D+ -> omega l N` while the production
list omits it. Higher-resonance and nonresonant semileptonic modes are not
accounted for. Weak cascades such as `B -> D/Ds X` followed by charm-hadron HNL
production are also absent.

**Required work:**

- generate a channel-closure table comparing implemented modes with all
  kinematically allowed modes in the selected rate model;
- add material missing exclusive, higher-resonance, and nonresonant channels,
  or bound their effect;
- include weak heavy-hadron cascades with decay vertices and avoid overlap
  with prompt heavy-flavor production;
- perform the same closure audit for tau-producing parent decays beyond the
  selected 11 sources.

**Completion test:** every allowed parent fraction and decay class is
implemented or assigned a documented negligible bound, with cascade
double-counting tests.

## P1: decay, lifetime, and polarization model

### 13. Replace or pin the HNLCalc rate model

**Current code:** `vendored/HNLCalc` is now pinned to our fork
`leoredi/HNLCalc` `physics-fixes` @ `e292cef9` (upstream `laroccod/HNLCalc`
`main` @ `07f84728` plus four local-fix commits; see `vendored/PROVENANCE.md`).
Direct meson and baryon channels use its embedded form factors and
absolute branching ratios. Several parameterizations are old or weakly sourced
-- e.g. several channels reuse the pion form factor as a proxy, and the CKM
elements and decay constants are hardcoded without uncertainties. These are
HNLCalc's own community-standard choices, inherited by citing the tool rather
than independently revalidated; for an exclusion *contour* that is accepted
practice and sub-leading (a P1 normalization effect on a log-scale reach
boundary -- see the curve-impact audit -- not a P0 result-definition gap), so
they are flagged here for a precision upgrade, not as bespoke errors. Manual
`m_N -> 0` spot checks found substantial differences between selected HNLCalc
integrals and external tau-mode branching fractions, but the comparison script
and benchmark table are not yet committed. The exception is the four local fixes
the pinned fork carries over upstream (including a misplaced Källén-factor
parenthesis in `tau -> P N` and an `|U|^4` mixing factor in the baryon channel;
see `vendored/PROVENANCE.md`): these correct genuine upstream bugs, so the
vendored copy is a deliberate deviation from the citable release until the fixes
are upstreamed.

**Required work:**

- upstream the four local fixes, and either track a maintained HNLCalc release
  or replace it with a local rate library (the fork is pinned but unmaintained);
- update decay constants, CKM inputs, hadron masses, and semileptonic form
  factors to a coherent current set;
- include form-factor covariance and normalization uncertainties for every
  material direct channel, including baryons and Bc;
- validate integrated rates and differential shapes against independent
  calculations and SM-limit branching fractions for all three flavors;
- validate and assign uncertainty to the abrupt exclusive-hadron versus
  inclusive-parton width matching around 1 GeV, including continuity and
  quark-hadron-duality assumptions *(partially done: item 15's seam-derived
  width band `delta(m)` propagates this duality for the lifetime / upper edge;
  the lower-edge form-factor + absolute visible-BR normalization remains open)*;
- make the direct-HNL and induced-tau channels use compatible rate inputs.

**Completion test:** channel-by-channel benchmark tables and regression tests
cover integrated rates and selected differential distributions across mass.

### 14. Propagate parent-to-tau branching-fraction uncertainties

**Current code:** the 11 induced-tau sources use manually reviewed central
values for `Ds`, `D+`, `B+`, `B`, `Bs`, `Bc`, and `Lambda_b` decays. The review
is not encoded as a reproducible source/version table, and no
experimental/theory uncertainty or covariance is propagated. The `Bs` and
`Bc` entries are theory-derived, and the `Lambda_b` measurement contains an
external-normalization uncertainty.

**Required work:**

- store source, date/version, uncertainty components, and correlations beside
  every central value;
- propagate HFLAV/PDG experimental averages and lattice/theory uncertainties;
- correlate common inputs such as `Vcb`, decay constants, normalization modes,
  and form factors;
- define an update procedure rather than manually changing constants.

**Completion test:** a machine-readable input table generates the constants
and produces correlated induced-tau normalization variations.

### 15. Validate the FairShip lifetime and visible-fraction model

**Current code:** the HNL lifetime `ctau(U^2 = 1)` and the visible final states
now come from the FairShip decay templates (`analysis/generate_decay_templates.py`):
the lifetime is `HNLbranchings.computeNLifetime`, and the visible fraction is
implicit in the per-flavor template multiplicity. The previous
`data/ctau/*.dat` lifetime / `BR_vis` tables and the `BR_vis = 1` fallback are
gone. Templates are generated on a fixed mass grid and matched per mass label
(no interpolation between masses).

**Progress (2026-06-25):**

- *Majorana/Dirac factor-2 gate — verified.* The Majorana x2 lives once, in the
  decay (`NDecayWidth` x2/channel -> `ctau`; visible BRs split /2 into
  CP-conjugate final states); production carries only the particle+antiparticle
  charge factor (FONLL `2*sigma`, MadGraph `SM_HeavyN` N1 self-conjugate summing
  both W charges). No double-count on the `U^2 * P_decay` product.
- *Lifetime-uncertainty band — done (upper edge).* The HNL total-width / lifetime
  quark-hadron duality is propagated as a seam-derived band `delta(m)`
  (`analysis/width_band.py`: floor 5% / cap 20% / per-flavor, read off FairShip's
  own `max(meson, quark)` disagreement with a transition envelope so the
  crossover does not notch the band). Driven coherently through both the lifetime
  leg (`ctau`) and the composition leg (`vis_frac`) in
  `analysis/decay_model_band.py`; the composition leg self-cancels to ~1%
  (`Gamma_had` is in both `Gamma_vis` and `Gamma_tot`), so the band is
  lifetime-dominated (~0.08 dex on the dome / upper edge, where it co-sets the
  closure mass).
- *Still open — absolute visible-BR normalization.* A SEPARATE in-scope
  decay-side nuisance (~+/-10%; audit `visible_branching_fraction` = up to
  0.13 dex on the LOWER edge, from the HNLCalc/FairShip decay BRs + form factors)
  is NOT covered by the duality band above. It groups with the lower-edge
  normalizations (item 13 form factors, item 9 Bc). Also open: pin the FairShip
  revision + manifest, the GRENDEL visible definition, finite-template statistics.

**Known result feature (document for referees, 2026-06-25; literature
cross-check 2026-07-17):** the `Umu` limit degrades sharply below
`m_N ~= 0.25 GeV` (e.g. `u2_min` at 0.2 GeV is ~10x weaker for `Umu` than
`Ue`). This is **correct physics, not a bug**: GRENDEL acceptance requires
`>= 2` charged stable daughters (`analysis/decay_reco_acceptance.py`); the
dominant visible 2-track mode `N -> mu+- pi-+` closes below
`m_mu + m_pi = 0.245 GeV`, so the `Umu` visible fraction collapses. Soft
3-body (`N -> mu e nu`, NC `N -> nu e e`) remain at the few-percent level;
`N -> nu pi0 (-> gamma gamma)` dominates below threshold but has **zero**
charged tracks and is invisible to GRENDEL. `N -> e+- pi-+` keeps `Ue` open
down to ~0.14 GeV. Expect a visible step in the low-mass `Umu` (and
analogously `Utau`) contour at the corresponding meson threshold; state it
explicitly so it is not read as an artifact.

**Template measurement (FairShip caches, 2026-07-17):** rest-frame fraction
of templates with `>= 2` charged stable daughters:

| m_N (GeV) | Ue | Umu |
|---|---|---|
| 0.200 | ~61% | ~4% |
| 0.245 | ~63% | ~5% |
| 0.260 | ~63% | ~30% |

`Umu` jumps ~6x when `N -> mu pi` opens; published `u2_min` moves
`1.16e-7 -> 2.17e-8` over the same step. Production is kaon-dominated on
both flavors in this window (~100% `Kmeson` at 0.305 GeV) -- the cliff is
**decay visibility**, not production.

**Why HNLimits / community U_mu plots do not show this cliff:**

1. **Existing low-mass exclusions are mostly peak searches, not displaced
   LLP.** E949 (~0.18-0.30 GeV) and NA62 (~0.20-0.38 GeV) constrain
   `K+ -> mu+ N` via missing mass and assume the HNL escapes undecayed.
   They never require `N -> mu pi`, so they cannot exhibit a 0.245 GeV
   visibility cliff. Those fill the Hostert grey band below ~0.4 GeV
   ([mhostert/Heavy-Neutrino-Limits](https://github.com/mhostert/Heavy-Neutrino-Limits)).

2. **SHiP's published sub-kaon reach omits kaons.** arXiv:1811.00930:
   for `M_N <~ 500 MeV` kaons dominate production, but most stop in the
   target/hadron stopper, so SHiP includes **only charm and beauty**; the
   Fig. 5 dashed extension below the kaon mass is **D-only**. SHiP also
   requires `>= 2` charged tracks and would cliff if it claimed
   kaon-sourced `Umu` there -- it simply does not claim that shelf.

3. **GRENDEL does include the kaon channel** (Pythia 8.315 SoftQCD spectrum plus
   the charged-kaon transport proxy; see item 10). That populates the kaon
   window and therefore *exposes*
   the `Umu` 2-track cliff. Comparing GRENDEL's displaced contour to the
   Hostert composite below 0.25 GeV is apples-to-oranges.

Separate from this cliff: absolute strength of the whole BC6/BC7 shelf
below ~0.5 GeV is still limited by **kaon transport** (item 10), which is
a modeling uncertainty, not the reason community plots lack the 0.245 GeV
step.

**Required work:**

- pin the FairShip module revision and write a manifest (mass grid,
  `n_templates`, seed, couplings, decay-selection config);
- define exactly which modes count as visible for GRENDEL and confirm the
  template-derived visible fraction matches that definition;
- propagate lifetime and visible-fraction uncertainties and correlate them with
  production;
- quantify finite-template statistics and confirm the mass-grid spacing needs
  no interpolation (or add it);
- cross-check the FairShip total width / lifetime and visible fraction against
  an independent calculation.

**Completion test:** a pinned, manifested template generation, strict loading,
threshold-aware mass coverage, and regression benchmarks against an independent
lifetime/visible-fraction calculation.

### 16. Complete tau decay modes and implement realistic spin correlations

**Current code:** semileptonic energy and invariant-mass distributions are
HNLCalc-weighted, but remaining angles are sampled isotropically. Two-body
tau-to-HNL decays use fixed unit-analyzing-power longitudinal asymmetries
with the origin-dependent sign (`+1` for W-origin taus, `-1` for the
helicity-suppressed heavy-meson leptonic sources; see
`production/decay_engine/tau_decay.py`). The magnitude is still an
approximation: the true analyzing power is mass-dependent and below one for
the vector modes, the asymmetry is applied about the lab tau direction
(no Wigner rotation, which matters most for `Ds -> tau nu` where the tau is
nearly at rest in the parent frame), and tau charge is not retained in the
prompt-tau CSV. Semitauonic, baryonic, and Drell-Yan tau sources are treated
as unpolarized, although their polarization depends on charge, phase space,
and production kinematics. Tau-to-HNL decay generation includes only `pi`,
`K`, `rho`, `K*`, and leptonic three-body modes; multi-hadron spectral modes
such as `a1/3pi` are absent.

**Partial bound on the analyzing-power error (2026-07-17): a candidate BC8
(Utau) mover, not dismissible.** The unit-analyzing-power approximation is
**exact** for the pseudoscalar 2-body modes (`pi`, `K`): a spin-0 daughter
carries the full tau polarization, so `asym = +-1` is correct there. The error
lives only in the **vector** modes (`rho`, `K*`), whose true analyzing power is
`alpha_V = (m_tau^2 - 2 m_V^2)/(m_tau^2 + 2 m_V^2) = 0.45` (`rho`), `0.33`
(`K*`), not 1. From `compute_tau_production_br_components`, the vector modes are
**33-39% of total tau -> N production for `m_N < 1 GeV`** (falling to ~0 above
1.2 GeV as the 2-body vector channels close). That gives an upper bound on the
fractional yield perturbation of `f_vector x (1 - alpha_V) ~ 18-22%` at low
`m_N`, i.e. an edge shift bound of very roughly `<= 10%` -- **above the 4.1% MC
noise floor**, and consistent with the audit `production_spin_and_residual_angles`
(median 7.1%, max 16.8%).

Unlike `Lambda_c` (item 7, dismissed below noise), tau spin **cannot be
dismissed from the branching ratios alone**. The bound is an overestimate,
because the tau boost from the `D`/`B` parent partly washes the rest-frame
angle out (the asymmetry mostly shifts the HNL *energy* spectrum, hence the
decay length and fiducial fraction, rather than its direction). Pinning the
actual BC8 curve shift needs the controlled acceptance measurement below --
generating the induced-tau HNL 4-vectors under (A) the current `asym = +-1` and
(B) per-mode analyzing power (`+-1` for `pi`/`K`, `alpha_V x sign` for `rho`/`K*`)
and re-scanning. This is Utau-only; BC6/BC7 (`Ue`/`Umu`) are unaffected.

**Required work:**

- use full matrix-element event generation or helicity amplitudes for the
  relevant parent-to-tau-to-HNL and direct-HNL chains;
- retain tau spin density matrices through production and decay;
- include Drell-Yan and W tau polarization with their kinematic dependence;
- add all material multi-hadron tau decay modes using validated spectral
  functions or a spin-aware generator such as TAUOLA;
- demonstrate branching-fraction closure as a function of HNL mass;
- validate angular and energy distributions against an independent generator.

**Completion test:** spin-correlated benchmark distributions and a quantified
curve shift relative to the current approximation.

## P1/P2: numerical and workflow uncertainty

### 17. Quantify finite Monte Carlo and integration errors

**Current code:** most sampled pools default to 100,000 rows, HNLCalc
three-body integrals commonly use 500 samples, and rejection-sampling ceilings
are numerical approximations (the weighted samplers now detect a proposal
exceeding the ceiling and restart unbiased with a raised ceiling, so an
under-estimated grid scan can no longer silently truncate the density, but
grid resolution still sets the proposal efficiency). The curve contains no
finite-sample uncertainty.

**Required work:**

- run independent seeds and increasing pool sizes by channel and mass;
- propagate weighted-event variance through geometry and the U-squared scan;
- validate rejection envelopes and retry limits near thresholds;
- converge the HNLCalc integrations and differential samplers;
- allocate events adaptively to rare sources and low-acceptance regions.

**Completion test:** per-mass MC/integration uncertainty is below a declared
target or appears in the final band.

### 18. Converge the acceptance, serialization, and exclusion scans

**Current code:** the accepted decay probability is a Monte Carlo estimate over
`DECAY_SAMPLES = 100` decay vertices sampled uniformly along the in-volume path,
with the reconstruction + selection re-evaluated per vertex
(`analysis/decay_reco_acceptance.py::build_event_mc`/`scan_u2`). The mixing scan
has 200 points from `1e-12` to `1e-1`, followed by local interpolation. The mass
grid has a minimum spacing of 15 MeV and becomes much coarser at high mass;
plotted lines simply connect calculated points. `run_sensitivity.py` now exposes
`--decay-samples`, `--max-hit-events`, and `--mass-stride` (2026-06-23) so
`DECAY_SAMPLES` and the hit-event count can be varied for convergence/approximate
scans without code edits.

**Progress (2026-07-16): exact-hit convergence controls recorded.** Exact-50
versus the published exact-100 central anchors differs by at most 0.0444 dex on
the lower edge and 0.0331 dex on finite upper edges.  Chunked versus unchunked
exact evaluation differs by at most 0.00537/0.00103 dex (lower/upper), and the
tested one-worker/two-worker point is bit-for-bit identical.  Independent
exact-200 repeats validate the largest FONLL upper shifts.  Three independent
exact-400 Bc endpoint repeats agree on island topology; their finite boundary
spread is at most 4.5%.  Machine-readable results and hashes are in
`data/published/bundle/NUMERICAL_CONTROLS.json`.

**Required work:**

- increase or replace the eight-significant-digit CSV serialization and
  quantify its effect on reconstructed mass, boost, geometry, and acceptance,
  especially for ultra-boosted low-mass HNLs;
- demonstrate the decay-vertex Monte Carlo is converged (vary `DECAY_SAMPLES`
  and seeds) against an exact quadrature of the decay-in-volume integral;
- replace or validate the fixed mixing scan with bracketed root finding for
  both exclusion boundaries;
- demonstrate stability under mass-grid refinement, especially at production
  thresholds and where the dominant channel changes;
- do not claim sub-grid precision: the current grid cannot support estimates
  every 0.3 MeV.

**Completion test:** curve boundaries are stable under all three refinements
within a declared numerical tolerance.

### 19. Make caches and channel completeness self-validating

**Current code:** geometry caches are named only by flavor and mass and do not
yet embed a checksum of the input vectors, geometry version, origin, or cuts.
Since 2026-06-16 a cache is invalidated when its source 4-vector CSV is newer
than the cached NPZ (mtime guard), so regenerating a combined CSV forces a fresh
ray-cast; a full content/geometry checksum is still outstanding (use
`--force-geometry` after geometry-code edits). `run_sensitivity.py` now records
every skipped `(flavor, mass)` with its reason in `run_metadata.json` and exits
non-zero when no point is processed, so dropped masses are no longer silent.
Channel combination is strict by default everywhere since 2026-06-11:
`combine_channels.py` fails on missing or malformed channels while accepting
zero-byte closed-channel sentinels, partial runs must opt in via
`--allow-missing`, and `run_full_all.sh` therefore fails on an incomplete run.

**Negative optimization result (2026-06-25):** do **not** reuse the central
monolithic geometry cache row-for-row across FONLL variation runs. The audit in
`geometry_reuse_study/` compared exact geometry caches for two scale variations
against the central cache and found large hit-mask disagreement:
`xor_hits / variation_hits ~= 0.90` median for the full combined samples. The
regenerated FONLL channels (`Bmeson`, `Dmeson`, `Bbaryon`) are the failure mode,
with median mismatch around `1.8--2.0` per variation hit; row `i` in a varied
grid is not the same trajectory as row `i` in the central grid. Blind cache
copying would therefore invent and miss detector crossings. The hardlinked
FONLL-independent channels checked in the audit (`Bc`, `WZ`) matched exactly,
so a future speedup may be channel-aware: reuse central geometry only for
hardlinked frozen channels (`Bc`, `Kmeson`, `tau`, `WZ`) and ray-cast
regenerated channels normally. That is a useful geometry-stage optimization,
not a replacement for exact per-variation geometry.

**Required work:**

- key caches by input checksum and geometry/analysis configuration;
- if geometry reuse is introduced, make it channel-aware and prove row identity
  for each reused channel;
- write a run manifest with code revision, environment, seeds, input
  checksums, channel row counts, and summed weights;
- verify that no stale channel or cache from another run tag is consumed.

**Completion test:** deliberately modifying an input invalidates the cache,
and incomplete runs fail before producing a combined curve.

### 20. Pin the software and external tools

**Current code:** the Python environment, MadGraph, LHAPDF, PDF data,
compiler/runtime, and some vendored sources are not captured in one immutable
run description. (`vendored/HNLCalc` is now pinned to `leoredi/HNLCalc@e292cef9`;
the broader environment lock remains outstanding.)

**Required work:**

- commit an environment lock and record exact external executable versions;
- checksum vendored tables, UFO files, and generated process cards;
- preserve MG run cards, param cards, logs, and random seeds for final samples;
- add a one-command clean-room reproduction of representative production and
  analysis points.

**Completion test:** a fresh checkout can reproduce benchmark weights and
curve points within the stated stochastic tolerance.

### 21. Expand regression coverage for the new production paths

**Current code:** `tests/test_production_channels.py` (2026-06-11) covers the
11-source induced-tau table (source count, polarization tags, kinematic
thresholds), the induced-tau and `Bbaryon` drivers through their production
interfaces (on-shell HNLs, weight normalization, closed-mass sentinels, the
`dq2dm122` sampler), and `tests/test_combine_channels.py` distinguishes
missing channels from zero-byte sentinels under the strict default. The
stale future-dated test metadata has been removed.

**Required work:**

- tighten the 30% normalization tolerances once item 13 provides
  deterministic benchmark rates (HNLCalc's integrator is Monte Carlo at
  `nsample=500`, so the current checks pin factor-2/|U|^4-class bugs only).

### 22. Propagate all variations into final plots and tables

**Current code:** `plot_exclusion.py` draws the central connected curves with
honest open-boundary handling (fills to the axis edge, no false closing line,
marks the open direction, never bridges insensitive masses). The GRENDEL "money
plot" (`analysis/plot_money.py`) builds on it and overlays the in-scope
production+decay bands: FONLL theory (item 5), direct-Bc normalization
(`bc_nuisance.py`, item 9), and the decay-model width band (`decay_model_band.py`,
item 15). Detector, background, numerical, and luminosity bands are still absent
(P0 items 1-3). External constraints are not overlaid by this package-local
diagnostic; the central publication comparison in `shared/curves_PBC` does
include versioned existing bounds and proposal curves.

**Progress (2026-07-16): topology-safe diagnostic deliverable.** `analysis/plot_money.py`
produces the single-flavor (m_N, |U|^2) projection with: the FONLL/Bc/decay-model
bands above; a `band_registry`-driven combination (`combine_band.py`); a
real-point high-mass closure; a
hypothesis/scope `run_metadata.json` (Majorana, single-flavor, N>=3, idealized
partner handoff); and a machine-readable bundle (central + per-band CSVs). It is
reproducible from the tracked `data/published/bundle/` on a clean clone.  Ribbon
interpolation is limited to contiguous finite anchor segments and stops before
nuisance-induced topology changes.  The publication comparison figure remains
central-only; these are theory/model diagnostics, not confidence bands. Still
open: detector/background/luminosity/numerical bands, a formal correlated
combination across all axes, and external-constraint overlays.

**Required work:**

- add the detector, background, numerical, and luminosity bands once those
  inputs exist (P0 items 1-3, items 17-18);
- include external constraints with source/version tracking if required for the
  intended figure.

**Completion test:** every visible band can be traced to stored variation
outputs and regenerated from the run manifest.

## Work that can be done in the NNPDF40 workspace

The active FONLL grid-generation workspace is external to this repository:
`/Volumes/sandbox/projects/aaaPHYSICSaaa/shared/NNPDF40/fonll-local`. Do not maintain
a second active copy under this HNL tree; use the external workspace for item 5
and part of item 6. Status as of 2026-06-23:

1. **Done.** `scripts/generate_meson_grids.py --campaign` accepts
   `(mu_R, mu_F)` (verified against FONLL's `read ffact,fren` order), PDF
   members, and `m_b`/`m_c` variations, with per-grid tags, headers, and a
   SHA-256 manifest. `validate_variation_gate.py` is the hard-abort
   central-reproduction gate.
2. **Done.** The full campaign completed: `output/variation_manifest.json`
   indexes 218 coherent grids (100 NNPDF4.0 replicas per quark + 7-point scale
   + `m_b`/`m_c`), and the regenerated centrals reproduce the committed grids.
3. **Done with caveat.** `combine_variations.py` writes scale/PDF/mass/
   combined envelopes plus a manifest. The PDF statistics use replica
   members >= 1 only (member 0 excluded; an earlier revision included it and
   biased the band low). Envelope grids are pointwise constructions marked
   `envelope_band` and are refused by `production/fonll/fonll_parser.py`;
   only coherent individual grids may be sampled.
4. **Known artifact (restated 2026-07-17; the earlier text misattributed it).**
   Both charm `mu_F = 0.5` scale points request the proton PDF below the
   NNPDF4.0 grid minimum at low `pT`. The resulting suppression is an LHAPDF
   extrapolation, not a scale uncertainty. `NNPDF40_nlo_as_01180` has
   `QMin = 1.65 GeV`; the charm central scale is `mu = sqrt(m_c^2 + pT^2)` with
   `m_c = 1.5 GeV`, so `mu_F = 0.5 mu` starts at 0.75 GeV at `pT = 0` and stays
   below `QMin` for all `pT < 2.94 GeV`. Bottom is structurally immune:
   `mu_F = 0.5 sqrt(m_b^2 + pT^2) >= 2.375 GeV` always exceeds `QMin`, which is
   why its scale spread is well behaved.

   An earlier revision of this item named only the `(0.5, 0.5)` point and
   quoted "suppression x14-20 for pT <~ 2 GeV". That figure is *accurate* for
   the point it names (at y=0 it is x14 at pT=0.51 and x20 at pT=1.01); the
   defect is that it is **incomplete**. Two points are affected, and the
   omitted one -- `(1, 0.5)` -- is roughly twice as severe and is the one that
   actually sets the band edge. Ratios to central, integrated over rapidity:

   | pT (GeV) | `muR1_muF0p5` | `muR0p5_muF0p5` |
   |---|---|---|
   | 0.51 | 0.035 | 0.076 |
   | 1.01 | 0.023 | 0.051 |
   | 2.02 | 0.093 | 0.220 |
   | 4.04 | 0.440 | 0.819 |

   The suppression switches off at the `QMin` crossing, as an extrapolation
   artifact should. This matters because `(1, 0.5)` sets the *lower edge* of
   the charm scale envelope: its trapezoid-integrated ratio is 0.236, i.e.
   `-76%`, the single largest term in the charm production budget (PDF 1sigma
   is 14.8%; alpha_s 0.44%, item 5). The published charm scale band's lower
   edge is therefore likely too wide, and a correct prescription should
   *narrow* the band rather than widen it. See
   `fonll-local/PUBLICATION_BASELINE.md` for the prescription caveat before
   quoting the low-mass charm scale band.

   **This is a band-edge problem, not a central-curve problem.** The central
   charm grid also dips below `QMin` -- `mu_F = sqrt(m_c^2 + pT^2) < 1.65 GeV`
   for `pT < 0.687 GeV` -- but at the 0.505 GeV node spacing that is a single
   populated node carrying 2.8% of the charm trapezoid integral, below the 4.1%
   median production-MC noise measured in `audits/curve_impact_20260610/`. The
   published central exclusion curve is therefore not materially affected; only
   the scale band's lower edge is.

   **Required work:** regenerate the two charm `mu_F = 0.5` grids under a
   defensible prescription (freeze `mu_F` at `QMin`, or restrict the variation
   range and document the restriction), then rerun those two coherent chains --
   the same two-chain cost as the alpha_s companions in item 5, for a far
   larger effect. Before committing to it, compute the GRENDEL-accepted charm
   yield below `pT ~ 3 GeV`: the artifact is confined to that region, so the
   accepted fraction determines whether the contour moves or only the grid
   does. 60% of the charm trapezoid integral lies below `pT = 3 GeV`, but that
   is a production-level figure with no acceptance folded in, and it is not a
   substitute for the accepted-yield calculation.

   This item is charm-only, and so applies to the HNL benchmark alone: the BC4
   scalar and BC10 fermiophilic-ALP benchmarks sample bottom exclusively and
   are immune by construction.
5. **Alpha_s grids and combiner done; propagated companion runs outstanding.** The alpha_s companion
   grids are generated -- `as_01170`/`as_01190` central grids for both quarks --
   and `combine_alphas.py` writes the PDF4LHC alpha_s envelope to
   `output/envelopes/alphas_manifest.json` (diagnostic-only, refused by the
   sampler; only the three coherent central grids are sampled). On the HNL side
   tracked code can fold an alpha_s term:
   `analysis/combine_band.py --alphas-lo/--alphas-hi` adds the PDF4LHC
   half-difference of supplied `as_01170`/`as_01190` companion curves in
   quadrature. The current published bundle does not supply or include those
   curves. Outstanding: a tracked driver for the companion production+analysis chains
   (currently run by hand via the `run_variation_band` pattern); extended
   `pT`/rapidity coverage to bound the GRENDEL-accepted tail (item 6); and
   per-species fragmentation outputs/variations (item 8).

   **Expected effect (measured 2026-07-17): negligible.** The PDF4LHC alpha_s
   half-spread on the trapezoid-integrated grids is 0.44% (charm) and 1.49%
   (bottom). Folded in quadrature with the terms it joins -- bottom scale
   `+-42%`, charm scale `-76%/+112%`, PDF 1sigma 14.8% (charm) / 3.9% (bottom)
   -- it moves the bottom band by `sqrt(42^2 + 1.49^2) - 42 = 0.03` percentage
   points, below the campaign's own MC noise and invisible on any plot. Run it
   for PDF4LHC prescription completeness and provenance, not for physics; it
   should not be scheduled ahead of the charm `mu_F` artifact in item 4, which
   costs the same two chains and is ~100x larger.

   **Mechanical blocker (why it is still "by hand").** Two things, both small.
   `run_variation_band.py:67` hardcodes
   `GRID_STEM = "..._as_01180_..."`, and `_env_for` (line 178) rebuilds each
   grid path from that stem rather than from the manifest's `path` field -- the
   alpha_s companions differ *in the stem*, not the variation tag, so they are
   unreachable. And `variation_manifest.json` indexes only `as_01180` grids
   (2 central + 12 scale + 200 pdf + 4 mass = 218), so `discover_variations`
   cannot see the companions or their checksums. Fixing means parametrizing the
   stem per entry and synthesizing two coherent variation dicts; alpha_s is
   coherent across both quarks, so it is structurally a scale point.

A full replica campaign takes days and produces large logs unless run with
`--compress-logs`; it is resumable via `--reuse-existing-grids`.

## Recommended execution order

1. Define the detector signal, geometry, backgrounds, and statistical model
   (items 1--4).
2. Pin/regenerate HNL rates, lifetime, and visible fractions (items 13--15).
3. Produce FONLL uncertainty grids and fix the largest production models:
   kaons, charm baryons, Bc, fragmentation, and channel closure (items 5--12).
4. Upgrade HNL/tau decay completeness and spin correlations (items 13--16).
5. Complete numerical convergence, strict provenance, regression coverage,
   and final uncertainty propagation (items 17--22).

Until the P0 items are resolved, curve-shift estimates for the remaining items
are diagnostics of the current simplified analysis, not uncertainties on a
fully defined experimental exclusion.

## Appended note: production-definition mismatch in comparison plots

The current published GRENDEL HNL curve is inclusive in production: it combines
mesons, kaons, baryons, induced taus, prompt taus, and direct electroweak `WZ`
samples. External proposal curves are not necessarily inclusive in the same
way, so overlays against ANUBIS, FASER/FASER2, SHiP, CODEX-b, or similar
proposals can become a result-definition problem rather than a simple
sensitivity comparison.

The immediate ANUBIS issue is important. The current GRENDEL `WZ` channel
contains explicit MadGraph decay-chain samples `pp -> W -> ell N` and
`pp -> Z -> nu N`, merged under one channel. The newer SET-ANUBIS HNL study
uses direct electroweak production through charged-current Drell-Yan,
neutral-current Drell-Yan, and `W gamma` fusion, and leaves hadronic HNL
production for future work. Therefore the current GRENDEL `WZ` sample overlaps
with the ANUBIS electroweak category, but it is not the same production model:
it is missing `W gamma` fusion and does not keep per-row W/Z/process tags.

Comparison plots should be split by production class before making strong
claims:

- meson/tau-only GRENDEL for forward or beam-dump-like comparisons such as
  FASER/FASER2 and SHiP;
- electroweak-only GRENDEL for ANUBIS/ATLAS/CMS displaced-style comparisons;
- inclusive GRENDEL as a separately labeled "all production modes included"
  result;
- a channel-dominance diagnostic along the exclusion boundary, so each part of
  the contour can be interpreted by its controlling production mode.

Required follow-up:

- publish separate inclusive, meson/tau-only, and electroweak-only sensitivity
  CSVs, or add a channel-filter option to the sensitivity and plotting paths;
- add `W gamma` fusion to the direct-electroweak production set, or quantify
  its absence before using ANUBIS as an electroweak benchmark;
- carry direct-electroweak process identity through LHE-to-CSV conversion so
  W, Z, and future `W gamma` rows can be reweighted and plotted separately;
- update comparison captions so they state the production classes being
  compared, and avoid claiming "GRENDEL competes with ANUBIS/FASER" from the
  inclusive envelope alone.

**Completion test:** the comparison repository can produce committed overlays
for inclusive, meson/tau-only, and electroweak-only GRENDEL curves, plus a
channel-dominance diagnostic. Any ANUBIS comparison either includes
CCDY/NCDY/`W gamma` in the GRENDEL electroweak sample, or carries an explicit
caveat quantifying the missing `W gamma` component.

### Resolution (2026-06-29): PBC conformance, Majorana convention, CL, EW K-factor

Cross-checked against the PBC summary report (arXiv:2505.00947, Fig. 23, BC7
muon-coupled HNL) and the SET-ANUBIS HNL study. Three points that looked like
open conformance gaps are now settled; one production gap (`W gamma` fusion,
per-row W/Z tags, mass-range) remains as written above.

- **Majorana/Dirac (settled: Majorana, self-consistent).** GRENDEL is already
  Majorana on both sides, so this is not a free choice. Production: the HeavyN
  UFO declares `N1` self-conjugate (`name == antiname == 'N1'` in
  `vendored/SM_HeavyN_CKM_AllMasses_LO/particles.py`). Decay/lifetime: HNLCalc
  sums each charged-current mode together with its charge conjugate
  (`HNLCalc.py` ~L1690-1762: `lP`, `lV`, `lud`, `lhad` each append both `mode`
  and `conjugate(mode)`; `llnu` lists both orderings), i.e. a single `N` decays
  to lepton-number +1 and -1 final states -- the Majorana ~2x-width convention,
  not Dirac. This matches BC7 (defined as one Majorana HNL) and ANUBIS's
  minimal-Majorana BC7. P0 item 4 therefore reduces from "determine the
  convention" to "declare Majorana in the result metadata and audit the
  remaining factors of two under it"; the convention itself is no longer open.

- **Confidence level (settled: no change; keep zero-background `N >= 3`).** In
  PBC Fig. 23 the *existing* upper limits carry mixed CLs (90% for the
  beam-dump/fixed-target set, 95% for ATLAS/CMS), but the *projection* curves we
  sit among (ANUBIS, CODEX-b, FLArE, FASER2, SHiP) are not held to a common CL:
  the caption uses line style for background-estimate maturity (solid = data
  extrapolation, dashed = full MC, dotted = toy MC / negligible background), not
  a confidence level. A zero-background `N >= 3` (95%) contour is a legitimate
  member of that set, and `N >= 3` vs `N >= 2.3` (90%) is a ~30% yield shift --
  invisible on a log-log reach boundary. Conformance work is therefore labeling,
  not recomputation: state our contour definition explicitly, and tag GRENDEL's
  curve at the **dotted** maturity tier (zero background assumed; the cosmic
  decay-in-flight is modelled in `../higgs/` but not yet folded into the limit,
  P0 item 3), promoting to dashed only once that background enters the limit.

- **Electroweak K-factor (settled: flat 1.3 for W and Z is sound).**
  `K_FACTOR_EW = 1.3` is an upper-but-in-band NLO single-boson on-peak value;
  W and Z K-factors are driven by the same `q qbar -> V` QCD and differ by only
  a few percent, so a common value is well justified (the W/Z gap is far below
  the order/PDF ambiguity). The mass grid stops at 10 GeV, far below
  `m_W ~ 80.4` / `m_Z ~ 91.2`, so the W/Z are on-shell across the whole grid and
  the on-peak 1.3 is mass-independent here (the off-shell/high-mass K-factor
  rise only matters if the grid is extended toward and above the boson pole).
  Consequently the "unsplit Z scaled by the W K-factor" caveat above is a
  sub-few-percent normalization effect, not a normalization error: the per-row
  W/Z tag is needed for the EW-only-vs-inclusive comparison split and for future
  `W gamma` tagging, not to correct the inclusive normalization.
