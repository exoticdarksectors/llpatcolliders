# Remaining work — BC4 dark scalar (`scalar/`)

This is the authoritative inventory of work still needed before the generated
GRENDEL BC4 (Higgs-portal dark scalar) exclusion island can be presented as a
publication-grade physics result. It is derived from the code paths that produce
the curve, not from comments or earlier planning notes.

Scope is the `scalar/` package and the shared code it imports:
`hnl/analysis/decay_reco_acceptance.py`, `hnl/analysis/exclusion.py`,
`hnl/analysis/_engine.py`, `hnl/production/fonll/`,
`hnl/production/decay_engine/kinematics.py`, `higgs/grendel_geometry.py`, and
`higgs/reco_common.py` — the same ten files `uncertainty_band.py:113-124`
checksums as `CODE_INPUTS`. Items are classified as:

- **P0 — result definition:** missing detector, background, statistical, or
  production-scope inputs that change the *meaning* of the exclusion.
- **P1 — uncertainty & completeness:** missing normalization terms, channels,
  and uncertainty propagation needed for a defensible central curve and band.
- **P2 — validation & reproducibility:** provenance, convergence, and workflow
  checks needed to make the result repeatable and auditable.

An item is not complete merely because a central constant exists in the code. It
is complete when its source, uncertainty, propagation, and regression test are
present.

## What the code currently produces

`run_sensitivity.py` scans an 86-point mass grid (`production.py:46-59`;
`m_S = 0.14`–`4.70 GeV`) against 200 log-spaced `sin^2 theta` points from `1e-12`
to `1e-2` (`run_sensitivity.py:43`), requiring `N_signal >= 3`
(`hnl/analysis/constants.py:21`) at `3000 fb^-1`
(`hnl/analysis/constants.py:12`). The published curve and its provenance are
`data/published/bc4_island.csv` + `MANIFEST.json`; the theory-variation
diagnostic is `data/published/bundle/`.

Production (`production.py:76-122`) samples one shared FONLL bottom `(pT, y, phi)`
pool, rebuilds `B+`/`B0`/`Bs` four-vectors, performs a two-body
`B -> X_s(m_K proxy) S`, and weights by `2 * sigma_FONLL * f_frag * BR(b -> X_s S)`
at `sin^2 theta = 1`. Decay/acceptance (`acceptance.py`) samples a decay mode by
branching ratio, generates a two-body charged pair, and runs the shared GRENDEL
reconstruction and PR #13 selection. The coupling is factored out of both legs and
re-applied by the shared `scan_u2`.

## Resolved / not applicable — do not re-import from the HNL inventory

These are recorded because the HNL `REMAINING_WORK.md` carries items that look
transferable and are not. Verified against this package's code.

**The charm `mu_F` LHAPDF-extrapolation artifact does not affect BC4.** BC4
production is bottom-only and samples charm nowhere: `uncertainty_band.py:289`
filters `entry["quark"] == "bottom"` when discovering the FONLL campaign,
`:236` builds only `..._bottom.dat` grid paths, and `production.py:92` /
`run_sensitivity.py:96` call `get_sigma_total("bottom")` and
`sample_meson_4vectors(..., "bottom", ...)`. The HNL artifact is that
`NNPDF40_nlo_as_01180` has `QMin = 1.65 GeV` while charm's
`mu_F = 0.5*sqrt(m_c^2 + pT^2)` starts at `0.75 GeV` and stays below `QMin` for
all `pT < 2.94 GeV`, so the charm `mu_F = 0.5` scale points are an LHAPDF
extrapolation rather than a perturbative uncertainty. Bottom is structurally
immune: `mu_F = 0.5*sqrt(m_b^2 + pT^2) >= 2.375 GeV` always clears `QMin`. The
BC4 scale envelope is therefore well behaved by construction, and the HNL's
charm-`mu_F` item must **not** be copied here.

**BC4 has no tau production chain.** `scalar/` contains zero references to
`induced_tau` or `tau_decay` (verified by grep over the package). The tau
machinery vendored under `hnl/` is shared infrastructure this package never
calls; production is `b -> X_s S` only. The HNL's "tau spin correlations" item
does not apply. (This is distinct from the `S -> tau tau` *decay* mode, which
does exist and is covered by item 9 below.)

**The `sin^2 theta` factorization is exact, not an approximation.** Every partial
width and `BR(b -> X_s S)` is linear in `sin^2 theta`, so visible branching
ratios are ratios of widths and are coupling-independent
(`model.py:295-301`); the decay MC is therefore built once per mass and
reweighted (`acceptance.py:22-26`). No factor-of-two or convention audit of the
kind the HNL needs (Majorana vs Dirac) applies to BC4: `S` is a real scalar and
the only charge-counting factor is the `2 *` for `b`/`bbar` hadronization
(`production.py:114-115`), matching the FONLL convention.

## P0 — define the detector-level result

Items 1–3 are **shared with every GRENDEL benchmark**, not BC4-specific: BC4
imports the same reconstruction, geometry, selection, and threshold as the HNL
and BC5 analyses rather than reimplementing them. Verified by reading the BC4
call path (`acceptance.py:46-56`), not assumed from the HNL document.

### 1. Add detector response to the decay/reconstruction acceptance

**Current code:** `acceptance.py:96-163` runs a per-event Monte Carlo: for each
`S` four-vector crossing the fiducial volume it samples decay vertices uniformly
along the in-volume path, samples a decay mode by branching ratio, generates
isotropic back-to-back two-body daughters, boosts to the lab, requires both
daughters to clear `P_CUT`, and calls `reconstruct_decays` + `selection_mask`
imported from `hnl/analysis/decay_reco_acceptance.py` (`acceptance.py:46-49`),
which are the `higgs/reco_common.py` single source. Daughter speeds use the true
`beta = p/E` for the timing model (`acceptance.py:148-149`).

Detector response is **geometric wall hits with Gaussian position and time
smearing only**. Grepping `higgs/reco_common.py`, `higgs/grendel_geometry.py`,
`hnl/analysis/decay_reco_acceptance.py`, and `hnl/analysis/exclusion.py` for
material / tracking or vertexing efficiency / trigger / readout / pileup / dead
regions / occupancy returns no implementation and no efficiency constants
anywhere in the chain: there is no per-track or per-channel efficiency factor to
vary. The only response parameters are `HIT_RESOLUTION`, `SIGMA_T_DEFAULT`
(`higgs/reco_common.py:30`), and the PR #13 selection thresholds.

**Required work:**

- implement detector response beyond geometric hits + smearing: tracking,
  vertexing, particle-identification and momentum thresholds, material
  interactions, trigger/readout, and event-selection efficiencies;
- include pileup, dead regions, and occupancy if they affect the proposed
  detector;
- assign per-channel efficiency uncertainties — note that BC4's final states
  span `e e`, `mu mu`, `pi pi`, `K K`, `tau tau`, and hadronic continuum, so a
  single flat efficiency is not adequate;
- charged hadrons (`pi`, `K`) can interact or decay in flight before the outer
  tracker layer; neither is modelled.

**Completion test:** an efficiency map or detector simulation with versioned
inputs reproduces benchmark samples, and `acceptance.py` consumes per-channel
efficiencies with uncertainty variations.

### 2. Validate the geometry and detector configuration

**Current code:** `higgs/grendel_geometry.py` builds a closed tunnel mesh from
one survey polyline and a fixed cross-section (`TUNNEL_ALPHA/BETA/GAMMA/DELTA`,
lines 29-33) with a single `DETECTOR_THICKNESS = 0.24` m shell (line 40) that
simultaneously sets the fiducial inset and the tracker-layer spacing. Rays
originate at `CMS_ORIGIN` (`acceptance.py:49, 176`). The mesh contains no
supports, services, inactive regions, material, alignment uncertainty, or
configurable layout; the parameters are import-time constants.

**Required work:**

- confirm the survey coordinates, CMS IP transform, tunnel cross-section, and
  proposed active volume with the detector/design owners;
- version the geometry and expose layout parameters rather than relying on
  import-time constants;
- model inactive and inaccessible regions and evaluate alignment/survey
  variations;
- independently compare ray intersections against a reference geometry
  implementation. (Partially covered for the *backend*, not the geometry: the
  campaign cross-checked 25,000 Embree rays against the triangle intersector
  with zero hit-mask mismatches — `data/published/bundle/README.md`. That
  validates the intersector, not the mesh.)

**Completion test:** signed-off geometry inputs, reference intersection tests,
and a geometry-variation envelope propagated to the island.

### 3. Supply backgrounds and a statistical model

**Current code:** exclusion is `N_signal >= 3` (`hnl/analysis/constants.py:21`,
consumed via `acceptance.py:54` and `find_exclusion_band_refined`,
`hnl/analysis/exclusion.py:77-181`) at `L_INT_PB = 3e6`
(`hnl/analysis/constants.py:12`) — a zero-background 95% CL approximation. There
is no background count in the statistical model, no control region, no nuisance
parameter, no systematic uncertainty, and no coverage calculation. Unlike the
HNL chain, the BC4 path does not even model a cosmic decay-in-flight background
for reference: `acceptance.py` imports only the signal side. `STATUS.md` records
zero background as an explicit paper working assumption.

**Required work:**

- estimate beam-, collision-, cosmic-, neutrino-, and detector-induced
  backgrounds after the PR #13 selection for BC4's final states;
- define control samples and background uncertainties;
- include luminosity, reconstruction, geometry, signal-model, and finite-MC
  nuisance parameters;
- choose and implement the statistical prescription (e.g. CLs or a documented
  Poisson construction) and validate expected coverage;
- define how the open upper edge in the `0.14–0.20 GeV` electron-only interval
  is reported (`MANIFEST.json` `topology.upper_open_grid_run_GeV`).

**Completion test:** the final likelihood/configuration is versioned, has toy or
asymptotic validation, and produces expected limits with uncertainty bands.

### 4. Fix the production scope: the result is B-produced BC4, not inclusive BC4

**Current code:** `K -> pi S` **is implemented in the model layer but is not in
the central result.** `model.br_K_to_pi_S` (`model.py:389-402`) computes the
Winkler eq. A8/A9 branching ratio and is gated to `m_S < m_K - m_pi`
(`model.py:400`), and `tests/test_model.py::test_K_to_pi_S_only_below_threshold`
exercises it. But no production or analysis code calls it: `production.py`
imports only `br_B_to_Xs_S` (`production.py:108`), `B_SPECIES` covers the
b-hadron pool `521`/`511`/`531`/`5122` (`production.py:67-71`) but no light kaon,
and `generate_scalar_4vectors` returns
empty above the `B -> K S` ceiling (`production.py:89`). Its own docstring states
"Sampling a kaon flux at the LHC IP is deferred" (`model.py:393`). The published
`MANIFEST.json` `physics_inputs.production` says "inclusive `b -> X_s S` ...
summed over `B+`, `B0`, and `Bs`" and does not mention kaons.

The published curve is therefore **B-produced BC4**. This matches the present
comparison convention (`EXTERNAL_INPUTS_NEEDED.md`: the FASER2 curve is included
precisely because "its production is B-driven and applicable to BC4"), so the
current result is self-consistent — but the scope is implicit rather than
declared, and a reader may assume the low-mass end is inclusive.

A complete kaon contribution is **not** a total-rate correction. Kaons are
long-lived (`TAU_KPLUS`, `model.py:86`), so both the `S` momentum spectrum and
the `S` production *vertex* depend on where the kaon decays; the analysis rays
every scalar from `CMS_ORIGIN` (`acceptance.py:176`), which is only valid for
prompt parents. Adding kaons needs: a validated 14 TeV charged **and** neutral
kaon spectrum, kaon transport (magnetic bending, material, interaction losses)
with a sampled decay point, species and branching fractions, and transverse
acceptance from the displaced origin.

**Scope note — what this cannot change.** The affected region is only
`m_S < m_K - m_pi ~= 0.354 GeV`. It does **not** touch the headline reach
(`sin^2 theta = 6.5193e-12` at `m_S = 0.975 GeV`) or the `~3.825 GeV` high-mass
closure (`MANIFEST.json` `headline_reach`). It can only extend or deepen the
low-mass end.

**Required work — choose one and state it:**

- **(a) Keep the result defined as B-produced BC4** (recommended, matching the
  current comparison set): declare the production scope in `MANIFEST.json` as a
  machine-readable field, label the figure and paper text accordingly, and state
  that `m_S <~ 0.354 GeV` carries an unquantified kaon contribution that would
  only strengthen the limit; **or**
- **(b) Add kaons as a dedicated low-mass extension** with the full flux,
  transport, decay-location, and acceptance chain above, published as a
  separately labelled curve so the B-produced comparison stays intact.

Do not "correct the rate" by scaling: without transport and decay-location
modelling the added yield would be assigned the wrong origin and the wrong
acceptance.

**Completion test:** the published manifest and figure caption state the
production scope unambiguously, and either (a) the kaon-affected mass interval
is explicitly bounded, or (b) a transport-validated kaon channel is published
with its own uncertainty envelope.

## P1 — uncertainty and completeness

### 5. Propagate the `b -> X_s S` production-normalization uncertainty

**Current code:** the entire production normalization is
`model.br_B_to_Xs_S` (`model.py:366-379`), the Winkler eq. A7 **spectator
estimate** `Gamma = g^2 (m_B^2 - m_S^2)^2 / (32 pi m_B^3)`, built from
`g_phisb` (`model.py:323-328`). Every input to that coupling is a hardcoded
constant with no uncertainty and no source table: `M_TOP = 172.76` (`:57`),
`V_TB, V_TS = 0.99915, 0.0404` (`:58`), `M_B_QUARK = 4.18` (`:63`), and
`V_HIGGS` from `G_F`. The only test is
`test_production_br_matches_literature` (`tests/test_model.py`), which pins
`br_B_to_Xs_S(0.5, "B+") ~= 5.3 +- 0.8` — a ~15% tolerance that catches factor-2
errors, not a propagated uncertainty.

The 109-variation campaign (`uncertainty_band.py:1298`) varies the FONLL
**parent spectrum and cross-section** only. It does not vary this branching
ratio at all. Because `BR(b -> X_s S)` is a multiplicative normalization on
every event weight (`production.py:115`), its uncertainty maps directly onto the
lower island edge and is currently absent from the band.

**Partial-coherence finding (verify intent before acting):** the campaign's
bottom-mass axis (`mb_dn`/`mb_up`, `uncertainty_band.py:308-310`) swaps only the
FONLL grid via the `_bottom_grid` env-var context manager
(`uncertainty_band.py:372-382`). `model.M_B_QUARK`, which enters `g_phisb`
linearly and hence the yield quadratically, is untouched — so the `m_b` axis is
a production-kinematics variation only. This may be deliberate (the FONLL grids
use a pole mass; `model.py:62` documents `M_B_QUARK` as MSbar, and the two
schemes should not be varied with a common value), but the code records no such
statement, and the model-side `m_b` currently carries no uncertainty in any
direction. **Cannot be determined from code which is intended** — resolve with
the author before either coupling the axes or documenting the separation.

**Required work:**

- store source, version, uncertainty, and correlations for `m_t`, `V_ts`,
  `V_tb`, and the coupling-side `m_b`/renormalization scheme beside each value;
- assign and propagate a theory uncertainty on the eq. A7 spectator inclusive
  rate itself (its accuracy relative to a summed exclusive calculation is not
  assessed anywhere in the package);
- resolve the `m_b` scheme question above and either add a coherent coupling-side
  `m_b` axis or document why the axes are separate;
- add the resulting normalization interval to the published envelope as its own
  source.

**Completion test:** a machine-readable input table generates the constants, and
a production-normalization interval appears in
`bc4_single_source_variation_envelope.csv` alongside `scale`, `pdf`,
`bottom_mass`, and `decay_model`.

### 6. Pin and checksum the Winkler width table; bound the fixed decay constants

**Current code:** the central decay model is a 6,010-row digitization of Winkler
arXiv:1809.01876 Fig. 4 in `data/winkler_widths.csv`, loaded by
`model._load_winkler` (`model.py:116-134`) and interpolated linearly in
`log(m)`–`log(Gamma)` (`model.py:137-144`), returning `0.0` outside each curve's
digitized range. It is used below `M_SPECTATOR = 2.0 GeV` for `pipi`/`KK`/`4pi`
and above it for `ss`/`cc`/`gg` (`model.py:277-286`).

Three concrete provenance gaps, all verified:

1. **The referenced extractor documentation does not exist.** `model.py:106`
   says "see `scalar/data/winkler_widths.csv` + its extractor README". There is
   no README in `scalar/data/`; the directory holds only `winkler_widths.csv`,
   `published/`, and (in the llpatcolliders worktree only) `competitors/`. No
   digitization tool, figure revision, extraction procedure, or estimated
   digitization error is recorded anywhere.
2. **The table is not checksummed by anything.** `CODE_INPUTS`
   (`uncertainty_band.py:113-124`) hashes ten `.py`/`.yml` files and
   `_validate_recorded_code_state` (`:255-276`) proves them against the recorded
   commit — but `scalar/data/winkler_widths.csv` is not in that tuple. It is
   also absent from `data/published/MANIFEST.json` and from
   `data/published/bundle/UNCERTAINTY_MANIFEST.json` (which references only the
   string `"width_scheme": "winkler"`), and no test hashes it. **Editing the
   central decay input would invalidate no recorded checksum and fail no test**,
   while `producer_code_sha256` would still validate.
3. **No digitization uncertainty is propagated.** The decay-model source in the
   published envelope is the *interval to the `chpt_spectator` alternate*
   (`uncertainty_band.py:1318-1322`) — a model-choice envelope. Read-off error
   on the Fig. 4 curves is a separate, unpropagated term.

Additionally, fixed constants in the analytic sector carry no variation:
`ALPHA_S = 0.30` (`model.py:101`, a single value across the whole 0.5–5 GeV
window, entering `width_gg` quadratically — its own comment concedes only "~10%
level"), and `C_4PI = 5.1e-9` (`model.py:96`), a constant *tuned* to soften the
2 GeV seam. The `chpt_spectator` alternate exercises these jointly with the
scheme change, so their individual effect is not separable from the current band.

**Required work:**

- write the extractor README `model.py:106` promises, or remove the reference:
  source figure, revision, digitization tool, per-curve mass range, and
  estimated read-off error;
- add `scalar/data/winkler_widths.csv` to `CODE_INPUTS` and to the published
  manifests, and add a regression test pinning its SHA-256;
- propagate a digitization uncertainty as its own envelope source, separate from
  the `chpt_spectator` model interval;
- vary `ALPHA_S` (or run it to `m_S`) and `C_4PI` and bound their contributions;
- cross-check the digitized widths against an independent evaluation of the
  Winkler dispersive result rather than against our own re-derivation;
- make `_winkler_width`'s silent `return 0.0` outside the digitized range
  (`model.py:142-143`) either an error or an audited, tested extrapolation — a
  channel that falls off the table currently vanishes without a warning.

**Completion test:** the width table is hash-pinned in the manifests and tests,
its provenance is machine-readable, and a digitization interval appears in the
published envelope.

### 7. Bound or extend the FONLL `pT` and rapidity truncation

**Current code:** every grid, central and varied, stops at `pT = 50 GeV` and
`|y| = 3`. This is fixed in the filename stem consumed by both the central path
(`hnl/production/fonll/fonll_parser.py:24`,
`..._dsdpTdy_pt0-50_y-3to3_central_bottom.dat`) and the campaign
(`uncertainty_band.py:100`, `GRID_STEM = "..._pt0-50_y-3to3"`). The sampler
draws only from within the grid, so the omitted tail contributes zero and no
bound on it exists anywhere in the package.

This matters more for BC4 than a pure rate argument suggests: GRENDEL sits above
IP5 at large angle, so acceptance is not monotonic in parent `pT`, and the
canonical run deliberately importance-samples *toward* high `pT`
(`--high-pt-tilt-scale 5`, `README.md`; `production.py:157-161`) — i.e. the
configuration itself indicates the high-`pT` region carries acceptance weight.
The tilt reweights within the grid; it cannot recover the region beyond it.

**Required work:**

- extend the `pT` and rapidity grids, or place a quantitative bound on the
  omitted contribution **after GRENDEL acceptance** (not at production level);
- because the tilt already concentrates statistics near the `pT` ceiling, report
  the accepted-yield fraction in the top grid bins as the diagnostic that decides
  whether extension is needed.

Use the same three-column rectangular format as the committed grids:

```text
# pT_GeV  rapidity  d2sigma_dpTdy_pb_per_GeV
```

**Completion test:** either extended grids are in use, or a committed
acceptance-folded bound shows the omitted tail is below a declared tolerance on
both island edges.

### 8. b-baryons: DONE (measured 2026-07-17). Shared shape / static fragmentation: residual

**Status:** the b-baryon pool is now included. `production.py` sums `B+`, `B0`,
`Bs`, **and `Lambda_b`** (PDG 5122, the `FRAG_LAMBDA_B = 0.18755` fraction that
lumps the whole b-baryon pool). Previously `production.py` summed only the three
mesons (`0.36206 + 0.36206 + 0.08834 = 0.81246`), dropping 18.75% of the
b-hadron pool. Because `b -> s S` is spectator-independent
(`model.br_B_to_Xs_S` takes only `(m_parent, tau_parent)`), the baryon enters
the rate on the same footing as a meson; only the recoil mass and lifetime
differ. `br_B_to_Xs_S` now accepts `parent="Lambda_b"`, and `B_SPECIES[5122]`
uses the **Lambda recoil** `m_Lambda = 1.1157 GeV` (the lightest strange baryon,
the baryonic analogue of the kaon recoil for mesons), **not** `M_KPLUS` — with
`m_Lambda` the Lambda_b closes at `5.62 - 1.12 = 4.50 GeV`; a kaon recoil would
wrongly leak it to `5.13 GeV`.

**Measured impact (controlled A/B/C scan, `n_pool=1e5`, `n_samples=50`; config A
reproduces the published `u2_min` to 0.4%, so the shifts are trustworthy):**

| m_S (GeV) | Lambda_b shift in `u2_min` | dex | recoil error (m_Lambda vs m_K) |
|---|---|---|---|
| 0.975 (deep reach) | **-11.3%** | -0.052 | +0.3% |
| 2.000 (mid) | **-9.9%** | -0.045 | -0.2% |
| 3.500 (near closure) | **-13.9%** | -0.065 | -0.8% |

So Lambda_b **strengthens** the lower edge by 10-14% (`-0.045` to `-0.065` dex)
across the sensitive range, growing toward closure (the heavier parent has more
phase space at high `m_S`: `(m_Lambda_b^2 - m_S^2)^2` beats `(m_B^2 - m_S^2)^2`).
This is well above the `~0.007` dex finite-edge MC noise and is the one omission
in this document that moves the central curve. It is a strengthening: the
published limit was conservative in normalization by this amount.

**On the recoil-mass choice:** using `m_Lambda` is physically correct, but the
measurement shows the choice is numerically **immaterial** for this island
(`<= 0.8%`, `<= 0.004` dex, even at 3.5 GeV). The region where `m_Lambda` and
`m_K` diverge (the Lambda_b kinematic ceiling, 4.5-5.1 GeV) sits far above the
island's 3.7 GeV closure, and at lower `m_S` the S spectrum is recoil-insensitive
at a ~5.6 GeV parent. The earlier concern that a kaon recoil would make the
`-0.045` dex optimistic is therefore not borne out numerically — but `m_Lambda`
remains the correct implementation.

**Full-consistency republication: DONE (2026-07-18).** The published set is now
`Lambda_b`-consistent and hash-linked: the 400k central island (deepest reach
`6.519e-12` vs the meson-only `7.399e-12`, closure `3.70 -> 3.80` GeV, +1
sensitive point), plus a matching **109-variation + 2 numerical-control
uncertainty bundle** at 100k parents per species (`data/published/bundle/`), a
`Lambda_b` main `MANIFEST.json`, and updated `README`/`STATUS`.
`test_published_curve.py` passes (central/bundle hash linkage + exact closure
interpolation). Reproduce with:

```sh
# central (400k):
python -m scalar.run_sensitivity --n-pool 400000 --n-samples 100 --seed 42 \
    --high-pt-tilt-scale 5.0 --nominal-mixture-fraction 0.5 --force-produce \
    --output tmp/bc4_island_lambdab.csv --vector-dir tmp/llp_4vectors_lambdab
# bundle (100k, needs ~270 GB scratch):
python -m scalar.uncertainty_band run  --n-pool 100000 --n-samples 100 --seed 42 --workers 12 --scratch-dir <SCRATCH>
python -m scalar.uncertainty_band collect --n-pool 100000 --n-samples 100 --seed 42 --scratch-dir <SCRATCH>
```

Note the bundle uses 100k parents per variation (vs the 400k central) because
the full campaign at 200k overruns scratch disk; the lower edge stays fully
physics-dominated, but on the noisier upper edge the numerical control is
marginal at 2.0 and 2.3 GeV (recorded in the bundle manifest and the test). Two
small refreshes remain: the six-million-event high-mass closure control and the
`paper/figures/bc4_exclusion.pdf` adapter output are still the meson-only
baseline (they do not affect the published central or band).

**Remaining (item 8 residual):** Lambda_b lumps `Xi_b`/`Omega_b` (different mass
and lifetime, a small fraction of b-baryons); a per-species baryon split and a
`pT`-dependent fragmentation model are still open, and the shared-shape and
static-fraction approximations below persist.

Two further approximations sit on the same code path:

- **Shared shape.** All three species draw from one bottom `(pT, y)` pool and
  differ only by the mass used to rebuild `E`/`pz`
  (`hnl/production/fonll/meson_sampler.py:45-58`, "the shared-shape
  approximation"). `Bs` in particular does not have its own spectrum.
- **Static fractions.** `FRAG_B` values are constants "measured in
  acceptance-specific inputs, applied over the full FONLL grid as constants"
  (`hnl/production/constants.py:64-65`), extrapolated over the whole grid with
  no `pT`/`y` dependence and no covariance.

`hnl/production/constants.py:82-85` further notes the remainder is not a measured
b-baryon fraction: `Xi_b`/`Omega_b` are neither separated nor validated as
`Lambda_b` equivalents.

**Required work:**

- add the b-baryon contribution (at minimum a `Lambda_b`-like entry using its own
  `m` and `tau`), or bound its omission quantitatively on both island edges;
- provide differential production fractions/spectra for `Lambda_b`, `Xi_b`, and
  `Omega_b`, or document the closure convention;
- provide species-specific `d2sigma/dpTdy` grids or differential fragmentation
  fractions;
- propagate the experimental fraction covariance and its `pT`/rapidity
  dependence.

**Completion test:** bottom species close without an unidentified remainder, the
fractions reproduce their source measurements in the quoted acceptance, and the
fraction covariance produces a correlated island variation.

### 9. Complete the scalar decay final-state model in the acceptance

**Current code:** `acceptance.py:66-76` maps each decay mode to a single
`(daughter mass, charged-track fraction)` pair and generates one isotropic
two-body pair. The approximations, in the code's own terms
(`acceptance.py:61-65`):

- `pipi` -> charged fraction `2/3` (removes `pi0 pi0`); `KK` -> `0.5` (removes
  `K0 K0bar`). These are isospin counting factors, applied as a random accept
  (`acceptance.py:119`) rather than by generating the neutral final state.
- `ss`, `cc`, `gg`, and `4pi` are all modelled as **two leading charged pions**
  with charged fraction `1.0`. Real multi-body hadronic final states have >2
  tracks, softer per-track momenta (the `P_CUT` requirement at
  `acceptance.py:153` bites differently), and a different opening-angle
  distribution — all of which feed the PR #13 selection. The code flags this as
  "a documented approximation; these channels are sub-dominant in the < 2 GeV
  region that drives the reach"; that is a correct statement about the *deepest*
  point (`0.975 GeV`) but these channels are exactly the ones that set the
  high-mass end, where the island closes (`~3.825 GeV`).
- `tautau` -> two taus at `M_TAU` with charged fraction `1.0`
  (`acceptance.py:69`, "tau directions proxy the leading tracks"). The taus'
  own decays are not generated: `BR(tau -> 1 charged prong) ~ 0.85`, so
  `1.0` overstates the charged yield, and the tau flight distance and the
  daughter-vs-tau direction difference are both ignored. `tautau` opens above
  `2 m_tau ~= 3.55 GeV` (`tests/test_model.py::test_spectator_region_opens_gg_cc_tautau`),
  i.e. inside the high-mass closure region.
- Charged `pi`/`K` daughters are stable in the model; decay-in-flight and
  hadronic interaction before the outer layer are not modelled (see item 1).

**Required work:**

- generate real hadronic multi-body final states (or a validated spectral/shower
  model) for `ss`, `cc`, `gg`, and `4pi`, and feed the actual track multiplicity
  and momenta into the best-two-track selection;
- decay the taus (validated generator or spectral functions) instead of using
  them as track proxies;
- generate neutral sub-modes explicitly rather than applying a scalar
  charged-fraction accept, so the selection sees the true final state;
- quantify the resulting shift, prioritizing `m_S > 2 GeV` where these channels
  dominate and where the high-mass closure is defined.

**Completion test:** the visible-final-state model is validated against an
independent generator, and the closure mass is stable under the replacement
within a declared tolerance.

### 10. Add the FONLL alpha_s companion (completeness only — measured negligible)

**Current code:** `uncertainty_band.py:1369` records the limitation verbatim:
`"No FONLL alpha_s companion grids are included."` Unlike the HNL chain — which
has a supported (if hand-driven) `combine_band.py --alphas-lo/--alphas-hi` fold
— the BC4 package has **no code path of any kind** to consume alpha_s companion
curves: `discover_variations` (`uncertainty_band.py:284-339`) accepts only
`variation_kind` in `{scale, pdf, mass}` and hard-asserts
`{"central": 1, "scale": 6, "pdf": 100, "mass": 2}` (`:337-339`), and
`_grid_path` (`:236`) rebuilds every path from the fixed
`GRID_STEM = "..._as_01180_..."` (`:100`) — the `as_01170`/`as_01190` companions
differ **in the stem**, so they are unreachable even if present. Adding them
means parametrizing the stem per entry and relaxing the count assertion.

**Do not over-prioritize this.** The PDF4LHC alpha_s half-spread on the
trapezoid-integrated **bottom** grids is **1.49%**. The terms it would join in
the BC4 envelope are the bottom scale envelope (`~+-42%`) and the bottom PDF
1sigma (`3.9%`). In quadrature against `42%` it moves the band by
`sqrt(42^2 + 1.49^2) - 42 ~= 0.03` percentage points — below the campaign's own
numerical-control scatter and invisible on any plot. Run it for PDF4LHC
prescription completeness and provenance, not for physics. It should be
scheduled behind items 5, 6, and 8, each of which is a larger effect on the same
edge. (Note: the HNL's charm alpha_s figure of 0.44% is irrelevant here — BC4
never samples charm; see the resolved section.)

**Required work:**

- parametrize the grid stem per manifest entry rather than from the module-level
  `GRID_STEM`, and source paths from the manifest's `path` field;
- relax the campaign completeness assertion to accept the two coherent alpha_s
  points (alpha_s is coherent across quarks, so it is structurally a scale
  point);
- run the two companion chains and add the PDF4LHC half-difference as an
  envelope source.

**Completion test:** `UNCERTAINTY_MANIFEST.json` no longer lists the alpha_s
limitation, and the envelope carries an alpha_s source with its measured
(negligible) size.

## P2 — validation and reproducibility

### 11. Make the geometry cache self-validating

**Current code:** `acceptance._geometry` (`acceptance.py:166-184`) keys the
ray-cast cache **only by the four-vector CSV stem** (`geom_{csv_path.stem}.npz`)
and invalidates it only when the CSV is newer than the NPZ (mtime guard,
`acceptance.py:172`). The write is atomic (tmp + `fsync` + `os.replace`,
`:177-183`), which is good. But the cache embeds no checksum of the input
vectors, no geometry version, no origin, and no selection configuration:
**editing `higgs/grendel_geometry.py` does not invalidate any cache**, and a
re-scan will silently reuse ray-casts from the previous mesh. `run_sensitivity.py`
exposes no `--force-geometry` equivalent (`run_sensitivity.py:155-170`); the only
recourse is deleting the cache directory by hand.

The campaign path is safer by construction — each variation gets its own
`runs/<variation>/geometry_cache/` (`uncertainty_band.py:23-28`) and is never
shared across variations (`independence_policy`, `:1292-1296`) — so this is an
interactive/dev-run hazard, not a defect in the published campaign.

**Required work:**

- key caches by input checksum plus a geometry/analysis configuration hash
  (`CODE_INPUTS` already hashes `higgs/grendel_geometry.py` — reuse it);
- add a force-regenerate flag to `run_sensitivity.py`;
- verify no stale cache from another run tag can be consumed.

**Completion test:** editing the geometry module or the input vectors
invalidates the cache automatically.

### 12. Quantify the remaining numerical convergence axes

**Current code — what is already done (do not redo):** the published curve
records at least 2,249 effective events on every finite sensitive lower edge and
at least 2,049 at every sensitive peak; independent six-million-event **meson-only**
controls at 3.75/3.80/3.85 GeV place the meson-only closure within 0.006 GeV of the
meson-only central 3.7975 GeV, so that endpoint is **not** a finite-pool artifact
(`MANIFEST.json` `numerical_validation`; `data/published/README.md`). The published
Lambda_b closure (3.825 GeV) rests on the same >=2,049 effective events per peak; a
matching Lambda_b 6M control is a pending refresh. Two fresh-seed same-physics
central repeats are carried as numerical controls, reported separately and excluded
from the envelope (`uncertainty_band.py:1331-1334`); `test_published_curve.py`
asserts they are subdominant to the physical sources at every mass except `u2_max`
at `2.0` and `2.3 GeV`.

**What is still open:**

- the **upper**-edge statistics are much weaker than the lower edge: minimum
  sensitive upper-edge event ESS is **201** (`MANIFEST.json`
  `numerical_validation.minimum_sensitive_upper_edge_event_ess`), an order of
  magnitude below the lower edge, and the one flagged control mass is an
  `u2_max` point;
- `n_samples = 100` decay vertices per entering scalar
  (`run_sensitivity.py:159`) is not demonstrated converged against an exact
  quadrature of the decay-in-volume integral; there is no published
  `n_samples` scan (contrast the HNL chain's exact-50 vs exact-100 controls);
- the `sin^2 theta` scan is 200 fixed log points (`run_sensitivity.py:43`)
  followed by refinement in `find_exclusion_band_refined`
  (`hnl/analysis/exclusion.py:77`); the refinement is present, but its tolerance
  is not reported in the manifest;
- mass-grid stability: spacing ranges from 20 MeV to 100 MeV
  (`production.py:46-59`) and the published closure at 3.825 GeV is a
  **log-yield interpolation between the 3.80 and 3.90 GeV rows**
  (`MANIFEST.json` `closure_note`, reproduced in
  `test_published_curve.py`). Do not quote sub-grid precision beyond what the
  six-million-event controls independently support;
- four-vector CSVs are written at eight significant digits
  (`production.py:139-140`); the effect on reconstructed geometry for
  ultra-boosted low-mass scalars is unquantified.

**Required work:** raise or bound the upper-edge statistics; publish an
`n_samples` and seed scan; record the refinement tolerance; demonstrate
stability under mass-grid refinement near 2 GeV (the scheme handover) and near
the closure.

**Completion test:** both island edges are stable under sample-count, seed, and
mass-grid refinement within a declared numerical tolerance.

### 13. Pin the software environment and external tools

**Current code:** `uncertainty_band.py` records a great deal already — git HEAD
(`:243-247`), SHA-256 of the ten `CODE_INPUTS` files validated against the named
commit (`:255-276`), python/numpy/pandas versions, the ray backend, and the
platform (`:1361-1367`), plus per-grid FONLL checksums verified at discovery
(`:317-320`). What is **not** pinned: the FONLL grids' own generation
provenance lives in an external workspace referenced only by
`--grid-dir` and a manifest hash (`:1308`); `hnl/environment.yml` is hashed but a
full environment lock is not committed; and the central `run_sensitivity.py`
path (as opposed to the campaign) records no manifest at all — it writes
`bc4_island.csv` with no run metadata, and publishing is a manual
`cp` + hand-edit of `MANIFEST.json` (`data/published/README.md`).

The campaign requires `embreex==4.4.0` (`data/published/bundle/README.md`); that
pin lives in prose, not in a lock file.

**Required work:**

- commit an environment lock (including `embreex`) and record exact external
  tool versions;
- have `run_sensitivity.py` emit a run manifest (code revision, seeds, input
  checksums, config) so the publish step can be mechanical rather than a manual
  copy + hand-edited JSON;
- add a one-command clean-room reproduction of a representative mass point.

**Completion test:** a fresh checkout reproduces benchmark weights and island
points within the stated stochastic tolerance, and publishing does not require
hand-editing a manifest.

### 14. Extend regression coverage

**Current code:** `tests/` covers more than the HNL package does at the
equivalent stage. `test_model.py` pins the vev, threshold behaviour, BR
normalization for both schemes, the `f0(980)` model difference at 0.975 GeV
(`central > 5x alternate`), coupling linearity/inverse-linearity, the
`B -> K S` kinematic close, and `K -> pi S` gating. `test_published_curve.py`
hash-links the published CSV to its manifest, verifies the mass grid matches
`production.MASS_GRID` exactly, re-derives the closure interpolation, and
asserts the bundle's axis counts (`central 1, scale 6, pdf 100, mass 2,
decay_model 1, numerical_control 2`) and hash links.

Gaps, all verified:

- **no test hashes `data/winkler_widths.csv`** (item 6) — the central decay
  input is unpinned by the test suite;
- `test_production_br_matches_literature` uses a ~15% tolerance; it pins
  factor-2 class bugs only, and cannot tighten until item 5 supplies a
  benchmark rate with a stated uncertainty;
- no test covers `production.generate_scalar_4vectors` weight normalization
  (that `sum(w)` reproduces `2 * sigma * sum(f_frag) * BR`), which is the single
  factor that makes `N_signal` absolute;
- no test covers the `acceptance._MODE_DAUGHTER` charged-fraction bookkeeping
  against `model.branching_ratios`.

**Required work:** add the width-table hash test; add a production weight-closure
test; tighten the normalization tolerance once item 5 lands.

### 15. Propagate all variations into the final plots and tables

**Current code:** `plot_exclusion.py` draws the island from the CSV with
honest segment handling — `_segments` yields only contiguous sensitive runs and
never bridges insensitive masses (`plot_exclusion.py:35-40`). It is explicitly
**not** the paper figure source (`plot_exclusion.py:1-8`): the publication
comparison is produced elsewhere (`shared/curves_PBC` in the llpatcolliders
worktree; `competitor-curves` in the standalone bundle — see the divergence note
below), and the package-local competitor loader is a labelled legacy diagnostic.

The published envelope
(`data/published/bundle/bc4_single_source_variation_envelope.csv`) carries
`scale`, `pdf`, `bottom_mass`, and `decay_model` as **one-source-at-a-time**
intervals, outermost-endpoint display, rebased onto the canonical central curve,
with nothing added in quadrature (`uncertainty_band.py:1310-1336`). The bundle
README and the manifest `limitations` both state it is not a confidence band.
That honesty is correct and should be preserved. Absent: detector, background,
luminosity, and numerical bands (P0 items 1–3, item 12), the production
normalization (item 5), the width digitization (item 6), and alpha_s (item 10).

The topology handling is already conservative: variations that restore/remove
sensitivity or flip an open-edge state are listed by variation and mass rather
than converted into a finite ribbon (`bundle/README.md`;
`physical_variation_topology_summary` flags `3.8 GeV`). Preserve that property
when adding sources.

**Required work:**

- add the detector, background, luminosity, and numerical bands once those inputs
  exist;
- add the item 5/6/10 sources to the envelope registry;
- if a combined band is ever wanted, define and validate a correlated
  combination — do not quadrature the current single-source intervals;
- keep the primary proposed-experiment plot central-only until a defensible
  probabilistic band exists.

**Completion test:** every displayed band traces to stored variation outputs and
regenerates from the manifest, and the topology guard still refuses to ribbon
across a sensitivity change.

## Divergence between the two `scalar/` packages

Verified by `diff -rq` on 2026-07-17. The two copies —
`llpatcolliders_BC4_PR17/scalar/` and `grendel_scalar/scalar/` — are **identical
in every physics-bearing file**: `model.py`, `production.py`, `acceptance.py`,
`run_sensitivity.py`, `uncertainty_band.py`, all of `tests/`,
`data/winkler_widths.csv`, and all of `data/published/` (curve, manifests,
bundle). Every claim in this document holds equally in both.

The differences are path/prose only, and none affects an item above:

- `plot_exclusion.py`, `STATUS.md`, `EXTERNAL_INPUTS_NEEDED.md`: the
  authoritative comparison renderer is named `shared/curves_PBC` in the
  llpatcolliders worktree and `competitor-curves` in the standalone bundle.
- `README.md`: the llpatcolliders copy hardcodes a local conda interpreter path;
  the bundle copy uses a placeholder.
- `data/competitors/` and `tmp/` exist only in the llpatcolliders worktree. The
  bundle copy has no `data/competitors/`, so `plot_exclusion.py`'s legacy loader
  finds nothing there and skips with a note — harmless, since it is not the paper
  figure source.

## Recommended execution order

1. **Declare the production scope (item 4a)** — cheap, and it is the only P0
   item that is a decision rather than a campaign. Doing this first prevents the
   B-produced result from being read as inclusive at low mass.
2. **Pin the decay input (item 6, provenance half)** — write the extractor
   README, add `winkler_widths.csv` to `CODE_INPUTS` and the manifests, add the
   hash test. Hours of work; it closes the largest audit hole in the package
   (the central decay input is currently editable without tripping any check).
3. **Define the detector signal, geometry, backgrounds, and statistical model
   (items 1–3)** — shared across GRENDEL benchmarks; coordinate with the HNL and
   BC5 owners rather than solving locally. These gate the *meaning* of the
   result.
4. **Fix the production normalization and completeness (items 5, 8)** — the
   `b -> X_s S` uncertainty and the omitted 18.75% b-baryon fraction are the two
   largest un-propagated terms on the lower edge, and item 8's `Lambda_b` entry
   is a small, well-defined code change.
5. **Complete the decay final-state model (item 9)** and bound the FONLL
   truncation (item 7) — both concentrate on `m_S > 2 GeV`, i.e. the high-mass
   closure.
6. **Finish uncertainty and digitization propagation (items 6 remainder, 10)** —
   alpha_s last; it is measured at 1.49% against a 42% scale envelope.
7. **Close out reproducibility (items 11–15)** — cache keying, the remaining
   convergence axes (upper-edge statistics first), environment lock, tests, and
   envelope propagation.

Until the P0 items are resolved, the variation envelope in
`data/published/bundle/` is a theory/model diagnostic of the current simplified
analysis, not an uncertainty on a fully defined experimental exclusion. The
package's own manifests say so; keep it that way.
