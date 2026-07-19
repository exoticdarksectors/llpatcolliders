# `hnl/` - Heavy Neutral Lepton production and sensitivity

This package runs the full GRENDEL HNL chain for proton-proton collisions at
14 TeV:

1. produce HNL four-vectors from mesons, kaons, induced taus, prompt taus,
   and electroweak W/Z processes;
2. combine the production channels for each flavor and mass;
3. ray-cast the HNL trajectories through the GRENDEL detector geometry;
4. for the events that cross the fiducial volume, sample flavor-dependent
   FairShip decays, reconstruct the visible tracks with the shared GRENDEL
   reconstruction, and apply the PR #13 selection to get the acceptance;
5. scan the active-flavor mixing and produce sensitivity tables and plots.

The supported flavor hypotheses are `Ue`, `Umu`, and `Utau`. The committed
mass grid contains 123 points from 0.2 to 10 GeV (incl. a 3.62-3.70 GeV
refinement that resolves where the exclusion islands pinch shut).

## Layout

```text
hnl/
|-- README.md
|-- environment.yml
|-- config_mass_grid.py
|-- run_all.py                 production orchestrator
|-- run_analysis.py            sensitivity entry point
|-- run_full_all.sh            complete three-flavor workflow
|-- analysis/                  FairShip decay templates, GRENDEL acceptance,
|                              exclusion band, plots
|-- data/
|   `-- production/fonll/central/   committed central FONLL grids (consumed)
|-- production/
|   |-- combine_channels.py
|   |-- decay_engine/          meson, baryon, kaon, and tau decays
|   |-- fonll/                 heavy-meson spectrum parser and sampler
|   `-- madgraph/
|       |-- run_tau_production.py
|       |-- run_wz_production.py
|       `-- run_wz_sharded.py
|-- tests/
|-- vendored/
|   |-- HNLCalc/               production BRs, differential rates, lifetimes
|   |-- fairship/              FairShip HNL decay modules (rest-frame sampling)
|   `-- SM_HeavyN_CKM_AllMasses_LO/   MadGraph HeavyN model
`-- tmp/                       generated runs, caches, decay templates; git-ignored
```

The GRENDEL detector geometry and 4-hit reconstruction are not duplicated
here: the analysis imports `grendel_geometry` and `reco_common` directly from
the sibling `../higgs/` package, the single source on
exoticdarksectors/llpatcolliders `main` (PR #13).

## Physics Chain

Production includes:

- charged and semileptonic decays of `K`, `D`, `Ds`, `B`, `Bs`, and `Bc`;
- b-baryons through `Lambda_b -> Lambda_c l N` (the "Bbaryon" channel);
- induced taus from the leptonic `Ds/D+/B+/Bc -> tau nu` modes and the
  semitauonic `B/Bs -> D(*)/Ds(*) tau nu` and `Lambda_b -> Lambda_c tau nu`
  modes (polarization treatment under "Production Decay Kinematics" below);
- prompt taus from `W -> tau nu` and `gamma*/Z -> tau tau`;
- `W/Z -> ell N` through MadGraph.

The channel labels combined per mass are `Bmeson`, `Dmeson`, `Bc`,
`Bbaryon`, `tau`, `induced_tau`, `Kmeson`, and `WZ`
(`production/combine_channels.py::CHANNELS`).

The meson spectra come from the committed 14 TeV FONLL tables. HNLCalc
provides production branching ratios, differential three-body rates, and
tau-decay rates. The HNL lifetime and visible final states used downstream come
from the FairShip decay modules instead (see "Signal Decay and Acceptance").

For a production row generated at unit mixing, the downstream yield is

```text
N_signal = L_int * sum_i [
    weight_i * U_alpha^2 * P_decay_i(U_alpha^2) * acceptance_i
]
```

The analysis implements this scan at `3000 fb^-1` and uses
`N_signal >= 3` as its default threshold. Here `acceptance_i` is the GRENDEL
reconstruction + selection efficiency described next, and the lifetime
reweighting plus decay-in-volume probability `P_decay_i(U^2)` are applied in
`analysis/decay_reco_acceptance.py::scan_u2`.

### Signal Decay and Acceptance

The HNL's own decay is flavor dependent -- the visible final states differ for
`Ue`, `Umu`, and `Utau` -- so a fixed `e+e-` final state is wrong for the muon
and tau scenarios. Decays are sampled with the vendored FairShip HNL modules
(`vendored/fairship/`, driven through `ROOT.TPythia8`), which also supply the
flavor-aware lifetime `ctau(U^2 = 1)`. This is a separate template-generation
stage (`analysis/generate_decay_templates.py`) that caches rest-frame decays to
`tmp/decay_templates/<flavor>/templates_<mass>.npz`. It is the only stage that
needs PyROOT + Pythia8, and runs under the Homebrew-ROOT venv, not the conda
`hnl` env.

For each four-vector that crosses the fiducial volume, the analysis samples a
decay vertex along the flight path, boosts a FairShip template to the lab,
keeps the two highest-momentum charged stable daughters (best-two-track),
ray-casts them to the tracker walls, and runs the shared
`../higgs/reco_common` bounded 4-hit reconstruction followed by the PR #13
selection (gate / pointing / collinearity / timing). The surviving fraction is
`acceptance_i`. This is the same reconstruction and selection the cosmic-decay
background uses, so signal and background share one definition.

### Production Decay Kinematics

These are the parent-to-HNL production decays (the HNL's own decay is covered
above). Two-body decays use exact rest-frame kinematics followed by a Lorentz
boost.
The polarized two-body tau modes use a fixed unit longitudinal analyzing
power whose sign follows the tau origin: `+1` (HNL forward) for W-origin
taus, `-1` (HNL backward) for the helicity-suppressed heavy-meson leptonic
sources (`production/decay_engine/tau_decay.py`).

Meson three-body decays are sampled with
`decay_3body_weighted_dq2dE`, using the HNLCalc differential
`dBR/(dq2 dE)` expression. The b-baryon mode `Lambda_b -> Lambda_c l N` uses
`decay_3body_weighted_dq2dm122` with the HNLCalc `dBR/(dq2 dm12^2)`
expression. Leptonic tau three-body decays use
`decay_3body_weighted_dE` with the HNLCalc `dBR/dE` expression. These are not
flat phase-space samplers. The remaining approximation is the treatment of
angular and spin correlations that are not contained in those one- or
two-dimensional differential rates.

### Stored Weights

All production CSVs are headerless:

```text
weight, E, px, py, pz
```

`weight` is in pb and factors out the active mixing, so it corresponds to
`U_alpha^2 = 1`. Four-momenta are in GeV in the CMS lab frame.

Direct meson weights use

```text
2 * sigma_FONLL * fragmentation_fraction * BR(parent -> N + X) / N_sampled
```

where the factor of two converts the FONLL quark-plus-antiquark convention to
the total rate. `Bc` uses its separately configured cross section.
`sample_meson_4vectors` also supports an optional high-`pT` proposal mixture;
in that mode its `sampling_weight` is the exact nominal/proposal probability
ratio and must multiply the ordinary event weight. The default sampler remains
the nominal FONLL distribution with unit sampling weights.

Induced-tau weights additionally include the parent-to-tau branching ratios
(`parent -> tau nu` for the two-body modes, `B -> D(*) tau nu` for the
semitauonic ones). The b-baryon weight carries the bottom-fragmentation closure
remainder (`FRAG_LAMBDA_B`) and `BR(Lambda_b -> Lambda_c l N)`. Prompt-tau and
W/Z rows use their MadGraph event weights and the electroweak K-factor, which
is now keyed per process (`K_FACTOR_EW_BY_PROCESS` in
`production/constants.py`: the W value for `W/Z -> ell N`, and the W vs
Drell-Yan values per tau origin) rather than a single flat constant in the
drivers. The charged-kaon channel uses the Pythia 8.315 SoftQCD spectrum and
transport weight documented in `production/constants.py`.

## Setup

Create and activate the committed environment:

```bash
conda env create -f hnl/environment.yml
conda activate hnl
```

The direct runtime dependencies include `numpy`, `pandas`, `scipy`, `sympy`,
`mpmath`, `particle`, `matplotlib`, `tqdm`, `trimesh`, `rtree`, `numba`, and
`cycler`. The frozen campaigns pin `embreex==4.4.0` for deterministic,
accelerated ray intersections. `pytest` is included for tests. The full runner
checks these imports before starting expensive production.

MadGraph5_aMC@NLO v3.6.6 is required for prompt-tau and W/Z production. Its
executable is resolved in this order:

1. `$HNL_MG5_EXE`;
2. `hnl/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC`;
3. the shared project-level install;
4. the legacy sibling-checkout install.

LHAPDF (binary + library) is supplied by the conda environment. Provision the
PDF data set once into the env with `lhapdf install NNPDF40_nlo_as_01180`, or
point `$HNL_LHAPDF_DATA` at an existing LHAPDF data directory. Resolution is
conda-only (no sibling-checkout fallback); a missing PDF backend fails with a
precise error before MadGraph launches.

The signal decay-template stage (`analysis/generate_decay_templates.py`)
additionally needs PyROOT + Pythia8 (`ROOT.TPythia8`), which the conda `hnl`
env does not provide. Run it under a Python that has both; on this development
machine that is the venv at
`/Volumes/sandbox/projects/aaaPHYSICSaaa/.venvs/fairship` (Homebrew ROOT with
`--with-pythia`, plus `numpy`/`scipy`). Everything else -- production,
geometry, and the exclusion scan -- runs entirely in the conda env.

## Full Run

From the repository root:

```bash
conda activate hnl
cd hnl
./run_full_all.sh
```

The script runs non-W/Z production, W/Z production in independent mass shards,
channel combination, and analysis. FairShip template generation uses a
separate interpreter: set `HNL_TEMPLATE_PYTHON` to a Python with
`ROOT.TPythia8`, or pre-generate the templates as described under Manual
Operation. If the main `HNL_PYTHON` already provides `ROOT.TPythia8`, the
runner detects and uses it. Otherwise template generation is skipped and the
analysis consumes existing files; a clean run with no templates fails rather
than publishing an empty result.

By default it uses:

- 6 non-W/Z production workers;
- 12 cores for prompt-tau MadGraph generation;
- 8 simultaneous W/Z shards with one MG5 core each;
- 3 analysis workers.

Analysis workers load large combined CSVs and geometry arrays. On an 18 GB
development machine, 12 analysis workers caused severe memory pressure and a
userspace watchdog restart. Three workers completed reliably. This limit is
about memory, not CPU utilization.

The runner creates a unique tag such as `full_20260608_231500`, preventing a
new invocation from silently overwriting an older run. To resume a run, reuse
its tag:

```bash
HNL_RUN_TAG=full_20260608_231500 ./run_full_all.sh
```

W/Z production passes `--skip-existing`, so non-empty completed mass points
are retained. Other useful overrides are:

```bash
HNL_PYTHON=/path/to/python \
HNL_TEMPLATE_PYTHON=/path/to/root-enabled-python \
HNL_PRODUCTION_WORKERS=6 \
HNL_PROMPT_TAU_CORES=12 \
HNL_WZ_JOBS=8 \
HNL_WZ_CORES_PER_JOB=1 \
HNL_ANALYSIS_WORKERS=3 \
./run_full_all.sh
```

The target W/Z CPU use is `HNL_WZ_JOBS * HNL_WZ_CORES_PER_JOB`. Avoid setting
both values high enough to oversubscribe the machine.

Generated files are written below:

```text
tmp/runs/<tag>/llp_4vectors/
tmp/runs/<tag>/analysis/
tmp/cache/
tmp/wz_shards/
```

Top-level stage logs are written as `tmp/<tag>_*.log`.

## Manual Operation

### Production without W/Z

```bash
python -u run_all.py \
  --flavor Ue Umu Utau \
  --no-wz \
  --skip-combine \
  --workers 6 \
  --prompt-tau-nb-core 12
```

`run_all.py` parallelizes channel/flavor jobs. Its built-in W/Z path processes
the masses serially inside each flavor job; the full runner deliberately uses
the sharded W/Z driver instead.

### Sharded W/Z Production

```bash
python -u production/madgraph/run_wz_sharded.py \
  --flavor Umu \
  --jobs 8 \
  --nb-core 1 \
  --skip-existing
```

Run this once per flavor. Every shard has an isolated MadGraph work directory,
which avoids process-tree collisions.

### Combine

```bash
python -m production.combine_channels --flavor Ue Umu Utau
```

### Decay templates

Generate the FairShip rest-frame decay templates once per flavor (this is the
only stage that needs PyROOT + Pythia8, so use the ROOT venv, not conda):

```bash
/Volumes/sandbox/projects/aaaPHYSICSaaa/.venvs/fairship/bin/python \
  analysis/generate_decay_templates.py \
  --flavor Ue Umu Utau --n-templates 20000 --out tmp/decay_templates
```

The analysis reads `tmp/decay_templates/<flavor>/templates_<mass>.npz`; a mass
point with no template file is skipped with a note. `--skip-existing` reuses
already-generated files.

### Analysis

```bash
python -u run_analysis.py --flavor Ue Umu Utau --workers 3
```

The analysis writes:

```text
tmp/runs/<tag>/analysis/hnl_sensitivity.csv
tmp/runs/<tag>/analysis/hnl_exclusion.pdf
tmp/runs/<tag>/analysis/hnl_exclusion.png
tmp/runs/<tag>/analysis/scan_status.json
tmp/runs/<tag>/analysis/run_metadata.json
```

If an exclusion boundary is not crossed inside the configured `U^2` scan, the
CSV stores a `NaN` boundary and sets `u2_min_open` or `u2_max_open`. The plot
fills to the corresponding axis edge, omits a false closing line, and marks
the open direction.

`analysis.exclusion.find_exclusion_band_refined` can refine a bracketed island
against an already-built deterministic yield evaluator without resampling the
Monte Carlo. `analysis.decay_reco_acceptance.signal_contribution_diagnostics`
reports sample- and event-level effective statistics for the same frozen
weighted estimator. These shared helpers are used by the BC4/BC10 convergence
campaigns and do not alter the default HNL scan unless called explicitly.

Use `--plot-only` to regenerate the plot from an existing sensitivity CSV and
`--force-geometry` to rebuild cached ray intersections.

## Downstream / handoff (who consumes the curve)

This repository **only produces the GRENDEL curve.** The final comparison
figures -- GRENDEL overlaid on the PBC BC7 contours and the HNLimits community
compilation of existing exclusions + competitor projections -- are made in the
workspace repository `shared/curves_PBC`, which reads our curve as input.

The contract is a single CSV with (at minimum) the columns
`mass_GeV, flavor, u2_min, u2_max` (plus `has_sensitivity`, `u2_max_open`),
exactly the schema `run_sensitivity.py` writes. `u2_min`/`u2_max` are the lower
and upper mixing edges of the excluded island; a blank/`NaN` `u2_max` with
`u2_max_open = True` means the island is open upward.

Because each analysis run lands in a per-run, git-ignored
`tmp/runs/<tag>/analysis/hnl_sensitivity.csv` (the tag changes every run and the
directory is eventually cleaned up), the canonical curve is **published to a
stable, committed path**:

```text
data/published/grendel_hnl_sensitivity.csv   # the curve consumers read
data/published/MANIFEST.json                 # provenance: run, sha, cut, reach
data/published/README.md                     # how to re-publish
```

`shared/curves_PBC` points its default HNL input at this file. Override it with
`HNL_GRENDEL_CSV`; `CURVES_PBC_GRENDEL_CSV` remains a legacy alias. After
producing a better run, re-publish by
copying its `hnl_sensitivity.csv` into `data/published/` and refreshing the
manifest (see `data/published/README.md`); do **not** point the consumer at a
`tmp/runs/<tag>` path.

The current published curve uses the
`central_newgrids_20260623/analysis_exact_100_betafix`
run (`P > 100 MeV` track cut) with the low-mass BC6/BC7 charged-kaon channel
replaced by the 2026-07-18 Pythia/transport rerun. Note for the figure caption:
the high-mass island closes near `m_N ~ 3.6 GeV` because of the `ctau ~ 1/m_N^5`
lifetime law (peak yield `~ sigma * beta*gamma / m_N^5`), **not** a B-meson /
`|V_ub|^2` production cutoff -- `W/Z` produce `N` out to 10 GeV but those are too
short-lived to reach a 22 m displaced detector (that regime belongs to prompt /
near-IP displaced-vertex searches).

### Why the island closes at 3.6 GeV (quantitative)

The published run itself demonstrates the mechanism; all numbers below are
read directly from `data/published/`:

- **It is not missing production.** At the closure point the electroweak
  channel already dominates: `bundle/channel_breakdown_u2min.csv` gives the
  `WZ` channel **89.5% of the peak yield at 3.6 GeV (Umu)** (Bmeson 5.9%,
  Bc 4.5%). More W/Z Monte Carlo statistics cannot reopen the island. A
  genuinely missing production channel would be a model change and would need
  enough accepted yield to overcome the steep lifetime suppression discussed
  below.
- **The collapse is the m^-6 lifetime wall.** `peak_N` in
  `grendel_hnl_sensitivity.csv` is the yield at the *optimal* mixing, i.e.
  the best GRENDEL can do at that mass: 3.02 at 3.62 GeV (last point with
  `N >= 3`), 1.63 at 4 GeV, 0.40 at 5 GeV, 0.022 at 8 GeV, and 0.0056 at
  10 GeV. From 4 to 8 GeV both `peak_N` and `peak_u2` fall by about
  `70 ~= 2^6`: five powers of mass from
  `Gamma_N ~ G_F^2 U^2 m^5`, and one from the boost
  `beta*gamma ~ E/m` at the approximately mass-independent
  `E ~ m_W/2` of on-shell `W -> l N`.
- **Analytic form.** For a shell detector at distance `d` with fiducial depth
  `dL`, writing `lambda_1 = beta*gamma * ctau(U^2=1) ~ m^-6`, the scan
  `N(U^2) = sigma_1 U^2 L * (dL U^2/lambda_1) exp(-d U^2/lambda_1)` peaks at
  `U^2_opt = 2 lambda_1/d` with
  `N_max ~ sigma_W L dL lambda_1 / d^2`. At 10 GeV even the optimal coupling
  (`peak_u2 = 9.7e-10`) yields only about 60 produced HNLs in all of
  `3 ab^-1`, before acceptance.
- **Geometry enters only through a sixth root.** The closure mass scales as
  `(sigma L dL / d^2)^(1/6)`, so large geometric or rate changes move it
  modestly. Digitized projections in `shared/curves_PBC` show the same family
  pattern for Umu: CODEX-b about 3.0--3.3 GeV, GRENDEL 3.63--3.69 GeV,
  MATHUSLA 4.05 GeV, SHiP 5.1--5.8 GeV, and ANUBIS 6.8--7.7 GeV. Moving the
  GRENDEL closure from 3.63 to 10 GeV would require roughly
  `(10/3.63)^6 ~= 440` times more peak yield. The higher-mass on-shell W/Z
  regime therefore belongs to near-IP searches rather than a detector about
  20 m from the IP.

## Paths and Files

Path defaults are centralized in `production/paths.py`. Supported overrides
include:

```text
HNL_RUN_TAG
HNL_RUN_DIR
HNL_TMP_DIR
HNL_MG5_WORK_DIR
HNL_TAU_POOL_CSV
HNL_LLP_VECTORS_DIR
HNL_ANALYSIS_DIR
```

Mass filenames use three decimal places with `p` replacing the decimal point:
`1.025 GeV -> mN_1p025.csv`. Use
`config_mass_grid.parse_mass_from_filename` rather than parsing labels by
hand.

Zero-byte channel files are intentional threshold sentinels. They mean that
the channel is closed at that mass. Non-empty files near a threshold can have
fewer rows than the requested pool because only a subset of parent species is
kinematically allowed.

## Tests

From the repository root:

```bash
conda activate hnl
python -P -m pytest hnl/tests/ -q
```

The tests cover FONLL parsing, nominal and importance sampling, mass labels,
channel combination,
two- and three-body kinematics, tau and kaon production, W/Z CSV conversion,
open and refined contour handling, weighted-signal effective statistics, the
FairShip->GRENDEL acceptance core (boost, best-two tracks,
reconstruction/selection, and the interior lifetime peak of the `U^2` scan),
and a meson-production smoke path. `-P` prevents any
stale ignored MG5 parser cache named `py.py` from shadowing pytest's
compatibility module.

## Physics Inputs and Limitations

- The heavy-meson backend is the committed central NNPDF4.0 NLO FONLL grid in
  `data/production/fonll/central/`, generated in the external workspace
  `/Volumes/sandbox/projects/aaaPHYSICSaaa/shared/NNPDF40/fonll-local`.
- The FONLL tables stop at `pT = 50 GeV`. The central curve uses one central
  scale/PDF choice; scale/PDF/mass variation grids are produced in that external
  FONLL workspace and propagated into a band via `run_variation_band.py`
  (`--grid-dir`, default `tmp/fonll/output`) + `analysis/combine_band.py`.
  `combine_band.py` can optionally fold a PDF4LHC alpha_s term from companion
  curves (`--alphas-lo/--alphas-hi`), but those companion analysis chains are
  not part of the published bundle. The post-`beta=p/E` exact-hit campaign
  propagated 111 coherent curves (central plus 110 non-central scale, PDF, and
  heavy-quark-mass members) with no topology changes among the non-central
  members.
  Its median combined half-widths are `-0.101/+0.131` dex on the lower edge
  (scale dominated) and `-0.0066/+0.0081` dex on finite upper edges.  The largest
  upper shift, 0.114 dex, survives independent exact-200 controls; the previous
  order-one upper band was an artifact of a 4,000-hit resampling cap.  The full
  raw table, manifests, controls, and combined result are tracked in
  `data/published/bundle/`.  `analysis/plot_money.py` renders these as a separate
  theory/model diagnostic; the paper comparison remains central-only.
- Charm and bottom species share one heavy-flavor shape per table; species
  fractions are applied in the event weights. `Bc` reuses the bottom shape at
  the Bc mass. Per-species and dedicated-Bc shapes are planned
  (see `REMAINING_WORK.md`).
- b-baryons are now included through `Lambda_b -> Lambda_c l N`, which carries
  the closure remainder of the LHCb bottom-fragmentation ratios;
  `Xi_b`/`Omega_b` are not separated and are approximated as Lambda_b-like.
  The Lambda_b pT-y shape reuses the bottom FONLL grid at the Lambda_b mass.
- Bottom and charm fragmentation inputs are measured in restricted,
  different acceptances and are extrapolated as constants over the FONLL
  phase space; their uncertainties are not propagated
  (see `REMAINING_WORK.md`).
- The charged-kaon flux uses a Pythia 8.315 SoftQCD spectrum and a
  decay-before-absorption transport weight (see `REMAINING_WORK.md`); it
  dominates the normalization uncertainty at the lightest masses. The BC6/BC7
  low-mass curves were rerun and republished with this model (2026-07-18,
  `d_esc = 1.5 m`; the transport band is `data/published/bundle/kaon_desc_band.csv`).
  The transport escape length remains a proxy pending a real material map.
- Meson, baryon, and tau three-body energy/`q^2` distributions are
  HNLCalc-weighted, but complete multidimensional matrix-element spin
  correlations are not modeled (see `REMAINING_WORK.md`).
- The direct-HNL rates use HNLCalc's own form factors, CKM elements, and decay
  constants -- community-standard parameterizations (some are proxies, e.g.
  several channels reuse the pion form factor; CKM/decay constants are hardcoded
  without uncertainties), inherited by citing HNLCalc rather than independently
  revalidated -- adequate for an exclusion contour but not a precision input.
  Induced-tau modes use externally normalized branching fractions; direct-HNL
  channels retain HNLCalc's absolute normalization. The vendored fork also
  corrects four upstream HNLCalc bugs, so it deviates from the citable release
  until those are upstreamed (see `vendored/PROVENANCE.md`, `REMAINING_WORK.md`).
- The HNL total-width / lifetime quark-hadron duality is propagated as a separate
  seam-derived decay-model band `delta(m)` (`analysis/width_band.py` +
  `analysis/decay_model_band.py`), driven coherently through `ctau` and
  `vis_frac` (the composition leg self-cancels to ~1%); stable finite upper-edge
  shifts have median half-widths near `-0.094/+0.115` dex.  One nuisance
  direction removes the island near high-mass closure, so those endpoints are
  recorded as topology changes rather than drawn as a smooth ribbon. The
  *absolute* visible-BR normalization is a distinct,
  still-unpropagated decay nuisance (see `REMAINING_WORK.md` item 15).
- The electroweak K-factor is a per-process table whose entries all currently
  hold the inclusive `1.3` constant; differential NLO/LO values are an optional
  upgrade (see `REMAINING_WORK.md`).
- The signal acceptance is the GRENDEL 4-hit reconstruction + PR #13 selection
  efficiency (best-two-track) applied to the FairShip decays, the same
  selection the cosmic-decay background uses. The exclusion itself is still a
  zero-background `N_signal >= 3` threshold: the cosmic background is computed
  by the shared `../higgs/` machinery but is not yet folded into the limit, and
  no systematic uncertainties are included (see `REMAINING_WORK.md`).

The FONLL and HNLCalc source details, versions, and local modifications are
recorded in `vendored/PROVENANCE.md`. Physics, detector, numerical, and
reproducibility work still needed is tracked in `REMAINING_WORK.md`.
