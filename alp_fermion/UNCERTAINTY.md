# BC10 theory/model variation diagnostics

The publication campaign propagates every nuisance through a complete signal
run. No event-level importance reweighting is used. The compact result is a
non-statistical, one-source-at-a-time variation envelope, not a confidence
interval or a simultaneous-combination coverage statement. It is retained as
an audit diagnostic and is not drawn on the primary comparison plots.

## Inputs and combination

- **FONLL scales:** independent 600,000-parent runs for the six non-central
  coherent 7-point `(mu_R, mu_F)` grids. The contour uncertainty is their
  asymmetric envelope around the independently regenerated central run.
- **FONLL PDF:** independent 600,000-parent runs for all 100 NNPDF4.0 NLO
  replicas. Their 16th and 84th percentiles in `log10(1/f)` define the PDF
  interval; the sample standard deviation is retained only as an audit field.
- **FONLL bottom mass:** independent 600,000-parent runs at `m_b = 4.5` and
  `5.0 GeV`, around the `4.75 GeV` central grid. The maximum absolute contour
  displacement is used.
- **Two-gluon decay surrogate:** full-branching Pythia samples in which the
  imported `a -> gg` branching fraction is represented by pure `u ubar`,
  `d dbar`, or `s sbar`. The central templates use an equal `u/d/s` mixture.
  Each alternative is propagated through its own geometry, reconstruction,
  and sensitivity run. Widths and branching fractions remain fixed to the
  decoded arXiv:2501.04525/SensCalc v1.3.3 tables.
- **`C_bs` scheme:** `+/-20%` amplitude variations around the ALPINIST
  one-loop/RG value, implemented as full production runs with rates scaled by
  `0.8^2` and `1.2^2`, followed by independent geometry/reconstruction scans.
- **Heavy-pseudoscalar decay structure:** the exact arXiv:2310.03524 widths,
  all 32 exclusive branching channels, and all 18 squared matrix elements
  shipped in the same pinned SensCalc v1.3.3 release are propagated through
  independent 20,000-event Pythia templates and a complete reconstruction
  scan. Production physics is unchanged, so this run reuses the central
  production vectors exactly. This is a one-sided named model comparison, not
  a calibrated uncertainty. It is published as a separate named contour data
  product and does **not** enter the pointwise envelope.
- **Numerical controls:** two same-physics central repeats use fresh 600,000-
  parent pools and distinct production and reconstruction RNG seeds with the
  central templates. They are excluded from the theory/model envelope.

All shifts are evaluated in `log10(1/f)`. Scale, bottom-mass,
gluon-surrogate, and `C_bs` sources retain their named extrema around the
central contour. The overall display envelope is the outermost boundary among
those four intervals and the PDF 16th/84th-percentile interval, with only one
source varied at a time. Sources are not added in quadrature. If a variation
removes an island boundary, `*_variation_missing` is set rather than treating
the missing boundary as zero displacement. Central insensitive rows and the
excluded eta/eta-prime pole rows remain explicit gaps with NaN envelope edges.
Named physical variations that restore or remove sensitivity, or change an
open-edge state, are listed separately by mass rather than converted into a
finite displacement around a missing boundary.
The structural table separately records whether the 2310 model restores or
removes sensitivity and whether an open-boundary state changes. In particular,
central gaps restored by the alternative remain explicit topology changes;
they are never converted into an envelope displacement around a nonexistent edge.
The compact table also reports each repeat boundary, median/max absolute repeat
shift in dex and fractional coupling, and flags a mass when the largest repeat
shift is at least as large as the physical envelope shift on that boundary.
Sensitivity and open-edge topology disagreements are recorded independently,
including the responsible repeat names and masses, even where the central
contour has no boundary to compare. These flags never restore a gap in the
central contour and never enter the physical envelope.

## Run layout and restart policy

`run_uncertainty_campaign.py` writes heavy artifacts only below the explicit
`--scratch-root` (or `$ALP_UNCERTAINTY_SCRATCH`) path. Each variation owns a
run directory containing:

```text
runs/<variation>/
  llp_4vectors/            # fresh 600k production output while running
  analysis/                # checkpoint CSV, geometry while running, plot
  production.log
  sensitivity.log
  production.complete.json
  variation.complete.json
```

Completion markers are written atomically only after a stage validates its
expected outputs. `alp_production --resume` and `sensitivity --resume` make an
interrupted variation restartable. The campaign can be partitioned over
several processes with `--worker-index` and `--worker-count`; the assignment is
stable because it follows the manifest order.

After a non-central variation validates completely, the completion marker is
written with vector and geometry tree hashes, exact grid/seed/code/template
provenance, commands, and log hashes. The runner then removes that variation's
raw vectors and geometry cache and atomically records `storage_state=compacted`.
An interruption during reclamation resumes from the `compacting` state. The
central vectors and geometry are always retained because the three gluon
surrogate variations and the 2310 structural comparison reuse central
production exactly. Pass
`--keep-intermediates` to suppress reclamation for new full runs.

Template bundles are shared through `$ALP_TEMPLATE_DIR`; run trees and source
directories are selected through environment variables, not repository
symlinks. The raw contour table, compact
`bc10_single_source_variation_envelope.csv`, dedicated
`bc10_decay_2310_structural_alternative.csv`, and exact registry/provenance
manifest are promoted to `alp_fermion/data/published/bundle/`. The definitive
refined-grid comparison additionally includes the byte-identical raw contour
`bc10_decay_2310_structural_curve_dense.csv`, its derived central-vs-structural
table, and `DENSE_STRUCTURAL_MANIFEST.json`, so the bundle is self-contained.
The manifest requires one consistent code state for the pointwise campaign and
one for the structural add-on, and records both. This permits the structural
implementation commit to be applied after already validated pointwise runs without
weakening reproducibility within either group.

## Reproduction

Generate the three pure-flavor alternatives in the ROOT/Pythia environment:

```bash
PYROOT=/path/to/fairship/bin/python
for spec in u:2234 d:3234 s:4234; do
  q=${spec%:*}; seed=${spec#*:}
  $PYROOT -m alp_fermion.generate_decay_templates_pythia \
    --mass 2.20 2.25 2.30 2.35 2.40 2.45 2.50 2.55 2.60 2.65 \
           2.70 2.75 2.80 2.85 2.90 2.95 3.00 3.10 3.20 \
           3.30 3.40 3.50 3.60 3.70 3.80 3.90 4.00 4.10 4.20 \
           4.30 4.40 4.50 4.60 4.65 4.70 4.75 \
    --n-templates 20000 --seed "$seed" --gluon-surrogate "$q" --resume \
    --out /scratch/bc10_uncertainty/templates/gg_$q
done

$PYROOT -m alp_fermion.generate_decay_templates_pythia \
  --decay-model 2310_structural --mass-grid-file /scratch/final_mass_grid.csv \
  --n-templates 20000 --seed 5234 --resume \
  --out /scratch/bc10_uncertainty/templates/decay_2310_structural
```

Start one or more FONLL workers, then the three gluon variants, two `C_bs`
variants, and the structural comparison. Run numerical controls after the
physics variations:

```bash
PY=/path/to/llpatcolliders_FONLL/bin/python
for worker in 0 1 2; do
  $PY -m alp_fermion.run_uncertainty_campaign \
    --scratch-root /scratch/bc10_uncertainty \
    --grid-dir /path/to/fonll-nnpdf40/output \
    --mass-grid-file /scratch/final_mass_grid.csv \
    --worker-index "$worker" --worker-count 3 --resume
done

$PY -m alp_fermion.run_uncertainty_campaign \
  --scratch-root /scratch/bc10_uncertainty \
  --grid-dir /path/to/fonll-nnpdf40/output \
  --mass-grid-file /scratch/final_mass_grid.csv \
  --axes decay_gg cbs --resume

$PY -m alp_fermion.run_uncertainty_campaign \
  --scratch-root /scratch/bc10_uncertainty \
  --grid-dir /path/to/fonll-nnpdf40/output \
  --mass-grid-file /scratch/final_mass_grid.csv \
  --axes decay_structure \
  --structural-template-dir \
    /scratch/bc10_uncertainty/templates/decay_2310_structural --resume

$PY -m alp_fermion.run_uncertainty_campaign \
  --scratch-root /scratch/bc10_uncertainty \
  --grid-dir /path/to/fonll-nnpdf40/output \
  --mass-grid-file /scratch/final_mass_grid.csv \
  --axes numerical_control --resume

# Rerun the exact 2310 contour directly on the refined high-statistics central
# grid, then pin its vectors, geometry, templates, code, and output hashes.
ALP_LLP_VECTORS_DIR=/scratch/bc10_final/llp_4vectors \
ALP_TEMPLATE_DIR=/scratch/bc10_final/decay_templates_2310 \
ALP_ANALYSIS_DIR=/scratch/bc10_final/analysis_2310 \
ALP_GEOM_CACHE_DIR=/scratch/bc10_final/geometry_cache \
$PY -m alp_fermion.sensitivity \
  --output /scratch/bc10_final/analysis_2310/bc10_sensitivity.csv \
  --decay-samples 60 --resume

RUN_COMMIT=0123456789abcdef0123456789abcdef01234567
RUN_COMMAND='ALP_* paths as above; python -m alp_fermion.sensitivity --output ... --decay-samples 60 --resume'
$PY -m alp_fermion.record_dense_structural \
  --central-curve /scratch/bc10_final/analysis/bc10_sensitivity.csv \
  --central-manifest alp_fermion/data/published/MANIFEST.json \
  --structural-curve /scratch/bc10_final/analysis_2310/bc10_sensitivity.csv \
  --mass-grid /scratch/bc10_final/mass_grid.csv \
  --vector-dir /scratch/bc10_final/llp_4vectors \
  --geometry-dir /scratch/bc10_final/geometry_cache \
  --template-dir /scratch/bc10_final/decay_templates_2310 \
  --producer-git-sha "$RUN_COMMIT" \
  --command "$RUN_COMMAND" \
  --output /scratch/bc10_final/analysis_2310/DENSE_STRUCTURAL_MANIFEST.json

$PY -m alp_fermion.combine_uncertainty_band \
  --scratch-root /scratch/bc10_uncertainty \
  --grid-dir /path/to/fonll-nnpdf40/output \
  --dense-structural-curve \
    /scratch/bc10_final/analysis_2310/bc10_sensitivity.csv \
  --dense-structural-manifest \
    /scratch/bc10_final/analysis_2310/DENSE_STRUCTURAL_MANIFEST.json
```

`final_mass_grid.csv` must contain a strictly increasing `mass_GeV` column.
Replace `RUN_COMMIT` and `RUN_COMMAND` with the exact committed code state and
fully expanded launch command used for the direct refined-grid run.
The runner records both its file hash and a canonical hash of the parsed mass
list in every production and variation marker, and refuses to reuse central
vectors generated on a different grid. Omit `--mass-grid-file` only when
reproducing the package's built-in grid exactly.
