# GRENDEL BC10 variation-envelope bundle

Version-controlled theory/model variation inputs for the canonical BC10
sensitivity curve in `../bc10_sensitivity.csv`.

Files:

- `bc10_uncertainty_variations.csv` - all 117 exact full-run contours: central,
  six scale choices, 100 NNPDF replicas, two bottom-mass choices, three gluon
  decay surrogates, two `C_bs` scheme choices, one exact 2310 decay-structure
  comparison, and two same-physics numerical controls.
- `bc10_single_source_variation_envelope.csv` - the 99-row plot input with the
  central contour, one-source-at-a-time display envelope, open/missing flags,
  named source extrema, PDF percentiles, and PDF log-coupling standard
  deviation for audit. It also contains the repeat-control boundary shifts and
  numerical-spread flags.
- `bc10_decay_2310_structural_alternative.csv` - the matched 99-point campaign
  audit of the one-sided structural contour, including explicit
  restored/removed sensitivity and open-boundary topology flags.
- `bc10_decay_2310_structural_alternative_dense.csv` - the definitive
  145-point structural-topology diagnostic, rerun directly with the refined
  canonical grid and high-statistics central production vectors. It is the
  dashed-contour input and is excluded from the pointwise halo.
- `DENSE_STRUCTURAL_MANIFEST.json` - hashes, code state, templates, production
  vectors, and command provenance for the refined structural run.
- `UNCERTAINTY_MANIFEST.json` - exact variation registry, input and output
  checksums, code/template provenance, combination prescription, and headline
  shifts.

The display object is explicitly a
`single_source_variation_envelope`, not a confidence interval. Scale,
bottom-mass, `C_bs`, and `a -> gg` surrogate alternatives use their named
pointwise extrema. The 100 NNPDF replicas use their 16th and 84th percentiles
in `log10(1/f)`; their raw extrema do not define the display envelope. The
overall edge is the outermost of these source intervals with one source varied
at a time. Sources are not added in quadrature. Both 2310 structural files are
named model comparisons, not interval endpoints, and never enter this envelope;
use the dense file for topology interpretation and the matched-grid file only
to audit the 117-run campaign.

Two same-physics central repeats use fresh 600,000-parent pools, independent
production seeds, distinct reconstruction RNG offsets, and the central
templates. They do not enter the physical envelope. For each boundary the
bundle reports their median and maximum absolute shifts in dex and fractional
coupling, and flags masses where the maximum repeat shift is not smaller than
the physical one-source envelope shift.

Every FONLL and `C_bs` variation uses a fresh 600,000-parent production pool
and full geometry, reconstruction, and sensitivity scan. The `u`, `d`, and
`s` gluon-surrogate templates and the exact 2310 structural templates use
central production-vector reuse because their production physics is identical,
followed by independent downstream simulation. No importance reweighting or
detector-outcome reuse is used.

Central insensitive rows and the excluded eta/eta-prime pole windows remain
explicit NaN gaps. Heavy resumable run trees stay in external scratch. After a
non-central run validates and records its raw vector/geometry hashes, those
regenerable trees are compacted; the sensitivity CSV, markers, provenance, and
logs remain. Central raw output is retained for exact gluon-variant reuse.

Reproduce in the `llpatcolliders_FONLL` environment using the commands in
`alp_fermion/UNCERTAINTY.md`, then collect with:

```bash
python -m alp_fermion.combine_uncertainty_band \
  --scratch-root /scratch/bc10_uncertainty \
  --grid-dir /path/to/fonll-nnpdf40/output \
  --dense-structural-curve /scratch/bc10_2310_dense.csv \
  --dense-structural-manifest /scratch/bc10_2310_dense.manifest.json
```
