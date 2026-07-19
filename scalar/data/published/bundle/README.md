# GRENDEL BC4 uncertainty bundle

Version-controlled theory uncertainty inputs for the canonical BC4 sensitivity
curve in `../bc4_island.csv`.

Files:

- `bc4_uncertainty_variations.csv` — all 110 independently simulated physics
  curves plus two same-physics numerical-control repeats.
- `bc4_single_source_variation_envelope.csv` — canonical central contour,
  per-source intervals, their outer display envelope, open-edge flags, source
  names, and audit shifts in dex for both contour edges.
- `UNCERTAINTY_MANIFEST.json` — pinned FONLL grid checksums, campaign
  configuration, combination prescription, output checksum, and limitations.

The FONLL contribution uses all 109 coherent bottom grids from the NNPDF4.0 NLO
campaign: central, the six non-central points of the standard seven-point scale
set, 100 Monte Carlo PDF replicas, and `m_b = 4.5, 5.0 GeV`. Every variation
uses a fresh 200,000-event parent pool, producing 600,000 scalar events per
mass after the three B-species contributions. Scalar kinematics, geometry,
decay/reconstruction MC, and the sensitivity scan are all independent between
variations; no importance reweighting or detector-outcome reuse is used.

Two additional central-FONLL/Winkler runs use fresh production and
reconstruction seeds. They quantify finite-simulation scatter only and are
excluded from the scale, PDF, bottom-mass, decay-model, and headline envelopes.
For each boundary the compact table reports the repeat median/max absolute dex
shift, maximum fractional shift, ratios to every physical source displacement,
and a flag when the repeat maximum is not smaller than the largest physical
source displacement at that mass. The manifest summarizes all flagged masses.
Sensitivity and open-edge topology disagreements are recorded separately
against both the independent campaign central and the canonical published
curve, including the responsible repeat names and masses. They remain
numerical diagnostics only and cannot restore a gap or enter an envelope.

Combination is in `log10(sin^2 theta)`: the scale source uses the extrema of
the coherent seven-point set, the PDF source uses the 16th and 84th percentiles
of the 100 replicas (with the sample standard deviation retained for audit),
and the bottom-mass source uses the extrema of the central, 4.5, and 5.0 GeV
curves.
The decay-model contribution comes from a separate fresh central-FONLL
production and full detector simulation using analytic LO-ChPT widths below
2 GeV plus spectator widths above 2 GeV, compared with the matched/dispersive
Winkler central run. It is an alternate-model envelope, not a Gaussian error.
The displayed `single_source_variation_envelope` is the outermost boundary of
these one-source-at-a-time intervals, rebased onto the canonical central curve.
Nothing is added in quadrature; it is not a confidence band and does not claim
simultaneous-source coverage.
Physical variations that restore or remove sensitivity, or change an open-edge
state, are listed by variation and mass against both the independent campaign
central and the canonical contour. They are not converted into a finite band
width where the corresponding reference boundary does not exist.

Reproduce from the repository root in the `llpatcolliders_FONLL` environment:

```bash
python -m scalar.uncertainty_band run --workers 2 \
  --grid-dir /Volumes/sandbox/projects/aaaPHYSICSaaa/shared/NNPDF40/fonll-local/output \
  --scratch-dir /Volumes/GRENDEL/extra_space/bc4_uncertainty
python -m scalar.uncertainty_band status
python -m scalar.uncertainty_band collect
```

The heavy resumable run tree stays in external scratch and is not published.
Each vector, geometry cache, result, and aggregate curve is written atomically
and pinned by checksums. After each variation is validated, raw vectors and
geometry are reclaimed; per-stage tree hashes, compact results, seeds,
provenance, and logs remain. The campaign requires `embreex==4.4.0`; 25,000
actual scalar rays were cross-checked against the triangle backend with zero
hit-mask mismatches and entry/exit differences below `2.2e-14 m`. This envelope
does not include FONLL alpha-s companions, detector-response systematics, or
background uncertainty.
