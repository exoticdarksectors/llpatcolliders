# GRENDEL HNL diagnostic bundle

Version-controlled inputs for the HNL production/decay uncertainty diagnostic.
These variations are **not statistical confidence bands** and are not drawn on
the central-only paper comparison figure.  `analysis/plot_money.py` renders them
as a separate diagnostic from a clean clone, with no `tmp/` run tree required.

## Contents

- `../grendel_hnl_sensitivity.csv` -- canonical beta-fixed central contour with
  the 2026-07-18 charged-kaon low-mass republication.
- `hnl_band_fonll.csv` -- combined FONLL scale/PDF/heavy-quark-mass result.
- `hnl_band_fonll_raw.csv` -- all 111 x 54 = 5,994 variation results.
- `FONLL_MANIFEST.json` -- base-campaign configuration and per-variation CSV
  hashes plus the structured 0.305 GeV partial-recompute contract.
- `FONLL_GRID_MANIFEST.json` -- sanitized 218-grid source manifest; original
  source-manifest SHA-256 is embedded.
- `decay_model_band.csv` -- coherent total-width/lifetime and visible-fraction
  variation, including dense closure anchors.
- `bc_nuisance_band.csv` -- direct-Bc +/-40% normalization diagnostic.  Each row
  records whether it came from the exact-50 broad scan, exact-200 closure scan,
  or exact-400 endpoint controls.
- `channel_breakdown_u2min.csv` -- exact-hit, eight-channel composition at the
  lower edge, peak, and finite upper edge.
- `kaon_desc_band.csv` -- charged-kaon `d_esc=1/1.5/3 m` transport diagnostic
  on the low-mass Ue/Umu lower edge.
- `NUMERICAL_CONTROLS.json`, `fonll_outlier_controls.csv`, and
  `bc_nuisance_endpoint_controls.csv` -- convergence and independent-seed checks.
- `MANIFEST.json` -- bundle hashes, schemas, summary statistics, and topology.

## FONLL campaign

The 111 coherent variations (central + six scale + 100 NNPDF4.0 replicas + four
`m_b`/`m_c` variations) were rerun after the timing `beta=p/E` correction using
all detector hits, 50 decay samples per hit, `event_chunk=1000`, two analysis
workers, and 100,000 generated events per production channel.  Every variation
has all 54 requested anchors.  All 110 non-central members agree with central on
finite/open topology at every lower and upper boundary.

The tracked raw table is hybrid only at `m_N=0.305 GeV`: its Ue/Umu rows were
recomputed for all 111 members after the Pythia-plus-transport kaon publication;
the other 5,772 rows retain the original campaign.  At this 99.85%-kaon-dominated
anchor the apparent member spread is acceptance-MC noise, so the combined FONLL
lower-edge ribbon is deliberately blank.  The raw rows and dex components remain
auditable, while `kaon_desc_band.csv` carries the relevant low-mass systematic.

In `log10(U^2)`, the combined FONLL lower-edge diagnostic has median
`-0.101/+0.131` dex half-widths (scale dominated; largest upward shift 0.301
dex).  Finite upper edges are much smaller: median `-0.0066/+0.0081` dex, with a
validated maximum upward shift of 0.114 dex.  The previous order-one upper-edge
structure came from the 4,000-hit weighted-resampling cap and is retired.

The decay-model variation has median stable finite-upper shifts of roughly
`-0.094/+0.115` dex.  At the high-mass closure, one nuisance direction removes
the island; those points are recorded as topology changes, not converted into a
smooth ribbon.  The Bc diagnostic similarly records three controlled topology
changes: Bc-up creates narrow Ue/Umu islands at 3.63 GeV, while Bc-down removes
the narrow central Umu island at 3.62 GeV.

## Numerical controls

- Exact-50 versus the published exact-100 central anchors differs by at most
  0.0444 dex on the lower edge and 0.0331 dex on finite upper edges.
- Chunked versus unchunked exact evaluation differs by at most 0.00537 dex
  (lower) and 0.00103 dex (upper).
- One-worker and two-worker results are bit-for-bit identical at the tested
  low-mass point.
- Independent exact-200 repeats validate the two largest FONLL upper-edge
  shifts to within 0.0032 dex.
- Three independent exact-400 Bc endpoint runs agree on topology; finite
  boundary spreads are at most 4.5% at the threshold-sensitive points.

Full definitions and source hashes are in `NUMERICAL_CONTROLS.json`.

## Reproduce the diagnostic

From `hnl/` in the `hnl` conda environment:

```bash
python -m analysis.plot_money --out-dir tmp/money_repro
```

This writes `money_plot.pdf`, `money_plot.png`, all plotted input CSVs (including
the kaon transport band), and `run_metadata.json`.  Ribbons are interpolated only within contiguous finite
anchor segments, so the renderer never extrapolates through a topology change.
