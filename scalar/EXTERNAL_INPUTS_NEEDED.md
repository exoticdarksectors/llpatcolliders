# External inputs -- BC4 dark scalar (`scalar/`)

_Reviewed 2026-07-17. The central decay input and the publication comparison
curves are both in place. This file records their provenance and scope._

## 1. Winkler dispersive hadronic widths -- satisfied

`scalar/data/winkler_widths.csv` contains the digitized dispersive result from
Winkler, arXiv:1809.01876, Fig. 4. `model.partial_widths` reads it below the
2 GeV handover and uses the perturbative spectator calculation above. The
table stores `Gamma_xx / sin^2(theta)` for the hadronic channels; provenance is
retained with the data.

The former analytic LO-ChPT calculation remains available as the
`chpt_spectator` alternate scheme. `scalar/uncertainty_band.py` propagates it
through a separate full simulation as a named decay-model diagnostic. It is
not used for the canonical curve and is not interpreted as a confidence
interval.

## 2. Existing bounds and proposed-experiment curves -- satisfied centrally

The authoritative comparison renderer and external data live in
`shared/curves_PBC`, not in `scalar/data/competitors/`. The BC4 publication
plot currently includes:

- the existing excluded region assembled from the documented PBC sources;
- the official SHiP projection from arXiv:2504.06692;
- the 300/fb, mixing-only CODEX-b projection from arXiv:1911.00481;
- the FASER2 curve from the 2025 PBC report, whose production is B-driven and
  applicable to BC4.

MATHUSLA and ANUBIS are deliberately omitted because no current-geometry,
mixing-only curves suitable for this comparison were identified. By user
policy, no ATLAS, CMS, or LHCb *projection* curves are included. Existing LHCb
constraints are not projections and remain part of the excluded landscape.

`scalar/plot_exclusion.py` retains a legacy optional CSV loader for local
diagnostics, but it is not the paper-figure source. Source files,
transformations, and citations for the actual comparison are documented under
`shared/curves_PBC/scalar/data/reference_curves/` and in that repository's
README.
