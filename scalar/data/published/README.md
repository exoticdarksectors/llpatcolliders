# Published GRENDEL BC4 curve (canonical handoff)

This is the **stable, canonical** GRENDEL BC4 (Higgs-portal dark scalar)
sensitivity curve, following the same publish pattern as
`hnl/data/published/`: a *copy* of one chosen analysis run promoted to a
fixed, committed path so consumers (the paper, comparison plots) never read a
git-ignored `tmp/` file.

## Files

- `bc4_island.csv` — the curve. Same schema as the analysis output
  (`scalar/run_sensitivity.py`); the columns a consumer requires are
  `mass_GeV, u2_min, u2_max` (plus `has_sensitivity`), where `u2` is
  `sin^2 theta`.
- `MANIFEST.json` — provenance: source run, git sha, csv sha256, inputs,
  luminosity, threshold, headline reach, and the closure note.

## To re-publish (after a better run)

```bash
cp scalar/tmp/bc4_island.csv scalar/data/published/bc4_island.csv
# then refresh MANIFEST.json (sha, date, csv_sha256, headline_reach)
```

## Current contents

- Source (Lambda_b republication, 2026-07-18): the central island uses a
  1.6M-scalar high-pT importance production pool per mass (400000 parents x
  B+, B0, Bs, Lambda_b), 100 decay/reconstruction samples per entering scalar,
  and the timing-chi2 daughter `beta = p/E` correction. The uncertainty bundle
  is a matching 109-variation + 2 numerical-control campaign at 100000 parents
  per species; central and bundle are hash-linked.
- Reach: the sensitive grid spans `m_S = 0.14-3.80 GeV`; the upper boundary is
  open above the scan ceiling for the electron-only `0.14-0.20 GeV` extension.
  The deepest lower edge is `sin^2 theta = 6.5193e-12` at 0.975 GeV, near the
  f0(980) enhancement. Log-yield interpolation places the final closure at
  approximately 3.825 GeV. This supersedes the 2026-07-15 meson-only curve
  (deepest `7.3988e-12`, closure `3.798`): adding the b-baryon pool (`Lambda_b`,
  physical `Lambda` recoil) strengthens the lower edge by 8-18% (growing toward
  closure). The lower edge is fully physics-dominated; on the noisier upper edge
  the 100k-parent bundle leaves the numerical control marginal at 2.0 and
  2.3 GeV. The six-million-event closure control and the paper figure remain the
  meson-only baseline (small pending refreshes; see `MANIFEST.json`).

## Numerical convergence

The importance proposal gives at least 2249 effective events on every finite
sensitive lower edge and at least 2049 at every sensitive peak. Independent
six-million-event controls at 3.75, 3.80, and 3.85 GeV give peak yields of
3.925, 2.858, and 2.069, respectively; their interpolated closure is 3.792 GeV.
Exact values, seeds, hashes, and the production configuration are pinned in
`MANIFEST.json`.
