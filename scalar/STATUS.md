# BC4 (light dark scalar) -- publication status

_Last reviewed 2026-07-17. Active worktree:
`signal-models/bc4-scalar/llpatcolliders_BC4_PR17`, branch `bc4-scalar`,
package `scalar/`._

## Canonical result

The stable curve is `scalar/data/published/bc4_island.csv`; its exact source,
configuration, checksum, topology, and convergence controls are pinned in the
adjacent `MANIFEST.json` and described in `data/published/README.md`.

- 86-point grid, `m_S = 0.14--4.70 GeV`.
- Sensitive grid run `0.14--3.80 GeV` (Lambda_b, 2026-07-18); log-yield
  interpolation places the final closure at `~3.825 GeV`.
- Deepest lower edge `sin^2(theta) = 6.5193e-12` at `m_S = 0.975 GeV`.
- The 2026-07-18 curve adds the b-baryon pool (`Lambda_b`, physical `Lambda`
  recoil) to `b -> X_s S`; it strengthens the lower edge by 8-18% over the
  meson-only curve (deepest was `7.3988e-12`, closure `3.7975`). The central and
  its 109-variation uncertainty bundle are regenerated together and hash-linked.
- In the electron-only `0.14--0.20 GeV` interval, the upper edge lies above
  the configured coupling scan and is stored as open, not artificially closed.
- Central production uses 400,000 importance-sampled FONLL parents for each of
  `B+`, `B0`, `Bs`, and `Lambda_b` (1.6 million scalar events per mass), with a
  5 GeV high-pT tilt, a 50% nominal mixture, and 100 decay/reconstruction samples
  per scalar entering the detector. The mesons recoil against a kaon; `Lambda_b`
  recoils against a physical `Lambda`.
- The four-hit timing calculation uses each daughter track's true
  `beta = p/E`.

Independent six-million-event **meson-only** controls at 3.75, 3.80, and 3.85 GeV
reproduce the meson-only closure (3.792 GeV) within 0.006 GeV, so that endpoint is
not a finite-pool artifact. The published `Lambda_b` closure (3.825 GeV) instead
rests on the at-least-2049 effective events at every sensitive peak; a matching
`Lambda_b` six-million-event control is a small pending refresh (see
`MANIFEST.json`).

## Physics definition

BC4 is the minimal Higgs-mixing benchmark: one parameter, `sin^2(theta)`,
controls production, lifetime, and all decay branching fractions. The signal
criterion is `N_signal >= 3` at `3000 fb^-1` under the paper's explicit
zero-background working assumption.

Production is normalized inclusively with `b -> X_s S` (Winkler
arXiv:1809.01876, Eq. A7), summed over the `B+`, `B0`, `Bs`, and `Lambda_b`
b-hadron pool (the rate is spectator-independent). The unobserved strange system
is represented kinematically by two-body recoil -- a kaon for the mesons and a
physical `Lambda` for `Lambda_b`; this proxy does not change the inclusive
normalization. Direct `gg -> S` production belongs outside the minimal b-hadron
BC4 definition used here.

The central decay model uses the digitized Winkler dispersive hadronic widths
below 2 GeV and the perturbative spectator treatment above 2 GeV, together
with analytic leptonic widths. The alternate LO-ChPT/spectator model is used
only as a named model diagnostic.

## Variation diagnostics and figures

`scalar/data/published/bundle/` contains 109 coherent FONLL-grid curves, one
alternate decay-model curve, and two excluded same-physics numerical repeats.
The derived outer envelope is a one-source-at-a-time diagnostic in
`log10(sin^2(theta))`; it is not a confidence interval and sources are not
combined in quadrature.

The publication comparison is produced in `shared/curves_PBC`, not by the
package-local legacy overlay loader. It includes the current excluded region
and the in-scope SHiP, CODEX-b, and FASER2 projections. MATHUSLA and ANUBIS are
omitted because no current-geometry mixing-only curves are available. By user
policy, no ATLAS, CMS, or LHCb projection curve is used; existing LHCb bounds
remain in scope.

## Remaining limitations

- Detector response and backgrounds are not yet validated at publication
  fidelity; the paper labels zero background as a working assumption.
- `K -> pi S` exists in the model layer, but a realistic LHC kaon-flux sample
  is not included in the central production model.
- The strange recoil system and exclusive decay topology are approximations
  documented in the package and manifests.
- The named theory/model curves are diagnostics rather than probabilistic
  uncertainty bands, so the primary proposed-experiment plot remains
  central-only.

For rerun commands use `scalar/README.md`; for external-input history and plot
policy use `scalar/EXTERNAL_INPUTS_NEEDED.md`.
