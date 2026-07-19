# External inputs — BC10 fermiophilic ALP

_Updated 2026-07-13._ The production normalization and the arXiv:2310.03524
decay baseline have been replaced by external data products. The
arXiv:2501.04525 widths and exclusive branching ratios are now the default,
and its production modes have been audited against the GRENDEL acceptance.

## 1. Absolute B → K a production normalization — SATISFIED

The leading-log placeholder (`C_TOP = 1`) was replaced by the one-loop RG
coefficient of GKOZ (arXiv:2310.03524), evaluated with ALPINIST's
implementation (`tools/compute_cbs_alpinist.py`, BSD-3 port; pinned commit in
`data/alpinist/COMMIT_SHA.txt`):

    g_bs = CBS_EFF * (1/f),  CBS_EFF = 3.518383e-4  at Lambda_UV = 1 TeV.

The finite one-loop + RG terms partially cancel the small log at 1 TeV: the
proper coefficient is ~2.5x *smaller* in amplitude (6.3x in rate) than the
leading-log placeholder. Production now also spans the **full kaon tower**
(K, K0*(700/1430), K*(892/1410/1680), K1(1270/1400), K2*(1430)) for B+ and B0
with the Boiarska et al. (arXiv:1904.10447) form factors — GKOZ's "almost 4x"
uplift over K + K*(892) (we get 3.9-4.1x depending on m_a).

**Residual:** cross-checking CBS_EFF against GKOZ Table 1 (|C_bs| = 1.8e-3 at
f_GKOZ = 1 GeV) agrees to 6-20% in amplitude depending on the m_b scheme used
to unfold their C_bs = c_bs m_b/2f definition. The island's lower
(production-limited) edge inherits <~20% in 1/f.

`B_s -> phi a` is not included. It is also absent from the pinned ALPINIST
implementation and from the exclusive channel set in GKOZ. Boiarska et al.
provides the kaon-tower form factors for `B+`/`B0`, but not a corresponding
`B_s` tower. Therefore an `O(10%)` uplift is not a sourced result of those
references. Using `f_s/(f_u+f_d) = 0.122` and a ground-state `B_s -> phi` rate
comparable to `B -> K*(892)` suggests only a few-percent yield correction;
adding it requires a separately pinned `B_s -> phi` form-factor calculation.

## 2. Data-driven decay model — 2501 TABLE UPGRADE SATISFIED

`model.alp_partial_widths` and `model.alp_total_width` read the exact exported
SensCalc v1.3.3 arXiv:2501.04525 tables over 0.01--10 GeV. The export includes
the widths, exclusive branching ratios, and squared matrix elements, with
source hashes and conversion checks in `data/senscalc_2501/PROVENANCE.md`.
The exact 1 GeV anchor is `BR(a -> mu mu) = 0.2022433833`.

Publication templates sample the full exclusive branching mixture, let
Pythia decay unstable daughters and hadronize partonic modes, and leave the
actual two-track decision to the GRENDEL reconstruction. The old visible-only
two-track proxy is gated behind `templates.py --legacy-proxy`.

Three-body primary decays are generated with Pythia's flat phase-space mode
and reweighted event by event with the exact exported squared matrix elements.
The weights are normalized within every exclusive channel, preserving the
imported branching mixture while replacing its Dalitz shape. Numerical
evaluation of the exported `CForm` expressions is cross-checked against
Wolfram Engine anchors in `tests/test_decay_matrix_elements.py`.

**Residual:** the two-gluon mode uses an equal `u/d/s` jet surrogate because
the Pythia external-decay interface cannot fragment a bare colour-singlet
gluon pair. This affects decay acceptance, not the imported total width or
branching ratios, and remains a decay-model systematic rather than exact 2501
partonic kinematics.

**Convention note (for future overlays / axis labels):** our 1/f axis is the
BNT (arXiv:1708.00443) convention g_aff = c_f m_f / f with c_f = 1. GKOZ
normalise with 1/(2 f_GKOZ), so 1/f_here = 2/f_GKOZ. Any comparison curve
digitized from GKOZ/PBC BC10 figures must be mapped accordingly.

## 3. Production modes discussed in arXiv:2501.04525 — AUDITED

The current GRENDEL signal contains the exclusive `B+`/`B0` kaon tower only.
ArXiv:2501.04525 discusses a broader proton-collision production policy:
exclusive B decays, Drell-Yan/gluon fusion, proton bremsstrahlung, light-meson
decays, and quark fragmentation using generalized rather than bare meson-ALP
mixing angles. The latter modes cannot be declared covered by the current
code. Their relevance must be evaluated with the GRENDEL geometry and LHC
kinematics; an inclusive production-probability comparison alone is not an
acceptance calculation. SensCalc `v.1.3.3` is the pinned source for this audit.

The source-level findings and per-channel decisions are recorded in
`data/senscalc_2501/PRODUCTION_AUDIT.md`. Drell-Yan is negligible throughout
the current island compared with the open B tower, and the obsolete
flux-times-mixing shortcut is explicitly rejected. The decoded generalized
fragmentation contribution is at most 0.15% of the B tower outside the
excluded light-meson pole windows. Light-meson decay rates are over five
orders smaller than accepted B production at low mass. Their omission is
therefore documented and quantified rather than assumed.

## 4. Competitor + existing-bound curves — SATISFIED OUTSIDE THIS PACKAGE

The signal package publishes the GRENDEL contour without duplicating external
limits. The sibling `curves_PBC` repository builds the paper comparison from
source-pinned GKOZ LHCb, arXiv:2501.04525 beam-dump, ALPINIST NA62, and selected
FIPs-2022 existing bounds, plus the published SHiP and DarkQuest projections.
ATLAS, CMS, and LHCb projection curves are excluded by analysis policy. The
local `plot.py` still accepts an optional long-format overlay CSV
(`curve, m_a_GeV, invf_GeV_inv`) for focused validation plots.
