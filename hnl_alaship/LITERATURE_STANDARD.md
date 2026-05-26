# HNL Exclusion Curves — Literature Standard for Transverse LLP Detectors

Reference for the GARGOYLE detector proposal (Citron & Redi). Captures the
conventions used by MATHUSLA, ANUBIS, CODEX-b, FASER2, SHiP, and the PBC
working group, so that GARGOYLE's HNL plots match the field's expectations.
Companion to `README.md` (architecture) and `HOWTO_RUN.md` (commands).

## 1. Single-flavor benchmarks (PBC BC6 / BC7 / BC8)

The community standard is to scan **one mixing element at a time**. Three
benchmarks, three plots:

| Benchmark | $U_e^2 : U_\mu^2 : U_\tau^2$ | Flavor tag |
|-----------|-------------------------------|------------|
| BC6       | 1 : 0 : 0                     | `Ue`       |
| BC7       | 0 : 1 : 0                     | `Umu`      |
| BC8       | 0 : 0 : 1                     | `Utau`     |

Reference: Antusch et al., *New Benchmark Models for Heavy Neutral Lepton
Searches*, Eur. Phys. J. C 83 (2023) 1110 (arXiv:2207.02742). Adopted by the
Physics Beyond Colliders BSM working group report (J. Phys. G 47 (2020)
010501, arXiv:1901.09966).

Combined-coupling scenarios exist in the literature but are non-standard for
new-detector proposals; do not lead with them.

## 2. Plot conventions

- Axes: $m_N$ [GeV] on x (linear or log), $|U_\alpha|^2$ on y, log scale.
- One figure per flavor (BC6, BC7, BC8).
- 95% CL exclusion contour, closed curve in $(m_N, |U|^2)$.
- Mass range: $0.1$–$15$ GeV (paper abstract). Realistic reach for a
  transverse LHC detector closes well below 15 GeV; do not draw past where
  $N_s$ falls below threshold.
- Standard overlays:
  - **BBN floor** (lower bound on $|U|^2$ from cosmology, typically below
    $\sim 10^{-9}$).
  - **Type-I seesaw band** between active-neutrino mass scales
    $0.05$ eV $\le m_\nu \le 0.12$ eV (light grey band).
  - **Existing experimental constraints**: PS191, CHARM, NA62, T2K, Belle,
    BEBC, DELPHI, ATLAS/CMS/LHCb prompt searches (whichever are relevant
    per flavor).
  - **Competitor proposed detectors**: MATHUSLA, ANUBIS, CODEX-b, FASER2,
    SHiP. SHiP is beam-dump (different $\sqrt{s}$, different luminosity);
    include it for completeness but label clearly.
- Luminosity tag in caption: $3\,\mathrm{ab}^{-1}$, $\sqrt{s}=14$ TeV.

## 3. Statistics

- 95% CL exclusion ⇔ Poisson upper limit of 3 expected signal events under
  zero-background hypothesis. This is the convention used by MATHUSLA,
  CODEX-b, FASER2, SHiP, and most PBC HNL projections.
- ANUBIS additionally publishes a $B = 90$ conservative-background curve;
  if GARGOYLE backgrounds are non-negligible after vetoes, mirror this with
  a second contour at the relevant $B$.
- Threshold convention is set in `hnl_alaship/analysis/exclusion.py` —
  current value $N_s \ge 3$ matches the standard.

## 4. Production channels by mass region

| Mass region                   | Dominant channels                                | Tools                              |
|-------------------------------|--------------------------------------------------|------------------------------------|
| $m_N \lesssim m_D$ (~1.87 GeV)| $D, D_s$ leptonic + semileptonic                 | FONLL grids + HNLCalc prod BRs     |
| $m_D \lesssim m_N \lesssim m_B$ | $B, B_s$ leptonic + semileptonic                 | FONLL grids + HNLCalc prod BRs     |
| $m_N \lesssim m_\tau$ (1.78 GeV), BC8 only | $\tau \to N \ell \nu$ via $W^*$        | MadGraph + HNLCalc tau BRs         |
| $m_N \gtrsim 5$ GeV           | $W^\pm \to \ell^\pm N$, $Z \to \nu N$ direct      | MadGraph (`SM_HeavyN_CKM_AllMasses_LO`) |
| Above $m_B$ ~5.3 GeV          | $B_c$ via $b\bar{c}$                              | FONLL (small contribution)         |

Tau-parent matters specifically for BC8 below the $\tau$ threshold —
omitting it loses the BC8 reach in that window. Bc and W/Z extend the
high-mass reach above the meson threshold.

## 5. Acceptance methodology — two rigor levels in the literature

| Level | Used by | Decay model | Daughter check |
|-------|---------|-------------|----------------|
| Analytic | MATHUSLA, ANUBIS, CODEX-b, DDC (arXiv:2512.13011) | Exponential decay-in-volume formula: $\varepsilon = (\delta\phi/2\pi)\,e^{-D/\beta\gamma c\tau}\,(1-e^{-L/\beta\gamma c\tau})$ | None — assumes 100% daughter detection if HNL decays inside |
| Full MC | SHiP/FairShip, **GARGOYLE** | Pythia/FairShip rest-frame templates + 3D Lorentz boost to lab | Daughter ray-cast against detector mesh, $p$-cut + separation cuts |

GARGOYLE uses the **full MC + ray-cast** approach (`hnl_alaship/analysis/`),
which is more conservative at high $m_N$ where daughters can fly out of the
~3 m × 3 m tunnel cross-section. This means **GARGOYLE-vs-MATHUSLA/ANUBIS/
CODEX-b on the same plot is not apples-to-apples**. See `README.md` §
"Acceptance model" for the impact (15–30% tighter contour above $m_N \sim
1.5$ GeV; reach closes ~3.8 GeV instead of ~5 GeV in BC7).

Two acceptable presentations:
1. **Single curve, footnoted.** Show the rigorous GARGOYLE contour and note
   in the caption that competitor curves use analytic acceptance.
2. **Dual curve.** Add a dashed "analytic-acceptance" GARGOYLE contour for
   apples-to-apples comparison, alongside the rigorous solid contour.

Option 2 is the strongest defense against "you're cherry-picking" objections
in review — the dashed curve will be tighter and visually closer to the
competitor curves at high mass; the gap between dashed and solid is the
honest impact of full daughter containment.

## 6. Decay BRs and ctau

- ctau($U^2 = 1$) and visible-decay BRs come from FairShip's `hnl.py`
  (`production/generate_ctau_tables.py`).
- Production BRs (meson → HNL) come from HNLCalc's `branchingratios.dat`
  (correct for production; the issue documented in `README.md` is on the
  *decay* side, not production).
- Caveat: HNLCalc decay tables miss 17/38 three-body channels — that's why
  decay-side calculations use FairShip instead. Do not regress this choice.

## 7. Required overlay curves — sourcing

Reference curves live in
`/Volumes/sandbox/projects/aaaPHYSICSaaa/curves_from_others/hnl_matrix/`
(separate repo, treat as the canonical curve library). Per-flavor coverage
expected:

- **BC6 (Ue)**: MATHUSLA, ANUBIS, CODEX-b, SHiP, existing exclusions (NA62,
  T2K, PS191, CHARM, BEBC, DELPHI).
- **BC7 (Umu)**: MATHUSLA, ANUBIS, CODEX-b, SHiP, existing (NA62, T2K,
  PS191, BEBC, DELPHI).
- **BC8 (Utau)**: MATHUSLA, ANUBIS (if available), SHiP, existing (CHARM,
  DELPHI). Curves are sparser for BC8; if a competitor is not published,
  omit and note in caption.

Plotting code: `analysis/reference_curves.py` and
`analysis/plot_exclusion.py`. Add new curves under
`curves_from_others/hnl_matrix/published/` and register in the source
registry there, not inline in `hnl_alaship/`.

## 8. Pipeline-vs-standard checklist (audit 2026-04-30)

Verified:
- [x] All three flavors implemented end-to-end (production, decay, sensitivity,
      plot). `production/decay_engine/generate_meson_csvs.py:69`,
      `production/madgraph/run_tau_production.py:79`,
      `production/madgraph/run_wz_production.py:90`,
      `analysis/run_sensitivity.py:302`.
- [x] Single-coupling scans, no combined mixing. Each flavor produces its own
      production sample and ctau table; `run_full_pipeline.sh:69-116`.
- [x] $N_s \ge 3$ threshold. `analysis/constants.py:13`
      (`N_THRESHOLD = 3.0`), used in `analysis/exclusion.py:16-57`.
- [x] Log–log axes, separate panels per flavor.
      `analysis/plot_exclusion.py:206-209`, x-limit [0.15, 6.0] GeV,
      y-limit [1e-12, 1e-1].
- [x] Tau-parent for BC8 — `production/madgraph/run_tau_production.py`,
      no artificial cutoff below $m_\tau$.
- [x] $B_c$ for high-mass reach — handled in
      `production/decay_engine/generate_meson_csvs.py:45-49,111-114` via
      fragmentation fractions.
- [x] HL-LHC luminosity 3 ab$^{-1}$ — `analysis/constants.py`.
- [x] **Pipeline already runs end-to-end and produces
      `output/analysis/gargoyle_hnl_exclusion.{pdf,png}`** (single figure
      with three panels, one per flavor; last full run 2026-03-13,
      4207 s wall time, 5000 mother samples / 12 position bins / 20000
      templates / 10 decays per bin). All caches under
      `output/analysis/{geometry_cache,decay_cache}/{Ue,Umu,Utau}/` valid.
      To refresh from cache: `conda run -n llpatcolliders python analysis/run_sensitivity.py --plot-only`.

Resolved (this commit):
- [x] **All 12 reference curves staged** for ANUBIS, CODEX-b, MATHUSLA, and
      PastExclusion across BC6/BC7/BC8. Sourced via the
      `curves_from_others/hnl_matrix/published/manifest.csv` provenance map.
      File format: 2-column `mass_GeV  u2` with a `# kind: contour` or
      `# kind: envelope` header indicating whether the file traces a closed
      island (drawn as a single connected line) or a "ceiling" curve (the
      shaded region above is excluded). Loader: `analysis/reference_curves.py`.
- [x] **Type-I seesaw band overlay** added to all three panels via
      `_plot_seesaw_band()` in `analysis/plot_exclusion.py`. The band uses
      $|U|^2 = m_\nu / m_N$ for $m_\nu \in [0.05, 0.12]$ eV.
- [x] **Past-experiment exclusion overlay** wired through as the
      `PastExclusion` curve set (HNLimits compilation). Plotted as a shaded
      "everything above is excluded" region; values smoothed with a rolling
      log-mass minimum to suppress source-handoff discontinuities.

Still open:
- [ ] **No BBN floor overlay.** Independent of past-collider exclusions; add
      a horizontal/curved line at the cosmologically allowed $|U|^2$ ceiling.
- [ ] **No SHiP curves.** No clean manifest source available; the
      `SHiP_Umu.dat` previously checked in was hand-digitized and has been
      removed pending a proper conversion.
- [ ] **Optional dual-curve presentation.** A dashed analytic-acceptance
      contour alongside the solid full-MC contour would defend against
      "not apples-to-apples" review objections (see §5).
- [ ] **Captions** must state $\sqrt{s} = 14$ TeV, $\mathcal{L} = 3\,\mathrm{ab}^{-1}$,
      acceptance model (full MC + daughter ray-cast), BR/ctau source
      (FairShip), $N_s \ge 3$ at 95% CL.

## 9. Paper-text alignment (Redi & Citron draft)

The draft `\paragraph{Heavy Neutral Leptons}` in §Sensitivity is currently
empty. When written it should:
- Reference BC6/BC7/BC8 explicitly by name and cite arXiv:2207.02742.
- State acceptance methodology (full 3D boost + daughter ray-cast against
  the PX56 tunnel mesh) and acknowledge the apples-to-apples caveat.
- Cite competitor curves: MATHUSLA (arXiv:1811.00927, 2308.05860), ANUBIS
  (arXiv:1909.13022, 2512.14942), CODEX-b (arXiv:1911.00481), SHiP
  (arXiv:1504.04855, 2112.01487), FASER2 (arXiv:1811.12522). Latest
  comparison study: arXiv:2512.13011.
- State $\mathcal{L} = 3\,\mathrm{ab}^{-1}$, 95% CL, $N_s \ge 3$.

## 10. Sources

- Antusch et al., *New Benchmark Models for HNL Searches*, EPJC 83 (2023)
  1110 — arXiv:2207.02742.
- PBC BSM working group, *Physics Beyond Colliders at CERN*, J. Phys. G 47
  (2020) 010501 — arXiv:1901.09966.
- DDC 2026 update — arXiv:2512.13011.
- ANUBIS sensitivity & proANUBIS — arXiv:2512.14942.
- HNLs at ANUBIS (original) — Phys. Rev. D 101 (2020) 055034.
- MATHUSLA HNL public files — github.com/davidrcurtin/MATHUSLA_LLPfiles_RHN_{Ue,Umu,Utau}.
- Sensitivity of the intensity frontier — arXiv:1902.06240.
