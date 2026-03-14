# hnl_alaship — FairShip-based HNL Pipeline

FairShip-backed GARGOYLE HNL sensitivity pipeline, replacing:
- **MATHUSLA decay templates** with **TPythia8** (via ROOT, in conda env)
- **HNLCalc decay BRs/ctau** with **FairShip hnl.py** (fixes double-counting bug above 1 GeV)
- **z-only boost** with **full 3D Lorentz boost** along actual parent 3-momentum (15-30% more accurate above 1.5 GeV)

HNLCalc is **kept for production BRs only** (correct there; `branchingratios.dat` is missing 17/38 three-body channels accounting for 28-45% of yield at 1.0-1.25 GeV).

## Flavors

All three flavors are supported from the start:
- `Ue` — electron coupling `[1, 0, 0]`
- `Umu` — muon coupling `[0, 1, 0]`
- `Utau` — tau coupling `[0, 0, 1]`

## Architecture

```
PRODUCTION (HNL 4-vectors + weights)
├── Heavy-flavor (Bmeson, Dmeson, Bc):  FONLL grids + HNLCalc (production BRs)
├── Tau-parent:                          MadGraph (Docker) + HNLCalc (tau BRs)
└── W/Z direct:                          MadGraph (Docker)

GEOMETRY (which mothers hit the detector)
└── Ray-cast mothers from IP to mesh:    gargoyle_geometry.py (shared)

DECAY + ACCEPTANCE (HNL decays inside detector)
├── ctau(U²=1):                          FairShip hnl.py
├── Rest-frame decay templates:          FairShip TPythia8/ROOT
├── 3D boost to lab frame:               fairship_decay.py (general Lorentz boost)
└── Daughter ray-cast + P_CUT + sep:     decay_acceptance.py

SENSITIVITY
├── U² scan:                             sensitivity.py
├── Exclusion band:                      exclusion.py
└── Plot + overlays:                     plot_exclusion.py + reference_curves
```

## Boost method

The `boost_decay_to_lab` function performs a **full 3D Lorentz boost** along
the actual parent 3-momentum direction (not a z-only approximation). This
gives 15-30% tighter exclusion contours above m_N ~ 1.5 GeV.

## Acceptance model: daughter containment

This pipeline uses **full daughter ray-casting** against the tunnel mesh,
which is more rigorous than the analytic 2-body formula used by the previous
`hnl/` pipeline and by most PBC comparison experiments (ANUBIS, MATHUSLA,
CODEX-b). The key differences and their impact:

### What changed vs. the old `hnl/` (z-boost) pipeline

| Feature | Old `hnl/` | This pipeline (`hnl_alaship/`) |
|---------|-----------|-------------------------------|
| Daughter acceptance | Analytic 2-body: fraction of cos θ* range passing cuts (infinite-plane assumption) | MC decay templates + ray-cast each daughter against the finite tunnel mesh |
| Daughter mass | M_electron (0.5 MeV) for all — most forgiving for P_CUT | Actual masses from TPythia8 (pions, muons, etc.) |
| Decay topology | 2-body only | Multi-body (all channels from FairShip) |
| Daughter containment | Not checked — assumes daughters hit detector if HNL decays inside | Explicitly required — daughters must intersect the mesh |

### Why this reduces high-mass reach

The GARGOYLE tunnel cross-section is ~2.9 m × 3.15 m. At high m_N the HNL
rest-frame decay products carry more energy and are emitted more isotropically.
After boosting to the lab, a significant fraction of daughters fly **out of
the tunnel** even when the HNL decayed inside. The old analysis missed this
because it never checked whether daughters actually hit the detector walls.

Quantitative impact at Umu:

| m_N (GeV) | Old peak_N | New peak_N | Notes |
|-----------|-----------|-----------|-------|
| 3.0 | ~60 | ~15 | Meson-dominated, moderate reduction |
| 4.0 | ~10 | ~2.2 | Falls below threshold (3.0) |
| 5.0 | ~4 | ~0.5 | Only Bc + WZ contribute |

The sensitivity contour closes at **~3.8 GeV** (new) vs **~5 GeV** (old).

### Comparison with other experiments

Published sensitivity curves for ANUBIS, MATHUSLA, CODEX-b, and SHiP
typically use analytic acceptances without full daughter containment in a
finite detector volume. Overlaying GARGOYLE (rigorous) against these
(optimistic) is **not an apples-to-apples comparison**.

Options for the paper:
1. Show the rigorous result as the primary GARGOYLE curve, and note
   that reference curves do not include daughter containment.
2. Add a dashed "analytic acceptance" contour (no daughter ray-casting)
   for direct comparison with PBC experiments.
3. Both.

### W/Z channel cross-section (verified correct)

The W/Z production weight (~29,000 pb at U²=1, constant across m_N) has
been verified against the expected σ(pp→W)×BR(W→ℓN) ≈ 22.5 nb, consistent
with MadGraph LO × K_FACTOR_EW = 1.3. The weight is extracted correctly
from the LHE event header by `lhe_to_csv.py` (XWGTUP field = σ/N_events).
Cross-sections are NOT the cause of the reduced reach.

### Geometric acceptance (verified correct)

GARGOYLE geometric acceptance is ~0.7% (2100/300k events hit the tunnel).
Mean path through the detector is ~2.3 m, mean entry distance ~28 m.
Compared to ANUBIS (~10% acceptance, ~18 m path length), GARGOYLE has
roughly 100× less signal yield, which is a genuine geometric limitation
of the smaller detector — not a bug.

## Key constants

| Parameter | Value |
|-----------|-------|
| Luminosity | 3000 fb⁻¹ (HL-LHC) |
| Exclusion threshold | N_signal ≥ 3 (95% CL, zero background) |
| U² scan | 10⁻¹² – 10⁻¹, 200 log points |
| P_CUT | 0.600 GeV/c |
| SEP_MIN / SEP_MAX | 1 mm / 1 m |
| FONLL mass max | 10.0 GeV (WZ extends above meson threshold) |
| Mass grid | 116 points × 3 flavors = 348 mass points |

Verified result: **270/348 mass points have sensitivity** (full 3D daughter ray-casting).

## Dependencies

- `conda env: llpatcolliders` (includes ROOT 6.38+ with Pythia8)
- Shared geometry: `../geometry/gargoyle_geometry.py`
- Docker (for MadGraph tau/WZ production only)

## See also

- `HOWTO_RUN.md` for step-by-step instructions
- `../geometry/gargoyle_geometry.py` for the shared detector geometry
