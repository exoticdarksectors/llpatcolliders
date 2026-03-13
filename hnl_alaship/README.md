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

## Key constants

| Parameter | Value |
|-----------|-------|
| Luminosity | 3000 fb⁻¹ (HL-LHC) |
| Exclusion threshold | N_signal ≥ 3 (95% CL, zero background) |
| U² scan | 10⁻¹² – 10⁻¹, 200 log points |
| P_CUT | 0.600 GeV/c |
| SEP_MIN / SEP_MAX | 1 mm / 1 m |
| FONLL mass max | 5.0 GeV |
| Mass grid | 96 points × 3 flavors = 288 mass points |

Verified result: **268/288 mass points have sensitivity** (Ue: 89/96, Umu: 89/96, Utau: 90/96).

## Dependencies

- `conda env: llpatcolliders` (includes ROOT 6.38+ with Pythia8)
- Shared geometry: `../geometry/gargoyle_geometry.py`
- Docker (for MadGraph tau/WZ production only)

## See also

- `HOWTO_RUN.md` for step-by-step instructions
- `../geometry/gargoyle_geometry.py` for the shared detector geometry
