# Dark Photon Portal

Dark photon A' (PDG 6000115, vector) produced via h -> A'A'. R-ratio branching ratios from DeLiVeR VMD+PDG below 1.7 GeV, perturbative QCD above.

Mass benchmarks: 0.5, 1.0, 1.5, 2.0 GeV.

Also explored (zero PX56 sensitivity confirmed):
- Meson decay: eta -> A' gamma, omega -> A' pi0 (SoftQCD)
- Drell-Yan: qq -> A' (NewGaugeBoson)

## Cross-section

`--xsec 54700` (sigma(gg->h) = 54.7 pb, 14 TeV)

## Generating cmnd files

```bash
python dark_photon/make_dp_cmnd.py --mass 2.0 --outfile dark_photon/cmnd/heavy_dp_m2.cmnd
```

## Production

```bash
bash generator/produce.sh dark_photon/cmnd/heavy_dp_m2.cmnd 50000 dp_heavy_m2 4 --outdir output/dark_photon
```

## Analysis

```bash
python analysis/decayProbPerEvent_Ntrack.py output/dark_photon/data/dp_heavy_m2.csv \
    --xsec 54700 --outdir output/dark_photon --external-dir dark_photon/external

python analysis/signal_surface_hitmap.py output/dark_photon/data/dp_heavy_m2.csv \
    --outdir output/dark_photon
```

## BR validation

```bash
python dark_photon/validate_brs.py
```

## Note on BR table

The VMD BR table (`br_tables/dp_brs_deliver.csv`) has 22 mass points in 0.2–1.7 GeV
with only 3 points in the rho/omega resonance region (0.70, 0.77, 0.80 GeV). Linear
interpolation through resonances is coarse but validated against DeLiVeR direct
computation via `validate_brs.py` (< 20% in resonance windows, < 10% elsewhere).

Above 1.7 GeV, perturbative QCD R-ratio BRs are used. There is a ~15% BR
discontinuity at the VMD/pQCD boundary — this is physical (genuine mismatch
between the two regimes) and propagates as ~15% discontinuity in the exclusion
contour. Within systematic uncertainty for a sensitivity study.

Channels with BR < 0.5% (`MIN_BR = 0.005`) are merged into pi+pi-. This affects
minor hadronic channels and simplifies the Pythia8 decay table while preserving
the dominant kinematic features.

## Note on acceptance model

The current analysis uses the N-track acceptance model (`decayProbPerEvent_Ntrack.py`)
with analytic formulas. This does not ray-cast individual daughters to the
detector mesh. For the HNL portal (multi-body decays), a full 3D boost +
daughter ray-cast was found to give 15-30% tighter exclusion above m > 1.5 GeV
(see `hnl_alaship/`). For dark photon decays (mostly 2-body at low mass,
multi-hadron above 1 GeV via R-ratio), a cross-check with the 3D method is
recommended as a future validation, especially for the higher mass benchmarks.
