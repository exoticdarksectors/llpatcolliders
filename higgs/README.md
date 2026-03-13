# Higgs Portal: h -> SS

Dark scalar S (PDG 6000113) produced via SM Higgs decay. Default mass 15 GeV, decay S -> b bbar.

## Cross-section

`--xsec 54700` (sigma(gg->h) = 54.7 pb, N3LO, 14 TeV)

## Production

```bash
bash generator/produce.sh higgs/cmnd/higgs_h_to_ss.cmnd 50000 higgs_m15 4 --outdir output/higgs
```

## Analysis

```bash
python analysis/decayProbPerEvent_2body.py output/higgs/data/higgs_m15.csv \
    --xsec 54700 --outdir output/higgs --external-dir higgs/external

python analysis/signal_surface_hitmap.py output/higgs/data/higgs_m15.csv \
    --outdir output/higgs
```

## Background

```bash
python analysis/background_trident.py --outdir output/higgs
```

## Note on acceptance model

The current analysis uses the analytic 2-body acceptance formula (opening angle
from Lorentz kinematics, separation = theta × d_remaining). This is a good
approximation for 2-body decays at high gamma but does not ray-cast individual
daughters to the detector mesh. For the HNL portal (multi-body decays), a full
3D boost + daughter ray-cast was found to give 15-30% tighter exclusion above
m > 1.5 GeV (see `hnl_alaship/`). For 2-body scalar decays the effect is
expected to be small, but a cross-check with the 3D method at low gamma
(high mass) is recommended as a future validation.
