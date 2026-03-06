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
