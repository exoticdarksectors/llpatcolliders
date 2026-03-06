# ALP Portal

Two production channels:

- **Heavy ALP** (PDG 6000113): h -> aa, mass ~ 1-60 GeV, decay a -> mu+ mu-
- **Light ALP** (PDG 9000001): B -> K(*) a via FCNC, mass ~ 0.2-4 GeV

Both cmnd files use BR=1 (efficiency maps); rescale by true BR via `--xsec`.

## Cross-sections

- Heavy ALP: `--xsec 54700` (sigma(gg->h) = 54.7 pb, 14 TeV)
- Light ALP: `--xsec 373000000` (sigma(pp->bb) ~ 373 ub, inclusive)

## Production

```bash
bash generator/produce.sh alps/cmnd/heavy_alp.cmnd 50000 alp_heavy_m15 4 --outdir output/alps
bash generator/produce.sh alps/cmnd/light_alp.cmnd 50000 alp_light_m1 4 50000 --outdir output/alps
```

## Analysis

```bash
python analysis/decayProbPerEvent_2body.py output/alps/data/alp_heavy_m15.csv \
    --xsec 54700 --outdir output/alps --external-dir higgs/external

python analysis/decayProbPerEvent_2body.py output/alps/data/alp_light_m1.csv \
    --xsec 373000000 --outdir output/alps

python analysis/signal_surface_hitmap.py output/alps/data/alp_heavy_m15.csv \
    --outdir output/alps
```

**Note:** `alps/external/` is intentionally empty. The heavy ALP (h->aa) uses the
same PBC benchmark curves (h->SS) as the Higgs portal, located in `higgs/external/`.
