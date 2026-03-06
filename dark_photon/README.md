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
