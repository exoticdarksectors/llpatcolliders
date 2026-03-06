# LLP at Colliders — PX56 Tunnel Detector

Sensitivity studies for long-lived particle (LLP) searches at the CMS PX56 drainage tunnel, for HL-LHC (14 TeV, 3000 fb^-1).

## Repository Structure

```
generator/          Shared Pythia8 event generator (one binary for all models)
geometry/           Shared PX56 tunnel geometry (gargoyle_geometry.py)
analysis/           Shared analysis scripts (decay probability, surface hitmap)

higgs/              Higgs portal: h -> SS (dark scalar)
alps/               ALP portal: h -> aa (heavy) and B -> K(*)a (light)
dark_photon/        Dark photon: h -> A'A' (R-ratio BRs)
FCC/                FCC detector studies (notebooks)
```

Each model directory contains: `cmnd/` (Pythia8 configs), `external/` (comparison curves), and a `README.md` with run instructions.

## Quick Start

```bash
# 1. Build the generator (requires Pythia8)
cd generator && bash build.sh && cd ..

# 2. Generate events (example: heavy ALP, 50k events, 4 parallel jobs)
bash generator/produce.sh alps/cmnd/heavy_alp.cmnd 50000 alp_heavy_m15 4 --outdir output/alps

# 3. Run analysis
python analysis/decayProbPerEvent_2body.py output/alps/data/alp_heavy_m15.csv \
    --xsec 54700 --outdir output/alps --external-dir alps/external

# 4. Surface hitmap
python analysis/signal_surface_hitmap.py output/alps/data/alp_heavy_m15.csv \
    --outdir output/alps
```

## Analysis Scripts

| Script | Method | Use for |
|--------|--------|---------|
| `decayProbPerEvent_2body.py` | Analytical 2-body acceptance | Any 2-body decay (higgs, alps) |
| `decayProbPerEvent_Ntrack.py` | N-track daughter pairing | Multi-body decay (dark photon) |
| `signal_surface_hitmap.py` | Surface coverage optimization | Any model |
| `background_trident.py` | Beam muon trident background MC | Any model |

## Adding a New Model

1. Create `<model>/cmnd/` with Pythia8 config (set `LLP:pdgId`)
2. Add external comparison curves to `<model>/external/`
3. Use the shared generator and analysis scripts
