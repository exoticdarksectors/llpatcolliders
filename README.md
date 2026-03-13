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
hnl_alaship/        HNL portal: FairShip backend (full 3D boost, all 3 flavors) [ACTIVE]
hnl_legacy/         HNL portal: original z-boost pipeline (superseded, kept for reference)
FCC/                FCC detector studies (notebooks)
```

Each model directory contains comparison curves in `external/` and a `README.md` with run instructions. Higgs/ALP/dark photon use the shared Pythia8 generator (`cmnd/` configs). The active HNL pipeline is `hnl_alaship/` (FairShip-based, full 3D Lorentz boost, all three flavors Ue/Umu/Utau). The original z-boost pipeline is preserved in `hnl_legacy/` for reference.

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

### HNL Portal (FairShip / hnl_alaship)

Uses FairShip hnl.py for decay widths/lifetime, TPythia8 for rest-frame decay
templates, and a full 3D Lorentz boost for daughter acceptance. Supports all
three flavors (Ue, Umu, Utau). Requires the `llpatcolliders` conda env with ROOT.

```bash
cd hnl_alaship
conda run -n llpatcolliders python production/decay_engine/generate_meson_csvs.py --flavor Ue Umu Utau --channel all
conda run -n llpatcolliders python production/madgraph/run_tau_production.py --flavor Ue Umu Utau
conda run -n llpatcolliders python production/madgraph/run_wz_production.py --flavor Ue Umu Utau
conda run -n llpatcolliders python production/combine_channels.py --flavor Ue Umu Utau
conda run -n llpatcolliders python production/generate_ctau_tables.py --flavor Ue Umu Utau
conda run -n llpatcolliders python analysis/run_sensitivity.py --flavor Ue Umu Utau
```

See `hnl_alaship/HOWTO_RUN.md` for details.

## Analysis Scripts

| Script | Method | Use for |
|--------|--------|---------|
| `decayProbPerEvent_2body.py` | Analytical 2-body acceptance | Any 2-body decay (higgs, alps) |
| `decayProbPerEvent_Ntrack.py` | N-track daughter pairing | Multi-body decay (dark photon) |
| `signal_surface_hitmap.py` | Surface coverage optimization | Any model |
| `background_trident.py` | Beam muon trident background MC | Any model |

## Acceptance Models

Three distinct acceptance approaches are used across portals (must be documented
in the paper as a systematic):

| Model | Script | Portals | Method | Limitation |
|-------|--------|---------|--------|------------|
| Analytic 2-body | `decayProbPerEvent_2body.py` | Higgs, ALPs | Boost-invariant opening angle formula (massless daughters) | Only valid for stable 2-body final states (e+e-, mu+mu-). **Not valid for tau+tau-.** |
| N-track | `decayProbPerEvent_Ntrack.py` | Dark photon | Daughter-pair opening angles from generator output | Approximates detector as flat plane at fixed distance. No per-daughter ray-cast. |
| Full 3D MC | `hnl_alaship/analysis/` | HNL | Full 3D Lorentz boost + per-daughter ray-cast to mesh | Most accurate. 15-30% tighter exclusion than z-boost above m > 1.5 GeV. |

Cross-validation at shared mass points between models is recommended before
publication.

## Adding a New Model

1. Create `<model>/cmnd/` with Pythia8 config (set `LLP:pdgId`)
2. Add external comparison curves to `<model>/external/`
3. Use the shared generator and analysis scripts
