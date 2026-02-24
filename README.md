# Tunnel LLP Detector — Active Workflow

This repository is currently centered on a Pythia generation + signal acceptance pipeline for LLPs in the CMS PX56 tunnel.

## Current workflow

See `howto.md` for the full command sequence.

Pipeline stages:
1. **Production** — `pythiaStuff/parallel_produce.sh` runs `main144` in parallel batches until a target number of LLP rows is reached. Both cmnd files use BR=1.0 (efficiency maps); rescale by true BR in analysis via `--xsec`. Output goes to `output/`.
2. **Signal acceptance** — `decayProbPerEvent_2body.py` ray-casts, computes decay probabilities, produces exclusion plots.
3. **Surface hit-map** — `signal_surface_hitmap_v2.py` maps accepted decays to tunnel surface coordinates.

## Main scripts

### `decayProbPerEvent_2body.py`
- Builds tunnel fiducial mesh from survey centerline.
- Ray-casts LLP directions from the IP to get entry/exit distances.
- Computes lifetime-dependent decay probability with 2-body acceptance cuts.
- Writes `particle_decay_results_2body.csv`, `event_decay_statistics_2body.csv`, and exclusion/separation plots.

### `signal_surface_hitmap_v2.py`
- Uses the same CSV input format (`event, id, pt, eta, phi, momentum, mass`).
- Maps accepted decays to tunnel surface coordinates.
- Produces surface coverage and efficiency/cost optimization plots.

### `backgroundFullGeo.py` (optional)
- Full-geometry muon-trident background Monte Carlo.
- Not part of the default `howto.md` run path, but kept for dedicated background studies.

## Geometry and physics constants

Shared tunnel model parameters are defined in the Python scripts via:
- `TUNNEL_ALPHA`, `TUNNEL_BETA`, `TUNNEL_GAMMA`, `TUNNEL_DELTA`
- `DETECTOR_THICKNESS`
- Centerline list `correctedVert` (shifted so the IP is at `[0,0,0]`)

## Dependencies

Python packages:
- `numpy`
- `pandas`
- `scipy`
- `trimesh`
- `matplotlib`
- `tqdm`

Build/runtime:
- Pythia8 (used by `pythiaStuff/main144.cc` and `pythiaStuff/make.sh`)
- ROOT is not required for the default non-ROOT build path.
