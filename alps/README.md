# Tunnel LLP Detector — Active Workflow

This repository is currently centered on a Pythia generation + signal acceptance pipeline for LLPs in the CMS PX56 tunnel.

## Current workflow

See `howto.md` for the full command sequence.

Pipeline stages:
1. **Production** — `generator/produce.sh` runs `generator` in parallel batches until a target number of **generated events** is reached. Both cmnd files use BR=1.0 (efficiency maps); rescale by true BR in analysis and always pass explicit `--xsec` (do not use the script default). Output goes to `output/`.
2. **Signal acceptance** — `decayProbPerEvent_2body.py` ray-casts, computes decay probabilities with 2-body (a→μ⁺μ⁻) acceptance cuts, produces exclusion plots with mass-matched external comparison curves.
3. **Surface hit-map** — `signal_surface_hitmap_v2.py` maps accepted decays to tunnel surface coordinates with closed-form decay probability and multi-seed adaptive 2D region growing.

HL-LHC normalization notes (`sqrt(s)=14 TeV`, `L=3000 fb^-1`):
- Heavy (`h -> aa`): `--xsec 54700` (σ(gg→h) ≈ 54.7 pb, N3LO, 14 TeV). Exclusion curve gives BR(h→aa)_min.
- Light (`B -> K(*)a`): `--xsec 373000000` (σ(pp→bb̄, inclusive) ≈ 373 µb). Generator uses BR(B→Ka)=1, so the exclusion curve gives BR(B→Ka)_min directly.

## Main scripts

### `gargoyle_geometry.py`
- Single source of truth for all tunnel geometry (cross-section constants, centerline
  survey data, mesh builder, ray-casting utilities).
- Shared by all three analysis scripts — do not edit geometry in individual scripts.
- Builds the fiducial-volume trimesh on import; exports `mesh_fiducial` and
  `path_3d_fiducial`.

### `decayProbPerEvent_2body.py`
- Ray-casts LLP directions from the IP to get entry/exit distances.
- Computes lifetime-dependent decay probability with 2-body (μ⁺μ⁻) acceptance cuts.
- Extended lifetime scan: 10^-2 to 10^5.5 ns (80 points by default).
- Mass-matched external comparison curves via `overlay_mass_matched_external_curves()`:
  infers sample mass from CSV, loads only `external/<STEM>_m<tag>.csv` files at
  the matching mass.
- Writes `particle_decay_results_2body.csv`, `event_decay_statistics_2body.csv`, and
  exclusion/separation plots to `--outdir`.
- Key CLI args: `--xsec` (fb), `--lumi` (fb^-1), `--outdir`, `--n-events`,
  `--lifetime-min-ns`, `--lifetime-max-ns`, `--lifetime-points`.

### `signal_surface_hitmap_v2.py`
- Maps accepted decays to tunnel surface coordinates (arc-length s, profile angle θ).
- Uses vectorized closed-form decay probability (no 2-body acceptance model).
- Computes contiguous-arc efficiency **and** multi-seed adaptive 2D region-growing
  efficiency vs. instrumented area for detector optimization.
- Produces surface coverage, efficiency comparison (fixed arc vs adaptive 2D), and
  density cross-section plots.  Output saved to `--outdir`/images/.

### `backgroundFullGeo.py` (optional)
- Full-geometry muon-trident background Monte Carlo.
- Imports the same mesh and path from `gargoyle_geometry.py`.
- Not part of the default `howto.md` run path, but kept for dedicated background studies.

### `visualize_tunnel.py`
- Static 4-panel figure + interactive 3D plot of the tunnel geometry.
- Useful for geometry validation (run after any change to `gargoyle_geometry.py`).

## Geometry

All tunnel model parameters live in `gargoyle_geometry.py`:
- Cross-section constants: `TUNNEL_ALPHA/BETA/GAMMA/DELTA`, `DETECTOR_THICKNESS`
- Centerline: `correctedVertWithShift` (shifted so CMS IP5 is at `[0,0,0]`, Y=22 m above)
- Exports: `mesh_fiducial`, `path_3d_fiducial`

## Dependencies

Python packages:
- `numpy`
- `pandas`
- `scipy`
- `trimesh`
- `matplotlib`
- `tqdm`

Build/runtime:
- Pythia8 — `generator/build.sh` auto-resolves via `$PYTHIA8_DIR` or a sibling `<repo-root>/../pythia8315` directory.
- ROOT is not required for the default build path (CSV-only output).
