# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Physics analysis scripts for the **GARGOYLE** detector concept — a long-lived particle (LLP) detector proposed for the CMS PX56 drainage tunnel at CERN. The three main scripts evaluate signal acceptance, background rates, and tracking coverage optimization.

## Running the scripts

```bash
# Signal acceptance (lifetime scan, comparison to MATHUSLA/CODEX-b/ANUBIS)
python decayProbPerEvent_2body.py

# Background estimation (beam muon tridents)
python background_trident.py

# Tracking coverage / surface hit map
python signal_surface_hitmap.py LLP.csv              # path-length weighting
python signal_surface_hitmap.py LLP.csv 1e-7         # decay prob at τ = 100 ns
```

Dependencies: `numpy`, `pandas`, `scipy`, `trimesh`, `matplotlib`, `tqdm`

## CSV input format

All signal scripts read LLP events with columns:
```
event, id, pt, eta, phi, momentum, mass
```
Units: `pt` and `momentum` in GeV/c, `mass` in GeV/c², angles in radians (CMS convention).

## Architecture

### `gargoyle_geometry.py` — shared geometry module (imported by all three scripts)

The central module. Builds the 3D fiducial volume mesh on import and exposes it as `mesh_fiducial` and `path_3d_fiducial`. Key responsibilities:

- **Cross-section** (`tunnel_profile_points`): arch-on-wall profile parameterised by `TUNNEL_ALPHA/BETA/GAMMA/DELTA`. Inset by `DETECTOR_THICKNESS` (24 cm) to define the fiducial air gap.
- **Mesh extrusion** (`create_profile_mesh`): sweeps the 2D profile along the 3D centreline using Frenet-like frames (Y-up CMS convention).
- **Centreline** (`correctedVertWithShift`): 47 surveyed PX56 tunnel points converted from mm to metres and shifted so the CMS IP is at the origin. Tunnel runs at `Y_POSITION = 22 m` above the IP.
- **Ray-casting utilities**: `cache_geometry` (batch ray-cast + cache entry/exit distances), `eta_phi_to_direction` (η,φ → unit vector), `calculate_decay_length` (βγcτ).

Because the mesh is built at module import time (last line of the file), importing `gargoyle_geometry` takes a few seconds.

### `decayProbPerEvent_2body.py` — signal acceptance

Integrates decay-in-flight probability × two-body acceptance over each LLP's path through the fiducial volume, for a scan over proper lifetimes. The acceptance function enforces:
- minimum electron momentum (`P_CUT = 0.600 GeV`)
- minimum pair separation at the exit wall (`SEP_MIN = 1 mm`)
- maximum pair separation (`SEP_MAX = 1 m`)

Outputs `particle_decay_results_2body.csv` and comparison plots against external detector curves from `external/`.

**Key variable:** `outString` at the top of the file controls the output filename suffix.

### `background_trident.py` — muon trident background MC

Estimates beam-muon trident rate (μN → μNe⁺e⁻ in air) via forced production. Muon flux is anchored to the milliQan measurement (0.3 fb⁻¹ cm⁻² at 33 m, ∝ 1/r²). Energy spectrum is E⁻²·⁷ exp(−E/200 GeV) above 15 GeV. All geometry calls use `gargoyle_geometry` identically to the signal code.

### `signal_surface_hitmap.py` — tracking coverage optimisation

Maps exit points of e⁺e⁻ pairs onto (s, θ) coordinates (distance along centreline × angle around profile) and finds the minimum instrumented area needed to capture a target signal fraction. Two methods compared:
- **Fixed angular arc**: optimal contiguous arc of the cross-section profile (1D).
- **Adaptive 2D region**: greedy 8-connected region growing in (s, θ) space starting from the highest-signal bin.

Cost estimate: ~$7M for 500 m² with 4 layers of triangular scintillator strips (linear scaling).

## Coordinate convention

CMS standard throughout: X = horizontal transverse, Y = vertical (up), Z = beam axis. The tunnel centreline lies in the X-Z plane at fixed Y = 22 m. `eta_phi_to_direction` maps (η, φ) → (x, y, z) using the standard CMS θ = 2 arctan(e^{-η}) formula.

## External comparison files

`external/` contains sensitivity curves for MATHUSLA, CODEX-b, ANUBIS (nominal and optimistic) used in the lifetime scan plots of `decayProbPerEvent_2body.py`.
