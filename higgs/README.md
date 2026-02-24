# Tunnel LLP Detector — Analysis Scripts

Three scripts for evaluating the physics reach and background rates of a long-lived particle detector in the CMS PX56 drainage tunnel. All three share the same tunnel geometry framework: a 3D mesh built from surveyed centerline coordinates with an arch-and-wall cross-section profile, inset by 24 cm of detector thickness to define the fiducial air volume.

## Scripts

### `decayProbPerEvent_2body.py` — Signal acceptance

Computes the decay-in-flight probability and analysis acceptance for LLPs produced at the CMS IP and decaying to two-body final states (e⁺e⁻) inside the tunnel fiducial volume.

**Input:** A CSV file with columns `event, eta, phi, momentum, mass`, one row per LLP from whatever generator (Pythia, MadGraph, etc.).

**What it does:**

- Ray-casts each LLP direction from the IP against the fiducial mesh to find entry/exit distances.
- For a scan over proper lifetimes, integrates the differential decay probability along each LLP's path, weighted by the two-body acceptance at each decay point: minimum electron energy (600 MeV) and minimum pair separation (1 mm) at the exit wall.
- Produces a lifetime scan plot (decay probability × BR vs cτ) and comparison curves against MATHUSLA, CODEX-b, and ANUBIS from external CSV files.

**Usage:**
```
python decayProbPerEvent_2body.py
```
Edit `sample_csv` and the external file paths in `__main__` as needed. The lifetime scan range and number of points are configurable.

**Output:** `particle_decay_results_2body.csv` with per-lifetime acceptance, plus comparison plots.

---

### `background_trident.py` — Beam muon trident backgrounds

Monte Carlo estimation of the beam muon trident background rate (μN → μNe⁺e⁻ in air) using the full tunnel geometry.

**What it does:**

- Samples muon directions (η, φ) isotropically over the solid angle subtended by the tunnel as seen from the IP.
- Ray-casts each direction to find entry/exit points and true path length through the fiducial air volume.
- Weights each muon by the flux at its entry distance, anchored to the milliQan measurement (0.3 fb⁻¹ cm⁻² at 33 m, falling as 1/r²). Energy spectrum: E⁻²·⁷ exp(−E/200 GeV), threshold 15 GeV.
- Uses forced trident production: every muon that hits the tunnel produces a trident, weighted by the actual production probability P = n_air × σ_trident × L_path. This avoids the statistics problem of sampling a ~10⁻⁶ probability process.
- Samples pair kinematics (energy fraction, asymmetry, opening angle) and computes separation at the exit wall.
- Applies signal-like cuts: both electrons above 600 MeV, pair separation between 1 mm and SEP_MAX.

**Usage:**
```
python background_fullgeo_mc.py
```
Number of rays, random seed, and cuts are set at the top of `__main__`. Default is 500k rays.

**Output:** Diagnostic plots (angular distributions, energy spectra, separation distributions, cut flow) and printed summary of weighted background rate per fb⁻¹.

---

### `signal_surface_hitmap_v2.py` — Tracking coverage optimization

Determines which tunnel surfaces receive signal decay products, allowing cost-optimized placement of fine-grained tracking (triangular scintillator strips) vs simple veto panels.

**Input:** Same CSV format as the signal script, plus an optional lifetime argument.

**What it does:**

- Ray-casts LLPs from the IP and records the exit point on the fiducial boundary where the e⁺e⁻ pair hits the tunnel wall.
- Maps each exit point to (s, θ) coordinates: distance along the tunnel centerline and angle around the cross-section profile.
- Classifies exit points by surface (floor, right wall, arch/ceiling, left wall) using geometry-derived boundary angles at the springline and floor corners.
- Computes contiguous angular coverage efficiency: for each possible contiguous arc of the tunnel profile, finds the optimal arc that captures a target fraction of signal. This gives physically meaningful answers like "a 95° arc captures 80% of the signal" rather than cherry-picked non-contiguous bins.
- Estimates tracking cost based on instrumented area, scaled from the full-coverage estimate of ~$7M for 500 m² with 4 layers of triangular strips.

**Usage:**
```
python signal_surface_hitmap_v2.py LLP.csv              # path-length weighting
python signal_surface_hitmap_v2.py LLP.csv 1e-7         # decay probability at τ = 100 ns
```

**Output:** A multi-panel plot (`*_surface_hitmap_v2.png`) with signal exit and entry heat maps, the 80% efficiency contiguous region overlay, efficiency vs area curve, surface fraction bar chart, cross-section scatter, and 1D angular distribution. Printed tables give the area and cost for 50–99% efficiency targets.

## Shared geometry

All three scripts use the same tunnel definition:

- **Centerline:** 47 survey points (`correctedVert`) defining the PX56 drainage tunnel path in CMS coordinates, shifted so the IP is at the origin.
- **Cross-section:** Arch-on-wall profile with 2.90 m floor width, 1.25 m wall height, and 1.90 m arch height (3.15 m total). Parameterized by `TUNNEL_ALPHA/BETA/GAMMA/DELTA`.
- **Fiducial volume:** The cross-section inset by `DETECTOR_THICKNESS` (24 cm) on all sides except the floor, extruded along the full centerline. Built as a `trimesh` mesh for ray-casting.
- **IP position:** Origin (0, 0, 0) after the coordinate shift applied to the survey data.

## Dependencies

`numpy`, `pandas`, `scipy`, `trimesh`, `matplotlib`, `tqdm`
