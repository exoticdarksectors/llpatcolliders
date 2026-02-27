# Dark Photon (vector portal) — Active Workflow

Dark photon A': massive spin-1 vector boson with kinetic mixing epsilon with SM photon.

Two production channels:
- **Heavy A'**: exotic Higgs decay `pp → h(125) → A'A'` — mass benchmarks:
  - `generator/heavy_dp_m05.cmnd` — m_A' = 0.5 GeV
  - `generator/heavy_dp_m1.cmnd`  — m_A' = 1 GeV
  - `generator/heavy_dp_m15.cmnd` — m_A' = 15 GeV
- **Light A'** (0.2–5 GeV): B-meson FCNC `pp → bb̄ → B → K(*)A'` *(not yet implemented)*

Decay: kinetic-mixing BRs from the R-ratio (mass-dependent leptonic + hadronic channels).
Pythia8 handles full decay and hadronization. Generator writes charged daughters to `_daughters.csv`.

## Build

```fish
set -x PYTHIA8_DIR /Users/fredi/sandbox-offline/pythia8315
bash /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon/generator/build.sh
```

## Production

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon

# h→A'A', three mass points (σ(pp→h) ≈ 60 pb = 60000 fb)
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m05.cmnd 10000 dp_heavy_m05 4 1500
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m1.cmnd  10000 dp_heavy_m1  4 1500
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m15.cmnd 10000 dp_heavy_m15 4 1500
```

## Analysis (conda env llpatcolliders)

Always pass `--xsec` explicitly. Run from `dark_photon/` so `external/*.csv` resolves correctly.

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon
cd $DP

# HL-LHC baseline: sqrt(s)=14 TeV, L=3000 fb^-1
for tag in dp_heavy_m05 dp_heavy_m1 dp_heavy_m15
    # N-track analysis (≥2 charged tracks, matching MATHUSLA/ANUBIS/CODEX-b)
    conda run -n llpatcolliders python decayProbPerEvent_Ntrack.py \
        output/$tag.csv --xsec 60000 --lumi 3000 --outdir output/
    conda run -n llpatcolliders python signal_surface_hitmap_v2.py \
        output/$tag.csv --outdir output/
end

# Geometry validation (after any change to gargoyle_geometry.py):
conda run -n llpatcolliders python visualize_tunnel.py
```

## Key differences from ALP

| Property        | ALP (`6000113`)  | Dark photon (`6000115`) |
|-----------------|------------------|-------------------------|
| Spin            | 0 (pseudoscalar) | 1 (vector)              |
| Pythia spinType | 1                | 3                       |
| Particle name   | ALP              | Aprime                  |
| isResonance     | off              | off                     |
| Production      | same (h→XX)      | same                    |
| Decay           | μ⁺μ⁻ (BR=1)      | R-ratio BRs (full)      |
| Analysis        | 2-body analytical | ≥2 charged tracks (DV)  |

## External comparison curves

`external/` accepts digitised sensitivity curves from ANUBIS/CODEX-b/MATHUSLA dark photon papers.
The existing h→SS dark-Higgs curves (from the MATHUSLA/ANUBIS/CODEX-b proposals) are a valid
geometric proxy for h→A'A' at the same mass — same production mechanism, same boost distribution.
No dedicated h→A'A' curves at m=15 GeV have been published yet.
