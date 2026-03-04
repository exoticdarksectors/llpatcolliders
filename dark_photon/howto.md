# How-to: full run sequence

## One-time: build Pythia8

```fish
cd /Users/fredi/sandbox-offline/pythia8315
./configure                          # add --with-root=... if you want ROOT support
make -j(sysctl -n hw.logicalcpu)
```

## One-time: build generator binary

```fish
set -x PYTHIA8_DIR /Users/fredi/sandbox-offline/pythia8315
bash /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon/generator/build.sh
```

---

## Dark photon production & analysis

### Production

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon

# Optional reset:
bash $DP/generator/clean.sh --all

# h→A'A', three mass points (σ(gg→h) ≈ 54.7 pb N3LO 14 TeV = 54700 fb)
# Argument 2 is target generated events (not LLP rows).
# Argument 4 is max parallel generator jobs (roughly max CPU cores used).
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m05.cmnd 50000 dp_heavy_m05 4
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m1.cmnd  50000 dp_heavy_m1  4
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m15.cmnd 50000 dp_heavy_m15 4

# Light A' via meson decay (m=0.5 GeV):
bash $DP/generator/produce.sh $DP/generator/meson_dp_omega_m05.cmnd 100000 dp_meson_omega_m05 4
bash $DP/generator/produce.sh $DP/generator/meson_dp_eta_m05.cmnd   100000 dp_meson_eta_m05   4

# Light A' via Drell-Yan (m=1.0 GeV):
bash $DP/generator/produce.sh $DP/generator/dy_dp_m1.cmnd 10000 dp_dy_m1 4

# Output structure:
#   data   -> LLP CSV, daughters CSV, meta JSON
#   images -> plots from analysis scripts
```

### Analysis (meson production — ε² scan)

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon
cd $DP

# σ_inel ≈ 80 mb = 8e10 fb
for parent in omega eta
    conda run -n llpatcolliders python decayProbPerEvent_Ntrack.py \
        output/data/dp_meson_{$parent}_m05.csv \
        --xsec 80000000000 --lumi 3000 --outdir output/ \
        --sep-min 0.001 --sep-max 1.0 \
        --production meson --parent-meson $parent --dp-mass 0.5
    conda run -n llpatcolliders python signal_surface_hitmap_v2.py \
        output/data/dp_meson_{$parent}_m05.csv --outdir output/
end
```

### Analysis (Drell-Yan DP — ε² scan)

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon
cd $DP

# σ(pp→A') at ε²=1 from Pythia8 log (grep "sum" output/data/.dp_dy_m1_part_*.log)
# Convert mb → fb (1 mb = 1e12 fb).
conda run -n llpatcolliders python decayProbPerEvent_Ntrack.py \
    output/data/dp_dy_m1.csv \
    --xsec <sigma_fb_from_log> --lumi 3000 --outdir output/ \
    --production drell_yan --dp-mass 1.0
```

### Analysis (heavy DP — Higgs production)

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon
cd $DP

for tag in dp_heavy_m05 dp_heavy_m1 dp_heavy_m15
    # N-track analysis with pair+separation selection:
    # >=2 charged daughters with p>600 MeV, union of all valid pair separation windows [1 mm, 1 m]
    conda run -n llpatcolliders python decayProbPerEvent_Ntrack.py \
        output/data/$tag.csv --xsec 54700 --lumi 3000 --outdir output/ \
        --sep-min 0.001 --sep-max 1.0
    conda run -n llpatcolliders python signal_surface_hitmap_v2.py \
        output/data/$tag.csv --outdir output/
end

# Geometry validation (after any change to gargoyle_geometry.py):
conda run -n llpatcolliders python visualize_tunnel.py
```
