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

# h→A'A', three mass points (σ(pp→h) ≈ 60 pb = 60000 fb)
# Argument 2 is target generated events (not LLP rows).
# Argument 4 is max parallel generator jobs (roughly max CPU cores used).
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m05.cmnd 10000 dp_heavy_m05 4
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m1.cmnd  10000 dp_heavy_m1  4
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m15.cmnd 10000 dp_heavy_m15 4

# Output structure:
#   data   -> LLP CSV, daughters CSV, meta JSON
#   images -> plots from analysis scripts
```

### Analysis

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon
cd $DP

for tag in dp_heavy_m05 dp_heavy_m1 dp_heavy_m15
    # N-track analysis with pair+separation selection:
    # >=2 charged daughters with p>600 MeV, union of all valid pair separation windows [1 mm, 1 m]
    conda run -n llpatcolliders python decayProbPerEvent_Ntrack.py \
        output/data/$tag.csv --xsec 60000 --lumi 3000 --outdir output/ \
        --sep-min 0.001 --sep-max 1.0
    conda run -n llpatcolliders python signal_surface_hitmap_v2.py \
        output/data/$tag.csv --outdir output/
end

# Geometry validation (after any change to gargoyle_geometry.py):
conda run -n llpatcolliders python visualize_tunnel.py
```
