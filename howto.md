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
# bash $DP/generator/clean.sh --all

# h→A'A', three mass points (σ(pp→h) ≈ 60 pb = 60000 fb)
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m05.cmnd 10000 dp_heavy_m05 4 1500
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m1.cmnd  10000 dp_heavy_m1  4 1500
bash $DP/generator/produce.sh $DP/generator/heavy_dp_m15.cmnd 10000 dp_heavy_m15 4 1500
```

### Analysis

```fish
set DP /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/dark_photon
cd $DP

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
