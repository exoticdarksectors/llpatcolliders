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
bash /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/alps/generator/build.sh
```

---

## ALP production & analysis

### Production

```fish
set ALPS /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/alps

# Optional reset:
bash $ALPS/generator/clean.sh --all

# Heavy ALP (h→aa, m=15 GeV)
# Argument 2 is target generated events (not LLP rows).
# Argument 4 is max parallel generator jobs.
bash $ALPS/generator/produce.sh $ALPS/generator/heavy_alp.cmnd 50000 alp_heavy_m15 4

# Light ALP (B→K(*)a, m=1 GeV)
# Argument 5 is events per Pythia job (default auto).
bash $ALPS/generator/produce.sh $ALPS/generator/light_alp.cmnd 50000 alp_light_m1 4 50000

# Output structure:
#   data/   — LLP CSV, meta JSON
#   images/ — plots from analysis scripts
```

### Analysis

```fish
set ALPS /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/alps
cd $ALPS

# Heavy ALP (σ(gg→h) ≈ 54.7 pb N3LO 14 TeV = 54700 fb)
conda run -n llpatcolliders python decayProbPerEvent_2body.py \
    output/data/alp_heavy_m15.csv --xsec 54700 --lumi 3000 --outdir output
conda run -n llpatcolliders python signal_surface_hitmap_v2.py \
    output/data/alp_heavy_m15.csv --outdir output

# Light ALP (σ(pp→bb̄, inclusive) ≈ 373 µb = 3.73e8 fb)
conda run -n llpatcolliders python decayProbPerEvent_2body.py \
    output/data/alp_light_m1.csv --xsec 373000000 --lumi 3000 --outdir output
conda run -n llpatcolliders python signal_surface_hitmap_v2.py \
    output/data/alp_light_m1.csv --outdir output

# Geometry validation (after any change to gargoyle_geometry.py):
conda run -n llpatcolliders python visualize_tunnel.py
```

### Cross-sections (HL-LHC, sqrt(s) = 14 TeV, L = 3000 fb⁻¹)

| Channel | Cross-section | `--xsec` value |
|---------|--------------|-----------------|
| Heavy ALP (gg→h) | 54.7 pb (N3LO) | `54700` |
| Light ALP (pp→bb̄) | ~373 µb | `373000000` |

### Optional CLI args for `decayProbPerEvent_2body.py`

| Argument | Default | Description |
|----------|---------|-------------|
| `--lifetime-min-ns` | 0.01 | Minimum lifetime for scan (ns) |
| `--lifetime-max-ns` | 3.16e5 | Maximum lifetime for scan (ns) |
| `--lifetime-points` | 80 | Number of scan points |
| `--n-events` | auto | Total generated events (read from `_meta.json` if available) |
