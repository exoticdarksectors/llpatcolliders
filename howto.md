# How-to: full run sequence

## One-time: build Pythia8 and main144

```bash
# Build Pythia8 (adjust path and ROOT prefix as needed)
cd /path/to/pythia8315
./configure --with-root=$(root-config --prefix)
make -j$(sysctl -n hw.logicalcpu)   # macOS; use nproc on Linux

# Build main144 against Pythia8
cd /path/to/llpatcolliders/pythiaStuff
bash make.sh
```

## Production (from pythiaStuff/)

```bash
bash parallel_produce.sh higgsLL.cmnd  10000 alp_heavy_m15
bash parallel_produce.sh alp_meson.cmnd 10000 alp_light_m1

# cleanup partial runs before restarting:
bash clean_production.sh --all
```

## Analysis (from repo root)

```bash
conda activate llpatcolliders

python decayProbPerEvent_2body.py output/alp_heavy_m15.csv --xsec 60000
python decayProbPerEvent_2body.py output/alp_light_m1.csv

python signal_surface_hitmap_v2.py output/alp_heavy_m15.csv
python signal_surface_hitmap_v2.py output/alp_light_m1.csv

# geometry validation (run once after any change to gargoyle_geometry.py):
python visualize_tunnel.py
```
