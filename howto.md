# How-to: full run sequence

## One-time: build Pythia8 and main144

```fish
cd /Users/fredi/sandbox-offline/pythia8315
./configure --with-root=/opt/homebrew/Cellar/root/6.38.00
make -j(sysctl -n hw.logicalcpu)

cd /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/pythiaStuff
bash make.sh
```

## Production (from pythiaStuff/)

```fish
# Optional reset BEFORE starting a fresh production (destructive: removes output/*.csv):
# bash clean_production.sh --all

bash parallel_produce.sh higgsLL.cmnd  10000 alp_heavy_m15
bash parallel_produce.sh alp_meson.cmnd 10000 alp_light_m1
```

## Analysis (from repo root)

```fish
conda activate llpatcolliders

# Always pass --xsec explicitly (do not rely on decayProbPerEvent_2body.py default).
# HL-LHC baseline: sqrt(s)=14 TeV, L=3000 fb^-1.
set XSEC_HEAVY_FB 60000
# Light sample (alp_meson.cmnd): set benchmark-specific effective rate
# XSEC_LIGHT_FB = sigma(pp->bb, generator cuts) * P(B->K(*)a), in fb.
set XSEC_LIGHT_FB 52000  # placeholder only; replace for your coupling benchmark

python decayProbPerEvent_2body.py output/alp_heavy_m15.csv --xsec $XSEC_HEAVY_FB --lumi 3000
python decayProbPerEvent_2body.py output/alp_light_m1.csv --xsec $XSEC_LIGHT_FB --lumi 3000

python signal_surface_hitmap_v2.py output/alp_heavy_m15.csv
python signal_surface_hitmap_v2.py output/alp_light_m1.csv

# geometry validation (run once after any change to gargoyle_geometry.py):
python visualize_tunnel.py
```
