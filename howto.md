cd /Users/fredi/sandbox-offline/pythia8315
./configure --with-root=/opt/homebrew/Cellar/root/6.38.00
make -j$(sysctl -n hw.logicalcpu)

cd /Users/fredi/sandbox-offline/llpatcolliders_MATT/llpatcolliders/pythiaStuff
bash make.sh

conda activate llpatcolliders

bash parallel_produce.sh higgsLL.cmnd 10000 alp_heavy_m15
bash parallel_produce.sh alp_meson.cmnd 10000 alp_light_m1

cd ..
python decayProbPerEvent_2body.py output/alp_heavy_m15.csv --xsec 60000
python decayProbPerEvent_2body.py output/alp_light_m1.csv
python signal_surface_hitmap_v2.py output/alp_heavy_m15.csv
python signal_surface_hitmap_v2.py output/alp_light_m1.csv

# cleanup if aborted / restarting:
bash clean_production.sh --all