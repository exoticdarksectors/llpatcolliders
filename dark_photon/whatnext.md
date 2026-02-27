 Ready to code?

 Here is Claude's plan:
╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
 Plan: Two new physics folders — dark Higgs & dark photon

 Context

 The repository currently covers ALPs (user, alps branch) and dark Higgs (Matt, higgs branch —
 old code, pre-modern pipeline). The four canonical portals benchmarked by transverse LLP
 detectors (ANUBIS, CODEX-b, MATHUSLA) at the HL-LHC are:

 ┌──────────────┬─────────────────┬────────────────────┐
 │    Portal    │      Model      │        Plan        │
 ├──────────────┼─────────────────┼────────────────────┤
 │ Scalar       │ Dark Higgs h→SS │ dark_higgs/ below  │
 ├──────────────┼─────────────────┼────────────────────┤
 │ Pseudoscalar │ ALP h→aa, B→Ka  │ Done (alps/)       │
 ├──────────────┼─────────────────┼────────────────────┤
 │ Vector       │ Dark photon A'  │ dark_photon/ below │
 ├──────────────┼─────────────────┼────────────────────┤
 │ Fermion      │ HNLs            │ Excluded           │
 └──────────────┴─────────────────┴────────────────────┘

 Both new folders mirror the alps/ structure exactly. Both reuse generator.cc unchanged
 (PDG ID is already configurable via LLP:pdgId in the cmnd file). All Python analysis
 scripts are reused as-is. Only new cmnd files and shell script copies are needed.

 ---
 ---
 Part 1 — dark_higgs/ (scalar portal)

 Physics

 Dark Higgs S: CP-even scalar mixing with SM Higgs via portal coupling sin θ.
 Two production channels mirror the ALP structure:
 - Heavy S (1–60 GeV): exotic Higgs decay pp → h(125) → S S
 - Light S (0.2–5 GeV): B-meson FCNC pp → bb̄ → B → K(*)

 Decay: S → μ⁺μ⁻ (BR=1 for efficiency maps, rescale in analysis).
 Cross-sections: same as corresponding ALP channels (60 pb heavy, 0.37 mb light).
 External comparison curves in alps/external/ are valid for this exact model
 (MATHUSLA/CODEX-b/ANUBIS for H(125)→SS at m_S = 1 and 15 GeV — Matt digitised them).

 Target layout

 dark_higgs/
 ├── generator/
 │   ├── heavy_dh.cmnd    ← NEW: h→SS + S→μ+μ-, PDG 6000114, spin-0 scalar
 │   ├── light_dh.cmnd    ← NEW: B→K(*)S + S→μ+μ-, PDG 6000114
 │   ├── build.sh         ← copy of alps/generator/build.sh
 │   ├── produce.sh       ← copy of alps/generator/produce.sh
 │   └── clean.sh         ← copy of alps/generator/clean.sh
 ├── output/
 │   └── .gitkeep
 └── README.md

 heavy_dh.cmnd — diff from alps/generator/heavy_alp.cmnd

 ┌─────────────┬──────────────────────────────────────────────┬──────────────────────────────────────────────────────┐
 │    Item     │                heavy_alp.cmnd                │                    heavy_dh.cmnd                     │
 ├─────────────┼──────────────────────────────────────────────┼──────────────────────────────────────────────────────┤
 │ PDG ID      │ 6000113                                      │ 6000114                                              │
 ├─────────────┼──────────────────────────────────────────────┼──────────────────────────────────────────────────────┤
 │ Particle    │ ALP ALP                                      │ S S                                                  │
 │ name        │                                              │                                                      │
 ├─────────────┼──────────────────────────────────────────────┼──────────────────────────────────────────────────────┤
 │ LLP:pdgId   │ missing (relies on default 6000113)          │ LLP:pdgId = 6000114                                  │
 ├─────────────┼──────────────────────────────────────────────┼──────────────────────────────────────────────────────┤
 │ h→SS        │ 25:addChannel = 1 0.0 0 6000113 6000113 +    │ 25:addChannel = 1 1.0 100 6000114 6000114 (BR and    │
 │ channel     │ fragile 25:76:bRatio index hack              │ meMode inline — avoids index fragility)              │
 ├─────────────┼──────────────────────────────────────────────┼──────────────────────────────────────────────────────┤
 │ 25:onIfAny  │ = 6000113                                    │ = 6000114                                            │
 ├─────────────┼──────────────────────────────────────────────┼──────────────────────────────────────────────────────┤
 │ Decay       │ 6000113:oneChannel = 1 1.0 0 13 -13          │ 6000114:oneChannel = 1 1.0 0 13 -13                  │
 └─────────────┴──────────────────────────────────────────────┴──────────────────────────────────────────────────────┘

 spinType = 1 (spin-0), same as ALP. CP-parity distinction does not affect Pythia8 kinematics.

 light_dh.cmnd — diff from alps/generator/light_alp.cmnd

 Replace every 9000001 with 6000114 throughout. Replace particle name ALP with S.
 All 8 B-meson channels (B±, B⁰, B̄⁰ → K, K* + S) identical otherwise

 Also: fix alps/generator/heavy_alp.cmnd

 Add LLP:pdgId = 6000113 explicitly (1 line — currently relies on default).

 Analysis invocation

 bash dark_higgs/generator/build.sh

 bash dark_higgs/generator/produce.sh heavy_dh.cmnd 10000 dh_heavy_m15
 bash dark_higgs/generator/produce.sh light_dh.cmnd 10000 dh_light_m1 8 50000

 # Reuse alps/ analysis scripts — external curves in alps/external/ are valid here
 conda run -n llpatcolliders python alps/decayProbPerEvent_2body.py \
     dark_higgs/output/dh_heavy_m15.csv \
     --xsec 60000 --lumi 3000 --outdir dark_higgs/output/

 conda run -n llpatcolliders python alps/decayProbPerEvent_2body.py \
     dark_higgs/output/dh_light_m1.csv \
     --xsec 373000000 --lumi 3000 --outdir dark_higgs/output/

 ---
 ---
 Part 2 — dark_photon/ (vector portal)

 Physics

 Dark photon A': massive spin-1 vector boson with kinetic mixing ε with SM photon.
 Minimal vector extension of the SM; benchmarked in every ANUBIS/CODEX-b/MATHUSLA paper.
 Two production channels mirror the ALP/dark Higgs structure:
 - Heavy A' (1–60 GeV): exotic Higgs decay pp → h(125) → A' A'
 - Light A' (0.2–5 GeV): B-meson FCNC pp → bb̄ → B → K(*) A

 Decay: A' → μ⁺μ⁻ (BR=1 for efficiency maps).
 Cross-sections: same as corresponding ALP channels.
 External comparison curves: need to be digitised from ANUBIS/CODEX-b/MATHUSLA dark
 photon sensitivity papers; placeholder external/.gitkeep created.

 Key difference from ALPs/dark Higgs: A' is spin-1 (spinType = 3 in Pythia8).
 The B→KA' matrix element differs (vector vs pseudoscalar recoil), but the A' kinematics
 (pt, eta, phi, |p|, m) captured in the CSV are unaffected by this at the acceptance level.

 Target layout

 dark_photon/
 ├── generator/
 │   ├── heavy_dp.cmnd    ← NEW: h→A'A', PDG 6000115, spinType=3
 │   ├── light_dp.cmnd    ← NEW: B→K(*)A', PDG 6000115
 │   ├── build.sh         ← copy of alps/generator/build.sh
 │   ├── produce.sh       ← copy of alps/generator/produce.sh
 │   └── clean.sh         ← copy of alps/generator/clean.sh
 ├── output/
 │   └── .gitkeep
 ├── external/
 │   └── .gitkeep
 └── README.md

 heavy_dp.cmnd — diff from alps/generator/heavy_alp.cmnd

 ┌────────────────┬─────────────────────────────────────┬───────────────────────────────────────────┐
 │      Item      │           heavy_alp.cmnd            │               heavy_dp.cmnd               │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ PDG ID         │ 6000113                             │ 6000115                                   │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ Particle name  │ ALP ALP                             │ Aprime Aprime                             │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ Spin           │ spinType = 1                        │ spinType = 3 (vector)                     │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ LLP:pdgId      │ missing                             │ LLP:pdgId = 6000115                       │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ h→A'A' channel │ 25:addChannel... + index hack       │ 25:addChannel = 1 1.0 100 6000115 6000115 │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ 25:onIfAny     │ = 6000113                           │ = 6000115                                 │
 ├────────────────┼─────────────────────────────────────┼───────────────────────────────────────────┤
 │ Decay          │ 6000113:oneChannel = 1 1.0 0 13 -13 │ 6000115:oneChannel = 1 1.0 0 13 -13       │
 └────────────────┴─────────────────────────────────────┴───────────────────────────────────────────┘

 light_dp.cmnd — diff from alps/generator/light_alp.cmnd

 Replace every 9000001 with 6000115. Replace particle name ALP with Aprime.
 Add spinType = 3 to particle definition (vector boson). All B-meson channels otherwise identical.

 Analysis invocation

 bash dark_photon/generator/build.sh

 bash dark_photon/generator/produce.sh heavy_dp.cmnd 10000 dp_heavy_m15
 bash dark_photon/generator/produce.sh light_dp.cmnd 10000 dp_light_m1 8 50000

 conda run -n llpatcolliders python alps/decayProbPerEvent_2body.py \
     dark_photon/output/dp_heavy_m15.csv \
     --xsec 60000 --lumi 3000 --outdir dark_photon/output/

 conda run -n llpatcolliders python alps/decayProbPerEvent_2body.py \
     dark_photon/output/dp_light_m1.csv \
     --xsec 373000000 --lumi 3000 --outdir dark_photon/output/

 ---
 ---
 Summary of all file changes

 ┌─────────────────────────────────────┬──────────────────────────────┐
 │                File                 │            Action            │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/generator/heavy_dh.cmnd  │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/generator/light_dh.cmnd  │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/generator/build.sh       │ NEW (copy)                   │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/generator/produce.sh     │ NEW (copy)                   │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/generator/clean.sh       │ NEW (copy)                   │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/output/.gitkeep          │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_higgs/README.md                │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/generator/heavy_dp.cmnd │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/generator/light_dp.cmnd │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/generator/build.sh      │ NEW (copy)                   │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/generator/produce.sh    │ NEW (copy)                   │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/generator/clean.sh      │ NEW (copy)                   │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/output/.gitkeep         │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/external/.gitkeep       │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ dark_photon/README.md               │ NEW                          │
 ├─────────────────────────────────────┼──────────────────────────────┤
 │ alps/generator/heavy_alp.cmnd       │ +1 line: LLP:pdgId = 6000113 │
 └─────────────────────────────────────┴──────────────────────────────┘

 generator.cc — no changes.
 All Python analysis scripts — no changes.

 ---
 Verification (both folders)

 # dark_higgs smoke test
 bash dark_higgs/generator/produce.sh heavy_dh.cmnd 200 test_dh 2 500
 cat dark_higgs/output/test_dh_meta.json      # llp_pdg_id must be 6000114

 # dark_photon smoke test
 bash dark_photon/generator/produce.sh heavy_dp.cmnd 200 test_dp 2 500
 cat dark_photon/output/test_dp_meta.json     # llp_pdg_id must be 6000115

 # Analysis on both
 conda run -n llpatcolliders python alps/decayProbPerEvent_2body.py \
     dark_higgs/output/test_dh.csv --xsec 60000 --lumi 3000 \
     --outdir dark_higgs/output/

 conda run -n llpatcolliders python alps/decayProbPerEvent_2body.py \
     dark_photon/output/test_dp.csv --xsec 60000 --lumi 3000 \
     --outdir dark_photon/output/