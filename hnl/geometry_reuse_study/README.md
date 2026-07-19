# Geometry Reuse Study

Read-only audit for the proposed FONLL-band shortcut: reuse the central
geometry cache for variation runs.

Conclusion: do not force row-for-row reuse of the monolithic central geometry
cache for FONLL variations. Two completed scale variations show large hit-mask
disagreement relative to their own exact ray-cast:

- `scale_muR0p5_muF0p5`: median `xor_hits / variation_hits = 0.896`
- `scale_muR0p5_muF1`: median `xor_hits / variation_hits = 0.910`

The regenerated FONLL channels are the problem:

- `Bmeson`, `Dmeson`, `Bbaryon`: median mismatch is about `1.8-2.0` per
  variation hit.
- Frozen/hardlinked channels checked here (`Bc`, `WZ`) match exactly.

Interpretation: the central and varied samples have similar aggregate
acceptance, but not the same row-by-row trajectories. Blindly copying central
geometry into variation runs would miss or invent many detector hits.

Safe optimization path: make geometry channel-aware. Reuse central geometry
only for hardlinked frozen channels (`Bc`, `Kmeson`, `tau`, `WZ`) and compute
geometry normally for regenerated FONLL channels (`Bmeson`, `Dmeson`,
`Bbaryon`, `induced_tau`). On the current band grid, frozen channels are about
half the rows overall (`Ue/Umu` median `0.50`, `Utau` median `0.67`), so this
is a useful geometry-stage speedup, but not a minutes-level replacement for the
exact band.

Files:

- `geometry_reuse_audit.py`: compares central and variation geometry caches.
- `scale_muR0p5_muF0p5_geometry_reuse.csv`
- `scale_muR0p5_muF1_geometry_reuse.csv`
- `frozen_channel_row_fractions.csv`
