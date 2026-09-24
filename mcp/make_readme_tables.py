#!/usr/bin/env python3
"""Regenerate the three result tables in README.md from the analysis output."""
import json

G   = {o["m"]: o for o in json.load(open("out/results_gam.json"))}
Z   = {o["m"]: o for o in json.load(open("out/results_gmZ.json"))}
peZ = {r["m"]: r for r in json.load(open("out/per_event_gmZ.json"))}
M   = sorted(Z)
REF = 0.02384

main = ("| m [GeV] | sigma gmZ [fb] | **R gmZ** | sigma gam [fb] | R gam | "
        "sigma ratio | R ratio |\n|---|---|---|---|---|---|---|\n" + "\n".join(
  f"| {m} | {Z[m]['sigma_fb']:.4g} | {Z[m]['R']:.4f} +- {Z[m]['R_err']:.4f} | "
  f"{G[m]['sigma_fb']:.4g} | {G[m]['R']:.4f} +- {G[m]['R_err']:.4f} | "
  f"{Z[m]['sigma_fb']/G[m]['sigma_fb']:.3f} | {Z[m]['R']/G[m]['R']:.3f} |" for m in M))

pe = ("| m [GeV] | R per particle | <N_acc>/event | P(>=1 chi) | P(both chi) | "
      "P(both)/R_pp^2 |\n|---|---|---|---|---|---|\n" + "\n".join(
  f"| {m} | {peZ[m]['R']:.4f} | {peZ[m]['Nacc']:.5f} | {peZ[m]['pAny']:.5f} | "
  f"{peZ[m]['pBoth']:.3g} | {peZ[m]['corr']:.2f} |" for m in M)
  + f"\n| *isotropic* | 1.0000 | {2*REF:.5f} | {2*REF-REF**2:.5f} | {REF**2:.3g} | 1.00 |")

vel = ("| m [GeV] | bg 2% | bg 16% | bg med | bg 84% | beta med | f(bg<1) | "
       "f(bg<2) | dt/10m med | dt/10m 98% |\n"
       "|---|---|---|---|---|---|---|---|---|---|\n" + "\n".join(
  f"| {m} | {Z[m]['bg_q'][0]:.2f} | {Z[m]['bg_q'][1]:.2f} | {Z[m]['bg_q'][2]:.2f} | "
  f"{Z[m]['bg_q'][3]:.2f} | {Z[m]['beta_q'][2]:.3f} | {Z[m]['f_bg_lt1']:.2f} | "
  f"{Z[m]['f_bg_lt2']:.2f} | {Z[m]['delay_med']:.1f} ns | {Z[m]['delay_98']:.0f} ns |"
  for m in M))

s = open("README.md").read()
def swap(text, anchor, table):
    i = text.index(anchor); j = text.index("|", i); k = text.index("\n\n", j)
    return text[:j] + table + text[k:]
s = swap(s, "## Results", main)
s = swap(s, "## Per-particle vs per-event", pe)
s = swap(s, "## Velocity / time of flight", vel)
open("README.md", "w").write(s)
print(f"README tables regenerated for {len(M)} masses")
