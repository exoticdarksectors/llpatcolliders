"""Sensitivity scan: horseshoe vs pentagon cross-section, 1 cm vs 4 cm bars.
Bar width -> SEP_MIN (min resolvable pair separation): 1 cm -> 0.01 m, 4 cm -> 0.04 m.
Output: excluded BR(H->LLP ee) vs proper decay length ctau, for 0.5 and 15 GeV."""
import numpy as np, trimesh
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

import decayProbPerEvent_2body as S
from grendel_geometry import (mesh_fiducial, path_3d_fiducial, create_profile_mesh,
    tunnel_profile_points, cache_geometry, DETECTOR_THICKNESS, TUNNEL_ALPHA,
    TUNNEL_WALL_HEIGHT)

C_LIGHT=2.99792458e8
ORIGIN=[0,0,0]

# ---- pentagon (inner/fiducial) mesh, inscribed in the arch at same location ----
fid=tunnel_profile_points(inset=DETECTOR_THICKNESS)
hw=TUNNEL_ALPHA/2-DETECTOR_THICKNESS
yf=fid[:,1].min(); yt=fid[:,1].max(); ysp=yf+TUNNEL_WALL_HEIGHT
pent=np.array([[-hw,yf],[hw,yf],[hw,ysp],[0,yt],[-hw,ysp]])
pv,pf=create_profile_mesh(path_3d_fiducial, pent)
mesh_pent=trimesh.Trimesh(vertices=pv,faces=pf); mesh_pent.fix_normals()
print("horseshoe volume %.1f m3 | pentagon volume %.1f m3" %
      (mesh_fiducial.volume, mesh_pent.volume))

MASSES={"0.5 GeV":"LLP0p5GeVSmall.csv", "15 GeV":"LLPSmall.csv"}
GEOMS ={"horseshoe":mesh_fiducial, "pentagon":mesh_pent}
SEPS  ={"1 cm bar":0.01, "4 cm bar":0.04}
lifetimes=np.logspace(-10.5,-3.5,16)
ctau=C_LIGHT*lifetimes

# cache ray-cast once per (mass,geom); scan both sep_min on each cache
results={}
for mlabel,csv in MASSES.items():
    for glabel,mesh in GEOMS.items():
        print("\n== caching %s / %s ==" % (mlabel,glabel))
        gc=cache_geometry(csv, mesh, ORIGIN)
        for slabel,sep in SEPS.items():
            r=S.analyze_decay_vs_lifetime(csv, gc, lifetimes,
                                          sep_min=sep, sep_max=S.SEP_MAX)
            results[(mlabel,glabel,slabel)]=np.array(r['exclusion'])
            print("   %s: min excl BR = %.2e" % (slabel, np.min(r['exclusion'])))

np.savez("sens_pentagon_scan.npz", ctau=ctau,
         **{f"{m}|{g}|{s}":results[(m,g,s)] for (m,g,s) in results})

# ---- plot: 2 panels (mass), 4 curves each (geom x bar) ----
styles={"horseshoe":"-","pentagon":"--"}
colors={"1 cm bar":"C0","4 cm bar":"C3"}
fig,axes=plt.subplots(1,2,figsize=(14,5.5),sharey=True)
for ax,mlabel in zip(axes,MASSES):
    for glabel in GEOMS:
        for slabel in SEPS:
            ax.plot(ctau, results[(mlabel,glabel,slabel)],
                    styles[glabel], color=colors[slabel], lw=2,
                    label=f"{glabel}, {slabel}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Proper decay length $c\tau$ (m)")
    ax.set_title(f"m(LLP) = {mlabel}")
    ax.grid(alpha=.3, which="both")
    ax.axhline(1.0, color="grey", ls=":", lw=1)
axes[0].set_ylabel(r"Excluded BR($H\to$ LLP $\to ee$)  [lower = better]")
axes[0].legend(fontsize=8, title="solid=horseshoe, dashed=pentagon")
fig.suptitle("GRENDEL sensitivity: horseshoe vs pentagon, 1 cm vs 4 cm bars "
             "(HL-LHC, 3 ab$^{-1}$, 3 events)")
fig.tight_layout()
fig.savefig("sens_pentagon_scan.png", dpi=130)
print("\nsaved sens_pentagon_scan.png")
