"""Horseshoe, 1 cm bars: full tunnel vs only X>-3 instrumented. Excluded BR vs ctau."""
import numpy as np, trimesh
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import decayProbPerEvent_2body as S
from grendel_geometry import (mesh_fiducial, path_3d_fiducial, create_profile_mesh,
    tunnel_profile_points, cache_geometry, DETECTOR_THICKNESS)
C=2.99792458e8; ORIGIN=[0,0,0]; SEP=0.01

# truncated path: interpolate exact X=-3 boundary then keep X>-3 nodes
p=path_3d_fiducial; X=p[:,0]
i0=np.where(X>-3.0)[0].min()                 # first node above -3 (node 24)
# interpolate boundary between i0-1 and i0
a,b=p[i0-1],p[i0]; t=(-3.0-a[0])/(b[0]-a[0])
bnd=a+t*(b-a)
trunc=np.vstack([bnd, p[i0:]])
prof=tunnel_profile_points(inset=DETECTOR_THICKNESS)
tv,tf=create_profile_mesh(trunc, prof); mesh_x=trimesh.Trimesh(vertices=tv,faces=tf); mesh_x.fix_normals()
print("full vol %.1f m3 | X>-3 vol %.1f m3 (%.0f%%)"%(mesh_fiducial.volume,mesh_x.volume,
      100*mesh_x.volume/mesh_fiducial.volume))

MASSES={"0.5 GeV":"LLP0p5GeVSmall.csv","15 GeV":"LLPSmall.csv","40 GeV":"LLP40GeVSmall.csv"}
GEOMS={"full":mesh_fiducial,"X>-3":mesh_x}
lifetimes=np.logspace(-10.5,-3.5,16); ctau=C*lifetimes
res={}
for ml,csv in MASSES.items():
    for gl,mesh in GEOMS.items():
        gc=cache_geometry(csv,mesh,ORIGIN)
        r=S.analyze_decay_vs_lifetime(csv,gc,lifetimes,sep_min=SEP,sep_max=S.SEP_MAX)
        res[(ml,gl)]=np.array(r['exclusion'])
        print("  %s %s: min excl BR %.2e"%(ml,gl,np.min(r['exclusion'])))

# retained fraction = full_reach / cut_reach at each mass optimum
fig,axes=plt.subplots(1,3,figsize=(16,5),sharey=True)
for ax,ml in zip(axes,MASSES):
    ax.plot(ctau,res[(ml,"full")],"-",color="C0",lw=2.2,label="full tunnel")
    ax.plot(ctau,res[(ml,"X>-3")],"--",color="C3",lw=2.2,label="X>-3 only (50% length)")
    r=np.min(res[(ml,"X>-3")])/np.min(res[(ml,"full")])
    ax.set_title("m=%s   (best-BR x%.2f worse)"%(ml,r))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=.3,which="both")
    ax.set_xlabel(r"$c\tau$ (m)"); ax.axhline(1,color="grey",ls=":")
axes[0].set_ylabel("Excluded BR(H->ee)  [lower=better]"); axes[0].legend(fontsize=9)
fig.suptitle("Horseshoe, 1 cm bars: full tunnel vs instrumenting only X>-3")
fig.tight_layout(); fig.savefig("sens_xcut.png",dpi=130)
print("saved sens_xcut.png")
