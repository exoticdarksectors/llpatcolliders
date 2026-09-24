#!/usr/bin/env python3
"""Figure: chi single-particle eta spectra vs A(eta), R(m), and the
acceptance-weighted beta*gamma distribution."""
import numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# All generated masses drive the R panel; a readable subset drives the
# spectrum and velocity panels.
MASSES = sorted(res_all_masses := [o["m"] for o in
                json.load(open("out/results_gmZ.json"))])
SHOW   = [10, 20, 30, 40, 50, 100, 200, 500, 1000]
a = np.load(os.environ.get("MCP_AETA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "external", "A_eta.npy")))
ctr, A = a[0], a[1]
dEta = ctr[1] - ctr[0]
edges = np.append(ctr - dEta / 2., ctr[-1] + dEta / 2.)

res  = {o["m"]: o for o in json.load(open("out/results_gmZ.json"))}
resG = {o["m"]: o for o in json.load(open("out/results_gam.json"))}

# Mass is an ordered magnitude -> one hue, light to dark (never a rainbow).
ramp = plt.get_cmap("Blues")(np.linspace(0.42, 0.97, len(SHOW)))
INK, INK2, MUTED = "#1c1f23", "#4a5057", "#9aa1a8"
REF_C = "#b4531f"

def read(fn):
    meta, rows, full = {}, [], None
    for line in open(fn):
        if line.startswith("#"):
            p = line[1:].split()
            if len(p) >= 2: meta[p[0]] = p[1]
        elif line.startswith("FULLHIST"):
            full = np.array(line.split()[1:], dtype=float)
        elif line[0].isalpha():
            continue
        else:
            rows.append(line.split())
    return meta, np.array(rows, dtype=float), full

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "axes.linewidth": 0.8, "figure.facecolor": "white"})

fig = plt.figure(figsize=(11.5, 6.6))
gs = GridSpec(2, 2, height_ratios=[2.6, 1.0], hspace=0.08, wspace=0.26,
              left=0.07, right=0.98, top=0.90, bottom=0.09)
axF = fig.add_subplot(gs[0, 0])
axA = fig.add_subplot(gs[1, 0], sharex=axF)
gsR = GridSpec(2, 2, hspace=0.42, wspace=0.26, left=0.07, right=0.98,
               top=0.90, bottom=0.09)
axR = fig.add_subplot(gsR[0, 1])
axB = fig.add_subplot(gsR[1, 1])

# ---- panel 1: single-particle eta pdf, normalised over ALL eta ----------
etaFine = np.linspace(-3, 3, 400)
axF.plot(etaFine, 0.5 / np.cosh(etaFine) ** 2, color=REF_C, lw=2, ls="--",
         label="isotropic reference", zorder=5)
for col, m in zip(ramp, SHOW):
    meta, d, full = read(f"out/z_m{m}.txt")
    nAll = float(meta["nChiAll"])
    fh = np.linspace(-6, 6, 241)
    ctrF = 0.5 * (fh[1:] + fh[:-1])
    pdf = full / nAll / (fh[1] - fh[0])
    sel = np.abs(ctrF) < 3
    axF.plot(ctrF[sel], pdf[sel], color=col, lw=2, label=f"{m} GeV")
axF.axvspan(edges[0], edges[-1], color="#eef2f6", zorder=0)
axF.set_ylabel(r"$f(\eta)$  (per $\chi$, normalised over all $\eta$)")
axF.set_xlim(-3, 3); axF.set_ylim(0, 0.56)
axF.grid(axis="y", color="#e6e9ec", lw=0.7); axF.set_axisbelow(True)
axF.spines[["top", "right"]].set_visible(False)
axF.tick_params(labelbottom=False)
axF.legend(frameon=False, ncol=3, fontsize=7, loc="upper left",
           handlelength=1.6, title=r"$m_\chi$", title_fontsize=8)
axF.set_title(r"$q\bar q\to\gamma^*/Z\to\chi\bar\chi$ at 14 TeV: "
              r"single-particle $\eta$", loc="left", fontsize=10, color=INK)

# ---- panel 2: the supplied acceptance -----------------------------------
axA.step(ctr, A, where="mid", color=INK2, lw=1.5)
axA.fill_between(ctr, A, step="mid", color=INK2, alpha=0.16)
axA.set_xlabel(r"$\eta$"); axA.set_ylabel(r"$A(\eta)$")
axA.grid(axis="y", color="#e6e9ec", lw=0.7); axA.set_axisbelow(True)
axA.spines[["top", "right"]].set_visible(False)
axA.set_ylim(0, 0.18)

# ---- panel 3: R(m) ------------------------------------------------------
ms = np.array(MASSES, float)
Rv = np.array([res[m]["R"] for m in MASSES])
Re = np.array([res[m]["R_err"] for m in MASSES])
RvG = np.array([resG[m]["R"] for m in MASSES])
ReG = np.array([resG[m]["R_err"] for m in MASSES])
axR.errorbar(ms, RvG, yerr=ReG, color="#b4531f", lw=2, ls="--", marker="s",
             ms=3.5, capsize=0, zorder=2, label=r"$\gamma^*$ only")
axR.errorbar(ms, Rv, yerr=Re, color="#1f5fa8", lw=2, marker="o", ms=4.5,
             capsize=0, zorder=3, label=r"$\gamma^*/Z$  (hypercharge mixing)")
axR.legend(frameon=False, fontsize=8, loc="upper left", handlelength=2.0)
axR.axvline(45.59, color=MUTED, lw=1, ls=":")
axR.annotate(r"$m_Z/2$", (45.59, 0.66), fontsize=7.5, color=INK2,
             ha="center", va="top", rotation=90,
             textcoords="offset points", xytext=(-4, 0))
for m, r in zip(ms, Rv):
    if m in (20, 45, 1000):
        axR.annotate(f"{r:.3f}", (m, r), textcoords="offset points",
                     xytext=(0, 10), ha="center", fontsize=8, color=INK)
axR.set_xscale("log"); axR.set_xlabel(r"$m_\chi$ [GeV]")
axR.set_ylabel(r"$R(m)=\int\! A f\,d\eta\,/\,0.02384$")
axR.set_ylim(0, 0.72)
axR.minorticks_off()
axR.grid(color="#e6e9ec", lw=0.7); axR.set_axisbelow(True)
axR.spines[["top", "right"]].set_visible(False)
TK = [10, 20, 50, 100, 200, 500, 1000]
axR.set_xticks(TK); axR.set_xticklabels([str(m) for m in TK])
axR.set_title("Acceptance relative to the isotropic reference, "
              "with and without the Z",
              loc="left", fontsize=10, color=INK)

# ---- panel 4: acceptance-weighted beta*gamma CDF ------------------------
for col, m in zip(ramp, SHOW):
    meta, d, full = read(f"out/z_m{m}.txt")
    eta, bg = d[:, 0], d[:, 1]
    sel = np.abs(eta) < edges[-1]
    w = A[np.digitize(eta[sel], edges) - 1]
    x, ww = bg[sel], w
    i = np.argsort(x); x, ww = x[i], ww[i]
    axB.plot(x, np.cumsum(ww) / ww.sum(), color=col, lw=2)
axB.axvline(1.0, color=MUTED, lw=1, ls=":")
axB.set_xscale("log"); axB.set_xlim(0.05, 50); axB.set_ylim(0, 1)
axB.set_xlabel(r"$\beta\gamma$")
axB.set_ylabel("cumulative fraction\n(weighted by $A$)")
axB.grid(color="#e6e9ec", lw=0.7); axB.set_axisbelow(True)
axB.spines[["top", "right"]].set_visible(False)
axB.annotate(r"$m_\chi$ increasing", xy=(0.35, 0.62), xycoords="data",
             fontsize=8, color=INK2)
axB.annotate("", xy=(0.28, 0.75), xytext=(0.5, 0.55),
             arrowprops=dict(arrowstyle="->", color=INK2, lw=1))
axB.set_title(r"Velocity of accepted $\chi$", loc="left", fontsize=10,
              color=INK)

fig.savefig("out/mcp_dy_summary.png", dpi=180, bbox_inches="tight")
print("wrote out/mcp_dy_summary.png")
