#!/usr/bin/env python3
"""Compare the reference DY cross sections with this generation."""
import numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ref = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), "external", "dy_xsecs_eta1_RUN3.txt"), delimiter=",")
rm, rs, ra = ref[:, 0], ref[:, 1], ref[:, 2]
Z = {o["m"]: o for o in json.load(open("out/results_gmZ.json"))}
M = [m for m in sorted(Z) if m <= 79]

INK, INK2, MUTED = "#1c1f23", "#4a5057", "#9aa1a8"
BLUE, RUST, TEAL = "#1f5fa8", "#b4531f", "#2c7a70"
MZ2 = 91.1876 / 2.

def load(fn):
    meta, pairs = {}, []
    for L in open(fn):
        if L.startswith("#"):
            p = L[1:].split()
            if len(p) >= 2: meta[p[0]] = p[1]
        elif L.startswith("PAIR"):
            pairs.append(L.split()[1:])
    return meta, np.array(pairs, dtype=float)

incl, cut1, refv, refa, accv = [], [], [], [], []
for m in M:
    meta, e = load(f"out/z_m{m}.txt")
    nEv = float(meta["nEventAccepted"]); s = float(meta["sigma_mb"]) * 1e9
    acc = (np.abs(e[:, 0]) < 1.) | (np.abs(e[:, 1]) < 1.)
    incl.append(s); cut1.append(s * acc.sum() / nEv)
    accv.append(acc.sum() / nEv)
    refv.append(np.interp(m, rm, rs)); refa.append(np.interp(m, rm, ra))
M = np.array(M, float); incl = np.array(incl); cut1 = np.array(cut1)
refv = np.array(refv); refa = np.array(refa); accv = np.array(accv)

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "axes.linewidth": 0.8, "figure.facecolor": "white"})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.2, 4.2))
fig.subplots_adjust(left=0.07, right=0.99, top=0.87, bottom=0.14, wspace=0.24)

sel = (rm >= 8) & (rm <= 79)
ax.semilogy(rm[sel], rs[sel], color=INK2, lw=2.2, label="reference file")
ax.semilogy(M, incl, color=BLUE, lw=2, ls="--", marker="o", ms=4,
            label=r"this work, inclusive")
ax.semilogy(M, cut1, color=RUST, lw=2, ls="--", marker="s", ms=4,
            label=r"this work, $\times$ P($\geq 1\ \chi$ at $|\eta|<1$)")
ax.semilogy(rm[sel], rs[sel] * np.interp(rm[sel], rm, ra), color=TEAL, lw=2,
            label=r"reference, col 2 $\times$ col 3")
ax.axvline(MZ2, color=MUTED, lw=1, ls=":")
ax.set_xlabel(r"$m_\chi$ [GeV]"); ax.set_ylabel(r"$\sigma$ [pb],  $Q_\chi = e$")
ax.set_xlim(8, 79); ax.set_ylim(0.3, 6e3)
ax.grid(color="#e6e9ec", lw=0.7); ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8.5, loc="lower left", handlelength=2.2)
ax.set_title("Cross sections", loc="left", fontsize=10, color=INK)

bx.plot(M, refv / incl, color=BLUE, lw=2, marker="o", ms=5,
        label=r"col 2 / my inclusive $\sigma$")
bx.plot(M, refa / accv, color=TEAL, lw=2, marker="D", ms=4.5,
        label=r"col 3 / my P($\geq 1\ \chi$ at $|\eta|<1$)")
bx.axvline(MZ2, color=MUTED, lw=1, ls=":")
bx.axhline(1.0, color=MUTED, lw=1, ls=":")
bx.set_ylim(0.85, 1.35)

bx.set_xlabel(r"$m_\chi$ [GeV]"); bx.set_ylabel("ratio")
bx.set_xlim(8, 55)
bx.grid(color="#e6e9ec", lw=0.7); bx.set_axisbelow(True)
bx.spines[["top", "right"]].set_visible(False)
bx.legend(frameon=False, fontsize=8.5, loc="upper left", handlelength=2.2)
bx.annotate("normalisation only (LO vs NLO)", (32, 1.20),
            textcoords="offset points", xytext=(0, 8), ha="center",
            fontsize=8, color=BLUE)
bx.annotate(r"acceptance agrees to 2%", (32, 1.0),
            textcoords="offset points", xytext=(0, -14), ha="center",
            fontsize=8, color=TEAL)
bx.set_title("Column 3 is P($\\geq$1 $\\chi$ in $|\\eta|<1$)", loc="left",
             fontsize=10, color=INK)
fig.savefig("out/mcp_dy_ref_compare.png", dpi=180, bbox_inches="tight")
print("wrote out/mcp_dy_ref_compare.png")
