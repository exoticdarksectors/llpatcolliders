#!/usr/bin/env python3
"""Zoom on the 35-50 GeV scan, where the chi chi threshold crosses the Z pole."""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Z = {o["m"]: o for o in json.load(open("out/results_gmZ.json"))}
G = {o["m"]: o for o in json.load(open("out/results_gam.json"))}
FINE = [30] + list(range(35, 51))
INK, INK2, MUTED = "#1c1f23", "#4a5057", "#9aa1a8"
BLUE, RUST = "#1f5fa8", "#b4531f"
MZ2 = 91.1876 / 2.

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "axes.linewidth": 0.8, "figure.facecolor": "white"})

fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.7))
fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.16, wspace=0.30)

m = np.array(FINE, float)
fine = m >= 35

def frame(ax, title, ylab):
    ax.axvline(MZ2, color=MUTED, lw=1, ls=":")
    ax.annotate(r"$2m_\chi = m_Z$", (MZ2, ax.get_ylim()[1]), fontsize=7.5,
                color=INK2, ha="right", va="top", rotation=90,
                textcoords="offset points", xytext=(-3, -4))
    ax.set_xlabel(r"$m_\chi$ [GeV]"); ax.set_ylabel(ylab)
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.grid(color="#e6e9ec", lw=0.7); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(28, 52)

ax = axes[0]
ax.errorbar(m, [Z[k]["R"] for k in FINE], yerr=[Z[k]["R_err"] for k in FINE],
            color=BLUE, lw=1.2, ls="-", marker="o", ms=4, capsize=0,
            alpha=0.45, zorder=2)
ax.errorbar(m[fine], [Z[k]["R"] for k in np.array(FINE)[fine]],
            yerr=[Z[k]["R_err"] for k in np.array(FINE)[fine]],
            color=BLUE, lw=2, marker="o", ms=5.5, capsize=0, zorder=3)
ax.set_ylim(0.08, 0.23); frame(ax, "Acceptance", r"$R(m)$")

ax = axes[1]
ax.semilogy(m, [Z[k]["sigma_fb"] for k in FINE], color=BLUE, lw=2, marker="o",
            ms=4.5, label=r"$\gamma^*/Z$")
ax.semilogy(m, [G[k]["sigma_fb"] for k in FINE], color=RUST, lw=2, ls="--",
            marker="s", ms=3.5, label=r"$\gamma^*$ only")
ax.set_ylim(5e3, 5e6)
frame(ax, r"Cross section, $Q_\chi = e$", r"$\sigma$ [fb]")
ax.legend(frameon=False, fontsize=8, loc="lower left", handlelength=2.0)

ax = axes[2]
ax.plot(m, [Z[k]["bg_q"][2] for k in FINE], color=BLUE, lw=2, marker="o", ms=4.5)
ax.fill_between(m, [Z[k]["bg_q"][1] for k in FINE],
                [Z[k]["bg_q"][3] for k in FINE], color=BLUE, alpha=0.14, lw=0)
ax.axhline(1.0, color=MUTED, lw=1, ls=":")
ax.set_ylim(0, 4.2)
frame(ax, r"Velocity of accepted $\chi$", r"$\beta\gamma$  (median, 16-84%)")

fig.savefig("out/mcp_dy_zoom.png", dpi=180, bbox_inches="tight")
print("wrote out/mcp_dy_zoom.png")
