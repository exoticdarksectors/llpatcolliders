#!/usr/bin/env python3
"""Fold the Pythia chi eta spectra through A(eta) and report R(m) and the
beta*gamma / time-of-flight information."""
import numpy as np, glob, os, json, sys, re

TAG = sys.argv[1] if len(sys.argv) > 1 else "gam"
PRE = {"gam": "mcp_m", "gmZ": "z_m"}[TAG]

AETA = os.environ.get("MCP_AETA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "external", "A_eta.npy"))
# R is the acceptance of the real DY spectrum relative to an ISOTROPIC single
# particle (flat in cos(theta)), both normalised over all eta:
#     REF = int A(eta) * (1/2) sech^2(eta) d(eta).
# REF is derived from the A(eta) actually in use, not hard-coded, so R stays a
# pure shape ratio and is comparable across geometry definitions.  (The
# superseded hand-made array gave REF = 0.02384, which is where that number in
# the git history comes from.)
# Discover the masses actually generated for this variant.
MASSES = sorted(int(re.search(rf'{PRE}(\d+)\.txt$', f).group(1))
                for f in glob.glob(f'out/{PRE}*.txt'))
C_M_PER_NS = 0.299792458

a = np.load(AETA)
ctr, A = a[0], a[1]
dEta = ctr[1] - ctr[0]
edges = np.append(ctr - dEta / 2., ctr[-1] + dEta / 2.)
ETAMAX = edges[-1]
REF = float((A * 0.5 / np.cosh(ctr) ** 2 * dEta).sum())


def read(fn):
    meta, rows = {}, []
    full = None
    with open(fn) as f:
        for line in f:
            if line.startswith("#"):
                p = line[1:].split()
                if len(p) >= 2:
                    meta[p[0]] = p[1]
            elif line.startswith("FULLHIST"):
                full = np.array(line.split()[1:], dtype=float)
            elif line[0].isalpha():
                continue
            else:
                rows.append(line.split())
    d = np.array(rows, dtype=float)
    return meta, d, full


def wquant(x, w, q):
    i = np.argsort(x)
    x, w = x[i], w[i]
    c = np.cumsum(w) - 0.5 * w
    return np.interp(np.array(q) * w.sum(), c, x)


out = []
for m in MASSES:
    meta, d, full = read(f"out/{PRE}{m}.txt")
    eta, bg, pt = d[:, 0], d[:, 1], d[:, 2]
    nAll = float(meta["nChiAll"])
    sig_mb = float(meta["sigma_mb"])
    sig_err = float(meta["sigmaErr_mb"])

    # Per-particle acceptance: A(eta) from the supplied binning, 0 outside.
    inwin = np.abs(eta) < ETAMAX
    idx = np.digitize(eta[inwin], edges) - 1
    Aper = A[idx]
    nWin = inwin.sum()

    # Convention 1: f normalised inside |eta| < 1.2 (shape-only ratio).
    I_win = Aper.mean()
    I_win_err = Aper.std(ddof=1) / np.sqrt(nWin)

    # Convention 2: f normalised over all eta (absolute per-particle accept.).
    I_all = Aper.sum() / nAll
    I_all_err = np.sqrt((Aper ** 2).sum()) / nAll

    frac_win = nWin / nAll

    # beta*gamma, weighted by the acceptance (i.e. for the particles that
    # actually land on the detector).
    w = Aper
    bgw = bg[inwin]
    beta = bgw / np.sqrt(1. + bgw ** 2)
    qs = [0.02, 0.16, 0.50, 0.84, 0.98]
    bgq = wquant(bgw, w, qs)
    betaq = wquant(beta, w, qs)
    # Time-of-flight delay relative to a massless particle, per 10 m of path.
    delay = (1. / beta - 1.) * 10. / C_M_PER_NS
    dq = wquant(delay, w, [0.50, 0.84, 0.98])
    fslow = lambda t: (w[bgw < t].sum() / w.sum())

    out.append(dict(
        m=m, sigma_fb=sig_mb * 1e12, sigma_fb_err=sig_err * 1e12,
        I_win=I_win, I_win_err=I_win_err,
        R_shape=I_win / REF, R_shape_err=I_win_err / REF,
        I_all=I_all, I_all_err=I_all_err,
        R=I_all / REF, R_err=I_all_err / REF,
        frac_win=frac_win, nWin=int(nWin),
        bg_med=bgq[2], bg_q=list(bgq), beta_q=list(betaq),
        delay_med=dq[0], delay_84=dq[1], delay_98=dq[2],
        f_bg_lt1=fslow(1.), f_bg_lt2=fslow(2.), f_bg_lt3=fslow(3.),
        pt_med=float(np.median(pt[inwin])),
        a_iso=REF, a_eta=os.path.basename(AETA),
    ))

print("=" * 108)
print(f"A(eta) = {os.path.basename(AETA)}   isotropic closure REF = {REF:.5f}")
print("q qbar -> gamma* -> chi chibar, 14 TeV, pure photon exchange, "
      "Pythia 8.315 (NNPDF2.3 LO), Q = 1 (scale sigma by eps^2)")
print("=" * 108)
hdr = (f"{'m [GeV]':>8} {'sigma(eps=1) [fb]':>18} {'int A f deta':>13}"
       f" {'R = /0.02384':>18} {'frac |eta|<1.2':>15} {'<A>|window':>11}")
print(hdr); print("-" * 108)
for o in out:
    print(f"{o['m']:>8} {o['sigma_fb']:>18.4g} "
          f"{o['I_all']:>13.5f} {o['R']:>11.4f} +- {o['R_err']:<5.4f}"
          f" {o['frac_win']:>15.4f} {o['I_win']:>11.5f}")

print()
print("=" * 108)
print("Acceptance-weighted velocity of the chi that hit the detector "
      "(weights = A(eta))")
print("=" * 108)
print(f"{'m [GeV]':>8} {'bg 2%':>8} {'bg 16%':>8} {'bg med':>8} {'bg 84%':>8}"
      f" {'beta med':>9} {'beta 2%':>9} {'f(bg<1)':>9} {'f(bg<2)':>9}"
      f" {'f(bg<3)':>9} {'dt/10m med':>11} {'dt/10m 98%':>11}")
print("-" * 108)
for o in out:
    print(f"{o['m']:>8} {o['bg_q'][0]:>8.2f} {o['bg_q'][1]:>8.2f} "
          f"{o['bg_med']:>8.2f} {o['bg_q'][3]:>8.2f} {o['beta_q'][2]:>9.4f} "
          f"{o['beta_q'][0]:>9.4f} {o['f_bg_lt1']:>9.4f} {o['f_bg_lt2']:>9.4f} "
          f"{o['f_bg_lt3']:>9.4f} {o['delay_med']:>9.3f}ns {o['delay_98']:>9.3f}ns")

json.dump(out, open(f"out/results_{TAG}.json", "w"), indent=1)

import csv
with open(f"out/results_{TAG}.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["m_GeV", "sigma_eps1_fb", "sigma_err_fb", "int_A_f_deta",
                "R", "R_err", "frac_eta_lt_1p2", "bg_p02", "bg_p16", "bg_med",
                "bg_p84", "bg_p98", "beta_med", "f_bg_lt1", "f_bg_lt2",
                "f_bg_lt3", "dt_per10m_med_ns", "dt_per10m_p98_ns"])
    for r in out:
        w.writerow([r["m"], f"{r['sigma_fb']:.5g}", f"{r['sigma_fb_err']:.3g}",
                    f"{r['I_all']:.5f}", f"{r['R']:.4f}", f"{r['R_err']:.4f}",
                    f"{r['frac_win']:.4f}"] + [f"{v:.4g}" for v in r["bg_q"]] +
                   [f"{r['beta_q'][2]:.4f}", f"{r['f_bg_lt1']:.4f}",
                    f"{r['f_bg_lt2']:.4f}", f"{r['f_bg_lt3']:.4f}",
                    f"{r['delay_med']:.3f}", f"{r['delay_98']:.3f}"])
print(f"\nwrote out/results_{TAG}.json, out/results_{TAG}.csv")
