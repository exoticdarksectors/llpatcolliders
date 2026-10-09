#!/usr/bin/env python3
"""Per-event acceptance: A(eta) treated as a per-particle probability, with
the real eta1-eta2 correlation from the pair boost."""
import numpy as np, json, sys, glob, re, os

TAG = sys.argv[1] if len(sys.argv) > 1 else "gam"
PRE = {"gam": "mcp_m", "gmZ": "z_m"}[TAG]

a = np.load(os.environ.get("MCP_AETA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "external", "A_eta.npy")))
ctr, A = a[0], a[1]
dEta = ctr[1] - ctr[0]
edges = np.append(ctr - dEta / 2., ctr[-1] + dEta / 2.)
REF = 0.02384
# Discover the masses actually generated for this variant.
MASSES = sorted(int(re.search(rf'{PRE}(\d+)\.txt$', f).group(1))
                for f in glob.glob(f'out/{PRE}*.txt'))

def Aof(eta):
    eta = np.asarray(eta)
    out = np.zeros_like(eta)
    m = np.abs(eta) < edges[-1]
    out[m] = A[np.digitize(eta[m], edges) - 1]
    return out

rows = []
for m in MASSES:
    meta, pairs = {}, []
    for line in open(f"out/{PRE}{m}.txt"):
        if line.startswith("#"):
            p = line[1:].split()
            if len(p) >= 2: meta[p[0]] = p[1]
        elif line.startswith("PAIR"):
            pairs.append(line.split()[1:])
    nEv = float(meta["nEventAccepted"])
    nChi = float(meta["nChiAll"])
    e = np.array(pairs, dtype=float)
    A1, A2 = Aof(e[:, 0]), Aof(e[:, 1])

    perPart = Aof(np.concatenate([e[:, 0], e[:, 1]])).sum() / nChi
    Nacc = (A1 + A2).sum() / nEv                 # mean accepted chi per event
    pAny = (1. - (1. - A1) * (1. - A2)).sum() / nEv
    pBoth = (A1 * A2).sum() / nEv
    # Errors (Poisson-like, from the event-level spread).
    def sem(x):
        s = np.zeros(int(nEv)); s[:len(x)] = x
        return s.std(ddof=1) / np.sqrt(nEv)
    rows.append(dict(m=m, perPart=perPart, R=perPart / REF,
        Nacc=Nacc, Nacc_err=sem(A1 + A2),
        pAny=pAny, pAny_err=sem(1. - (1. - A1) * (1. - A2)),
        pBoth=pBoth, pBoth_err=sem(A1 * A2),
        corr=pBoth / perPart ** 2))

iso_any = 2 * REF - REF ** 2
iso_both = REF ** 2
print("=" * 104)
print("Per-event acceptance (A(eta) as a per-particle probability; eta1,eta2 "
      "correlated by the pair boost)")
print("=" * 104)
print(f"{'m [GeV]':>8} {'R per particle':>15} {'<N_acc>/evt':>13} "
      f"{'P(>=1 chi)':>13} {'P(both chi)':>15} {'P(both)/R_pp^2':>16}")
print("-" * 104)
for r in rows:
    print(f"{r['m']:>8} {r['R']:>15.4f} {r['Nacc']:>13.5f} "
          f"{r['pAny']:>13.5f} {r['pBoth']:>10.3e}+-{r['pBoth_err']:<.0e} "
          f"{r['corr']:>16.2f}")
print("-" * 104)
print(f"{'isotropic':>8} {1.0:>15.4f} {2*REF:>13.5f} {iso_any:>13.5f} "
      f"{iso_both:>10.3e}{'':>7} {1.00:>16.2f}")
print()
print("Ratios to the isotropic, uncorrelated reference:")
print(f"{'m [GeV]':>8} {'<N_acc> ratio':>15} {'P(>=1) ratio':>14} "
      f"{'P(both) ratio':>15}")
for r in rows:
    print(f"{r['m']:>8} {r['Nacc']/(2*REF):>15.4f} {r['pAny']/iso_any:>14.4f} "
          f"{r['pBoth']/iso_both:>15.4f}")
json.dump(rows, open(f"out/per_event_{TAG}.json", "w"), indent=1)
