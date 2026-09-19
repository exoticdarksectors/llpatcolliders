"""
Standalone exclusion-limit panel (the bottom-right panel of the full exclusion
figure) for given masses, with the current default selection. One figure per mass.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial

LIFETIMES = np.logspace(-10.5, -3.5, 20)

# Signal-yield thresholds to draw. The exclusion normalisation is
# N_sig / (P * L * sigma), so the curve for a given N_sig is just the
# baseline (sig.N_SIG_EXCL-event) curve scaled by N_sig / N_SIG_EXCL.
N_SIG_LEVELS = (3, 10)
N_SIG_STYLES = ('-', '--')
GRENDEL_COLOR = 'blue'


def panel(csv, mass, out, mini=False, n_sig_levels=N_SIG_LEVELS):
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    scan = sig.analyze_decay_vs_lifetime(csv, geo, LIFETIMES)
    mc = sig.sample_separations(geo, 1e-6, n_samples_per_particle=200)
    mc_scan = sig.mc_exclusion_vs_lifetime(mc, LIFETIMES, scan['total_events'])

    fig, ax = plt.subplots(figsize=(7, 6))
    ctm = LIFETIMES * sig.SPEED_OF_LIGHT
    scale = 100 if mini else 1
    base = np.asarray(mc_scan['exclusion'], dtype=float) * scale
    acc = np.asarray(scan['exclusion'], dtype=float) * scale
    # N_sig thresholds share one colour, differentiated by linestyle.
    for n_sig, ls in zip(n_sig_levels, N_SIG_STYLES):
        ax.loglog(ctm, base * (n_sig / sig.N_SIG_EXCL), color=GRENDEL_COLOR, linewidth=2,
                  linestyle=ls,
                  label=fr'GRENDEL ($N_\mathrm{{sig}} > {n_sig}$)')
    # Acceptance-only (N_sig > 3), same colour, lightened.
    ax.loglog(ctm, acc, color=GRENDEL_COLOR, linewidth=2, linestyle=':',
              alpha=0.5, label='GRENDEL (acceptance only)')

    def overlay(path, color, label, ls='-',scale=1):
        if os.path.exists(path):
            d = np.loadtxt(path, delimiter=',')
            ax.loglog(d[:, 0], d[:, 1]*scale, color=color, linewidth=2, ls=ls, label=label)
    if not mini:
        if mass == 15:
            overlay('external/MATHUSLA.csv', 'green', r'MATHUSLA 3000 fb$^{-1}$','--')
            overlay('external/CODEX.csv', 'cyan', r'CODEX-b 300 fb$^{-1}$','--')
            # overlay('external/CMS.csv', 'purple', 'CMS')
            # overlay('external/ANUBISOpt_update.csv', 'magenta', r'ANUBIS (Zero background) 3000 fb$^{-1}$', '--')
            overlay('external/ANUBIS1DV.csv', 'magenta', r'ANUBIS 1DV 3000 fb$^{-1}$', '--')
            overlay('external/ANUBIS2DV.csv', 'magenta', r'ANUBIS 2DV 3000 fb$^{-1}$', ':')
            ax.set_ylim([1e-5, 0.1])
        elif mass == 0.5:
            overlay('external/CODEX0p5_2.csv', 'cyan', r'CODEX-b 300 fb$^{-1}$',"--")
            ax.set_ylim([3e-5, 0.1])
    if mini:
        if mass == 15:
            overlay('external/CMS_current.csv', 'cyan', 'CMS current','-')
            overlay('external/ATLAS_current.csv', 'cyan', 'ATLAS current','-')

    ax.set_xlabel(r'$c\tau$ (m)')
    ax.set_ylabel('BR')
    ax.set_title(f'$m = {mass}$ GeV')
    ax.grid(True, which='both', ls='-', alpha=0.2)
    ax.legend(fontsize=9, loc='lower right')
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('wrote', out, ' best excl BR =', f"{np.nanmin(mc_scan['exclusion']):.2e}")


def panel_vtx_scan(csv, mass, out, vtx_mins=(0.0, 0.30, 0.60, 1.00),
                   n_sig=3, mini=False):
    """Overlay exclusion curves for several reconstructed vertex stand-off
    requirements (d_implied >= vtx_min). The MC pass is stand-off-independent,
    so it runs once and only the cheap reweighting/selection re-runs per
    threshold.
    """
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    scan = sig.analyze_decay_vs_lifetime(csv, geo, LIFETIMES)
    mc = sig.sample_separations(geo, 1e-6, n_samples_per_particle=200)

    fig, ax = plt.subplots(figsize=(7, 6))
    ctm = LIFETIMES * sig.SPEED_OF_LIGHT
    escale = 100 if mini else 1
    cmap = plt.cm.viridis(np.linspace(0, 0.85, len(vtx_mins)))
    for vmin, color in zip(vtx_mins, cmap):
        ms = sig.mc_exclusion_vs_lifetime(mc, LIFETIMES, scan['total_events'],
                                          vtx_inner_min=vmin)
        excl = np.asarray(ms['exclusion'], dtype=float) * escale * (n_sig / sig.N_SIG_EXCL)
        label = ('no stand-off' if vmin <= 0
                 else fr'$\geq {vmin*100:.0f}$ cm from inner layer')
        ax.loglog(ctm, excl, color=color, linewidth=2, label=label)
        print(f'  vtx_min={vmin*100:.0f}cm  best excl BR = {np.nanmin(excl):.2e}')

    ax.set_xlabel(r'$c\tau$ (m)')
    ax.set_ylabel('BR')
    ax.set_title(fr'$m = {mass}$ GeV, $N_\mathrm{{sig}} > {n_sig}$')
    ax.grid(True, which='both', ls='-', alpha=0.2)
    ax.legend(fontsize=9, loc='lower right', title='reco vertex stand-off')
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('wrote', out)


def plot_vtx_variable(csv, mass, out, ctau_ref=10.0,
                      vtx_marks=(0.30, 0.60, 1.00)):
    """Diagnostic of the cut variable itself: the reconstructed
    vertex->inner-layer distance d_implied = sep_inner / open_angle, for signal
    that passes every OTHER cut. Left: decay-weighted distribution (reweighted
    to ctau_ref, near the sensitivity peak) with the candidate thresholds.
    Right: signal efficiency retained vs stand-off threshold.
    """
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    mc = sig.sample_separations(geo, 1e-6, n_samples_per_particle=200)

    # Selection with the stand-off OFF (all other cuts), so we see the marginal
    # effect of the stand-off requirement alone.
    sel = sig.selection_mask(mc, vtx_inner_min=0.0)

    # Decay weight reweighted to a representative lifetime near the peak.
    lam = mc['betagamma'] * ctau_ref            # decay length (m)
    w = (mc['path_len'] / mc['n_per']) * (1.0 / lam) * np.exp(-mc['d'] / lam)

    di = np.asarray(mc['d_implied'], dtype=float)
    ws = w[sel]
    dis = di[sel]
    finite = np.isfinite(dis)

    fig, (axh, axe) = plt.subplots(1, 2, figsize=(12, 5))

    # --- Left: weighted distribution of d_implied ---
    xmax = 2.0
    bins = np.linspace(0, xmax, 60)
    axh.hist(np.clip(dis[finite], 0, xmax), bins=bins, weights=ws[finite],
             color='steelblue', alpha=0.8)
    for vm in vtx_marks:
        axh.axvline(vm, color='crimson', ls='--', lw=1.5)
        axh.text(vm, axh.get_ylim()[1]*0.92, f'{vm*100:.0f} cm',
                 rotation=90, va='top', ha='right', color='crimson', fontsize=9)
    axh.set_xlabel(r'reco vertex stand-off  $d_\mathrm{impl}=s_\mathrm{in}/\theta_\mathrm{open}$ (m)')
    axh.set_ylabel('decay-weighted signal (a.u.)')
    axh.set_title(fr'$m={mass}$ GeV, $c\tau={ctau_ref:.0f}$ m selected signal')
    axh.set_xlim(0, xmax)

    # --- Right: efficiency retained vs threshold ---
    thr = np.linspace(0, xmax, 200)
    tot = ws.sum()
    eff = np.array([ws[dis >= t].sum() / tot for t in thr])
    axe.plot(thr, eff, color='navy', lw=2)
    for vm in vtx_marks:
        e = ws[dis >= vm].sum() / tot
        axe.axvline(vm, color='crimson', ls='--', lw=1)
        axe.plot([vm], [e], 'o', color='crimson')
        axe.annotate(f'{vm*100:.0f} cm: {e:.2f}', (vm, e),
                     textcoords='offset points', xytext=(6, 6), fontsize=9,
                     color='crimson')
    axe.set_xlabel('reco vertex stand-off requirement (m)')
    axe.set_ylabel('signal efficiency retained')
    axe.set_title('Marginal impact of the stand-off cut')
    axe.set_xlim(0, xmax)
    axe.set_ylim(0, 1.02)
    axe.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print('wrote', out, '| retained:',
          ', '.join(f'{v*100:.0f}cm={ws[dis>=v].sum()/tot:.3f}'
                    for v in vtx_marks))


panel('LLP0p5GeV.csv', 0.5, 'exclusion_panel_0p5GeV_fixCodex.png')
panel('LLPSmall.csv', 15, 'exclusion_panel_15GeV.png')
# panel('LLPSmall.csv', 15, 'exclusion_panel_15GeV_mini.png', True)
# panel_vtx_scan('LLPSmall.csv', 15, 'exclusion_panel_15GeV_vtxscan.png')
# plot_vtx_variable('LLPSmall.csv', 15, 'vtx_standoff_variable_15GeV.png')
