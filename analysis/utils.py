import os
import glob

import numpy as np


def infer_sample_mass(mass_values):
    """Infer the LLP benchmark mass from the sample payload."""
    masses = np.asarray(mass_values, dtype=float)
    masses = masses[np.isfinite(masses)]
    if masses.size == 0:
        return None
    return float(np.median(masses))


def candidate_mass_tags(mass):
    """Return filename tags commonly used for mass-stamped external curves."""
    if mass is None or not np.isfinite(mass):
        return []

    mass = float(mass)
    tags = []

    if mass < 1.0:
        hundredths = int(round(mass * 10))
        tags.append(f"m0{hundredths}")

    rounded_int = int(round(mass))
    if np.isclose(mass, rounded_int):
        tags.extend([f"m{rounded_int}", f"m{rounded_int}p0", f"m{rounded_int}.0"])
    else:
        mass_text = f"{mass:g}"
        tags.extend([f"m{mass_text.replace('.', 'p')}", f"m{mass_text}"])

    deduped = []
    for tag in tags:
        if tag not in deduped:
            deduped.append(tag)
    return deduped


def overlay_mass_matched_external_curves(ax, sample_mass, external_dir,
                                         llp_pdg_id=None):
    """
    Overlay only external curves that are explicitly mass-tagged.

    Parameters
    ----------
    ax : matplotlib Axes
    sample_mass : float or None
        LLP mass in GeV inferred from the sample.
    external_dir : str
        Path to the directory containing external CSV curves.
    llp_pdg_id : int or None
        PDG ID of the LLP.  9000001 selects BKS curves for light ALP;
        anything else selects h->SS dark-Higgs curves.
    """
    LIGHT_ALP_PDG = 9000001

    if sample_mass is None:
        print("  Note: could not infer LLP mass; skipping external curves.")
        return 0

    tags = candidate_mass_tags(sample_mass)
    loaded = 0

    if llp_pdg_id == LIGHT_ALP_PDG:
        # Light ALP: use BKS curves only
        for tag in tags:
            path = os.path.join(external_dir, "BKS", f"CODEX_BKS_{tag}.csv")
            if os.path.isfile(path):
                data = np.loadtxt(path, delimiter=",")
                ax.loglog(data[:, 0], data[:, 1],
                          color="cyan", linewidth=2, linestyle="-",
                          label=f"CODEX-b B→KS (m={sample_mass:.1f} GeV)")
                loaded += 1
                break
        ax.set_ylabel(r'BR$(B \to K^{(*)} a)_{\min}$')
    else:
        # Heavy ALP / dark photon: h->SS dark-Higgs curves (mass-matched)
        stems = [
            ("MATHUSLA", "MATHUSLA", "green", "-"),
            ("CODEX-b", "CODEX", "cyan", "-"),
            ("ANUBIS", "ANUBIS", "purple", "-"),
            ("ANUBIS Opt", "ANUBISOpt", "purple", "--"),
            ("ANUBIS Cons", "ANUBISUpdateCons", "magenta", "--"),
        ]

        # Draw excluded region (grey fill) if available
        for tag in tags:
            excl_path = os.path.join(external_dir, f"excluded_{tag}.csv")
            if os.path.isfile(excl_path):
                excl = np.loadtxt(excl_path, delimiter=",")
                ax.fill_between(excl[:, 0], excl[:, 1], 1.0,
                                color='gray', alpha=0.35, zorder=0)
                ax.plot(excl[:, 0], excl[:, 1], color='gray',
                        linewidth=1.2, zorder=0)
                ax.text(np.sqrt(excl[:, 0].min() * excl[:, 0].max()), 0.4,
                        'excluded', color='gray', fontsize=14, alpha=0.9,
                        ha='center')
                break

        for label, stem, color, ls in stems:
            for tag in tags:
                path = os.path.join(external_dir, f"{stem}_{tag}.csv")
                if not os.path.isfile(path):
                    continue
                data = np.loadtxt(path, delimiter=",")
                ax.loglog(data[:, 0], data[:, 1],
                          color=color, linewidth=2, linestyle=ls, label=label)
                loaded += 1
                break

        ax.set_ylabel(r'BR$(h \to aa)_{\min}$')

    if loaded == 0:
        available = sorted(os.path.basename(p) for p in glob.glob(
            os.path.join(external_dir, "*.csv")))
        print(f"  Note: no mass-matched external curves found for "
              f"m={sample_mass:.3g} GeV.")
        if available:
            print(f"  Available external files: {', '.join(available)}")

    return loaded
