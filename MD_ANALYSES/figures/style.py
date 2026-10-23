#!/usr/bin/env python3
"""
style.py -- the house figure style, so every panel here looks like the others.

Taken from the two style modules already in the workspace rather than invented:
`EPITOPE_FOCUS/abstyle.py` (publication geometry -- 183 mm column pair, 7 pt
type, `pdf.fonttype 42` so text stays editable in Illustrator, PNG and vector
PDF written together) and `analysis_W395/scripts/plot_style.py` (one colour per
entity, one line style per replica, so a reader can tell three replicas of one
thing from three different things at a glance).

Acronyms keep their capitals: RMSD, RMSF, SASA, Rg. They are names, not words.
"""

MM = 1 / 25.4
WIDTH_1COL = 89 * MM
WIDTH_2COL = 183 * MM

#: One colour per measured subset; one dash per replica.
COLOR = {"receptor": "#1b6ca8", "ligand": "#d1495b", "glycans": "#3f8f4a",
         "solute": "#7a5195", "protein": "#8a8f98", "epitope": "#e07b39",
         "water": "#3d9dc4", "ions": "#c99700"}
DASHES = ["-", "--", ":", "-."]

LABEL = {"rmsd": "RMSD", "rmsf": "RMSF", "rg": "Rg", "sasa": "SASA",
         "receptor": "Receptor", "ligand": "Ligand", "glycans": "Glycans",
         "epitope": "Epitope", "solute": "Solute", "protein": "Protein"}

RC = {
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "lines.linewidth": 0.8,
    "grid.color": "#e6e6e6",
    "grid.linewidth": 0.4,
    # Editable text in the vector output, which is the whole reason for it.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.dpi": 300,
    "figure.dpi": 110,
}


def use():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(RC)
    return plt


def label(key):
    return LABEL.get(key, str(key)[:1].upper() + str(key)[1:])


def save(fig, stem):
    """PNG at 300 dpi and a vector PDF, always as a pair."""
    import os
    os.makedirs(os.path.dirname(stem), exist_ok=True)
    out = []
    for ext in ("png", "pdf"):
        p = f"{stem}.{ext}"
        fig.savefig(p)
        out.append(os.path.basename(p))
    return out
