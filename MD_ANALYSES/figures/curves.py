#!/usr/bin/env python3
"""
curves.py -- the static figures, from the tables the stages already wrote.

Reads `data/*.tsv`, never a trajectory: a figure can always be redrawn in
seconds, and redrawing one cannot change the numbers. This module must not
import MDAnalysis, and the page modules must not import this -- three ways to
look at one set of rows, one source of truth.
"""

import os
from collections import defaultdict

from core import io as _io

from . import style

#: What the RMSD/RMSF panels show. The rest stays in data/*.tsv.
PLOT_SUBSETS = ("receptor", "ligand", "epitope")


def _read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def load(data_dir, stage=None):
    rows = []
    for name in sorted(os.listdir(data_dir)):
        if not name.endswith(".tsv"):
            continue
        if stage and f"__{stage}." not in name:
            continue
        rows += list(_read_tsv(os.path.join(data_dir, name)))
    return rows


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def timeseries(rows, out_dir, observable="rmsd", title=None):
    """One panel per subset, replicas overlaid: the shape of the run."""
    plt = style.use()
    series = defaultdict(lambda: defaultdict(lambda: ([], [])))
    for r in rows:
        if r.get("scope") != "frame" or r.get("observable") != observable:
            continue
        if not _io.absent(r.get("protomer")):
            continue          # per-protomer rows are their own observable
        x, y = series[r["subset"]][r["replica"]]
        t, v = _num(r.get("time_ns")), _num(r.get("value"))
        if t is not None and v is not None:
            x.append(t)
            y.append(v)
    if not series:
        return []
    # Same choice as the page: glycans and the whole-solute trace are in the
    # tables, not in these panels.
    subsets = [s for s in PLOT_SUBSETS if s in series] or sorted(series)
    fig, axes = plt.subplots(len(subsets), 1, sharex=True, squeeze=False,
                             figsize=(style.WIDTH_2COL,
                                      1.05 * len(subsets) + 0.5))
    for ax, sub in zip(axes[:, 0], subsets):
        reps = sorted(series[sub])
        for i, rep in enumerate(reps):
            x, y = series[sub][rep]
            ax.plot(x, y, style.DASHES[i % len(style.DASHES)],
                    color=style.COLOR.get(sub, "#555"), label=rep)
        ax.set_ylabel(f"{style.label(observable)}\n({style.label(sub)}, Å)")
        ax.grid(True, axis="y")
        ax.margins(x=0.01)
        if ax is axes[0, 0]:
            ax.legend(ncol=len(reps), loc="upper left")
    axes[-1, 0].set_xlabel("time (ns)")
    if title:
        fig.suptitle(title, y=1.0, fontsize=8, ha="left", x=0.0)
    fig.align_ylabels(axes[:, 0])
    out = style.save(fig, os.path.join(out_dir, f"{observable}_timeseries"))
    plt.close(fig)
    return out


def per_residue(rows, out_dir, observable="rmsf", subsets=("receptor", "ligand")):
    """RMSF along the chain, replicas overlaid, one panel per subset."""
    plt = style.use()
    series = defaultdict(lambda: defaultdict(list))
    labels = defaultdict(list)
    for r in rows:
        if r.get("scope") != "residue" or r.get("observable") != observable:
            continue
        if not _io.absent(r.get("protomer")):
            continue
        if subsets and r["subset"] not in subsets:
            continue
        v = _num(r.get("value"))
        if v is not None:
            series[r["subset"]][r["replica"]].append(v)
            if r["replica"] == sorted(series[r["subset"]])[0]:
                labels[r["subset"]].append(r.get("segid") or "")
    if not series:
        return []
    subs = [s for s in subsets if s in series] or sorted(series)
    fig, axes = plt.subplots(len(subs), 1, squeeze=False,
                             figsize=(style.WIDTH_2COL, 1.15 * len(subs) + 0.5))
    for ax, sub in zip(axes[:, 0], subs):
        reps = sorted(series[sub])
        for i, rep in enumerate(reps):
            y = series[sub][rep]
            ax.plot(range(len(y)), y, style.DASHES[i % len(style.DASHES)],
                    color=style.COLOR.get(sub, "#555"), label=rep, linewidth=0.6)
        # Mark where one chain ends and the next begins, so a reader is not
        # looking at a concatenation of three copies as if it were one chain.
        segs = labels.get(sub) or []
        last, start = None, 0
        for j, s in enumerate(segs + [None]):
            if s != last:
                if last is not None and j - start > 3:
                    ax.axvline(j - 0.5, color="#bbb", linewidth=0.4)
                    ax.annotate(last, (start + (j - start) / 2, 1.0),
                                xycoords=("data", "axes fraction"),
                                ha="center", va="bottom", fontsize=5.5,
                                color="#777", annotation_clip=False)
                last, start = s, j
        ax.set_ylabel(f"{style.label(observable)}\n({style.label(sub)}, Å)")
        ax.set_xlabel("residue (in selection order)")
        ax.grid(True, axis="y")
        ax.margins(x=0.005)
        if ax is axes[0, 0]:
            ax.legend(ncol=len(reps), loc="upper right")
    out = style.save(fig, os.path.join(out_dir, f"{observable}_per_residue"))
    plt.close(fig)
    return out


def composition_card(rows, out_dir):
    """The composition table as a figure, for a slide that needs one."""
    plt = style.use()
    want = ["water_molecules", "ions_total", "salt_concentration",
            "box_volume", "solute_atoms", "glycan_residues"]
    pretty = {"water_molecules": "Water molecules", "ions_total": "Ions",
              "salt_concentration": "Salt (mol/L)", "box_volume": "Box (Å³)",
              "solute_atoms": "Solute atoms", "glycan_residues": "Glycan residues"}
    got, reps = defaultdict(dict), []
    for r in rows:
        if r.get("stage") != "composition" or r.get("scope") != "run":
            continue
        key = r["observable"]
        if key not in want:
            continue
        got[key][r["replica"]] = _num(r.get("value"))
        if r["replica"] not in reps:
            reps.append(r["replica"])
    if not got:
        return []
    reps = sorted(reps)
    fig, ax = plt.subplots(figsize=(style.WIDTH_1COL, 0.28 * len(want) + 0.5))
    ax.axis("off")
    cells = [[f"{got[k].get(rep):,.4g}" if got[k].get(rep) is not None else "—"
              for rep in reps] for k in want if k in got]
    tab = ax.table(cellText=cells, rowLabels=[pretty[k] for k in want if k in got],
                   colLabels=reps, loc="center", cellLoc="right")
    tab.auto_set_font_size(False)
    tab.set_fontsize(6.5)
    tab.scale(1, 1.25)
    for (_r, _c), cell in tab.get_celld().items():
        cell.set_linewidth(0.3)
        cell.set_edgecolor("#dddddd")
    out = style.save(fig, os.path.join(out_dir, "composition"))
    plt.close(fig)
    return out


def all_for(bundle_dir):
    """Every figure this bundle's data supports."""
    data = os.path.join(bundle_dir, "data")
    figs = os.path.join(bundle_dir, "figures")
    rows = load(data)
    out = []
    out += timeseries(rows, figs, "rmsd")
    out += timeseries(rows, figs, "rg")
    out += per_residue(rows, figs, "rmsf")
    out += composition_card(rows, figs)
    return out
