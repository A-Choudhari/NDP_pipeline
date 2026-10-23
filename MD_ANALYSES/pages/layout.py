#!/usr/bin/env python3
"""
layout.py -- what a page shows, and in which column.

THE SHAPE THE OWNER ASKED FOR (2026-09-24): "left side the protein analyses and
on the right the ligands, if present... when we do the analyses on only the Apo
protein complex in water the graphs and maps can be in a landscape format,
while with ligands the formattation need to be splitted into the 2 main
sections."

So the page has one layout with two states, not two layouts:

    with a ligand   |  PROTEIN  |  LIGAND  |      two columns, paired rows
    apo             |        PROTEIN       |     one landscape column

and the decision is made from the DATA, not from a flag: a system whose rows
carry no `ligand` subset has no ligand half, because there is nothing to put
in it. That also means a monomer, a nucleic complex or a small-molecule
system gets the right page without anyone configuring it.

COLOUR. Cyan for the protein, red/orange for the ligand, each as a gradient so
three replicas or three protomers are told apart within one family while the
two families stay unmistakable at a glance. Taken from the owner's own
preference, and close to the lab's existing `plot_style.py` chain palette.

THE MEAN AND ITS SHADE. RMSD and RMSF are drawn as the mean over three
protomers and three replicas with a band at +/- one standard deviation --
"a shade to show the standard deviation is a must". The band is the spread
across those nine units, not a confidence interval on the mean: it says how
much the protomers and replicas disagree, which is the honest thing a reader
wants from replicate simulations.
"""

#: Cyan for protein, red/orange for ligand; light to dark within each family.
PROTEIN_GRADIENT = ["#0b6e8f", "#1596b8", "#35b8d6", "#7fd8e8"]
LIGAND_GRADIENT = ["#a01526", "#d1495b", "#e8743c", "#f2a65a"]
NEUTRAL = "#8a8f98"

#: Which family a measured subset belongs to. Anything not named here is
#: neutral and lands in the protein column, because that is where a reader
#: looks for "the rest of the system".
FAMILY = {
    "receptor": "protein", "protein": "protein", "solute": "protein",
    "glycans": "protein", "epitope": "protein", "engaged": "protein",
    "ligand": "ligand", "antibody": "ligand", "paratope": "ligand",
    "interface": "ligand",
}

#: The rows each column draws, in order. A panel whose rows are absent is
#: skipped silently -- a page for a system that has not had a stage run must
#: not show an empty frame with a title.
PANELS = {
    "protein": [
        ("rmsd", "receptor", "RMSD", "Å",
         "RMSD — backbone-fitted, per protomer, mean ± SD"),
        ("rg", "receptor", "Rg", "Å",
         "Rg — radius of gyration"),
        ("rmsf", "receptor", "RMSF", "Å",
         "RMSF — per residue, mean over protomers and replicas"),
    ],
    "ligand": [
        ("rmsd", "ligand", "RMSD", "Å",
         "RMSD — ligand, fitted on its receptor protomer"),
        ("rmsf", "ligand", "RMSF", "Å",
         "RMSF — per residue of the ligand"),
    ],
}

#: Composition rows, in the order a person reads them, with their labels.
COMPOSITION = [
    ("atoms_total", "Atoms, total"),
    ("solute_atoms", "Solute atoms"),
    ("water_molecules", "Water molecules (TIP3P)"),
    ("ions_SOD", "Na⁺ ions"),
    ("ions_CLA", "Cl⁻ ions"),
    ("ions_total", "Ions, total"),
    ("salt_concentration", "Salt (mol/L)"),
    # X x Y x Z on one line. The cell vectors are called a, b, c in the .xst
    # the numbers come from, but the box here is orthorhombic -- a is along x,
    # b along y, c along z -- so the axes are what a reader means.
    ("box_edges", "Box X × Y × Z (Å)"),
    ("box_volume", "Box volume (Å³)"),
    ("water_density", "Water density (Å⁻³)"),
]


def family_of(subset):
    return FAMILY.get(str(subset), "protein")


def has_ligand(rows):
    """True when anything was measured on a ligand.

    Read from the rows rather than from the spec: the page is built from what
    was MEASURED, so a system whose ligand stage has not run yet gets the
    landscape layout and gains the second column when the numbers arrive --
    instead of showing an empty half that looks like a failure.
    """
    for r in rows:
        if family_of(r.get("subset")) == "ligand" and r.get("value") is not None:
            return True
    return False


def gradient(family, i, n):
    """The i-th of n colours within a family's gradient."""
    g = LIGAND_GRADIENT if family == "ligand" else PROTEIN_GRADIENT
    if n <= 1:
        return g[1]
    step = (len(g) - 1) / max(1, n - 1)
    return g[min(len(g) - 1, int(round(i * step)))]


def composition_single_column(rows):
    """The composition as ONE column of values, not one per replica.

    The three replicas are the same system built once and run three times, so
    three identical columns of numbers said nothing three times. Where they do
    differ -- the box breathes independently in each -- the value is given as
    the mean and the spread is shown beside it, which is the only part of that
    table where the replicas carry information.
    """
    got = {}
    for r in rows:
        if r.get("stage") != "composition" or r.get("scope") != "run":
            continue
        key = r["observable"] + (f"_{r['subset']}" if r.get("subset") else "")
        v = r.get("value")
        if v is not None:
            got.setdefault(key, []).append(float(v))
    out = []
    for key, label in COMPOSITION:
        if key == "box_edges":
            edges = [got.get(f"box_edge_{ax}") for ax in ("a", "b", "c")]
            if not all(edges):
                continue
            means = [sum(e) / len(e) for e in edges]
            out.append((label,
                        " × ".join(f"{m:.1f}" for m in means),
                        None, len(edges[0])))
            continue
        vals = got.get(key)
        if not vals:
            continue
        mean = sum(vals) / len(vals)
        # A count is an integer. It arrives as a float because the tidy table
        # has one `value` column for every observable, and printing "668.0"
        # sodium ions makes a counted thing look measured.
        if all(float(v).is_integer() for v in vals):
            mean = int(round(mean))
        spread = max(vals) - min(vals)
        # A quantity identical in every replica is stated once. One that
        # varies gets its variation shown, because that is the only thing
        # three copies of a number can tell you.
        varies = spread > abs(mean) * 1e-6
        out.append((label, mean, spread if varies else None, len(vals)))
    return out


def read_credits(figures_dir):
    """What `credits.yaml` says about each imported picture.

    Kept here, beside the other page-shape rules, because BOTH pages show
    these figures and a second copy of the parser is how the two would start
    disagreeing about a caption. Deliberately a few lines of hand parsing
    rather than a YAML dependency: the file is written by `import_slides.py`
    in exactly this shape, and the pages must build with nothing installed
    beyond what the analyses already need.
    """
    import os
    meta, path, cur = {}, os.path.join(figures_dir, "credits.yaml"), None
    if not os.path.exists(path):
        return meta
    for line in open(path):
        line = line.strip()
        if line.startswith("- file:"):
            cur = line.split(":", 1)[1].strip()
            meta[cur] = {}
        elif cur and ":" in line and not line.startswith("#"):
            k, v = line.split(":", 1)
            meta[cur][k.strip()] = v.strip().strip('"')
    return meta
