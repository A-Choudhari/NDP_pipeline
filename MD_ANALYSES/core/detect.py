#!/usr/bin/env python3
"""
detect.py -- read a topology, say what the system is, propose the selections.

The owner confirms; nothing here is treated as fact. That is deliberate: a
mis-detected component silently produces the wrong analysis set, so the
proposal is written into settings.yaml where it can be read and corrected
before anything runs.

HOW IT DECIDES. Segments are classified by the residues they hold, and protein
segments that share a residue sequence are recognised as copies of one entity
(core/topology.entities). Then:

  * the largest protein entity and anything covalently part of the same
    assembly is the RECEPTOR;
  * a smaller protein entity present in the same copy number, or in pairs, is a
    LIGAND -- and if its chains look like an immunoglobulin fold pair (two
    chains of 100-250 residues, same copy count as the receptor) it is called
    an ANTIBODY, which is only a label for choosing analyses;
  * a short chain (< 50 residues) is a PEPTIDE, a non-standard entity under
    200 atoms is a SMALL-MOLECULE, nucleic entities make it a NUCLEIC complex;
  * glycans, water, ions and lipids are bulk classes, detected by residue name.

On the CDR trimers this yields: gp120 x3 (479 res) + gp41 x3 (153 res) as the
receptor, three Fab copies (heavy 119 res + light ~110 res) as the antibody,
735 glycan residues, 236,838 waters and 1,365 ions -- which is exactly what
those systems are.
"""

import os

from . import topology as T
from .spec import Spec

#: An immunoglobulin variable+constant domain pair sits in this range. Used as
#: a LABEL for bundle selection, never as a claim about the molecule.
AB_RES_RANGE = (90, 260)
PEPTIDE_MAX_RES = 50
SMALL_MOLECULE_MAX_ATOMS = 200


def describe(topology_path, limit=None):
    """Everything detected about one topology, as plain data."""
    atoms = T.read_topology(topology_path, limit=limit)
    segs = T.segments(atoms)
    ents = T.entities(segs)
    comp = T.composition(atoms, segs)

    proteins = [e for e in ents if e.kind == "protein"]
    nucleics = [e for e in ents if e.kind == "nucleic"]
    others = [e for e in ents if e.kind == "other"]

    # WHO CARRIES THE GLYCANS. The PSF's bond list says which protein segments
    # N-glycans are attached to, and that is the discriminator a size rule
    # cannot be: gp41 is 153 residues -- squarely inside the immunoglobulin
    # range -- and was proposed as a LIGAND until this evidence was used. It
    # carries two glycans per protomer; the Fab chains carry none. A
    # glycoprotein antigen and an aglycosylated antibody differ in chemistry.
    links = T.bonded_segments(topology_path, atoms)
    hosts = T.glycan_hosts(links, segs)
    glycosylated = set(hosts)

    receptor, ligands = [], []
    if proteins:
        main = proteins[0]
        receptor = [main]
        for e in proteins[1:]:
            carries = any(sg in glycosylated for sg in e.segids)
            same_assembly = e.copies == main.copies and (
                carries or e.n_residues > AB_RES_RANGE[1])
            (receptor if same_assembly else ligands).append(e)

    components = []
    if proteins:
        components.append("protein")
    if comp.get("glycan_residues"):
        components.append("glycans")
    if nucleics:
        components.append("nucleic")
    if comp.get("lipid_residues"):
        components.append("membrane")
    if comp.get("water_molecules"):
        components.append("water")
    if comp.get("ions_total"):
        components.append("ions")

    kind_of_ligand = None
    if ligands:
        ig_like = [e for e in ligands
                   if AB_RES_RANGE[0] <= e.n_residues <= AB_RES_RANGE[1]]
        if len(ig_like) >= 2:
            kind_of_ligand = "antibody"
        elif any(e.n_residues <= PEPTIDE_MAX_RES for e in ligands):
            kind_of_ligand = "peptide"
        else:
            kind_of_ligand = "protein-partner"
        components.append(kind_of_ligand if kind_of_ligand != "protein-partner"
                          else "partner")
    for e in others:
        if e.n_atoms // max(1, e.copies) <= SMALL_MOLECULE_MAX_ATOMS:
            components.append("small-molecule")
            ligands.append(e)
            kind_of_ligand = kind_of_ligand or "small-molecule"
            break

    return {
        "topology": topology_path,
        "composition": comp,
        "components": components,
        "receptor": [e.as_dict() for e in receptor],
        "ligands": [e.as_dict() for e in ligands],
        "ligand_kind": kind_of_ligand,
        "entities": [e.as_dict() for e in ents],
        "glycan_segids": sorted(s.segid for s in segs.values()
                                if s.kind == "glycan"),
        "glycan_hosts": {k: v for k, v in sorted(hosts.items())},
        "ion_names": list(comp.get("ions_by_name", {})),
    }


def _segid_selection(segids):
    return "segid " + " ".join(segids) if segids else ""


def selections_from(found):
    """Concrete selections for what was detected, in VMD syntax."""
    sel = {}
    rec = [s for e in found["receptor"] for s in e["segids"]]
    lig = [s for e in found["ligands"] for s in e["segids"]]
    if rec:
        sel["receptor"] = _segid_selection(sorted(rec))
    if lig:
        sel["ligand"] = _segid_selection(sorted(lig))
    if found["glycan_segids"]:
        # Written as a range where the segids are a numbered family, because
        # `segid GA1 to GA23 GAE GAH` is what a person would type and read.
        sel["glycans"] = _segid_selection(_compress(found["glycan_segids"]))
    if found["composition"].get("water_molecules"):
        sel["water"] = "water"
    if found["ion_names"]:
        sel["ions"] = "resname " + " ".join(found["ion_names"])
    sel["solute"] = "not water and not ions"
    return sel


def _compress(segids):
    """`GA1..GA23, GAE, GAH` -> `GA1 to GA23 GAE GAH`, which round-trips."""
    import re
    fams, loose = {}, []
    for s in segids:
        m = re.match(r"^([A-Za-z]+)(\d+)$", s)
        if m:
            # Zero-padded numbers (G001) form their own family by width, and
            # the original text is kept: rebuilding it from int() wrote `G1`,
            # which matched no atom, and the whole system was skipped.
            digits = m.group(2)
            width = len(digits) if digits.startswith("0") else 0
            fams.setdefault((m.group(1), width), []).append((int(digits), s))
        else:
            loose.append(s)
    out = []
    for (stem, width), items in sorted(fams.items()):
        items.sort()
        nums = [n for n, _ in items]
        if len(nums) > 2 and nums == list(range(nums[0], nums[-1] + 1)):
            out.append(f"{items[0][1]} to {items[-1][1]}")
        else:
            out += [s for _, s in items]
    return out + sorted(loose)


def propose(record_id, mds, label=None, limit=None):
    """A Spec proposed from the first MD's topology.

    `mds` are core.spec.MD objects; the topology of the first is read, because
    replicas of one system share a topology by construction (that is what makes
    them replicas -- the metadata pipeline groups them by a byte-identical
    topology hash).
    """
    if not mds:
        raise ValueError(f"{record_id}: no MDs to detect from")
    top = mds[0].solvated_topology or mds[0].topology
    found = describe(top, limit=limit)
    spec = Spec(id=record_id, label=label,
                components=found["components"],
                selections=selections_from(found),
                mds=mds)
    spec.options["detected"] = {
        "topology": os.path.basename(top),
        "composition": found["composition"],
        "entities": found["entities"],
        "ligand_kind": found["ligand_kind"],
    }
    return spec, found


if __name__ == "__main__":                                 # pragma: no cover
    import json
    import sys
    for p in sys.argv[1:]:
        f = describe(p)
        print(json.dumps({k: f[k] for k in
                          ("components", "ligand_kind", "receptor", "ligands")},
                         indent=2))
        print(json.dumps(selections_from(f), indent=2))
