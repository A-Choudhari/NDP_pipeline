#!/usr/bin/env python3
"""
topology.py -- what a system is made of, read from its topology alone.

Pure standard library: a 105 MB solvated PSF is parsed in a few seconds and
nothing here needs MDAnalysis, so component detection runs before the analysis
environment is even activated, and the answer can be shown to a person for
confirmation before any trajectory is opened.

WHAT IT READS. The CHARMM PSF `!NATOM` block, which is the only section that
matters for deciding what a system contains:

    atomid  segid  resid  resname  atomname  type  charge  mass  ...

Atom ORDER is preserved exactly; every index this module reports is an index
into that order, which is also the order of coordinates in the DCD.
"""

import hashlib
import os
import re
from collections import Counter, OrderedDict

#: Residue names by class. Deliberately explicit rather than clever: a system
#: whose residues are not listed comes back as `other`, which is a question for
#: a person, not something to guess at.
PROTEIN = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # CHARMM protonation variants of histidine, and the termini/cysteine forms
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH",
    "LYN", "ACE", "CT3", "NME", "NMA",
}

NUCLEIC = {
    "ADE", "GUA", "CYT", "THY", "URA", "DA", "DC", "DG", "DT", "DU",
    "A", "C", "G", "U", "RA", "RC", "RG", "RU",
    "DA5", "DA3", "DC5", "DC3", "DG5", "DG3", "DT5", "DT3",
}

#: The glycan set is the VMD macro this workspace already relies on
#: (~/.vmdrc line 44), plus the CHARMM names seen in these trimers
#: (BGLC AMAN BMAN BGAL AFUC ANE5). Keeping the same list means a selection
#: written for VMD and one written here pick the same atoms.
GLYCAN = {
    "NAG", "BGL", "BGLN", "BGLCNA", "AGLCNA", "AFU", "FUC", "AFUC", "BFUC",
    "AGAL", "BGAL", "GAL", "AMA", "MAN", "AMAN", "BMA", "BMAN", "BCNA",
    "ANE", "ANE5", "ANE5AC", "BNE5AC", "ANE5GC", "BNE5GC", "AGAN", "AGALNA",
    "AGALNAC", "BGALNAC", "BGLC", "AGLC", "ASIA", "BSIA", "AXYL", "BXYL",
    "ARAM", "BRAM", "BGALNA",
}

WATER = {"TIP3", "TIP3P", "TIP4", "TP4E", "SPC", "SPCE", "WAT", "HOH", "SOL",
         "T3P", "TIP", "OPC"}

#: Monatomic ions. Anything polyatomic (phosphate, acetate) is a ligand.
ION = {"SOD", "CLA", "POT", "CES", "CAL", "MG", "MG2", "ZN", "ZN2", "LIT",
       "RUB", "BAR", "NA", "NA+", "CL", "CL-", "K", "K+", "CA", "CA2", "FE",
       "FE2", "FE3", "MN", "CU", "CU1", "CU2", "IOD", "BR"}

#: Lipids, from the same ~/.vmdrc macro (line 45) plus the common CHARMM36 set.
LIPID = {"DLPE", "DMPC", "DPPC", "GPC", "LPPC", "PALM", "PC", "PGCL", "POPC",
         "POPE", "POPS", "POPI", "CHL1", "PSM", "DOPC", "DOPE", "DOPS",
         "SAPI", "SOPC", "DPPE", "DLPC", "DMPE", "CHOL", "TIP3C"}

CLASSES = ("protein", "nucleic", "glycan", "water", "ion", "lipid", "other")


def classify_resname(resname):
    r = str(resname).strip().upper()
    if r in WATER:
        return "water"
    if r in ION:
        return "ion"
    if r in PROTEIN:
        return "protein"
    if r in GLYCAN:
        return "glycan"
    if r in NUCLEIC:
        return "nucleic"
    if r in LIPID:
        return "lipid"
    return "other"


class Atom:
    """One row of the !NATOM block. Slots, because there are 767,697 of them."""

    __slots__ = ("index", "segid", "resid", "resname", "name", "mass")

    def __init__(self, index, segid, resid, resname, name, mass):
        self.index = index
        self.segid = segid
        self.resid = resid
        self.resname = resname
        self.name = name
        self.mass = mass


_NATOM = re.compile(r"^\s*(\d+)\s*!NATOM")


def read_psf(path, limit=None):
    """Atoms of a CHARMM PSF, in file order.

    `limit` stops early, for a cheap look at a huge solvated topology when only
    the composition of the first segments is wanted.
    """
    atoms = []
    with open(path, errors="ignore") as fh:
        n = None
        for line in fh:
            m = _NATOM.match(line)
            if m:
                n = int(m.group(1))
                break
        if n is None:
            raise ValueError(f"no !NATOM record in {path}")
        for i in range(n if limit is None else min(n, limit)):
            line = fh.readline()
            if not line:
                raise ValueError(f"{path}: !NATOM says {n} atoms, file ends at {i}")
            f = line.split()
            if len(f) < 8:
                raise ValueError(f"{path}: short !NATOM line {i + 1}: {line!r}")
            # id segid resid resname name type charge mass
            try:
                mass = float(f[7])
            except ValueError:
                mass = 0.0
            atoms.append(Atom(i, f[1], f[2], f[3].upper(), f[4].upper(), mass))
    return atoms


def read_topology(path, limit=None):
    """Atoms of a topology. PSF today; the door is open for prmtop/PDB."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".psf":
        return read_psf(path, limit=limit)
    raise NotImplementedError(
        f"{ext or path} is not a topology this module reads yet (PSF only)")


class Segment:
    """One segid: its class, size, residues, and a fingerprint of its sequence."""

    __slots__ = ("segid", "kind", "n_atoms", "n_residues", "resnames",
                 "first_index", "last_index", "sequence_hash")

    def __init__(self, segid):
        self.segid = segid
        self.kind = None
        self.n_atoms = 0
        self.n_residues = 0
        self.resnames = Counter()
        self.first_index = None
        self.last_index = None
        self.sequence_hash = None

    def as_dict(self):
        return {"segid": self.segid, "kind": self.kind, "atoms": self.n_atoms,
                "residues": self.n_residues,
                "sequence_hash": self.sequence_hash,
                "range": [self.first_index, self.last_index]}

    def __repr__(self):                                    # pragma: no cover
        return f"<Segment {self.segid} {self.kind} {self.n_atoms} atoms>"


def segments(atoms):
    """Group atoms into segments, in first-appearance order.

    A segment's KIND is the class of the majority of its residues, because a
    protein segment may legitimately carry a few non-standard residues (a
    phosphoserine, a covalent ligand) without ceasing to be protein.

    The SEQUENCE HASH is the ordered residue-name sequence, so two segments
    that hold the same chain -- the three copies of gp120, the three copies of
    a Fab light chain -- hash identically and can be recognised as copies of one
    entity without any naming convention.
    """
    out = OrderedDict()
    seen = {}
    seq = {}
    for a in atoms:
        s = out.get(a.segid)
        if s is None:
            s = out[a.segid] = Segment(a.segid)
            seen[a.segid] = set()
            seq[a.segid] = []
            s.first_index = a.index
        s.n_atoms += 1
        s.last_index = a.index
        key = (a.resid, a.resname)
        if key not in seen[a.segid]:
            seen[a.segid].add(key)
            s.n_residues += 1
            s.resnames[a.resname] += 1
            seq[a.segid].append(a.resname)
    for segid, s in out.items():
        kinds = Counter(classify_resname(r) for r in seq[segid])
        s.kind = kinds.most_common(1)[0][0]
        s.sequence_hash = hashlib.sha1(
            "|".join(seq[segid]).encode()).hexdigest()[:12]
    return out


class Entity:
    """One distinct molecule, and every segment that is a copy of it."""

    __slots__ = ("kind", "sequence_hash", "segids", "n_residues", "n_atoms")

    def __init__(self, kind, sequence_hash):
        self.kind = kind
        self.sequence_hash = sequence_hash
        self.segids = []
        self.n_residues = 0
        self.n_atoms = 0

    @property
    def copies(self):
        return len(self.segids)

    def as_dict(self):
        return {"kind": self.kind, "copies": self.copies,
                "segids": list(self.segids), "residues": self.n_residues,
                "atoms_per_copy": self.n_atoms // max(1, self.copies)}

    def __repr__(self):                                    # pragma: no cover
        return (f"<Entity {self.kind} x{self.copies} "
                f"{self.n_residues} res {self.segids[:3]}>")


def entities(segs, kinds=("protein", "nucleic", "other")):
    """Distinct molecules among the segments, largest first.

    Water, ions, lipids and glycans are not entities in this sense: they are
    bulk, and counting 48 water segments as 48 molecules would say nothing.
    They are reported by `composition()` instead.
    """
    out = {}
    for s in segs.values():
        if s.kind not in kinds:
            continue
        e = out.get(s.sequence_hash)
        if e is None:
            e = out[s.sequence_hash] = Entity(s.kind, s.sequence_hash)
            e.n_residues = s.n_residues
        e.segids.append(s.segid)
        e.n_atoms += s.n_atoms
    return sorted(out.values(),
                  key=lambda e: (-(e.n_atoms // max(1, e.copies)), e.segids[0]))


def bonded_segments(path, atoms):
    """Which segments are covalently bonded to which, from the PSF !NBOND block.

    The evidence that settles what a receptor IS. On the CDR trimers there are
    no protein-protein bonds at all (gp120 and gp41 are separate chains, and
    CHARMM-GUI did not patch the interchain disulfides), but there ARE 75
    glycan-to-protein bonds -- and those say that gp120 and gp41 carry
    N-glycans while the Fab chains do not. A glycoprotein antigen and an
    aglycosylated antibody are told apart by their chemistry, not by a size
    heuristic that would have to guess.

    Returns {segid: {segid, ...}}, only for pairs in DIFFERENT segments, and
    only among non-water segments (a water-water bond says nothing).
    """
    seg = [a.segid for a in atoms]
    n_bonds = None
    with open(path, errors="ignore") as fh:
        for line in fh:
            if "!NBOND" in line:
                n_bonds = int(line.split()[0])
                break
        if not n_bonds:
            return {}
        need, toks = n_bonds * 2, []
        while len(toks) < need:
            line = fh.readline()
            if not line:
                break
            toks += line.split()
    out = {}
    water_or_ion = {"water", "ion"}
    for i in range(0, min(len(toks), need), 2):
        try:
            a, b = int(toks[i]) - 1, int(toks[i + 1]) - 1
        except ValueError:
            continue
        if not (0 <= a < len(seg) and 0 <= b < len(seg)) or seg[a] == seg[b]:
            continue
        if classify_resname(atoms[a].resname) in water_or_ion:
            continue
        if classify_resname(atoms[b].resname) in water_or_ion:
            continue
        out.setdefault(seg[a], set()).add(seg[b])
        out.setdefault(seg[b], set()).add(seg[a])
    return out


def glycan_hosts(links, segs):
    """{protein segid: [its glycan segids]} -- who carries which glycans."""
    out = {}
    for segid, partners in links.items():
        if segs.get(segid) is None or segs[segid].kind == "glycan":
            continue
        gl = sorted(p for p in partners
                    if segs.get(p) is not None and segs[p].kind == "glycan")
        if gl:
            out[segid] = gl
    return out


def composition(atoms, segs=None):
    """How much of each class the system holds: the numbers a record states.

    Water and ions are counted as MOLECULES, not atoms -- "710,514 atoms of
    TIP3" is true and useless; "236,838 water molecules" is what a reader of
    the record needs.
    """
    segs = segs or segments(atoms)
    per_class = Counter()
    for a in atoms:
        per_class[classify_resname(a.resname)] += 1
    waters, ions, by_ion = set(), set(), Counter()
    glycan_res, lipid_res = set(), set()
    for a in atoms:
        k = classify_resname(a.resname)
        key = (a.segid, a.resid)
        if k == "water":
            waters.add(key)
        elif k == "ion":
            ions.add(key)
            by_ion[a.resname] += 1
        elif k == "glycan":
            glycan_res.add(key)
        elif k == "lipid":
            lipid_res.add(key)
    return {
        "atoms_total": len(atoms),
        "atoms_by_class": {k: per_class.get(k, 0) for k in CLASSES
                           if per_class.get(k, 0)},
        "water_molecules": len(waters),
        "ions_total": len(ions),
        "ions_by_name": dict(sorted(by_ion.items())),
        "glycan_residues": len(glycan_res),
        "lipid_residues": len(lipid_res),
        "segments": len(segs),
        "solute_atoms": sum(n for k, n in per_class.items()
                            if k not in ("water", "ion")),
    }


def summary(path, limit=None):
    """Everything this module can say about one topology, ready to print."""
    atoms = read_topology(path, limit=limit)
    segs = segments(atoms)
    ents = entities(segs)
    return {"topology": path,
            "composition": composition(atoms, segs),
            "entities": [e.as_dict() for e in ents],
            "segment_kinds": dict(Counter(s.kind for s in segs.values()))}


if __name__ == "__main__":                                 # pragma: no cover
    import json
    import sys
    for p in sys.argv[1:]:
        print(json.dumps(summary(p), indent=2))
