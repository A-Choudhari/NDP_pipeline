#!/usr/bin/env python3
"""
protomers.py -- the repeating units of a symmetric assembly, named A, B, C.

WHY THIS EXISTS. A C3 trimer measured as one body reports the motion of the
whole assembly: three protomers breathing apart show up as an RMSD that keeps
climbing, and nothing in that number says whether any single protomer is
actually deforming. Fitting each protomer on itself separates the two -- and
having three of them in every box turns one trajectory into three samples of
the same observable, which is what makes a mean with a band honest rather than
decorative.

HOW A PROTOMER IS FOUND. From `core.topology.entities`, which groups segments
into molecules by their residue SEQUENCE: the three copies of gp120 hash
identically and land in one entity, and so do the three gp41, the three Fab
heavy chains and the three light chains. Protomer i is then the i-th copy of
every entity, in the order the segments appear in the topology -- which is the
order the system was built in.

NEVER FROM THE SPELLING. `segid.startswith("A")` works on 1F4 and gives a wrong
answer without complaining on JL427, whose six Fab chains are called JLA, JLB,
JLC, JLK, JLX and JLY. Naming is not evidence; the topology is. This is the
class of bug this codebase keeps finding, so it is excluded by construction:
nothing here ever looks at a character of a segid.

WHICH RULE WON, MEASURED ON THE REAL TOPOLOGIES (2026-09-24, frame 0 of R1).

  Receptor -> COPY INDEX. On 1F4 the copy index gives A120+A41, B120+B41,
  C120+C41, which is what the lab's own `4-OPEN/CDR/SCRIPTS/DGscript.tcl` means
  by "chain A B C -- gp120/gp41 protomers, 632 res each" (479 + 153 = 632).
  Proximity must NOT be used here and would be wrong if it were: in this Env
  trimer each gp41 packs mostly against the NEIGHBOURING gp120. A41 makes 2,475
  heavy-atom contacts within 6 A with C120 and 105 with A120, so a "nearest
  chain" rule would build a protomer of A120+C41 and quietly measure a unit
  that does not exist.

  Ligand -> COPY INDEX, CHECKED AGAINST CONTACT. On 1F4 the copy index already
  agrees with DGscript (IFH+IFL -> A, IFM+IFN -> B, IFQ+IFS -> C) and with the
  contacts, and on YZ23 likewise. On JL427 it does not: the heavy-chain entity
  is [JLK, JLX, JLY], and JLK makes 153 contacts with the protomer built THIRD
  and 0 with the first, so the copy index would have bolted it onto the wrong
  protomer. So where frame-0 coordinates are available the ligand copy is
  assigned to the receptor protomer it is in contact with, and a disagreement
  with the copy index is RECORDED in `notes` rather than hidden.

  (Note for whoever reads DGscript next: its JL427 line, "JLA+JLB -> A
  JLX+JLY -> B", cannot be right -- JLA and JLB are copies of one another, as
  are JLX and JLY, so neither pair can be a heavy+light Fab. Its 1F4 and YZ23
  lines are confirmed by both rules here.)

WHAT IS LEFT OUT, AND SAID SO. Only entities present in the modal copy number
take part. A single chain in a hetero-complex -- one CD4, one nanobody on a
trimer -- is not a third of anything, and forcing it into a protomer would
produce three "protomers" of which two are missing a piece. It is returned in
`excluded` so the stage can write a row about it.

A MONOMER IS A TRIMER WITH n=1. When everything is present once, one protomer
is returned holding the whole receptor and the whole ligand, so no caller ever
needs a special case for the non-symmetric systems.
"""

import itertools
from collections import Counter

#: Above this many copies the placement stops being brute-forced (n! grows),
#: and the copy index stands. Nothing in this workspace is above 3, and a
#: guessed answer for a 12-mer would be worse than an honest convention.
MAX_ASSIGNMENT = 8


def _label(i):
    """A, B, C ... then P27 onwards.

    Twenty-six letters is already far more protomers than anything in this
    workspace has; the numeric tail exists so the function cannot return a
    duplicate label and silently merge two units.
    """
    return chr(ord("A") + i) if i < 26 else f"P{i + 1}"


def _selection(segids):
    """Segids as the VMD/CHARMM string the rest of the toolkit speaks."""
    return ("segid " + " ".join(segids)) if segids else None


class Split:
    """The protomers of one assembly, plus everything that did not fit in one.

    Iterates as the list of protomer records, so a caller that only wants the
    units can write `for p in split(...)` and ignore the rest.
    """

    __slots__ = ("protomers", "excluded", "n", "rule", "notes", "entity_of")

    def __init__(self, protomers, excluded, n, rule, notes, entity_of):
        self.protomers = protomers
        self.excluded = excluded
        #: How many copies of the repeating unit there are (3 here).
        self.n = n
        #: Which rule placed the ligand copies: "copy_index" or "contact".
        self.rule = rule
        #: Anything a person should know about this split, in plain words.
        self.notes = notes
        #: segid -> the label of the entity it is a copy of. The label is the
        #: segid of the copy that belongs to protomer A, because there is no
        #: name in a topology for "gp120": it is a real chain, it is stable
        #: across replicas, and it is what makes (entity, resid) a key that
        #: identifies THE SAME residue of three copies.
        self.entity_of = entity_of

    def __iter__(self):
        return iter(self.protomers)

    def __len__(self):
        return len(self.protomers)

    def __getitem__(self, i):
        return self.protomers[i]

    def as_dict(self):
        return {"n": self.n, "rule": self.rule, "protomers": self.protomers,
                "excluded": self.excluded, "notes": list(self.notes)}

    def __repr__(self):                                    # pragma: no cover
        return (f"<Split n={self.n} rule={self.rule} "
                f"{[p['id'] for p in self.protomers]} "
                f"excluded={len(self.excluded)}>")


def _entities_for(entities, segids):
    """The entities whose segments all belong to one role, in entity order.

    An entity that straddles the two roles is a selection mistake -- the same
    molecule cannot be both the receptor and the ligand -- and is returned
    separately so it can be reported instead of being half-used.
    """
    mine, straddling = [], []
    wanted = set(segids)
    if not wanted:
        return mine, straddling
    for e in entities:
        hit = [s for s in e.segids if s in wanted]
        if not hit:
            continue
        if len(hit) == len(e.segids):
            mine.append(e)
        else:
            straddling.append((e, hit))
    return mine, straddling


def _modal_copies(entities, primary):
    """How many copies the repeating unit has: the modal copy count.

    COUNTED BY ATOMS, NOT BY ENTITIES, because an entity is not a vote. A
    trimer carrying two different stray chains is (3, 1, 1) by entity and would
    come out a MONOMER on a show of hands -- one repeating unit holding the
    whole assembly, every per-protomer number silently identical to the
    whole-assembly one. Weighted by atoms, the 3-copy subunit that is most of
    the mass wins, the strays are excluded, and the split says so.

    A residual tie goes to the PRIMARY entity -- the largest receptor entity,
    the thing the assembly is made of -- rather than to whichever key the
    counter happened to see first.
    """
    if not entities:
        return 1
    mass = Counter()
    for e in entities:
        mass[e.copies] += max(1, e.n_atoms)
    best = max(mass.values())
    tied = sorted(c for c, m in mass.items() if m == best)
    if len(tied) == 1:
        return tied[0]
    return primary.copies if primary is not None and primary.copies in tied \
        else max(tied)


def contact_counts(X, receptor_groups, ligand_groups, cutoff=6.0):
    """counts[j][i] = atom pairs within `cutoff` between ligand j and
    receptor protomer i.

    The coordinates must ALREADY be in one periodic image -- in the stage they
    are, because core.pbc has imaged every measured group against the receptor
    before anything is measured. Plain Euclidean distance on un-imaged
    coordinates would report zero contacts for a Fab that `wrapAll on` wrote a
    box length away, and the pairing would then be decided by which copy
    happened to be nearest to nothing.
    """
    from scipy.spatial import cKDTree
    rt = [cKDTree(X[idx]) for idx in receptor_groups]
    lt = [cKDTree(X[idx]) for idx in ligand_groups]
    return [[int(r.count_neighbors(lig, cutoff)) for r in rt] for lig in lt]


def _assign(rows, n):
    """Place n copies on n protomers so that total contact is greatest.

    THE WHOLE ASSIGNMENT AT ONCE, not one argmax per copy. Taking each copy's
    best protomer independently breaks on a chain that barely touches the
    receptor: on 1F4 replica 2-R2 frame 0 the light chains IFL and IFN both
    peak on protomer B (IFL makes only a handful of contacts anywhere), so the
    greedy answer is not even a permutation. Maximising the total instead gives
    IFL its own protomer by elimination, which is the only consistent reading
    of the same numbers.

    n! is brute-forced because n is 3. Returns (permutation, total, how many
    permutations tie for the best) -- a tie means the contacts do not decide,
    and the caller must not pretend they did.
    """
    best, total, ties = None, -1, 0
    for perm in itertools.permutations(range(n)):
        t = sum(rows[i][perm[i]] for i in range(n))
        if t > total:
            best, total, ties = perm, t, 1
        elif t == total:
            ties += 1
    return best, total, ties


def split(entities, receptor_segids, ligand_segids=(), contacts=None):
    """Split an assembly into its protomers.

    `entities`        the list from `core.topology.entities`.
    `receptor_segids` / `ligand_segids`
                      the segids the spec's selections actually resolved to.
                      Segids, not selection strings: this module must never
                      have to parse a selection, and a caller that has a
                      Universe already knows the answer exactly.
    `contacts`        optional. Either {ligand segid: [pairs to protomer 0,
                      1, ...]}, or a callable given the receptor protomers'
                      segid lists that returns such a mapping (see
                      `contact_counts`). Supplying it turns the ligand pairing
                      from a convention into a measurement.

    Returns a `Split`, which iterates as a list of
    {"id": "A", "receptor": "segid A120 A41", "ligand": "segid IFH IFL"}.
    """
    notes = []
    rec_ents, rec_straddle = _entities_for(entities, receptor_segids)
    lig_ents, lig_straddle = _entities_for(entities, ligand_segids or ())
    excluded = []
    for e, hit in rec_straddle + lig_straddle:
        excluded.append({"role": "mixed", "segids": list(e.segids),
                         "copies": e.copies, "residues": e.n_residues,
                         "reason": f"only {len(hit)} of this entity's "
                                   f"{e.copies} copies are inside one "
                                   f"selection; receptor and ligand must not "
                                   f"split a molecule"})
        notes.append(f"entity {e.segids} straddles the receptor/ligand "
                     f"selections and takes part in no protomer")

    primary = rec_ents[0] if rec_ents else (lig_ents[0] if lig_ents else None)
    n = _modal_copies(rec_ents + lig_ents, primary)

    taking, dropped = [], []
    for role, ents in (("receptor", rec_ents), ("ligand", lig_ents)):
        for e in ents:
            (taking if e.copies == n else dropped).append((role, e))
    for role, e in dropped:
        excluded.append({"role": role, "segids": list(e.segids),
                         "copies": e.copies, "residues": e.n_residues,
                         "reason": f"{e.copies} copy/copies where the "
                                   f"assembly "
                                   f"has {n}: not a repeating unit of it"})
        notes.append(f"{role} entity {e.segids} has {e.copies} copies against "
                     f"the assembly's {n}; measured only as part of the whole")

    entity_of = {}
    for _role, e in taking:
        label = e.segids[0]
        for s in e.segids:
            entity_of[s] = label

    # A MONOMER IS A TRIMER WITH n=1. Give back one unit holding everything, so
    # the stage's loop is the same code for a trimer and for a single chain.
    if n <= 1:
        rec = [s for _r, e in taking if _r == "receptor" for s in e.segids]
        lig = [s for _r, e in taking if _r == "ligand" for s in e.segids]
        rec = rec or sorted(receptor_segids)
        lig = lig or sorted(ligand_segids or ())
        p = {"id": _label(0), "index": 0,
             "receptor_segids": rec, "ligand_segids": lig,
             "receptor": _selection(rec), "ligand": _selection(lig)}
        notes.append("one copy of everything: the single protomer IS the "
                     "whole assembly")
        return Split([p], excluded, 1, "monomer", notes, entity_of)

    rec_taking = [e for role, e in taking if role == "receptor"]
    lig_taking = [e for role, e in taking if role == "ligand"]

    # The receptor defines the units, by build order. See the module docstring
    # for why contacts are not consulted here: gp41 touches the NEIGHBOURING
    # gp120 far more than its own, so a nearest-chain rule builds A120+C41.
    receptor_groups = [[e.segids[i] for e in rec_taking] for i in range(n)]

    rule = "copy_index"
    # Where each ligand copy goes. Default: same copy index. Then, if contacts
    # were offered, the measurement gets to overrule the convention.
    place = {e.segids[i]: i for e in lig_taking for i in range(n)}
    if lig_taking and contacts is not None:
        counts = contacts(receptor_groups) if callable(contacts) else contacts
        by_contact, ok = {}, True
        if n > MAX_ASSIGNMENT:                             # pragma: no cover
            ok = False
            notes.append(f"{n} copies is more than this module will "
                         f"brute-force "
                         f"({MAX_ASSIGNMENT}); the copy index stands")
        for e in (lig_taking if ok else ()):
            rows = []
            for sg in e.segids:
                row = list(counts.get(sg) or [])
                if len(row) != n:
                    ok = False
                    notes.append(f"no contact row of length {n} for {sg}: the "
                                 f"copy index is used for the whole ligand")
                    break
                rows.append(row)
            if not ok:
                break
            perm, total, ties = _assign(rows, n)
            if not total:
                ok = False
                notes.append(f"none of {', '.join(e.segids)} touches any "
                             f"receptor protomer within the cutoff, so the "
                             f"contacts cannot place them; the copy index "
                             f"stands")
                break
            if ties > 1:
                # Equally good placements are not a placement. Say so and keep
                # the convention rather than inventing a unit with two light
                # chains and one with none.
                ok = False
                notes.append(f"{ties} equally good placements of "
                             f"{', '.join(e.segids)} on the protomers; the "
                             f"contacts do not decide, so the copy index "
                             f"stands and this split needs a human look")
                break
            for i, sg in enumerate(e.segids):
                if not any(rows[i]):
                    # Touching nothing is not the same as touching everything.
                    # On 1F4 frame 0 the heavy chain IFQ makes no heavy-atom
                    # contact within 6 A with ANY gp120 -- it sits on the
                    # glycans of its protomer -- while its light chain IFS
                    # makes 352 with C120. The assignment places it on what is
                    # left, which is the only consistent reading, and says so.
                    notes.append(f"{sg} touches no receptor protomer within "
                                 f"the cutoff; the assignment leaves it "
                                 f"protomer {_label(perm[i])}")
                by_contact[sg] = perm[i]
        if ok:
            moved = sorted(sg for sg, j in by_contact.items()
                           if place.get(sg) != j)
            if moved:
                notes.append(
                    "contact pairing disagrees with the copy index for "
                    + ", ".join(f"{sg}: copy {place[sg]} -> protomer "
                                f"{_label(by_contact[sg])}" for sg in moved)
                    + " -- the contacts win, because a ligand belongs to the "
                      "protomer it is bound to")
            place = by_contact
            rule = "contact"

    protomers = []
    for i in range(n):
        lig = [s for e in lig_taking for s in e.segids if place.get(s) == i]
        # Keep the ligand's chains in entity order (heavy then light), not in
        # the order they happened to be matched, so the selection string for
        # protomer B reads the same way as the one for protomer A.
        order = {s: k for k, e in enumerate(lig_taking) for s in e.segids}
        lig.sort(key=lambda s: order[s])
        protomers.append({
            "id": _label(i), "index": i,
            "receptor_segids": list(receptor_groups[i]),
            "ligand_segids": lig,
            "receptor": _selection(receptor_groups[i]),
            "ligand": _selection(lig),
        })
    return Split(protomers, excluded, n, rule, notes, entity_of)


def from_topology(path, receptor_segids, ligand_segids=(), contacts=None):
    """`split` straight from a topology file, for a caller with no Universe."""
    from . import topology as T
    ents = T.entities(T.segments(T.read_topology(path)))
    return split(ents, receptor_segids, ligand_segids, contacts=contacts)
