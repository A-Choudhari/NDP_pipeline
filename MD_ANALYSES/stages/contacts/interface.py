#!/usr/bin/env python3
"""
interface.py -- which residues of the receptor and the ligand actually touch,
how often, and for how long.

WHY THIS STAGE EXISTS. `structure` says the Fab stayed on (RMSD after fitting
the receptor) and `composition` says how many chains there are. Neither can say
WHERE the antibody sits. For an epitope paper that is the whole claim: 1F4 is
described as binding gp120 around residues 384-401, and until something counts
heavy-atom contacts per residue that description rests on the model it was built
from, not on the microseconds of dynamics that followed. A contact map is also
the only observable here that can FALSIFY a selection -- if the
highest-occupancy receptor residues come out somewhere other than the epitope,
either the system is not what the record says or `receptor`/`ligand` picked the
wrong segments, and both are things a person must be told rather than have
averaged away.

WHAT IS MEASURED, IN ONE PASS OVER THE STRIPPED TRAJECTORY.

  contact_occupancy      per (receptor residue, ligand residue) pair, the
  (scope=pair)           fraction of sampled frames in which any heavy atom of
                         one is within 4.5 A of any heavy atom of the other.
                         The map itself.
  residue_occupancy      per residue, the fraction of frames it makes ANY
  partners_per_frame     interface contact, and the mean number of partner
  (scope=residue)        residues it holds at once. Two different questions: a
                         residue can be at the interface in every frame and
                         hold one partner (a pin) or six (a hub), and a
                         histogram of the first alone cannot tell them apart.
  contact_pairs          per protomer per replica: how big the interface is,
  contacts_per_frame     how many contacts it holds at a time, and how many
  interface_residues     residues each side ever commits to it. EVER, and that
                         word is load-bearing: on 1F4 at stride 5 about 60
                         receptor residues per protomer touch the Fv at some
                         point in 500 ns, while 11 of them are in contact more
                         than half the time. The first number is the area the
                         interface explores, not the epitope; the occupancies
                         are what separate the two, and a caption that quotes
                         `interface_residues` as an epitope size is wrong.
  contact_lifetime_*     how long one residue pair stays in contact, from
  (scope=run)            geom.run_lengths -- the same kernel as ion dwell and
                         water residence, so the three cannot disagree about
                         what a run is.
  pairing_evidence       the diagnostic that says how much contact the split
                         into protomers was measured on.

ONLY PAIRS THAT WERE EVER IN CONTACT GET A ROW. The full matrix for these
systems is 1,896 x 699 = 1.3 M pairs per replica, of which some two thousand are
ever within 4.5 A. Absence from the table IS zero, stated here once rather than
written 1.3 M times.

THE STRIPPED TRAJECTORY, DELIBERATELY. `md.inputs(solvated=False)`: the solute
is the whole of what an interface is made of, and the stripped file is 1.7 GB
against 23 GB solvated -- thirteen times less I/O for exactly the same geometry.
Water-mediated contacts are a different observable and belong to the water
stage, which reads the file that has water in it.

NOTHING IS WRAPPED AND NOTHING IS IMAGED HERE. `cellOrigin` for these systems is
(169.7, 169.0, 158.7), so MDAnalysis's `wrap()` moves atoms about the wrong
reference and silently tears the solute apart. Every distance below -- the
contact search AND the neighbour counts the protomer split is measured on -- is
a minimum-image distance computed by `capped_distance(..., box=ts.dimensions)`,
which applies the convention internally and never touches a coordinate;
`core.pbc.box_lengths` is the guard that the cell really is orthorhombic, which
that arithmetic assumes. It is also why this stage does not call
`core.protomers.contact_counts`, whose docstring requires coordinates already
brought into one image: handing `split` a minimum-image count instead needs no
imaging at all. Putting the assembly back together WHOLE -- for a render, or for
anything that needs a centre -- is the structure stage's job, and has already
happened by the time these per-protomer numbers are read.

WHAT A PROTOMER IS. `core.protomers` owns that, and this stage asks it rather
than deciding for itself. That matters more here than almost anywhere: on these
Env trimers gp120 of one protomer packs against gp41 of the NEXT (A41 makes
2,475 heavy-atom contacts within 6 A with C120 and 105 with A120), so a stage
that grouped receptor chains by proximity would build a unit that does not
exist, and one that grouped them by the letter in the segid would be making
exactly the naming assumption core/selections.py exists to refuse.
`core.protomers.split` takes the copy index from the residue SEQUENCE for the
receptor and the measured contact for the ligand, and says in `notes` when the
two disagree. Those notes are copied into this stage's manifest verbatim.

EQUIVALENCE IS BY RESID WITHIN ENTITY, NEVER BY SEGID SPELLING. Two segments
holding the same residue sequence are copies of one entity, and
`split.entity_of` labels them, so residue 384 of A120 and residue 384 of B120
are the same residue of the same molecule. An aggregate pair row pools

    (receptor entity, receptor resid, ligand entity, ligand resid, offset)

where `offset` is (receptor protomer - ligand protomer) mod n: 0 is the ligand
against the protomer it is bound to, 1 or 2 the neighbouring protomer it also
brushes. The offset is what keeps the two apart -- and under the three-fold
symmetry it is also what makes a neighbour contact in protomer A the same
measurement as the one in protomer B, which a plain cross/cognate flag would
not be.

A UNIT THAT NEVER SAW A PAIR CONTRIBUTES ITS ZERO. When an aggregate is built
every (replica, protomer) unit contributes a value, and a unit in which the pair
was never in contact contributes 0.0 rather than being left out. Averaging only
over the protomers that happened to show a contact turns a contact seen in one
protomer of nine into an occupancy of 1.0, which is how a figure comes to claim
a universal interaction from a single-protomer excursion.
"""

from core import threads                                   # noqa: F401

import numpy as np

from core import frames as F
from core import geom, io, pbc, protomers
from core import selections as S
from core import spec as SP

STAGE = "interface_contacts"

#: None, because everything this stage reports by default is a TIME AVERAGE,
#: and a time average has no minimum spacing -- one frame per 20 ns estimates an
#: occupancy less precisely, not differently. The lifetime observables are the
#: exception and carry their own floor; see LIFETIME_MIN_SPACING_NS.
MIN_SPACING_NS = None

#: The floor for `contact_lifetime_*` ONLY. A side-chain contact at a protein
#: interface breaks and re-forms on the 0.1-1 ns scale, so at coarser than
#: 0.5 ns per frame the mean lifetime is a statement about the stride: every
#: episode looks like exactly one frame. Violating it does NOT kill the stage --
#: that would throw away the occupancies, which are the point -- it marks the
#: lifetime rows `resolved=False` and copies check_spacing's own message into
#: `floor_note`. Note the consequence, deliberately: at DEFAULT_STRIDE, which is
#: chosen for the map, the lifetimes are NOT resolved. Lifetimes want
#: stride 1-2.
LIFETIME_MIN_SPACING_NS = 0.5

#: Frames, matching core.spec.DEFAULT_STRIDE["contacts"]. At 0.2 ns/frame that
#: is one frame per nanosecond: an occupancy converges long before a contact map
#: changes shape, and 505 frames per replica cost seconds of arithmetic against
#: minutes of I/O.
DEFAULT_STRIDE = 5

#: Heavy-atom to heavy-atom contact distance. 4.5 A is the convention this lab's
#: interface scripts already use: far enough to hold a hydrogen bond or a salt
#: bridge through their heavy atoms, short enough not to report a second shell
#: as an interaction.
CONTACT_A = 4.5

#: The radius the EVIDENCE for the protomer split is counted at, on the first
#: sampled frame. It is NOT a contact cutoff and must not be read as one: its
#: whole job is to give `core.protomers.split` enough signal to place each
#: ligand chain on the protomer it belongs to.
#:
#: 12 A RATHER THAN THE 6 A core.protomers.contact_counts DEFAULTS TO, AND THE
#: MEASUREMENT THAT DECIDED IT. Frame 0 of the three 1F4 replicas, where the
#: split puts each Fv:
#:
#:            6 A                              12 A
#:   1-R1     IFH+IFL / IFM+IFN / IFQ+IFS      same
#:   2-R2     IFH+IFN / IFM+IFL / IFQ+IFS      IFH+IFL / IFM+IFN / IFQ+IFS
#:   3-R3     IFH+IFL / IFM+IFN / IFQ+IFS      same
#:
#: The R2 answer at 6 A tears two Fvs in half. That is not a defect in `split`,
#: which places each entity by the contact it is handed and says so in `notes`:
#: on that frame the light chain IFL has lifted off and makes ZERO heavy-atom
#: contacts with any gp120 within 6 A, so the assignment can only place it by
#: elimination. At 12 A it makes 51 with A120 and none with anything else, all
#: three replicas agree with each other, and they agree with the Fv chain
#: pairing measured directly (IFH-IFL 511 atom pairs within 5 A on that frame,
#: IFH-IFN none). More evidence, same rule.
PAIRING_A = 12.0

#: Above this share of one-frame episodes the lifetime distribution is against
#: the sampling floor: half of everything seen was over before the next sampled
#: frame, so the mean is the frame spacing and not a contact lifetime.
SINGLE_FRAME_MAX = 0.5

#: The columns this stage adds to core.io's block, declared as DEFAULTS so every
#: row carries them and the header does not depend on which rows a particular
#: system happened to produce. Two tables whose columns differ cannot be
#: concatenated without a realignment nobody would notice going wrong.
EXTRA_COLUMNS = {
    "partner_segid": None, "partner_resid": None, "partner_resname": None,
    "entity": None, "partner_entity": None, "protomer": None,
    "protomer_offset": None, "cutoff_A": CONTACT_A,
}

#: What `subset` says about a pair, given its protomer offset. `unpaired` is for
#: a chain core.protomers left out of every protomer -- one CD4 on a trimer is
#: not a third of anything -- whose contacts are real, are written, and are not
#: pooled with a repeating unit.
COGNATE, CROSS, UNPAIRED = "cognate", "cross_protomer", "unpaired"


def _universe(topology, trajectories):
    import MDAnalysis as mda
    return mda.Universe(topology, *trajectories)


def _groups(spec, u):
    """The two sides of the interface, checked before a frame is read.

    An empty selection is an error, never an empty result: "the ligand never
    came within 4.5 A of the receptor" and "the ligand selection matched
    nothing" produce the same empty table and are entirely different facts --
    the first is a finding, the second is a typo.
    """
    want = {"receptor": spec.selection("receptor"),
            "ligand": spec.selection("ligand")}
    missing = [k for k, v in want.items() if not str(v or "").strip()]
    if missing:
        raise S.EmptySelection(
            f"{STAGE} needs both a `receptor` and a `ligand` selection, and "
            f"this spec has no {' or '.join(missing)}. An apo system has no "
            f"interface, so this stage refuses rather than reporting one with "
            f"no contacts in it.")
    resolved, problems = S.resolve_all(u, want)
    if problems:
        raise S.EmptySelection("; ".join(problems))
    return resolved


def neighbour_counts(a_xyz, a_code, n_a, b_xyz, b_code, n_b, box,
                     cutoff=PAIRING_A):
    """(n_a x n_b) heavy-atom pairs within `cutoff`, by group code.

    The evidence the protomer split is measured on. Atom PAIRS rather than a
    centre-to-centre distance, because a centre is the one thing that cannot be
    computed safely here: two chains living in different periodic images have a
    midpoint in empty space, and nothing in this stage is ever wrapped. A count
    of minimum-image neighbours has no such failure mode, which is why this is
    what gets handed to `core.protomers.split` rather than its own
    `contact_counts` -- that one takes plain Euclidean distances and says it
    needs coordinates already brought into one image.

    Integer codes rather than segid strings so the whole tally is one bincount:
    at PAIRING_A the receptor-ligand search returns tens of thousands of pairs
    and a Python loop over them would be the slowest thing this stage does.
    """
    from MDAnalysis.lib.distances import capped_distance
    out = np.zeros((n_a, n_b), dtype=np.int64)
    pairs = capped_distance(np.asarray(a_xyz, dtype=np.float32),
                            np.asarray(b_xyz, dtype=np.float32),
                            cutoff, box=box, return_distances=False)
    if not len(pairs):
        return out
    a = np.asarray(a_code, dtype=np.int64)[pairs[:, 0]]
    b = np.asarray(b_code, dtype=np.int64)[pairs[:, 1]]
    flat = np.bincount(a * n_b + b, minlength=n_a * n_b)
    return flat.reshape(n_a, n_b).astype(np.int64)


def split_assembly(topology, rec_seg_atom, rec_xyz, lig_seg_atom, lig_xyz,
                   box, cutoff=PAIRING_A):
    """The protomers, from `core.protomers`, measured on this one frame.

    This stage supplies the contact evidence and nothing else. Which chains are
    copies of which molecule, which copy indexes them, and how a ligand chain
    with no contacts at all is placed are all `core.protomers`' decisions, and
    keeping them there is what stops two stages from disagreeing about what
    protomer B is.

    Returns (Split, {ligand segid: [pairs to protomer 0, 1, ...]}).
    """
    rec_segids = sorted(set(str(s) for s in rec_seg_atom))
    lig_segids = sorted(set(str(s) for s in lig_seg_atom))
    r_code = {s: i for i, s in enumerate(rec_segids)}
    l_code = {s: i for i, s in enumerate(lig_segids)}
    per_segid = neighbour_counts(
        lig_xyz, [l_code[str(s)] for s in lig_seg_atom], len(lig_segids),
        rec_xyz, [r_code[str(s)] for s in rec_seg_atom], len(rec_segids), box,
        cutoff=cutoff)
    evidence = {}

    def contacts(receptor_groups):
        # `split` hands back the receptor protomers it built, so the per-segid
        # tally above is summed into per-protomer rows without this stage ever
        # having to decide for itself which segids make a protomer.
        for s in lig_segids:
            evidence[s] = [int(sum(per_segid[l_code[s], r_code[g]]
                                   for g in grp))
                           for grp in receptor_groups]
        return evidence

    return (protomers.from_topology(topology, rec_segids, lig_segids,
                                    contacts=contacts), evidence)


def protomer_offsets(protomer_of_residue, index_of, units, n):
    """offset[ligand unit, receptor residue]: how far round the ring it is.

    0 means the receptor residue belongs to the protomer this ligand is bound
    to, 1 and 2 mean the next protomers round, and -1 means a chain
    core.protomers left out of every protomer and that therefore has no
    position in the ring at all.

    THIS IS THE AGGREGATION KEY, which is why it is a difference and not a
    boolean. Under three-fold symmetry the Fv of protomer A brushing protomer B
    and the Fv of protomer B brushing protomer C are THE SAME measurement, and
    pooling them is what makes a neighbour contact a number with a spread
    instead of three unrelated rows. A cognate/cross flag would pool them with
    the A-to-C contact as well, which is a different piece of geometry.
    """
    out = np.full((len(units), len(protomer_of_residue)), -1, dtype=np.int16)
    for ui, pid in enumerate(units):
        for j, p in enumerate(protomer_of_residue):
            out[ui, j] = -1 if p is None else (index_of[p] - index_of[pid]) % n
    return out


def residue_pairs(rec_xyz, lig_xyz, res_of_rec_atom, res_of_lig_atom,
                  n_lig_res, box, cutoff=CONTACT_A):
    """The distinct (receptor residue, ligand residue) pairs in contact.

    ONE box-aware neighbour search per frame, then a residue-level unique: a
    pair counts ONCE however many of its atoms are inside the cutoff, because
    the observable is "are these two residues touching", not "how many of their
    atoms are". Without the unique, a large residue pair with twelve close atoms
    would outweigh a hydrogen bond twelvefold and the map would be a map of
    residue size.

    `box` goes straight to `capped_distance`, which applies the minimum image
    itself and leaves every coordinate where it was -- the whole reason nothing
    here has to be wrapped.
    """
    from MDAnalysis.lib.distances import capped_distance
    empty = np.zeros(0, dtype=np.int64)
    pairs = capped_distance(np.asarray(rec_xyz, dtype=np.float32),
                            np.asarray(lig_xyz, dtype=np.float32),
                            cutoff, box=box, return_distances=False)
    if not len(pairs):
        return empty, empty
    n_lig_res = np.int64(n_lig_res)
    r = np.asarray(res_of_rec_atom, dtype=np.int64)[pairs[:, 0]]
    l = np.asarray(res_of_lig_atom, dtype=np.int64)[pairs[:, 1]]
    key = np.unique(r * n_lig_res + l)
    return key // n_lig_res, key % n_lig_res


def lifetime_stats(lengths):
    """mean, median and the one-frame share of a set of episode lengths."""
    a = np.asarray(lengths, dtype=np.float64)
    if not len(a):
        return None
    return {"n": int(len(a)), "mean": float(a.mean()),
            "median": float(np.median(a)),
            "single_frame_fraction": float((a == 1).mean())}


def lifetime_floor(stats, spacing_ns, prior_note=None):
    """Whether a contact lifetime is resolved, and the note that says why not.

    TWO INDEPENDENT GUARDS, because either alone can be fooled. The a-priori one
    is core.spec.check_spacing's verdict on the frame spacing, passed in as
    `prior_note`. The measured one is here: a distribution in which most
    episodes last exactly one sampled frame has been cut off by the sampling and
    its mean is the frame spacing, not a lifetime -- which happens at a spacing
    check_spacing would pass, whenever the contact turns over faster than
    anyone expected.
    """
    note = prior_note
    if stats["single_frame_fraction"] >= SINGLE_FRAME_MAX:
        extra = (f"{stats['single_frame_fraction']:.0%} of {stats['n']} "
                 f"episodes last a single sampled frame"
                 + (f" ({spacing_ns:g} ns)" if spacing_ns else "")
                 + "; the mean lifetime is then the frame spacing, not the "
                   "contact's")
        note = f"{note} {extra}" if note else extra
    return note is None, note


def pool(values, replicas):
    """An aggregate over the pooled units: mean, sd, n, s.e.m. over REPLICAS.

    The three protomers of one trimer share a box, a thermostat and a starting
    structure -- they are pseudoreplicates, not independent samples -- so the
    spread across all nine units is reported as `sd` (what the protomers did)
    and the standard error is taken over the three REPLICA means only (what an
    independent repeat would do). Quoting sd/sqrt(9) instead would claim three
    times the precision the sampling supports, which is the mistake core/io.py's
    SCOPES note exists to prevent.
    """
    a = np.asarray(values, dtype=np.float64)
    out = {"value": float(a.mean()), "n": int(len(a)),
           "sd": float(a.std(ddof=1)) if len(a) > 1 else None, "sem": None}
    per_replica = {}
    for r, v in zip(replicas, values):
        per_replica.setdefault(r, []).append(float(v))
    means = np.asarray([float(np.mean(v))
                        for _, v in sorted(per_replica.items())])
    if len(means) > 1:
        out["sem"] = float(means.std(ddof=1) / np.sqrt(len(means)))
    out["n_replicas"] = int(len(means))
    return out


class Pooled:
    """Every unit's value for every key, so a missing unit contributes its zero.

    Keys live in GROUPS, and each group carries its own list of units. The
    groups all end up the same size here -- core.protomers guarantees that every
    participating entity has exactly one copy per protomer -- but they are kept
    apart because they need not be: an entity excluded from the split has no
    protomer at all, and zero-filling it against nine units would divide by nine
    something that never had nine chances to happen.
    """

    def __init__(self):
        self.units = {}          # group -> [(replica, protomer), ...] in order
        self.values = {}         # group -> {key -> {(replica, protomer): v}}
        self.meta = {}           # group -> {key -> the row's display fields}

    # `protomer`, not `unit`: `unit` is a core.io COLUMN (fraction, count, ns)
    # and it arrives in **meta, so the two must not share a parameter name.
    def register(self, group, replica, protomer):
        seen = self.units.setdefault(group, [])
        if (replica, protomer) not in seen:
            seen.append((replica, protomer))

    def add(self, group, key, replica, protomer, value, **meta):
        self.register(group, replica, protomer)
        self.values.setdefault(group, {}).setdefault(key, {})
        self.values[group][key][(replica, protomer)] = float(value)
        self.meta.setdefault(group, {}).setdefault(key, meta)

    def items(self, group):
        """(key, meta, values, replicas) per key, zero-filled, in key order."""
        units = self.units.get(group, [])
        for key in sorted(self.values.get(group, {})):
            got = self.values[group][key]
            yield (key, self.meta[group][key],
                   [got.get(u, 0.0) for u in units],
                   [r for r, _ in units])


def pooled_groups(rec_seg, lig_seg, entity_of):
    """Every Pooled group one replica will write into.

    Registered up front for every protomer, so that a residue seen in one
    protomer and nowhere else is still averaged over all of them. Without this
    the zero-fill would have nothing to fill against and that residue would
    average to itself -- an occupancy of 1.0 out of a single-protomer excursion.
    """
    out = ["pair", "run"]
    for side, segs in (("receptor", rec_seg), ("ligand", lig_seg)):
        for entity in sorted({entity_of[s] for s in segs if s in entity_of}):
            for obs in ("residue_occupancy", "partners_per_frame"):
                out.append(f"residue:{side}:{entity}:{obs}")
    return out


def _sample(u, picked, stride):
    """Iterate the wanted frames the cheaper of the two available ways.

    Above about stride 3 a seek on this NFS share is cheaper than reading the
    frames in between, so the wanted frames are reached by random access; at
    stride 1 or 2 the sequential walk wins and the unwanted frames are simply
    not used. Either way the walk STOPS at the last wanted frame -- a
    `max_frames` that skipped the remainder instead would still read the whole
    1.7 GB, which is this stage's entire cost. Same rule as
    stages/solvent/ions.py: it is a fact about this filesystem, not about the
    analysis.
    """
    if stride >= 3:
        for ts in u.trajectory[list(picked)]:
            yield ts
        return
    want, last = set(picked), max(picked)
    for ts in u.trajectory:
        if ts.frame in want:
            yield ts
        if ts.frame >= last:
            break


def _residue_table(ag):
    """Atom -> local residue index, and (segid, resid, resname) per residue.

    Computed ONCE per replica. Rebuilding it from the AtomGroup inside the frame
    loop is what made the first draft of the ion stage unusable, and there are
    2,525 frames here too.
    """
    res_ids, res_of_atom = np.unique(ag.resindices, return_inverse=True)
    res_of_atom = np.ravel(res_of_atom)
    residues = ag.universe.residues[res_ids]
    return (res_of_atom,
            [str(s) for s in residues.segids],
            [int(i) for i in residues.resids],
            [str(r) for r in residues.resnames])


def run(spec, out_dir, stride=None, write=True, max_frames=None):
    stride = int(stride or DEFAULT_STRIDE)
    rows = io.Rows(system=spec.id, stage=STAGE, stride=stride, resolved=True,
                   **EXTRA_COLUMNS)
    agg = Pooled()
    inputs = []
    splits = {}

    for md in spec.mds:
        # The STRIPPED inputs, and only those: an interface is made of solute,
        # and this is thirteen times less I/O for exactly the same geometry.
        top, trajs = md.inputs(solvated=False)
        if not trajs:
            continue
        inputs += [top] + list(trajs)
        u = _universe(top, trajs)
        sel = _groups(spec, u)

        rec = u.select_atoms(sel["receptor"].mdanalysis + " and " + S.HEAVY)
        lig = u.select_atoms(sel["ligand"].mdanalysis + " and " + S.HEAVY)
        for name, ag in (("receptor", rec), ("ligand", lig)):
            if not len(ag):                                # pragma: no cover
                raise S.EmptySelection(
                    f"{md.id}: the {name} has no heavy atoms by name. Every "
                    f"contact here is measured against `{S.HEAVY}` -- heavy "
                    f"atoms by NAME: the masses in these topologies are "
                    f"guessed (MDAnalysis reads the atom name CA as carbon) "
                    f"and a mass-based filter quietly returns a different set.")

        rec_of_atom, rec_seg, rec_resid, rec_resname = _residue_table(rec)
        lig_of_atom, lig_seg, lig_resid, lig_resname = _residue_table(lig)
        n_rec_res, n_lig_res = len(rec_seg), len(lig_seg)
        rec_seg_atom = np.asarray([str(s) for s in rec.segids])
        lig_seg_atom = np.asarray([str(s) for s in lig.segids])

        # WHICH FRAMES. Seams first: a restart re-emits a step with DIFFERENT
        # coordinates, and a duplicated frame silently lengthens every contact
        # episode that spans it by one. core/frames.py finds them from the DCD
        # headers rather than from the coordinates.
        timeline = F.build(md.id, list(md.timeline_inputs()),
                           framestep_ns=md.framestep_ns)
        seams = set(timeline.seams)
        n_raw = len(u.trajectory)
        if n_raw != len(timeline) + len(seams):
            seams = set()          # already deduplicated: drop nothing
        picked = [i for i in range(n_raw) if i not in seams][::stride]
        if max_frames:
            picked = picked[:int(max_frames)]
        if not picked:                                     # pragma: no cover
            continue

        step_ns = md.framestep_ns or timeline.framestep_ns
        spacing_ns = (step_ns * stride) if step_ns else None
        # The a-priori spacing verdict is CARRIED, not raised: the occupancies
        # are time averages a coarse stride does not harm, and discarding them
        # to protect the lifetimes would lose the cheap answers for the
        # expensive one. The message goes into the lifetime rows verbatim.
        life_note = None
        try:
            SP.check_spacing(f"{STAGE}:contact_lifetime", step_ns, stride,
                             LIFETIME_MIN_SPACING_NS)
        except SP.SpacingTooCoarse as exc:
            life_note = str(exc)

        split = evidence = entity_of = units = unit_index = None
        rec_prot = lig_prot = lig_unit_i = offset_of = pairs_per_frame = None
        # {pair key -> the sampled frames it was on}, kept sparse. A dense
        # (1,896 x 699) mask is 1.3 M bits per frame and about 0.15% of it is
        # ever True.
        seen = {}
        rec_any = np.zeros(n_rec_res, dtype=np.int64)
        lig_any = np.zeros(n_lig_res, dtype=np.int64)
        rec_partners = np.zeros(n_rec_res, dtype=np.int64)
        lig_partners = np.zeros(n_lig_res, dtype=np.int64)

        n_used = 0
        for ts in _sample(u, picked, stride):
            # The guard, not a correction: this stage's distances and
            # capped_distance's `box=` both assume 90/90/90, which these runs
            # are (`useFlexibleCell no`). NonOrthorhombic propagates on purpose.
            L = pbc.box_lengths(ts.dimensions)
            if L is None:
                raise ValueError(
                    f"{md.id} frame {picked[n_used]}: this frame carries no "
                    f"periodic cell. Every contact here is a minimum-image "
                    f"contact and there is nothing honest to do without one.")
            box = ts.dimensions
            rec_xyz = rec.positions
            lig_xyz = lig.positions

            if split is None:
                # The split, from the FIRST SAMPLED frame -- the first frame
                # this stage actually reads, so a short test and a production
                # run derive it the same way.
                split, evidence = split_assembly(top, rec_seg_atom, rec_xyz,
                                                 lig_seg_atom, lig_xyz, box)
                entity_of = split.entity_of
                units = [p["id"] for p in split]
                unit_index = {p["id"]: p["index"] for p in split}
                prot_of_seg = {s: p["id"] for p in split
                               for s in (list(p["receptor_segids"])
                                         + list(p["ligand_segids"]))}
                rec_prot = [prot_of_seg.get(s) for s in rec_seg]
                lig_prot = [prot_of_seg.get(s) for s in lig_seg]
                lig_unit_i = np.asarray(
                    [unit_index[p] if p else -1 for p in lig_prot])
                offset_of = protomer_offsets(rec_prot, unit_index, units,
                                             split.n)
                pairs_per_frame = np.zeros(len(units), dtype=np.int64)
                splits[md.id] = dict(split.as_dict(),
                                     pairing_evidence=dict(evidence),
                                     pairing_cutoff_A=PAIRING_A)
                for group in pooled_groups(rec_seg, lig_seg, entity_of):
                    for pid in units:
                        agg.register(group, md.id, pid)

            r_idx, l_idx = residue_pairs(rec_xyz, lig_xyz, rec_of_atom,
                                         lig_of_atom, n_lig_res, box)
            if len(r_idx):
                # The same pairs answer both histograms: how many partners a
                # residue holds (bincount over the pairs) and whether it is at
                # the interface at all (unique residues). One search, two
                # questions, and they cannot disagree.
                rec_partners += np.bincount(r_idx, minlength=n_rec_res)
                lig_partners += np.bincount(l_idx, minlength=n_lig_res)
                rec_any[np.unique(r_idx)] += 1
                lig_any[np.unique(l_idx)] += 1
                placed = lig_unit_i[l_idx]
                pairs_per_frame += np.bincount(placed[placed >= 0],
                                               minlength=len(units))
                for k in (r_idx * np.int64(n_lig_res) + l_idx).tolist():
                    seen.setdefault(k, []).append(n_used)
            n_used += 1

        if not n_used:                                     # pragma: no cover
            continue
        ns_sampled = round(n_used * spacing_ns, 3) if spacing_ns else None
        base = dict(replica=md.id, ns_per_frame=step_ns, n_frames=n_used,
                    ns_sampled=ns_sampled)

        # ----------------------------------------------------- the contact map
        # Only pairs seen at least once; absence IS zero (see the docstring).
        # Sorted by key, so the same trajectory writes the same bytes twice.
        per_unit_pairs = {name: 0 for name in units}
        per_unit_res = {name: (set(), set()) for name in units}
        lengths = {name: [] for name in units}
        censored = {name: 0 for name in units}
        for key in sorted(seen):
            j, k = divmod(key, n_lig_res)
            name = lig_prot[k]
            ui = int(lig_unit_i[k])
            off = int(offset_of[ui, j]) if ui >= 0 else -1
            frames_on = seen[key]
            occ = len(frames_on) / n_used
            subset = UNPAIRED if off < 0 else (COGNATE if off == 0 else CROSS)
            e_rec = entity_of.get(rec_seg[j], rec_seg[j])
            e_lig = entity_of.get(lig_seg[k], lig_seg[k])
            side = dict(segid=rec_seg[j], resid=rec_resid[j],
                        resname=rec_resname[j], entity=e_rec,
                        partner_resid=lig_resid[k],
                        partner_resname=lig_resname[k], partner_entity=e_lig,
                        protomer_offset=(off if off >= 0 else None),
                        subset=subset)
            rows.add(observable="contact_occupancy", scope="pair",
                     scope_id=f"{rec_seg[j]}:{rec_resid[j]}|"
                              f"{lig_seg[k]}:{lig_resid[k]}",
                     value=occ, unit="fraction", protomer=name,
                     partner_segid=lig_seg[k], **side, **base)
            if subset == UNPAIRED:
                # Not a repeating unit of this assembly, so there is nothing to
                # pool it with. The row above still says that it happened.
                continue
            per_unit_pairs[name] += 1
            per_unit_res[name][0].add(j)
            per_unit_res[name][1].add(k)
            # Pooled by (entity, resid, entity, resid, offset) -- resid WITHIN
            # ENTITY, never segid spelling -- and the offset keeps a
            # neighbouring protomer's contact out of the cognate map's average
            # while still pooling it with the same neighbour contact in the
            # other two protomers. BOTH segids are replaced by their entity
            # labels on the way in: an aggregate row that kept the segid of
            # whichever protomer happened to be written first would read as a
            # statement about that one copy, and it is about all of them.
            agg.add("pair", (e_rec, rec_resid[j], e_lig, lig_resid[k], off),
                    md.id, name, occ,
                    observable="contact_occupancy", unit="fraction",
                    scope_id=f"{e_rec}:{rec_resid[j]}|{e_lig}:{lig_resid[k]}"
                             f"@offset{off}",
                    **dict(side, segid=e_rec, partner_segid=e_lig))
            # Episodes, for the lifetime. The gaps are the point: a pair that is
            # on, off for one frame, then on again is TWO contacts, and
            # geom.run_lengths is the kernel every lifetime in this toolkit
            # uses, so contact lifetime, ion dwell and water residence cannot
            # come to different conclusions about one mask.
            first = frames_on[0]
            mask = np.zeros(frames_on[-1] - first + 1, dtype=bool)
            mask[np.asarray(frames_on, dtype=np.int64) - first] = True
            for start, length in geom.run_lengths(mask):
                if first + start == 0 or first + start + length == n_used:
                    censored[name] += 1
                lengths[name].append(int(length))

        # ------------------------------------------------- the two histograms
        for side_name, segs, resids, resnames, any_, partners, prot in (
                ("receptor", rec_seg, rec_resid, rec_resname, rec_any,
                 rec_partners, rec_prot),
                ("ligand", lig_seg, lig_resid, lig_resname, lig_any,
                 lig_partners, lig_prot)):
            for j in np.flatnonzero(any_).tolist():
                segid = segs[j]
                entity = entity_of.get(segid, segid)
                where = dict(scope="residue", scope_id=f"{segid}:{resids[j]}",
                             subset=side_name, segid=segid, resid=resids[j],
                             resname=resnames[j], entity=entity,
                             protomer=prot[j], **base)
                for obs, val, un in (("residue_occupancy", any_[j] / n_used,
                                      "fraction"),
                                     ("partners_per_frame",
                                      partners[j] / n_used, "count")):
                    rows.add(observable=obs, value=val, unit=un, **where)
                    if prot[j] is None:
                        continue          # excluded from the split: see above
                    # Pooled over the COPIES OF THIS ENTITY, one per protomer:
                    # a gp120 residue is averaged over the three gp120 copies,
                    # not over every receptor segment in the system.
                    agg.add(f"residue:{side_name}:{entity}:{obs}",
                            (entity, resids[j]), md.id, prot[j], val,
                            observable=obs, unit=un, subset=side_name,
                            scope_id=f"{entity}:{resids[j]}", segid=entity,
                            resid=resids[j], resname=resnames[j],
                            entity=entity)

        # ---------------------------------------------- per protomer, per run
        scale = spacing_ns or 1.0
        life_unit = "ns" if spacing_ns else "frame"
        for p in split:
            name = p["id"]
            run_at = dict(scope="run", scope_id=f"{md.id}:{name}",
                          protomer=name, **base)
            for obs, val, un, sub in (
                    ("contact_pairs", per_unit_pairs[name], "count", None),
                    ("contacts_per_frame",
                     pairs_per_frame[unit_index[name]] / n_used, "count", None),
                    ("interface_residues", len(per_unit_res[name][0]), "count",
                     "receptor"),
                    ("interface_residues", len(per_unit_res[name][1]), "count",
                     "ligand")):
                rows.add(observable=obs, value=val, unit=un, subset=sub,
                         **run_at)
                agg.add("run", (obs, sub or ""), md.id, name, val,
                        observable=obs, unit=un, subset=sub,
                        scope_id=f"{obs}:{sub}" if sub else obs)
            st = lifetime_stats(lengths[name])
            if st is not None:
                ok, note = lifetime_floor(st, spacing_ns, life_note)
                for obs, val in (("contact_lifetime_mean", st["mean"] * scale),
                                 ("contact_lifetime_median",
                                  st["median"] * scale)):
                    rows.add(observable=obs, value=val, unit=life_unit,
                             resolved=ok, floor_note=note, **run_at)
                    agg.add("run", (obs, ""), md.id, name, val,
                            observable=obs, unit=life_unit, resolved=ok,
                            floor_note=note, scope_id=obs)
                # The guard's own evidence, always resolved: this is the number
                # that decides the two above, so it has to be readable alone.
                for obs, val, un in (
                        ("contact_lifetime_single_frame_fraction",
                         st["single_frame_fraction"], "fraction"),
                        ("contact_lifetime_episodes", st["n"], "count"),
                        # Episodes still open at the first or the last sampled
                        # frame are censored: their true length is longer than
                        # measured, so the mean above is a LOWER bound. Counted
                        # rather than dropped -- dropping them removes the long
                        # ones and biases the mean up.
                        ("contact_lifetime_censored", censored[name],
                         "count")):
                    rows.add(observable=obs, value=val, unit=un, **run_at)
                    agg.add("run", (obs, ""), md.id, name, val,
                            observable=obs, unit=un, scope_id=obs)
            # What the split was measured on. NOT aggregated -- it is a
            # diagnostic about this replica's first frame, and an average of it
            # would mean nothing -- but written, so "which chains are protomer
            # B, and on how much evidence" is answerable from the table alone.
            # A ZERO HERE IS THE INTERESTING CASE: it means core.protomers
            # placed that chain by elimination rather than by contact, and
            # `split.rule` and the notes in the manifest say so.
            for segid in p["ligand_segids"]:
                rows.add(observable="pairing_evidence", subset=split.rule,
                         label=segid,
                         value=evidence.get(segid, [0] * split.n)[
                             unit_index[name]],
                         unit="count", partner_segid=segid,
                         partner_entity=entity_of.get(segid, segid),
                         cutoff_A=PAIRING_A, param_name="pairing_rule",
                         param_value=split.rule, **run_at)

    # ------------------------------------------------------------- aggregates
    # Pooled over protomers AND replicas: `value` the mean, `sd` the spread
    # across the units, `n` their count, `sem` over REPLICAS only. `param_name`
    # says that a row pools protomers, which core/io.py requires of any row
    # that does.
    for group in sorted(agg.values):
        for key, meta, values, replicas in agg.items(group):
            st = pool(values, replicas)
            meta = dict(meta)
            note = meta.pop("floor_note", None)
            ok = bool(meta.pop("resolved", True))
            scope_id = meta.pop("scope_id",
                                "|".join(str(k) for k in key if str(k) != ""))
            n_units = max(1, st["n"] // max(1, st["n_replicas"]))
            if st["sem"] is None:
                extra = ("one replica: there is no standard error over "
                         "replicas, and `sd` is over protomers that share a "
                         "box and a thermostat")
                note = f"{note} {extra}" if note else extra
                ok = False
            rows.add(scope="aggregate", scope_id=scope_id, value=st["value"],
                     sd=st["sd"], sem=st["sem"], n=st["n"], resolved=ok,
                     floor_note=note, param_name="pooled_over",
                     param_value=f"{n_units} protomers x "
                                 f"{st['n_replicas']} replicas", **meta)

    # The cutoffs and the split travel with the result, so a table can never be
    # read without the geometry that produced it -- and core.protomers' own
    # `notes` are copied verbatim, because "the contact pairing disagrees with
    # the copy index" is exactly the sentence a person needs and exactly the one
    # a summary would drop. io.needs_run() compares stage, spec_hash, stride and
    # inputs only, so editing a cutoff here does NOT invalidate an existing
    # result: delete the stem's files to force the re-run.
    man = io.manifest(STAGE, spec, inputs, stride=stride,
                      extra={"cutoffs_A": {"contact": CONTACT_A,
                                           "pairing": PAIRING_A},
                             "max_frames": max_frames,
                             "protomer_split": splits,
                             "pairing_source":
                                 "core.protomers.split, on the first sampled "
                                 "frame, from minimum-image contact counts"})
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows
