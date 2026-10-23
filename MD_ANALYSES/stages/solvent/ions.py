#!/usr/bin/env python3
"""
ions.py -- where the counterions sit, how they are hydrated, how long they stay.

WHY THIS STAGE EXISTS. Nothing in this workspace measures ions. `composition`
counts them out of the PSF -- 668 SOD and 697 CLA for 1F4 -- and that is the
whole of what is known about them: a 0.15 M salt bath is declared once and
never looked at again. These are glycosylated antigens whose shield is dense
and carries anionic sialic acids, and what a charged epitope presents to a
partner depends on which counterions are camped on it and for how long. A Na+
that holds one aspartate for 40 ns is part of the structure the reader is
being shown; one that brushes it for a single frame is noise. Without this
stage the two are indistinguishable, because both are absent.

FOUR OBSERVABLES, ONE PASS.

  ion_contacts         per residue, the fraction of sampled frames with an ion
                       of that species within 4.0 A of its heavy atoms: the
                       map of where the salt actually goes.
  coordination_number  mean water oxygens within the first (3.2 A) and second
                       (5.6 A) shell radii of each species. Both are
                       CUMULATIVE counts from the ion, so the second includes
                       the first. This one has a known answer, which is the
                       point of measuring it: Na+ comes out near 5-6 in bulk
                       TIP3P, and on 1F4/1-R1 it is 5.79. A first shell far
                       from that means the water selection, the cutoff or the
                       minimum image is wrong -- and it would be wrong for the
                       other three observables too. This is the stage's own
                       control, not decoration.
  ion_partition        every ion classified each frame as glycan_shell /
                       protein_shell / both / bulk, as mean fractions AND as
                       enrichments.
  ion_dwell            how long an ion stays at one residue, from
                       geom.run_lengths on the per-ion/per-residue contact
                       mask -- the same kernel as every other lifetime here,
                       so ion dwell and water residence cannot end up
                       disagreeing about what a run is.

COUNTS ARE NOT COMPARABLE BETWEEN SYSTEMS; ENRICHMENTS ARE. "410 ions in the
glycan shell" is mostly a statement about box size and salt count -- the OPEN,
OCCLUDED-OPEN and CLOSED boxes differ by tens of thousands of waters -- so
each class also gets its count divided by what the bulk would put in the same
shell volume. THE SHELL VOLUME IS ESTIMATED WITH WATER AS THE PROBE: the water
oxygens in that same 4.0 A shell, divided by the bulk water number density
(N_water / V_box). Its one assumption is that shell water sits at bulk
density, which is wrong by the 10-20% that interfacial water is denser -- so
the enrichment of a tight first shell is an under-estimate. That assumption is
stated here and recorded in every row (`volume_probe =
water_at_bulk_density`), rather than left for a reader to rediscover. The
estimate's virtue is that V_box cancels algebraically:

    enrichment = (ions in class / ions total) / (waters in class / waters total)

so nothing depends on a box volume the `.xst` and the DCD headers might
disagree about, and a value near 1.0 means "as much salt here as anywhere".

CONVENTIONS THAT ARE NOT NEGOTIABLE, each one a bug this lab has paid for:

  * water is counted by its OXYGEN, `name OH2`. `name O` in a CHARMM topology
    returns every backbone carbonyl, and a coordination number built from
    those is plausible, wrong and silent.
  * heavy atoms come from `not name H*` (core.selections.HEAVY), by NAME. The
    masses in these topologies are broken.
  * NOTHING IS WRAPPED. `cellOrigin` for 1F4 is (169.7, 169.0, 158.7), so
    MDAnalysis's `wrap()` and `pack_into_box()` move atoms about the wrong
    reference and tear the solute apart. Every distance here is a minimum-image
    distance computed by `capped_distance(..., box=ts.dimensions)`, which
    applies the convention internally and never touches a coordinate;
    core.pbc.box_lengths is the guard that the cell really is orthorhombic,
    which that arithmetic assumes.
  * `segid`, never `chainID`: chainID collapses the three antibody copies and
    all 75 glycan segments onto one letter each.
  * NO MASS IS EVER READ. MDAnalysis's guesser maps the atom name `CA` to
    carbon, so a calcium ion would weigh 12. Every observable here is pure
    geometry and there is nothing for a wrong mass to corrupt -- which is a
    property worth stating, because a centre-of-mass version of any of these
    numbers would be quietly wrong on these topologies.
  * the SOLVATED trajectory, via `md.inputs(solvated=True)`. A stripped file
    has no ions at all, and that method raises rather than letting this stage
    report 1,365 confident zeros.

COST. The solvated replicas are 23-30 GB on a ~100 MB/s share, so the
trajectory is read ONCE and all four observables share the walk. Three
neighbour searches per frame cost ~0.5 s against ~0.15 s of I/O, and the
expensive one is the water probe -- it buys the only figure that can be
compared between two systems with different boxes, so it is paid for
knowingly. `max_frames` STOPS the read rather than skipping the rest of it,
which is the difference between a 15-frame test costing seconds and costing
four minutes.
"""

from core import threads                                   # noqa: F401

import numpy as np

from core import frames as F
from core import geom, io, pbc
from core import selections as S
from core import spec as SP

STAGE = "ions"

#: Declared for `ion_dwell`, the only observable here with dynamic content: an
#: ion-residue episode is over in a few ns, so at coarser than 0.5 ns per frame
#: the mean dwell is a statement about the stride and not about the ion. The
#: other three are time averages and need no particular spacing, which is why a
#: stride that violates this does NOT kill the stage -- it marks the dwell rows
#: `resolved=False` and copies check_spacing's own message into `floor_note`.
#: Note the consequence, deliberately: at DEFAULT_STRIDE, which is chosen for
#: the occupancy observables, dwell is NOT resolved. Dwell wants stride 1-2.
MIN_SPACING_NS = 0.5
DEFAULT_STRIDE = 10

#: Heavy-atom contact distance. The 4.0 A this workspace's TCL scripts already
#: use for contacts, and wide enough to hold a contact ion pair without
#: reaching into the solvent-separated one.
CONTACT_A = 4.0

#: The two shell radii, CUMULATIVE: `second_shell` counts everything inside
#: 5.6 A, first shell included.
#:
#: 3.2 A is the first minimum of the Na+-O RDF and the right boundary for a
#: cation. IT IS NOT THE RIGHT ONE FOR Cl-, whose ionic radius is larger: on
#: 1F4/1-R1, frame 0, the cumulative count per ion runs
#:
#:            2.8 A    3.2 A    3.5 A    3.8 A    4.0 A
#:     SOD     5.52     5.75     5.95     6.50     7.34
#:     CLA     0.00     3.36     6.50     7.59     8.25
#:
#: -- Na+ is flat from 2.8 to 3.5 (a real shell, 5.5-6.0), while Cl- has no
#: water at all inside 2.8 A and its own first minimum near 3.8-3.9 A. So the
#: CLA first-shell number under this COMMON cutoff is a definition ("water
#: oxygens within 3.2 A"), not the whole shell, and it halves it. The cutoff
#: therefore travels with the value in `param_value` and can never be read
#: apart from it. Per-species cutoffs were the alternative and are worse: two
#: numbers in one column measured at two radii, comparable to nothing.
#: First-shell radius PER SPECIES, because one number is wrong for both. The
#: cumulative water count against cutoff, measured on 1F4/1-R1 frame 0:
#:
#:     cutoff   2.8   3.2   3.5   3.8   3.9   4.0   5.6
#:     SOD     5.52  5.75  5.95  6.50  6.89  7.34 25.52
#:     CLA     0.00  3.36  6.50  7.59  7.88  8.25 23.29
#:
#: Na+ is flat from 2.8 to 3.5 A -- a real shell, and 3.2 sits inside it.
#: Cl- has NO water inside 2.8 A and its own first minimum is near 3.8-3.9 A,
#: so 3.2 cuts its first peak in half and reports 3.4 where the shell holds
#: about 7.6. Reporting both species at one radius made them look comparable
#: while making one of them wrong; each now gets the radius its own g(r)
#: minimum defines, and the radius travels with every row it produced.
#: (CHARMM Cl- in TIP3P is known to be somewhat over-coordinated, so ~7.6-7.9
#: is the expected value for this force field, not a defect here.)
FIRST_SHELL_BY_SPECIES = {"SOD": 3.2, "NA": 3.2, "POT": 3.5, "K": 3.5,
                          "CAL": 3.1, "CA": 3.1, "MG": 2.6, "ZN": 2.4,
                          "CLA": 3.9, "CL": 3.9, "BR": 4.0, "IOD": 4.3}
#: Used when a species is not in that table: a person must then read the
#: `first_shell_A` column rather than assume.
FIRST_SHELL_A = 3.2
SECOND_SHELL_A = 5.6


def first_shell_for(species):
    """The first-shell radius for one ion species, and whether it was known."""
    key = str(species).strip().upper()
    if key in FIRST_SHELL_BY_SPECIES:
        return FIRST_SHELL_BY_SPECIES[key], True
    return FIRST_SHELL_A, False

#: The water oxygen by a name that can only be one (CHARMM OH2, GROMACS OW),
#: spelled out because `name O` is the trap.
WATER_OXYGEN = "water and name " + " ".join(S.WATER_O_NAMES)

#: What an ion can be next to. Mutually exclusive and exhaustive, so the four
#: fractions sum to exactly 1 and `bulk` is a measurement rather than a
#: leftover.
PARTITION_CLASSES = ("glycan_shell", "protein_shell", "both", "bulk")

#: Above this share of one-frame episodes the dwell distribution is against the
#: sampling floor: half of everything you saw was over before the next frame,
#: so the mean is the frame spacing and not a residence time.
SINGLE_FRAME_MAX = 0.5


def _universe(topology, trajectories):
    import MDAnalysis as mda
    return mda.Universe(topology, *trajectories)


def _groups(spec, u):
    """Every selection this stage needs, checked before a frame is read.

    An empty selection is an error, never an empty result: "no ion came within
    4 A of anything" and "the ion selection matched nothing" produce the same
    table and are different facts.
    """
    want = {"ions": spec.selection("ions"),
            "solute": spec.selection("solute"),
            "glycans": spec.selection("glycans")}
    res, problems = S.resolve_all(u, want, allow_empty=("glycans",))
    if problems:
        raise S.EmptySelection("; ".join(problems))
    try:
        res["water_oxygen"] = S.resolve(u, "water_oxygen", WATER_OXYGEN)
    except S.EmptySelection as exc:
        raise S.EmptySelection(
            f"{exc} This stage counts water by its OXYGEN, `name OH2`, which "
            f"is what CHARMM calls it. Do not relax that to `name O`: in a "
            f"CHARMM topology it returns every backbone carbonyl instead, and "
            f"the coordination numbers built from them look entirely "
            f"reasonable.") from exc
    return res


def shell_classes(xyz, solute_xyz, is_glycan, box, cutoff=CONTACT_A):
    """Which shell each point of `xyz` is in, and the pairs that decided it.

    ONE box-aware neighbour search answers three questions at once -- which
    residues a point touches, whether it touches glycan, whether it touches
    the rest of the solute -- so the residue contacts, the partition and the
    dwell mask all come from the same pairs and cannot disagree.

    `box` goes straight to `capped_distance`, which applies the minimum image
    itself and leaves every coordinate exactly where it was; that is the whole
    reason nothing here has to be wrapped.

    Returns (class code per point, point indices, solute atom indices).
    """
    from MDAnalysis.lib.distances import capped_distance
    pairs = capped_distance(np.asarray(xyz, dtype=np.float32),
                            np.asarray(solute_xyz, dtype=np.float32),
                            cutoff, box=box, return_distances=False)
    i = pairs[:, 0] if len(pairs) else np.zeros(0, dtype=np.int64)
    j = pairs[:, 1] if len(pairs) else np.zeros(0, dtype=np.int64)
    on_glycan = np.zeros(len(xyz), dtype=bool)
    on_other = np.zeros(len(xyz), dtype=bool)
    if len(i):
        gly = is_glycan[j]
        on_glycan[i[gly]] = True
        on_other[i[~gly]] = True
    code = np.full(len(xyz), PARTITION_CLASSES.index("bulk"), dtype=np.int8)
    code[on_glycan & ~on_other] = PARTITION_CLASSES.index("glycan_shell")
    code[on_other & ~on_glycan] = PARTITION_CLASSES.index("protein_shell")
    code[on_glycan & on_other] = PARTITION_CLASSES.index("both")
    return code, i, j


def coordination(ion_xyz, water_o_xyz, box, first=FIRST_SHELL_A,
                 second=SECOND_SHELL_A):
    """Water oxygens per ion within the first and the second shell radius.

    Both counts are CUMULATIVE -- the second includes the first -- because
    that is what a running coordination number is and what the RDF minima that
    set the radii refer to. An annulus would need saying so in the row.

    One search at the OUTER radius, then a threshold on the distances it has
    already returned: searching twice would double the per-frame cost for an
    answer the first search already contains, and two searches could drift
    apart into a first shell larger than the second.
    """
    from MDAnalysis.lib.distances import capped_distance
    n = len(ion_xyz)
    pairs, d = capped_distance(np.asarray(ion_xyz, dtype=np.float32),
                               np.asarray(water_o_xyz, dtype=np.float32),
                               second, box=box, return_distances=True)
    if not len(pairs):
        return np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
    idx = pairs[:, 0]
    # `first` may be one radius or one PER ION, so that Na+ and Cl- can each
    # be counted inside their own g(r) minimum in a single search: comparing
    # every pair against its own ion's radius costs one extra fancy-index.
    first_arr = np.broadcast_to(np.asarray(first, dtype=np.float64), (n,))
    return (np.bincount(idx[d <= first_arr[idx]], minlength=n).astype(np.int64),
            np.bincount(idx, minlength=n).astype(np.int64))


def enrichment(n_class, n_total, water_class, water_total):
    """Ions in a shell against what bulk density would put in that shell.

    The shell volume is estimated from the water it holds, so the box volume
    cancels and what is left is the ions' share of the class over water's
    share of it. Returns None when the probe is empty: a shell with no water
    in it has no volume estimate, and a ratio over zero has put an `inf` on a
    figure in this workspace before.
    """
    if not n_total or not water_class or not water_total:
        return None
    return (n_class / n_total) / (water_class / water_total)


def dwell_stats(lengths):
    """mean, median and the one-frame share of a set of episode lengths."""
    a = np.asarray(lengths, dtype=np.float64)
    if not len(a):
        return None
    return {"n": int(len(a)), "mean": float(a.mean()),
            "median": float(np.median(a)),
            "single_frame_fraction": float((a == 1).mean())}


def dwell_floor(stats, spacing_ns, prior_note=None):
    """Whether a dwell time is resolved, and the note that says why it is not.

    TWO INDEPENDENT GUARDS, because either alone can be fooled. The a-priori
    one is core.spec.check_spacing's verdict on the frame spacing, passed in as
    `prior_note`. The measured one is here: a distribution in which most
    episodes last exactly one frame has been cut off by the sampling, and its
    mean is the frame spacing rather than a residence time -- which happens at
    a spacing check_spacing would pass, whenever the real exchange is faster
    than anyone expected.
    """
    note = prior_note
    if stats["single_frame_fraction"] >= SINGLE_FRAME_MAX:
        extra = (f"{stats['single_frame_fraction']:.0%} of {stats['n']} "
                 f"episodes last a single frame"
                 + (f" ({spacing_ns:g} ns)" if spacing_ns else "")
                 + "; the mean dwell is then the frame spacing, not the ion's "
                   "residence time")
        note = f"{note} {extra}" if note else extra
    return note is None, note


def _sample(u, picked, stride):
    """Iterate the wanted frames the cheaper of the two available ways.

    Above about stride 3 a seek on this NFS share (~215 ms) is cheaper than
    reading the frames in between (~92 ms each), so the frames are reached by
    random access; at stride 1 or 2 the sequential walk wins and the unwanted
    frames are simply not used. Either way the walk STOPS at the last wanted
    frame -- a `max_frames` that skipped the remainder instead would still
    read all 23 GB, which is the whole cost of this stage.
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


def run(spec, out_dir, stride=None, write=True, max_frames=None):
    stride = int(stride or DEFAULT_STRIDE)
    rows = io.Rows(system=spec.id, stage=STAGE, stride=stride, resolved=True)
    inputs = []

    for md in spec.mds:
        # The solvated inputs and nothing else. A stripped trajectory holds no
        # ions, and this raises instead of reporting confident zeros.
        top, trajs = md.inputs(solvated=True)
        if not trajs:
            continue
        inputs += [top] + list(trajs)
        u = _universe(top, trajs)
        sel = _groups(spec, u)

        # One atom per ion: core.topology.ION is the monatomic set by
        # construction ("anything polyatomic is a ligand"), so an ion's atom IS
        # the ion and no residue grouping is needed. Its species is its resname.
        ion_ag = u.select_atoms(sel["ions"].mdanalysis)
        water_o = u.select_atoms(sel["water_oxygen"].mdanalysis)
        solute = u.select_atoms(sel["solute"].mdanalysis + " and " + S.HEAVY)
        if not len(solute):                                # pragma: no cover
            raise S.EmptySelection(
                f"{md.id}: the solute has no heavy atoms by name; every "
                f"contact here is measured against `{S.HEAVY}`")

        species = np.asarray([str(r) for r in ion_ag.resnames])
        names = sorted(set(species.tolist()))
        sp_of_ion = np.searchsorted(np.asarray(names), species)
        species_mask = [species == n for n in names]
        n_per_species = np.asarray([int(m.sum()) for m in species_mask])
        n_sp = len(names)

        # `protein_shell` means "the solute that is not glycan": for these
        # systems that is protein, and for one with a small-molecule ligand it
        # would include the ligand. The two sets partition the solute exactly,
        # which is what makes `both` and `bulk` mean anything.
        if "glycans" in sel and sel["glycans"].n_atoms:
            gly = u.select_atoms(sel["glycans"].mdanalysis + " and " + S.HEAVY)
            is_glycan = np.isin(solute.indices, gly.indices)
        else:
            is_glycan = np.zeros(len(solute), dtype=bool)

        # Atom -> residue as a compact index, computed once. Doing it per frame
        # from the AtomGroup was the cost that made the first draft unusable.
        res_ids, res_of_atom = np.unique(solute.resindices, return_inverse=True)
        res_of_atom = np.ravel(res_of_atom)
        n_res = len(res_ids)
        residues = u.residues[res_ids]
        seg_of = [str(s) for s in getattr(residues, "segids", [""] * n_res)]
        resid_of = [int(i) for i in
                    getattr(residues, "resids", list(range(n_res)))]
        resname_of = [str(r) for r in
                      getattr(residues, "resnames", [""] * n_res)]

        # WHICH FRAMES. Seams first: a restart re-emits a step with DIFFERENT
        # coordinates, and a duplicated frame silently lengthens every dwell
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
        # The a-priori spacing verdict is CARRIED, not raised: three of the four
        # observables here are time averages a coarse stride does not harm, and
        # discarding them to protect the fourth would lose the cheap answers for
        # the expensive one. The message goes into the dwell rows verbatim.
        dwell_note = None
        try:
            SP.check_spacing("ion_dwell", step_ns, stride, MIN_SPACING_NS)
        except SP.SpacingTooCoarse as exc:
            dwell_note = str(exc)

        contact_frames = np.zeros((n_sp, n_res), dtype=np.int64)
        hit = np.zeros((n_sp, n_res), dtype=bool)
        coord_sum = np.zeros((n_sp, 2))
        coord_sq = np.zeros((n_sp, 2))
        coord_n = np.zeros(n_sp)
        part_ions = np.zeros((n_sp, len(PARTITION_CLASSES)), dtype=np.int64)
        part_water = np.zeros(len(PARTITION_CLASSES), dtype=np.int64)
        vol_shell = np.zeros(len(PARTITION_CLASSES))
        # The dwell mask, kept sparse: {(ion, residue) key -> sampled frames}.
        # A dense (ion x residue x frame) array is 1,365 x 3,330 x 2,520 bits,
        # and about 0.2% of it would be True.
        episodes = {}

        # Each ion carries the radius its own species' first minimum defines.
        first_by_species = [first_shell_for(nm) for nm in names]
        first_per_ion = np.empty(len(sp_of_ion), dtype=np.float64)
        for s_i, (radius, _known) in enumerate(first_by_species):
            first_per_ion[sp_of_ion == s_i] = radius

        n_used = 0
        for ts in _sample(u, picked, stride):
            # The guard, not a correction: this stage's arithmetic and
            # capped_distance's `box=` both assume 90/90/90, which these runs
            # are (`useFlexibleCell no`). NonOrthorhombic propagates on purpose.
            L = pbc.box_lengths(ts.dimensions)
            if L is None:
                raise ValueError(
                    f"{md.id} frame {picked[n_used]}: this frame carries no "
                    f"periodic cell. Every distance here is a minimum-image "
                    f"distance and there is nothing honest to do without one.")
            box = ts.dimensions
            ion_xyz = ion_ag.positions
            solute_xyz = solute.positions
            water_xyz = water_o.positions

            code_i, pi, pj = shell_classes(ion_xyz, solute_xyz, is_glycan, box)
            # A residue counts ONCE per frame however many of its atoms an ion
            # is near, and once however many ions of the species are near it:
            # this observable is "was there a SOD here", not "how many".
            hit[:] = False
            if len(pi):
                res = res_of_atom[pj]
                hit[sp_of_ion[pi], res] = True
                for key in np.unique(pi.astype(np.int64) * n_res + res).tolist():
                    episodes.setdefault(key, []).append(n_used)
            contact_frames += hit
            np.add.at(part_ions, (sp_of_ion, code_i), 1)

            # The volume probe: the SAME classification applied to water, so
            # the shell whose volume is estimated is by construction the shell
            # the ions were counted in.
            code_w, _, _ = shell_classes(water_xyz, solute_xyz, is_glycan, box)
            per_class = np.bincount(code_w, minlength=len(PARTITION_CLASSES))
            part_water += per_class
            vol_shell += per_class / (len(water_o) / float(L.prod()))

            n1, n2 = coordination(ion_xyz, water_xyz, box,
                                  first=first_per_ion)
            for s, m in enumerate(species_mask):
                coord_sum[s] += (n1[m].sum(), n2[m].sum())
                coord_sq[s] += ((n1[m] ** 2).sum(), (n2[m] ** 2).sum())
            coord_n += n_per_species
            n_used += 1

        base = dict(replica=md.id, scope="run", scope_id=md.id,
                    ns_per_frame=step_ns, n_frames=n_used,
                    ns_sampled=(round(n_used * spacing_ns, 3)
                                if spacing_ns else None))

        # ------------------------------------------------------- ion_contacts
        # Only residues an ion reached at least once get a row. A residue
        # absent from this table was never within CONTACT_A of an ion of that
        # species in the sampled frames -- i.e. absence IS zero, stated here
        # because 3,330 residues x 2 species of mostly-zero rows would bury
        # the few hundred that carry the answer.
        for s, name in enumerate(names):
            for j in np.flatnonzero(contact_frames[s]).tolist():
                rows.add(observable="ion_contacts", subset=name,
                         value=contact_frames[s, j] / n_used, unit="fraction",
                         replica=md.id, scope="residue",
                         scope_id=int(res_ids[j]), segid=seg_of[j],
                         resid=resid_of[j], resname=resname_of[j],
                         ns_per_frame=step_ns, n_frames=n_used,
                         param_name="cutoff_A", param_value=CONTACT_A)

        # -------------------------------------------------- coordination_number
        for s, name in enumerate(names):
            if not coord_n[s]:
                continue                                   # pragma: no cover
            radius, known = first_by_species[s]
            for k, (label, cut) in enumerate((("first_shell", radius),
                                              ("second_shell", SECOND_SHELL_A))):
                mean = coord_sum[s, k] / coord_n[s]
                var = coord_sq[s, k] / coord_n[s] - mean * mean
                rows.add(observable="coordination_number", subset=name,
                         label=label, value=mean, unit="count",
                         param_name="cutoff_A", param_value=cut, **base)
                # The spread over ions and frames: a first shell with a large
                # sd is not a shell, and the mean should not be quoted as one.
                rows.add(observable="coordination_number_sd", subset=name,
                         label=label, value=float(np.sqrt(max(0.0, var))),
                         unit="count", param_name="cutoff_A", param_value=cut,
                         **base)

        # ------------------------------------------------------ ion_partition
        water_total = int(part_water.sum())
        for s, name in enumerate(names):
            total = int(part_ions[s].sum())
            if not total:
                continue                                   # pragma: no cover
            for c, cls in enumerate(PARTITION_CLASSES):
                rows.add(observable="ion_partition", subset=name, label=cls,
                         value=part_ions[s, c] / total, unit="fraction",
                         param_name="cutoff_A", param_value=CONTACT_A, **base)
                e = enrichment(int(part_ions[s, c]), total,
                               int(part_water[c]), water_total)
                rows.add(observable="ion_partition_enrichment", subset=name,
                         label=cls, value=e, unit="ratio",
                         param_name="volume_probe",
                         param_value="water_at_bulk_density",
                         resolved=e is not None,
                         floor_note=(None if e is not None else
                                     "no water in this shell to estimate its "
                                     "volume from, so there is no density to "
                                     "compare against"),
                         **base)
        # The volume the enrichment divided by, written down so the assumption
        # can be checked against a geometric estimate rather than believed.
        for c, cls in enumerate(PARTITION_CLASSES):
            rows.add(observable="shell_volume", label=cls,
                     value=vol_shell[c] / n_used, unit="A^3",
                     param_name="volume_probe",
                     param_value="water_at_bulk_density", **base)

        # ---------------------------------------------------------- ion_dwell
        lengths = {s: [] for s in range(n_sp)}
        censored = np.zeros(n_sp, dtype=np.int64)
        for key, fr in episodes.items():
            s = int(sp_of_ion[key // n_res])
            first = fr[0]
            mask = np.zeros(fr[-1] - first + 1, dtype=bool)
            mask[np.asarray(fr, dtype=np.int64) - first] = True
            # THE kernel. Every lifetime in this toolkit comes through
            # geom.run_lengths, so contact lifetimes, water residence and ion
            # dwell cannot come to different conclusions about one mask.
            for start, length in geom.run_lengths(mask):
                if first + start == 0 or first + start + length == n_used:
                    censored[s] += 1
                lengths[s].append(int(length))

        scale = spacing_ns or 1.0
        unit = "ns" if spacing_ns else "frame"
        for s, name in enumerate(names):
            st = dwell_stats(lengths[s])
            if st is None:
                continue                                   # pragma: no cover
            ok, note = dwell_floor(st, spacing_ns, dwell_note)
            for obs, val in (("ion_dwell_mean", st["mean"] * scale),
                             ("ion_dwell_median", st["median"] * scale)):
                rows.add(observable=obs, subset=name, value=val, unit=unit,
                         resolved=ok, floor_note=note, param_name="cutoff_A",
                         param_value=CONTACT_A, **base)
            # The guard's own evidence, always resolved: this number is exactly
            # what decides the two above, so it must be readable on its own.
            rows.add(observable="ion_dwell_single_frame_fraction", subset=name,
                     value=st["single_frame_fraction"], unit="fraction", **base)
            rows.add(observable="ion_dwell_episodes", subset=name,
                     value=st["n"], unit="count", **base)
            # Episodes still open at the first or the last sampled frame are
            # censored: their true length is unknown and longer than measured,
            # so the mean above is a LOWER bound. Counted rather than dropped,
            # because dropping them removes the long ones and biases it up.
            rows.add(observable="ion_dwell_censored", subset=name,
                     value=int(censored[s]), unit="count", **base)

    # The cutoffs are recorded so a table can never be read without the
    # geometry that produced it. io.needs_run() compares stage, spec_hash,
    # stride and inputs only, so editing a cutoff here does NOT invalidate an
    # existing result: delete the stem's files to force the re-run.
    man = io.manifest(STAGE, spec, inputs, stride=stride,
                      extra={"cutoffs_A": {"contact": CONTACT_A,
                                           "first_shell": FIRST_SHELL_A,
                                           "second_shell": SECOND_SHELL_A},
                             "max_frames": max_frames,
                             "volume_probe": "water_at_bulk_density"})
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows
