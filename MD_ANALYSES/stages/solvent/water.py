#!/usr/bin/env python3
"""
water.py -- the interfacial water shell, the bridges it makes, and how long it
stays. Ported from ABZYME/3-OCCLUDED-OPEN/CDR/ANALYSIS/70_water_shell.py,
71_water_bridges.py and 72_water_residence.py, which are correct and were
verified against an independent brute force; what changed here is only what made
them unportable.

WHY THESE THREE TOGETHER. They are one trajectory pass. A solvated replica is
23-30 GB on a ~100 MB/s share, so reading it is the entire budget: ~4 minutes a
replica, against which every distance search below is noise. Splitting the three
observables into three modules would have meant three passes and three times the
wall clock for nothing. `run()` walks the file once and writes three tidy
tables; `run_shell` / `run_bridges` / `run_residence` exist for a driver that
wants one of them, and each of those costs its own pass -- which is why the
driver should call `run()`.

WHAT IS MEASURED
----------------
An INTERFACIAL WATER at a frame is a water OXYGEN that lies within
D_WATER_SHELL = 3.5 A, BY MINIMUM IMAGE, of a heavy atom on BOTH sides of the
interface the spec declares:

    receptor side   spec.selections["receptor"], plus spec.selections["glycans"]
                    (see RECEPTOR SIDE below)
    ligand side     spec.selections["ligand"]

That is the full, unrestricted two-sided definition. It is COMPUTED, exactly, on
a reduced pair of atom sets, because an interfacial water's receptor partner
must lie within 3.5 + 3.5 = 7.0 A of some ligand atom and vice versa:

    rec_if = receptor heavy atoms within D_INTERFACE = 7.0 A of the ligand
    lig_if = ligand heavy atoms within D_INTERFACE = 7.0 A of the receptor

Restricting the two-sided search to rec_if / lig_if changes nothing -- no
interfacial water can be missed -- while cutting 23,064 x 5,376 heavy atoms
(1F4 OPEN, trimer + glycans against three Fvs) down to ~1,020 and ~1,130, a few
hundred per protomer. `n_rif_atoms` / `n_aif_atoms` are written out at EVERY
analysed frame, unthinned, so the reduction is auditable: they are the only
evidence a reader has that the thing that was searched was big enough.

TWO SEARCH MODES, and the output says which ran. The original's fast path -- one
bounding-box water window plus two small KD-trees -- needs the interface patch
to be small against L/2, which holds for one Fab against one protomer and need
not hold for a whole assembly, a filament or a membrane. When `patch_ok` fails,
the frame falls back to a box-aware `capped_distance` over every water oxygen:
slower, exact, and the same code path the original's own verification test used
as its INDEPENDENT recomputation. `n_frames_patch_search` and
`n_frames_exact_search` count them, so the mode is never a hidden assumption.

THE TWO CONTEXT COUNTS. `n_rec_only` and `n_lig_only` are the ONE-SIDED counts
on those same interface faces. They are the denominators that separate "this
interface is wet" from "these two proteins are of course solvated". Because
rec_if and lig_if are supersets of the atoms any interfacial water touches,

    n_interfacial <= min(n_rec_only, n_lig_only)            BY CONSTRUCTION

which this module asserts on every frame, as the original's test did. Note the
sense of the names: these are the ONE-SIDED counts, not the exclusive "near this
side and not the other" -- the exclusive counts are written too, as
`n_rec_excl` / `n_lig_excl`, and the three partition cleanly.

PBC -- THE TWO RULES THIS MODULE OBEYS
--------------------------------------
NAMD runs these systems with `wrapAll on` / `wrapWater on`, so water IS wrapped
and every shell search must be minimum-image. But `cellOrigin` is NOT the
origin -- 1F4's is (169.7, 169.0, 158.7) -- so the cell spans origin +- L/2 and
not [0, L). MDAnalysis's `dimensions` implicitly assumes [0, L), so
`AtomGroup.wrap()` / `pack_into_box()` would translate every atom by
-cellOrigin + L/2 and silently tear a solute that is intact.

NOTHING HERE IS EVER WRAPPED. Not an AtomGroup, not a Timestep. The only
periodic operation is an explicit minimum-image displacement applied to a COPY,

    d = x - c;  d -= L * rint(d / L)

which is `core.pbc.minimum_image`. (MDAnalysis's box-aware `capped_distance` /
`calc_angles` are used on the slow path and for the hydrogen-bond geometry: they
perform that same arithmetic internally on a copy of the coordinates they were
handed, and the Universe's own positions are untouched. The original's stage 71
and its verification test both relied on exactly that.)

TWO KINDS OF IMAGING, AND THEY ARE NOT INTERCHANGEABLE. `core.pbc` also offers
`shift_into_image`, which moves whole molecules -- that is the right tool for an
RMSD, a centre of mass or a render, where tearing a molecule in half is the
failure. A distance search needs the opposite: every ATOM placed in the image
nearest the patch being searched, molecule or no molecule. This module images
atoms about a patch centre, and never computes a molecular centre from the
result.

RECEPTOR SIDE, AND WHY GLYCANS ARE ON IT
----------------------------------------
The original's receptor side was "protomer P's gp120 + gp41 + protomer P's 25
glycans", because glycan burial IS part of this interface -- 1F4 buries 732 A^2
of glycan, 62% of its total. Dropping the glycans would not make the number
approximate, it would make it a different number.

Generically, the receptor side is the union of the spec selections named in
`spec.options["receptor_side"]`, default `("receptor", "glycans")`, and the
ligand side those in `spec.options["ligand_side"]`, default `("ligand",)`. A
system whose glycans sit on the LIGAND says so there. Atoms that land on both
sides are removed from the receptor side and counted (`n_side_overlap_atoms`):
an atom cannot be on both sides of an interface, and silently leaving it on both
would let a water bridge a side to itself.

The receptor-partner CLASSES (which part of the receptor a water or a bridge
touches) are likewise the spec's own named selections. The original detected
them with `np.char.startswith(seg, "G")` for glycans and `endswith("41")` for
gp41 -- segid-prefix guessing, which is exactly right for one project's naming
convention and silently wrong for every other.

TIME
----
There is deliberately NO module-level frame count and NO module-level ns/frame.
The projects differ -- this OPEN set is 0.2 ns/frame, OCCLUDED-OPEN and CLOSED
are 0.1 -- and a constant here is how a table ends up pooling them. Spacing comes
from `md.framestep_ns`, per replica, always, and every threshold in this file is
declared in NANOSECONDS and converted per replica. The original expressed
`FRAME_ROW_STRIDE`, `HYD_STRIDE`, `LONG_FRAMES` and `MAX_TAU_FRAMES` in frames,
which silently meant different durations the moment the project changed.

THE SAMPLING FLOOR, WHICH IS NOT OPTIONAL READING
-------------------------------------------------
Interfacial water exchanges on ps to ns. At the original's 0.1 ns/frame a large
majority of residences were already SINGLE FRAMES
(`residence_frac_runs_single_frame` near 0.9), which means 0.1 ns is at or below
the floor, not comfortably above it. At this project's 0.2 ns/frame, residence
times and bridge lifetimes are simply not resolved: a one-frame residence means
"present at this sample and not the next", which is anything from a fleeting
excursion to just under 0.4 ns.

So `MIN_SPACING_NS = 0.02` for the residence observable, and every row of the
residence table carries the verdict of `core.spec.check_spacing`: at 0.2 ns the
stage still RUNS -- the tail is real and is the thing worth comparing between
systems -- but the rows are marked `resolved=False` and carry a `floor_note`
saying what spacing would be needed. A plausible number with no warning is worse
than no number.

Shell counts and bridge counts are instantaneous, and a time average of an
instantaneous property is unbiased at any spacing, so those rows are resolved.

WHAT THIS PORT DOES NOT DO
--------------------------
The original's cross-protomer secondary interfaces (antibody copy P against a
non-cognate protomer Q) have no generic analogue: a spec declares ONE receptor
and ONE ligand, not a set of cognate pairs. A system that needs them declares
the second interface as its own spec.
"""

from core import threads                                   # noqa: F401
import collections

import numpy as np

from core import frames as F
from core import geom, io, pbc
from core import selections as S
from core import spec as SP

#: The stage family. Each observable writes its own file under its own name,
#: because `core.spec.BUNDLES` lists them separately and a reader asking for the
#: shell should not have to read a residence table to find it.
STAGE = "water"
SHELL = "water_shell"
BRIDGES = "water_bridges"
RESIDENCE = "water_residence"
STAGES = (SHELL, BRIDGES, RESIDENCE)

#: The coarsest spacing at which the finest observable here still means
#: something. See THE SAMPLING FLOOR above: at 0.1 ns/frame ~90% of residences
#: were single frames, so the floor is below 0.1, and resolving a distribution
#: whose mode is tens of ps needs ~0.02 ns.
MIN_SPACING_NS = 0.02
DEFAULT_STRIDE = 1

#: Per observable, because they do not share a floor. A count is a count at any
#: spacing; a lifetime is not.
MIN_SPACING_NS_FOR = {SHELL: None, BRIDGES: None, RESIDENCE: MIN_SPACING_NS}

# ---------------------------------------------------------------- geometry, A
D_WATER_SHELL = 3.5      # water O to a heavy atom, both sides
D_INTERFACE = 2.0 * D_WATER_SHELL   # 7.0, the EXACT reduction radius
D_HBOND = 3.5            # donor-acceptor heavy-atom distance, upper bound
A_HBOND = 150.0          # D-H...A angle, strictly greater than
D_H_COVALENT = 1.2       # donor-hydrogen pairing distance (HBA's default)
MIN_DA = 1.0             # HBA's own min_cutoff: an atom cannot bond to itself

#: A patch (interface plus its water window) must be small enough that a single
#: minimum-image displacement about one centre is unambiguous. half + pad must
#: stay below this fraction of L/2; measured worst case is ~0.6 for one Fab
#: against one protomer. When it is not met the frame falls back to a box-aware
#: neighbour search, which is slower and exact, and the two modes are counted.
PATCH_MAX_FRAC = 0.90

# ---------------------------------------------------------------- durations, ns
LONG_NS = 1.0            # "long-resident" threshold
MAX_TAU_NS = 100.0       # S(tau) is numerical noise long before this
HYDRATION_EVERY_NS = 5.0     # whole-side hydration reference
FRAME_ROW_EVERY_NS = 1.0     # thinning for the per-frame TSV rows
STABLE_FRAC = 0.70       # one uninterrupted episode covering this much of the window
FLICKER_MEAN_NS = 2.0    # mean episode shorter than this
HIST_EDGES_NS = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0,
                          20.0, 50.0, 100.0, 1e9])

N_BLOCKS = 5
PERSIST = 0.50           # occupancy threshold for a "persistent" bridged pair
OCC_THRESHOLDS = (0.05, 0.10, 0.25, 0.50)

#: Reading is the budget: a forward walk streams, a seek costs ~0.2 s on this
#: NFS mount. Walk forwards when most of what is read is used; seek when more
#: than this fraction of it would be thrown away. A pilot (`max_frames=10,
#: stride=100`) must cost seconds, and a forward walk to frame 900 of a 23 GB
#: file costs 90.
SEEK_WHEN_DISCARDING = 3

#: O(n) build, because the tree is rebuilt every frame and never reused.
KDT = dict(compact_nodes=False, balanced_tree=False)

#: Water oxygen by NAME, in preference order. CHARMM calls it OH2. `name O` is
#: tried only INSIDE the water group and only when it picks exactly one atom per
#: water residue: selected globally it returns backbone carbonyls instead and
#: yields a plausible, wrong answer, which is the mistake this list exists to
#: prevent.
WATER_O_NAMES = S.WATER_O_NAMES + ("O",)

#: Spec selections that may name a part of the receptor worth reporting
#: separately. Whichever of these the spec declares and which intersect the
#: receptor side become the `subset` of the per-class rows, and a spec can name
#: others through `options["receptor_classes"]`. This is the replacement for the
#: original's segid-prefix guessing; the classes may overlap, exactly as
#: "glycan" and "gp120" did, so they need not sum to the total.
CLASS_NAMES = ("receptor", "glycans", "epitope", "engaged")

DEFAULT_RECEPTOR_SIDE = ("receptor", "glycans")
DEFAULT_LIGAND_SIDE = ("ligand",)


# =========================================================================
# selections and topology
# =========================================================================
def _universe(topology, trajectories):
    import MDAnalysis as mda
    return mda.Universe(topology, *trajectories)


def _side(u, spec, names, label):
    """One side of the interface: the union of the named spec selections.

    Returns the FULL group and its heavy-atom subset, because the two answer
    different questions and conflating them is a silent failure: the shell is
    defined on HEAVY atoms by NAME (`not name H*`, never by mass -- these
    topologies have broken masses, and the lab's own contact script documents a
    mass test keeping 100 of 4,963 heavy atoms), while the hydrogen-bond
    geometry needs the hydrogens, and a donor set built from a heavy-only group
    is empty.
    """
    parts, used = [], []
    for name in names:
        text = spec.selections.get(name)
        if not str(text or "").strip():
            continue
        ag = u.select_atoms(S.translate(text))
        if not len(ag):
            continue
        parts.append(ag)
        used.append(name)
    if not parts:
        raise S.EmptySelection(
            f"{label}: none of {list(names)} is declared and non-empty in this "
            f"spec. A water shell has two sides; this one has none.")
    ag = parts[0]
    for extra in parts[1:]:
        ag = ag.union(extra)
    return ag, used


def _sides(u, spec):
    """Receptor and ligand groups, disjoint, with the overlap counted.

    An atom on both sides would let a water "bridge" one side to itself, and the
    bridge count would be quietly inflated by however many atoms the two
    selections share. So the overlap is removed from the receptor side -- the
    ligand is the smaller, more specific set -- and the number is written out.
    """
    rec_names = tuple(spec.options.get("receptor_side") or DEFAULT_RECEPTOR_SIDE)
    lig_names = tuple(spec.options.get("ligand_side") or DEFAULT_LIGAND_SIDE)
    lig, lig_used = _side(u, spec, lig_names, "ligand side")
    rec, rec_used = _side(u, spec, rec_names, "receptor side")
    overlap = np.intersect1d(rec.indices, lig.indices)
    if len(overlap):
        rec = rec[~np.isin(rec.indices, overlap)]
        if not len(rec):
            raise S.EmptySelection(
                "receptor side is entirely contained in the ligand side; these "
                "two selections describe one group, not an interface")
    rec_h = rec.select_atoms(S.HEAVY)
    lig_h = lig.select_atoms(S.HEAVY)
    if not len(rec_h) or not len(lig_h):
        raise S.EmptySelection(
            "one side has no heavy atoms; check that the selections name a "
            "solute and not, say, only hydrogens")
    return rec, lig, rec_h, lig_h, rec_used, lig_used, int(len(overlap))


def water_oxygens(u, water_sel):
    """The water oxygens, by NAME, with the count checked against the residues.

    The check is the point: a name that matches some-but-not-all waters gives a
    smaller shell and no error at all.
    """
    wat = u.select_atoms(S.translate(water_sel))
    if not len(wat):
        raise S.EmptySelection(
            f"water: {water_sel!r} matched no atoms. This stage needs the "
            f"SOLVATED trajectory.")
    n_res = wat.n_residues
    for nm in WATER_O_NAMES:
        og = wat.select_atoms(f"name {nm}")
        if len(og) == n_res:
            return og, nm
    raise ValueError(
        f"water: none of {WATER_O_NAMES} picks exactly one oxygen per water "
        f"residue ({n_res} residues; counts "
        f"{[len(wat.select_atoms('name ' + n)) for n in WATER_O_NAMES]}). Name "
        f"the oxygen explicitly rather than let this guess.")


def water_hydrogens(u, wo):
    """Global indices of each water's two hydrogens, in `wo`'s own order.

    TIP3 is written OH2, H1, H2, so the hydrogens are the next two atom indices
    and no second selection is needed. That is a shortcut in COST only, and it
    is checked rather than assumed -- the fallback asks the residues.
    """
    names = u.atoms.names.astype(str)
    idx = np.asarray(wo.indices, dtype=np.int64)
    if int(idx.max()) + 2 < u.atoms.n_atoms:
        h1, h2 = idx + 1, idx + 2
        if (np.char.startswith(names[h1], "H").all()
                and np.char.startswith(names[h2], "H").all()):
            return h1, h2
    hyd = wo.residues.atoms.select_atoms("name H*")
    if len(hyd) != 2 * len(idx) or not np.array_equal(
            np.asarray(hyd.resindices), np.repeat(wo.resindices, 2)):
        raise ValueError(
            "water hydrogens are not two per oxygen in residue order; the "
            "bridge search cannot pair them. Check the water model.")
    hh = np.asarray(hyd.indices, dtype=np.int64).reshape(len(idx), 2)
    return hh[:, 0], hh[:, 1]


def water_block(wo):
    """(offset, count, contiguous?) for the stride-3 water-oxygen slice.

    In a CHARMM-built solvated system the waters are one uniform OH2/H1/H2 block
    after the solute, so `pos[W0:W0+3*NW:3]` is a basic slice -- a VIEW, no copy
    -- and gathering 236,838 oxygen positions costs nothing. Fancy indexing the
    same atoms copies 2.8 MB per frame, which at 2,520 frames is 7 GB of pointless
    memory traffic. The optimisation is a real win and it holds on these builds.

    It is TESTED, not asserted: a build that interleaves ions with water, or
    orders a water H before its O, is unusual but not wrong, and the right
    response is an index array, not a crash.
    """
    idx = np.asarray(wo.indices, dtype=np.int64)
    W0, NW = int(idx[0]), len(idx)
    contiguous = bool(np.array_equal(idx, np.arange(W0, W0 + 3 * NW, 3)))
    return W0, NW, contiguous


def _class_masks(u, spec, rec_h):
    """Per-receptor-atom membership in each named sub-selection the spec declares.

    This is what replaces `np.char.startswith(seg, "G")` (glycan) and
    `np.char.endswith(seg, "41")` (gp41). Those work perfectly on one project's
    segid convention and return an empty mask, with no error, on any other --
    the per-class rows would all read zero and look like a finding.
    """
    out = {}
    names = tuple(spec.options.get("receptor_classes") or CLASS_NAMES)
    for name in names:
        text = spec.selections.get(name)
        if not str(text or "").strip():
            continue
        try:
            ag = u.select_atoms(S.translate(text))
        except Exception:                                  # noqa: BLE001
            continue                                       # not this topology's
        if not len(ag):
            continue
        m = np.isin(rec_h.indices, ag.indices)
        if m.any():
            out[name] = m
    return out


# =========================================================================
# geometry kernels -- explicit minimum image, on a copy, always
# =========================================================================
def patch_ok(half, pad, L):
    """Is a single-centre minimum image unambiguous for a patch this big?"""
    if L is None:
        return True
    return bool((np.asarray(half) + pad < PATCH_MAX_FRAC * 0.5 * np.asarray(L)).all())


def image_about(X, c, L):
    """Copy of X placed in the image nearest to c. X is never modified."""
    if L is None:
        return np.asarray(X, dtype=np.float64)
    return pbc.minimum_image(np.asarray(X, dtype=np.float64) - c, L) + c


def coherent(X, L):
    """(imaged coordinates, centre, half-extent) for one group, safely.

    Two passes on purpose. A group split across the boundary has a MEANINGLESS
    centroid -- it sits in the middle of the box, between the two halves -- so
    imaging about it would leave the group split. The first pass images about
    one arbitrary ATOM, which is always in some real image and therefore always
    a valid reference; the centroid of that result is meaningful, and the second
    pass re-centres on it to make the half-extent as small as it can be. On a
    trajectory that was never wrapped the second pass is a no-op.
    """
    X = np.asarray(X, dtype=np.float64)
    if L is None or not len(X):
        c = X.mean(0) if len(X) else np.zeros(3)
        return X, c, (np.abs(X - c).max(0) if len(X) else np.zeros(3))
    Y = image_about(X, X[0], L)
    c = Y.mean(0)
    Y = image_about(Y, c, L)
    c = Y.mean(0)
    return Y, c, np.abs(Y - c).max(0)


def near_atoms(X, cen, half, pad, L):
    """Indices into X of atoms within half+pad of `cen` on every axis, and their
    coordinates imaged about `cen`.

    A bounding-box prefilter, exact when half+pad < L/2 (which `patch_ok`
    guarantees before this is called): an atom further than half+pad from the
    centre on any axis cannot be within `pad` of any atom of the patch.
    """
    d = pbc.minimum_image(np.asarray(X, dtype=np.float64) - cen, L) \
        if L is not None else np.asarray(X, dtype=np.float64) - cen
    k = (np.abs(d) <= (half + pad)).all(1)
    idx = np.flatnonzero(k)
    return idx, np.ascontiguousarray(d[idx] + cen)


def water_window(wx, wxyz, cen, half, pad, L):
    """Water oxygens within half+pad of `cen` on every axis, imaged about `cen`.

    Exactness: a water within `pad` of a patch atom is within half+pad of `cen`
    on every axis, provided half+pad < L/2 (asserted by the caller), because the
    two displacements then add without folding. So this is a superset of every
    water the fine search can find -- a prefilter, never a decision.

    The x axis is filtered first on a contiguous copy of the x column, and y/z
    only on the survivors; that is what keeps an O(236k) pass off the critical
    path three times per frame.
    """
    R = np.asarray(half) + pad
    d = np.asarray(wx, dtype=np.float64) - cen[0]
    if L is not None:
        d = pbc.minimum_image(d, L[0])
    s = np.flatnonzero(np.abs(d) <= R[0])
    sub = np.asarray(wxyz[s], dtype=np.float64)
    for ax in (1, 2):
        d = sub[:, ax] - cen[ax]
        if L is not None:
            d = pbc.minimum_image(d, L[ax])
        k = np.abs(d) <= R[ax]
        s = s[k]
        sub = sub[k]
    return s, np.ascontiguousarray(image_about(sub, cen, L))


def within(tree, X, cut):
    """Boolean mask: rows of X within `cut` of any point in `tree`."""
    d, _ = tree.query(X, k=1, distance_upper_bound=cut)
    return d <= cut


# =========================================================================
# hydrogen-bond machinery (the bridge definition)
# =========================================================================
def _map_sel(m):
    """{resname -> names} as a selection string."""
    return " or ".join(f"(resname {k} and name {' '.join(v)})"
                       for k, v in sorted(m.items()))


def polar_maps(ag, box):
    """{resname: names} for N/O donors bearing a polar H, and for those hydrogens.

    Derived from GEOMETRY rather than a hand-written table because the table is
    always wrong somewhere: BGLC carries a free O4-HO4 only in some linkage
    contexts, and ANE5 puts HT1/HT2/HT3 on a methyl CARBON, which any name-only
    rule would wrongly call polar.

    The 1.2 A rule is not an approximation to the PSF's bonds; it IS the bonds --
    the original's test asserts set equality between the two over every solute
    it has. Using geometry keeps this stage's donor set identical to the direct
    hydrogen-bond stage's, so the water-mediated and direct counts are
    comparable cell for cell. That comparison is the whole reason for the stage.
    """
    from MDAnalysis.lib.distances import capped_distance
    NO = ag.select_atoms("name N* or name O*")
    H = ag.select_atoms("name H*")
    if not len(NO) or not len(H):
        raise ValueError("no N/O or no H atoms on this side; a water bridge "
                         "needs both, so the topology or the selection is wrong")
    pairs = capped_distance(NO.positions, H.positions, max_cutoff=D_H_COVALENT,
                            box=box, return_distances=False)
    if not len(pairs):
        raise ValueError(f"no donor-hydrogen pairs within {D_H_COVALENT} A")
    dmap = collections.defaultdict(set)
    hmap = collections.defaultdict(set)
    dn, dr = NO.names.astype(str), NO.resnames.astype(str)
    hn, hr = H.names.astype(str), H.resnames.astype(str)
    for i, j in pairs:
        dmap[dr[i]].add(dn[i])
        hmap[hr[j]].add(hn[j])
    return ({k: sorted(v) for k, v in dmap.items()},
            {k: sorted(v) for k, v in hmap.items()})


class SidePolar:
    """One side's donor / hydrogen / acceptor atom indices, in HBA's own layout.

    `don` and `hyd` are PARALLEL, one entry per D-H bond, so a donor carrying two
    hydrogens (ARG NH1, LYS NZ, water OH2) appears twice and gets two independent
    chances -- exactly what MDAnalysis's HydrogenBondAnalysis._get_dh_pairs
    returns, and the reason the counts here line up with a direct H-bond stage.
    """

    __slots__ = ("don", "hyd", "acc", "label")

    def __init__(self, ag, dmap, hmap, box, label):
        from MDAnalysis.lib.distances import capped_distance
        don_ag = ag.select_atoms(_map_sel(dmap))
        hyd_ag = ag.select_atoms(_map_sel(hmap))
        acc_ag = ag.select_atoms("name N* or name O*")
        if not len(acc_ag):
            raise ValueError(f"{label}: no acceptor (N/O) atoms")
        p = capped_distance(don_ag.positions, hyd_ag.positions,
                            max_cutoff=D_H_COVALENT, box=box,
                            return_distances=False)
        if not len(p):
            raise ValueError(f"{label}: no donor-hydrogen pairs")
        self.don = don_ag.indices[p[:, 0]].astype(np.int64)
        self.hyd = hyd_ag.indices[p[:, 1]].astype(np.int64)
        self.acc = acc_ag.indices.astype(np.int64)
        self.label = label


def side_hbonds(pos, box, Ow, wH1, wH2, side):
    """Every hydrogen bond between a set of waters and one side of the interface.

    Water is both a donor (via H1/H2) and an acceptor (via its oxygen); both
    directions are searched on both sides, so all four combinations of
    (water donates / water accepts) x (receptor / ligand) are covered.

    The criteria are HBA's, exactly as the direct hydrogen-bond stage configures
    them: donor-acceptor heavy-atom distance in the half-open interval
    (1.0, 3.5] A -- the 1.0 A floor is HBA's own `min_cutoff`, and the interval
    is (min, max], not [min, max] -- and a D-H...A angle STRICTLY GREATER than
    150 deg, because HBA uses `angles > cut`, not `>=`. At thousands of frames
    the difference is immaterial; nothing here is left to chance.

    Angles come from `calc_angles` WITH the box, i.e. from minimum-image
    displacement vectors, never from raw coordinate differences.

    Returns (water_slot, side_atom_index), deduplicated: a water bonded to the
    same side atom through both of its hydrogens counts once.
    """
    from MDAnalysis.lib.distances import capped_distance, calc_angles
    if not len(Ow):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    ws, as_ = [], []

    # --- water donates: D = water O, H = H1/H2, A = side N/O ---------------
    p = capped_distance(Ow, pos[side.acc], max_cutoff=D_HBOND,
                        min_cutoff=MIN_DA, box=box, return_distances=False)
    if len(p):
        wi, ai = p[:, 0], p[:, 1]
        A = pos[side.acc[ai]]
        D = Ow[wi]
        for wH in (wH1, wH2):
            ang = np.rad2deg(calc_angles(D, pos[wH[wi]], A, box=box))
            m = ang > A_HBOND
            if m.any():
                ws.append(wi[m])
                as_.append(side.acc[ai[m]])

    # --- side donates: D = side N/O, H = its H, A = water O ----------------
    p = capped_distance(pos[side.don], Ow, max_cutoff=D_HBOND,
                        min_cutoff=MIN_DA, box=box, return_distances=False)
    if len(p):
        di, wi = p[:, 0], p[:, 1]
        ang = np.rad2deg(calc_angles(pos[side.don[di]], pos[side.hyd[di]],
                                     Ow[wi], box=box))
        m = ang > A_HBOND
        if m.any():
            ws.append(wi[m])
            as_.append(side.don[di[m]])

    if not ws:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    w = np.concatenate(ws).astype(np.int64)
    a = np.concatenate(as_).astype(np.int64)
    uw = np.unique(np.stack([w, a], axis=1), axis=0)
    return uw[:, 0], uw[:, 1]


def check_water_geometry(pos, box, oidx, h1, h2):
    """Assert the structural O-H pairing really is the covalent one."""
    from MDAnalysis.lib.distances import calc_bonds
    for h in (h1, h2):
        d = calc_bonds(pos[oidx], pos[h], box=box)
        if d.size and float(d.max()) > D_H_COVALENT:
            raise ValueError(
                f"water O-H distance up to {d.max():.3f} A exceeds the "
                f"{D_H_COVALENT} A pairing cutoff: the oxygens and hydrogens "
                f"this stage paired are not the same molecule")


# =========================================================================
# the per-frame shell search
# =========================================================================
def _shell_patch(wx, wxyz, RIF, AIF, L):
    """The fast path: one water window, two small KD-trees.

    Valid only when the interface patch is small against L/2, which the caller
    has checked. RIF and AIF are already imaged into one coherent patch, so
    plain Euclidean distances inside it ARE the minimum-image distances.
    """
    from scipy.spatial import cKDTree
    IF = np.vstack([RIF, AIF])
    cif = IF.mean(0)
    half = np.abs(IF - cif).max(0)
    if not patch_ok(half, D_WATER_SHELL, L):
        return None
    s, Wc = water_window(wx, wxyz, cif, half, D_WATER_SHELL, L)
    if not len(s):
        return s, np.zeros(0, bool), np.zeros(0, bool), Wc
    mr = within(cKDTree(RIF, **KDT), Wc, D_WATER_SHELL)
    ma = within(cKDTree(AIF, **KDT), Wc, D_WATER_SHELL)
    return s, mr, ma, Wc


def _shell_exact(W, RIF_raw, AIF_raw, box):
    """The slow path: a box-aware neighbour search over every water oxygen.

    Used when the interface is too spread out for a single-centre minimum image
    to be unambiguous -- three Fvs around a trimer, a filament, a membrane. It
    is not an approximation of the fast path, it is the same answer computed
    without the bounding-box prefilter, and it is what the original's
    verification test used as its INDEPENDENT recomputation.
    """
    from MDAnalysis.lib.distances import capped_distance
    pr = capped_distance(W, RIF_raw, max_cutoff=D_WATER_SHELL, box=box,
                         return_distances=False)
    pa = capped_distance(W, AIF_raw, max_cutoff=D_WATER_SHELL, box=box,
                         return_distances=False)
    Sr = np.unique(pr[:, 0]) if len(pr) else np.zeros(0, np.int64)
    Sa = np.unique(pa[:, 0]) if len(pa) else np.zeros(0, np.int64)
    s = np.union1d(Sr, Sa).astype(np.int64)
    return s, np.isin(s, Sr), np.isin(s, Sa)


def _pairs_patch(Wi, XIF):
    """(water slot, side atom) pairs within the shell, inside one imaged patch."""
    if not len(Wi) or not len(XIF):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    dd = ((Wi[:, None, :] - XIF[None, :, :]) ** 2).sum(-1)
    iw, ia = np.nonzero(dd <= D_WATER_SHELL * D_WATER_SHELL)
    return iw.astype(np.int64), ia.astype(np.int64)


def _pairs_exact(Wpos, XIF_raw, box):
    from MDAnalysis.lib.distances import capped_distance
    if not len(Wpos) or not len(XIF_raw):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    p = capped_distance(Wpos, XIF_raw, max_cutoff=D_WATER_SHELL, box=box,
                        return_distances=False)
    if not len(p):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return p[:, 0].astype(np.int64), p[:, 1].astype(np.int64)


def _hydration(wx, wxyz, X, cen, half, L):
    """First-shell water count of one already-imaged, compact atom set.

    -1 means "the group is too big for single-centre imaging to be safe", which
    is information, not a failure: a ratio built from it would be wrong, and a
    -1 is filtered out downstream where a 0 would not be.
    """
    from scipy.spatial import cKDTree
    if not patch_ok(half, D_WATER_SHELL, L):
        return -1
    s, Wc = water_window(wx, wxyz, cen, half, D_WATER_SHELL, L)
    if not len(s):
        return 0
    return int(within(cKDTree(X, **KDT), Wc, D_WATER_SHELL).sum())


# =========================================================================
# the accumulators
# =========================================================================
class Scan:
    """Everything one replica's single pass produced."""

    def __init__(self, md, stride, n, n_rec_res, n_lig_res):
        self.md_id = md.id
        self.framestep_ns = md.framestep_ns
        self.stride = stride
        self.eff_ns = float(md.framestep_ns or 0.0) * stride
        self.n = n
        z = lambda: np.zeros(n, np.int32)                   # noqa: E731
        self.n_if, self.n_rec, self.n_lig = z(), z(), z()
        self.n_rif, self.n_aif, self.n_cand = z(), z(), z()
        self.time_ns = np.zeros(n, np.float64)
        self.raw_frame = np.zeros(n, np.int64)
        self.occ_f, self.occ_w = [], []
        self.rec_wet = np.zeros(n_rec_res, np.float64)
        self.lig_wet = np.zeros(n_lig_res, np.float64)
        self.rec_ever = np.zeros(n_rec_res, bool)
        self.lig_ever = np.zeros(n_lig_res, bool)
        self.rec_if = np.zeros(n_rec_res, np.float64)
        self.lig_if = np.zeros(n_lig_res, np.float64)
        self.hyd_k, self.hyd_rec, self.hyd_lig = [], [], []
        self.n_patch = 0
        self.n_exact = 0
        self.n_overlap = 0
        self.n_prod = 0
        self.n_analysed = 0
        self.bridges = None
        self.occ = None
        self.water_ordinal = None
        self.water_resid = None
        self.water_segid = None
        self.water_atomindex = None
        self.notes = []


class Bridges:
    """Per-frame bridge counts and the counters the occupancies come from.

    Counts only -- never per-frame Python objects -- so memory is O(distinct
    pairs), not O(frames x pairs).
    """

    def __init__(self, n, class_names):
        z = lambda: np.zeros(n, np.int32)                   # noqa: E731
        self.n_water = z()          # distinct bridging waters
        self.n_pairs = z()          # distinct (receptor residue, ligand residue)
        self.n_triples = z()        # distinct (rec atom, water, lig atom)
        self.n_flagged = z()        # interfacial waters offered to the search
        self.by_class = {c: z() for c in class_names}
        self.pair_frames = collections.Counter()
        self.pair_waters = collections.Counter()
        self.rec_frames = collections.Counter()
        self.rec_waters = collections.Counter()
        self.lig_frames = collections.Counter()
        self.lig_waters = collections.Counter()
        self.wat_bridge = collections.defaultdict(list)


# =========================================================================
# one replica, one pass
# =========================================================================
def _frame_plan(md, u, stride, max_frames):
    """Which raw frames to analyse, what time each one is, and stream or seek.

    SEAM FRAMES ARE DROPPED. A restart re-emits a step with different
    coordinates -- the window resumed from an earlier checkpoint and
    re-simulated up to it -- and core/frames.py finds those raw indices from the
    DCD headers, never from the coordinates. On 1F4 that is 5 frames of 2,525.

    TIME COMES FROM THE HEADERS, not from the raw index. `raw_index x
    framestep_ns` counts a dropped seam as advancing the clock, which it does
    not: it would put every frame after the last seam 1 ns late on a 504 ns run,
    and the error grows with the number of restarts. `Timeline.time_ns` reads the
    step counter instead.
    """
    timeline = F.build(md.id, list(md.timeline_inputs()),
                       framestep_ns=md.framestep_ns)
    n_raw = len(u.trajectory)
    seams = set(timeline.seams)
    if n_raw != len(timeline) + len(seams):
        # A merged file may already have been deduplicated; say so rather than
        # dropping frames that are not there.
        seams = set()
    kept = np.array([i for i in range(n_raw) if i not in seams], dtype=np.int64)
    take = np.arange(len(kept))[::max(1, stride)]
    if max_frames:
        take = take[:int(max_frames)]
    want = kept[take]
    if len(kept) == len(timeline):
        times = np.array([timeline.time_ns(int(j)) for j in take], np.float64)
    else:
        times = take.astype(np.float64) * float(md.framestep_ns or 0.0)
    seek = bool(len(want) * SEEK_WHEN_DISCARDING < n_raw)
    return timeline, want, times, seek, len(seams)


def _walk(u, want, seek):
    """(position in `want`, Timestep) for every frame to analyse.

    Streaming is the production path: a 23 GB file is not something to seek
    around in. Seeking is for a capped or heavily strided run, where a forward
    walk would read gigabytes only to discard them.
    """
    if seek:
        for k, ts in enumerate(u.trajectory[want.tolist()]):
            yield k, ts
        return
    pos = {int(f): k for k, f in enumerate(want)}
    last = int(want[-1])
    for ts in u.trajectory:
        k = pos.get(ts.frame)
        if k is not None:
            yield k, ts
        if ts.frame >= last:
            break


def scan_md(spec, md, stride, max_frames=None, do_bridges=True, progress=0):
    """One replica: one forward pass, shell and bridges together."""
    top, trajs = md.inputs(solvated=True)
    u = _universe(top, trajs)

    rec_all, lig_all, rec, lig, rec_used, lig_used, n_overlap = _sides(u, spec)
    wo, o_name = water_oxygens(u, spec.selection("water"))
    W0, NW, contiguous = water_block(wo)
    h1, h2 = water_hydrogens(u, wo)
    w_index = np.asarray(wo.indices, dtype=np.int64)

    rec_idx = np.asarray(rec.indices, dtype=np.int64)
    lig_idx = np.asarray(lig.indices, dtype=np.int64)

    # Per-residue keying is on RESINDEX, never on (segid, resid): a topology
    # may reuse a resid across segments, and a shell that keys on the pair
    # quietly merges two residues into one row.
    rec_res = np.asarray(rec.residues.resindices, dtype=np.int64)
    lig_res = np.asarray(lig.residues.resindices, dtype=np.int64)
    nres = int(u.atoms.n_residues)
    m = np.full(nres, -1, np.int64)
    m[rec_res] = np.arange(len(rec_res))
    rec_row = m[rec.resindices]          # per receptor heavy ATOM -> residue row
    m = np.full(nres, -1, np.int64)
    m[lig_res] = np.arange(len(lig_res))
    lig_row = m[lig.resindices]

    classes = _class_masks(u, spec, rec)
    # global atom index -> row in the receptor heavy side, for the class
    # breakdown. Donors and acceptors are N/O, so a bridge's receptor partner is
    # always in this set.
    rec_pos = np.full(u.atoms.n_atoms, -1, np.int64)
    rec_pos[rec_idx] = np.arange(len(rec_idx))

    timeline, want, want_ns, seek, n_seams = _frame_plan(md, u, stride,
                                                        max_frames)
    n = len(want)
    if not n:
        raise ValueError(f"{md.id}: no frames to analyse")

    sc = Scan(md, stride, n, len(rec_res), len(lig_res))
    sc.n_analysed = n
    sc.n_overlap = n_overlap
    sc.notes = [f"receptor side = {'+'.join(rec_used)} ({len(rec_idx)} heavy "
                f"of {len(rec_all)})",
                f"ligand side = {'+'.join(lig_used)} ({len(lig_idx)} heavy "
                f"of {len(lig_all)})",
                f"water oxygen = name {o_name} ({NW})",
                f"water block contiguous = {contiguous}",
                f"seam frames dropped = {n_seams}",
                f"side overlap atoms removed = {n_overlap}",
                f"frame access = {'seek' if seek else 'stream'}",
                f"analysed {n} of {len(timeline)} deduplicated frames, "
                f"{want_ns[0]:g}-{want_ns[-1]:g} ns"]
    # A hole in the time axis means a segment was never copied or the run
    # restarted from an earlier checkpoint. An average across it is an average
    # across time that was never simulated, so it is stated, not hidden.
    for a, b, ns in F.gaps(timeline):
        sc.notes.append(f"GAP {ns} ns in the time axis between steps {a} "
                        f"and {b}")

    eq_ns = float(spec.options.get("equilibration_ns") or 0.0)
    hyd_every = max(1, int(round(HYDRATION_EVERY_NS / sc.eff_ns))) \
        if sc.eff_ns else 1

    sides = None
    if do_bridges:
        sc.bridges = Bridges(n, tuple(classes))

    for k, ts in _walk(u, want, seek):
        if progress and k and k % progress == 0:
            print(f"[{STAGE}] {md.id} frame {k}/{n}", flush=True)
        pos = ts.positions                      # a view, not a copy
        box = ts.dimensions
        L = pbc.box_lengths(box)                # raises on a triclinic cell
        sc.time_ns[k] = float(want_ns[k])
        sc.raw_frame[k] = int(want[k])
        prod = sc.time_ns[k] >= eq_ns
        if prod:
            sc.n_prod += 1

        # Water oxygen positions. The contiguous case is a basic slice and
        # therefore a view; the fallback copies, which is correct and slower.
        wxyz = pos[W0:W0 + 3 * NW:3] if contiguous else pos[w_index]
        wx = np.ascontiguousarray(wxyz[:, 0])

        # ---- the ligand side, made coherent ------------------------------
        A, cA, halfA = coherent(pos[lig_idx], L)

        # ---- the whole-side hydration reference, on its own schedule -----
        # BEFORE the interface work, and on the schedule regardless of what
        # the interface is doing: this is the bulk DENOMINATOR, and a frame
        # where the two sides happen not to touch is exactly a frame whose
        # hydration a ratio needs. Sampling it only when the interface is wet
        # would bias it upward. It is also the most expensive thing in the
        # loop -- a first-shell count over a whole side is a window over every
        # water -- which is why it runs every HYDRATION_EVERY_NS and not every
        # frame.
        if k % hyd_every == 0:
            Rc, cRc, halfRc = coherent(pos[rec_idx], L)
            sc.hyd_k.append(k)
            sc.hyd_lig.append(_hydration(wx, wxyz, A, cA, halfA, L))
            sc.hyd_rec.append(_hydration(wx, wxyz, Rc, cRc, halfRc, L))

        # ---- the provably-lossless reduction to interface atoms ----------
        # Any water within 3.5 A of BOTH sides has a receptor partner within
        # 3.5 + 3.5 = 7.0 A of a ligand atom, and vice versa. Restricting the
        # two-sided search to these faces cannot miss a single water.
        # `rif_loc` / `aif_loc` are positions WITHIN the side, not global atom
        # indices, so every downstream lookup (residue row, class mask) is a
        # plain gather. The original mixed the two spaces and needed a
        # searchsorted to get back.
        patch = patch_ok(halfA, D_INTERFACE, L)
        if patch:
            from scipy.spatial import cKDTree
            rk, Rn = near_atoms(pos[rec_idx], cA, halfA, D_INTERFACE, L)
            if not len(rk):
                continue
            tA = cKDTree(A, **KDT)
            rif = np.flatnonzero(within(tA, Rn, D_INTERFACE))
            if not len(rif):
                continue
            rif_loc = rk[rif]
            RIF = np.ascontiguousarray(Rn[rif])
            aif_loc = np.flatnonzero(within(cKDTree(RIF, **KDT), A, D_INTERFACE))
            AIF = np.ascontiguousarray(A[aif_loc])
        else:
            from MDAnalysis.lib.distances import capped_distance
            p = capped_distance(np.ascontiguousarray(pos[rec_idx]),
                                np.ascontiguousarray(pos[lig_idx]),
                                max_cutoff=D_INTERFACE, box=box,
                                return_distances=False)
            if not len(p):
                continue
            rif_loc = np.unique(p[:, 0])
            aif_loc = np.unique(p[:, 1])
            RIF = np.ascontiguousarray(pos[rec_idx[rif_loc]])
            AIF = np.ascontiguousarray(pos[lig_idx[aif_loc]])

        sc.n_rif[k] = len(rif_loc)
        sc.n_aif[k] = len(aif_loc)
        if prod:
            np.add.at(sc.rec_if, np.unique(rec_row[rif_loc]), 1.0)
            np.add.at(sc.lig_if, np.unique(lig_row[aif_loc]), 1.0)

        # ---- the two one-sided criteria, then their intersection ---------
        got = _shell_patch(wx, wxyz, RIF, AIF, L) if patch else None
        if got is None:
            sc.n_exact += 1
            s, mr, ma = _shell_exact(np.ascontiguousarray(wxyz), RIF, AIF, box)
            Wc = None
        else:
            sc.n_patch += 1
            s, mr, ma, Wc = got
        sc.n_cand[k] = len(s)
        if not len(s):
            continue
        mi = mr & ma
        sc.n_rec[k] = int(mr.sum())
        sc.n_lig[k] = int(ma.sum())
        ni = int(mi.sum())
        sc.n_if[k] = ni

        if ni == 0:
            continue
        wsel = np.flatnonzero(mi)
        sc.occ_f.append(np.full(ni, k, np.int32))
        sc.occ_w.append(s[wsel].astype(np.int32))

        # ---- per-residue wetness, from the interfacial waters only -------
        if prod:
            if Wc is not None:
                Wi = Wc[wsel]
                iw_r, ia_r = _pairs_patch(Wi, RIF)
                iw_a, ia_a = _pairs_patch(Wi, AIF)
            else:
                Wp = np.ascontiguousarray(wxyz[s[wsel]])
                iw_r, ia_r = _pairs_exact(Wp, RIF, box)
                iw_a, ia_a = _pairs_exact(Wp, AIF, box)
            for side_loc, rows_of, ia, iw, wet, ever in (
                    (rif_loc, rec_row, ia_r, iw_r, sc.rec_wet, sc.rec_ever),
                    (aif_loc, lig_row, ia_a, iw_a, sc.lig_wet, sc.lig_ever)):
                if not len(ia):
                    continue
                rows = rows_of[side_loc[ia]]
                # One count per (residue, water) pair, not per atom pair: a
                # water touching four atoms of one residue is ONE water on that
                # residue. The `rows * ni + iw` key deduplicates the pair in a
                # single `unique`, and `// ni` recovers the residue row.
                np.add.at(wet, np.unique(rows * ni + iw) // ni, 1.0)
                ever[np.unique(rows)] = True

        # ---- bridges, on the interfacial waters of this frame ------------
        if do_bridges:
            ordinals = s[wsel].astype(np.int64)
            o_idx = w_index[ordinals]
            if sides is None:
                # The polar sets are built ONCE, from the full groups (the
                # hydrogens are the point), at the first frame that has an
                # interfacial water -- they are a property of the topology, not
                # of the frame.
                dmap, hmap = polar_maps(rec_all.union(lig_all), box)
                sides = {"rec": SidePolar(rec_all, dmap, hmap, box, "receptor"),
                         "lig": SidePolar(lig_all, dmap, hmap, box, "ligand")}
                check_water_geometry(pos, box, o_idx[:512],
                                     h1[ordinals][:512], h2[ordinals][:512])
            _bridge_frame(sc.bridges, k, pos, box, o_idx, h1[ordinals],
                          h2[ordinals], sides, rec_pos, classes,
                          u.atoms.resindices, prod)

    _finish(sc, u, n)
    return sc, rec, lig, classes


def _bridge_frame(B, k, pos, box, o_idx, wH1, wH2, sides, rec_pos, classes,
                  resindex, prod):
    """One frame of the bridge search.

    A BRIDGE is one water simultaneously hydrogen-bonded to a receptor-side
    donor/acceptor and to a ligand-side one. `n_water` -- the number of DISTINCT
    BRIDGING WATERS -- is the primary observable, because a water is the physical
    object a bridge is made of. `n_pairs` (distinct residue pairs so connected)
    and `n_triples` (distinct rec-atom/water/lig-atom triples, the closest
    analogue to a direct H-bond record count) are reported beside it and are NOT
    interchangeable.
    """
    B.n_flagged[k] = len(o_idx)
    if not len(o_idx):
        return
    Ow = pos[o_idx]
    rw, ra = side_hbonds(pos, box, Ow, wH1, wH2, sides["rec"])
    aw, aa = side_hbonds(pos, box, Ow, wH1, wH2, sides["lig"])
    if not len(rw) or not len(aw):
        return
    both = np.intersect1d(np.unique(rw), np.unique(aw))
    nb = len(both)
    B.n_water[k] = nb
    if nb == 0:
        return
    for g in o_idx[both]:
        B.wat_bridge[int(g)].append(k)

    kr, ka = np.isin(rw, both), np.isin(aw, both)
    rw, ra = rw[kr], ra[kr]
    aw, aa = aw[ka], aa[ka]
    ri = np.searchsorted(both, rw)              # water slot within `both`
    ai = np.searchsorted(both, aw)

    n_rec_at = np.bincount(ri, minlength=nb).astype(np.int64)
    n_lig_at = np.bincount(ai, minlength=nb).astype(np.int64)
    B.n_triples[k] = int((n_rec_at * n_lig_at).sum())

    # Receptor-partner classes, from the spec's selections. A water touching
    # both a glycan and the protein is counted in BOTH, so the classes need not
    # sum to n_water -- which is the honest reading of "this water touches a
    # glycan", and is what the original did.
    rows = rec_pos[ra]
    for name, mask in classes.items():
        ok = rows >= 0
        if not ok.any():
            continue
        hit = np.zeros(len(rows), bool)
        hit[ok] = mask[rows[ok]]
        B.by_class[name][k] = int(np.unique(ri[hit]).size)

    R = np.unique(np.stack([ri, resindex[ra]], axis=1), axis=0)
    Aq = np.unique(np.stack([ai, resindex[aa]], axis=1), axis=0)
    grid = np.arange(nb)
    r0 = np.searchsorted(R[:, 0], grid, "left")
    r1 = np.searchsorted(R[:, 0], grid, "right")
    a0 = np.searchsorted(Aq[:, 0], grid, "left")
    a1 = np.searchsorted(Aq[:, 0], grid, "right")

    pairs_here = set()
    for w in range(nb):
        rs = R[r0[w]:r1[w], 1]
        as_ = Aq[a0[w]:a1[w], 1]
        if prod:
            for y in as_:
                B.lig_waters[int(y)] += 1
        for x in rs:
            xi = int(x)
            if prod:
                B.rec_waters[xi] += 1
            for y in as_:
                key = (xi, int(y))
                pairs_here.add(key)
                if prod:
                    B.pair_waters[key] += 1
    B.n_pairs[k] = len(pairs_here)
    if prod:
        for key in pairs_here:
            B.pair_frames[key] += 1
        for x in np.unique(R[:, 1]):
            B.rec_frames[int(x)] += 1
        for y in np.unique(Aq[:, 1]):
            B.lig_frames[int(y)] += 1


def _finish(sc, u, n):
    """Assemble the occupancy matrix and normalise the per-residue counters."""
    if sc.occ_f:
        ff = np.concatenate(sc.occ_f)
        ww = np.concatenate(sc.occ_w)
    else:
        ff = np.zeros(0, np.int32)
        ww = np.zeros(0, np.int32)
    uw, inv = np.unique(ww, return_inverse=True)
    occ = np.zeros((len(uw), n), np.uint8)
    if len(uw):
        occ[inv, ff] = 1
    sc.occ = occ
    sc.water_ordinal = uw.astype(np.int64)
    sc.occ_f = sc.occ_w = None
    denom = max(1, sc.n_prod)
    for a in (sc.rec_wet, sc.lig_wet, sc.rec_if, sc.lig_if):
        a /= denom

    # THE construction invariant, asserted on every frame, as the original's
    # test did. rec_if and lig_if are supersets of the atoms any interfacial
    # water touches, so the two-sided count cannot exceed either one-sided one.
    # If it does, the shell construction is wrong and nothing below is worth
    # writing.
    bad = int((sc.n_if > np.minimum(sc.n_rec, sc.n_lig)).sum())
    if bad:
        raise AssertionError(
            f"{sc.md_id}: n_interfacial exceeds a one-sided count at {bad} "
            f"frame(s). By construction that is impossible; the shell "
            f"construction is wrong.")


# =========================================================================
# residence kernels
# =========================================================================
def all_runs(occ, lo=0, hi=None):
    """Run-length encode every row of a (n_water, n_frame) occupancy matrix.

    Uses `core.geom.run_lengths` per row -- the contact-lifetime routine --
    deliberately, so "long-resident water" and "stable contact" mean the same
    thing, measured the same way, on the same grid. No second definition is
    invented here.
    """
    occ = np.asarray(occ)
    hi = occ.shape[1] if hi is None else hi
    rows, starts, lens = [], [], []
    for i in range(occ.shape[0]):
        rl = geom.run_lengths(occ[i, lo:hi].astype(bool))
        if not len(rl):
            continue
        rows.append(np.full(len(rl), i, np.int64))
        starts.append(rl[:, 0].astype(np.int64))
        lens.append(rl[:, 1].astype(np.int64))
    if not rows:
        z = np.zeros(0, np.int64)
        return z, z, z
    return np.concatenate(rows), np.concatenate(starts), np.concatenate(lens)


def survival(lens, n_win, eff_ns, max_tau_frames):
    """S(tau) from run lengths alone.

    S(tau) = the probability that a water present at time t is still
    CONTINUOUSLY present at t + tau. Counting every occupied frame as a possible
    origin, S(tau) = sum_runs max(L - tau, 0) / sum_runs L, so it follows from
    the same run lengths and needs no separate pass. S(0) = 1 by construction.
    """
    lens = np.asarray(lens, dtype=np.int64)
    tau = np.arange(0, min(max_tau_frames, n_win) + 1, dtype=np.int64)
    if lens.size == 0 or lens.sum() == 0:
        return tau * eff_ns, np.full(len(tau), np.nan)
    Lm = int(lens.max())
    c = np.bincount(lens, minlength=Lm + 1).astype(np.float64)
    ell = np.arange(Lm + 1, dtype=np.float64)
    S_c = np.concatenate([np.cumsum(c[::-1])[::-1], [0.0]])
    S_lc = np.concatenate([np.cumsum((c * ell)[::-1])[::-1], [0.0]])
    t = np.minimum(tau, Lm)      # tau > L contributes nothing; keeps the index in range
    num = S_lc[t + 1] - t.astype(np.float64) * S_c[t + 1]
    den = float((c * ell).sum())
    return tau * eff_ns, num / den


def _fit_exp(t, y):
    """Single exponential A*exp(-t/tau). Returns (tau_ns, A, r2) or nans."""
    from scipy.optimize import curve_fit
    m = np.isfinite(y) & (y > 1e-4) & (t > 0)
    if m.sum() < 5:
        return float("nan"), float("nan"), float("nan")
    try:
        # np.trapezoid, not np.trapz: np.trapz was REMOVED in NumPy 2 and this
        # environment is NumPy 2.4, where the original raises AttributeError.
        p, _ = curve_fit(lambda x, A, tau: A * np.exp(-x / tau), t[m], y[m],
                         p0=[1.0, max(float(np.trapezoid(y[m], t[m])), 0.1)],
                         bounds=([0.0, 1e-3], [10.0, 1e4]), maxfev=20000)
    except Exception:                                      # noqa: BLE001
        return float("nan"), float("nan"), float("nan")
    pred = p[0] * np.exp(-t[m] / p[1])
    ss = ((y[m] - pred) ** 2).sum()
    st = ((y[m] - y[m].mean()) ** 2).sum()
    return float(p[1]), float(p[0]), float(1 - ss / st) if st > 0 else float("nan")


def _fit_biexp(t, y):
    """A1*exp(-t/tau1) + A2*exp(-t/tau2), tau1 < tau2.

    Returns a dict of nans on failure rather than raising: a fit that will not
    converge is information, not an error.
    """
    from scipy.optimize import curve_fit
    nanout = dict(A1=float("nan"), tau1=float("nan"), A2=float("nan"),
                  tau2=float("nan"), r2=float("nan"), tau_amp=float("nan"))
    m = np.isfinite(y) & (y > 1e-4) & (t > 0)
    if m.sum() < 8:
        return nanout
    guess = max(float(np.trapezoid(y[m], t[m])), 0.1)
    try:
        p, _ = curve_fit(
            lambda x, A1, t1, A2, t2: A1 * np.exp(-x / t1) + A2 * np.exp(-x / t2),
            t[m], y[m], p0=[0.8, 0.2 * guess, 0.2, 3.0 * guess],
            bounds=([0.0, 1e-3, 0.0, 1e-3], [10.0, 1e4, 10.0, 1e4]),
            maxfev=40000)
    except Exception:                                      # noqa: BLE001
        return nanout
    A1, t1, A2, t2 = p
    if t1 > t2:
        A1, t1, A2, t2 = A2, t2, A1, t1
    pred = A1 * np.exp(-t[m] / t1) + A2 * np.exp(-t[m] / t2)
    ss = ((y[m] - pred) ** 2).sum()
    st = ((y[m] - y[m].mean()) ** 2).sum()
    tot = A1 + A2
    return dict(A1=float(A1), tau1=float(t1), A2=float(A2), tau2=float(t2),
                r2=float(1 - ss / st) if st > 0 else float("nan"),
                tau_amp=float((A1 * t1 + A2 * t2) / tot) if tot > 0
                else float("nan"))


# =========================================================================
# statistics helpers
# =========================================================================
def _drift(t, y):
    """Least-squares slope in units of y per 100 ns, plus R^2.

    Lifted from ANALYSIS/geom.drift_slope, which core/geom.py has not ported.
    A drift on an equilibrium observable is the cheapest convergence check there
    is, and leaving it out would silently drop it from every water table.
    """
    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return float("nan"), float("nan")
    t, y = t[ok], y[ok]
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    ss_res = ((y - A @ coef) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    return float(coef[0] * 100.0), \
        float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def _blocks(x, n=N_BLOCKS):
    if len(x) < n:
        return []
    e = np.linspace(0, len(x), n + 1).astype(int)
    return [(i, float(x[e[i]:e[i + 1]].mean())) for i in range(n)]


def _stats(rows, base, obs, subset, x, unit, times=None):
    """mean / sd / correlation-corrected sem / drift, and the block means."""
    x = np.asarray(x, dtype=np.float64)
    nan = float("nan")
    rows.add(observable=obs + "_mean", subset=subset,
             value=float(x.mean()) if len(x) else nan, unit=unit, **base)
    rows.add(observable=obs + "_sd", subset=subset,
             value=float(x.std(ddof=1)) if len(x) > 1 else nan, unit=unit, **base)
    rows.add(observable=obs + "_sem", subset=subset,
             value=geom.plateau_sem(x) if len(x) > 3 else nan, unit=unit, **base)
    if times is not None and len(x) > 2:
        sl, r2 = _drift(times, x)
        rows.add(observable=obs + "_drift", subset=subset, value=sl,
                 unit=unit + "/100ns", r2=r2, **base)
    blk = dict(base)
    blk["scope"] = "block"
    for bi, bv in _blocks(x):
        blk["scope_id"] = bi
        rows.add(observable=obs + "_mean", subset=f"block{bi}", value=bv,
                 unit=unit, **blk)


def _verdict(stage, framestep_ns, stride, floor):
    """(effective spacing, resolved?, floor_note) from core.spec.check_spacing.

    The stage still RUNS when the sampling cannot support the observable -- the
    tail of the residence distribution is real and is the thing worth comparing
    between systems -- but the rows say so, rather than printing a plausible
    number with no warning.
    """
    eff = float(framestep_ns or 0.0) * max(1, int(stride))
    if not framestep_ns:
        # `check_spacing` returns None for an unknown spacing, which would read
        # as "no floor to violate". It is the opposite: a duration computed from
        # an unknown frame spacing is not a duration.
        return eff, False, (f"{stage}: this replica declares no framestep_ns, "
                            f"so a residence time cannot be put in nanoseconds "
                            f"at all. Set framestep_ns on the MD.")
    try:
        got = SP.check_spacing(stage, framestep_ns, stride, floor)
        return (got if got is not None else eff), True, None
    except SP.SpacingTooCoarse as exc:
        note = str(exc)
        if int(stride) <= 1:
            note += (" At stride 1 there is no stride left to lower: this needs "
                     f"a trajectory WRITTEN every {floor:g} ns or finer.")
        return eff, False, note


# =========================================================================
# emit
# =========================================================================
def _run_base(sc, md, extra=None):
    d = dict(replica=md.id, scope="run", scope_id=md.id,
             ns_per_frame=md.framestep_ns, n_frames=sc.n_analysed,
             ns_sampled=round(sc.n_analysed * sc.eff_ns, 3))
    d.update(extra or {})
    return d


def _emit_shell(rows, spec, sc, md, rec, lig, classes, resolved, note):
    """The shell table: per frame, per run, per residue, per class.

    `resolved`/`note` are the residence verdict, carried here for the handful
    of shell rows that are secretly durations -- see the `dur` block below.
    """
    nif, nrec, nlig = sc.n_if, sc.n_rec, sc.n_lig
    prod = sc.time_ns >= float(spec.options.get("equilibration_ns") or 0.0)
    if not prod.any():
        prod = np.ones(sc.n, bool)
    times = sc.time_ns
    base = _run_base(sc, md)

    # --- per-frame -------------------------------------------------------
    # THE AUDIT, unthinned. n_rif/n_aif are the only evidence a reader has that
    # the set the two-sided search ran on was a superset of the answer. Thinning
    # them to save rows would remove exactly the frames where it was not.
    for k in range(sc.n):
        fb = dict(replica=md.id, scope="frame", scope_id=k,
                  ns_per_frame=md.framestep_ns, n_frames=1,
                  frame=int(sc.raw_frame[k]), time_ns=float(times[k]))
        rows.add(observable="n_rif_atoms", subset="all",
                 value=int(sc.n_rif[k]), unit="count", **fb)
        rows.add(observable="n_aif_atoms", subset="all",
                 value=int(sc.n_aif[k]), unit="count", **fb)

    every = max(1, int(round(FRAME_ROW_EVERY_NS / sc.eff_ns))) if sc.eff_ns else 1
    for k in range(0, sc.n, every):
        fb = dict(replica=md.id, scope="frame", scope_id=k,
                  ns_per_frame=md.framestep_ns, n_frames=1,
                  frame=int(sc.raw_frame[k]), time_ns=float(times[k]))
        for obs, x in (("n_interfacial", nif), ("n_rec_only", nrec),
                       ("n_lig_only", nlig), ("n_candidate_waters", sc.n_cand)):
            rows.add(observable=obs, subset="all", value=int(x[k]),
                     unit="count", **fb)

    # --- per run ---------------------------------------------------------
    nrec_x = (nrec - nif).astype(np.int32)
    nlig_x = (nlig - nif).astype(np.int32)
    for obs, x in (("n_interfacial", nif), ("n_rec_only", nrec),
                   ("n_lig_only", nlig), ("n_rec_excl", nrec_x),
                   ("n_lig_excl", nlig_x), ("n_rif_atoms", sc.n_rif),
                   ("n_aif_atoms", sc.n_aif),
                   ("n_candidate_waters", sc.n_cand)):
        _stats(rows, base, obs, "all", x[prod].astype(np.float64), "count",
               times=times[prod])

    with np.errstate(invalid="ignore", divide="ignore"):
        f_rec = np.where(nrec[prod] > 0, nif[prod] / np.maximum(nrec[prod], 1),
                         np.nan)
        f_lig = np.where(nlig[prod] > 0, nif[prod] / np.maximum(nlig[prod], 1),
                         np.nan)
        per_atom = np.where(sc.n_aif[prod] > 0,
                            nif[prod] / np.maximum(sc.n_aif[prod], 1), np.nan)
    for obs, x in (("frac_shared_rec", f_rec), ("frac_shared_lig", f_lig),
                   ("waters_per_lig_if_atom", per_atom)):
        good = x[np.isfinite(x)]
        _stats(rows, base, obs, "all",
               good if len(good) else np.array([np.nan]), "fraction")

    occ_prod = sc.occ[:, prod]
    n_ever = int(occ_prod.any(1).sum()) if len(occ_prod) else 0
    frames_each = occ_prod.sum(1).astype(np.float64) if len(occ_prod) \
        else np.zeros(0)
    long_res = int((frames_each * sc.eff_ns >= LONG_NS).sum())
    nan = float("nan")
    for obs, val, unit in (
            ("n_frames_no_contact", float((sc.n_rif[prod] == 0).sum()), "count"),
            ("n_frames_dry_interface", float((nif[prod] == 0).sum()), "count"),
            ("frac_frames_dry", float((nif[prod] == 0).mean()), "fraction"),
            ("n_waters_ever_interfacial", float(n_ever), "count"),
            ("n_receptor_residues_wet", float((sc.rec_wet > 0).sum()), "count"),
            ("n_ligand_residues_wet", float((sc.lig_wet > 0).sum()), "count"),
            ("n_frames_analysed", float(sc.n_analysed), "count"),
            ("n_frames_production", float(sc.n_prod), "count"),
            ("n_frames_patch_search", float(sc.n_patch), "count"),
            ("n_frames_exact_search", float(sc.n_exact), "count"),
            ("n_side_overlap_atoms", float(sc.n_overlap), "count")):
        rows.add(observable=obs, subset="all", value=val, unit=unit, **base)

    # These three are DURATION claims wearing a count's clothes: "long
    # resident" is `frames x spacing >= 1 ns`, and at a spacing of 20 ns every
    # single-frame visit clears that bar and the answer is a tautology. They
    # carry the residence verdict, not the shell's.
    dur = dict(base)
    dur["resolved"] = resolved
    dur["floor_note"] = note
    for obs, val, unit in (
            (f"n_waters_long_resident_ge_{LONG_NS:g}ns", float(long_res),
             "count"),
            (f"frac_waters_long_resident_ge_{LONG_NS:g}ns",
             float(long_res / n_ever) if n_ever else nan, "fraction"),
            ("mean_frames_per_water",
             float(frames_each[frames_each > 0].mean()) if n_ever else nan,
             "count")):
        rows.add(observable=obs, subset="all", value=val, unit=unit, **dur)

    for name, val, unit in (("d_water_shell", D_WATER_SHELL, "A"),
                            ("d_interface_reduction", D_INTERFACE, "A"),
                            ("effective_spacing", sc.eff_ns, "ns")):
        rows.add(observable=name, subset="all", value=val, unit=unit,
                 param_name=name, param_value=val, **base)
    # `text`, not `note`: `note` is the floor_note parameter, and shadowing it
    # here would silently retire the sampling-floor warning.
    for i, text in enumerate(sc.notes):
        rows.add(observable="scan_note", subset="all", value=i, unit="index",
                 label=text, **base)

    # whole-side hydration reference: how much of a side's water is doing
    # interface work, against the side's own first shell rather than a guess.
    hk = np.asarray(sc.hyd_k, np.int64)
    for tag, arr in (("rec", np.asarray(sc.hyd_rec, np.int64)),
                     ("lig", np.asarray(sc.hyd_lig, np.int64))):
        if not len(hk):
            continue
        ok = (arr >= 0) & prod[hk]
        # A -1 is a REFUSAL: that side was too big for a single-centre minimum
        # image on that frame, so no honest first-shell count exists for it.
        # Said out loud, because otherwise the only clue is a smaller n_frames
        # on the rows below and a reader would have to infer it.
        rows.add(observable=f"n_{tag}_hydration_samples_refused", subset="all",
                 value=float(((arr < 0) & prod[hk]).sum()), unit="count",
                 **dict(base, n_frames=int(prod[hk].sum())))
        if not ok.any():
            continue
        hb = dict(base)
        hb["n_frames"] = int(ok.sum())
        _stats(rows, hb, f"n_{tag}_hydration", "all", arr[ok].astype(float),
               "count")
        rows.add(observable=f"frac_{tag}_water_interfacial", subset="all",
                 value=float(nif[hk[ok]].mean() / max(1e-9, arr[ok].mean())),
                 unit="fraction", **hb)

    # --- per residue -----------------------------------------------------
    for tag, ag, wet, ever, iff in (
            ("receptor", rec, sc.rec_wet, sc.rec_ever, sc.rec_if),
            ("ligand", lig, sc.lig_wet, sc.lig_ever, sc.lig_if)):
        residues = ag.residues
        segids = np.asarray(residues.segids).astype(str)
        resids = np.asarray(residues.resids)
        resnames = np.asarray(residues.resnames).astype(str)
        show = np.flatnonzero(ever | (iff > 0))
        for j in show:
            rows.add(observable="residue_wetness", subset=tag, scope="residue",
                     scope_id=int(residues.resindices[j]), value=float(wet[j]),
                     unit="count", replica=md.id, ns_per_frame=md.framestep_ns,
                     n_frames=sc.n_prod, segid=segids[j], resid=int(resids[j]),
                     resname=resnames[j],
                     frac_frames_in_interface=float(iff[j]))

    # --- per class, from the spec's own selections -----------------------
    # `rec.residues.resindices` is the SAME order the wetness rows were
    # accumulated in (both come from `ag.residues`), so isin against it aligns;
    # sorting it separately would silently shuffle the rows.
    rec_res_order = np.asarray(rec.residues.resindices)
    for name, mask in classes.items():
        # a residue is in the class when any of its atoms is
        idx = np.asarray(rec[mask].residues.resindices)
        rmask = np.isin(rec_res_order, idx)
        rows.add(observable="receptor_wetness", subset=name, scope="region",
                 scope_id=name, value=float(sc.rec_wet[rmask].sum()),
                 unit="count", replica=md.id, ns_per_frame=md.framestep_ns,
                 n_frames=sc.n_prod, n_residues=int(rmask.sum()),
                 n_wet_residues=int((sc.rec_wet[rmask] > 0).sum()))


def _emit_bridges(rows, spec, sc, md, resid_meta):
    B = sc.bridges
    prod = sc.time_ns >= float(spec.options.get("equilibration_ns") or 0.0)
    if not prod.any():
        prod = np.ones(sc.n, bool)
    times = sc.time_ns
    base = _run_base(sc, md)
    n_prod = int(prod.sum())

    # A bridging water is hydrogen bonded to an N/O within 3.5 A on BOTH sides,
    # so it is necessarily inside the 3.5 A shell on both sides: the bridge
    # count can never exceed the interfacial water count. If it does, the search
    # is broken and nothing below is worth writing.
    if (B.n_water > sc.n_if).any():
        bad = int((B.n_water > sc.n_if).sum())
        raise AssertionError(
            f"{md.id}: bridging waters outnumber interfacial waters at {bad} "
            f"frame(s); that is impossible")

    series = [("wbridge_count", B.n_water),
              ("wbridge_respair_count", B.n_pairs),
              ("wbridge_triple_count", B.n_triples),
              ("n_interfacial_water", sc.n_if)]
    series += [(f"wbridge_count_{c}", v) for c, v in sorted(B.by_class.items())]
    for obs, x in series:
        _stats(rows, base, obs, "all", x[prod].astype(np.float64), "count",
               times=times[prod])

    ni = sc.n_if[prod].astype(np.float64)
    frac = np.divide(B.n_water[prod], ni, out=np.full(ni.shape, np.nan),
                     where=ni > 0)
    rows.add(observable="wbridge_fraction_of_interfacial_mean", subset="all",
             value=float(np.nanmean(frac)) if np.isfinite(frac).any()
             else float("nan"), unit="fraction", **base)
    rows.add(observable="wbridge_n_respairs", subset="distinct",
             value=float(len(B.pair_frames)), unit="count", **base)

    every = max(1, int(round(FRAME_ROW_EVERY_NS / sc.eff_ns))) if sc.eff_ns else 1
    for k in range(0, sc.n, every):
        fb = dict(replica=md.id, scope="frame", scope_id=k,
                  ns_per_frame=md.framestep_ns, n_frames=1,
                  frame=int(sc.raw_frame[k]), time_ns=float(times[k]))
        for obs, x in (("wbridge_count", B.n_water),
                       ("wbridge_respair_count", B.n_pairs),
                       ("wbridge_triple_count", B.n_triples)):
            rows.add(observable=obs, subset="all", value=int(x[k]),
                     unit="count", **fb)

    # --- per residue pair -------------------------------------------------
    nfp = max(1, n_prod)
    n_persist = 0
    for (rr, aa), cnt in sorted(B.pair_frames.items(), key=lambda kv: -kv[1]):
        rm, am = resid_meta.get(rr), resid_meta.get(aa)
        if rm is None or am is None:
            continue
        occ = cnt / nfp
        n_persist += int(occ >= PERSIST)
        rows.add(observable="wbridge_occupancy", subset="respair", scope="pair",
                 scope_id=f"{rm[0]}:{rm[1]}|{am[0]}:{am[1]}", value=occ,
                 unit="fraction", replica=md.id, ns_per_frame=md.framestep_ns,
                 n_frames=n_prod, segid=am[0], resid=int(am[1]),
                 resname=am[2], partner_segid=rm[0], partner_resid=int(rm[1]),
                 partner_resname=rm[2],
                 n_wbridge_mean=B.pair_waters[(rr, aa)] / nfp,
                 n_frames_bridged=int(cnt))
    rows.add(observable="wbridge_n_respairs", subset=f"occ_ge_{PERSIST}",
             value=float(n_persist), unit="count", **base)

    for tag, fr, wa in (("receptor", B.rec_frames, B.rec_waters),
                        ("ligand", B.lig_frames, B.lig_waters)):
        for rr, cnt in sorted(fr.items(), key=lambda kv: -kv[1]):
            meta = resid_meta.get(rr)
            if meta is None:
                continue
            rows.add(observable=f"wbridge_occupancy_{tag}", subset="residue",
                     scope="residue", scope_id=int(rr), value=cnt / nfp,
                     unit="fraction", replica=md.id,
                     ns_per_frame=md.framestep_ns, n_frames=n_prod,
                     segid=meta[0], resid=int(meta[1]), resname=meta[2],
                     n_wbridge_mean=wa[rr] / nfp)


def _emit_residence(rows, spec, sc, md, tag, occ, resolved, note):
    """Residence statistics for one occupancy matrix, reported three ways."""
    eq = float(spec.options.get("equilibration_ns") or 0.0)
    lo = int(np.searchsorted(sc.time_ns, eq)) if eq else 0
    hi = sc.n
    n_win = hi - lo
    eff = sc.eff_ns
    if n_win < 2 or not len(occ):
        return
    long_frames = max(1, int(round(LONG_NS / eff))) if eff else 1
    max_tau = max(1, int(round(MAX_TAU_NS / eff))) if eff else n_win

    row, start, ln = all_runs(occ, lo, hi)
    nsr = ln * eff
    total = int(ln.sum())
    edge = (start == 0) | (start + ln == n_win)
    unc = ~edge
    nan = float("nan")
    base = dict(replica=md.id, scope="run", scope_id=md.id,
                ns_per_frame=md.framestep_ns, n_frames=n_win,
                ns_sampled=round(n_win * eff, 3), resolved=resolved,
                floor_note=note)

    def add(obs, val, unit, **kw):
        rows.add(observable=obs, subset=tag, value=val, unit=unit, **base, **kw)

    add("residence_n_runs", float(len(ln)), "count")
    add("residence_n_waters", float(len(np.unique(row))) if len(row) else 0.0,
        "count")
    add("residence_water_frames_total", float(total), "count")
    add("residence_mean_ns", float(nsr.mean()) if len(nsr) else nan, "ns")
    add("residence_median_ns", float(np.median(nsr)) if len(nsr) else nan, "ns")
    add("residence_p90_ns",
        float(np.percentile(nsr, 90)) if len(nsr) else nan, "ns")
    add("residence_max_ns", float(nsr.max()) if len(nsr) else nan, "ns")
    # CENSORING. Runs touching either edge of the window are truncated, biasing
    # every mean downward. Both the censored and uncensored statistics are
    # reported, plus the share of water-time in an edge-touching run, so the
    # size of the bias is visible rather than assumed.
    add("residence_mean_uncensored_ns",
        float(nsr[unc].mean()) if unc.any() else nan, "ns")
    add("residence_median_uncensored_ns",
        float(np.median(nsr[unc])) if unc.any() else nan, "ns")
    add("residence_censored_run_fraction",
        float(edge.mean()) if len(ln) else nan, "fraction")
    add("residence_censored_time_fraction",
        float(ln[edge].sum() / total) if total else nan, "fraction")

    long_run = ln >= long_frames
    add(f"residence_n_runs_ge_{LONG_NS:g}ns", float(long_run.sum()), "count")
    add(f"residence_n_waters_ge_{LONG_NS:g}ns",
        float(len(np.unique(row[long_run]))) if long_run.any() else 0.0, "count")
    add(f"residence_frac_time_ge_{LONG_NS:g}ns",
        float(ln[long_run].sum() / total) if total else nan, "fraction")
    add(f"residence_frac_runs_ge_{LONG_NS:g}ns",
        float(long_run.mean()) if len(ln) else nan, "fraction")

    # THE SAMPLING-FLOOR GUARD. When this sits near 0.9 the mean and median
    # residence are pinned at the sampling floor and must not be quoted as
    # physical exchange times -- which is why it is reported next to them and
    # not buried in a log.
    short = ln <= 1
    add("residence_frac_runs_single_frame",
        float(short.mean()) if len(ln) else nan, "fraction")
    add("residence_frac_time_single_frame",
        float(ln[short].sum() / total) if total else nan, "fraction")
    add("residence_sampling_floor_ns", float(eff), "ns",
        param_name="one_frame_is", param_value=eff)

    n_wat = occ.shape[0]
    longest = np.zeros(n_wat, np.int64)
    nruns = np.zeros(n_wat, np.int64)
    tot = np.zeros(n_wat, np.int64)
    if len(row):
        np.maximum.at(longest, row, ln)
        np.add.at(nruns, row, 1)
        np.add.at(tot, row, ln)
    seen = tot > 0
    add("residence_waters_seen", float(seen.sum()), "count")
    add("residence_longest_mean_ns",
        float((longest[seen] * eff).mean()) if seen.any() else nan, "ns")
    add(f"residence_frac_waters_longest_ge_{LONG_NS:g}ns",
        float((longest[seen] >= long_frames).mean()) if seen.any() else nan,
        "fraction")
    add("residence_visits_per_water",
        float(nruns[seen].mean()) if seen.any() else nan, "count")
    add("water_turnover_per_ns",
        float(len(ln) / (n_win * eff)) if eff else nan, "count/ns")

    # The occupancy-versus-stability distinction, applied to water exactly as it
    # is to contacts: stable = ONE uninterrupted episode covering STABLE_FRAC of
    # the window; flicker = mean episode shorter than FLICKER_MEAN_NS. For
    # individual waters the first is a near-impossible bar -- bulk water
    # exchanges on ps-ns -- and saying so is the point of reporting it.
    mean_ep = np.divide(tot, np.maximum(nruns, 1)) * eff
    add("residence_n_waters_stable",
        float(((longest[seen] / n_win) >= STABLE_FRAC).sum()) if seen.any()
        else 0.0, "count", param_name="stable_frac", param_value=STABLE_FRAC)
    add("residence_n_waters_flicker",
        float((mean_ep[seen] < FLICKER_MEAN_NS).sum()) if seen.any() else 0.0,
        "count", param_name="flicker_mean_ns", param_value=FLICKER_MEAN_NS)
    add("residence_frac_waters_flicker",
        float((mean_ep[seen] < FLICKER_MEAN_NS).mean()) if seen.any() else nan,
        "fraction")

    occ_frac = occ[:, lo:hi].mean(axis=1)
    for thr in OCC_THRESHOLDS:
        hi_occ = (occ_frac > 0) & (occ_frac >= thr)
        stable = hi_occ & (longest >= long_frames)
        rows.add(observable="waters_occ_ge", subset=f"{tag}_occ{thr}",
                 value=float(hi_occ.sum()), unit="count",
                 n_stable=float(stable.sum()),
                 frac_stable=float(stable.sum() / hi_occ.sum())
                 if hi_occ.any() else nan, **base)

    h, _ = np.histogram(nsr, bins=HIST_EDGES_NS)
    hb = dict(base)
    hb["scope"] = "bin"
    for i in range(len(h)):
        hb["scope_id"] = i
        rows.add(observable="residence_hist", subset=tag, value=float(h[i]),
                 unit="count", bin_lo_ns=float(HIST_EDGES_NS[i]),
                 bin_hi_ns=float(HIST_EDGES_NS[i + 1])
                 if HIST_EDGES_NS[i + 1] < 1e8 else float("inf"),
                 bin_index=i, **hb)

    # SURVIVAL, reported three ways, because no one of them is trustworthy
    # alone: the model-free integral is robust and has no fit to fail; the
    # mono-exponential is systematically wrong when the decay is not
    # single-exponential, which for interfacial water it is not; the
    # bi-exponential splits it, and the amplitude split IS the answer -- the
    # fast component is the flicker and the slow one is the structural water.
    tau, Sv = survival(ln, n_win, eff, max_tau)
    ok = np.isfinite(Sv)
    add("survival_tau_integral_ns",
        float(np.trapezoid(Sv[ok], tau[ok])) if ok.sum() > 1 else nan, "ns")
    t_mono, a_mono, r2_mono = _fit_exp(tau, Sv)
    add("survival_tau_mono_ns", t_mono, "ns", r2=r2_mono, amplitude=a_mono)
    bi = _fit_biexp(tau, Sv)
    add("survival_tau_fast_ns", bi["tau1"], "ns", amplitude=bi["A1"],
        r2=bi["r2"])
    add("survival_tau_slow_ns", bi["tau2"], "ns", amplitude=bi["A2"],
        r2=bi["r2"])
    add("survival_tau_amp_weighted_ns", bi["tau_amp"], "ns", r2=bi["r2"])
    tot_amp = bi["A1"] + bi["A2"]
    add("survival_slow_amplitude_fraction",
        float(bi["A2"] / tot_amp) if tot_amp and np.isfinite(tot_amp)
        and tot_amp > 0 else nan, "fraction")
    sb = dict(base)
    sb["scope"] = "tau"
    step = max(1, len(tau) // 200)
    for i in range(0, len(tau), step):
        sb["scope_id"] = int(i)
        rows.add(observable="survival", subset=tag, value=float(Sv[i]),
                 unit="probability", tau_ns=float(tau[i]), **sb)


# =========================================================================
# driver
# =========================================================================
def run(spec, out_dir, stride=None, write=True, max_frames=None, stage=None,
        progress=0):
    """Everything water, in ONE pass per replica.

    `stage` selects which of `STAGES` to produce; the default is all three,
    which is the only arrangement that reads the trajectory once. Returns
    {stage name: Rows}.
    """
    stride = int(stride or DEFAULT_STRIDE)
    wanted = tuple(STAGES) if stage is None else \
        tuple(s for s in STAGES if s in ([stage] if isinstance(stage, str)
                                         else list(stage)))
    if not wanted:
        raise ValueError(f"stage must be one of {STAGES}")
    do_bridges = BRIDGES in wanted or RESIDENCE in wanted

    out = {s: io.Rows(system=spec.id, stage=s, stride=stride, resolved=True)
           for s in wanted}
    inputs = []

    for md in spec.mds:
        top, trajs = md.inputs(solvated=True)
        inputs += [top] + list(trajs)
        sc, rec, lig, classes = scan_md(spec, md, stride,
                                        max_frames=max_frames,
                                        do_bridges=do_bridges,
                                        progress=progress)
        # ONE verdict per replica, shared by every table: the residence floor is
        # a property of the sampling, not of which file a row lands in.
        _, resolved, note = _verdict(RESIDENCE, md.framestep_ns, stride,
                                     MIN_SPACING_NS_FOR[RESIDENCE])
        if SHELL in out:
            _emit_shell(out[SHELL], spec, sc, md, rec, lig, classes,
                        resolved, note)
        if BRIDGES in out and sc.bridges is not None:
            meta = _residue_meta(rec, lig)
            _emit_bridges(out[BRIDGES], spec, sc, md, meta)
        if RESIDENCE in out:
            _emit_residence(out[RESIDENCE], spec, sc, md, "shell",
                            sc.occ.astype(bool), resolved, note)
            if sc.bridges is not None:
                bocc, _ = _bridge_matrix(sc)
                if len(bocc):
                    # The bridged state is the one that carries a mechanism: a
                    # long-lived bridge is a structural water doing work, a
                    # long-lived neighbour may just be parked.
                    _emit_residence(out[RESIDENCE], spec, sc, md, "bridge",
                                    bocc, resolved, note)

    for s, rows in out.items():
        man = io.manifest(s, spec, inputs, stride=stride, extra={
            "family": STAGE,
            "d_water_shell": D_WATER_SHELL,
            "d_interface_reduction": D_INTERFACE,
            "d_hbond": D_HBOND, "a_hbond": A_HBOND,
            "min_spacing_ns": MIN_SPACING_NS_FOR[s],
            "max_frames": max_frames})
        if write:
            io.write(out_dir, f"{spec.label}__{s}", rows, manifest=man)
    return out


def _residue_meta(rec, lig):
    """resindex -> (segid, resid, resname) for both sides.

    Straight from the topology. The original joined a FROZEN residue index
    instead, because one of its systems carried Kabat insertion codes that its
    solvated PSF had renumbered (52A -> 901) and a label rebuilt from that
    topology moved two paratope residues out of their CDRs. That repair is
    specific to a topology this toolkit does not own; here the labels are
    whatever the topology says, and a topology that renumbers residues will be
    reported with the renumbered ids.
    """
    meta = {}
    for ag in (rec, lig):
        res = ag.residues
        for ridx, seg, rid, rn in zip(res.resindices, res.segids, res.resids,
                                      res.resnames):
            meta[int(ridx)] = (str(seg), int(rid), str(rn))
    return meta


def _bridge_matrix(sc):
    """(n_bridging_water, n_frame) occupancy of the BRIDGED state."""
    B = sc.bridges
    keys = sorted(B.wat_bridge)
    mat = np.zeros((len(keys), sc.n), bool)
    for i, g in enumerate(keys):
        mat[i, B.wat_bridge[g]] = True
    return mat, np.array(keys, dtype=np.int64)


def run_shell(spec, out_dir, **kw):
    """The shell alone. Costs its own pass -- prefer `run()`."""
    return run(spec, out_dir, stage=SHELL, **kw)[SHELL]


def run_bridges(spec, out_dir, **kw):
    """The bridges alone. Costs its own pass -- prefer `run()`."""
    return run(spec, out_dir, stage=BRIDGES, **kw)[BRIDGES]


def run_residence(spec, out_dir, **kw):
    """Residence alone. Costs its own pass -- prefer `run()`."""
    return run(spec, out_dir, stage=RESIDENCE, **kw)[RESIDENCE]
