#!/usr/bin/env python3
"""
rdf.py -- g(r): how solvent is structured around the solute, as one curve.

WHY THIS STAGE EXISTS. `composition.py` says the box holds 236,838 waters and
1,365 ions; `grids.py` says where a few of them sit. Neither says whether the
solvent is ORDERED, and that is the question a hydration argument stands on. A
radial distribution function answers it in a form that is comparable between
systems, between replicas and against experiment: the first peak of the
protein-polar/water-oxygen g(r) is the hydrogen-bonded first shell, and it is at
2.8 A in every folded protein ever measured. If a run puts it somewhere else,
something is wrong with the run, the force field, or the analysis -- and that is
worth knowing before any of the pretty numbers are quoted.

WHAT IS PAIRED WITH WHAT, AND WHY MORE THAN WAS ASKED.

    waterO_around_protein         all protein heavy atoms
    waterO_around_protein_polar   protein N* and O* only
    waterO_around_glycan          all glycan heavy atoms
    waterO_around_glycan_polar    glycan N* and O* only
    <ion>_around_protein          one curve per ion resname
    <ion>_around_glycan
    <anion>_around_<cation>       the ion-pairing curve

The polar curves are not a refinement, they are a different question. Measured
here on 1F4 1-R1, both put the first maximum at 2.78 A -- it is the same
hydrogen bond -- but against ALL protein heavy atoms it stands at g = 0.21 with
a prominence of 0.33, because two thirds of those atoms are carbon and the
carbon-water shell at 3.6-3.8 A is the TALLEST thing inside 6 A and buries it.
Against N*/O* the same frames give g = 0.58 at a prominence of 0.58, and the
coordination number that comes with it is waters per polar atom, which is a
quantity somebody can use. So both are written, named apart by their `subset`,
and neither can be quoted as the other. See `first_shell()` for why "the first
peak" here means the first local maximum and never the tallest one.

WHY g(r) DOES NOT GO TO 1 HERE, AND WHY THAT IS NOT A BUG. The normalisation is
cpptraj's `volume`: the expected count in a shell is N1*N2*dV/<V>, measured from
this run's own cell volume rather than assumed. For a small solute in water that
sends g(r) to 1.0 by ~10 A. This solute is a glycosylated trimer with three Fabs
bound and it fills most of its box, so most protein atoms are BURIED: at 15 A
from a buried atom there is still protein, not water. The measured tail is
~0.52, and it is climbing. That is the excluded volume of a large solute, it is
reported as `g_tail` rather than hidden, and it is the reason the first-shell
COORDINATION NUMBER (from the integral, a count of real neighbours) is the
quantity to compare between systems -- not the height of the peak.

COST, WHICH IS WHAT SETS THE STRIDE. `radial` is O(N1 x N2) per frame with no
pair list: 236,838 water oxygens against 20,265 protein heavy atoms is 4.8e9
distances for ONE frame, measured at 45 s on the serial build (8.6 s at 8
OpenMP threads, 4.7 s at 32 -- this is the stage the `omp` flag was added for).
Statistics are not the problem; a single frame already averages 4.8e9 pairs and
gives a smooth curve. So the stride is chosen by cost, not by convergence, and
the default is deliberately coarse. Frames are strided AT READ TIME by cpptraj
(`trajin <file> 1 last <stride>`), so a stride-250 pass over a 23 GB replica
transfers ~100 MB.

RESTART SEAMS ARE NOT DROPPED HERE, and that is a deliberate difference from
stages/structure/basic.py. There a re-simulated step is a visible spike on a
time axis; here every frame is poured into the same histogram, and 5 frames in
2,525 shift an equilibrium average by 0.2% of nothing. Correcting for them would
mean giving up cpptraj's read-time stride, which is the entire reason this stage
is fast.
"""

from core import threads                                   # noqa: F401
import os

import numpy as np

from core import cpptraj as C
from core import io
from core import selections as S
from core import spec as SP

STAGE = "rdf"

#: g(r) is an equilibrium average with no time axis in it, so no frame spacing
#: is too coarse -- widely spaced frames are MORE independent, which helps.
#: What it needs is pairs, and it has 4.8e9 of them per frame.
MIN_SPACING_NS = None

#: Set by cost, not by statistics: at 45 s/frame on the serial build this is
#: ~11 frames and ~16 minutes for a 504 ns replica, or ~2 minutes with
#: `omp=True`. Lowering it buys almost nothing; the curve is already smooth.
DEFAULT_STRIDE = 250

#: Angstrom. Far enough to reach the second shell and the plateau (or, on a
#: solute this size, to show that there is no plateau yet).
DEFAULT_RMAX = 15.0

#: Angstrom. Fine enough to place a 2.8 A peak to better than +/-0.05 A.
DEFAULT_DR = 0.05

#: Where a FIRST shell may be. Outside this there is no first shell, there is a
#: second one, and a "peak" found at the edge of the window is reported as the
#: non-answer it is.
FIRST_SHELL = (1.5, 6.0)

#: Bins of smoothing used to LOCATE the peak and the following minimum (0.15 A
#: at the default dr). The value reported at that r is the raw one: smoothing
#: is for finding, never for quoting.
SMOOTH_BINS = 3

#: A peak bin whose bulk expectation is below this many pairs is counting
#: coincidences. 25 pairs is ~20% shot noise on that bin.
MIN_PAIRS_PER_BIN = 25.0

#: Which ion resnames are cations, so a cation-anion curve can be paired up.
#: The same list core/spec.py's sibling stage (composition.py) uses, kept here
#: rather than imported because neither stage owns the other.
CATIONS = ("SOD", "POT", "NA", "K", "LIT", "CES", "RUB", "CAL", "MG", "MG2",
           "ZN", "ZN2", "CA", "CA2", "BAR", "FE", "FE2", "FE3", "MN", "CU",
           "CU1", "CU2")


class Sel:
    """One side of a pair: what to call it, its cpptraj mask, its atom count."""

    __slots__ = ("name", "mask", "n_atoms")

    def __init__(self, name, mask, n_atoms):
        self.name = name
        self.mask = mask
        self.n_atoms = n_atoms

    def __repr__(self):                                    # pragma: no cover
        return f"<{self.name} {self.mask} {self.n_atoms} atoms>"


def _universe(topology):
    """Topology only. The frames belong to cpptraj; MDAnalysis is here to turn
    this lab's VMD-syntax selections into atom numbers and to COUNT them, which
    is what makes an empty mask an error before cpptraj gets a chance to shrug
    it off with exit 0."""
    import MDAnalysis as mda
    return mda.Universe(topology)


def _group(u, name, selection, polar=False):
    """A solute side of a pair, as an atom-number mask restricted by name.

    Heavy atoms BY NAME (`!@H*`), never by mass -- these topologies carry broken
    masses, which is the reason core/selections.py defines HEAVY the way it
    does. `polar` keeps only N* and O*: the atoms that can donate or accept the
    hydrogen bond whose length the first peak IS.
    """
    res = S.resolve(u, name, selection)
    ag = u.select_atoms(res.mdanalysis)
    base = C.mask_from_indices(ag.indices)
    if polar:
        sub = u.select_atoms(res.mdanalysis + " and (name N* or name O*)")
        if not len(sub):
            return None
        return Sel(name + "_polar", f"({base})&(@N*,O*)", len(sub))
    return Sel(name, C.heavy(base), len(ag.select_atoms(S.HEAVY)))


def water_oxygen_sel(u):
    """Water oxygen as `:TIP3@OH2`, with the name checked rather than assumed.

    `OH2` is the convention in every CHARMM36 system here. `O` is a protein
    backbone carbonyl: a mask of `@O` intended for water silently selects the
    backbone instead -- right count, smooth curve, wrong molecule -- so the name
    is looked up and anything else is refused by name.
    """
    wat = u.select_atoms(S.translate("water"))
    if not len(wat):
        raise S.EmptySelection(
            "no water in this topology: a solvent stage was pointed at a "
            "stripped system")
    resnames = sorted({str(r) for r in wat.atoms.resnames})
    names = {str(n) for n in wat.atoms.names}
    oxy = S.water_oxygen_name(wat)
    if not oxy:
        raise S.EmptySelection(
            f"water here is {resnames} with atom names {sorted(names)}: no "
            f"`OH2` or `OW` on every water. Say which name is the oxygen "
            f"rather than letting a stage guess.")
    ox = u.select_atoms(S.translate("water") + f" and name {oxy}")
    return Sel("waterO", ":" + ",".join(resnames) + "@" + oxy, len(ox))


def ion_sels(u):
    """One Sel per ion resname present, never one lumped `ions` mask.

    Sodium and chloride sit in opposite places; their sum is a map of where the
    protein is charged at all, which answers nothing.
    """
    ions = u.select_atoms(S.translate("ions"))
    out = []
    for rn in sorted({str(r) for r in ions.atoms.resnames}):
        out.append(Sel(rn, f":{rn}", len(u.select_atoms(f"resname {rn}"))))
    return out


def pairs_for(u, spec):
    """Every (around, of) pair this system earns, in a stable order.

    `around` is cpptraj's mask1 -- the species whose distribution is being
    measured -- and `of` is mask2, the reference it is measured from. The
    subset name says it in that order (`waterO_around_protein`) so a row can
    never be read backwards.
    """
    water = water_oxygen_sel(u)
    ions = ion_sels(u)
    solutes = []
    for key in ("protein", "glycans"):
        if not spec.selection(key):
            continue
        try:
            heavy = _group(u, key, spec.selection(key))
        except S.EmptySelection:
            continue                       # an aglycosylated system has no
        solutes.append(heavy)              # glycans; that is not an error
        polar = _group(u, key, spec.selection(key), polar=True)
        if polar is not None:
            solutes.append(polar)

    out = []
    for sol in solutes:
        out.append((water, sol))
    for ion in ions:
        for sol in solutes:
            if sol.name.endswith("_polar"):
                continue                   # the ion curves are already sparse
            out.append((ion, sol))
    cations = [i for i in ions if i.name.upper() in CATIONS]
    anions = [i for i in ions if i.name.upper() not in CATIONS]
    for cat in cations:
        for an in anions:
            out.append((an, cat))
    return out, water, ions


def _smooth(y, n=SMOOTH_BINS):
    """A TRIANGULAR moving average, not a boxcar, and the difference decides
    whether a sharp peak can be found at all.

    A boxcar turns a single-bin spike into three EQUAL bins, and a local maximum
    that needs both neighbours to be strictly lower then does not exist -- so a
    curve with one sharp shell (an ion pair at a fixed distance, which is
    exactly what the test fixture is) reports no first shell. A triangular
    kernel leaves the maximum where the spike was.
    """
    n = int(n)
    if n <= 1 or y.size < n:
        return y
    if n % 2 == 0:
        n += 1
    c = n // 2
    k = np.array([c + 1 - abs(i - c) for i in range(n)], dtype=float)
    return np.convolve(y, k / k.sum(), mode="same")


def first_shell(r, g, nint, window=FIRST_SHELL, smooth=SMOOTH_BINS):
    """(r_peak, g_peak, r_min, coordination, prominence) for the FIRST shell.

    THE FIRST LOCAL MAXIMUM, NOT THE TALLEST. This distinction decides the
    answer on these systems. Against all protein heavy atoms the TALLEST feature
    inside 6 A is the carbon-water shell at 3.73 A, because two thirds of those
    atoms are carbon; the first maximum is the hydrogen bond at 2.78 A. "The
    first hydration shell" means the latter, and an `argmax` over the window
    reports the former -- with no warning, and off by a full angstrom.

    LOCATED ON A SMOOTHED COPY, QUOTED FROM THE RAW ONE. At a dozen frames an
    ion curve is spiky enough that its tallest raw bin is a coincidence while
    the shape under it is not, so the maximum is FOUND on a 3-bin average and
    then the raw maximum within that same width is what gets reported. Reporting
    the smoothed height instead would quote a number that is in no bin.

    The coordination number is the running integral at the first MINIMUM after
    the peak -- the shell boundary -- and not at a round number like 3.5 A,
    which is a different quantity for every pair and silently so.

    Returns None when there is no first shell to report. A curve that only rises
    across the window has no peak, and answering "6.0 A" because that is where
    the window ends would be inventing one.
    """
    lo, hi = window
    idx = np.nonzero((r >= lo) & (r <= hi))[0]
    if idx.size < 3 or not np.any(g[idx] > 0):
        return None
    s = _smooth(g, smooth)
    k = None
    for i in idx[1:-1]:
        # Strict on both sides, which is also what keeps a run of empty bins
        # (g = 0 below contact distance) from counting as a maximum.
        if s[i] > s[i - 1] and s[i] > s[i + 1]:
            k = int(i)
            break
    if k is None:
        return None
    w = max(1, int(smooth) // 2)
    a, b = max(0, k - w), min(g.size, k + w + 1)
    kk = int(a + np.argmax(g[a:b]))
    # Out to the first minimum of the smoothed curve: the shell boundary. The
    # comparison is STRICT, so a run of equal bins ends the walk. With `<=` a
    # curve that falls to zero and stays there -- one sharp shell and nothing
    # beyond it, which is what a two-species test case looks like -- walks to
    # r_max and reports no shell at all.
    j = k
    while j + 1 < s.size and s[j + 1] < s[j]:
        j += 1
    if j >= s.size - 1:
        return None                        # still falling at r_max: no shell
    coord = float(nint[j]) if nint is not None and j < nint.size else None
    prom = float(1.0 - s[j] / s[k]) if s[k] else 0.0
    return float(r[kk]), float(g[kk]), float(r[j]), coord, prom


def _flag(ok, note, cap_note=None):
    """`resolved` and `floor_note` in one place, so a capped probe run can never
    be read as a measurement even when its statistics look fine."""
    notes = [n for n in (None if ok else note, cap_note) if n]
    return {"resolved": bool(ok) and cap_note is None,
            "floor_note": "; ".join(notes) or None}


def run(spec, out_dir, stride=None, write=True, max_frames=None,
        r_max=DEFAULT_RMAX, dr=DEFAULT_DR, omp=True, threads=None):
    """g(r) and its integral for every pair, as tidy rows."""
    stride = int(stride or DEFAULT_STRIDE)
    rows = io.Rows(system=spec.id, stage=STAGE, stride=stride, resolved=True)
    inputs = []
    runs = {}
    cap_note = (None if not max_frames else
                f"a capped probe: only the first {max_frames} strided frame(s) "
                f"of the first segment were read, so this samples the head of "
                f"the run and not the run")

    for md in spec.mds:
        # Solvent stages must read the SOLVATED trajectory; asking for it here
        # means a spec that only declares the stripped one fails loudly instead
        # of returning a g(r) of nothing.
        top, trajs = md.inputs(solvated=True)
        if not trajs:
            continue
        inputs += [top] + list(trajs)
        SP.check_spacing(STAGE, md.framestep_ns, stride, MIN_SPACING_NS)

        u = _universe(top)
        pairs, water, ions = pairs_for(u, spec)
        anchor, _ = C.mask_for(u, "receptor", spec.selection("receptor"))
        work = C.workdir_for(out_dir, STAGE, md.id)
        vol_path = os.path.join(work, "volume.dat")

        # `autoimage` cannot change any number below -- `radial` takes
        # minimum-image distances itself, and r_max (15 A) is well under half
        # the shortest cell edge (176 A), so the convention is exact. It is
        # here because every pass over a solvated trajectory in this workspace
        # is imaged the same way, and because it costs ~0.05 s against the 45 s
        # a single `radial` frame costs.
        lines = [f"parm {top}"] + C.trajin_lines(trajs, stride, max_frames) + [
            f"autoimage anchor {anchor}",
            f"volume out {vol_path}",
        ]
        expect = [vol_path]
        files = {}
        for a, b in pairs:
            label = f"{a.name}_around_{b.name}"
            g_path = os.path.join(work, f"g_{label}.dat")
            n_path = os.path.join(work, f"n_{label}.dat")
            files[label] = (a, b, g_path, n_path)
            lines.append(
                f"radial out {g_path} {dr} {r_max} {a.mask} {b.mask} "
                f"volume intrdf {n_path}")
            expect += [g_path, n_path]

        res = C.run(lines, work, expect=expect, omp=omp, threads=threads,
                    name="radial")
        runs[md.id] = res

        v_mean, v_sd, n_vol = C.mean_volume(vol_path)
        n_frames = res.frames_processed or n_vol
        ns_sampled = (round(n_frames * (md.framestep_ns or 0) * stride, 3)
                      if md.framestep_ns else None)
        base = dict(replica=md.id, ns_per_frame=md.framestep_ns,
                    n_frames=n_frames, ns_sampled=ns_sampled)
        rbase = dict(scope="run", scope_id=md.id, **base)

        rows.add(observable="box_volume", value=v_mean, unit="A^3",
                 **_flag(True, None, cap_note), **rbase)
        rows.add(observable="box_volume_sd", value=v_sd, unit="A^3",
                 **_flag(True, None, cap_note), **rbase)
        # The bulk number density of each species, measured: atoms from the
        # topology over the cell volume these frames had. For TIP3P water this
        # is the 0.031-0.033 A^-3 that says the box is liquid water and not
        # something that went wrong in equilibration.
        for sel in [water] + ions:
            rows.add(observable="bulk_density", subset=sel.name,
                     value=sel.n_atoms / v_mean, unit="1/A^3",
                     param_name="atoms", param_value=str(sel.n_atoms),
                     **_flag(True, None, cap_note), **rbase)

        for label, (a, b, g_path, n_path) in files.items():
            _, gd = C.read_datafile(g_path)
            _, nd = C.read_datafile(n_path)
            r, g = gd[:, 0], gd[:, 1]
            nint = nd[:, 1] if nd.shape[1] > 1 else None
            if nint is not None and nint.size != r.size:    # pragma: no cover
                raise C.CpptrajError(
                    f"{label}: g(r) has {r.size} bins and its integral has "
                    f"{nint.size}; they cannot be the same curve")

            for i in range(r.size):
                brow = dict(subset=label, scope="bin", scope_id=str(i),
                            r_A=round(float(r[i]), 4), **base)
                rows.add(observable="g_r", value=float(g[i]), unit="ratio",
                         **_flag(True, None, cap_note), **brow)
                if nint is not None:
                    rows.add(observable="coordination",
                             value=float(nint[i]), unit="count",
                             **_flag(True, None, cap_note), **brow)

            tail = float(g[r >= (r_max - 1.0)].mean()) if r.size else None
            rows.add(observable="g_tail", subset=label, value=tail,
                     unit="ratio", param_name="over_last_A", param_value="1.0",
                     **_flag(True, None, cap_note), **rbase)

            shell = first_shell(r, g, nint)
            if shell is None:
                rows.add(observable="first_peak_r", subset=label, value=None,
                         unit="A",
                         **_flag(False,
                                 "no first-shell maximum between "
                                 f"{FIRST_SHELL[0]} and {FIRST_SHELL[1]} A: "
                                 f"this curve only rises, so there is nothing "
                                 f"here to call a shell", cap_note), **rbase)
                continue
            r_peak, g_peak, r_min, coord, prom = shell
            # The bulk expectation for the peak bin, which is the noise scale:
            # N1*N2*4*pi*r^2*dr/V per frame.
            per_frame = (a.n_atoms * b.n_atoms * 4.0 * np.pi * r_peak ** 2 *
                         dr / v_mean)
            n_pairs = per_frame * n_frames
            thin = (f"the peak bin expects only {n_pairs:.0f} pairs from bulk "
                    f"({MIN_PAIRS_PER_BIN:g} is the floor here), so its height "
                    f"and position are shot noise; read more frames")
            ok = n_pairs >= MIN_PAIRS_PER_BIN
            sb = dict(subset=label, label=os.path.basename(g_path), **rbase)
            rows.add(observable="first_peak_r", value=r_peak, unit="A",
                     param_name="pairs_expected_in_bin",
                     param_value=f"{n_pairs:.1f}",
                     **_flag(ok, thin, cap_note), **sb)
            rows.add(observable="first_peak_g", value=g_peak, unit="ratio",
                     param_name="prominence", param_value=f"{prom:.3f}",
                     **_flag(ok, thin, cap_note), **sb)
            rows.add(observable="first_shell_coordination", value=coord,
                     unit="count", param_name="first_min_r",
                     param_value=f"{r_min:.3f}",
                     **_flag(ok, thin, cap_note), **sb)

    man = io.manifest(STAGE, spec, inputs, stride=stride, extra={
        "cpptraj": {k: v.as_dict() for k, v in runs.items()},
        "r_max_A": float(r_max), "dr_A": float(dr),
        "max_frames": max_frames})
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows
