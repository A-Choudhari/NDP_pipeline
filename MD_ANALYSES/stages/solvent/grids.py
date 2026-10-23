#!/usr/bin/env python3
"""
grids.py -- where the water and the ions actually sit, as a map you can open.

WHY THIS STAGE EXISTS. Every other solvent number in this bundle is a scalar:
how many waters, what the salt concentration is, how far the first hydration
shell is. None of them can answer "is there a water bridging the CDR loop to
that glycan, and is it the SAME water site in every replica?" A 3D occupancy
grid can, because it keeps the position. It is also the only output here a
person can look at directly: `mol addfile <file>.dx` in VMD, Isosurface,
done.

THE GRID IS ONLY AS GOOD AS THE FRAME IT IS BINNED IN. Two things have to
happen before a single atom is binned, and skipping either produces a picture
that looks fine and means nothing:

  autoimage   NAMD runs these with `wrapAll on`, so molecules are wrapped into
              the cell independently and a water hydrogen-bonded to the Fab can
              be written a box length away. Binned as-is, it lands in a random
              voxel. `autoimage` with the receptor as anchor is the sanctioned
              imaging in this workspace (the same call the STRIP_*/strip.in
              files make), and it runs BEFORE anything else -- never a `wrap`,
              because `cellOrigin` here is not the origin.

  rms fit     The solute tumbles and drifts. Without fitting it to a fixed
              reference, a "density" is the time-average of a rotating object
              -- a spherically smeared blur of the whole box with no sites in
              it at all. The fit is on the SPEC'S `receptor` selection (CA
              atoms, as in stages/structure/basic.py), so the map is in the
              receptor's frame of reference and a peak beside the epitope is
              beside the epitope in every frame.

WHY THE VALUES ARE MULTIPLES OF BULK, NOT A^-3. An isosurface at "0.05 A^-3"
means nothing on its own: it is 1.6x bulk for water and 500x bulk for sodium,
and it is not comparable between two systems whose boxes differ by 3%. So the
grid is divided by the bulk number density of that species, MEASURED from this
run (atoms of the species from the topology, mean cell volume from the frames
that were read), and cpptraj's `normdensity density <rho>` writes the ratio
directly. An isovalue of 5 then means the same thing everywhere: five times as
much of this species as the bulk holds. The divisor is written out in a row,
because a normalised number whose divisor is not recorded cannot be undone.

WHAT A 0.5 A GRID COSTS, AND WHEN IT IS SHOT NOISE. A bulk voxel holds
`n_frames * v * rho` atoms on average -- 0.125 A^3 * 0.0312 A^-3 = 0.0039 water
oxygens PER FRAME. At stride 200 (12 frames) that is 0.05 expected counts, so
almost every voxel is 0 and the few that are 1 come out at 21x bulk: the map is
a list of places a water happened to be, not a density. Every threshold row
therefore carries `resolved` computed from that expectation, and the peak
carries `grid_peak_counts` -- how many atom-visits the tallest voxel actually
collected -- which is the number that tells you whether you have found a site.
For the ions this never goes away at 0.5 A (668 sodiums in 7.6e6 A^3 is
0.03 expected counts per voxel even at stride 1), and that is not a defect: a
bound ion that sits in one voxel for 30% of the run stands 10,000x above that
floor, which is exactly what makes it visible. The bulk of an ion grid is
noise; its peaks are sites.

COST. The extent has to cover the solute, and this trimer-plus-three-Fabs spans
177 x 170 x 152 A after padding -- 36 million voxels at 0.5 A, and a `.dx` of a
few hundred MB per species at production stride. That is what the requested
spacing costs on a solute this size; the size written is reported as a row
rather than discovered later, and `spacing=1.0` divides it by eight when the
map is for looking at rather than measuring.
"""

from core import threads                                   # noqa: F401
import os

import numpy as np

from core import cpptraj as C
from core import io
from core import selections as S
from core import spec as SP

STAGE = "density_grids"

#: A time-averaged density has no dynamic content, so no spacing is too coarse
#: -- coarser frames are MORE independent, which is better here, not worse.
#: What this stage needs is a NUMBER of frames, and that is checked per run
#: against the bulk occupancy each voxel expects (see `_resolution`).
MIN_SPACING_NS = None

#: Matches core.spec.DEFAULT_STRIDE["grids"]. At 0.2 ns/frame that is one frame
#: per 2 ns, 252 frames of a 504 ns replica -- about one expected count per
#: bulk water voxel at 0.5 A. Enough to see sites, not enough to call the bulk
#: converged, and the rows say which.
DEFAULT_STRIDE = 10

#: Angstrom. The requested default; small enough to separate two water sites
#: 2.8 A apart.
DEFAULT_SPACING = 0.5

#: Angstrom of empty space kept around the solute's own bounding box, so the
#: first two hydration shells and any ion site on the surface are inside the
#: map rather than clipped at its face.
DEFAULT_PADDING = 6.0

#: A bulk voxel must expect at least this many counts before the "voxels above
#: N x bulk" numbers mean anything: below it, a single visit by one atom is
#: already 20x bulk and the counts measure sampling, not density.
MIN_BULK_COUNTS = 5.0

#: A peak voxel must have collected at least this many atom-visits before it is
#: called a site rather than an accident.
MIN_PEAK_COUNTS = 3

#: Isovalues the summary reports, as multiples of bulk.
THRESHOLDS = (2.0, 5.0, 10.0)

#: Refuse rather than allocate: cpptraj holds the grid as doubles, so this is
#: ~1 GB per species, and three species at once is where this box starts to
#: swap. The message says what to do instead.
MAX_VOXELS = 120_000_000


def _universe(topology):
    """The topology alone. No trajectory is ever opened in Python here: the
    frames are cpptraj's job, and MDAnalysis is used only to turn this lab's
    VMD-syntax selections into atom numbers cpptraj understands."""
    import MDAnalysis as mda
    return mda.Universe(topology)


def water_oxygen(u):
    """(resnames, atom name) for water oxygen, and an error if it is not OH2.

    THE NAME IS OH2, NOT O. In CHARMM36 water the oxygen is `OH2`; `O` is a
    protein backbone carbonyl. A grid asked for `name O` inside a water mask
    comes back empty, and a mask of `name O` WITHOUT the water restriction
    comes back as the protein backbone -- 5,000 atoms, a perfectly smooth map,
    and completely the wrong molecule. So the name is looked up rather than
    assumed, and anything other than the convention is reported by name.
    """
    wat = u.select_atoms(S.translate("water"))
    if not len(wat):
        raise S.EmptySelection(
            "no water in this topology, so there is nothing to grid; a solvent "
            "stage was run on a stripped system")
    names = set(str(n) for n in wat.atoms.names)
    resnames = sorted({str(r) for r in wat.atoms.resnames})
    oxy = S.water_oxygen_name(wat)
    if oxy:
        return resnames, oxy
    heavy = sorted(n for n in names if not n.startswith("H"))
    raise S.EmptySelection(
        f"water here is {resnames} with atom names {sorted(names)}: no `OH2`. "
        f"Every stage in this workspace names water oxygen OH2; the heavy "
        f"name(s) present are {heavy}. Say which one is the oxygen rather than "
        f"letting a stage guess.")


def species_of(u):
    """The species to grid: water oxygen first, then one entry per ion resname.

    Ions are taken one resname at a time on purpose. A single `ions` grid
    superposes sodium and chloride, and the two sit in opposite places -- one
    on the acidic patches, one on the basic ones -- so their sum is a map of
    where the protein is charged at all, which is not a question anybody asked.
    """
    resnames, oxy = water_oxygen(u)
    out = [("water_" + oxy, ":" + ",".join(resnames) + "@" + oxy)]
    ions = u.select_atoms(S.translate("ions"))
    for rn in sorted({str(r) for r in ions.atoms.resnames}):
        out.append((rn, f":{rn}"))
    return out


def _bbox(pdb_path):
    """(low corner, high corner) over every model in a cpptraj-written PDB.

    Several frames spread across the run, not one: after the fit the receptor
    is fixed but the glycans and the Fab elbows are not, and a box drawn round
    frame 1 alone clips whatever moved.
    """
    lo = np.array([np.inf] * 3)
    hi = np.array([-np.inf] * 3)
    n = 0
    with open(pdb_path, errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]),
                       float(line[46:54]))
            except ValueError:                             # pragma: no cover
                continue
            lo = np.minimum(lo, xyz)
            hi = np.maximum(hi, xyz)
            n += 1
    if not n:
        raise C.CpptrajError(
            f"{os.path.basename(pdb_path)} holds no atoms: the probe pass "
            f"wrote a file but stripped everything, so the grid extent is "
            f"unknown")
    return lo, hi, n


def _counts(lo, hi, spacing, padding):
    """Bin counts and grid centre for a box that covers `lo..hi` plus padding.

    Counts are rounded UP to even numbers: cpptraj places the grid symmetrically
    about `gridcenter`, and an odd count puts the centre on a voxel face, which
    shifts the whole map by half a voxel against the coordinates the fit was
    done in.
    """
    span = (hi - lo) + 2.0 * float(padding)
    n = np.ceil(span / float(spacing)).astype(int)
    n = n + (n % 2)
    centre = (hi + lo) / 2.0
    return n, centre


def _resolution(n_frames, voxel_volume, bulk_density):
    """How many counts a BULK voxel expects. The honesty check for a grid."""
    return float(n_frames) * float(voxel_volume) * float(bulk_density)


def run(spec, out_dir, stride=None, write=True, max_frames=None,
        spacing=DEFAULT_SPACING, padding=DEFAULT_PADDING, omp=False,
        threads=None):
    """One `.dx` per species per replica, and one row per number about it."""
    stride = int(stride or DEFAULT_STRIDE)
    rows = io.Rows(system=spec.id, stage=STAGE, stride=stride, resolved=True)
    inputs = []
    runs = {}

    grids_dir = os.path.join(os.path.dirname(os.path.abspath(out_dir)), "grids")
    os.makedirs(grids_dir, exist_ok=True)

    for md in spec.mds:
        # `solvated=True` raises if only the stripped trajectory is declared --
        # gridding water on a file the water was removed from would otherwise
        # produce a perfectly empty map and no complaint.
        top, trajs = md.inputs(solvated=True)
        if not trajs:
            continue
        inputs += [top] + list(trajs)
        SP.check_spacing(STAGE, md.framestep_ns, stride, MIN_SPACING_NS)

        u = _universe(top)
        anchor, rec = C.mask_for(u, "receptor", spec.selection("receptor"))
        # Fit on CA atoms of the receptor, as stages/structure/basic.py does, so
        # the two stages share a frame of reference and an RMSD from one can be
        # read against a density from the other.
        ca = u.select_atoms(rec.mdanalysis + " and name CA")
        fit = f"({anchor})&(@CA)" if len(ca) else C.heavy(anchor)
        species = species_of(u)
        solvent_res = sorted({str(r) for r in
                              u.select_atoms(S.translate("water")).resnames} |
                             {str(r) for r in
                              u.select_atoms(S.translate("ions")).resnames})
        n_of = {name: len(u.select_atoms(_mda_for(mask)))
                for name, mask in species}

        work = C.workdir_for(out_dir, STAGE, md.id)

        # ---------------------------------------------------------- probe pass
        # One frame from each segment: enough to bound the solute over the whole
        # run and to measure the cell volume the bulk density is divided by,
        # for the price of ten frame reads (~90 MB) instead of a second pass.
        probe_pdb = os.path.join(work, "extent.pdb")
        probe_vol = os.path.join(work, "volume_probe.dat")
        probe = C.run(
            [f"parm {top}"] +
            [f"trajin {p} 1 1 1" for p in trajs] +
            [f"autoimage anchor {anchor}",
             f"rms first {fit} mass",
             f"volume out {probe_vol}",
             # Solvent and hydrogens in ONE strip, not two: a second `strip`
             # whose mask matches nothing (a system with no hydrogens) is a
             # skipped action, and a skipped action is a failed run here.
             "strip :" + ",".join(solvent_res) + "|@H*",
             f"trajout {probe_pdb} pdb"],
            work, expect=[probe_pdb, probe_vol], omp=False, name="extent")
        lo, hi, n_atoms = _bbox(probe_pdb)
        v_probe, _, _ = C.mean_volume(probe_vol)
        n, centre = _counts(lo, hi, spacing, padding)
        n_voxels = int(np.prod(n.astype(np.int64)))
        if n_voxels > MAX_VOXELS:
            raise C.CpptrajError(
                f"{md.id}: a {spacing} A grid over this solute is "
                f"{n[0]}x{n[1]}x{n[2]} = {n_voxels:,} voxels, past the "
                f"{MAX_VOXELS:,} this stage will allocate. Pass a coarser "
                f"`spacing` (1.0 A is one eighth of the voxels) or grid a "
                f"smaller system.")

        voxel_volume = float(spacing) ** 3
        # The divisor: this species' atom count (exact, from the topology) over
        # the cell volume this run measured. Not the setup script's target
        # concentration, and not a textbook 0.0334 A^-3.
        bulk = {name: n_of[name] / v_probe for name, _ in species}

        # ---------------------------------------------------------- main pass
        dx_of = {}
        lines = [f"parm {top}"] + C.trajin_lines(trajs, stride, max_frames) + [
            f"autoimage anchor {anchor}",
            f"rms first {fit} mass out {os.path.join(work, 'rmsfit.dat')}",
            f"volume out {os.path.join(work, 'volume.dat')}",
        ]
        for name, mask in species:
            dx = os.path.join(grids_dir,
                              f"{spec.label}__{md.id}__{name}.dx")
            dx_of[name] = dx
            lines.append(
                f"grid out {dx} {n[0]} {spacing} {n[1]} {spacing} "
                f"{n[2]} {spacing} gridcenter "
                f"{centre[0]:.3f} {centre[1]:.3f} {centre[2]:.3f} "
                f"{mask} normdensity density {bulk[name]:.8g}")
        res = C.run(lines, work,
                    expect=list(dx_of.values()) +
                    [os.path.join(work, "volume.dat")],
                    omp=omp, threads=threads, name="grid")
        # Both passes, because the extent the grid was drawn on is as much
        # a part of the result as the binning that followed it.
        runs[md.id] = {"grid": res.as_dict(), "extent": probe.as_dict()}

        v_mean, v_sd, n_frames = C.mean_volume(os.path.join(work,
                                                            "volume.dat"))
        if res.frames_processed and res.frames_processed != n_frames:
            n_frames = res.frames_processed
        ns_sampled = (round(n_frames * (md.framestep_ns or 0) * stride, 3)
                      if md.framestep_ns else None)
        base = dict(replica=md.id, scope="run", scope_id=md.id,
                    ns_per_frame=md.framestep_ns, n_frames=n_frames,
                    ns_sampled=ns_sampled)

        rows.add(observable="box_volume", value=v_mean, unit="A^3", **base)
        rows.add(observable="box_volume_sd", value=v_sd, unit="A^3", **base)
        rows.add(observable="grid_voxels", value=n_voxels, unit="count",
                 param_name="grid_shape",
                 param_value=f"{n[0]}x{n[1]}x{n[2]}", **base)
        rows.add(observable="grid_spacing", value=float(spacing), unit="A",
                 **base)
        rows.add(observable="solute_atoms_bounded", value=n_atoms,
                 unit="count", param_name="padding_A",
                 param_value=f"{padding:g}",
                 **base)
        for lbl, val in zip("xyz", (hi - lo)):
            rows.add(observable="solute_span", subset=lbl, value=float(val),
                     unit="A", **base)

        # The probe measured the cell on 10 frames; the run measured it on all
        # of them. The grid was already divided by the probe's number, so any
        # drift between the two is a scale error on every voxel -- stated, and
        # flagged when it is big enough to matter to an isovalue.
        drift = abs(v_mean - v_probe) / v_probe if v_probe else 0.0
        rows.add(observable="bulk_density_drift", value=float(drift),
                 unit="fraction", resolved=(drift < 0.02),
                 floor_note=(None if drift < 0.02 else
                             f"the cell volume over the analysed frames "
                             f"({v_mean:.4g} A^3) differs from the probe's "
                             f"({v_probe:.4g} A^3) by {drift:.1%}; every voxel "
                             f"value is scaled by that, so divide an isovalue "
                             f"by {v_mean / v_probe:.4f} before comparing runs"),
                 **base)

        for name, mask in species:
            g = C.read_dx(dx_of[name])
            if g.n_voxels != n_voxels:                     # pragma: no cover
                raise C.CpptrajError(
                    f"{dx_of[name]}: {g.n_voxels:,} voxels, the script asked "
                    f"for {n_voxels:,}")
            expect_bulk = _resolution(n_frames, voxel_volume, bulk[name])
            peak = float(g.values.max())
            # cpptraj wrote value = count / (n_frames * v * rho), so the raw
            # count comes straight back out. At 0.5 A a voxel holds at most one
            # atom per frame, so the count IS the number of frames that voxel
            # was occupied; at a coarse spacing it can exceed the frame count,
            # which is why this is named for counts and not for frames.
            peak_counts = peak * n_frames * voxel_volume * bulk[name]
            mb = os.path.getsize(dx_of[name]) / 1e6
            dxname = os.path.basename(dx_of[name])
            sbase = dict(subset=name, label=dxname, **base)
            need = (MIN_BULK_COUNTS / (voxel_volume * bulk[name])
                    if bulk[name] else float("inf"))
            thin = (f"a bulk voxel expects {expect_bulk:.3g} counts here; "
                    f"below {MIN_BULK_COUNTS:g} these numbers count visits, "
                    f"not density. {need:,.0f} frames at this spacing would "
                    f"be needed, against the {n_frames} read.")

            rows.add(observable="bulk_density", value=bulk[name],
                     unit="1/A^3", param_name="mask", param_value=mask,
                     **sbase)
            rows.add(observable="species_atoms", value=n_of[name],
                     unit="count", **sbase)
            rows.add(observable="grid_peak_density", value=peak,
                     unit="x_bulk", resolved=(peak_counts >= MIN_PEAK_COUNTS),
                     floor_note=(None if peak_counts >= MIN_PEAK_COUNTS else
                                 f"the tallest voxel collected only "
                                 f"{peak_counts:.1f} atom-visit(s) in "
                                 f"{n_frames} frames: that is an atom passing "
                                 f"through, not a site"),
                     **sbase)
            rows.add(observable="grid_peak_counts",
                     value=round(float(peak_counts), 2), unit="count",
                     **sbase)
            rows.add(observable="grid_peak_occupancy",
                     value=float(peak_counts / n_frames) if n_frames else None,
                     unit="atoms_per_frame", **sbase)
            rows.add(observable="grid_bulk_counts_expected",
                     value=expect_bulk, unit="count",
                     resolved=(expect_bulk >= MIN_BULK_COUNTS),
                     floor_note=(None if expect_bulk >= MIN_BULK_COUNTS
                                 else thin), **sbase)
            for t in THRESHOLDS:
                rows.add(observable="grid_voxels_above",
                         value=int((g.values >= t).sum()), unit="count",
                         param_name="threshold_x_bulk", param_value=f"{t:g}",
                         resolved=(expect_bulk >= MIN_BULK_COUNTS),
                         floor_note=(None if expect_bulk >= MIN_BULK_COUNTS
                                     else thin), **sbase)
            rows.add(observable="grid_mean_density",
                     value=float(g.values.mean()), unit="x_bulk", **sbase)
            rows.add(observable="grid_file_mb", value=round(mb, 1), unit="MB",
                     param_name="dx_file", param_value=dxname, **sbase)
            del g                                  # 36M doubles; do not keep 3

    man = io.manifest(STAGE, spec, inputs, stride=stride, extra={
        "cpptraj": runs,
        "grid_spacing_A": float(spacing), "grid_padding_A": float(padding),
        "max_frames": max_frames,
        "grids_dir": os.path.basename(grids_dir)})
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows


def _mda_for(cpptraj_mask):
    """The MDAnalysis equivalent of one of the simple masks built above.

    Only the two shapes this stage makes -- `:RES1,RES2@NAME` and `:RES` -- are
    understood, and anything else raises rather than being guessed at: the atom
    COUNT taken from here becomes the bulk density divisor, so a mask that is
    reinterpreted differently on the two sides would scale every voxel.
    """
    body = cpptraj_mask
    if not body.startswith(":"):
        raise ValueError(f"not a residue mask: {cpptraj_mask!r}")
    name = None
    if "@" in body:
        body, name = body[1:].split("@", 1)
    else:
        body = body[1:]
    sel = "resname " + " ".join(body.split(","))
    if name:
        sel += f" and name {name}"
    return sel
