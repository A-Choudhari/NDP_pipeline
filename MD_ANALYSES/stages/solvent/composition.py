#!/usr/bin/env python3
"""
composition.py -- how much water, how many ions, how big the box.

The first thing a reader of a record wants to know about a solvated system, and
until now nowhere in this workspace: `config2.py` carries N_WATER/N_SOD/N_CLA as
hand-typed constants for three systems, which is a table that can go out of date
without anyone noticing.

WHERE THE NUMBERS COME FROM.

  counts   the topology itself (core/topology.composition), so they are exact
           and cost one PSF parse -- 236,838 TIP3 waters, 668 SOD, 697 CLA for
           1F4, matching the PSF's own residue records.

  box      NAMD's `.xst` files, which record the periodic cell every
           `xstFreq` step. They are text, a few hundred kB, and denser in time
           than the trajectory itself -- so the box volume series is free,
           where reading it from the DCD frames would cost a seek per frame
           (~0.2 s each on this NFS mount: 8 minutes per replica for 56 bytes
           a time).

  salt     computed from the ion counts and the mean box volume, not from what
           the setup script was asked for. `-sc 0.15` is an instruction; this
           is the measurement.

The `.xst` step column is the same step counter the DCD headers use, so the box
series is trimmed to the production window rather than including heating and
equilibration, which have a different (and much larger) box.
"""

import glob
import os

from core import io
from core import topology as T

STAGE = "composition"
#: One frame per 20 ns is plenty for a mean box volume; the `.xst` is denser
#: anyway. Declared so the stage cannot be run at a spacing that says nothing.
MIN_SPACING_NS = None            # a time average has no dynamic content
DEFAULT_STRIDE = 100

#: Avogadro's number x 1e-27 (A^3 -> litres), so N / (V[A^3] * this) is molar.
_MOL_PER_A3 = 6.02214076e23 * 1e-27


def read_xst(path):
    """(step, a, b, c) for every line of a NAMD .xst file.

    Only the diagonal is read: these boxes are orthorhombic (`useFlexibleCell
    no`), and a triclinic cell would need the full nine components -- which is
    checked, not assumed.
    """
    out = []
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 10:
                continue
            try:
                step = int(float(f[0]))
                ax, ay, az = float(f[1]), float(f[2]), float(f[3])
                bx, by, bz = float(f[4]), float(f[5]), float(f[6])
                cx, cy, cz = float(f[7]), float(f[8]), float(f[9])
            except ValueError:
                continue
            off = max(abs(ay), abs(az), abs(bx), abs(bz), abs(cx), abs(cy))
            if off > 1e-6 * max(ax, by, cz):
                raise ValueError(
                    f"{os.path.basename(path)} step {step}: the cell is not "
                    f"orthorhombic (off-diagonal {off:.3g}); this stage reads "
                    f"only the diagonal and would be silently wrong")
            out.append((step, ax, by, cz))
    return out


def box_series(run_dir, first_step=None, last_step=None, pattern="pr*.xst"):
    """The production box, from every `.xst` beside the trajectory."""
    rows = []
    for p in sorted(glob.glob(os.path.join(run_dir, pattern))):
        if ".BAK" in p or "/old" in p.lower():
            continue
        rows += read_xst(p)
    rows.sort(key=lambda r: r[0])
    seen, out = set(), []
    for step, a, b, c in rows:
        if step in seen:
            continue                      # a restart re-emits its first step
        if first_step is not None and step < first_step:
            continue
        if last_step is not None and step > last_step:
            continue
        seen.add(step)
        out.append((step, a, b, c))
    return out


def _stats(values):
    n = len(values)
    if not n:
        return None, None
    mean = sum(values) / n
    if n < 2:
        return mean, 0.0
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return mean, var ** 0.5


def run(spec, out_dir, stride=None, write=True):
    """One row per number, per replica, plus the system's own summary."""
    stride = int(stride or DEFAULT_STRIDE)
    rows = io.Rows(system=spec.id, stage=STAGE, stride=stride, resolved=True)
    inputs, per_replica = [], []

    for md in spec.mds:
        top = md.solvated_topology or md.topology
        atoms = T.read_topology(top)
        segs = T.segments(atoms)
        comp = T.composition(atoms, segs)
        inputs.append(top)

        run_dir = os.path.dirname(md.solvated_trajectories[0]
                                  if md.solvated_trajectories
                                  else md.trajectories[0])
        box = box_series(run_dir)
        vols = [a * b * c for _, a, b, c in box]
        vmean, vsd = _stats(vols)
        ameans = [_stats([r[i] for r in box])[0] for i in (1, 2, 3)]

        base = dict(replica=md.id, scope="run", scope_id=md.id,
                    ns_per_frame=md.framestep_ns, n_frames=len(box),
                    ns_sampled=None)
        rows.add(observable="water_molecules", value=comp["water_molecules"],
                 unit="count", **base)
        for name, n in comp["ions_by_name"].items():
            rows.add(observable="ions", subset=name, value=n, unit="count",
                     **base)
        rows.add(observable="ions_total", value=comp["ions_total"],
                 unit="count", **base)
        rows.add(observable="atoms_total", value=comp["atoms_total"],
                 unit="count", **base)
        rows.add(observable="solute_atoms", value=comp["solute_atoms"],
                 unit="count", **base)
        rows.add(observable="glycan_residues", value=comp["glycan_residues"],
                 unit="count", **base)
        if vmean:
            for label, val in zip("abc", ameans):
                rows.add(observable="box_edge", subset=label, value=val,
                         unit="A", **base)
            rows.add(observable="box_volume", value=vmean, unit="A^3", **base)
            rows.add(observable="box_volume_sd", value=vsd, unit="A^3", **base)
            # Salt as measured: ions over the box they actually occupy.
            cations = sum(n for r, n in comp["ions_by_name"].items()
                          if r in ("SOD", "POT", "NA", "K", "LIT", "CES",
                                   "RUB", "CAL", "MG", "ZN", "ZN2", "CA"))
            anions = comp["ions_total"] - cations
            pairs = min(cations, anions)
            rows.add(observable="salt_concentration", value=
                     pairs / (vmean * _MOL_PER_A3), unit="mol/L", **base)
            rows.add(observable="water_density", value=
                     comp["water_molecules"] / vmean, unit="1/A^3", **base)
        per_replica.append((md.id, comp, vmean, vsd, len(box)))

    # The system's own line: replicas of one system share a topology, so a
    # disagreement in composition means they are not replicas.
    if per_replica:
        first = per_replica[0][1]
        for md_id, comp, *_ in per_replica[1:]:
            if comp["water_molecules"] != first["water_molecules"] or \
                    comp["ions_by_name"] != first["ions_by_name"]:
                rows.add(replica=md_id, scope="system", observable="mismatch",
                         value=1, unit="flag", resolved=False,
                         floor_note=("this replica's composition differs from "
                                     f"{per_replica[0][0]}: they cannot be "
                                     "replicas of one system"))

    man = io.manifest(STAGE, spec, inputs, stride=stride)
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows


def table(rows):
    """The composition as a small table, for the record and the page."""
    out = {}
    for r in rows.rows:
        if r.get("scope") != "run":
            continue
        key = r["observable"] + (f"_{r['subset']}" if r.get("subset") else "")
        out.setdefault(key, {})[r["replica"]] = r["value"]
    return out
