#!/usr/bin/env python3
"""
basic.py -- RMSD, radius of gyration and per-residue RMSF, in ONE pass.

THE POINT OF ONE PASS. Reading is the budget: this NFS mount delivers ~94 MB/s
cold, so a 1.7 GB stripped replica costs ~18 s to read and a 23 GB solvated one
costs ~4 minutes. Computing RMSD, Rg and RMSF are all cheap next to that, so
they share a single walk of the trajectory and every additional observable
added to this stage is nearly free. The per-protomer work below rides on that
same walk: three extra Kabsch fits of ~632 CA and a few thousand transformed
atoms per frame, against a frame that cost milliseconds to read.

WHAT IS MEASURED, AND AGAINST WHAT. Fitting and measuring are separate
selections, because they answer different questions: the ligand's RMSD after
fitting the RECEPTOR is how far the ligand moved on its target, while fitting
the ligand to itself would report only its internal deformation. Both are
written, named by their `subset`, so neither can be mistaken for the other.

THREE PROTOMERS AND THREE REPLICAS. These systems are C3 trimers, so one
trajectory holds three copies of every measurement. Whole-assembly rows are
still written -- they are what the current pages plot -- but they answer a
different question: an assembly RMSD that climbs to 5 A can be three rigid
protomers splaying apart, and nothing in that number says so. So each protomer
is ALSO fitted on its own CA and measured on its own heavy atoms (its internal
deformation), the ligand bound to it is measured after fitting THAT protomer
(how far the ligand moved on its target), and the nine curves are pooled into
`scope="aggregate"` rows carrying mean, `sd` and `n`: the band the page draws.

SD OVER NINE, S.E.M. OVER THREE. The three protomers of one trimer share a box,
a thermostat and a barostat. They are pseudoreplicates: their agreement says
the measurement is reproducible, not that the sampling is converged. So `sd` is
the honest spread over all nine protomer-replicas (that is what a shaded band
means: this is how much the value varies between equivalent units), while `sem`
is computed over REPLICAS only, n=3, because three independent boxes are the
only independent samples there are. core/io.py's SCOPES docstring states the
same rule, and every aggregate row says `sem_over = replicas` out loud.

REPLICAS OF DIFFERENT LENGTHS ARE POOLED OVER WHAT THEY SHARE. The apo runs do
not stop at the same frame, and a mean over a time point where only one replica
still has data is not a mean over three. Pooling uses the frames every replica
has, `n_frames` on the aggregate rows is that pooled length, and a
`frames_dropped` row per replica says how much of the longer runs was left out.
Nothing is padded and nothing is truncated quietly.

SEAM FRAMES ARE DROPPED. A restart re-emits a step with different coordinates
(the window resumed from an earlier checkpoint and re-simulated); core/frames.py
finds those raw indices from the DCD headers, and they are skipped here so the
time axis holds every simulated step exactly once. On 1F4 that is 5 frames of
2,525 per replica -- small, and exactly the kind of small that turns into an
unexplained 0.2 ns spike on a plot.
"""

from core import threads                                   # noqa: F401
import os

import numpy as np

from core import frames as F
from core import geom, io, pbc
from core import protomers as PRO
from core import selections as S
from core import topology as T

STAGE = "structure"
MIN_SPACING_NS = None            # a static average needs no particular spacing
DEFAULT_STRIDE = 1

#: What to fit on, when the spec does not say. CA atoms of the receptor: the
#: assembly whose frame of reference everything else is measured in.
FIT_DEFAULT = "receptor"

#: The subsets that get the per-protomer treatment. Glycans are deliberately
#: not among them: they are not part of a protomer in any structural sense, and
#: 735 sugar residues across 75 segments would triple the table to say nothing.
PROTOMER_SUBSETS = ("receptor", "ligand")

#: How close two heavy atoms have to be to count as evidence that a ligand copy
#: sits on a given receptor protomer. Only the RANKING matters, never the
#: number, so this is a generous contact distance rather than a tuned one.
CONTACT_A = 6.0


def _universe(topology, trajectories):
    import MDAnalysis as mda
    return mda.Universe(topology, *trajectories)


#: Subsets a system may simply not have. An apo protein has no ligand, a
#: peptide has no glycans, and neither is a mistake -- so when the selection
#: is this kit's GENERIC default and matches nothing, that subset is dropped
#: and the rest of the stage runs.
#:
#: The distinction is what a person wrote. `glycans: glycan` is the default
#: and matching nothing means "no glycans here". `glycans: segid GA1 to GA23`
#: is somebody's sentence, and matching nothing means they got it wrong --
#: which stays fatal, because the lab's own TCL scripts carry the warning that
#: the wrong selection line silently gives zero contacts.
#:
#: `protein` joined this list the first time the kit met a system that has
#: none: a free aptamer, all nucleic acid. The generic `protein` selection
#: matched nothing, the structure stage raised, and the whole FOLDER stopped --
#: taking with it the records that had already built cleanly. `receptor` and
#: `solute` are deliberately NOT here: they are what everything else is fitted
#: and measured against, and a system with neither is not one this kit can say
#: anything about.
OPTIONAL_SUBSETS = ("glycans", "epitope", "ligand", "protein")


def _groups(spec, u):
    """The selections this stage measures, resolved and reported."""
    wanted = ["receptor", "ligand", "glycans", "epitope", "solute", "protein"]
    sel = {n: spec.selection(n) for n in wanted if spec.selection(n)}
    resolved, problems = S.resolve_all(u, sel)
    if problems:
        declared = getattr(spec, "declared", set())
        fatal = []
        for problem in problems:
            name = str(problem).split(":", 1)[0].strip()
            if name in OPTIONAL_SUBSETS and name not in declared:
                print(f"        {name}: none in this topology; not measured")
                continue
            fatal.append(problem)
        if fatal:
            raise S.EmptySelection("; ".join(fatal))
    return resolved


def _segids(u, resolved):
    """Every segid a resolved selection really covers.

    NOT `Resolved.segids`, which is capped at 40 for display. Reading the cap
    as the answer would drop segments from a protomer without any error -- the
    same shape of bug the selection translator's range expansion once had.
    """
    return sorted({str(s) for s in u.select_atoms(resolved.mdanalysis).segids})


def _split_assembly(u, groups, topology_path, X):
    """The protomers of this system, paired on the coordinates in `X`.

    `X` must be the first measured frame AFTER core.pbc has imaged everything
    against the receptor: a Fab that `wrapAll on` wrote a box length away has
    no contacts with anything, and the pairing would then be decided by which
    protomer happened to be nearest to nothing.

    Returns None when the topology is not one core.topology can read; the
    caller writes a row saying so rather than pretending the split happened.
    """
    if "receptor" not in groups:
        return None
    try:
        ents = T.entities(T.segments(T.read_topology(topology_path)))
    except (NotImplementedError, ValueError, OSError):
        return None
    rec = _segids(u, groups["receptor"])
    lig = _segids(u, groups["ligand"]) if "ligand" in groups else []

    def contacts(receptor_groups):
        rg = [u.select_atoms(f"segid {' '.join(g)} and {S.HEAVY}").indices
              for g in receptor_groups]
        lg = [u.select_atoms(f"segid {s} and {S.HEAVY}").indices for s in lig]
        return dict(zip(lig, PRO.contact_counts(X, rg, lg, cutoff=CONTACT_A)))

    return PRO.split(ents, rec, lig, contacts=(contacts if lig else None))


def _units(u, split):
    """One measurement kit per protomer: what to fit on, what to measure.

    The fit group is the protomer's OWN CA. That is the whole point of doing
    this per protomer: fitting on the trimer leaves the inter-protomer opening
    inside every protomer's RMSD, which is why the assembly number is always
    the larger one.
    """
    out = []
    for p in split:
        fit = u.select_atoms(f"{p['receptor']} and name CA")
        if not len(fit):                                   # pragma: no cover
            fit = u.select_atoms(f"{p['receptor']} and {S.HEAVY}")
        if not len(fit):                                   # pragma: no cover
            continue
        measure = []
        for sub in PROTOMER_SUBSETS:
            if not p.get(sub):
                continue
            ag = u.select_atoms(f"{p[sub]} and {S.HEAVY}")
            if len(ag):
                measure.append((sub, ag, ag.indices,
                                geom.RunningRMSF(len(ag))))
        out.append({"id": p["id"], "fit": fit.indices, "measure": measure})
    return out


def _fit_note(protomer_id, subset):
    """What a per-protomer RMSD row was fitted on AND measured on.

    Both halves in one field, because `subset` alone cannot tell the two apart
    for a reader who meets the row on its own: "receptor" means this protomer
    fitted on itself, "ligand" means the same fit with the cognate ligand
    measured on top of it, and those are different observables.
    """
    what = "own_heavy" if subset == "receptor" else "cognate_ligand_heavy"
    return f"protomer_{protomer_id}_receptor_CA:{what}"


def _pooled(values_by_unit):
    """(mean, sd over every unit, sem over REPLICAS, n, n_protomers,
    n_replicas).

    `values_by_unit` is {(replica, protomer): value}. The split between the two
    spreads is the point: `sd` is over all nine protomer-replicas, because that
    is what a band around a mean means; `sem` is over the three replica means,
    because protomers in one box share a thermostat and are pseudoreplicates,
    not independent samples.
    """
    vals = np.asarray(list(values_by_unit.values()), dtype=np.float64)
    reps = sorted({k[0] for k in values_by_unit})
    prots = sorted({k[1] for k in values_by_unit})
    per_rep = np.asarray(
        [np.mean([v for (r, _p), v in values_by_unit.items() if r == rep])
         for rep in reps], dtype=np.float64)
    sd = float(vals.std(ddof=1)) if len(vals) > 1 else float("nan")
    sem = (float(per_rep.std(ddof=1) / np.sqrt(len(per_rep)))
           if len(per_rep) > 1 else float("nan"))
    return (float(vals.mean()), sd, sem, int(len(vals)),
            len(prots), len(reps))


def _aggregate_rmsd(rows, series, meta):
    """Pooled RMSD vs time: one row per time point per subset.

    `series` is {(subset, replica, protomer): {frame: value}}.
    """
    for sub in sorted({k[0] for k in series}):
        keys = [k for k in series if k[0] == sub]
        common = set.intersection(*[set(series[k]) for k in keys])
        if not common:                                     # pragma: no cover
            continue
        # How much of each replica was left out, said once per replica rather
        # than implied by a shorter curve.
        for rep in sorted({k[1] for k in keys}):
            own = set()
            for k in keys:
                if k[1] == rep:
                    own |= set(series[k])
            dropped = len(own) - len(common)
            rows.add(observable="frames_dropped", subset=sub, value=dropped,
                     unit="count", replica=rep, scope="run", scope_id=rep,
                     n_frames=len(common),
                     ns_per_frame=meta["ns_per_frame"].get(rep),
                     param_name="pooled_over",
                     param_value="protomers+replicas",
                     floor_note=(None if not dropped else
                                 f"{dropped} frame(s) of this replica lie "
                                 f"beyond the last frame every replica has, "
                                 f"and are not in the pooled curve"))
        nspf = meta["ns_per_frame"].get(sorted({k[1] for k in keys})[0])
        for f in sorted(common):
            mean, sd, sem, n, n_p, n_r = _pooled(
                {(k[1], k[2]): series[k][f] for k in keys})
            rows.add(observable="rmsd", subset=sub, scope="aggregate",
                     scope_id=f, value=mean, unit="A", sd=sd, sem=sem, n=n,
                     n_frames=len(common), frame=f,
                     time_ns=meta["time_ns"].get(f), ns_per_frame=nspf,
                     replica=None, protomer=None,
                     n_protomers=n_p, n_replicas=n_r,
                     param_name="sem_over", param_value="replicas")


def _aggregate_rmsf(rows, series, meta):
    """Pooled RMSF vs residue, indexed by (entity, resid).

    EQUIVALENCE IS WITHIN AN ENTITY, NOT BY POSITION IN A LIST. Residue 384 of
    A120, B120 and C120 is the same residue of three copies of one chain;
    residue 384 of gp120 and residue 384 of a Fab heavy chain are not the same
    residue of anything. Averaging a running index would have silently averaged
    the second pair, and the resulting trace would still have looked like an
    RMSF plot.
    """
    n_used = meta["rmsf_frames"]
    lo = min(n_used.values()) if n_used else None
    uneven = len(set(n_used.values())) > 1
    for (sub, entity, resid) in sorted(series,
                                      key=lambda k: (k[0], k[1], k[2])):
        rec = series[(sub, entity, resid)]
        mean, sd, sem, n, n_p, n_r = _pooled(rec["vals"])
        rows.add(observable="rmsf", subset=sub, scope="aggregate",
                 scope_id=f"{entity}:{resid}", value=mean, unit="A",
                 sd=sd, sem=sem, n=n, n_frames=lo,
                 replica=None, protomer=None, label=entity,
                 resid=resid, resname=rec["resname"],
                 n_protomers=n_p, n_replicas=n_r,
                 param_name="sem_over", param_value="replicas",
                 floor_note=(None if not uneven else
                             "pooled over replicas of different lengths "
                             f"({min(n_used.values())}-{max(n_used.values())} "
                             f"frames): an RMSF is an average over a whole "
                             f"run, so the longer runs weigh more"))


def run(spec, out_dir, stride=None, write=True, max_frames=None):
    stride = int(stride or DEFAULT_STRIDE)
    rows = io.Rows(system=spec.id, stage=STAGE, stride=stride, resolved=True)
    inputs = []

    # What the pooled rows are built from, filled as each replica is walked so
    # the trajectories are still read exactly once.
    series = {}              # (subset, replica, protomer) -> {frame: value}
    rmsf_series = {}             # (subset, entity, resid) -> {resname, vals}
    meta = {"ns_per_frame": {}, "time_ns": {}, "rmsf_frames": {}}

    # THE SPLIT IS DECIDED ONCE, FOR THE SYSTEM, NOT ONCE PER REPLICA.
    # Replicas of one system share a topology by construction -- that is what
    # makes them replicas -- so the protomers are the same chains in all of
    # them. Re-deciding per replica would let a marginal contact flip one
    # replica's pairing, and then "protomer B" would name different chains in
    # different boxes and the pooled mean would average unlike things: a wrong
    # number with a perfectly convincing error band.
    assembly = None

    for md in spec.mds:
        top, trajs = md.inputs(solvated=False)
        if not trajs:
            continue
        inputs += [top] + list(trajs)
        u = _universe(top, trajs)
        groups = _groups(spec, u)
        meta["ns_per_frame"][md.id] = md.framestep_ns

        fit_name = FIT_DEFAULT if FIT_DEFAULT in groups else \
            next(iter(groups))
        fit = u.select_atoms(groups[fit_name].mdanalysis + " and name CA")
        if not len(fit):                                   # pragma: no cover
            fit = u.select_atoms(groups[fit_name].mdanalysis + " and " + S.HEAVY)

        # Which raw frames are seams, from the headers rather than from the
        # coordinates: the merged trajectory holds every raw frame, including
        # the re-simulated ones.
        timeline = F.build(md.id, list(md.timeline_inputs()),
                           framestep_ns=md.framestep_ns)
        seams = set(timeline.seams)
        n_raw = len(u.trajectory)
        if n_raw != len(timeline) + len(seams):
            # A merged file may already have been deduplicated; say so rather
            # than dropping frames that are not there.
            seams = set()

        measured = {name: u.select_atoms(r.mdanalysis + " and " + S.HEAVY)
                    for name, r in groups.items()}
        rmsf_of = {name: geom.RunningRMSF(len(ag))
                   for name, ag in measured.items() if len(ag)}
        # Which atoms move together when a group has to be re-imaged: by
        # segment, so a chain is never torn in half.
        # Precomputed once per MD: where each segment starts inside the group
        # and how long it is. Recomputing this per frame was the whole cost.
        seg_of = {name: pbc.runs_of(np.asarray(ag.segids))
                  for name, ag in measured.items()}
        imaged = {name: 0 for name in measured}
        split_frames = {name: 0 for name in measured}

        ref = None
        units = []
        kept = 0
        for ts in u.trajectory:
            if ts.frame in seams:
                continue
            if (kept % stride) or (max_frames and kept >= max_frames * stride):
                kept += 1
                continue
            X = u.atoms.positions.astype(np.float64)
            # THE FAB THAT FLEW AWAY. `wrapAll on` wraps every molecule into
            # the cell independently, so a bound antibody can be written a box
            # length from the trimer it is bound to. Nothing moved; its image
            # did. Each measured group is put back into the receptor's image,
            # whole and by segment, before anything is measured -- otherwise
            # the RMSD of a ligand that crosses the boundary jumps by ~200 A
            # and the renders show it in orbit.
            try:
                L = pbc.box_lengths(ts.dimensions)
            except pbc.NonOrthorhombic:
                L = None
            if L is not None:
                anchor = pbc.assembly_center(X[fit.indices])
                for name, ag in measured.items():
                    idx = ag.indices
                    # Two reductions to ask whether anything is out of place;
                    # the per-segment work only happens when it is.
                    if not pbc.needs_imaging(X[idx], anchor, L):
                        continue
                    Y2, n = pbc.shift_into_image(X[idx], anchor, L,
                                                 groups=seg_of[name])
                    if n:
                        X[idx] = Y2
                        imaged[name] += n
                    if pbc.looks_split(X[idx], L):
                        split_frames[name] += 1
            if ref is None:
                ref = X.copy()
                ref_fit = fit.indices
                # The split is decided on the first MEASURED frame, i.e. after
                # the imaging above: contacts on un-imaged coordinates are not
                # contacts.
                if assembly is None:
                    assembly = _split_assembly(u, groups, top, X)
                    if assembly is not None:
                        assembly.notes.append(
                            f"decided on {md.id} frame {kept} and reused for "
                            f"every replica: replicas of one system share a "
                            f"topology, so the protomers are the same chains")
                units = _units(u, assembly) if assembly is not None else []
            R, t, _ = geom.kabsch(X[ref_fit], ref[ref_fit])
            Y = X @ R.T + t
            time_ns = round(kept * (md.framestep_ns or 0), 4)
            meta["time_ns"].setdefault(kept, time_ns)
            base = dict(replica=md.id, scope="frame", scope_id=kept,
                        ns_per_frame=md.framestep_ns, n_frames=1,
                        frame=kept, time_ns=time_ns)
            for name, ag in measured.items():
                if not len(ag):
                    continue
                idx = ag.indices
                d = Y[idx] - ref[idx]
                rows.add(observable="rmsd", subset=name,
                         value=float(np.sqrt((d * d).sum() / len(idx))),
                         unit="A", param_name="fit_on", param_value=fit_name,
                         **base)
                rows.add(observable="rg", subset=name,
                         value=geom.radius_of_gyration(X[idx]),
                         unit="A", **base)
                rmsf_of[name].add(Y[idx])
            # Per protomer: its own fit, then its own heavy atoms and the
            # ligand sitting on it, in that one frame of reference.
            for unit in units:
                Rp, tp, _ = geom.kabsch(X[unit["fit"]], ref[unit["fit"]])
                for sub, _ag, idx, acc in unit["measure"]:
                    Yp = X[idx] @ Rp.T + tp
                    d = Yp - ref[idx]
                    val = float(np.sqrt((d * d).sum() / len(idx)))
                    rows.add(observable="rmsd", subset=sub, value=val,
                             unit="A", protomer=unit["id"],
                             param_name="fit_on",
                             param_value=_fit_note(unit["id"], sub), **base)
                    acc.add(Yp)
                    series.setdefault((sub, md.id, unit["id"]), {})[kept] = val
            kept += 1

        # Per-residue RMSF, and the run-level summary.
        for name, ag in measured.items():
            if not len(ag) or rmsf_of[name].n < 2:
                continue
            per_atom = rmsf_of[name].rmsf()
            res = ag.residues
            start = 0
            for r in res:
                n = len(r.atoms.intersection(ag))
                if not n:
                    continue
                val = float(np.sqrt((per_atom[start:start + n] ** 2).mean()))
                start += n
                rows.add(observable="rmsf", subset=name, value=val, unit="A",
                         replica=md.id, scope="residue",
                         scope_id=int(r.resindex),
                         ns_per_frame=md.framestep_ns,
                         n_frames=rmsf_of[name].n,
                         segid=str(getattr(r, "segid", "")),
                         resid=int(r.resid), resname=str(r.resname))
        # ... and the same per residue per protomer, which is what the pooled
        # trace is built from.
        entity_of = assembly.entity_of if assembly is not None else {}
        for unit in units:
            for sub, ag, _idx, acc in unit["measure"]:
                if acc.n < 2:
                    continue
                meta["rmsf_frames"][md.id] = acc.n
                per_atom = acc.rmsf()
                start = 0
                for r in ag.residues:
                    n = len(r.atoms.intersection(ag))
                    if not n:
                        continue
                    val = float(np.sqrt(
                        (per_atom[start:start + n] ** 2).mean()))
                    start += n
                    segid = str(getattr(r, "segid", ""))
                    entity = entity_of.get(segid, segid)
                    rows.add(observable="rmsf", subset=sub, value=val,
                             unit="A",
                             replica=md.id, scope="residue",
                             scope_id=int(r.resindex), protomer=unit["id"],
                             ns_per_frame=md.framestep_ns, n_frames=acc.n,
                             segid=segid, resid=int(r.resid),
                             resname=str(r.resname), label=entity)
                    key = (sub, entity, int(r.resid))
                    rec = rmsf_series.setdefault(
                        key, {"resname": str(r.resname), "vals": {}})
                    rec["vals"][(md.id, unit["id"])] = val
        for name in measured:
            vals = [r["value"] for r in rows.rows
                    if r.get("observable") == "rmsd" and r.get("subset") == name
                    and r.get("replica") == md.id and r.get("scope") == "frame"
                    and r.get("protomer") is None]
            if vals:
                a = np.asarray(vals)
                rows.add(observable="rmsd_mean", subset=name,
                         value=float(a.mean()), unit="A", replica=md.id,
                         scope="run", scope_id=md.id, n_frames=len(a),
                         ns_per_frame=md.framestep_ns,
                         ns_sampled=round(len(a) * (md.framestep_ns or 0) *
                                          stride, 3),
                         param_name="sem_plateau",
                         param_value=round(geom.plateau_sem(a), 4))
                rows.add(observable="rmsd_final", subset=name,
                         value=float(a[-1]), unit="A", replica=md.id,
                         scope="run", scope_id=md.id, n_frames=len(a),
                         ns_per_frame=md.framestep_ns)
        # The same summary per protomer, so "is one protomer different from
        # the other two?" can be answered from the run rows alone.
        for unit in units:
            for sub, _ag, _idx, _acc in unit["measure"]:
                a = np.asarray(list(series[(sub, md.id, unit["id"])].values()))
                if not len(a):                             # pragma: no cover
                    continue
                rows.add(observable="rmsd_mean", subset=sub,
                         value=float(a.mean()), unit="A", replica=md.id,
                         scope="run", scope_id=md.id, n_frames=len(a),
                         protomer=unit["id"], ns_per_frame=md.framestep_ns,
                         ns_sampled=round(len(a) * (md.framestep_ns or 0) *
                                          stride, 3),
                         param_name="fit_on",
                         param_value=_fit_note(unit["id"], sub))
                rows.add(observable="rmsd_final", subset=sub,
                         value=float(a[-1]), unit="A", replica=md.id,
                         scope="run", scope_id=md.id, n_frames=len(a),
                         protomer=unit["id"], ns_per_frame=md.framestep_ns)
        # What the split decided, and what it could not place. A protomer
        # analysis that quietly dropped a chain is worse than none.
        if assembly is not None:
            rows.add(observable="protomers", subset="assembly",
                     value=assembly.n, unit="count", replica=md.id,
                     scope="run", scope_id=md.id,
                     n_frames=kept, ns_per_frame=md.framestep_ns,
                     label="|".join(p["id"] for p in assembly),
                     param_name="pairing_rule", param_value=assembly.rule,
                     floor_note=("; ".join(assembly.notes) or None))
            for e in assembly.excluded:
                rows.add(observable="protomer_excluded", subset=e["role"],
                         value=e["copies"], unit="count", replica=md.id,
                         scope="run", scope_id=md.id, n_frames=kept,
                         ns_per_frame=md.framestep_ns,
                         label=" ".join(e["segids"]), resolved=False,
                         floor_note=e["reason"])
        elif "receptor" in groups:
            rows.add(observable="protomers", subset="assembly", value=0,
                     unit="count", replica=md.id, scope="run", scope_id=md.id,
                     n_frames=kept, ns_per_frame=md.framestep_ns,
                     resolved=False,
                     floor_note=f"{os.path.basename(top)} is not a topology "
                                f"core.topology can read, so the assembly was "
                                f"not split into protomers; only the "
                                f"whole-assembly rows are present")
        # How much imaging this trajectory needed, stated rather than hidden:
        # a large count means the file is not picture-ready and should be
        # re-imaged (cpptraj `autoimage anchor <receptor>`) before rendering.
        for name in measured:
            rows.add(observable="segments_reimaged", subset=name,
                     value=imaged[name], unit="count", replica=md.id,
                     scope="run", scope_id=md.id, n_frames=kept,
                     ns_per_frame=md.framestep_ns,
                     resolved=(split_frames[name] == 0),
                     floor_note=(None if split_frames[name] == 0 else
                                 f"{split_frames[name]} frame(s) where this "
                                 f"group still spans more than half the box "
                                 f"after imaging: check the trajectory"))

    # The rows the page draws: one mean, one band, over protomers AND replicas.
    if series:
        _aggregate_rmsd(rows, series, meta)
    if rmsf_series:
        _aggregate_rmsf(rows, rmsf_series, meta)

    man = io.manifest(STAGE, spec, inputs, stride=stride)
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows
