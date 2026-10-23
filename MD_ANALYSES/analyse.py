#!/usr/bin/env python3
"""
analyse.py -- analyses for one folder's records, from the glossary it already has.

    python3 analyse.py <folder> [--system CDR_1-1F4_3x500ns] [--only structure]
                                [--stride 1] [--detect-only]

    <folder>/MD_RESULTS/analyses/<record>/data/*.tsv *.parquet   the numbers
    <folder>/MD_RESULTS/analyses/<record>/frames.json            the time axis
    <folder>/MD_RESULTS/analyses/<record>/page_plotly.html       the page
    <folder>/MD_RESULTS/analyses/<record>/figures_user/          YOUR pictures

WHERE THE INPUTS COME FROM. The metadata pipeline has already read every
trajectory header and written `MD_RESULTS/glossary.parquet`: the topology, the
production segments in time order, the frame spacing, the atom count. This
reads that, so the two never disagree about what an MD is -- and so no
trajectory is scanned twice.

STRIPPED VERSUS SOLVATED. Structure stages read the solute-only trajectory
where one exists (1.7 GB for 1F4, against 23 GB solvated): thirteen times less
I/O for the same geometry. Solvent stages must read the solvated one. Both are
recorded per MD, and a stage asks for what it needs -- it never guesses. The
stripped file is found by convention (`STRIP*/merged.dcd`, newest generation,
with `merged.psf` beside it) and NAMED in the spec, because two generations of
that directory exist side by side in this very folder.
"""

from core import threads                                   # noqa: F401
import argparse
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import detect, frames as F, spec as SP           # noqa: E402

RESULTS = "MD_RESULTS"
PLOT_SUBSETS = ("receptor", "ligand")
ANALYSES = "analyses"


def _settings(folder):
    """Where the person declares what only they know, for THIS folder.

    It had never been defined. The one line that uses it runs only when a
    record has no readable frame spacing, so the NameError sat here until a
    folder finally produced one -- and then it did not skip that record, it
    ended the whole run: HD4 stopped at the tenth of twenty-five, inside the
    message explaining the skip.
    """
    return os.path.join(_results(folder), "settings.yaml")


def _results(folder):
    return os.path.join(folder, RESULTS)


def load_glossary(folder):
    import pandas as pd
    p = os.path.join(_results(folder), "glossary.parquet")
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"no glossary at {p} -- run the NDP pipeline on {folder} first")
    return pd.read_parquet(p)


def load_records(folder):
    """The published records: file stem -> (record document)."""
    out = {}
    for d in (_results(folder), os.path.join(_results(folder), "NDP_records")):
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if not name.endswith(".json") or name == "catalog.json":
                continue
            with open(os.path.join(d, name)) as fh:
                try:
                    doc = json.load(fh)
                except Exception:                          # noqa: BLE001
                    continue
            if isinstance(doc, dict) and doc.get("lab_glossary"):
                out[name[:-5]] = doc
    return out


def stripped_for(run_dir, n_expected=None, root=None, n_solvated=None):
    """The solute-only trajectory beside a run, or (None, None).

    CHOSEN BY CONTENT, NOT BY DATE. A replica here can hold three generations
    of `STRIP*/merged.dcd`, and the older ones cover only part of the run --
    JL427/3R's held 1,293 frames of 2,526 until it was re-made. Picking "the
    newest directory" would analyse a third of a trajectory and say nothing.
    So the candidate whose frame count matches the raw segments wins; if none
    does, the longest one is returned together with a warning the caller can
    print, and nothing is chosen silently.
    """
    from core.dcdheader import read_dcd_header
    # ANY .dcd under a STRIP* directory, not just one called merged.dcd. The
    # occluded-open replicas call theirs `STRIP/OCCOPEN_1R_stripped.dcd`, and
    # a kit that only recognises one lab member's filename is a kit that works
    # on one folder.
    cands = sorted(glob.glob(os.path.join(run_dir, "STRIP*", "*.dcd")))
    if not cands:
        # A DIRECTORY IS NOT THE ONLY PLACE A STRIP LIVES. `8-IgG` and the
        # whole Aptamers corpus keep theirs loose in the replica as
        # `merged.dcd`, beside the solvated run it was made from -- and with
        # no STRIP*/ to glob, the structure stages were reading 581,016-atom
        # production files to measure a 20,862-atom solute: thirteen times the
        # I/O for the same geometry.
        #
        # The filter is the same one the metadata half already applies (its
        # `atoms` basis): a candidate qualifies only if it holds FEWER atoms
        # than the solvated run. That is what keeps `heating.dcd`, `npt2.dcd`
        # and `min1.dcd` out -- they carry the production atom count exactly,
        # so a name-based rule would let them through and a frame-count rule
        # would sometimes prefer them (X-MODELS' npt2.dcd holds 200 frames
        # where each production segment holds 100).
        if not n_solvated:
            return None, None
        loose = []
        for c in sorted(glob.glob(os.path.join(run_dir, "*.dcd"))):
            try:
                h = read_dcd_header(c)
            except (OSError, ValueError):
                continue
            if (h.get("natom") or 0) and h["natom"] < n_solvated:
                loose.append(c)
        if not loose:
            return None, None
        cands = loose
    best, best_n, best_h = None, -1, None
    for c in sorted(cands, key=lambda q: os.path.basename(os.path.dirname(q)),
                    reverse=True):
        h = read_dcd_header(c)
        n = h.get("n_frames") or 0
        if n_expected and n == n_expected:
            best, best_h = c, h
            break
        if n > best_n:
            best, best_n, best_h = c, n, h
    top = solute_topology(run_dir, (best_h or {}).get("natom"), root)
    return (top, best) if top else (None, None)


def solute_topology(run_dir, natom, root=None):
    """A topology with `natom` atoms, for a stripped trajectory.

    MATCHED BY ATOM COUNT, and searched beyond the replica's own directory.
    `<run_dir>/merged.psf` is the convention here and is tried first, but the
    occluded-open runs keep theirs in `MDs/1-INPUTs/merged.psf` -- two levels
    up and in a sibling of the replica -- because one solute topology serves
    all three replicas. Guessing by name would pick the 712k-atom
    `solvated_ions.psf` sitting right next to the trajectory; the atom count
    cannot be guessed wrong, and it is already in the DCD header.

    THE SEARCH STOPS ONE DIRECTORY ABOVE THE FOLDER YOU NAMED. Without that
    bound, a replica sitting directly in a project's own folder climbs to the
    shared parent -- `ABZYME/`, holding every project in the lab -- and an
    atom count that matches by coincidence would pair a trajectory with some
    other system's topology. Nothing would raise: the atom counts agree, so
    every analysis would run and every number would be wrong. A topology
    further away than that is not "beside the data", and belongs in
    `settings.yaml` instead.
    """
    if not natom:
        return None
    here = os.path.join(run_dir, "merged.psf")
    if os.path.exists(here) and psf_natoms(here) == natom:
        return here
    up1 = os.path.dirname(run_dir.rstrip("/"))
    up2 = os.path.dirname(up1)
    limit = os.path.dirname(os.path.abspath(root)) if root else None
    seen = set()
    for base in (run_dir, up1, up2):
        if not base or base in seen or not _within(base, limit):
            continue
        seen.add(base)
        for pattern in ("*.psf", os.path.join("*", "*.psf")):
            for cand in sorted(glob.glob(os.path.join(base, pattern))):
                # Every candidate is read, and reading one costs a single
                # 1 MB head even for the 96 MB solvated topology. An earlier
                # version skipped files below `natom * 60` bytes on the
                # grounds that a topology cannot be smaller than its atom
                # table -- true of the lab's files, false of anything with a
                # short header, and it made the rule impossible to test
                # without writing megabytes of fixture.
                if psf_natoms(cand) == natom:
                    return cand
    return None


def _within(path, limit):
    """Is `path` inside `limit`? No limit means anywhere."""
    if not limit:
        return True
    path, limit = os.path.abspath(path), os.path.abspath(limit)
    return path == limit or path.startswith(limit + os.sep)


def psf_natoms(path, budget=1 << 20):
    """The `!NATOM` count, from the head of the file and nothing more.

    The budget is a megabyte because a PSF header is not short: psfgen writes
    one REMARKS line per patch, and the CDR topologies carry 896 of them
    before `!NATOM`. A first attempt read 80 lines, found nothing, and quietly
    decided those systems had no solute topology at all -- so three of the
    four fell back to analysing the solvated trajectory.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(budget)
    except OSError:
        return None
    m = re.search(rb"(\d+)\s*!NATOM", head)
    return int(m.group(1)) if m else None


def solvated_topology_for(run_dir, n_published, root=None):
    """The topology of the solvated build, when the PUBLISHED run is the solute.

    A record can name the solute trajectory as its production data: 8-IgG's
    1F4 does, because one replica's raw segments no longer exist and the
    merged solute copy is the only complete record of it. The glossary then
    hands over `merged.psf`, 20,862 atoms, and every question about the
    solvent becomes unanswerable -- the system detects as `protein, glycans`,
    `composition` never runs, and the record would say nothing about the
    186,366 waters and 1,056 ions the box was built with.

    They are not lost, only unnamed. The solvated topology is still beside the
    run, and which one it is is EVIDENCE rather than a guess: it is the one
    whose atom count matches the largest trajectory in the replica. That rules
    out `1-INPUTs/solvated.psf` (583,128 atoms), a build that matches nothing
    here and would win any "take the biggest topology" rule.
    """
    from core.dcdheader import read_dcd_header
    biggest = 0
    for c in glob.glob(os.path.join(run_dir, "*.dcd")):
        try:
            biggest = max(biggest, read_dcd_header(c).get("natom") or 0)
        except (OSError, ValueError):
            continue
    if not biggest or (n_published and biggest <= n_published):
        return None
    return solute_topology(run_dir, biggest, root)


def specs_from(folder, only_system=None):
    """One Spec per published record, with its MDs resolved from the glossary."""
    g = load_glossary(folder)
    rows = {str(r.md_key): r for _, r in g.iterrows()}
    out = []
    for stem, doc in load_records(folder).items():
        if only_system and only_system not in (stem, doc.get("identifier")):
            continue
        mds = []
        for key in doc["lab_glossary"]["md_keys"]:
            r = rows.get(str(key))
            if r is None:
                continue
            root, rel = str(r.source_root), str(r.source_dir)
            run_dir = os.path.join(root, rel)
            segs = [os.path.join(run_dir, s)
                    for s in str(r.segments).split("; ") if s]
            top_solv = os.path.join(root, str(r.topology))
            # What the production segments really hold, so a stale merged copy
            # cannot be mistaken for the whole run. UNIQUE frames, from the
            # timeline, not the sum of the headers: a run that restarted
            # re-simulates the checkpoint frame, so the segments of the
            # occluded-open replicas sum to 10,081 where the run is 10,080
            # frames long -- and the strip, which drops that duplicate with
            # `trajin ... 2 last`, was being reported as a short strip.
            from core.dcdheader import read_dcd_header
            try:
                n_raw = len(F.build(str(r.md_id), segs).frames)
            except (ValueError, OSError):
                n_raw = sum((read_dcd_header(q).get("n_frames") or 0)
                            for q in segs)
            try:
                n_solv = read_dcd_header(segs[0]).get("natom") if segs else None
            except (OSError, ValueError):
                n_solv = None
            top_dry, dry = stripped_for(run_dir, n_expected=n_raw,
                                        root=folder, n_solvated=n_solv)
            if dry:
                h = read_dcd_header(dry)
                if (h.get("n_frames") or 0) != n_raw:
                    print(f"      WARNING {r.md_id}: "
                          f"{os.path.relpath(dry, run_dir)} holds "
                          f"{h.get('n_frames')} frames, the segments hold "
                          f"{n_raw} unique frames. Re-run the strip, or the "
                          f"analysis covers only part of this replica.")
            # THE FRAME SPACING, AND WHERE IT IS ALLOWED TO COME FROM.
            # The glossary states it when the trajectory headers do. When
            # the production files are post-processed copies -- the whole
            # FULL-GP120 folder is -- the headers carry a placeholder and
            # the glossary column is empty; the person then declares the
            # interval in settings.yaml, which reaches the RECORD. Reading
            # it from there too is not a second source of truth: it is the
            # same declaration, and without it the analysis would invent
            # 4.9e-05 ns per frame from a placeholder and put every point of
            # every curve at t = 0.
            step = (float(r.framestep_ns)
                    if r.framestep_ns == r.framestep_ns else None)
            declared = (doc.get("sampling") or {}).get("frame_interval_ps")
            if step is None and declared:
                step = float(declared) / 1000.0
            # When the published segments ARE the solute, the glossary's
            # topology is the solute one too, and `solvated` would point at
            # itself. Recover the real solvated topology from the replica.
            solv_top, solv_trajs = top_solv, segs
            if not dry:
                found_solv = solvated_topology_for(run_dir, n_solv, root=folder)
                if found_solv:
                    # The topology, but NOT the trajectories: the solvated run
                    # is not what this record publishes, and on 3R most of it
                    # no longer exists. `composition` needs only the topology
                    # and the .xst series beside the run, so the waters, the
                    # ions, the box and the salt concentration are all still
                    # measured. A stage that genuinely needs solvated frames
                    # says so and stops, which is the right answer here.
                    solv_top, solv_trajs = found_solv, []
            mds.append(SP.MD(
                id=str(r.md_id), md_key=str(key),
                topology=top_dry or top_solv,
                trajectories=[dry] if dry else segs,
                framestep_ns=step,
                solvated_topology=solv_top, solvated_trajectories=solv_trajs))
        if not mds:
            continue
        spec, found = detect.propose(doc["lab_glossary"]["record_id"], mds,
                                     label=stem)
        spec.options["record"] = stem
        spec.options["title"] = doc.get("title") or stem
        out.append((spec, found, doc))
    return out


class _Rows:
    """The little a page needs from a Rows object when reading tables back."""

    def __init__(self, rows):
        self.rows = rows

    def __len__(self):
        return len(self.rows)


def bundle_dir(folder, stem):
    return os.path.join(_results(folder), ANALYSES, stem)


def write_frames(spec, out_dir):
    """`frames.json` -- the time axis of every replica, for the pages."""
    tls = []
    for md in spec.mds:
        src = list(md.timeline_inputs())
        if not src:
            continue
        t = F.build(md.id, src, framestep_ns=md.framestep_ns)
        # What the analysis actually read, which may be the merged copy whose
        # frame i is raw frame i of the segments above.
        t.analysed = [os.path.basename(p) for p in md.trajectories]
        t.analysed_path = md.trajectories[0] if md.trajectories else None
        if t.analysed_path:
            from core.dcdheader import read_dcd_header
            t.n_analysed = read_dcd_header(t.analysed_path).get("n_frames")
        tls.append(t)
    path = os.path.join(out_dir, "frames.json")
    F.write(path, tls)
    return path, tls


def run_stages(spec, out_dir, names, stride=None, done=None):
    """Run the stages asked for, in increasing order of what they cost.

    The solvent stages read the SOLVATED trajectories -- 23-30 GB a replica
    against 1.7 GB stripped -- so they are never implied by `--only structure`
    and never run unless named.

    `done` is filled as each stage finishes, so a caller that catches a later
    stage's failure still knows which stages wrote their files.
    """
    from stages.solvent import composition
    from stages.structure import basic
    data = os.path.join(out_dir, "data")
    done = {} if done is None else done
    if "composition" in names and spec.has_solvent:
        done["composition"] = composition.run(spec, data)
    if "structure" in names:
        done["structure"] = basic.run(spec, data, stride=stride)
    if "ions" in names and spec.has_solvent:
        from stages.solvent import ions
        done["ions"] = ions.run(spec, data, stride=stride)
    # The water family shares ONE walk of the solvated trajectory: asking for
    # any of the three runs all three, because the second and third would
    # otherwise re-read 23 GB to compute what the first pass already had.
    if "interface_contacts" in names:
        from core import selections as _S2
        from stages.contacts import interface
        try:
            done["interface_contacts"] = interface.run(spec, data, stride=stride)
        except _S2.EmptySelection as exc:
            print(f"        interface_contacts: skipped -- {exc}")
    if "binding_energy" in names:
        # Reads results that already exist on disk -- no trajectory, so it
        # costs nothing and is safe to run every time. A system with no such
        # calculation produces no rows and no error.
        from stages.contacts import binding
        done["binding_energy"] = binding.run(spec, data)
    if "rdf" in names and spec.has_solvent:
        from stages.solvent import rdf
        done["rdf"] = rdf.run(spec, data, stride=stride)
    if "density_grids" in names and spec.has_solvent:
        from stages.solvent import grids
        # 1.0 A, not the stage's 0.5 A default. At 0.5 A a species costs
        # ~350 MB per replica -- about 12 GB for this folder -- and the bundle
        # is meant to be something a person can copy. One angstrom divides
        # that by eight and still resolves a binding site; the spacing is
        # written into every row, so nobody has to guess which was used.
        done["density_grids"] = grids.run(spec, data, stride=stride,
                                          spacing=float(os.environ.get(
                                              "GRID_SPACING", 1.0)))
    water_names = {"water", "water_shell", "water_bridges", "water_residence"}
    if (water_names & set(names)) and spec.has_solvent:
        from core import selections as _S
        from stages.solvent import water
        try:
            done.update(water.run(spec, data, stride=stride))
        except _S.EmptySelection as exc:
            # An apo system has no second side, so there is no INTERFACIAL
            # water to count -- the stage says so and refuses. That is right,
            # and it must not take the rest of the system's analyses with it:
            # a run over four systems where one is apo should still finish.
            print(f"        water: skipped -- {exc}")
    return done


def build_page(spec, out_dir, done, tls, doc, which="story"):
    """Build the page from EVERY table in data/, not just what has just run.

    Running one stage must not silently drop the others from the page: with
    `--only ions` the page would have been rebuilt from ion rows alone and the
    RMSD curves would have vanished, with nothing to say they had.
    """
    from core import io as _io
    from pages import report as P
    data = os.path.join(out_dir, "data")
    rows = _io.read(data) if os.path.isdir(data) else []
    if not rows:
        for r in done.values():
            rows += r.rows
    # Plotted subsets. Glycans are measured and kept in data/, but they are
    # not plotted here: at 12-20 A against the receptor's 4-6 A they flatten
    # the axis, and 735 sugar residues across 75 segments make an RMSF trace
    # nobody can read.
    with open(os.path.join(out_dir, "frames.json")) as fh:
        frames_doc = json.load(fh)
    samp = doc.get("sampling", {})
    # "3 replicas, 1512.0 ns in total, 200 ps per frame" told a reader three
    # things they have to multiply together. What they want is how long the
    # sampling is: 1.5 us, written with the sign, by the same rule the record
    # uses for its length tag so the two can never disagree.
    total = samp.get("total_ns")
    subtitle = (f"{samp.get('replicas', len(spec.mds))} replicas"
                + (f" \u00b7 {_io.human_ns(total)}" if total else "")
                + f" \u00b7 {spec.bundle}")
    title = doc.get("title") or spec.options.get("title") or spec.label
    topology = os.path.basename(spec.mds[0].topology)
    made = []
    if which in ("story", "both"):
        from pages import story as ST
        made.append(ST.build(os.path.join(out_dir, "page_story.html"),
                             title=title, subtitle=subtitle, rows=rows,
                             frames_doc=frames_doc, topology_name=topology,
                             doc=doc, bundle=spec.bundle,
                             components=spec.components))
    if which == "paper":
        # Asked for BY NAME only, never through `both`. Same protection the
        # working page has: `both` meant story + plotly before this page
        # existed, and quietly adding a third to it would rewrite a file
        # somebody did not ask to have rewritten.
        #
        # The same builder, the same rows, with the article's figures set
        # among them. A separate file, so publishing the combined page never
        # changes the plain one -- and so a bundle with no paper simply has
        # one page rather than a broken second one.
        from pages import paper as PAP
        from pages import story as ST
        art = PAP.read(out_dir)
        if not art:
            print("        no paper.yaml in this bundle; run import_paper.py "
                  "to attach one. Nothing was written.")
        else:
            made.append(ST.build(os.path.join(out_dir, "page_paper_story.html"),
                                 title=title, subtitle=subtitle, rows=rows,
                                 frames_doc=frames_doc, topology_name=topology,
                                 doc=doc, bundle=spec.bundle,
                                 components=spec.components, paper=art))
    if which in ("plotly", "both"):
        made.append(P.build(os.path.join(out_dir, "page_plotly.html"),
                            title=title, subtitle=subtitle, rows=rows,
                            frames_doc=frames_doc, topology_name=topology))
    return made


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--system", help="one record stem or identifier")
    ap.add_argument("--only", default="composition,structure",
                    help="comma-separated stage names")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--detect-only", action="store_true",
                    help="say what each system is, write nothing")
    # The working page is only rewritten when asked for. It is the one the
    # owner has been reading all week, and a new page must not be able to
    # replace it as a side effect of rebuilding.
    ap.add_argument("--page", choices=("story", "plotly", "paper", "both"),
                    default="story",
                    help="which page(s) to write (default: story only; "
                         "`paper` is by name only and needs paper.yaml from "
                         "import_paper.py; "
                         "page_plotly.html is left alone unless named)")
    ap.add_argument("--pages-only", action="store_true",
                    help="rebuild the page and the figures from the tables "
                         "already written, without touching a trajectory")
    a = ap.parse_args(argv)
    folder = os.path.abspath(a.folder)
    names = {s.strip() for s in a.only.split(",") if s.strip()}

    found_any = False
    failures = []
    for spec, found, doc in specs_from(folder, a.system):
        found_any = True
        comp = found["composition"]
        print(f"\n  {spec.options['record']}")
        print(f"    detected : {', '.join(spec.components)}  -> {spec.bundle}")
        print(f"    system   : {comp['solute_atoms']:,} solute atoms, "
              f"{comp['water_molecules']:,} waters, "
              f"{comp['ions_total']:,} ions, "
              f"{comp['glycan_residues']} glycan residues")
        for name in ("receptor", "ligand", "glycans"):
            if spec.selection(name):
                print(f"    {name:<9}: {spec.selection(name)[:88]}")
        print(f"    runs     : {', '.join(spec.stages(implemented_only=True))}")
        if spec.planned():
            print(f"    planned  : {', '.join(spec.planned())} "
                  f"(not written yet)")
        # Solute-only or not is whether a DIFFERENT topology was paired with
        # it, not whether the file happens to be called `merged.dcd`:
        # `OCCOPEN_1R_stripped.dcd` was being reported as solvated.
        print(f"    MDs      : " + ", ".join(
            f"{m.id}({'solute-only' if m.topology != m.solvated_topology else 'solvated'})"
            for m in spec.mds))
        if a.detect_only:
            continue

        out = bundle_dir(folder, spec.options["record"])
        os.makedirs(os.path.join(out, "figures_user"), exist_ok=True)
        _, tls = write_frames(spec, out)
        # NO TIME AXIS, NO ANALYSIS. Every curve is drawn against time and
        # every caption quotes a time; a run whose frame spacing nobody can
        # state produces tables full of t = 0 and a page that looks finished.
        # Refusing here costs a minute of reading; the alternative costs
        # whoever reads the page.
        unknown = [t.md_id for t in tls if not t.framestep_ns]
        if unknown:
            print(f"    SKIPPED  : {len(unknown)} replica(s) have no frame "
                  f"spacing ({', '.join(unknown[:4])}"
                  f"{', ...' if len(unknown) > 4 else ''}).")
            print(f"               The trajectory headers do not state one "
                  f"(post-processed copies say nsavc=1), and neither does the "
                  f"record.")
            print(f"               Declare it in {os.path.relpath(_settings(folder), folder)} "
                  f"as the dataset's frame interval and re-run "
                  f"`ndp_build.py`, then this folder.")
            continue
        if a.pages_only:
            from core import io as _io
            data = os.path.join(out, "data")
            rows = _io.read(data) if os.path.isdir(data) else []
            if not rows:
                print("    nothing to rebuild from: run the stages first")
                continue
            done = {"existing": _Rows(rows)}
            print(f"    tables   : {len(rows):,} rows read back from data/")
        else:
            # ONE SYSTEM'S FAILURE IS NOT THE FOLDER'S. B40T77 holds a record
            # whose topology has no protein in it at all -- a free aptamer,
            # all nucleic acid -- and the exception from that one ended the
            # whole run, discarding the two systems that had already built and
            # the four that had not been reached. A folder of twenty-five
            # records must not be all-or-nothing.
            partial = {}
            try:
                done = run_stages(spec, out, names, stride=a.stride, done=partial)
            except Exception as exc:                       # noqa: BLE001
                first = str(exc).splitlines()[0][:150] if str(exc) else ""
                print(f"    FAILED   : {type(exc).__name__}: {first}")
                # Stages before the failure HAVE written their files; saying
                # "nothing was written" left fresh and stale tables side by side
                # with no hint which were which.
                kept = ", ".join(partial) or "none"
                print(f"               the rest of this folder continues; "
                      f"stages that finished first and wrote data/: {kept}.")
                failures.append((spec.options.get("record", spec.id),
                                 f"{type(exc).__name__}: {first}"))
                continue
            for stage, rows in done.items():
                print(f"    {stage:<9}: {len(rows):,} rows")
        for page in build_page(spec, out, done, tls, doc, which=a.page):
            print(f"    page     : {os.path.relpath(page, folder)} "
                  f"({os.path.getsize(page) / 1e6:.1f} MB)")
        try:
            from figures import curves
            made = curves.all_for(out)
            print(f"    figures  : {len(made)} file(s) in figures/ "
                  f"({', '.join(sorted(set(f.rsplit('.', 1)[0] for f in made)))})")
        except Exception as exc:                           # noqa: BLE001
            print(f"    figures  : skipped ({type(exc).__name__}: {exc})")

    if not found_any:
        print(f"\n  No records found under {_results(folder)}.\n")
        return 1
    if failures:
        # Named, at the end, where a person reading a long run will see them.
        # A folder that half worked must say which half.
        print(f"\n  {len(failures)} system(s) failed and were skipped:")
        for record, why in failures:
            print(f"    {record:<34} {why}")
        print()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
