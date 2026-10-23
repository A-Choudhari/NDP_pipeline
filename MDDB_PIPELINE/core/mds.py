#!/usr/bin/env python3
"""
mds.py -- from indexed run directories to MDs: files in time order, the right
topology, and which runs are replicas of the same system.

Everything here used to live inside the MDDB exporter (mddb_lib/mwf.py), which
meant a metadata record could only be produced by running the deposition
workflow's assembler. It is generic: an MD is a topology plus its production
trajectory segments, whoever asks.

    production_files(...)   the production files of one run directory, in TIME
                            order, and how that order was established
    replica_split(...)      one directory holding several replicas side by side
    topology_by_atoms(...)  the topology whose atom count EXACTLY matches
    topology_groups(...)    runs whose topology files are byte-identical --
                            the same simulated system, whatever the paths say
    md_records(...)         all of the above, as one record per MD
"""

import hashlib
import os
import re

#: Topology extensions, in the order a tie is broken (see _topology_rank).
TOPO_ANY = (".psf", ".prmtop", ".parm7", ".top", ".tpr", ".gro", ".pdb")


def topology_atoms(path):
    """Atom count declared by a topology, read from its header -- never loads
    coordinates. PSF: the `!NATOM` line. prmtop: the first POINTERS integer."""
    try:
        if path.endswith((".prmtop", ".parm7")):
            with open(path, errors="replace") as fh:
                head = fh.read(20000)
            m = re.search(r"%FLAG POINTERS\s*\n%FORMAT\([^)]*\)\s*\n\s*(\d+)", head)
            return int(m.group(1)) if m else None
        if path.endswith(".psf"):
            with open(path, errors="replace") as fh:
                for ln in fh:
                    if "!NATOM" in ln:
                        return int(ln.split()[0])
    except (OSError, ValueError):
        pass
    return None


def _topology_rank(path):
    """Ordinary prmtop, then PSF, then a chamber prmtop last -- the same
    ranking the index uses, so both layers pick the same file."""
    try:
        from .trajscan import topology_rank
        return topology_rank(path)
    except Exception:                                        # noqa: BLE001
        return 5


def topology_by_atoms(run_dir, natom, max_up=2):
    """The topology whose atom count EXACTLY equals the trajectory's.

    Names mislead in both directions: `merged.psf` is the stripped solute
    (55,818 atoms) sitting next to a full-system `pr.dcd` (767,697), and
    `solvated.psf` is the pre-neutralisation box, a few thousand atoms short.
    In ABZYME/4-OPEN/CDR the staging directory's `solvated_ions.psf` is a
    DIFFERENT file from the one each replica ran; matching on the atom count,
    from the replica's own directory outwards, is what keeps them apart.
    """
    if not natom or natom != natom:
        return None
    natom = int(natom)
    d = run_dir
    for _ in range(max_up + 1):
        try:
            names = sorted(os.listdir(d))
        except OSError:
            names = []
        hits = [n for n in names if n.endswith(TOPO_ANY)
                and topology_atoms(os.path.join(d, n)) == natom]
        if hits:
            hits.sort(key=lambda n: _topology_rank(os.path.join(d, n)))
            return os.path.join(d, hits[0])
        d = os.path.dirname(d)
    return None


def _file_md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def topology_groups(topologies):
    """{path: 8-character group id} -- one group per distinct topology CONTENT.

    Two MDs in the same group were run from the same system file; two in
    different groups were not, whatever their directory names suggest. This is
    what makes "one record per simulated system" a measurement rather than a
    guess about folder layout.
    """
    by_hash, groups = {}, {}
    for path in sorted(set(p for p in topologies if p)):
        try:
            h = _file_md5(path)
        except OSError:
            h = path
        by_hash.setdefault(h, []).append(path)
    for h, paths in by_hash.items():
        for p in paths:
            groups[p] = str(h)[:8]
    return groups


def _natural_key(name):
    """Sort key that reads embedded integers as numbers, not as text.

    The EXTENSION IS KEPT OUT of the compared chunks. With it inside, `pr.dcd`
    keys as ['pr.dcd'] and `pr2.dcd` as ['pr', 2, '.dcd']; 'pr' sorts before
    'pr.dcd', so every numbered continuation went ahead of the unnumbered
    segment it continues -- the lab's restart convention (pr.dcd, pr.part1.dcd,
    pr.part2.dcd ...), and it put the FIRST segment LAST in 128 of 139 MDs.
    """
    stem, ext = os.path.splitext(str(name))
    return ([int(t) if t.isdigit() else t.lower()
             for t in re.split(r"(\d+)", stem)], ext.lower())


def dcd_start_step(path):
    """The step number of a DCD's first frame, or None. Header only."""
    try:
        from .dcdheader import read_dcd_header
        head = read_dcd_header(path)
    except Exception:                                        # noqa: BLE001
        return None
    if not isinstance(head, dict) or "error" in head:
        return None
    start = head.get("istart")
    return int(start) if isinstance(start, (int, float)) else None


def start_of(path):
    """Where a trajectory file begins on its run's time axis, or None.

    A DCD header's first step, or the first value of an AMBER NetCDF `time`
    variable. Continuation segments begin at increasing values; replicas all
    begin at the same one.
    """
    low = path.lower()
    if low.endswith(".dcd"):
        return dcd_start_step(path)
    if low.endswith((".nc", ".netcdf")):
        try:
            from scipy.io import netcdf_file
            with netcdf_file(path, "r", mmap=True) as nc:
                t = nc.variables.get("time")
                first = float(t[0]) if t is not None and t.shape[0] else None
                del t
            return first
        except Exception:                                    # noqa: BLE001
            return None
    return None


def production_files(root, rel_path, classify, traj_ext):
    """(filenames, evidence) -- the production files of one run, in TIME order.

    ORDER IS PART OF THE DATA: everything downstream concatenates these in the
    order given, so this list defines the time axis. A DCD header carries the
    step of its first frame, and ordering on that is the real chronology; only
    where no header can supply it do we fall back to a natural sort -- and then
    say so, rather than letting a caller assume the order was verified.

    Evidence: 'single' | 'header' | 'natural' | 'none'.
    """
    d = os.path.join(root, rel_path)
    try:
        names = sorted(os.listdir(d))
    except OSError:
        return [], "none"
    files = [n for n in names
             if n.endswith(traj_ext) and classify(n, d) == "production"]
    if len(files) < 2:
        return files, "single"
    starts = {n: dcd_start_step(os.path.join(d, n)) for n in files}
    if all(v is not None for v in starts.values()):
        if len(set(starts.values())) == len(files):
            return sorted(files, key=lambda n: starts[n]), "header"
        # TWO SEGMENTS BEGINNING AT THE SAME STEP. A window that died after a
        # frame or two, and the window that resumed at that same step: CDR
        # 3-YZ23/3R (pr.part10, 1 frame) and MD/1R (pr.part2, 1 frame). Their
        # order is still determinate -- at one save interval the short one's
        # frames are inside the long one's span, so either order gives the same
        # time axis -- and calling that "unverified" sent two finished 504 ns
        # replicas to review for nothing. It is only ambiguous when the tied
        # segments were saved at different intervals.
        sizes = {n: _dcd_shape(os.path.join(d, n)) for n in files}
        groups = {}
        for n in files:
            groups.setdefault(starts[n], []).append(n)
        nested = all(
            len({sizes[n][1] for n in g}) == 1 and all(sizes[n][0] for n in g)
            for g in groups.values() if len(g) > 1)
        if nested:
            return sorted(files, key=lambda n: (starts[n], sizes[n][0])), "header"
    return sorted(files, key=_natural_key), "natural"


def _dcd_shape(path):
    """(frames, save interval) of a DCD, or (None, None)."""
    try:
        from core.dcdheader import read_dcd_header
        h = read_dcd_header(path)
        if "error" in h:
            return None, None
        return h.get("n_frames"), h.get("nsavc")
    except Exception:                                      # noqa: BLE001
        return None, None


#: A replica label inside a FILE name: R1, 1R, rep1, replica_2, delimited by
#: `.`, `_` or `-`. Deliberately narrower than the index's REPLICA_RE, which
#: also reads `pr2` and `run3` as replicas: inside one directory
#: `pr.dcd, pr2.dcd, pr3.dcd` are restart segments of one run.
_FILE_REPLICA = re.compile(
    r"(?:^|[._-])(R(\d{1,2})|(\d{1,2})R|rep(?:lica)?[_-]?(\d{1,2}))(?=$|[._-])",
    re.I)


#: A bare number closing a file's stem: `spike_open_run_1.dcd` ...
#: `_6.dcd` (SPIKE_CLEAVED, 2026-09-30). Tried only when _FILE_REPLICA finds
#: nothing, and trusted only on POSITIVE evidence -- see _end_split -- because
#: `md_1.xtc, md_2.xtc` is just as often one run's restart segments.
_END_REPLICA = re.compile(r"[_-](\d{1,2})$")


def replica_split(root, rel_path, files, rules=()):
    """[(replica number, [file])] when a directory holds REPLICAS side by side.

    One directory is normally one MD, its production files the segments of one
    time axis. `ABZYME/1-FULL-GP120/MDs` breaks that: R1.dcd ... R5.dcd are five
    independent replicas, and concatenating them would publish one MD whose
    "time" jumps back to the start four times.

    Split only on evidence, all of it required:
      * every file carries exactly one replica label, and the labels differ;
      * with the label taken out, every name is the same;
      * the files do not continue one another: where a start can be read they
        must not all differ.

    When that finds no split, two narrower routes are tried, in order:
      * `rules` -- the author's `folder: replicas:` patterns (capture group =
        replica number). The author's word: splits even where no start can be
        read, but never where the starts show the files continue one another;
      * a bare `_<n>` closing every stem, which also needs every start READ and
        all of them EQUAL -- replicas begin together, restart segments do not.
    """
    if len(files) < 2:
        return None
    split = _label_split(root, rel_path, files)
    if split is None and rules:
        split = _rule_split(root, rel_path, files, rules)
    if split is None:
        split = _end_split(root, rel_path, files)
    return split


def _starts(root, rel_path, files):
    return [start_of(os.path.join(root, rel_path, n)) for n in files]


def _continues(starts):
    """True when every start was read and no two are the same: segments."""
    return all(v is not None for v in starts) and len(set(starts)) == len(starts)


def _by_number(files, number_of):
    """{number: file} when every file yields a distinct number and the names
    agree once it is taken out; else None."""
    by_rest = {}
    for name in files:
        got = number_of(name)
        if got is None:
            return None
        num, rest = got
        if num in by_rest.setdefault(rest, {}):
            return None
        by_rest[rest][num] = name
    if len(by_rest) != 1:
        return None
    (_, by_num), = by_rest.items()
    return by_num


def _rule_split(root, rel_path, files, rules):
    def number_of(name):
        for rx in rules:
            m = rx.search(name)
            if m and m.group(1) and m.group(1).isdigit():
                return int(m.group(1)), name[:m.start(1)] + "\0" + name[m.end(1):]
        return None
    by_num = _by_number(files, number_of)
    if by_num is None or _continues(_starts(root, rel_path, files)):
        return None
    return [(num, [by_num[num]]) for num in sorted(by_num)]


def _end_split(root, rel_path, files):
    def number_of(name):
        stem, ext = os.path.splitext(name)
        m = _END_REPLICA.search(stem)
        return (int(m.group(1)), stem[:m.start(1)] + "\0" + ext) if m else None
    by_num = _by_number(files, number_of)
    if by_num is None:
        return None
    starts = _starts(root, rel_path, files)
    if any(v is None for v in starts) or len(set(starts)) != 1:
        return None
    return [(num, [by_num[num]]) for num in sorted(by_num)]


def _label_split(root, rel_path, files):
    by_rest = {}
    for name in files:
        stem = os.path.splitext(name)[0]
        hits = list(_FILE_REPLICA.finditer(stem))
        if len(hits) != 1:
            return None
        m = hits[0]
        num = int(next(g for g in m.groups()[1:] if g))
        rest = stem[:m.start(1)] + "\0" + stem[m.end(1):]
        if num in by_rest.setdefault(rest, {}):
            return None
        by_rest[rest][num] = name
    if len(by_rest) != 1:
        return None
    starts = [start_of(os.path.join(root, rel_path, n)) for n in files]
    if all(v is not None for v in starts) and len(set(starts)) == len(starts):
        return None
    (_, by_num), = by_rest.items()
    return [(num, [by_num[num]]) for num in sorted(by_num)]


def md_records(root, runs, classify, traj_ext, md_ids, basis=None, replicas=()):
    """One record per MD, ready for the glossary.

    `runs` is the indexed run rows selected for a dataset; `md_ids` maps a run's
    relative path to its stable id. Returns (records, notes) where a record is:

        md_id, source_dir, trajectories (time order), ordering, topology,
        replica (when a directory was split), stage_assumed (files that are
        production only by default)
    """
    records, notes = [], []
    for r in runs.sort_values("rel_path").itertuples():
        rel = str(r.rel_path)
        run_dir = os.path.join(root, rel)
        files, how = production_files(root, rel, classify, traj_ext)
        if not files:
            notes.append(("mds", "note", f"{rel}: no file survives production "
                                         f"classification; not recorded"))
            continue
        natom = getattr(r, "natom", None)
        topology = topology_by_atoms(run_dir, natom)
        groups = replica_split(root, rel, files, rules=replicas)
        units = ([(f"{md_ids[rel]}_R{num}", group, "single", num)
                  for num, group in groups] if groups
                 else [(md_ids[rel], files, how, None)])
        for name, group, ordering, num in units:
            assumed = [f for f in group
                       if basis and basis(f, run_dir) == "default"]
            records.append({
                "md_id": name,
                "source_dir": rel,
                "trajectories": list(group),
                "ordering": ordering,
                "topology": (os.path.relpath(topology, root)
                             if topology else None),
                "replica": num,
                "stage_assumed": assumed,
            })
    return records, notes
