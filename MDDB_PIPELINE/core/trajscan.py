#!/usr/bin/env python3
"""
trajscan.py -- inventory MD replica-runs into runs.parquet, engine-agnostic.

Generalised from the ABZYME/Aptamers-specific scanner. Two rules carried over
because each was learned from a real bug:

1. FRAME SPACING DIFFERS BY STAGE. heating.dcd is written every 5,000 steps
   (10 ps/frame); pr.dcd every 50,000 (100 ps/frame). Taking one spacing per
   directory reported every 500 ns production run as 50 ns. Nanoseconds are
   therefore accumulated per file, from that file's own header.

2. EVIDENCE OR "unknown", NEVER A GUESS. The old scanner hardcoded this lab
   member's antibody and state vocabulary. Shipping that as a default would
   silently mislabel someone else's spike or APOBEC data, so only the
   near-universal replica patterns are inferred; everything else is "unknown"
   until the user declares a grammar in substrate.yaml.

Frame counts are only reported when they can be obtained CHEAPLY. `.mdcrd` is
ASCII with no header and `.xtc` needs an index scan, so those are null unless
--deep is passed. The n_frames_method column always says which happened.
"""

import os
import re
import time

import pandas as pd

from . import config
from .dcdheader import read_dcd_header

# Near-universal across engines and groups: 1R/R1/rep1/replica2/pr2/run3.
#
# Two defects fixed here, both of which produced wrong replica numbers rather
# than absent ones:
#
#   `3-R2` matched `(\d{1,2})-R(?:\d{1,2})?` and captured the PREFIX, so it was
#   recorded as replica 3 when the R-number says 2. The R-number is the replica;
#   the prefix is a sibling index. The R-number now wins, with a separate
#   alternative for a bare `3-R` that has none.
#
#   `pr2`, `PR3`, `prod4` matched nothing at all and fell through to `unknown`.
#   24 directories were affected.
REPLICA_RE = re.compile(
    r"^(?:"
    r"(?P<a>\d{1,2})R"                                # 1R, 2R
    r"|(?:\d{1,2}-)?R(?P<b>\d{1,2})"                  # R1, 1-R1, 3-R2
    r"|(?P<c>\d{1,2})-R"                              # 3-R, no R-number
    r"|rep(?:lica)?[_-]?(?P<d>\d{1,2})"               # rep1, replica_2
    r"|(?:pr|prod|production)[_-]?(?P<e>\d{1,2})"     # pr2, PR3, prod4
    r"|run[_-]?(?P<f>\d{1,2})"                        # run3
    r")$", re.I)

#: Directory names that denote a STAGE, not a replica. `pr`, `rep` and `md` are
#: where the trajectory lives, not which repeat it is. They were being written
#: into the same `replica` column as real indices -- 250 of 808 rows -- so a
#: consumer could not tell "replica 1" from "the directory called pr". They now
#: go in their own column.
BARE_REPLICA = {"pr", "prod", "production", "rep", "md"}


# --- which stage of a simulation a trajectory file belongs to ----------------
#
# THIS USED TO MATCH FILENAME PREFIXES, AND PREFIXES LIE IN BOTH DIRECTIONS.
#
#   `protein_sol_nacl_eq.mdcrd` starts with "pr" -- from *pr*otein -- and was
#   read as production. 468 entries across the exported MDDB catalogue were
#   equilibration files offered as the simulation itself. The whole LC3B
#   dataset ("157 production runs", 3.0 GB) is made of them: 157 of its
#   directories hold the production input MD1.in and exactly 2 hold its output.
#
#   `base_model.lipids_protein_glycans.R1.dcd` matched no prefix at all and
#   became "other", so its directory recorded n_production = 0 and disappeared.
#   That directory holds five 32 GB trajectories of a glycosylated membrane
#   system: 161 GB, the largest run directory in this workspace, invisible.
#
# So: match on TOKENS, not prefixes, and when nothing matches say `unresolved`
# rather than `other`. "Other" reads like a decision. It was an omission, and
# an omission that hides 161 GB has to be visible in the report.
_TOKEN_SPLIT = re.compile(r"[^a-z0-9]+")

#: Unambiguous. Checked before the stage vocabularies so `prod_npt` reads as
#: production instead of being captured by the ensemble name `npt`.
STRONG_PRODUCTION = {"prod", "production", "prd"}

#: Real stage names, but also common English fragments -- `pr` is the whole
#: reason `protein` was misread. Checked last, after every other vocabulary.
WEAK_PRODUCTION = {"pr", "md", "run", "sim"}

STAGE_TOKENS = (
    ("minimisation", {"min", "minim", "minimise", "minimize", "minimisation",
                      "minimization", "em", "emin", "steep"}),
    ("heating", {"heat", "heating", "anneal", "annealing", "warmup"}),
    # `npt`/`nvt` are ensemble names rather than stage names, but in every
    # corpus seen here they name the equilibration steps. STRONG_PRODUCTION
    # overrides them, for the case where production ran in NPT and said so.
    ("equilibration", {"eq", "equi", "equil", "equilibration", "equilibrate",
                       "npt", "nvt", "nve", "relax"}),
)

MERGED_TOKENS = {"merged", "merge", "cat", "concat", "combined", "joined"}

#: Substrings rather than tokens: strong enough signals to catch however the
#: name is punctuated. `not_water_not_ions` has to match here because its
#: tokens are the innocuous words "not", "water" and "ions".
STRIPPED_MARKERS = ("stripped", "strip", "nowat", "no_wat", "nowater",
                    "not_water", "notwater", "noions", "no_ions")

#: Directory names that mean "this is not simulation output": analysis,
#: scratch and helper trees. A file inside one is never production BY DEFAULT
#: -- ABZYME/4-OPEN/CDR/SCRIPTS/SMOKETEST holds `_mmpbsa_com_tmp.dcd`, and
#: `1-1F4/0-MD/1-R1/ANALYSIS/` holds a 10-frame cpptraj extract; both were read
#: as simulations by a rule that only looked at the file name.
NOT_SIMULATION = {"analysis", "analyses", "scripts", "script", "smoketest",
                  "mmpbsa", "mmgbsa", "test", "tests", "tmp", "scratch",
                  "docking", "render", "rendering", "figures", "movie",
                  "movies", "benchmark", "benchmarks"}

#: Superseded or derived output rather than the run itself.
SUPERSEDED_MARKERS = ("-old", "_old", ".old", "bak", "backup", "~",
                      "failed")    # failed_heating_nvt_2026-07-30/: a discarded attempt
DERIVED_MARKERS = ("frames", "stride", "strided", "skip", "thinned", "subset")


def _tokens(stem):
    """Words in a filename stem, with digit prefixes and suffixes also offered
    as bare words so `eq1`, `md1` and `npt2` match their vocabularies."""
    out = set()
    for part in _TOKEN_SPLIT.split(stem):
        if not part:
            continue
        out.add(part)
        alpha = part.rstrip("0123456789")
        if alpha:
            out.add(alpha)
        tail = part.lstrip("0123456789")
        if tail:
            out.add(tail)
    return out


def stage_of(fn, parent=None):
    """Which stage a trajectory file belongs to.

    One of: stripped | superseded | derived | merged | minimisation | heating |
    equilibration | production | unresolved.

    `parent` is the containing directory's name, when the caller knows it. THE
    DIRECTORY IS EVIDENCE THE FILENAME DOES NOT CARRY, and ignoring it is wrong
    in a way that matters: `.../SASA-384-401/stripped/` holds sixteen files
    named `1R.pr.dcd`, `2R.pr.dcd` and so on -- names that say production --
    which are in fact 45,210-atom solute-only extracts of 788,461-atom runs.
    Only the directory says so. A derivative deposited as the simulation is
    precisely the error this classifier exists to prevent, so a parent naming a
    derivative overrides the filename. The weaker stage names (`pr`, `md`) never
    override anything, because a directory called `pr` says where a trajectory
    lives, not what it is.

    `unresolved` is not a synonym for "not production". It means nothing
    carried evidence either way, and it is counted and reported rather than
    dropped -- silently dropping these is what made 161 GB of real MD vanish.
    """
    if parent:
        p = os.path.basename(str(parent).rstrip("/\\")).lower()
        if any(m in p for m in STRIPPED_MARKERS):
            return "stripped"
        if any(m in p for m in SUPERSEDED_MARKERS):
            return "superseded"
        if any(m in p for m in DERIVED_MARKERS):
            return "derived"

    b = str(fn).lower()
    stem = b.rsplit(".", 1)[0] if "." in b else b

    if any(m in b for m in STRIPPED_MARKERS):
        return "stripped"
    if any(m in b for m in SUPERSEDED_MARKERS):
        return "superseded"
    if any(m in b for m in DERIVED_MARKERS):
        return "derived"

    tok = _tokens(stem)
    if tok & MERGED_TOKENS:
        return "merged"
    if tok & STRONG_PRODUCTION:
        return "production"
    for stage, vocab in STAGE_TOKENS:
        if tok & vocab:
            return stage
    if tok & WEAK_PRODUCTION:
        return "production"
    return "unresolved"


#: Formats that are PRODUCTION BY DEFAULT when their name says nothing.
#:
#: Decided 2026-09-21: a .dcd is the simulation until its author says otherwise.
#: `unresolved` was the honest answer, but it had a cost nobody could see from
#: the glossary: a folder whose trajectories are named for their contents rather
#: than their stage (`base_model.lipids_protein_glycans.R1.dcd` ... `R5.dcd`,
#: 161 GB) produced no MD at all, and an empty glossary says "nothing here".
#:
#: Two limits keep the default from overruling evidence:
#:   * it applies only where the name AND the directory are silent. `eq.dcd`,
#:     `stripped/1R.pr.dcd` and `pr6900frames.dcd` keep what their names say,
#:     and nothing inside `07_EQUILIBRATION_1/`, `ANALYSIS/` or `SCRIPTS/`
#:     becomes production by default;
#:   * it applies only in a directory with NO other production file. Where
#:     pr.dcd ... pr7.dcd are already production, a silent `4092-5091.dcd` beside
#:     them is almost always a frame-range extract of them, and appending it would
#:     put the same frames on the time axis twice.
#: Everything counted this way is marked (basis `default`) and becomes `review`
#: in the MD glossary until the author confirms it under `stages:`.
#: Other formats keep `unresolved`: `250ns.mdcrd` beside MD1.in is as often an
#: equilibration output as a production one.
DEFAULT_PRODUCTION_EXT = (".dcd",)


def classify_dir(dirpath, names, stage_rules=(), relpath=None, natoms=None):
    """{name: (stage, basis)} for the trajectory files of ONE directory.

    basis is how the stage was established, strongest first:

      author    a `stages:` rule in substrate.yaml matched the file's path
      name      stage_of() read it from the file name or its directory
      atoms     named (or defaulted) production, but holds FEWER atoms than the
                directory's largest production file: a solute-only subset.
                Needs `natoms` ({name: atom count}, from the headers)
      default   a silent .dcd in a directory with no other production file

    `relpath` is the directory relative to root; rules match `relpath/name`.
    Only a directory's own files are compared, which is why this takes the
    whole listing rather than one name: whether the default applies depends
    on the siblings.
    """
    base = str(relpath or "").replace(os.sep, "/").strip("/")
    out = {}
    for name in names:
        path = f"{base}/{name}" if base else name
        declared = next((stage for stage, pats in stage_rules
                         if any(p.search(path) for p in pats)), None)
        out[name] = ((declared, "author") if declared
                     else (stage_of(name, parent=dirpath), "name"))
    # The directory's own name is evidence too: `07_EQUILIBRATION_1/` or
    # `heating/` says what a silent file inside it is not.
    dtok = _tokens(os.path.basename(str(dirpath).rstrip("/\\")).lower())
    dir_says_other = (bool(dtok & MERGED_TOKENS) or bool(dtok & NOT_SIMULATION)
                      or any(dtok & vocab for _, vocab in STAGE_TOKENS))
    silent = [] if dir_says_other else [
        n for n, (st, _) in out.items()
        if st == "unresolved" and n.lower().endswith(DEFAULT_PRODUCTION_EXT)]

    # A SUBSET IS NOT THE RUN. A production run writes the whole system; a file
    # beside it with fewer atoms was cut from it. ABZYME/3-OCCLUDED-OPEN/MDs/
    # 2-REPLICA/1R holds md.dcd (45,210 atoms, a solute-only extract, written by
    # a tool) beside pr.dcd and pr_r1.dcd (711,846). `md` reads as production,
    # its header starts at step 1, so it was ordered FIRST and three published
    # MDs began with 1,105 frames of a different atom set. Only the author can
    # overrule this, under `stages:`.
    if natoms:
        cands = [n for n, (st, _) in out.items() if st == "production"] + silent
        known = {n: natoms[n] for n in cands if natoms.get(n)}
        if len(set(known.values())) > 1:
            top = max(known.values())
            for n, count in known.items():
                if count < top and out[n][1] != "author":
                    out[n] = ("stripped", "atoms")
    if not any(st == "production" for st, _ in out.values()):
        for name in silent:
            if out[name][0] == "unresolved":
                out[name] = ("production", "default")
    return out


def engine_of(filenames):
    exts = {os.path.splitext(f)[1].lower() for f in filenames}
    if ".dcd" in exts:
        return "NAMD/CHARMM"
    if exts & {".nc", ".netcdf", ".mdcrd", ".crd"}:
        return "AMBER"
    if exts & {".xtc", ".trr"}:
        return "GROMACS"
    return "unknown"


def sniff_format(path):
    """A trajectory's real format from its first bytes: netcdf | dcd | ascii | None.

    EXTENSIONS LIE. Both LC3B production trajectories are named `MD1.mdcrd` and
    are NetCDF (magic `CDF\\x02`). Tools that trust the extension read them as
    ASCII and fail -- MDAnalysis with a UnicodeDecodeError, and mwf maps `mdcrd`
    to its ASCII `crd` format (constants.py:305-339). Eight bytes settle it.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError:
        return None
    if head[:3] == b"CDF" or head[:4] == b"\x89HDF":
        return "netcdf"
    if len(head) >= 8 and head[4:8] == b"CORD":
        return "dcd"
    try:
        head.decode("ascii")
        return "ascii" if head else None
    except UnicodeDecodeError:
        return None


#: Below this, a stated frame spacing is not a save interval. Production
#: output in this workspace is saved every 50.4 ps at the finest; one
#: picosecond is two orders of magnitude below that and still far above any
#: integration step, so nothing real is excluded and every placeholder is.
MIN_CREDIBLE_SPACING_PS = 1.0


def dcd_spacing_ps(header):
    """Frame spacing a DCD header states, or None when it states a placeholder.

    A DCD written by a trajectory TOOL rather than by the engine does not know
    the original save interval and writes NSAVC = 1 -- one frame per step --
    with whatever DELTA it has to hand. Two of them are in this workspace and
    they do not look alike:

      VMD's DCD plugin        NSAVC=1, DELTA=1.0 AKMA  -> 0.0489 ps/frame
      cpptraj (`merged.dcd`)  NSAVC=1, DELTA=0.0205    -> 0.001  ps/frame

    Read literally the first turns 1-FULL-GP120's five 10,500-frame replicas
    into 0.5 ns each; the second turns 8-IgG/1-IF4's 3R, a 1,008 ns run, into
    5 picoseconds. Neither number has anything behind it, and both would have
    gone into the glossary and every record built on it.

    Matching VMD's exact pair caught the first and not the second, so the rule
    is the one the reasoning always rested on: NSAVC = 1 means the writer was
    not the engine, and a spacing finer than a picosecond is not a save
    interval. No all-atom engine steps 48.9 fs, and none saves every 1 fs step
    for a microsecond. Left unknown, which the owner can then DECLARE -- and
    `framestep_basis` says which of the two it was.
    """
    if header.get("nsavc") != 1:
        return header.get("dt_ps_header")
    if header.get("delta_akma") == 1.0:
        return None
    dt = header.get("dt_ps_header")
    if dt is not None and dt < MIN_CREDIBLE_SPACING_PS:
        return None
    return dt


def frame_info(path, deep=False):
    """Return (n_frames, natom, dt_ps, method). Never raises.

    method: 'header'      exact, read from the file header (cheap)
            'scan'        exact, obtained by scanning (only with --deep)
            'unavailable' the format has no cheap frame count
            'error:...'   the file could not be read
    """
    ext = os.path.splitext(path)[1].lower()
    # An AMBER-named file that is really NetCDF gets the NetCDF reader: its
    # header has the frame count and its time variable has the spacing.
    if ext in (".mdcrd", ".crd", ".trj") and sniff_format(path) == "netcdf":
        ext = ".nc"

    if ext == ".dcd":
        h = read_dcd_header(path)
        if "error" in h:
            return None, None, None, f"error:{h['error']}"
        return h["n_frames"], h["natom"], dcd_spacing_ps(h), "header"

    if ext in (".nc", ".netcdf"):
        # AMBER NetCDF. `frame` is the UNLIMITED dimension, and scipy reports an
        # unlimited dimension as None -- so `dimensions.get("frame")` was always
        # None here and every .nc run recorded a null frame count while
        # n_frames_method still claimed "header". That contradicted the contract
        # in AGENTS.md ("header is exact and cheap"), and the suite could not see
        # it because its invariant only checks unavailable => null, never
        # header => not-null. The real length is the coordinates array's first
        # axis, which costs nothing to read: NetCDF stores the shape, not the data.
        #
        # The sampling interval genuinely is absent from the header -- AMBER
        # records the time UNIT but never dt or ntwx -- so it is read from the
        # first two values of the `time` variable instead. That is two floats,
        # not a scan, and it is the difference between "sampling not readable
        # from headers" for every AMBER dataset and an exact answer.
        #
        # MEMORY-MAPPED, NOT READ. With mmap=False scipy copies every variable
        # into memory when the file is opened -- the whole coordinate array --
        # so "reading the header" of a 342 MB trajectory peaked at 737 MB, and
        # a colleague's 20 GB trajectory would have needed about 40 GB per file
        # on a shared machine. Mapped, the same call peaks at 70 MB, returns
        # the same values and is five times faster. A mapped file only closes
        # cleanly once nothing refers to its data, so the two variables are
        # released before the `with` block ends and plain numbers are returned.
        # NetCDF4 is HDF5 underneath, and scipy reads only the classic format.
        # Read its shape through h5py when that is installed; otherwise say what
        # is missing instead of a bare `error:TypeError`.
        try:
            with open(path, "rb") as fh:
                if fh.read(4) == b"\x89HDF":
                    return _netcdf4_header(path)
        except OSError:
            pass
        try:
            from scipy.io import netcdf_file
        except ImportError:
            return (None, None, None,
                    "error:scipy is not installed; it is needed to read "
                    "NetCDF headers")
        try:
            n = natom = dt_ps = None
            with netcdf_file(path, "r", mmap=True) as nc:
                natom = nc.dimensions.get("atom")
                n = nc.dimensions.get("frame")
                coords = nc.variables.get("coordinates")
                if n is None and coords is not None:
                    n = int(coords.shape[0])
                times = nc.variables.get("time")
                if times is not None and getattr(times, "shape", (0,))[0] >= 2:
                    try:
                        step = float(times[1]) - float(times[0])
                        # AMBER writes picoseconds and says so in time.units;
                        # anything else is not something to silently rescale.
                        units = getattr(times, "units", b"")
                        if isinstance(units, bytes):
                            units = units.decode("ascii", "replace")
                        if step > 0 and str(units).lower().startswith("pico"):
                            dt_ps = step
                    except (IndexError, TypeError, ValueError):
                        dt_ps = None
                del coords, times
            return (int(n) if n else None,
                    int(natom) if natom else None, dt_ps, "header")
        except Exception as exc:                             # noqa: BLE001
            return None, None, None, f"error:{type(exc).__name__}"

    if ext in (".xtc", ".trr", ".mdcrd", ".crd"):
        if not deep:
            # No cheap total. Scanning an .xtc or an ASCII .mdcrd means reading
            # the whole file; on 6 TB of data that is hours for a number nobody
            # asked for yet.
            return None, None, None, "unavailable"
        try:
            import mdtraj.formats as F
            cls = {".xtc": F.XTCTrajectoryFile, ".trr": F.TRRTrajectoryFile}.get(ext)
            if cls is None:
                return None, None, None, "unavailable"
            with cls(path) as t:
                return len(t), None, None, "scan"
        except Exception as exc:                             # noqa: BLE001
            return None, None, None, f"error:{type(exc).__name__}"

    return None, None, None, "unavailable"


def _netcdf4_header(path):
    """(n_frames, natom, dt_ps, method) for a NetCDF4/HDF5 trajectory."""
    try:
        import h5py
    except ImportError:
        return (None, None, None,
                "error:NetCDF4 (HDF5) file; install h5py to read its header")
    try:
        with h5py.File(path, "r") as f:
            coords = f.get("coordinates")
            if coords is None:
                return None, None, None, "error:no coordinates variable"
            n, natom = int(coords.shape[0]), int(coords.shape[1])
            dt_ps = None
            times = f.get("time")
            if times is not None and times.shape and times.shape[0] >= 2:
                units = times.attrs.get("units", b"")
                if isinstance(units, bytes):
                    units = units.decode("ascii", "replace")
                step = float(times[1]) - float(times[0])
                if step > 0 and str(units).lower().startswith("pico"):
                    dt_ps = step
            return n or None, natom or None, dt_ps, "header"
    except Exception as exc:                                 # noqa: BLE001
        return None, None, None, f"error:{type(exc).__name__}"


def _first_label(rules, relpath):
    """First declared label whose pattern matches, else 'unknown'.

    Declaration order in substrate.yaml decides ties, so the result is
    predictable and re-orderable by the person who wrote the rules. No match is
    reported as `unknown` -- absence, not a guess.
    """
    for label, pats in rules:
        for pat in pats:
            if pat.search(relpath):
                return label
    return "unknown"


def _captured_label(rules, relpath):
    """Label taken from the path by the first capturing rule that matches."""
    for rx in rules:
        m = rx.search(relpath)
        if m:
            for g in m.groups():
                if g:
                    return g
    return "unknown"


def parse_path(parts, relpath=None, system_rules=(), condition_rules=(),
               system_captures=(), condition_captures=()):
    """Structure extraction from a run's path.

    `replica` is inferred generically (1R/2R/3R, rep1, ...) because replica
    naming is near-universal. `system` and `condition` are NOT inferred: they
    come only from rules the user declared in substrate.yaml, because what
    counts as a "system" is a modelling decision, not a string pattern. With no
    rules declared both stay `unknown`, which is what every build did before
    this was configurable.
    """
    # Two separate questions, previously answered into one column: WHICH REPEAT
    # is this (a replica index), and WHAT IS THE DIRECTORY CALLED (a stage name
    # like pr/rep/md). Mixing them meant `rep` and `2R` sat side by side as if
    # they were the same kind of value, and it is why `1R` (39 rows) was so much
    # rarer than `2R` (169): the first replica is usually the directory named
    # `pr` or `rep`, so it never got an index at all.
    #
    # Both branches now take the FIRST match walking down the path. The bare
    # branch used to have no `break` while the regex branch did, so a regex hit
    # took the first component and a bare hit took the last -- for
    # `.../Top5/1/pr/md/PR4` that discarded `PR4`, the one component naming the
    # replica, and recorded `md`.
    replica_index = stage_dir = None
    for p in parts:
        if replica_index is None:
            m = REPLICA_RE.match(p)
            if m:
                num = next((g for g in m.groups() if g), None)
                if num:
                    replica_index = f"{int(num)}R"
        if stage_dir is None and p.lower() in BARE_REPLICA:
            stage_dir = p.lower()

    # `replica` keeps its old meaning and its old values so every existing
    # consumer and ground-truth check still works; the two new columns are what
    # you read when you need to know which kind of answer you got.
    replica = replica_index or stage_dir or "unknown"

    if relpath is None:
        relpath = os.sep.join(parts)
    # Match on "/"-separated form so a rule written with "/" works regardless
    # of the platform separator.
    norm = relpath.replace(os.sep, "/")

    # An explicit `systems:` label wins over a captured one: it is the more
    # specific statement, and it is how you override a capture rule for the one
    # directory whose name lies.
    system = _first_label(system_rules, norm)
    if system == "unknown":
        system = _captured_label(system_captures, norm)
    condition = _first_label(condition_rules, norm)
    if condition == "unknown":
        condition = _captured_label(condition_captures, norm)
    return {"replica": replica, "replica_index": replica_index,
            "stage_dir": stage_dir, "system": system, "condition": condition}


#: Topology formats that state their own atom count cheaply, in a header.
TOPO_COUNTED = (".psf", ".prmtop", ".parm7")


def topology_atoms(path):
    """Atom count a topology declares, from its header. Never loads coordinates.

    PSF: the `!NATOM` line. prmtop/parm7: the first POINTERS integer.
    Returns None when the format does not say, or the file cannot be read.
    """
    try:
        low = path.lower()
        if low.endswith((".prmtop", ".parm7")):
            with open(path, errors="replace") as fh:
                head = fh.read(20000)
            m = re.search(r"%FLAG POINTERS\s*\n%FORMAT\([^)]*\)\s*\n\s*(\d+)",
                          head)
            return int(m.group(1)) if m else None
        if low.endswith(".psf"):
            with open(path, errors="replace") as fh:
                for line in fh:
                    if "!NATOM" in line:
                        return int(line.split()[0])
    except (OSError, ValueError, IndexError):
        pass
    return None


def is_chamber_prmtop(path):
    """True for a CHAMBER-style prmtop: CHARMM parameters written in AMBER format.

    ParmEd/chamber files carry `%FLAG CTITLE` where an ordinary prmtop has
    `%FLAG TITLE`. They describe the same atoms, but not every tool reads them:
    MDAnalysis refuses them outright ("Chamber-style TOP file ... not
    supported"), which is how this was found -- by handing a glossary entry to
    an analysis. They are also evidence about the science: the run used a CHARMM
    force field even though AMBER wrote the trajectory.
    """
    if not str(path).lower().endswith((".prmtop", ".parm7")):
        return False
    try:
        with open(path, errors="replace") as fh:
            return "%FLAG CTITLE" in fh.read(4000)
    except OSError:
        return False


def topology_rank(path):
    """Sort key among topologies that all match the trajectory's atom count.

    An ordinary prmtop first (engine-native, readable everywhere), then a PSF,
    and a chamber prmtop last -- chosen only when nothing else fits. All fifteen
    DESIGN-V3 runs were given a chamber `cx.prmtop` while `com.psf` describing
    the same 15,997 atoms sat beside it.
    """
    low = str(path).lower()
    if low.endswith((".prmtop", ".parm7")):
        kind = 2 if is_chamber_prmtop(path) else 0
    elif low.endswith(".psf"):
        kind = 1
    else:
        kind = 3
    return (kind, len(os.path.basename(low)))


def find_topology(dirpath, max_up=2, natom=None):
    """The topology describing this run. Returns (path, kind, atom_match).

    `atom_match` is True when the chosen topology's declared atom count equals
    the trajectory's, False when it demonstrably does not, and None when the
    question could not be asked -- so a mismatch becomes queryable instead of
    invisible.

    WHY THE NAME HEURISTIC WAS NOT ENOUGH. This used to prefer a filename
    containing "merged", then the shortest name. Both mislead, and measurably:

      - All three aptasuite-models runs were given `merged.psf`, which declares
        46,633 atoms, for a trajectory of 637,634. A 13.7x mismatch, recorded as
        if it were the answer.
      - Six of fifteen DESIGN-V3 runs were given `lig.psf` -- a 694-atom
        LIGAND-ONLY topology -- while the matching `com.psf` sat in the same
        directory. `com.psf` and `lig.psf` are both 7 characters, and `sort` is
        stable, so which one won came down to os.listdir order. Identical
        directories got different topologies.

    An exact atom-count match at this scale is a strong identity signal, so it
    is tried first and the name heuristic is kept only as the fallback for
    formats that declare no count (.pdb, .crd, .gro, .tpr).
    """
    d = dirpath
    for _ in range(max_up + 1):
        try:
            names = sorted(os.listdir(d))
        except OSError:
            return None, None, None
        # Symlinks are skipped for the same reason the trajectory scan skips
        # them (config.skip_symlink): a link can point into another user's tree,
        # and the index should describe this folder, not follow it elsewhere.
        # find_topology was the one place missing this guard, and the invariant
        # that was supposed to catch it greps the source for the string
        # "skip_symlink" rather than testing behaviour, so it passed anyway.
        cands = [f for f in names
                 if f.endswith(config.TOPO_EXT)
                 and not config.skip_symlink(os.path.join(d, f))]
        if cands:
            if natom:
                exact = [f for f in cands
                         if f.lower().endswith(TOPO_COUNTED)
                         and topology_atoms(os.path.join(d, f)) == int(natom)]
                if exact:
                    exact.sort(key=lambda f: topology_rank(os.path.join(d, f)))
                    return (os.path.join(d, exact[0]),
                            exact[0].rsplit(".", 1)[-1], True)
            real = [f for f in cands if not f.endswith(".pdb")] or cands
            real.sort(key=lambda f: (0 if "merged" in f.lower() else 1, len(f)))
            pick = os.path.join(d, real[0])
            declared = topology_atoms(pick) if real[0].lower().endswith(
                TOPO_COUNTED) else None
            match = None if (declared is None or not natom) \
                else (declared == int(natom))
            return pick, real[0].rsplit(".", 1)[-1], match
        d = os.path.dirname(d)
    return None, None, None


def scan(cfg, verbose=True):
    rows = []
    deep = cfg.get("deep_scan", False)
    # Compile the user's label rules ONCE, before the walk. A bad regex should
    # stop the build immediately with the offending pattern named, not fail
    # 8 TB into a scan.
    system_rules = config.compile_label_rules(cfg, "systems")
    condition_rules = config.compile_label_rules(cfg, "conditions")
    system_captures = config.compile_capture_rules(cfg, "systems_from")
    condition_captures = config.compile_capture_rules(cfg, "conditions_from")
    stage_rules = config.compile_stage_rules(cfg)
    if verbose and (system_rules or condition_rules
                    or system_captures or condition_captures):
        print(f"    label rules: {len(system_rules)} system + "
              f"{len(system_captures)} captured, "
              f"{len(condition_rules)} condition + "
              f"{len(condition_captures)} captured", flush=True)
    # Declared corpora are DIRECTORIES under root, so trajectories lying loose
    # in root itself were never scanned: ABZYME/1-FULL-GP120/MDs keeps its five
    # replicas at the top, beside the analysis folders that became its corpora,
    # and its index held every derivative and not the simulation. Root's own
    # files are scanned too, without descending (the corpora cover the rest).
    roots = [(name, base, None) for name, base in config.corpus_roots(cfg)]
    top = os.path.realpath(cfg["root"])
    if cfg.get("corpora") and all(os.path.realpath(b) != top for _, b, _ in roots):
        roots.append((os.path.basename(top) or "root", cfg["root"], 0))
    for corpus_name, base, depth in roots:
        if verbose:
            print(f"    corpus {corpus_name}: {base}"
                  + (" (its own files only)" if depth == 0 else ""), flush=True)
        t0, n_dirs = time.time(), 0
        for dirpath, _dirnames, filenames in config.walk(base, cfg["prune"],
                                     max_depth=depth,
                                     exclude=cfg.get("exclude_paths", ())):
            trajs = [f for f in filenames
                     if f.endswith(config.TRAJ_EXT)
                     and not config.skip_symlink(os.path.join(dirpath, f))]
            if not trajs:
                continue
            n_dirs += 1
            relp = config.rel(dirpath, cfg)
            meta = parse_path(relp.split(os.sep), relp,
                              system_rules, condition_rules,
                              system_captures, condition_captures)

            # Headers first: whether a file is a solute-only subset of its
            # neighbours is a question about atom counts, which only they give.
            info = {}
            for fn in sorted(trajs):
                full = os.path.join(dirpath, fn)
                try:
                    sz = os.path.getsize(full)
                except OSError:
                    continue
                info[fn] = (sz,) + tuple(frame_info(full, deep=deep))
            staged = classify_dir(dirpath, list(info), stage_rules, relp,
                                  natoms={fn: v[2] for fn, v in info.items()})
            by_stage, errs = {}, []
            n_default = n_author = n_atoms = 0
            prod_natoms = set()
            n_prod = n_merged = 0
            prod_ns = 0.0
            bytes_total = bytes_prod = 0
            natom = dt_prod = natom_source = None
            natom_bytes = 0
            natom_conflict = False
            dt_mixed = False
            methods = set()

            for fn in sorted(info):
                sz, nf, na, dtp, method = info[fn]
                st, basis = staged[fn]
                by_stage[st] = by_stage.get(st, 0) + 1
                if st == "production":
                    n_default += basis == "default"
                n_author += basis == "author"
                n_atoms += basis == "atoms"
                bytes_total += sz
                if st == "production":
                    bytes_prod += sz

                methods.add(method.split(":")[0])
                if method.startswith("error"):
                    errs.append(f"{fn}: {method}")
                    continue
                # ATOM COUNT MUST COME FROM A PRODUCTION FILE. This used to take
                # whichever file sorted first and yielded a count, so
                # ABZYME/3-OCCLUDED-OPEN/MDs/1R recorded 711,846 atoms read from
                # heating.dcd while its production md.dcd has 45,210. That is
                # not a cosmetic error: the deposition layer feeds this number
                # into its exact-atom-count topology matcher, so a wrong natom
                # silently selects a topology describing a different molecule.
                # A non-production count is still better than nothing, so it is
                # kept as a fallback and overwritten the moment a production
                # file supplies one.
                if na:
                    if st == "production":
                        # Among production files, the LARGEST one. Not the first
                        # alphabetically: ABZYME/3-OCCLUDED-OPEN/MDs/1R holds
                        # both pr.dcd (711,846 atoms, 45 GB) and md.dcd (45,210
                        # atoms, a solute-only derivative that the filename does
                        # not admit to being), and `md` sorts first. Size is the
                        # one signal that reliably picks the primary trajectory.
                        if natom_source != "production" or sz > natom_bytes:
                            natom, natom_source, natom_bytes = na, "production", sz
                        # Every count, not a comparison with the running one:
                        # that missed a smaller file sorted BEFORE the larger.
                        prod_natoms.add(na)
                        natom_conflict = len(prod_natoms) > 1
                    elif natom is None:
                        natom, natom_source = na, st
                if nf is None:
                    continue
                if st == "production":
                    n_prod += nf
                    if dtp:
                        prod_ns += nf * dtp / 1000.0
                        if dt_prod is None:
                            dt_prod = dtp
                        elif abs(dt_prod - dtp) > 1e-6:
                            dt_mixed = True
                elif st == "merged":
                    n_merged += nf

            topo, kind, topo_match = find_topology(dirpath, natom=natom)
            # A spacing that is only true of the first file is not a property of
            # the directory. MDs/1R published dt_ps_production = 4.9e-05 ps next
            # to prod_ns_header = 1008.1 ns, because the first production file
            # was a post-processed md.dcd whose header says NSAVC=1. The total
            # is still right (nanoseconds are accumulated per file, from each
            # file's own header) but the single dt is not, so it is withheld.
            # dt_note already explains why, and mwf's per-dataset framestep
            # drops nulls -- which is what it should have been doing all along.
            dt_out = None if dt_mixed else dt_prod
            # `round(x, 2) if x else None` stored a real 0.0 for a tiny nonzero
            # total, making "measured as zero" indistinguishable from "could not
            # be read". Below the rounding resolution we have no number to give.
            ns_out = round(prod_ns, 2) if prod_ns else None
            if prod_ns and ns_out == 0.0:
                ns_out = None
            rows.append({
                "run_id": relp.replace(os.sep, "__"),
                "owner": cfg["owner"],
                "corpus": corpus_name,
                "rel_path": relp,
                "replica": meta["replica"],
                # The two halves of what `replica` used to conflate: which
                # repeat this is, and what the directory is called.
                "replica_index": meta.get("replica_index"),
                "stage_dir": meta.get("stage_dir"),
                "system": meta["system"],
                "condition": meta["condition"],
                "engine": engine_of(trajs),
                "n_traj_files": len(trajs),
                "traj_gb": round(bytes_total / 1e9, 3),
                "traj_gb_production": round(bytes_prod / 1e9, 3),
                "n_frames_production": n_prod or None,
                "n_frames_merged": n_merged or None,
                "prod_ns_header": ns_out,
                "dt_ps_production": dt_out,
                "dt_mixed": dt_mixed,
                "natom": natom,
                # Which stage the atom count came from. `production` is the only
                # value the deposition layer should trust for topology matching.
                "natom_source": natom_source,
                # Production files in this directory disagree about how many
                # atoms they contain, which means one of them is a stripped or
                # solute-only derivative whose filename does not say so. `natom`
                # is the largest file's; treat the directory with suspicion.
                "natom_conflict": natom_conflict,
                # How the frame numbers above were obtained. Never silent.
                "n_frames_method": ("mixed" if len(methods - {"error"}) > 1
                                    else (sorted(methods - {"error"}) or ["none"])[0]),
                "n_minimisation": by_stage.get("minimisation", 0),
                "n_heating": by_stage.get("heating", 0),
                "n_equilibration": by_stage.get("equilibration", 0),
                "n_production": by_stage.get("production", 0),
                "n_merged_files": by_stage.get("merged", 0),
                "n_stripped": by_stage.get("stripped", 0),
                # Superseded (`pr6-old.dcd`) and derived (`pr6900frames.dcd`)
                # used to be counted as production; `unresolved` is the file
                # whose name says nothing, which used to be silently "other".
                "n_superseded": by_stage.get("superseded", 0),
                "n_derived": by_stage.get("derived", 0),
                "n_unresolved": by_stage.get("unresolved", 0),
                # Of n_production, how many are production only BY DEFAULT: a
                # .dcd whose name says nothing, in a directory with nothing else
                # that is production. Not established -- assumed, and flagged
                # until the author confirms it under `stages:`.
                "n_production_default": n_default,
                # Files whose stage the author declared under `stages:`.
                "n_stage_declared": n_author,
                # Files named as production that hold FEWER atoms than the
                # directory's largest production file: solute-only subsets whose
                # names do not say so (md.dcd, 45,210 atoms, beside pr.dcd,
                # 711,846). Counted in n_stripped, not n_production.
                "n_stripped_by_atoms": n_atoms,
                # Kept so existing consumers still read a number; it now means
                # "matched no stage vocabulary", which is what it always
                # implied. n_unresolved is the one to look at.
                "n_other": (by_stage.get("other", 0)
                            + by_stage.get("unresolved", 0)
                            + by_stage.get("superseded", 0)
                            + by_stage.get("derived", 0)),
                "topology_path": config.rel(topo, cfg) if topo else None,
                "topology_kind": kind,
                # True/False/None: does the chosen topology declare the same
                # atom count as the trajectory? None means the format does not
                # say. mwf rejects a mismatch at its first step, so a False here
                # is a run that cannot be processed until a topology is supplied.
                "topology_natom_match": topo_match,
                "header_errors": "; ".join(errs) if errs else None,
                "dt_note": ("production files here have DIFFERENT frame spacings; "
                            "prod_ns sums each with its own dt, and the single "
                            "dt_ps_production is withheld as unrepresentative"
                            if dt_mixed else None),
                "status": "current",
            })
            if verbose and n_dirs % 50 == 0:
                print(f"      {n_dirs} run dirs ({time.time()-t0:.0f}s)", flush=True)
        if verbose:
            print(f"      {n_dirs} run dirs in {time.time()-t0:.0f}s", flush=True)

    # A substrate with NO trajectories is legitimate -- a QM project, or someone
    # who only has analysis output. Returning a column-less empty frame made
    # substrate_build skip writing runs.parquet at all, and the acceptance suite
    # then aborted with "cannot continue without runs.parquet": a hard red
    # failure for a colleague who had done nothing wrong. Empty but correctly
    # shaped is the honest answer.
    RUN_COLUMNS = {
        "run_id": "object", "owner": "object", "corpus": "object",
        "rel_path": "object", "replica": "object", "replica_index": "object",
        "stage_dir": "object", "system": "object", "condition": "object",
        "engine": "object", "n_traj_files": "int64", "traj_gb": "float64",
        "traj_gb_production": "float64", "n_frames_production": "float64",
        "n_frames_merged": "float64", "prod_ns_header": "float64",
        "dt_ps_production": "float64", "dt_mixed": "bool", "natom": "float64",
        "natom_source": "object", "natom_conflict": "bool",
        "n_frames_method": "object",
        "n_minimisation": "int64", "n_heating": "int64",
        "n_equilibration": "int64", "n_production": "int64",
        "n_merged_files": "int64", "n_stripped": "int64",
        "n_superseded": "int64", "n_derived": "int64", "n_unresolved": "int64",
        "n_production_default": "int64", "n_stage_declared": "int64",
        "n_stripped_by_atoms": "int64",
        "n_other": "int64", "topology_path": "object", "topology_kind": "object",
        "topology_natom_match": "object", "header_errors": "object",
        "dt_note": "object", "status": "object",
    }
    df = pd.DataFrame(rows) if rows else config.empty_frame(RUN_COLUMNS)
    # This dict is a SECOND, hand-maintained copy of the column list built above,
    # and two copies of anything drift. A column added to the row but not here
    # would give an empty substrate a different shape from a populated one --
    # which is exactly the class of bug the empty frame was added to fix. Catch
    # it at the point of divergence rather than in someone's parquet a month on.
    if len(df):
        missing = set(df.columns) - set(RUN_COLUMNS)
        extra = set(RUN_COLUMNS) - set(df.columns)
        if missing or extra:
            raise RuntimeError(
                "trajscan RUN_COLUMNS is out of step with the rows it builds: "
                f"rows have {sorted(missing)} that the schema lacks; schema has "
                f"{sorted(extra)} that the rows lack. Update RUN_COLUMNS.")
        for c in ("owner", "corpus", "replica", "system", "condition", "engine",
                  "topology_kind", "n_frames_method", "natom_source",
                  "stage_dir", "status"):
            df[c] = df[c].astype("category")
    return df
