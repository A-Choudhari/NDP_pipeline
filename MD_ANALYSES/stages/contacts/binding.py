#!/usr/bin/env python3
"""
binding.py -- the binding free energies that ALREADY exist, read off disk.

WHY THIS STAGE READS INSTEAD OF COMPUTES. An end-point free energy is not an
observable of a trajectory the way an RMSD is: it is the output of a separate
program, run with its own topologies, its own solvent model and its own frame
selection, and re-running it here would take a GPU-day and produce a second set
of numbers to disagree with the first. What this workspace is missing is not
another MM/GBSA run -- 1F4's finished months ago -- it is a path by which those
nine numbers reach the same tidy table as everything else, carrying the
conditions they were produced under. So this stage reads no trajectory, has no
stride and no sampling floor. Its whole job is transport plus caveat.

THE METHOD IS PART OF THE ANSWER, NEVER AN ASSUMPTION. "Binding energy" is a
family, not a quantity: MM/GBSA, MM/PBSA and FEP/TI give numbers in the same
units that mean different things and are not interchangeable. This stage
therefore DISCOVERS which of them was run, records the method and the directory
it was found in, and refuses to print a value without them. A table that says
`-5.0 kcal/mol` and not `MMPBSA.py, generalized Born, igb=2, saltcon=0.15 M`
has already lost the only information that makes the number checkable.

AN ABSENT METHOD IS NOT AN ERROR. Three of the four CDR systems have no
end-point calculation at all. For them this stage finds nothing, writes an
empty table and raises nothing -- but "found nothing" and "could not look" are
different, and a directory the search could not read produces a row saying so.

WHERE IT LOOKS. From each MD's run directory, upward `LEVELS_UP` levels, and
down to `MAX_DEPTH` below each of those. For 1F4 that reaches
`0-MD/ANALYSIS/MMGBSA/<protomer>/<replica>/FINAL_MMGBSA.dat`, which is where
these results live. The radius is deliberately short. One level further up is
`CDR/`, whose shared `SCRIPTS/` directory holds an abandoned CaFE/NAMD MM-PBSA
attempt on 1F4 replica 1 -- and a search wide enough to find it for 1F4 is wide
enough to attribute it to JL427, YZ23 and the apo trimer too, which have
nothing to do with it. (That attempt produced no binding free energy in any
case: APBS died on a linker error, its PB log is empty, and only the gas-phase
MM term was ever written. The CaFE reader below therefore treats a working
directory with no `Delta:` block as an attempt, not a result, and says so in a
row rather than ingesting the MM term as if it were a ΔG.)

POOLING, AND WHY `n=9` IS NOT NINE SAMPLES. The nine values are three protomers
x three replicas. The three protomers of one trimer share a box, a thermostat
and a barostat -- `core/io.py`'s SCOPES docstring names them pseudoreplicates --
so the aggregate row reports `sd` over all nine pooled units but `sem` over the
three REPLICA means only. Dividing the nine-unit spread by sqrt(9) would shrink
the error bar by a factor this data does not earn.

WHAT THE NUMBERS ARE NOT. MM/GBSA here has no configurational entropy term, is
computed from a single trajectory of the complex (so every bonded term in the
difference is identically zero by construction), and depends strongly on the
implicit-solvent model it was run with. It ranks; it does not predict. That
sentence is not left to this docstring: it is written into the table as a
`resolved=False` row, beside the numbers, so that a reader who converts
-5 kcal/mol into a Kd has had to scroll past it.
"""

import glob
import math
import os
import re

from core import io

STAGE = "binding_energy"
#: No trajectory is read, so there is no frame spacing to be too coarse.
MIN_SPACING_NS = None
DEFAULT_STRIDE = None
UNIT = "kcal/mol"

#: How far up from a run directory to start looking, and how far down from
#: there. See WHERE IT LOOKS above: widening either of these starts attributing
#: one system's results to its neighbours.
LEVELS_UP = 2
MAX_DEPTH = 5

#: Directories that hold trajectories, caches and force fields. Skipped so the
#: search over an NFS mount stays in the tens of milliseconds it measured at.
#: `namd*` is deliberately NOT here: a FEP set is often written into a
#: directory named after the engine, and pruning it would hide exactly the
#: output this stage exists to find.
_PRUNE = re.compile(r"^(STRIP.*|__pycache__|\.git|\.ipynb_checkpoints|"
                    r"CHARMM36m|toppar|par_.*)$", re.I)

#: method key -> (human name, result globs, attempt-marker globs, root markers).
#:
#: ATTEMPT MARKERS exist so a run that started and produced nothing is
#: distinguishable from a method that was never tried.
#:
#: ROOT MARKERS name the directory that owns a set of results -- the one
#: holding the inputs and the summary tables. It matters because MMPBSA.py
#: writes one result file per subsystem: 1F4 has nine of them under
#: `MMGBSA/<protomer>/<replica>/`, and a discovery that returned nine separate
#: findings would pool each of them with itself and report nine aggregates of
#: n=1. The root is what makes the nine one calculation.
METHODS = {
    "mmpbsa.py": (
        "MM/GBSA or MM/PBSA (AmberTools MMPBSA.py)",
        ("FINAL_MMGBSA.dat", "FINAL_MMPBSA.dat", "FINAL_RESULTS_MMPBSA.dat",
         "FINAL_RESULTS_MMGBSA.dat"),
        ("_MMPBSA_*", "mmpbsa.log", "mmgbsa.log"),
        ("*mmgbsa*.in", "*mmpbsa*.in", "MMGBSA_SUMMARY.tsv",
         "MMPBSA_SUMMARY.tsv", "BINDING_TOTALS.tsv"),
    ),
    "cafe": (
        "MM/PBSA (CaFE plugin, NAMD + APBS)",
        ("result_*.log", "*_result.log"),
        ("com_mm.log", "rec_mm.log", "lig_mm.log", "_mmpbsa_com_tmp.dcd"),
        ("run.tcl", "average_mmpbsa.py"),
    ),
    "fep": (
        "alchemical free energy (FEP / TI)",
        ("*.fepout", "ti.out", "bar.dat", "mbar*.txt", "*.mbar"),
        ("fep*.conf", "ti*.conf"),
        ("fep*.conf", "ti*.conf", "*.fepin"),
    ),
}

#: `param_value` for the caveat rows. Written once, read by the page, and the
#: only place in this workspace where what MM/GBSA is NOT is stated next to
#: what it is.
CAVEAT = (
    "This is an end-point estimate, not a binding affinity. There is no "
    "configurational entropy term (-T dS is not computed), the three states "
    "come from ONE trajectory of the complex so every bonded term in the "
    "difference is exactly zero by construction, and the value depends "
    "strongly on the implicit-solvent model it was run with. Values computed "
    "with a different igb or saltcon are not comparable with these, and none "
    "of them converts to a Kd. Use them to RANK interfaces against each "
    "other; a -5 kcal/mol here does not mean micromolar."
)


# ---------------------------------------------------------------- arithmetic
def _mean(vals):
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def _sd(vals):
    """Sample standard deviation. Returns None for n<2, never 0.0.

    A single value has no spread; reporting 0.0 for it would draw a band of
    zero width around a number with no error bar at all.
    """
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return None
    mu = sum(vals) / len(vals)
    return math.sqrt(sum((v - mu) ** 2 for v in vals) / (len(vals) - 1))


# ------------------------------------------------------------------ discovery
def search_roots(spec, levels_up=LEVELS_UP):
    """Directories to look in: each run directory and `levels_up` ancestors.

    The SOLVATED trajectory is preferred for locating the run, because the
    stripped copy lives one level deeper (`<run>/STRIP*/merged.dcd`) and taking
    its dirname would start the search inside the strip directory and miss the
    `ANALYSIS/` sibling entirely.
    """
    roots = []
    for md in spec.mds:
        traj = (list(md.solvated_trajectories) or list(md.trajectories)
                or [None])[0]
        d = os.path.abspath(os.path.dirname(traj)) if traj else None
        if not d:
            continue
        for _ in range(levels_up + 1):
            if d not in roots:
                roots.append(d)
            parent = os.path.dirname(d)
            if parent == d:
                break
            d = parent
    # Deepest first, so a subtree is walked from the root CLOSEST to it and the
    # depth limit is measured from there. Walking 1-1F4/ first would consume
    # the whole budget on the way down and stop above the result files.
    return sorted(roots, key=lambda p: -p.count(os.sep))


def _walk(roots, max_depth=MAX_DEPTH):
    """(directory, [entry names]) under every root, bounded and de-duplicated.

    Entries are files AND subdirectory names, because `_MMPBSA_*` -- the
    fingerprint of an AmberTools run that did not clean up after itself -- is a
    directory. One listing per directory: testing for it with a second
    `os.listdir` would triple the syscalls of the whole search on an NFS mount.

    Yields the unreadable directories too, as (path, None): a search that could
    not read a directory must not be reported the same way as a search that
    read it and found nothing.
    """
    seen, failed = set(), []
    for root in roots:
        if not os.path.isdir(root):
            continue
        base = root.rstrip(os.sep).count(os.sep)
        for dirpath, dirnames, filenames in os.walk(
                root, onerror=lambda e: failed.append(getattr(e, "filename",
                                                              str(e)))):
            real = os.path.realpath(dirpath)
            if real in seen:
                dirnames[:] = []
                continue
            seen.add(real)
            entries = list(filenames) + list(dirnames)
            if dirpath.rstrip(os.sep).count(os.sep) - base >= max_depth:
                dirnames[:] = []
            dirnames[:] = [d for d in dirnames if not _PRUNE.match(d)]
            yield dirpath, entries
    for p in failed:
        yield p, None


def _match(name, pattern):
    import fnmatch
    return fnmatch.fnmatch(name, pattern)


def method_root(dirpath, method, levels=4):
    """The directory that OWNS a set of results: the one holding the inputs.

    Found by looking for the method's own root markers, not by counting two
    levels up. The `<protomer>/<replica>/` layout under `MMGBSA/` is a
    convention of this project, and a convention that is assumed rather than
    checked is one that breaks silently on the next system.
    """
    markers = METHODS[method][3]
    d = dirpath
    for _ in range(levels):
        if any(glob.glob(os.path.join(d, pat)) for pat in markers):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return dirpath


def discover(spec, roots=None):
    """What end-point method was run near these MDs, and where.

    Returns (found, unreadable). `found` is a list of dicts, ONE PER
    (method, root):

        {"method", "name", "root", "dirs": [...], "results": [paths],
         "attempt": bool}

    Grouping by root, not by directory, is the whole point: MMPBSA.py writes
    one result file per subsystem, 1F4 has nine of them under
    `MMGBSA/<protomer>/<replica>/`, and a discovery that returned nine separate
    findings would pool each with itself and report nine aggregates of n=1.

    `attempt` marks a root carrying a method's fingerprints but none of its
    result files: a run that started and produced nothing. That is a finding,
    not an absence, and it is reported as one.
    """
    roots = roots if roots is not None else search_roots(spec)
    groups, unreadable = {}, []
    for dirpath, entries in _walk(roots):
        if entries is None:
            unreadable.append(dirpath)
            continue
        for key, (label, result_pats, attempt_pats, _) in METHODS.items():
            got = sorted({n for pat in result_pats
                          for n in entries if _match(n, pat)})
            att = (not got) and any(_match(n, pat) for pat in attempt_pats
                                    for n in entries)
            if not got and not att:
                continue
            root = method_root(dirpath, key)
            rec = groups.setdefault((key, root), {
                "method": key, "name": label, "root": root,
                "dirs": [], "results": []})
            rec["dirs"].append(dirpath)
            rec["results"] += [os.path.join(dirpath, n) for n in got]
    found = []
    for _, rec in sorted(groups.items()):
        rec["results"] = sorted(set(rec["results"]))
        rec["dirs"] = sorted(set(rec["dirs"]))
        rec["attempt"] = not rec["results"]
        found.append(rec)
    return found, unreadable


# -------------------------------------------------------------- MMPBSA.py
#: `DELTA TOTAL  -5.0019  11.1274  0.7010` -- and also `1-4 VDW`, `DELTA G
#: gas`, `G solv`: the component name contains spaces, so it is whatever is
#: left once the three trailing numbers are taken off the end.
_NUM = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")

#: The block whose numbers are the answer. MMPBSA.py writes the absolute
#: energies of complex, receptor and ligand above it; those are five-figure
#: numbers that differ by the binding energy, and reading the wrong block gives
#: a value 2000x too large that is nonetheless perfectly well formed.
_DIFF_HEAD = re.compile(r"^Differences\s*\(Complex\s*-\s*Receptor\s*-\s*"
                        r"Ligand\)", re.I)
_MODEL_HEAD = re.compile(r"^(GENERALIZED BORN|POISSON BOLTZMANN|"
                         r"POISSON-BOLTZMANN)\s*:?\s*$", re.I)
_FRAMES = re.compile(r"using\s+([\d.]+)\s+complex frames", re.I)

#: GENERALIZED BORN -> GB, POISSON BOLTZMANN -> PB. The distinction is the
#: difference between MM/GBSA and MM/PBSA, which is the difference between two
#: methods whose numbers must never be averaged together.
_MODELS = {"GENERALIZED BORN": "GB", "POISSON BOLTZMANN": "PB",
           "POISSON-BOLTZMANN": "PB"}


def subset_name(component):
    """`DELTA G gas` -> `DELTA_G_GAS`, so a subset is one token in a table."""
    return re.sub(r"\s+", "_", component.strip()).upper()


def parse_mmpbsa(path):
    """The Differences block of a FINAL_*.dat, per solvent model.

    Returns {"n_frames": int|None, "models": {"GB": {COMPONENT: (mean, sd,
    sem)}}, "recomputed": False}. The file's own numbers are used: MMPBSA.py
    already computed the mean, the standard deviation and the standard error of
    the mean over the frames it used, and recomputing them from PER_FRAME.csv
    would give a second answer that can silently differ (it strides, it drops,
    and its DELTA block is written at a different precision).

    An empty models dict is a legitimate return -- a truncated or still-running
    file has no Differences block -- and the caller decides what to do about
    it. Raising here would make one bad file take a whole system's table with
    it.
    """
    out = {"n_frames": None, "models": {}, "recomputed": False,
           "n_used": None}
    model, in_diff = None, False
    try:
        with open(path, errors="ignore") as fh:
            lines = fh.readlines()
    except OSError:
        return out
    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()
        if out["n_frames"] is None:
            m = _FRAMES.search(stripped)
            if m:
                out["n_frames"] = int(float(m.group(1)))
        head = _MODEL_HEAD.match(stripped)
        if head:
            model = _MODELS.get(head.group(1).upper(), head.group(1).upper())
            in_diff = False
            continue
        if _DIFF_HEAD.match(stripped):
            in_diff = True
            continue
        # A bare `Complex:` / `Receptor:` / `Ligand:` heading ENDS the
        # differences block of the previous model. Without this, a
        # FINAL_RESULTS file that carries GB then PB would append PB's absolute
        # complex energies to GB's differences.
        if stripped.endswith(":") and not stripped.startswith("DELTA"):
            in_diff = False
            continue
        if not in_diff or not stripped or stripped.startswith("-"):
            continue
        f = stripped.split()
        if len(f) < 4 or not all(_NUM.match(x) for x in f[-3:]):
            continue
        name = " ".join(f[:-3])
        if not name:
            continue
        mean, sd, sem = (float(x) for x in f[-3:])
        out["models"].setdefault(model or "GB", {})[subset_name(name)] = \
            (mean, sd, sem)
    return out


def parse_per_frame(path):
    """The DELTA block of a PER_FRAME.csv, as a LAST resort.

    Only reached when the summary line is missing, and the caller marks every
    row built from it so a reader can tell the two apart. The s.e.m. here is
    sd/sqrt(n) over frames, which treats consecutive MD frames as independent
    and therefore UNDERSTATES the true error -- one more reason this is the
    fallback and not the path.
    """
    cols, rows = None, []
    try:
        with open(path, errors="ignore") as fh:
            in_delta = False
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                if s.lower().startswith("delta energy terms"):
                    in_delta, cols = True, None
                    continue
                if s.lower().endswith("energy terms"):
                    in_delta = False
                    continue
                if not in_delta:
                    continue
                f = s.split(",")
                if cols is None:
                    cols = [c.strip() for c in f]
                    continue
                try:
                    rows.append([float(x) for x in f])
                except ValueError:
                    continue
    except OSError:
        return {}
    if not cols or not rows:
        return {}
    out = {}
    for i, name in enumerate(cols):
        if i == 0 or i >= len(rows[0]):
            continue
        vals = [r[i] for r in rows if len(r) > i]
        mu, sd = _mean(vals), _sd(vals)
        sem = (sd / math.sqrt(len(vals))) if sd is not None else None
        out[subset_name(name)] = (mu, sd, sem)
    return {"n_frames": len(rows), "models": {"GB": out}, "recomputed": True,
            "n_used": len(rows)}


# ------------------------------------------------------------------ CaFE
_CAFE_DELTA = re.compile(r"^\s*Delta\s*:\s*$", re.I)


def parse_cafe(path):
    """The Delta block of a CaFE result log.

    CaFE reports `name: nframes mean sd` and no s.e.m., so the s.e.m. is
    derived as sd/sqrt(n) and the row says where it came from. A log with no
    Delta block means the run stopped before the PB/SA term -- which is exactly
    what happened to the one attempt in this workspace -- and returns nothing,
    so the gas-phase MM term cannot be mistaken for a binding free energy.
    """
    out, in_delta = {}, False
    try:
        with open(path, errors="ignore") as fh:
            for line in fh:
                if _CAFE_DELTA.match(line):
                    in_delta = True
                    continue
                if not in_delta:
                    continue
                s = line.strip()
                if not s or s.startswith("-") or s.startswith("*"):
                    if out:
                        break
                    continue
                f = s.replace(":", " ").split()
                if len(f) < 4 or not all(_NUM.match(x) for x in f[-3:]):
                    continue
                n, mean, sd = int(float(f[-3])), float(f[-2]), float(f[-1])
                sem = sd / math.sqrt(n) if n > 0 else None
                out[subset_name(" ".join(f[:-3]))] = (mean, sd, sem, n)
    except OSError:
        return {}
    if not out:
        return {}
    n = max(v[3] for v in out.values())
    return {"n_frames": n, "recomputed": False, "n_used": n,
            "models": {"PB": {k: v[:3] for k, v in out.items()}}}


# ------------------------------------------------------------------- FEP/TI
_FEP_TOTAL = re.compile(r"#Free energy change for (?:the )?whole window.*?"
                        r"=\s*([-\d.eE+]+)", re.I)
_FEP_LAST = re.compile(r"#Free energy change for last step.*?"
                       r"=\s*([-\d.eE+]+)", re.I)


def parse_fepout(path):
    """The terminal free-energy line of a NAMD `.fepout`, UNVERIFIED.

    NOTHING IN THIS WORKSPACE HAS FEP OUTPUT. This reader follows the
    documented NAMD format, but it has never been run against a real file here,
    so every row built from it is written `resolved=False` with that said in
    the note. A number from an untested parser that looks like every other
    number in the table is worse than no number.
    """
    total = None
    try:
        with open(path, errors="ignore") as fh:
            for line in fh:
                m = _FEP_TOTAL.search(line) or _FEP_LAST.search(line)
                if m:
                    total = float(m.group(1))
    except (OSError, ValueError):
        return {}
    if total is None:
        return {}
    return {"n_frames": None, "recomputed": False, "n_used": None,
            "models": {"FEP": {"DELTA_TOTAL": (total, None, None)}}}


# ------------------------------------------------------------------- inputs
_NAMELIST = re.compile(r"^&(\w+)")
_KV = re.compile(r"(\w+)\s*=\s*(\"[^\"]*\"|'[^']*'|[^,\s]+)")


def read_namelists(path):
    """`{namelist: {key: value}}` from an MMPBSA.py `.in` file.

    The method's identity lives here and nowhere else: `igb`, `saltcon` and
    `idecomp` are echoed into the top of every FINAL_*.dat, but only as the
    comment block, and a reader who quotes a value without them has quoted a
    number whose solvent model is unknown.
    """
    out, current = {}, None
    try:
        with open(path, errors="ignore") as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("/"):
                    current = None
                    continue
                m = _NAMELIST.match(s)
                if m:
                    current = m.group(1).lower()
                    out.setdefault(current, {})
                    s = s[m.end():]
                if current is None:
                    continue
                for k, v in _KV.findall(s):
                    out[current][k.lower()] = v.strip("\"'")
    except OSError:
        return {}
    return out


def method_parameters(root):
    """The parameters every result under `root` was produced with.

    A value per key, plus the set of files that disagreed about it. Two
    protomers run at different `igb` are two methods wearing one name, and the
    aggregate over them would be meaningless -- so a disagreement is returned,
    not resolved.

    MMPBSA.py's `.in` namelists only. CaFE keeps its settings in a `.tcl` run
    script and NAMD FEP in a `.conf`, and neither is parsed here: those methods
    therefore fall back to reporting their name alone in the caveat row, which
    is less than this gives for MM/GBSA and is said rather than papered over.
    """
    files = sorted(set(glob.glob(os.path.join(root, "*mmgbsa*.in")) +
                       glob.glob(os.path.join(root, "*mmpbsa*.in"))))
    values, conflicts = {}, {}
    for p in files:
        for group, kv in read_namelists(p).items():
            for k, v in kv.items():
                if k == "print_res":
                    continue       # a per-protomer residue list, not a method
                key = f"{group}.{k}"
                if key in values and values[key] != v:
                    conflicts.setdefault(key, set()).add(values[key])
                    conflicts[key].add(v)
                values[key] = v
    return values, {k: sorted(v) for k, v in conflicts.items()}, files


# ------------------------------------------------------------------- locate
def locate(path, md_ids):
    """(protomer, replica) for a result file, from the directory names.

    The replica is the first path component that matches one of the spec's own
    MD ids, and the protomer is the component above it. Matching against the
    spec rather than counting two levels up is what makes this survive a layout
    where the replica is not where it is here -- and when nothing matches, the
    caller is told (replica=None), rather than being handed the name of
    whatever directory happened to be there.
    """
    parts = os.path.abspath(path).split(os.sep)
    ids = {str(i) for i in md_ids}
    for i in range(len(parts) - 2, 0, -1):
        if parts[i] in ids:
            return (parts[i - 1] if i >= 1 else None), parts[i]
    return None, None


# --------------------------------------------------------------------- run
def run(spec, out_dir, write=True, roots=None):
    """Ingest every end-point free energy that exists near this system's MDs.

    Reads no trajectory. Returns the rows; writes nothing but an empty table
    when the system has no such results, which is the common case.
    """
    rows = io.Rows(system=spec.id, stage=STAGE, stride=None, resolved=True,
                   ns_per_frame=None)
    found, unreadable = discover(spec, roots=roots)
    md_ids = [m.id for m in spec.mds]
    inputs = []

    if unreadable:
        # "Nothing found" and "could not look" must not be the same row.
        rows.add(scope="system", observable="search_incomplete",
                 value=len(unreadable), unit="count", resolved=False,
                 floor_note=("the search could not read "
                             f"{len(unreadable)} director(ies) "
                             f"(first: {os.path.basename(unreadable[0])}); "
                             "an absence below may be a permission problem, "
                             "not an absent calculation"))

    for hit in found:
        root = hit["root"]
        params, conflicts, in_files = method_parameters(root)
        inputs += in_files

        rows.add(scope="system", observable="method_found", value=1,
                 unit="flag", subset=hit["method"], label=hit["name"],
                 method=hit["method"],
                 param_name="results_dir",
                 param_value=f"{os.path.basename(root)} "
                             f"({len(hit['results'])} result file(s) in "
                             f"{len(hit['dirs'])} director(ies))",
                 resolved=not hit["attempt"],
                 floor_note=None if not hit["attempt"] else
                 (f"a {hit['name']} working directory with no result file: "
                  "the run started and produced no free energy. Nothing from "
                  "it is ingested."),
                 source=os.path.basename(root))

        if hit["attempt"]:
            continue

        for key, val in sorted(params.items()):
            rows.add(scope="system", observable="method_parameter",
                     subset=hit["method"], param_name=key, param_value=val,
                     value=None, unit=None, method=hit["method"],
                     resolved=key not in conflicts,
                     floor_note=None if key not in conflicts else
                     (f"the input files under {os.path.basename(root)} "
                      f"disagree about {key} ({', '.join(conflicts[key])}): "
                      "results produced with different settings are not "
                      "comparable and must not be pooled"))

        units = _ingest(rows, hit, md_ids, root)
        _aggregate(rows, hit, units)
        _caveats(rows, hit, units, params)

    man = io.manifest(STAGE, spec, sorted(set(inputs)), stride=None, extra={
        # The nine result files share ONE basename, and io.manifest keys its
        # inputs by basename -- so passing them there would collapse nine files
        # into one entry and leave `needs_run` blind to eight of them changing.
        # Keyed by path relative to the method root instead, which is both
        # unique and free of the absolute paths io.py keeps out of the data.
        "binding_results": {
            f"{os.path.basename(h['root'])}/{os.path.relpath(p, h['root'])}":
                _stat(p)
            for h in found for p in h["results"]},
        "binding_methods": sorted({h["method"] for h in found}),
        "search_roots": [os.path.basename(r)
                         for r in (roots if roots is not None
                                   else search_roots(spec))],
        "unreadable": len(unreadable),
    })
    if write:
        io.write(out_dir, f"{spec.label}__{STAGE}", rows, manifest=man)
    return rows


def _stat(path):
    try:
        st = os.stat(path)
        return [st.st_size, int(st.st_mtime)]
    except OSError:
        return None


def _reader(method):
    return {"mmpbsa.py": parse_mmpbsa, "cafe": parse_cafe,
            "fep": parse_fepout}[method]


def _ingest(rows, hit, md_ids, root):
    """One `scope="run"` row per (protomer, replica, model, component).

    Returns {(model, component): [(protomer, replica, value), ...]} for the
    aggregate to pool.
    """
    units, read = {}, _reader(hit["method"])
    for path in sorted(hit["results"]):
        parsed = read(path)
        # A FINAL_*.dat with no Differences block is truncated or still being
        # written. PER_FRAME.csv beside it holds the same DELTA series, so the
        # run is recoverable -- but every row built that way SAYS SO, because
        # its s.e.m. is over frames, not over the file's own estimate.
        if not parsed.get("models") and hit["method"] == "mmpbsa.py":
            alt = os.path.join(os.path.dirname(path), "PER_FRAME.csv")
            if os.path.exists(alt):
                parsed = parse_per_frame(alt)
        if not parsed.get("models"):
            protomer, replica = locate(path, md_ids)
            rows.add(scope="run", replica=replica, scope_id=replica,
                     protomer=protomer, observable="binding_energy",
                     subset="DELTA_TOTAL", value=None, unit=UNIT,
                     method=hit["method"], resolved=False,
                     source=os.path.relpath(path, root),
                     floor_note=("no summary block and no usable PER_FRAME.csv "
                                 "in this result file: it is truncated, "
                                 "malformed or still being written"))
            continue

        protomer, replica = locate(path, md_ids)
        if replica is None:
            # Better an unlabelled row than a row labelled with the name of
            # whatever directory the file happened to sit in.
            rows.add(scope="system", observable="unplaced_result", value=1,
                     unit="flag", method=hit["method"], resolved=False,
                     source=os.path.relpath(path, root),
                     floor_note=("no directory on this path matches an MD id "
                                 f"of this system ({', '.join(md_ids)}), so "
                                 "the replica it belongs to is unknown and it "
                                 "is not pooled"))
            continue

        n = parsed.get("n_frames")
        for model, comps in sorted(parsed["models"].items()):
            for comp, (mean, sd, sem) in sorted(comps.items()):
                rows.add(
                    scope="run", replica=replica,
                    scope_id=f"{protomer}/{replica}" if protomer else replica,
                    protomer=protomer, observable="binding_energy",
                    subset=comp, value=mean, unit=UNIT, sd=sd, sem=sem, n=n,
                    n_frames=n, ns_sampled=None, label=model,
                    method=hit["method"], source=os.path.relpath(path, root),
                    resolved=not parsed.get("recomputed")
                    and hit["method"] != "fep",
                    floor_note=_row_note(parsed, hit))
                units.setdefault((model, comp), []).append(
                    (protomer, replica, mean))
    return units


def _row_note(parsed, hit):
    if hit["method"] == "fep":
        return ("parsed by a reader that has never been run against a real "
                "FEP output in this workspace: check it against the file "
                "before quoting it")
    if parsed.get("recomputed"):
        return ("recomputed from PER_FRAME.csv because the summary block was "
                "missing; the s.e.m. is sd/sqrt(n) over frames, which treats "
                "consecutive frames as independent and understates the error")
    return None


def _aggregate(rows, hit, units):
    """Pool protomers and replicas: sd over all of them, s.e.m. over replicas.

    The three protomers of one trimer share a box, a thermostat and a barostat,
    so they are pseudoreplicates (core/io.py, SCOPES). `n` is the number of
    pooled units and `sd` their spread -- but `sem` is computed from the
    REPLICA means alone, because dividing the nine-unit spread by sqrt(9) would
    claim an independence this design does not have.
    """
    for (model, comp), triples in sorted(units.items()):
        vals = [v for _, _, v in triples if v is not None]
        if not vals:
            continue
        by_replica = {}
        for _, replica, v in triples:
            if v is not None:
                by_replica.setdefault(replica, []).append(v)
        replica_means = [_mean(v) for v in by_replica.values()]
        rsd = _sd(replica_means)
        sem = (rsd / math.sqrt(len(replica_means))) if rsd is not None else None
        protomers = sorted({p for p, _, v in triples if p and v is not None})
        rows.add(
            scope="aggregate", scope_id="pooled", replica=None, subset=comp,
            observable="binding_energy", value=_mean(vals), unit=UNIT,
            sd=_sd(vals), sem=sem, n=len(vals), label=model,
            method=hit["method"], resolved=True,
            param_name="pooled_over",
            param_value=(f"{len(protomers)} protomer(s) x "
                         f"{len(by_replica)} replica(s); sem over replicas "
                         "only -- protomers of one trimer share a box and a "
                         "thermostat and are pseudoreplicates"),
            floor_note=None if len(by_replica) > 1 else
            ("one replica only: the s.e.m. over replicas cannot be estimated, "
             "and the spread shown is between protomers of the same box"))


def _caveats(rows, hit, units, params):
    """The rows a reader must scroll past to misuse the numbers above.

    Written as `resolved=False` so a figure greys them and a caption refuses
    to quote them as a value -- they are not measurements, they are the terms
    on which the measurements may be read.
    """
    model = "PB" if hit["method"] == "cafe" else (
        "FEP" if hit["method"] == "fep" else None)
    settings = ", ".join(f"{k.split('.')[-1]}={v}" for k, v in
                         sorted(params.items())
                         if k.split(".")[-1] in ("igb", "saltcon", "istrng",
                                                 "idecomp", "inp", "pbtemp"))
    rows.add(scope="system", observable="interpretation",
             subset=hit["method"], value=None, unit=UNIT, resolved=False,
             method=hit["method"], label=model,
             param_name="what_this_is_not",
             param_value=(settings or hit["name"]),
             floor_note=CAVEAT)

    # WHY `ns_per_frame` IS BLANK ON EVERY ROW THIS STAGE WRITES. io.py calls
    # that column never-optional, and it is right to -- but the honest value
    # here is unknown, and the MD's own 0.2 ns/frame is the one value it
    # certainly is not. Writing the MD's spacing would be wrong by whatever
    # stride the end-point trajectory was prepared at, and wrong in a way that
    # reads as a fact.
    rows.add(scope="system", observable="frame_spacing_unknown",
             subset=hit["method"], value=None, unit="ns", resolved=False,
             method=hit["method"],
             param_name="general.interval",
             param_value=params.get("general.interval", "unknown"),
             floor_note=(
                 "the end-point calculation ran over a trajectory that was "
                 "strided and equilibration-trimmed during preparation, and "
                 "its output records only the interval it applied to THAT "
                 "file (not the MD's). ns_per_frame is therefore left blank on "
                 "every row here rather than filled with the MD's own frame "
                 "spacing, which would be wrong by the preparation stride. "
                 "n_frames is what the calculation actually averaged over."))

    # What the per-run spread means for pooling, in the data rather than only
    # in a report: the within-run frame spread is not an error bar on the mean,
    # and it dwarfs the between-interface differences these numbers are used to
    # rank.
    key = next((k for k in sorted(units) if k[1] == "DELTA_TOTAL"), None)
    if not key:
        return
    run_rows = [r for r in rows.rows
                if r.get("scope") == "run" and r.get("subset") == "DELTA_TOTAL"
                and r.get("method") == hit["method"]
                and r.get("sd") is not None]
    if not run_rows:
        return
    within = _mean([r["sd"] for r in run_rows])
    vals = [v for _, _, v in units[key] if v is not None]
    between = _sd(vals)
    mu = _mean(vals)
    agg = next((r for r in rows.rows if r.get("scope") == "aggregate"
                and r.get("subset") == "DELTA_TOTAL"
                and r.get("method") == hit["method"]), {})
    sem = agg.get("sem")
    n_rep = len({r["replica"] for r in run_rows})

    note = (f"the frame-to-frame spread WITHIN a run is {within:.1f} kcal/mol "
            f"against a pooled mean of {mu:.1f}: it is LARGER THAN THE MEAN "
            "ITSELF, so no single frame and no short window carries this "
            "number -- only the run mean does, and the distribution it comes "
            "from straddles zero. ")
    note += (f"Pooling does not divide that spread by sqrt({len(vals)}): the "
             f"{len(vals)} units are {n_rep} replica(s) worth of "
             "pseudoreplicated protomers. ")
    if sem is not None:
        note += (f"The honest error bar on the pooled mean is the replica "
                 f"s.e.m., {sem:.1f} kcal/mol. ")
    if between is not None:
        note += (f"The between-unit spread is {between:.1f} kcal/mol, and that "
                 "is the most these values can resolve: a difference between "
                 "protomers smaller than it is not a difference.")
    rows.add(scope="system", observable="spread_warning", subset="DELTA_TOTAL",
             value=within, unit=UNIT, resolved=False, method=hit["method"],
             n=len(vals), sd=between, sem=sem,
             param_name="within_run_sd_over_pooled_sd",
             param_value=(f"{within:.2f} vs {between:.2f}"
                          if between is not None else f"{within:.2f}"),
             floor_note=note)


def table(rows):
    """{(protomer, replica): DELTA TOTAL} plus the pooled value, for a page."""
    out, pooled = {}, {}
    for r in rows.rows:
        if r.get("subset") != "DELTA_TOTAL" or r.get("value") is None:
            continue
        if r.get("scope") == "run":
            out[(r.get("protomer"), r.get("replica"))] = r["value"]
        elif r.get("scope") == "aggregate":
            pooled = {"value": r["value"], "sd": r.get("sd"),
                      "sem": r.get("sem"), "n": r.get("n")}
    return {"per_run": out, "pooled": pooled}
