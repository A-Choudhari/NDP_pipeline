#!/usr/bin/env python3
"""
cpptraj.py -- drive cpptraj, and refuse to believe it worked until it proves it.

WHY ANY STAGE SHELLS OUT AT ALL. Two of the solvent observables are bound by
reading and by pair counting at the same time: a 3D occupancy grid bins 236,838
water oxygens per frame, and one RDF walks 236,838 x 20,265 pairs per frame.
cpptraj does both in C++ and -- the part that decides it on this mount --
strides at READ time (`trajin <file> 1 last <stride>`), so a stride-200 pass
over a 23 GB solvated replica touches about 120 MB instead of all of it.
Measured on 1F4 1-R1, that is 15-30x faster than the same arithmetic in Python
for the same answer, which is why `stages/solvent/grids.py` and
`stages/solvent/rdf.py` are the two stages that do not walk frames themselves.

WHY THIS MODULE IS MOSTLY CHECKING. cpptraj reports most setup failures as
WARNINGS and still exits 0. Ask it for a mask that matches nothing and it says

    Warning: First mask has no atoms.
    Warning: Setup incomplete for [radial]: Skipping
    Warning: File 'rdf.dat' has no sets containing data.

writes no data file at all -- or, for `grid`, an all-zero `.dx` -- and returns
success. That is precisely the failure this workspace spends its effort on: a
plausible empty result with no error. (Verified here, 2026-09-23: exit 0, no
`rdf.dat`, a 2,257-byte `.dx` of zeros.) So every run is judged on four
independent things, and any one of them fails the run:

    1. the process exit status;
    2. every `Error:` line in the log;
    3. the warnings that MEAN NOTHING HAPPENED -- an empty mask, a skipped
       setup, a data set with no data -- promoted to errors (`FATAL_WARNINGS`);
    4. each expected output file exists AND is non-empty.

WHICH BINARY RAN IS A RESULT, NOT A SETTING. Two builds live on this machine and
they are not the same program: V6.24.0 serial (the mwf_env conda env) and V6.4.4
with OpenMP (AmberTools 2022, which needs its `amber.sh` sourced first). The
version recorded in a manifest is read out of cpptraj's own startup banner in
the log of the run that did the work -- never inferred from the path, because a
path can be a symlink to anything and a manifest that names the wrong program
is worse than one that names none.

THE OpenMP FLAG AND core/threads.py. threads.py pins OMP_NUM_THREADS=1 for every
entry point, which is right for BLAS inside a process pool and exactly wrong for
the cpptraj child: it would run one-threaded while the manifest claimed OpenMP.
So the child's thread count is set explicitly here, and the banner line
`| N OpenMP threads available.` is read back as proof. Worth it where it is
worth it: `radial` on 1F4 costs 45 s/frame serial, 8.6 s at 8 threads and 4.7 s
at 32 -- the same numbers, ~10x sooner. `grid` is not threaded, so a grid-only
run gains nothing from the OpenMP build and pays an older version for it.
"""

import os
import re
import shutil
import subprocess
import time

import numpy as np

#: Where cpptraj comes from, in order: an environment variable, then the PATH.
#: The serial build is the default (`grid` gains nothing from threads).
#:   export CPPTRAJ=/path/to/cpptraj                 # serial
#:   export CPPTRAJ_OMP=/path/to/cpptraj.OMP         # optional OpenMP build
#:   export AMBER_SH=/path/to/amber/amber.sh         # sourced before the OpenMP build
SERIAL_BIN = os.environ.get("CPPTRAJ") or shutil.which("cpptraj") or ""

#: The OpenMP build, and the script that has to be sourced before it will run
#: (it needs AMBERHOME and its own libraries on the path). Without one, runs
#: that ask for OpenMP use the serial build.
OMP_BIN = os.environ.get("CPPTRAJ_OMP") or shutil.which("cpptraj.OMP") or ""
OMP_SETUP = os.environ.get("AMBER_SH") or (
    os.path.join(os.environ["AMBERHOME"], "amber.sh")
    if os.environ.get("AMBERHOME") else "")

#: Threads for the OpenMP build when a caller does not say: half the cores of
#: whatever machine this is (it was a fixed 24, half of one workstation's 48),
#: so two stages can run at once without oversubscribing it.
OMP_THREADS = int(os.environ.get("CPPTRAJ_OMP_THREADS")
                  or max(1, (os.cpu_count() or 2) // 2))

#: Warnings that are normal for a CHARMM PSF plus a NAMD DCD, and say nothing
#: about whether the analysis is right. Each one is here because it appears on
#: EVERY run of these systems; surfacing them would train a person to ignore
#: the list, which is how a real warning gets missed.
BENIGN_WARNINGS = (
    # PSF carries no force-field parameters, so cpptraj guesses bond lengths
    # from elements. Nothing here uses bonds.
    "Determining bond length parameters from element types",
    # The six trimer chains are not contiguous in the PSF. True, and harmless:
    # every mask in these stages is by atom number or residue name.
    "molecules have non-contiguous segments of atoms",
    "The 'fixatomorder' command can be used to reorder",
    "associated coordinates.",
    # A PSF has no box; the DCD does. Taking it from the trajectory is correct.
    "but topology box type is 'None'",
    "Setting topology box information from trajectory",
    # A stride longer than a short restart segment yields one frame from it.
    # Expected with stride 200 over an 18-frame pr.part9.dcd -- and the frames
    # actually processed are counted separately, so nothing hides here.
    "is so large that only 1 set will be processed",
)

#: Warnings that mean the run did nothing, and are therefore errors here. This
#: is the list that makes "cpptraj exited 0" untrustworthy on its own.
FATAL_WARNINGS = (
    "has no atoms",                 # `radial`: first/second mask empty
    "corresponds to 0 atoms",       # `grid` and friends
    "No atoms selected",
    "Setup incomplete for",         # the action was silently skipped
    "contains no data",
    "has no sets containing data",
    "Could not set up",
    "Error setting up",
)

_BANNER = re.compile(r"CPPTRAJ:.*?V(\S+)\s*\(([^)]*)\)\s*(\S*)")
_THREADS = re.compile(r"\|\s*(\d+)\s+OpenMP threads available")
_FRAMES = re.compile(r"Read\s+(\d+)\s+frames and processed\s+(\d+)\s+frames")


class CpptrajError(RuntimeError):
    """cpptraj did not do what it was asked. Never swallowed, never retried."""


class Result:
    """What one cpptraj run produced, and the evidence that it did."""

    __slots__ = ("binary", "version", "flavour", "threads", "script", "log",
                 "outputs", "errors", "warnings", "frames_read",
                 "frames_processed", "seconds", "returncode")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    def as_dict(self):
        """The block a manifest carries: which program ran, and what it saw.

        Names, not paths: core/io.py's contract is that nothing written carries
        an absolute path, and the script and log live in a known place anyway
        (`<bundle>/cpptraj/<stage>/<md>/`, see `workdir_for`). The BINARY is the
        exception and stays absolute, because which of the two builds ran is the
        whole point of recording it.
        """
        return {"binary": self.binary, "version": self.version,
                "build": self.flavour, "threads": self.threads,
                "frames_read": self.frames_read,
                "frames_processed": self.frames_processed,
                "seconds": round(self.seconds or 0.0, 1),
                "script": os.path.basename(self.script or ""),
                "log": os.path.basename(self.log or ""),
                "outputs": list(self.outputs or ()),
                "warnings": list(self.warnings or ())}

    def __repr__(self):                                    # pragma: no cover
        return (f"<cpptraj {self.version} {self.flavour or 'serial'} "
                f"{self.frames_processed} frames {self.seconds:.0f}s>")


# ---------------------------------------------------------------------- masks
def mask_from_indices(indices, max_ranges=200):
    """A cpptraj atom mask for a set of 0-based atom indices.

    WHY THIS EXISTS. cpptraj masks have no `segid`, and every selection in this
    workspace is written with one (`segid A120 B120 C120 A41 B41 C41`) because
    `chainID` collapses the three antibody copies. So the selection is resolved
    against the real topology by core/selections.py -- which is where an empty
    one becomes an error -- and handed to cpptraj as explicit 1-based atom
    numbers. The same atoms, by index, with nothing reinterpreted on the way.

    Atom order in a PSF is the atom order in the DCD (core/topology.py), so the
    numbers mean the same thing to cpptraj as they do to MDAnalysis.

    `max_ranges` refuses a mask that would be a mile of text: `name OH2` over
    236,838 waters compresses to 236,838 one-atom ranges, which is not a mask,
    it is a denial-of-service on cpptraj's parser. Such a selection has to be
    expressed by NAME (`:TIP3@OH2`), and the error says so.
    """
    idx = np.unique(np.asarray(indices, dtype=np.int64)) + 1
    if not idx.size:
        raise CpptrajError("cannot build a cpptraj mask from an empty "
                           "selection: an empty selection is a mistake, not a "
                           "result")
    brk = np.nonzero(np.diff(idx) != 1)[0]
    starts = np.concatenate(([idx[0]], idx[brk + 1]))
    ends = np.concatenate((idx[brk], [idx[-1]]))
    if len(starts) > max_ranges:
        raise CpptrajError(
            f"this selection is {len(starts):,} disjoint index ranges over "
            f"{idx.size:,} atoms; an atom-number mask that long is unusable. "
            f"Express it by residue/atom NAME instead (e.g. ':TIP3@OH2' for "
            f"water oxygen) -- see stages/solvent/rdf.py for the named masks.")
    parts = [f"{a}" if a == b else f"{a}-{b}"
             for a, b in zip(starts.tolist(), ends.tolist())]
    return "@" + ",".join(parts)


def heavy(mask):
    """`mask` restricted to heavy atoms, BY NAME.

    Never by mass: these topologies carry broken masses, which is why
    core/selections.py defines HEAVY as `not name H*` and not `mass > 2`. Every
    hydrogen in CHARMM36 is named H*, and no heavy atom in these systems is.
    """
    return f"({mask})&(!@H*)"


def mask_for(universe, name, selection, max_ranges=200):
    """(cpptraj mask, Resolved) for one of the spec's VMD-syntax selections.

    Goes through core/selections.resolve, so an empty selection raises
    EmptySelection here -- before a trajectory is opened -- instead of becoming
    a `Warning: ... has no atoms` that cpptraj shrugs off with exit 0. The
    indices are then taken from the TRANSLATED string that resolve() reports, so
    the atoms cpptraj is given are the atoms that were checked and counted, not
    a second interpretation of the same text.
    """
    from . import selections as S
    res = S.resolve(universe, name, selection)
    ag = universe.select_atoms(res.mdanalysis)
    return mask_from_indices(ag.indices, max_ranges=max_ranges), res


# ------------------------------------------------------------------- trajin
def trajin_lines(paths, stride=1, max_frames=None):
    """`trajin <file> 1 last <stride>` per segment -- the whole point.

    The stride is applied by the READER, so the frames that are skipped are
    never transferred: the cost of a pass is (frames kept) x (frame size), not
    the file size. On a 23 GB solvated replica at stride 200 that is the
    difference between ~120 MB and all of it.

    THE STRIDE RESTARTS AT EVERY SEGMENT, and it is worth knowing before
    counting frames. Each `trajin` line strides its own file from ITS frame 1,
    so ten restart segments at stride 200 give 17 frames from 1F4's 2,525, not
    the 13 a single uniform stride over the concatenation would give -- the
    short segments each contribute their first frame. The extra frames are real
    frames at known times, so nothing is wrong with them; they simply mean the
    sample is slightly denser near a restart. The frame count that ends up in
    every row is the one cpptraj reports having processed, never this estimate.

    `max_frames` caps the FIRST segment only, and deliberately: it exists so a
    test or a pilot costs seconds, and a cap spread over ten segments would
    quietly change which part of the run was sampled. A capped run is a probe,
    not a measurement, and the row it writes says so.
    """
    stride = max(1, int(stride or 1))
    out = []
    for i, p in enumerate(paths):
        if max_frames and i == 0:
            last = 1 + (int(max_frames) - 1) * stride
            out.append(f"trajin {p} 1 {last} {stride}")
        elif max_frames:
            break
        else:
            out.append(f"trajin {p} 1 last {stride}")
    return out


# ---------------------------------------------------------------------- run
def _scan(log_text):
    """Errors, warnings, the banner and the frame counts, out of one log.

    Progress is written with carriage returns (` 0% 7% ...`), so the log is
    split on both line endings or the whole bar arrives as one enormous line
    and every `grep` on it matches.
    """
    lines = log_text.replace("\r", "\n").split("\n")
    errors, warnings = [], []
    version = flavour = None
    threads = frames_read = frames_processed = None
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if version is None:
            m = _BANNER.search(s)
            if m:
                version = "V" + m.group(1)
                flavour = (m.group(3) or "serial").strip() or "serial"
        m = _THREADS.search(s)
        if m:
            threads = int(m.group(1))
        m = _FRAMES.search(s)
        if m:
            frames_read = int(m.group(1))
            frames_processed = int(m.group(2))
        if s.startswith("Error:"):
            errors.append(s)
        elif s.startswith("Warning:") or "Warning:" in s[:12]:
            if any(b in s for b in BENIGN_WARNINGS):
                continue
            if any(f in s for f in FATAL_WARNINGS):
                errors.append(s)
            else:
                warnings.append(s)
        elif any(f in s for f in ("corresponds to 0 atoms",)):
            # Not prefixed with `Warning:` at all -- `grid` prints it plain.
            errors.append(s)
    return (errors, warnings, version, flavour, threads,
            frames_read, frames_processed)


def run(lines, workdir, expect=(), omp=False, threads=None, name="cpptraj",
        timeout=None):
    """Write a cpptraj script, run it, and return a Result -- or raise.

    `expect` is the list of files the script is supposed to have written. They
    are checked for existence and for a non-zero size, because a cpptraj that
    skipped an action leaves either no file or an empty one and still exits 0.

    The script and the log are kept beside the outputs, not in a temporary
    directory: when a number looks wrong, the first question is what cpptraj was
    actually asked, and the second is what it said.
    """
    os.makedirs(workdir, exist_ok=True)
    script = os.path.join(workdir, f"{name}.in")
    log = os.path.join(workdir, f"{name}.log")
    body = "\n".join(list(lines) + ["run", "quit", ""])
    with open(script, "w") as fh:
        fh.write(body)

    if omp and not (OMP_BIN and os.path.exists(OMP_BIN)):
        omp = False                      # no OpenMP build here: the serial one
    binary = OMP_BIN if omp else SERIAL_BIN
    if not binary or not os.path.exists(binary):
        raise CpptrajError(
            f"no cpptraj found{' at ' + binary if binary else ''}: install "
            f"AmberTools (e.g. `conda install -c conda-forge ambertools`) or "
            f"export CPPTRAJ=/path/to/cpptraj")

    env = dict(os.environ)
    if omp:
        # threads.py pinned this to 1 for BLAS. The child needs its own value,
        # or the OpenMP build runs single-threaded and the manifest lies.
        env["OMP_NUM_THREADS"] = str(int(threads or OMP_THREADS))
        cmd = (["bash", "-c",
                f'source "{OMP_SETUP}" && exec "{binary}" -i "{script}"']
               if OMP_SETUP and os.path.exists(OMP_SETUP)
               else [binary, "-i", script])
    else:
        env["OMP_NUM_THREADS"] = "1"
        cmd = [binary, "-i", script]

    # A hung NFS read blocks a cpptraj forever, and with no ceiling the whole
    # folder's run waits with it. CPPTRAJ_TIMEOUT_S sets one for every call
    # that does not pass its own; unset, there is none (an rdf can legitimately
    # take days).
    if timeout is None and os.environ.get("CPPTRAJ_TIMEOUT_S"):
        timeout = float(os.environ["CPPTRAJ_TIMEOUT_S"])
    t0 = time.time()
    try:
        with open(log, "w") as fh:
            proc = subprocess.run(cmd, cwd=workdir, stdout=fh,
                                  stderr=subprocess.STDOUT, env=env,
                                  timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise CpptrajError(
            f"cpptraj did not finish {os.path.basename(script)} within "
            f"{timeout:.0f} s and was stopped.\n  script: {script}\n  log:    {log}"
        ) from exc
    seconds = time.time() - t0

    with open(log, errors="ignore") as fh:
        text = fh.read()
    errors, warnings, version, flavour, nthreads, nread, ndone = _scan(text)

    missing = []
    for p in expect:
        full = p if os.path.isabs(p) else os.path.join(workdir, p)
        if not os.path.exists(full):
            missing.append(f"{p}: never written")
        elif os.path.getsize(full) == 0:
            missing.append(f"{p}: written empty")

    if proc.returncode != 0 or errors or missing:
        detail = "\n  ".join(errors + missing) or "(no message in the log)"
        raise CpptrajError(
            f"cpptraj {version or '?'} ({flavour or 'serial'}) failed on "
            f"{os.path.basename(script)} (exit {proc.returncode}):\n  "
            f"{detail}\n  script: {script}\n  log:    {log}")
    if version is None:
        # No banner means this was not cpptraj, or its output was swallowed.
        raise CpptrajError(
            f"no cpptraj banner in {log}: the version that ran cannot be "
            f"established, so the result cannot be recorded")

    return Result(binary=binary, version=version, flavour=flavour,
                  threads=nthreads, script=script, log=log,
                  outputs=[os.path.basename(p) for p in expect],
                  errors=errors, warnings=warnings, frames_read=nread,
                  frames_processed=ndone, seconds=seconds,
                  returncode=proc.returncode)


def workdir_for(out_dir, stage, md_id):
    """Where a stage's cpptraj scripts, logs and raw outputs live.

    Beside `data/`, not inside it: `data/` is the tidy contract that pages and
    aggregation read, and a directory of cpptraj `.dat` files in it would be
    picked up by core.io.read as if it were results.
    """
    return os.path.join(os.path.dirname(os.path.abspath(out_dir)),
                        "cpptraj", stage, str(md_id))


# ----------------------------------------------------------- reading it back
class Grid:
    """One OpenDX volumetric grid, as cpptraj wrote it."""

    __slots__ = ("path", "counts", "origin", "spacing", "values")

    def __init__(self, path, counts, origin, spacing, values):
        self.path = path
        self.counts = counts            # (nx, ny, nz)
        self.origin = origin            # lower corner, A
        self.spacing = spacing          # (dx, dy, dz), A
        self.values = values            # flat, in the units cpptraj wrote

    @property
    def voxel_volume(self):
        return float(self.spacing[0] * self.spacing[1] * self.spacing[2])

    @property
    def n_voxels(self):
        return int(self.values.size)


def read_dx(path):
    """Read an OpenDX grid without holding three copies of it in memory.

    The header is read line by line and the numeric body in one slice, so a
    39.5-million-voxel grid (which is what 0.5 A over this trimer costs) is one
    string and one array rather than 39.5 million Python objects.
    """
    counts = origin = None
    spacing = []
    items = None
    with open(path, errors="ignore") as fh:
        while True:
            line = fh.readline()
            if not line:
                raise CpptrajError(f"{path}: no 'data follows' -- not an "
                                   f"OpenDX grid, or truncated")
            s = line.strip()
            if s.startswith("object 1") and "counts" in s:
                counts = tuple(int(x) for x in s.split()[-3:])
            elif s.startswith("origin"):
                origin = tuple(float(x) for x in s.split()[1:4])
            elif s.startswith("delta"):
                spacing.append([float(x) for x in s.split()[1:4]])
            if "data follows" in s:
                m = re.search(r"items\s+(\d+)", s)
                items = int(m.group(1)) if m else None
                break
        body = fh.read()
    cut = body.find("\nobject")           # the trailing `object "..." class field`
    values = np.fromstring(body if cut < 0 else body[:cut], sep=" ")
    if counts is None or origin is None or len(spacing) != 3:
        raise CpptrajError(f"{path}: incomplete OpenDX header "
                           f"(counts={counts} origin={origin} "
                           f"deltas={len(spacing)})")
    if items is not None and values.size != items:
        raise CpptrajError(
            f"{path}: header promises {items:,} voxels, the file holds "
            f"{values.size:,}. A truncated grid averages to the wrong density, "
            f"so it is refused rather than used.")
    diag = (spacing[0][0], spacing[1][1], spacing[2][2])
    off = max(abs(spacing[0][1]), abs(spacing[0][2]), abs(spacing[1][0]),
              abs(spacing[1][2]), abs(spacing[2][0]), abs(spacing[2][1]))
    if off > 1e-9:
        raise CpptrajError(f"{path}: the grid axes are not orthogonal "
                           f"(off-diagonal delta {off:g}); nothing here can "
                           f"read a sheared grid")
    return Grid(path, counts, origin, diag, values)


def read_datafile(path):
    """(column names, array) for a cpptraj 'Standard Data File'.

    The header is `#Distance_(Ang) <set name>`; a set name can contain spaces
    (`:TIP3@OH2 => @1-40140&!@H*`), so the names are split on two-or-more
    spaces, never on one.
    """
    names = []
    with open(path, errors="ignore") as fh:
        first = fh.readline()
    if first.startswith("#"):
        names = [c.strip() for c in re.split(r"\s{2,}", first[1:].strip())
                 if c.strip()]
    data = np.loadtxt(path, comments="#", ndmin=2)
    if not data.size:
        raise CpptrajError(
            f"{os.path.basename(path)} holds no rows. cpptraj writes an empty "
            f"data file when an action was skipped and still exits 0, so this "
            f"is a failure, not an empty result.")
    return names, data


def mean_volume(path):
    """The mean unit-cell volume over the frames of a run, from `volume out`.

    Measured, not taken from the setup script's `-sc` argument or from a
    `.xst` beside the trajectory: the bulk density every grid and every g(r) is
    divided by comes from the frames that were actually read.
    """
    _, data = read_datafile(path)
    vols = data[:, -1]
    return float(vols.mean()), float(vols.std()), int(vols.size)


if __name__ == "__main__":                                 # pragma: no cover
    import sys
    for p in sys.argv[1:]:
        if p.endswith(".dx"):
            g = read_dx(p)
            print(f"{p}: {g.counts} voxels of {g.spacing} A at {g.origin}, "
                  f"max {g.values.max():.3f}, mean {g.values.mean():.4f}")
        else:
            names, d = read_datafile(p)
            print(f"{p}: {d.shape} {names}")
