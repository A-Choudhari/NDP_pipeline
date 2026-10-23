#!/usr/bin/env python3
"""
config.py -- substrate.yaml loading, bounded directory walking, and the safety
probes every other module depends on.

This is the module that makes the toolkit portable. The previous version of
this code hardcoded ROOT/ABZYME/APTAMERS; everything configurable now lives in
a per-user substrate.yaml.
"""

import os
import sys

# Vendor, cache and environment directories. Walking these on a 628 T NFS mount
# is slow and returns nothing anyone wants indexed.
DEFAULT_PRUNE = [
    "conda_pkgs", "site-packages", "__pycache__", ".git", "node_modules",
    ".ipynb_checkpoints", ".cache", ".conda", ".local", "envs",
    "env", "venv", ".venv", "miniconda3", "anaconda3",
    "lib", "lib64", "include", "share", "bin", "man",
    ".Trash", "toppar",
    # Vendored libraries and their bundled sample data. PyContact's
    # exampleData appeared inside an analysis tree mid-session and was counted
    # as a 10th CDR "MD run" (0.0 GB, replica unknown). 13 such directories
    # exist across this workspace alone; every colleague will have some.
    "vendor", "vendored", "third_party", "thirdparty", "external",
    "exampleData", "example_data", "examples", "sample_data", "sampledata",
    # `examples` did not catch `example-output`: 7 NAMD tutorial runs shipped
    # inside namd-tutorial-files/*/example-output were being inventoried as the
    # user's own simulations. Same defect as exampleData, one hyphen away.
    "example-output", "example_output", "exampleOutput",
    "test_data", "testdata", "tests", "testing", "demo", "demos",
    # Sibling tooling built on top of this index. MDDB_substrate is another
    # tool's OUTPUT and the pipelines are its code; indexing either
    # would make this catalogue describe its own ecosystem rather than the
    # science. Pruned in the shared default rather than per-user because every
    # colleague who adopts the MDDB layer would otherwise get an UNCOVERED
    # entry and a failing coverage check in their own folder.
    "MDDB_substrate", "NDP_PIPELINE", "MDDB_PIPELINE",
    # Our own visible output, which lives INSIDE the folder being indexed.
    # Left in, a second run would inventory the first run's records.
    "MD_RESULTS",
]

# Trajectory and topology extensions we understand.
# NOTE: `.crd` is deliberately NOT a trajectory extension. In CHARMM/CHARMM-GUI
# output it is a SINGLE-STRUCTURE coordinate file, and including it turned 220
# CHARMM-GUI build directories into phantom "AMBER runs" holding 9.7 GB between
# them. It belongs with the structure formats below.
TRAJ_EXT = (".dcd", ".nc", ".netcdf", ".xtc", ".trr", ".mdcrd")
TOPO_EXT = (".psf", ".prmtop", ".parm7", ".top", ".tpr", ".gro", ".crd", ".pdb")

TABLE_EXT = (".tsv", ".csv", ".dat")

# scipy reads AMBER NetCDF headers. It was not listed, so a Python without it
# indexed every .nc run with no atom or frame count, and only a per-file
# `error:ModuleNotFoundError` said so. Found by running the pipeline under a
# freshly built environment on 2026-09-16.
REQUIRED = ["pandas", "pyarrow", "yaml", "scipy"]
INSTALL_HINT = {
    "pandas": "pandas", "pyarrow": "pyarrow", "yaml": "pyyaml", "scipy": "scipy",
}


def check_dependencies(extra=()):
    """Fail early with an actionable message instead of a traceback.

    Colleagues run their own Python; none of this is guaranteed to be present.
    """
    missing = []
    for mod in list(REQUIRED) + list(extra):
        try:
            __import__(mod)
        except ImportError:
            missing.append(INSTALL_HINT.get(mod, mod))
    if missing:
        print("\n  Missing required packages: " + ", ".join(missing))
        print("\n  Install them with:\n")
        print(f"      {sys.executable} -m pip install " + " ".join(missing))
        print("\n  or, with conda:\n")
        print("      conda install -c conda-forge " + " ".join(missing) + "\n")
        return False
    return True


def skip_symlink(path):
    """True when `path` is a symlink and must not be read.

    `os.walk(followlinks=False)` refuses to *traverse* linked directories, which
    is often mistaken for symlink safety. A linked FILE is still listed and
    opened. On shared storage `notes.md -> /shared/<someone>/private.md`
    would otherwise be read and indexed.

    Lives here, not in one scanner, because it was originally applied only to
    the table scanner and silently missing from the document and trajectory
    scanners for a day.
    """
    try:
        return os.path.islink(path)
    except OSError:
        return True


def is_writable(path):
    """True when we could create a file in `path`.

    Used to decide whether a converted Parquet may sit next to its source.
    On this filesystem only ~15 of 67 top-level directories are writable by any
    given user, so this is the difference between working and crashing.
    """
    d = path if os.path.isdir(path) else os.path.dirname(path)
    return os.access(d, os.W_OK | os.X_OK)


def current_user():
    """The account running this, by uid.

    Not `$USER`: cron and `env -i` shells leave it unset, and a folder
    initialised that way was recorded as `owner: unknown` for good, so every
    attribution downstream fell back to nobody.
    """
    try:
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name
    except (KeyError, ImportError):
        return os.environ.get("USER") or "unknown"


def default_config(root, owner=None):
    return {
        "owner": owner or current_user(),
        "root": os.path.abspath(root),
        "topics": [],
        "corpora": [],
        "prune": list(DEFAULT_PRUNE),
        "min_table_bytes": 50_000,
        "deep_scan": False,
        # Never write into the data directories. Outputs go to SUBSTRATE/tables/.
        "write_in_place": False,
    }


def load(path):
    """Load and validate an old-layout substrate.yaml. Returns the config dict.

    The current layout keeps the same mapping in `settings.yaml`'s `folder:`
    section; `lib.settings.config_of` validates it through `validate()` below.
    This entry point stays for reading a folder that has not been migrated yet.
    """
    import yaml
    with open(path) as fh:
        cfg = yaml.safe_load(fh) or {}
    return validate(cfg, path)


def validate(cfg, path, folder=None):
    """Fill in the defaults, check what must be right, and mark what not to scan.

    `path` only names the file in error messages. `folder` is the data folder
    whose output directories must be excluded; it defaults to the directory the
    settings live under, which is where the old layout kept them.
    """
    cfg = dict(cfg or {})
    if not cfg.get("root"):
        raise ValueError(f"{path}: 'root' is required")
    cfg.setdefault("owner", current_user())
    cfg.setdefault("topics", [])
    cfg.setdefault("corpora", [])
    cfg.setdefault("prune", list(DEFAULT_PRUNE))
    cfg.setdefault("min_table_bytes", 50_000)
    cfg.setdefault("deep_scan", False)
    cfg.setdefault("write_in_place", False)
    cfg.setdefault("systems", {})
    cfg.setdefault("conditions", {})
    cfg.setdefault("systems_from", [])
    cfg.setdefault("conditions_from", [])
    cfg.setdefault("stages", {})
    for key in ("systems", "conditions", "stages"):
        if not isinstance(cfg[key], dict):
            raise ValueError(f"{path}: '{key}' must be a mapping of "
                             f"label -> pattern (or list of patterns)")
    for key in ("systems_from", "conditions_from"):
        if isinstance(cfg[key], str):
            cfg[key] = [cfg[key]]
        if not isinstance(cfg[key], list):
            raise ValueError(f"{path}: '{key}' must be a regex, or a list of "
                             f"regexes, each with one capture group")
    cfg["root"] = os.path.abspath(os.path.expanduser(cfg["root"]))
    if not os.path.isdir(cfg["root"]):
        raise ValueError(f"{path}: root {cfg['root']} is not a directory")

    # NEVER SCAN OUR OWN OUTPUT. The substrate lives inside the root it
    # describes, so a whole-root scan would index AGENTS.md, README.md,
    # LEGEND.md and the reports as if they were the user's science -- and the
    # document count would then grow on every build, from its own output. The
    # toolkit directory is excluded for the same reason.
    #
    # This is a PATH exclusion, deliberately not a name one: `ABZYME/SUBSTRATE`
    # is a real hand-built claims ledger in this workspace, and pruning every
    # directory called SUBSTRATE would silently drop it.
    from . import layout                                   # noqa: PLC0415
    toolkit = os.path.realpath(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    own = folder or os.path.dirname(os.path.abspath(path))
    cfg["exclude_paths"] = [
        os.path.realpath(os.path.dirname(os.path.abspath(path))),
        toolkit,
        # MD_RESULTS/ and .md_pipeline/, and the SUBSTRATE/ and MDDB_substrate/
        # a folder may still hold from the layout before them. An explicit
        # `prune:` REPLACES the defaults, so naming them by PATH here is what
        # keeps a colleague's older settings from indexing our own output --
        # and, deliberately not by NAME, because `ABZYME/SUBSTRATE` is a real
        # hand-built ledger in this workspace.
        *[os.path.realpath(p) for p in layout.excluded(own)],
        *[os.path.realpath(p) for p in layout.excluded(cfg["root"])],
        # The other pipeline, wherever it is installed beside this one.
        os.path.join(os.path.dirname(toolkit), "MDDB_PIPELINE"),
    ]
    return cfg


LABEL_SCAFFOLD = """
# --- what each run IS (optional) --------------------------------------------
# `system` and `condition` are `unknown` for every run until you say what your
# paths mean. Nothing is ever guessed: a label appears only where a rule you
# wrote matches. Uncomment and edit.
#
# Matched case-insensitively against each run's path relative to `root`, so a
# rule can key on a parent directory -- which is usually where the system lives.
# The FIRST declared match wins, so put specific rules above general ones.
#
# systems:
#   gp120-open:   ["3-OCCLUDED-OPEN", "4-OPEN"]   # one label, several patterns
#   gp120-deglyc: "2-NO-GLYC"
#
# conditions:
#   apo:  "/apo(/|$)"
#   holo: "/holo(/|$)"
#
# When the directory NAME is the thing simulated -- one folder per sequence,
# say -- take the label from the path instead of writing one rule per folder.
# The capture group becomes the label, used verbatim:
#
# systems_from:
#   - "GENTILE/[0-9]+-([ACGUT]{6})"     # 10-TGGGGC -> TGGGGC
#
# An explicit `systems:` entry beats a captured one, which is how you correct a
# single directory whose name lies. A pattern that cannot compile stops the
# build and names itself, rather than silently labelling nothing.
#
# After editing, re-run substrate_build.py and check it took:
#   python3 .../ask.py sql "SELECT system, COUNT(*) FROM runs GROUP BY system"
"""


STAGE_SCAFFOLD = """
# --- what each trajectory file IS (optional) ---------------------------------
# Every trajectory file is given a stage from its name and its directory:
# `eq`, `heat`, `min` -> not production; `prod`, `pr`, `md` -> production;
# `stripped`, `-old`, `frames` -> not the simulation itself.
#
# A .dcd whose name says nothing (`base_model.R1.dcd`, `all.dcd`, `4092-5091.dcd`)
# is counted as PRODUCTION BY DEFAULT when its directory holds no other
# production file. Until you confirm it, every MD built from such a file is in
# state `review` with the issue `stage_assumed`, and SUBSTRATE/reports/
# inventory.md lists each one. Other formats whose name says nothing stay
# `unresolved` and are not counted.
#
# Say what a file really is here, and that wins over both. Same shape as
# `systems:`: a stage, then a pattern (or a list), matched case-insensitively
# against the FILE's path relative to `root`. The first declared match wins.
# Stages: production, equilibration, heating, minimisation, merged, stripped,
# superseded, derived.
#
# stages:
#   derived:    ['/RMSF/all\\.dcd$', '/benchmark/', '_movie\\.dcd$']
#   stripped:   '_prot\\.dcd$'
#   production: 'base_model\\..*\\.R[0-9]+\\.dcd$'
#
# Write patterns in single quotes: YAML reads a backslash inside double quotes
# as an escape, and '\\.' (a literal dot) is then refused.
#
# After editing, re-run substrate_build.py, then mddb_build.py.
"""

#: What an author may declare a trajectory to be under `stages:`.
DECLARABLE_STAGES = ("production", "equilibration", "heating", "minimisation",
                     "merged", "stripped", "superseded", "derived")


def compile_stage_rules(cfg):
    """Compile the optional `stages:` rules: the author saying what a file is.

    Returns [(stage, [compiled, ...]), ...] in declaration order. Same shape
    as `systems:`, but matched against each trajectory FILE's path relative to
    root, not the directory's, because two files side by side can be different
    things. A stage name the classifier does not know stops the build: a typo
    (`production `, `derivd`) must not quietly leave a file production.
    """
    import re
    rules = []
    for stage, pats in (cfg.get("stages") or {}).items():
        stage = str(stage).strip()
        if stage not in DECLARABLE_STAGES:
            raise ValueError(
                f"substrate.yaml: stages: {stage!r} is not a stage. Use one of: "
                f"{', '.join(DECLARABLE_STAGES)}")
        if isinstance(pats, str):
            pats = [pats]
        compiled = []
        for pat in pats or []:
            try:
                compiled.append(re.compile(str(pat), re.I))
            except re.error as exc:
                raise ValueError(
                    f"substrate.yaml: stages: {stage}: bad regex {pat!r}: {exc}")
        if compiled:
            rules.append((stage, compiled))
    return rules


def folder_section(cfg, indent="  "):
    """The `folder:` section of settings.yaml, as text.

    Same mapping `save()` has always written, one level in, with the same
    scaffolds: a person reads this file before any documentation, so the
    optional features have to be visible in it.
    """
    import yaml
    body = {k: v for k, v in cfg.items()
            if k not in ("systems", "conditions", "systems_from",
                         "conditions_from", "stages", "exclude_paths")
            or (v and k != "exclude_paths")}
    dumped = yaml.safe_dump(body, sort_keys=False, default_flow_style=False)
    text = "folder:\n" + "".join(f"{indent}{line}\n" if line.strip() else "\n"
                                 for line in dumped.splitlines())
    scaffold = ""
    if not (cfg.get("systems") or cfg.get("systems_from")):
        scaffold += LABEL_SCAFFOLD
    if not cfg.get("stages"):
        scaffold += STAGE_SCAFFOLD
    # The scaffolds are comment blocks; indent them so they sit inside the
    # section they describe and a person can uncomment a line in place.
    text += "".join(f"{indent}{line}\n" if line.strip() else "\n"
                    for line in scaffold.splitlines())
    return text


def save(cfg, path):
    import yaml
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # Write the labelling scaffold as COMMENTS rather than empty keys. This file
    # is where people look before they read any documentation, and a generated
    # config that never mentions `systems:` is a feature nobody discovers -- the
    # first version of this toolkit documented the feature in AGENTS.md while
    # the config said nothing and the code did not implement it at all.
    body = {k: v for k, v in cfg.items()
            if k not in ("systems", "conditions", "systems_from",
                         "conditions_from", "stages", "exclude_paths")
            or (v and k != "exclude_paths")}
    with open(path, "w") as fh:
        fh.write("# substrate.yaml -- generated by substrate_init.py\n"
                 "# Edit freely; re-run substrate_build.py after any change.\n\n")
        yaml.safe_dump(body, fh, sort_keys=False, default_flow_style=False)
        if not (cfg.get("systems") or cfg.get("systems_from")):
            fh.write(LABEL_SCAFFOLD)
        if not cfg.get("stages"):
            fh.write(STAGE_SCAFFOLD)


def compile_label_rules(cfg, key):
    """Compile the optional `systems:` / `conditions:` rules from substrate.yaml.

    Returns [(label, [compiled, ...]), ...] in declaration order, or [] when the
    user declared none -- in which case the field stays `unknown` for every run.

    The shape is deliberately the simplest thing that describes a directory
    layout:

        systems:
          gp120-open:      "3-OCCLUDED-OPEN|4-OPEN"
          gp120-glycan-free: "2-NO-GLYC"
        conditions:
          apo: "/apo(/|$)"

    A label maps to one regex or a list of them, matched case-insensitively
    with `search` against the run's path RELATIVE TO root -- so a rule can key
    on a parent directory, which is usually where the system is encoded.

    Nothing is inferred without a rule. This function only makes the user's own
    declaration executable; it never proposes labels of its own, because a
    guessed `system` is exactly the kind of value that gets quoted in a paper.
    """
    import re
    rules = []
    for label, pats in (cfg.get(key) or {}).items():
        if isinstance(pats, str):
            pats = [pats]
        compiled = []
        for pat in pats:
            try:
                compiled.append(re.compile(pat, re.I))
            except re.error as exc:
                # A bad regex must name itself. Silently dropping it would
                # leave every run `unknown` with no indication why.
                raise ValueError(
                    f"substrate.yaml: {key}: {label}: bad regex {pat!r}: {exc}")
        if compiled:
            rules.append((str(label), compiled))
    return rules


def compile_capture_rules(cfg, key):
    """Compile the optional `systems_from:` / `conditions_from:` rules.

    Where `systems:` maps a label you wrote to a pattern, these take the label
    FROM THE PATH ITSELF via a capture group -- for a corpus whose directories
    are named after the thing they simulate, and where writing out one rule per
    directory would be fifty lines that rot the moment a directory is added:

        systems_from: "GENTILE/\\d+-([ACGUT]{6})"

    Returns a list of compiled patterns in declaration order. The captured text
    is used verbatim: no case-folding, no cleanup. A transformed label would no
    longer match the directory it came from, which is the one property that
    makes a captured label checkable.
    """
    import re
    pats = cfg.get(key) or []
    if isinstance(pats, str):
        pats = [pats]
    compiled = []
    for pat in pats:
        try:
            rx = re.compile(pat, re.I)
        except re.error as exc:
            raise ValueError(
                f"substrate.yaml: {key}: bad regex {pat!r}: {exc}")
        if rx.groups < 1:
            # Without a group there is nothing to name the label with, and the
            # rule would silently label nothing at all.
            raise ValueError(
                f"substrate.yaml: {key}: {pat!r} has no capture group -- "
                f"put ( ) around the part that names the system, or use "
                f"`{key.replace('_from', '')}:` to map a fixed label instead")
        compiled.append(rx)
    return compiled


def compile_replica_rules(cfg):
    """Compile the optional `folder: replicas:` rules: the author saying which
    trajectory FILES in one directory are side-by-side replicas.

        replicas:
          - '_run_(\\d+)\\.dcd$'     # the capture group is the replica number

    Matched against each production file's NAME. Absent, this returns [] and
    replica splitting is exactly what it was before the key existed. Read here
    rather than in validate(), so a folder that never declares it never has an
    empty `replicas:` written back into its settings.yaml.
    """
    import re
    pats = cfg.get("replicas") or []
    if isinstance(pats, str):
        pats = [pats]
    if not isinstance(pats, list):
        raise ValueError("settings.yaml: folder: replicas: must be a regex, or a "
                         "list of regexes, each with one capture group")
    compiled = []
    for pat in pats:
        try:
            rx = re.compile(str(pat), re.I)
        except re.error as exc:
            raise ValueError(
                f"settings.yaml: folder: replicas: bad regex {pat!r}: {exc}")
        if rx.groups < 1:
            raise ValueError(
                f"settings.yaml: folder: replicas: {pat!r} has no capture group -- "
                f"put ( ) around the replica number, e.g. '_(\\d+)\\.dcd$'")
        compiled.append(rx)
    return compiled


def corpus_roots(cfg):
    """Absolute (name, path) pairs to scan.

    With no `corpora` declared the whole root is scanned as one corpus -- so the
    tool does something useful on the very first run with zero configuration.
    """
    root = cfg["root"]
    if not cfg.get("corpora"):
        return [(os.path.basename(root) or "root", root)]
    out = []
    for c in cfg["corpora"]:
        if isinstance(c, str):
            out.append((c, os.path.join(root, c)))
        else:
            p = c.get("path", "")
            out.append((c.get("name") or os.path.basename(p) or "root",
                        p if os.path.isabs(p) else os.path.join(root, p)))
    return [(n, p) for n, p in out if os.path.isdir(p)]



def empty_frame(spec):
    """An empty DataFrame with the right columns AND the right dtypes.

    Columns alone are not enough, and the failure is silent. On an empty frame
    built from column names only, a bool column arrives as dtype `object`, so
    `df[~df.curated]` is no longer a boolean mask -- pandas reads the object
    Series as a list of COLUMN LABELS and returns a zero-column frame. The next
    attribute access then dies with "no attribute 'status'", nowhere near the
    cause. Fixing this for runs.parquet and not for docs, doc_refs and tables is
    what broke a colleague's first build.
    """
    import pandas as _pd
    return _pd.DataFrame({c: _pd.Series(dtype=t) for c, t in spec.items()})

def coverage(cfg):
    """Classify every top-level entry under `root`: is it actually indexed?

    Declared corpora give good per-project attribution, but they fail SILENTLY:
    a new project directory is simply absent from the registry, and nothing in
    the output says so. That is how AROMATIC_RINGS_ENERGIES_PROJECT (37 docs)
    and a 250 ns trajectory under DOWNLOAD/ stayed invisible -- the substrate
    looked complete because it only ever reported what it had scanned.

    Returns a list of (name, status, detail), status being:
      indexed  -- inside a declared corpus
      pruned   -- excluded by the prune list, by name
      output   -- the substrate's or toolkit's own directory
      file     -- a loose file at the root; corpora are directories
      UNCOVERED-- none of the above. This is the one that matters.
    """
    root = cfg["root"]
    prune = set(cfg.get("prune") or [])
    excl = {os.path.realpath(e) for e in (cfg.get("exclude_paths") or ())}
    tops = []
    for c in (cfg.get("corpora") or []):
        pth = c if isinstance(c, str) else (c.get("path") or "")
        first = str(pth).strip("/").split("/")[0]
        if first:
            tops.append(first)
    covered = set(tops)

    rows = []
    try:
        entries = sorted(os.listdir(root))
    except OSError:
        return rows
    for name in entries:
        full = os.path.join(root, name)
        if not os.path.isdir(full):
            rows.append((name, "file", "loose trajectory at the root: indexed"
                         if name.endswith(TRAJ_EXT) else "not a directory"))
        elif os.path.realpath(full) in excl:
            rows.append((name, "output", "this tool's own output"))
        elif name in covered:
            rows.append((name, "indexed", ""))
        elif name in prune or name.startswith("."):
            rows.append((name, "pruned", "excluded by `prune:`"))
        else:
            rows.append((name, "UNCOVERED", "not in any corpus and not pruned"))
    return rows


def walk(base, prune, follow_symlinks=False, max_depth=None, exclude=()):
    """Bounded os.walk that prunes vendor trees and never follows symlinks.

    Symlinks are not followed on purpose: this filesystem has cross-user links,
    and following them would index someone else's data into your registry.

    `exclude` is a list of absolute paths -- the substrate's own directory and
    the toolkit's -- matched by resolved path, never by name. Comparing names
    would drop `ABZYME/SUBSTRATE`, a real ledger that happens to share a name
    with the output directory. The realpath call is made only for a directory
    whose basename could match, so the common case costs one set lookup.
    """
    pruneset = set(prune)
    excl = {os.path.realpath(e) for e in exclude}
    excl_names = {os.path.basename(e) for e in excl}
    for dirpath, dirnames, filenames in os.walk(base, topdown=True,
                                                followlinks=follow_symlinks):
        dirnames[:] = [d for d in dirnames
                       if d not in pruneset and not d.startswith(".")
                       and not (d in excl_names
                                and os.path.realpath(os.path.join(dirpath, d)) in excl)]
        # Depth by path components, not by counting "/": on Windows the
        # separator is "\" and the limit never applied.
        rel = os.path.relpath(dirpath, base)
        depth = 0 if rel == os.curdir else len(rel.split(os.sep))
        if max_depth is not None and depth >= max_depth:
            dirnames[:] = []
        yield dirpath, dirnames, filenames


def rel(path, cfg):
    """Path relative to the configured root, for stable IDs across machines."""
    try:
        return os.path.relpath(path, cfg["root"])
    except ValueError:
        return path
