#!/usr/bin/env python3
"""
workdir.py -- give the MDDB workflow somewhere to write that is not your data.

THE PROBLEM. The installed workflow (mwf 0.1.10) treats each MD's `mdir` as a
scratch directory. It creates it, writes `.mwf_cache.json`, `.register.json`,
`structure.pdb`, `trajectory.xtc`, `metadata.json`, an `inpro/` directory of
converted trajectories and every `mda.*.json` analysis into it
(mwf.py:146-215, constants.py:132-228) -- and at the end of the run it DELETES
`time.txt`, `#*` and `f*.trajectory.xtc` from it (tools/remove_trash.py:6-16).
Every export used to set `mdir` to the production run directory itself, so
processing any dataset as generated would have written into, and deleted from,
raw simulation data. `MDDB_substrate/smoketest/` shows the artefacts a real run
leaves behind; that smoke test only stayed clean because someone built its
directory by hand out of symlinks. This module builds that layout for every MD.

THE LAYOUT.

    MDDB_substrate/
      exports/<dataset_id>/inputs.yaml      pristine, generated, never run in place
      exports/<dataset_id>/provenance.yaml  where every input really lives
      work/<dataset_id>/                    mwf -dir points here
        inputs.yaml                         the copy mwf is allowed to rewrite
        src.<topology>          -> source   (project-level topology, if shared)
        <md_id>/                            one mdir per MD, a single component
          src.<trajectory>      -> source
          src.<topology>        -> source   (per-MD topology, if they differ)

WHY EACH DETAIL IS THE WAY IT IS.

  `mdir` is a single relative component. The workflow creates it with `mkdir`,
  not `makedirs` (mwf.py:154), and an absolute `mdir` is silently re-rooted by
  string concatenation in Project.pathify (mwf.py:1688-1690).

  Inputs are SYMLINKS, never copies. There are ~6.4 TB of production trajectory
  here; copying to satisfy a relative-path convention would double it. The
  upstream docs recommend exactly this (docs/source/usage.rst:24-26).

  Every link is named `src.<original name>`. Writing to a path that is a
  symlink writes THROUGH the link into its target. If a source topology happened
  to be called `structure.pdb` or `topology.psf`, a link with that bare name
  would sit exactly where mwf writes its own output, and the raw file would be
  overwritten. No mwf output name begins with `src.` (constants.py:125-232), and
  RESERVED_NAMES below is checked anyway.

  The workflow rewrites inputs.yaml in place and strips its comments
  (core/inputs.py:172-193). So the export stays pristine and the work directory
  gets a copy to consume.

WHAT THIS DOES NOT PROMISE. Every check here is application-level: it stops
this toolkit from pointing the workflow at raw data. It is not filesystem
enforcement. Only read-only permissions or a read-only mount on the source tree
can guarantee that nothing else writes there.
"""

import hashlib
import os
import shutil

#: Prefix on every symlink this module creates. It is also how pruning knows a
#: file is ours to remove.
LINK_PREFIX = "src."

#: Names the workflow writes into a project or MD directory. A link must never
#: take one of these names, because writing to it would write through the link.
RESERVED_NAMES = frozenset({
    "inputs.yaml", "inputs.yml", "inputs.json",
    ".mwf_cache.json", ".aux.mwf_cache.json", ".register.json",
    "structure.pdb", "trajectory.xtc", "trajectory_SAMPLE.xtc",
    "topology.json", "topology.psf", "topology.prmtop", "topology.top",
    "topology.tpr", "metadata.json", "chains.json", "charges.txt",
    "resorted_charges.json", "resorted_bonds.json", "first_frame.pdb",
    "average.pdb", "pdb_references.json", "protein_references.json",
    "inchikey_references.json", "time.txt", "atommass.dat", "filter.ndx",
})

#: Name of the fingerprint file recording which export the work copy came from.
EXPORT_STAMP = ".export_fingerprint"


class UnsafeWorkdir(Exception):
    """A working directory would resolve somewhere the workflow must not write."""


def within(path, root):
    """True when `path` is `root` or lies beneath it, component by component.

    String-prefix tests are wrong here: '/data/run' is a prefix of
    '/data/run2', and '/data/run2' is not inside '/data/run'. Both sides are
    compared after full symlink resolution, so a path that reaches the root
    through a link is judged by where it actually lands.
    """
    p = os.path.realpath(path)
    r = os.path.realpath(root)
    try:
        return os.path.commonpath([p, r]) == r
    except ValueError:            # different drives, or one relative and one not
        return False


def authorised_root(mddb):
    """The one directory the workflow may write under: <MDDB_substrate>/work.

    The MDDB_substrate directory is resolved, but `work` is appended LEXICALLY.
    If `work` were a symlink into raw data, resolving it too would make every
    path beneath it look authorised. Appending the name without resolving it
    means such a link is caught by `within` instead.
    """
    return os.path.join(os.path.realpath(mddb), "work")


def assert_safe(md_dir, mddb, protected):
    """Raise UnsafeWorkdir unless `md_dir` is inside the authorised root and
    outside every protected (source) directory.

    Checked on the resolved path, so an escape through a symlink anywhere in the
    chain -- `work` itself, the dataset directory, the MD directory -- is caught.
    """
    root = authorised_root(mddb)
    if os.path.islink(root):
        raise UnsafeWorkdir(
            f"{root} is a symlink (to {os.path.realpath(root)}). The workflow's "
            f"working area must be a real directory inside MDDB_substrate.")
    if not within(md_dir, root):
        raise UnsafeWorkdir(
            f"working directory {md_dir} resolves to {os.path.realpath(md_dir)}, "
            f"which is outside {root}.")
    for src in protected:
        if src and within(md_dir, src):
            # A run lying loose at the folder root has the ROOT as its source
            # directory, and MDDB_substrate/ lives in the root, so every working
            # directory is "inside" it. That nesting is the only way this can
            # hold while md_dir is under the (real, unlinked) authorised root:
            # the work area is then our own output, which the index never scans,
            # and none of the run's files can be in it. Anything else still stops.
            if within(root, src) and not within(src, root):
                continue
            raise UnsafeWorkdir(
                f"working directory {md_dir} resolves inside source directory "
                f"{src}. The workflow writes into and deletes from its working "
                f"directories; that must never be simulation data.")


#: Extensions that tools read as NetCDF. A NetCDF file named anything else gets
#: a link ending in `.nc`, because mwf and MDAnalysis both choose a reader from
#: the extension.
NETCDF_EXTS = (".nc", ".netcdf", ".cdf", ".ncdf")


def corrected_ext(source_path):
    """The extension a link should carry so tools read the file as what it is,
    or None when the source name is already truthful."""
    if str(source_path).lower().endswith(NETCDF_EXTS):
        return None
    if not str(source_path).lower().endswith((".mdcrd", ".crd", ".trj")):
        return None
    try:
        from core.trajscan import sniff_format
    except ImportError:
        return None
    return ".nc" if sniff_format(source_path) == "netcdf" else None


def link_name(source_path, ext=None):
    """The link name for a source file: `src.<basename>`, never a reserved name.

    With `ext`, the source extension is replaced -- used only when the file's
    content contradicts its name (see corrected_ext).
    """
    base = os.path.basename(str(source_path))
    if ext:
        base = os.path.splitext(base)[0] + ext
    name = LINK_PREFIX + base
    if name in RESERVED_NAMES:            # impossible with the prefix; belt and braces
        raise UnsafeWorkdir(f"link name {name!r} collides with a workflow output")
    return name


def place_link(directory, source_path, ext=None):
    """Create or refresh one symlink. Returns the link name.

    Idempotent: an existing link to the same target is left alone, and a link
    to a different target is repointed (the source moved, or a better topology
    was chosen). A REAL file with that name is never replaced -- it could be
    workflow output, or something a person put there.
    """
    name = link_name(source_path, ext)
    path = os.path.join(directory, name)
    target = os.path.abspath(source_path)
    if os.path.islink(path):
        if os.readlink(path) == target:
            return name
        os.unlink(path)
    elif os.path.lexists(path):
        raise UnsafeWorkdir(
            f"{path} exists and is not a symlink; refusing to replace it.")
    os.symlink(target, path)
    return name


def _fingerprint(path):
    with open(path, "rb") as fh:
        return hashlib.blake2b(fh.read(), digest_size=16).hexdigest()


def materialise(mddb, dataset_id, doc, source_root):
    """Build work/<dataset_id>/ for one exported dataset.

    `doc` is the SOURCE view produced by mwf.build_dataset: each MD's `mdir` and
    topology are paths relative to `source_root`. Returns
    (runnable_doc, provenance) -- the doc rewritten to point at the links, and a
    record of what every link stands for.
    """
    project = os.path.join(authorised_root(mddb), dataset_id)
    mds = doc.get("mds") or []
    protected = {os.path.join(source_root, str(m["mdir"])) for m in mds}

    # Check every destination BEFORE creating anything, so an unsafe layout
    # fails with nothing half-built. realpath resolves whatever part of a path
    # already exists, so this also catches a pre-existing link on the way.
    assert_safe(project, mddb, protected)
    for m in mds:
        assert_safe(os.path.join(project, str(m["name"])), mddb, protected)
    os.makedirs(project, exist_ok=True)

    runnable = {k: v for k, v in doc.items() if not k.startswith("_")}
    provenance = {"source_root": os.path.realpath(source_root),
                  "work_dir": os.path.relpath(project, os.path.dirname(mddb)),
                  "mds": []}

    top = doc.get("input_topology_filepath")
    if top:
        src_top = os.path.join(source_root, top)
        runnable["input_topology_filepath"] = place_link(project, src_top)
        provenance["topology"] = top

    ordering = doc.get("_ordering") or {}
    split = doc.get("_split") or {}
    assumed = doc.get("_assumed") or {}
    out_mds = []
    for m in mds:
        md_id, rel = str(m["name"]), str(m["mdir"])
        md_dir = os.path.join(project, md_id)
        os.makedirs(md_dir, exist_ok=True)
        assert_safe(md_dir, mddb, protected)

        trajs, relabelled = [], {}
        for t in m.get("input_trajectory_filepaths") or []:
            src = os.path.join(source_root, rel, t)
            ext = corrected_ext(src)
            name = place_link(md_dir, src, ext)
            trajs.append(name)
            if ext:
                relabelled[name] = t
        # Only the path-bearing keys are rewritten to name links. Every OTHER
        # key on the entry is per-MD metadata the workflow uploads as that MD's
        # own record (framestep, temp, version, ...). Rebuilding the entry from
        # scratch dropped it, silently undoing the distribution of any field
        # whose replicas disagreed.
        entry = {k: v for k, v in m.items()
                 if k not in ("name", "mdir", "input_trajectory_filepaths",
                              "input_topology_filepath")}
        entry = {"name": md_id, "mdir": md_id,
                 "input_trajectory_filepaths": trajs, **entry}
        record = {"md_id": md_id, "source_dir": rel,
                  "trajectories": list(m.get("input_trajectory_filepaths") or []),
                  "links": trajs,
                  # By MD name; by directory for a doc built before replicas
                  # could be split out of one.
                  "ordering": ordering.get(md_id, ordering.get(rel, "unknown"))}
        if md_id in split:
            record["replica"] = int(split[md_id])
        if assumed.get(md_id):
            record["stage_assumed"] = list(assumed[md_id])
        if relabelled:
            # The file's content contradicts its extension; the link says what
            # it really is. The source filename is kept here, untouched.
            record["relabelled"] = relabelled
        if m.get("input_topology_filepath"):
            md_top = str(m["input_topology_filepath"])
            entry["input_topology_filepath"] = place_link(
                md_dir, os.path.join(source_root, md_top))
            record["topology"] = md_top
            record["topology_link"] = os.path.join(md_id,
                                                   entry["input_topology_filepath"])
        elif top:
            record["topology"] = top
            record["topology_link"] = runnable["input_topology_filepath"]
        out_mds.append(entry)
        provenance["mds"].append(record)

    runnable["mds"] = out_mds

    # Links from a previous build that no longer name an input -- a topology the
    # matcher has since ranked below a better one, a segment reclassified as
    # derived -- are removed. They are ours by construction (src.-prefixed
    # symlinks), and unlink removes a link without touching its target. Left in
    # place they would be harmless to mwf but misleading to a person.
    wanted = {project: {runnable.get("input_topology_filepath")}}
    for e in out_mds:
        wanted[os.path.join(project, e["name"])] = set(
            e["input_trajectory_filepaths"]) | {e.get("input_topology_filepath")}
    for directory, names in wanted.items():
        for n in os.listdir(directory):
            full = os.path.join(directory, n)
            if n.startswith(LINK_PREFIX) and os.path.islink(full) and n not in names:
                os.unlink(full)
    return runnable, provenance


def stage_inputs(mddb, dataset_id, export_path):
    """Copy the pristine export into the work directory for mwf to consume.

    Three cases, because mwf rewrites its inputs file in place:
      - no work copy yet: copy it.
      - the export has not changed since the last copy: leave the work copy
        alone, including whatever mwf legitimately wrote into it.
      - the export changed: the regenerated export is the truth, so it replaces
        the work copy -- after keeping the old one as inputs.yaml.previous.
    Returns one of 'created' | 'unchanged' | 'refreshed'.
    """
    project = os.path.join(authorised_root(mddb), dataset_id)
    work_inputs = os.path.join(project, "inputs.yaml")
    stamp = os.path.join(project, EXPORT_STAMP)
    new = _fingerprint(export_path)
    if not os.path.exists(work_inputs):
        shutil.copyfile(export_path, work_inputs)
        state = "created"
    else:
        old = None
        if os.path.exists(stamp):
            with open(stamp) as fh:
                old = fh.read().strip()
        if old == new:
            return "unchanged"
        # `.previous` deliberately does not end in `inputs.yaml`: mwf globs
        # `*inputs.yaml` in the project directory and stops on ambiguity.
        shutil.copyfile(work_inputs, work_inputs + ".previous")
        shutil.copyfile(export_path, work_inputs)
        state = "refreshed"
    from .atomic import write_text                           # noqa: PLC0415
    write_text(stamp, new + "\n")
    return state


def _only_our_links(directory):
    """True when a directory tree holds nothing but our links and empty dirs."""
    for dp, dns, fns in os.walk(directory):
        for fn in fns:
            full = os.path.join(dp, fn)
            ours = (fn.startswith(LINK_PREFIX) and os.path.islink(full)) \
                or (dp == directory and fn in ("inputs.yaml",
                                               "inputs.yaml.previous",
                                               EXPORT_STAMP))
            if not ours:
                return False
        for dn in dns:
            if os.path.islink(os.path.join(dp, dn)):
                return False
    return True


def prune(mddb, declared):
    """Remove working areas for datasets or MDs that are no longer exported.

    `declared` maps dataset_id -> set of md_ids. A stale area is removed only if
    it holds nothing but our own links and the staged inputs file. Anything else
    -- above all, workflow output -- is kept and reported, because processed
    results are expensive and a vanished MD id is a reason to look, not to
    delete. Links are removed with unlink, which never touches their targets.
    Returns (removed, kept).
    """
    root = authorised_root(mddb)
    removed, kept = [], []
    if not os.path.isdir(root) or os.path.islink(root):
        return removed, kept
    for did in sorted(os.listdir(root)):
        project = os.path.join(root, did)
        if not os.path.isdir(project) or os.path.islink(project):
            continue
        if did not in declared:
            if _only_our_links(project):
                shutil.rmtree(project)
                removed.append(did)
            else:
                kept.append(did)
            continue
        for md in sorted(os.listdir(project)):
            md_dir = os.path.join(project, md)
            if not os.path.isdir(md_dir) or os.path.islink(md_dir) \
                    or md in declared[did]:
                continue
            if _only_our_links(md_dir):
                shutil.rmtree(md_dir)
                removed.append(f"{did}/{md}")
            else:
                kept.append(f"{did}/{md}")
    return removed, kept


def verify(mddb):
    """Check every working area on disk. Returns a list of problems (empty = ok).

    Independent of how the areas were built, so it can be run before any
    workflow job as a last gate.
    """
    problems = []
    root = authorised_root(mddb)
    if not os.path.isdir(root):
        return [f"no working area at {root}"]
    if os.path.islink(root):
        return [f"{root} is a symlink"]
    for dp, dns, fns in os.walk(root):
        for name in list(fns) + list(dns):
            full = os.path.join(dp, name)
            if not os.path.islink(full):
                continue
            if not name.startswith(LINK_PREFIX):
                problems.append(f"unexpected symlink {full}")
                continue
            if name in RESERVED_NAMES:
                problems.append(f"link uses a workflow output name: {full}")
            target = os.path.realpath(full)
            if not os.path.exists(target):
                problems.append(f"dangling link {full} -> {os.readlink(full)}")
            elif within(target, root):
                problems.append(f"link points back into the work area: {full}")
        for dn in dns:
            full = os.path.join(dp, dn)
            if os.path.islink(full):
                problems.append(f"directory symlink inside work area: {full}")
    return problems
