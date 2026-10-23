#!/usr/bin/env python3
"""
identity.py -- deterministic, collision-free identities for datasets and MDs.

WHY THIS EXISTS. The exporter used to name each MD from the indexed replica
label (`1R`, `2R`, `rep`, `pr`) and fall back to the run directory's basename.
Both are reused across sibling systems, so three different DNA aptamers under
`Top1/`, `Top4/` and `Top5/` all became `1R`. Eleven of thirteen exported
datasets carried duplicate MD names; `gp120-gp41-glyco-HD4` had 93 MDs sharing
19 names.

That is not a cosmetic problem. The installed MDDB workflow (v0.1.10) refuses
the whole project before any processing starts:

    core/inputs.py:358
        raise InputError('Duplicated values in MD inputs (see warnings above).'
            ' All MD names and directories must be unique.')

and it raises that error only AFTER rewriting inputs.yaml in place.

THE FOUR PROPERTIES AN ID MUST HAVE, and how each is obtained here:

  deterministic   Derived only from the run's own path. No timestamps, no
                  random ids, no counters, no dependence on the order
                  os.listdir happened to return.

  stable          Derived from the run's own path AND a dataset root taken from
                  DECLARED configuration, never from the set of runs that
                  matched. Computing a common ancestor of the matched runs would
                  be shorter and prettier and would renumber every existing MD
                  the moment a colleague adds a run higher up the tree. A
                  published `<dataset_id>/<md_id>` must not move.

  unique          The path suffix is unique by construction, because two runs
                  cannot occupy one directory. Only sanitisation can collide.

  readable        Once a dataset is scoped to one system, the suffix is short:
                  for `.../top100/Top1/5/MODEL35` the id is simply `REP1`.

SANITISATION MUST MATCH THE WORKFLOW'S. mwf turns an MD name into a directory
by replacing spaces with '_' and DELETING '. , ; : º /'
(core/inputs.py:293-300, constants.py:232). So `rep 1` and `rep_1` both become
`rep_1`, and `A.B` and `AB` both become `AB` -- collisions invisible in our YAML
that surface only as an upstream error. We therefore collide-check on the
SANITISED form, not the pretty one, and disambiguate with a short hash of the
full source path when two distinct runs would otherwise fold together.
"""

import hashlib
import os
import re

#: Deleted outright by mwf when deriving a directory from a name
#: (mddb_workflow/utils/constants.py:232). 'º' is there because MDtraj cannot
#: read a path containing it.
MWF_FORBIDDEN = (".", ",", ";", ":", "º", "/")

#: Anything outside this set becomes '_'. Deliberately narrower than mwf's own
#: rule: mwf only strips the six characters above, which would let '(' or '#'
#: through into a directory name and then into a shell command somewhere
#: downstream.
_SAFE = re.compile(r"[^A-Za-z0-9_-]")


def sanitize(component):
    """One path component, reduced to a portable directory-safe token.

    Applies mwf's own transformation first so that our collision check sees
    exactly what the workflow will see, then narrows further.
    """
    text = str(component).replace(" ", "_")
    for ch in MWF_FORBIDDEN:
        text = text.replace(ch, "")
    text = _SAFE.sub("_", text)
    text = re.sub(r"_{2,}", "_", text).strip("_-")
    return text


def short_hash(text, length=6):
    """A stable, deterministic suffix. blake2b, not hash(), because Python's
    built-in hash is randomised per process and would produce a different id on
    every run."""
    return hashlib.blake2b(str(text).encode("utf-8"),
                           digest_size=8).hexdigest()[:length]


def dataset_root(spec, systems=None):
    """The declared path prefix a dataset's MD ids are measured from.

    Taken from configuration, never from the runs that matched -- see the
    `stable` property above. Returns "" when no prefix can be established, which
    yields ids built from the full relative path: long, but never wrong.
    """
    explicit = spec.get("md_root")
    if explicit:
        return str(explicit).strip("/")

    paths = spec.get("paths")
    if paths:
        cleaned = [str(p).strip("/") for p in paths if str(p).strip("/")]
        if not cleaned:
            return ""
        if len(cleaned) == 1:
            return os.path.dirname(cleaned[0])
        try:
            # commonpath over DECLARED paths is stable: the declaration is the
            # configuration, so it only changes when a person edits it.
            return os.path.commonpath(cleaned)
        except ValueError:
            return ""

    label = spec.get("system")
    if label and systems:
        declared = systems.get(str(label))
        # A system may be declared as one path or a list of them; a list means
        # one campaign living in several places, so there is no single prefix.
        if isinstance(declared, str):
            return declared.strip("/")
    return ""


def md_id(rel_path, root=""):
    """The identity of one MD within its dataset.

    `rel_path` is the run directory relative to the data root; `root` is the
    dataset's declared prefix. Components between the two become the id.
    """
    rel = str(rel_path).strip("/")
    root = str(root or "").strip("/")
    if root and (rel == root or rel.startswith(root + "/")):
        suffix = rel[len(root):].strip("/")
    else:
        suffix = rel
    parts = [sanitize(p) for p in suffix.split("/") if p]
    parts = [p for p in parts if p]
    if not parts:
        # The run directory IS the dataset root. Name it for its own basename
        # rather than emitting an empty id.
        parts = [sanitize(os.path.basename(rel)) or "md"]
    return "_".join(parts)


def assign_md_ids(rel_paths, root=""):
    """{rel_path: md_id} for one dataset, guaranteed collision-free.

    Two distinct runs can only collide after sanitisation. When they do, BOTH
    get a hash suffix rather than one keeping the bare name: picking a winner
    would depend on iteration order, and the loser's id would change if the
    winner were later removed.
    """
    proposed = {}
    for rel in rel_paths:
        proposed.setdefault(md_id(rel, root), []).append(str(rel))

    assigned = {}
    for ident, rels in proposed.items():
        if len(rels) == 1:
            assigned[rels[0]] = ident
            continue
        for rel in rels:
            assigned[rel] = f"{ident}_{short_hash(rel)}"
    return assigned


class IdentityCollision(Exception):
    """Raised before anything is written, naming every conflicting source path."""


def check_unique(entries, what="MD"):
    """Fail loudly on duplicates, naming the real paths behind each clash.

    `entries` is an iterable of (identity, source_path). The workflow's own
    error reports only the repeated name; a person then has to work out which
    directories produced it. This reports both sides.
    """
    seen = {}
    for ident, source in entries:
        seen.setdefault(ident, []).append(source)
    clashes = {k: v for k, v in seen.items() if len(v) > 1}
    if not clashes:
        return
    lines = [f"{len(clashes)} duplicated {what} identit"
             + ("ies" if len(clashes) > 1 else "y") + ":"]
    for ident, sources in sorted(clashes.items()):
        lines.append(f"  {ident!r} is claimed by {len(sources)} entries:")
        lines.extend(f"    - {s}" for s in sorted(sources))
    lines.append("")
    lines.append("The MDDB workflow rejects the whole project for this "
                 "(core/inputs.py:358), so the export is stopped here instead.")
    raise IdentityCollision("\n".join(lines))
