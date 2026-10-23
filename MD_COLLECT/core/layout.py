#!/usr/bin/env python3
"""
layout.py -- every path this kit writes, in one place.

The same rule the metadata pipeline's own `layout.py` follows: exactly one
module knows where things go, so "what does this write, and where" is answered
by reading one short file rather than by auditing every `open(..., "w")`.

THE SHARE IS PRISTINE. `<PREFIX>_MD_RESULTS/` is a real copy of somebody's
MD_RESULTS and NOTHING of ours goes inside it -- no manifest, no marker, no
log. A folder taken out of the share and handed to someone is byte for byte
what the source was. The cost is that such a folder, alone, no longer says
where it came from, which is why `index.json` carries the provenance of all of
them and why the README says so in its first lines.

DRY RUN IS A PROPERTY OF THIS MODULE. Every path-creating function takes the
run's `dry` flag and returns the path without creating anything when it is set,
so `--dry-run` cannot create the destination, the lock, a temp directory or a
log by accident. The alternative -- checking `if dry` at each of thirty call
sites -- is the shape that eventually misses one.
"""

import os

#: Where collections are copied. Set MD_COLLECT_DEST, or pass --dest.
SHARE_DEFAULT = os.environ.get("MD_COLLECT_DEST") or os.path.expanduser("~/MD_NDP_DATA")

#: What a collected folder is called. The prefix is the only variable part.
SUFFIX = "_MD_RESULTS"

#: Ours, hidden, mirroring `<folder>/.md_pipeline/` in the other kits.
INTERNAL = ".collect"

#: Written at the top of the share, for a person who opens it cold.
README = "README.md"
INDEX_MD = "INDEX.md"
INDEX_JSON = "index.json"

#: Never swept, never collected, never counted as a collection.
RESERVED = {README, INDEX_MD, INDEX_JSON, INTERNAL}


def ensure(path):
    """The ONLY thing in this module that creates a directory.

    Every other function returns a path and touches nothing. `--dry-run`
    created MD_NDP_DATA on its first outing because a path helper made the
    directory it was asked to name -- so the guarantee "writes nothing at all"
    is kept by having almost nothing here able to write.
    """
    os.makedirs(path, exist_ok=True)
    return path


def share(root=None):
    return os.path.abspath(root or SHARE_DEFAULT)


def collected(root, prefix):
    """Where one source's copy lives. Never created here: the copy is built
    under a temporary name and renamed into place, so this path existing is
    what `finished` means."""
    return os.path.join(os.path.abspath(root), f"{prefix}{SUFFIX}")


def internal(root):
    return os.path.join(os.path.abspath(root), INTERNAL)


def manifests(root):
    return os.path.join(internal(root), "manifests")


def manifest(root, prefix):
    return os.path.join(manifests(root), f"{prefix}.json")


def tmp(root):
    return os.path.join(internal(root), "tmp")


def partial(root, prefix):
    return os.path.join(tmp(root), f"{prefix}.{os.getpid()}.partial")


def displaced(root, prefix, stamp):
    """Where the previous copy is moved to, between the two renames."""
    return os.path.join(tmp(root), f"{prefix}.{stamp}.old")


def logs(root):
    return os.path.join(internal(root), "logs")


def lock(root):
    return os.path.join(internal(root), ".collect.lock")


def prefix_of(name):
    """`IgG-1F4_MD_RESULTS` -> `IgG-1F4`, or None if it is not one of ours."""
    return name[:-len(SUFFIX)] if name.endswith(SUFFIX) else None


def collections_in(root):
    """Every `<PREFIX>_MD_RESULTS` in the share, by prefix."""
    out = {}
    try:
        names = sorted(os.listdir(root))
    except OSError:
        return out
    for n in names:
        p = prefix_of(n)
        if p and os.path.isdir(os.path.join(root, n)):
            out[p] = os.path.join(root, n)
    return out
