#!/usr/bin/env python3
"""
prefix.py -- what a collected folder is called, and who decides.

Two sources, and which one was used is always reported:

  declared  `collection_prefix:` at the top of the source's own
            `MD_RESULTS/settings.yaml`. The owner's word.
  derived   built from the path, because nobody said. Deterministic, and it
            prints the line telling you to declare one instead.

WHY A NEW KEY IS SAFE TO PUT IN THAT FILE. The metadata pipeline never
rewrites settings.yaml wholesale: its editors change the text under a NAMED
section and leave unknown top-level keys alone, and the file is only created
when it is absent. A key at column zero, above `folder:`, is therefore
untouched by any number of rebuilds. It must go ABOVE `folder:`, because the
editors append at the end of the section they are given and stop at the next
column-zero key.

COLLISIONS ARE REFUSED, NOT RESOLVED. Two sources reaching one prefix is a
question about what these datasets ARE, and the wrong answer overwrites
somebody's collection with somebody else's. So every prefix for every source
is resolved before a single byte is copied, and a clash stops the whole run --
including the folders that did not clash, because half a shared folder is
harder to reason about than an unchanged one.
"""

import os
import re

#: Directory names that say where a run sits rather than what it is. Taken from
#: the same instinct as the metadata pipeline's `stage_dir` notion: `MDs`,
#: `2-REPLICA` and `runs` are carriers, not names.
CARRIER = {"md", "mds", "md_results", "runs", "run", "replica", "replicas",
           "data", "prod", "production", "analysis"}

_CARRIER_RE = re.compile(r"^\d+[-_]?replica$|^replica\d*$|^\d+r$", re.I)


def _is_carrier(name):
    low = name.lower()
    return low in CARRIER or low.isdigit() or bool(_CARRIER_RE.match(low))


def sanitize(text):
    """Letters, digits, dot, dash and underscore. Nothing that needs quoting
    in a shell, a URL or a file name."""
    out = re.sub(r"[^A-Za-z0-9._-]+", "-", str(text)).strip("-._")
    return re.sub(r"-{2,}", "-", out)


def declared(folder):
    """`collection_prefix:` from the source's settings.yaml, or None.

    Read with a two-line parser rather than a YAML load: this runs before
    anything is copied, on folders whose settings file may be mid-edit, and a
    parse error in an unrelated part of the file must not stop a collection.
    """
    p = os.path.join(folder, "MD_RESULTS", "settings.yaml")
    try:
        with open(p) as fh:
            for line in fh:
                m = re.match(r"^collection_prefix:\s*(.+?)\s*$", line)
                if m:
                    return sanitize(m.group(1).strip("\"'")) or None
                if re.match(r"^[A-Za-z_]+:", line) and not line.startswith(
                        "collection_prefix"):
                    continue
    except OSError:
        return None
    return None


def derive(folder, parts=2):
    """A prefix from the path, when nobody declared one.

    Takes the last components that NAME something, skipping the carriers, and
    writes them in PATH order so they read the way the tree does:
    `ABZYME/8-IgG/1-IF4` gives `8-IgG-1-IF4`, and `ABZYME/1-FULL-GP120/MDs`
    gives `1-FULL-GP120`, because `MDs` says where and not what.

    It is deliberately not clever. A derived prefix is reported as derived
    every single run, with the line that says to declare one, because only the
    owner knows whether `CDR` means this folder's CDR or another's.
    """
    keep = []
    for name in reversed(os.path.abspath(folder).strip("/").split("/")):
        if _is_carrier(name):
            continue
        keep.append(sanitize(name))
        if len(keep) >= parts:
            break
    keep.reverse()                       # path order: outermost first
    return "-".join(keep) or sanitize(os.path.basename(folder)) or "unnamed"


def resolve(folders):
    """[(folder, prefix, basis)], or raise on a collision.

    Declared prefixes are resolved first, so a derived one can never take a
    name somebody asked for.
    """
    out, owner = [], {}
    pending = []
    for f in folders:
        d = declared(f)
        (out if d else pending).append((f, d, "declared") if d else (f, None, None))
    for f, p, basis in list(out):
        if p in owner:
            raise Collision(p, owner[p], f, "both declared")
        owner[p] = f
    for f, _, _ in pending:
        p = derive(f)
        if p in owner:
            raise Collision(p, owner[p], f,
                            "one derived from its path; declare "
                            "`collection_prefix:` in one of them")
        owner[p] = f
        out.append((f, p, "derived"))
    return sorted(out, key=lambda t: t[1])


def partition(folders):
    """([(folder, prefix, basis)], [Collision]) -- collect what is unambiguous.

    `resolve` stops at the first clash, which at a few thousand folders means
    one pair of colliding names blocks every collection. Here the clashing
    folders are set aside and everything else goes ahead. Both folders of a
    clash are set aside, never one: which of two folders claiming a name should
    get it is a person's decision. A derived prefix that would take a declared
    one is the exception -- the declared one was chosen, so it keeps the name.
    """
    claims = {}
    for f in folders:
        d = declared(f)
        if d:
            claims.setdefault(d, []).append((f, "declared"))
        else:
            claims.setdefault(derive(f), []).append((f, "derived"))
    ok, clashes = [], []
    for p, owners in claims.items():
        decl = [o for o in owners if o[1] == "declared"]
        if len(owners) == 1:
            ok.append((owners[0][0], p, owners[0][1]))
            continue
        if len(decl) == 1:
            ok.append((decl[0][0], p, "declared"))
            for f, _ in owners:
                if f != decl[0][0]:
                    clashes.append(Collision(p, decl[0][0], f,
                                             "one derived from its path; declare "
                                             "`collection_prefix:` in one of them"))
            continue
        first = owners[0][0]
        for f, _ in owners[1:]:
            clashes.append(Collision(p, first, f, "both declared" if decl
                                     else "both derived from their paths"))
    return sorted(ok, key=lambda t: t[1]), clashes


class Collision(Exception):
    def __init__(self, prefix, first, second, why):
        self.prefix, self.first, self.second, self.why = prefix, first, second, why
        super().__init__(
            f"two sources want the prefix {prefix!r} ({why}):\n"
            f"      {first}\n"
            f"      {second}\n"
            f"  Nothing was copied. Declare `collection_prefix:` at the top of "
            f"one of their MD_RESULTS/settings.yaml, above `folder:`.")
