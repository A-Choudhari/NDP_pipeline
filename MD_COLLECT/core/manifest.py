#!/usr/bin/env python3
"""
manifest.py -- what this collector put in the share, and nothing else.

The manifest is the single authority for the sweep. A path listed in the
previous manifest and absent from the new one is removed; anything not listed
was not put there by us and is never touched. That rule exists because the
metadata pipeline once came one build away from deleting 39 files by sweeping a
tree it had not written -- every `frames.json` among them. *A sweep that
removes what it did not write is not a sweep.*

MTIME IN NANOSECONDS, not the whole seconds the other kits' manifests carry.
On NFS a file rewritten within the same second, to the same length, is exactly
the change a second-granularity comparison cannot see. `--verify-hashes`
remains the only thing that catches a same-size same-mtime rewrite, and the
README says so plainly rather than implying the skip is proof.

THE FINGERPRINT IS COMPUTED FROM `lstat` ALONE. No file is opened to decide
whether a source changed: hashing 2.4 GB to answer "did anything move" costs
half a minute of a shared export that forty people are reading, every run, to
tell you nothing new. Bytes are hashed while they are being copied, on both
streams, which is where the hash is free.
"""

import hashlib
import json
import os
import tempfile
import time

VERSION = "MD_COLLECT/1"

#: Directory entries we never carry: the other kits' machinery, and ours.
SKIP_DIRS = {".md_pipeline", ".collect", "__pycache__"}


def _atomic(path, text):
    d = os.path.dirname(path)
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".manifest.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def walk(root):
    """(dirs, files, links) below `root`, relative and sorted.

    A symlinked DIRECTORY is recorded as a link and never descended: following
    one would copy whatever it points at, which in `production_files/` is a
    32 GB trajectory, and could lead straight out of the tree.
    """
    dirs, files, links = [], {}, {}
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        keep = []
        for d in sorted(dirnames):
            full = os.path.join(dirpath, d)
            if d in SKIP_DIRS:
                continue
            if os.path.islink(full):
                links[os.path.relpath(full, root)] = os.readlink(full)
                continue
            keep.append(d)
        dirnames[:] = keep
        rel_dir = os.path.relpath(dirpath, root)
        if rel_dir != ".":
            dirs.append(rel_dir)
        for name in sorted(filenames):
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, root)
            if os.path.islink(full):
                links[rel] = os.readlink(full)
                continue
            try:
                st = os.lstat(full)
            except OSError:
                continue
            files[rel] = [st.st_size, st.st_mtime_ns]
    return sorted(set(dirs)), files, links


def fingerprint(dirs, files, links):
    """One digest over the shape of a tree, from lstat data only."""
    h = hashlib.sha256()
    for d in sorted(dirs):
        h.update(f"d\0{d}\0\n".encode())
    for rel in sorted(files):
        size, mtime = files[rel]
        h.update(f"f\0{rel}\0{size}\0{mtime}\n".encode())
    for rel in sorted(links):
        h.update(f"l\0{rel}\0{links[rel]}\n".encode())
    return h.hexdigest()[:16]


def read(path):
    try:
        with open(path) as fh:
            doc = json.load(fh)
    except (OSError, ValueError):
        return None
    return doc if isinstance(doc, dict) and doc.get("collector") else None


def write(path, doc):
    _atomic(path, json.dumps(doc, indent=1, sort_keys=True) + "\n")
    return path


def build(prefix, basis, source, folder, dirs, files, links, hashes,
          elapsed=0.0, skipped=None):
    return {
        "collector": VERSION,
        "prefix": prefix,
        "prefix_basis": basis,
        "source": source,
        "source_folder": folder,
        "destination": f"{prefix}_MD_RESULTS",
        "collected": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "collected_by": f"{os.environ.get('USER', 'unknown')}@"
                        f"{os.uname().nodename.split('.')[0]}",
        "elapsed_s": round(elapsed, 1),
        "bytes": sum(v[0] for v in files.values()),
        "counts": {"dirs": len(dirs), "files": len(files), "links": len(links)},
        "source_fingerprint": fingerprint(dirs, files, links),
        "dirs": sorted(dirs),
        "files": {k: [v[0], v[1], hashes.get(k, "")] for k, v in sorted(files.items())},
        "links": {k: links[k] for k in sorted(links)},
        "skipped": skipped or {},
    }


def diff(old, dirs, files, links):
    """(added, replaced, unchanged, swept) between a manifest and a fresh walk."""
    if not old:
        return sorted(files) + sorted(links), [], [], []
    was_f = {k: (v[0], v[1]) for k, v in (old.get("files") or {}).items()}
    was_l = dict(old.get("links") or {})
    added, replaced, unchanged = [], [], []
    for rel, (size, mtime) in files.items():
        if rel not in was_f:
            added.append(rel)
        elif was_f[rel] != (size, mtime):
            replaced.append(rel)
        else:
            unchanged.append(rel)
    for rel, target in links.items():
        if rel not in was_l:
            added.append(rel)
        elif was_l[rel] != target:
            replaced.append(rel)
        else:
            unchanged.append(rel)
    now = set(files) | set(links)
    swept = sorted((set(was_f) | set(was_l)) - now)
    return sorted(added), sorted(replaced), sorted(unchanged), swept
