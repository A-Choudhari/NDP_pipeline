#!/usr/bin/env python3
"""
copy.py -- a real copy, links kept as links, swapped into place atomically.

THREE PROPERTIES, EACH FOR A REASON THAT BIT SOMEBODY.

*Links stay links.* `MD_RESULTS/production_files/` holds symlinks pointing at
the raw trajectories with absolute paths. Dereferencing them turns a 2.4 GB
share into terabytes and, worse, produces a copy whose "trajectories" are
frozen duplicates rather than pointers at the run. They are recreated with
`os.symlink(os.readlink(src), dst)` -- target string unchanged, so they resolve
for everyone on this export. `rsync -L`, `cp -L` and `tar -h` follow symlinks;
the README says so, because a sync done that way is how the share becomes a
copy of the whole cluster.

*Nothing half-copied is ever named like a finished collection.* The tree is
built under `.collect/tmp/<PREFIX>.<pid>.partial/`, each file size-checked
and synced to disk, and only then renamed into place -- two renames on one filesystem, with the previous copy
displaced first so it survives until the new one is in position. The manifest
is written LAST, so "finished" means *manifest present and its fingerprint
matches the source*. A crash leaves either the old complete folder or a
`.partial` nobody claims, and the next run redoes the work.

*The sweep removes only what the manifest says we put there.* Files a person
added inside a collected folder are not in the manifest, so they are carried
forward into the new tree and reported as kept. The metadata pipeline once came
one build away from deleting 39 files by sweeping a tree it had not written.
"""

import hashlib
import os
import shutil
import time

from . import layout, manifest as M

CHUNK = 4 << 20


def copy_file(src, dst):
    """Copy, hashing the bytes as they pass, and return (digest, size).

    The hash is of the SOURCE stream, taken while the bytes are already moving,
    which is the only place it is free. The copy is checked by size and synced
    to disk before it can be renamed into place; re-reading it to compare
    hashes is what `--verify-hashes` does, at the cost of a second read.

    The full SHA-256 is kept. It was cut to 64 bits, which is a fingerprint,
    not an integrity check, for a collection meant to hold millions of files.
    """
    h = hashlib.sha256()
    n = 0
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        while True:
            b = fi.read(CHUNK)
            if not b:
                break
            h.update(b)
            fo.write(b)
            n += len(b)
        fo.flush()
        os.fsync(fo.fileno())
    shutil.copystat(src, dst, follow_symlinks=False)
    got = os.path.getsize(dst)
    if got != n:
        raise IOError(f"{dst}: wrote {n} bytes, file holds {got}")
    return h.hexdigest(), n


def build_tree(source, dest, dirs, files, links, carry=(), on_file=None):
    """Assemble the copy under `dest`. Returns {relpath: digest}."""
    os.makedirs(dest, exist_ok=True)
    for d in dirs:
        os.makedirs(os.path.join(dest, d), exist_ok=True)
    hashes = {}
    for rel in sorted(files):
        digest, _ = copy_file(os.path.join(source, rel), os.path.join(dest, rel))
        hashes[rel] = digest
        if on_file:
            on_file(rel)
    for rel in sorted(links):
        target = links[rel]
        full = os.path.join(dest, rel)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        if os.path.lexists(full):
            os.unlink(full)
        os.symlink(target, full)
    # Anything a person put inside the previous copy: brought across, never
    # eaten. It is reported as `kept (foreign)` and stays out of the manifest,
    # so we do not claim to own it on the next run either.
    for rel, src in carry:
        full = os.path.join(dest, rel)
        if os.path.lexists(full):
            continue
        os.makedirs(os.path.dirname(full), exist_ok=True)
        if os.path.islink(src):
            os.symlink(os.readlink(src), full)
        elif os.path.isfile(src):
            shutil.copy2(src, full)
    return hashes


def foreign(previous, old_manifest):
    """[(relpath, abspath)] inside a previous copy that we did not write."""
    if not previous or not os.path.isdir(previous) or not old_manifest:
        return []
    ours = set(old_manifest.get("files") or {}) | set(old_manifest.get("links") or {})
    out = []
    for dirpath, dirnames, filenames in os.walk(previous, followlinks=False):
        dirnames[:] = [d for d in dirnames if d not in M.SKIP_DIRS]
        for name in filenames + [d for d in list(dirnames)
                                 if os.path.islink(os.path.join(dirpath, d))]:
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, previous)
            if rel not in ours:
                out.append((rel, full))
    return sorted(out)


def swap(root, prefix, staged):
    """Put `staged` in place, keeping the previous copy until it is there."""
    final = layout.collected(root, prefix)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    old = None
    if os.path.isdir(final):
        old = layout.displaced(root, prefix, stamp)
        os.rename(final, old)
    try:
        os.rename(staged, final)
    except BaseException:
        if old and not os.path.isdir(final):
            os.rename(old, final)        # put it back; nothing was lost
        raise
    if old:
        shutil.rmtree(old, ignore_errors=True)
    return final


def verify(dest, files_manifest):
    """Re-hash what is in the share against the manifest. Returns the bad ones."""
    bad = []
    for rel, meta in (files_manifest or {}).items():
        want = meta[2] if len(meta) > 2 else ""
        full = os.path.join(dest, rel)
        if not want:
            continue
        try:
            h = hashlib.sha256()
            with open(full, "rb") as fh:
                while True:
                    b = fh.read(CHUNK)
                    if not b:
                        break
                    h.update(b)
            # Manifests written before the full digest was kept hold its first
            # 16 characters; compare at whatever length was recorded.
            if h.hexdigest()[:len(want)] != want:
                bad.append(rel)
        except OSError:
            bad.append(rel)
    return bad
