#!/usr/bin/env python3
"""
collect.py -- gather finished MD_RESULTS folders into one shared place.

    collect.py [SOURCE ...] [--all] [--dest DIR] [--dry-run]

Each source is a project's `MD_RESULTS/` directory. Each is copied, whole and
byte for byte, to `<dest>/<PREFIX>_MD_RESULTS/`, where the prefix comes from
`collection_prefix:` in that project's settings.yaml or, failing that, from its
path.

This is a SEPARATE kit. It does not build records and it does not analyse
anything; it reads finished output and puts it where the lab can see it. The
two pipelines that produce that output are untouched by it.

WHAT IT WILL NOT DO
  * follow a symlink out of a tree, or dereference one into the share
  * write anything outside the destination
  * copy when a source is unchanged
  * overwrite one project's collection with another's: two sources wanting one
    prefix stops the whole run before a byte moves
  * remove anything it did not put there, including a file a person added
    inside a collection, which is carried forward and reported
"""

import argparse
import contextlib
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import copy as C            # noqa: E402
from core import index as IDX         # noqa: E402
from core import layout as L          # noqa: E402
from core import manifest as M        # noqa: E402
from core import prefix as P          # noqa: E402
from core import safety as S          # noqa: E402

#: How deep `--all` looks below each root.
MAX_DEPTH = 7


def default_roots():
    """Where `--all` looks: MD_COLLECT_ROOTS (separated by `:`), else the
    `roots:` list in sources.yaml beside this script. Read without a YAML
    library, so the collector keeps needing nothing but the standard library."""
    env = os.environ.get("MD_COLLECT_ROOTS")
    if env:
        return [r for r in env.split(os.pathsep) if r]
    roots, inside = [], False
    try:
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "sources.yaml")) as fh:
            for line in fh:
                text = line.split("#", 1)[0].rstrip()
                if not text.strip():
                    continue
                if not line.startswith((" ", "-")):
                    inside = text.strip() == "roots:"
                    continue
                if inside and text.strip().startswith("- "):
                    roots.append(os.path.expanduser(text.strip()[2:].strip().strip("'\"")))
    except OSError:
        pass
    return roots


def discover(roots, depth=MAX_DEPTH):
    """Every MD_RESULTS under these roots, never following a symlink."""
    found = []
    for root in roots:
        root = os.path.abspath(root)
        base = root.rstrip("/").count("/")
        for dirpath, dirnames, _ in os.walk(root, followlinks=False):
            if dirpath.count("/") - base >= depth:
                dirnames[:] = []
                continue
            dirnames[:] = [d for d in sorted(dirnames)
                           if not d.startswith(".")
                           and d not in {"conda_pkgs", "site-packages", "envs",
                                         "env", "venv", "node_modules",
                                         "__pycache__", "af3_py312", "toppar"}
                           and not os.path.islink(os.path.join(dirpath, d))]
            if "MD_RESULTS" in dirnames:
                cand = os.path.join(dirpath, "MD_RESULTS")
                if os.path.exists(os.path.join(cand, "settings.yaml")):
                    found.append(cand)
                dirnames.remove("MD_RESULTS")
    return sorted(set(found))


def plan_one(root, source, prefix, basis, force=False, no_links=False):
    """What collecting this source would do. Reads lstat only."""
    dirs, files, links = M.walk(source)
    skipped = {}
    if no_links:
        skipped["production_files"] = f"{len(links)} symlink(s) omitted (--no-production-files)"
        links = {}
    old = M.read(L.manifest(root, prefix))
    dest = L.collected(root, prefix)
    fresh = M.fingerprint(dirs, files, links)
    unchanged = (not force and old
                 and old.get("source_fingerprint") == fresh
                 and os.path.isdir(dest))
    added, replaced, same, swept = M.diff(old, dirs, files, links)
    return {
        "source": source, "prefix": prefix, "basis": basis, "dest": dest,
        "dirs": dirs, "files": files, "links": links, "skipped": skipped,
        "old": old, "unchanged": bool(unchanged),
        "added": added, "replaced": replaced, "same": same, "swept": swept,
        "bytes": sum(files[r][0] for r in (added + replaced) if r in files),
        "total_bytes": sum(v[0] for v in files.values()),
    }


def do_one(root, p, dry=False):
    """Copy one planned source into the share. Returns a row for the index."""
    source, prefix = p["source"], p["prefix"]
    folder = os.path.dirname(source)
    if p["unchanged"]:
        return None
    if dry:
        return None
    start = time.time()
    L.ensure(L.tmp(root))
    staged = L.partial(root, prefix)
    if os.path.isdir(staged):
        shutil.rmtree(staged, ignore_errors=True)
    S.check_write_target(root, staged)
    carry = C.foreign(p["dest"], p["old"])
    hashes = C.build_tree(source, staged, p["dirs"], p["files"], p["links"],
                          carry=carry)
    C.swap(root, prefix, staged)
    doc = M.build(prefix, p["basis"], source, folder, p["dirs"], p["files"],
                  p["links"], hashes, elapsed=time.time() - start,
                  skipped=p["skipped"])
    doc["carried_foreign"] = [rel for rel, _ in carry]
    L.ensure(L.manifests(root))
    M.write(L.manifest(root, prefix), doc)
    return doc


def row_for(root, prefix, doc):
    dest = L.collected(root, prefix)
    records, total_ns, ids = IDX.catalogue_of(dest)
    unresolved = sum(1 for rel in (doc.get("links") or {})
                     if not os.path.exists(os.path.join(dest, rel)))
    return {"prefix": prefix, "prefix_basis": doc.get("prefix_basis", "?"),
            "source_folder": doc.get("source_folder", ""),
            "collected": doc.get("collected", ""),
            "records": records, "total_ns": total_ns, "identifiers": ids,
            "files": (doc.get("counts") or {}).get("files", 0),
            "links": (doc.get("counts") or {}).get("links", 0),
            "bytes": doc.get("bytes", 0), "unresolved_links": unresolved}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="*", help="MD_RESULTS directories")
    ap.add_argument("--all", action="store_true", help="find them under the roots")
    ap.add_argument("--root", action="append", default=[],
                    help="a root for --all (repeatable)")
    ap.add_argument("--dest", default=L.SHARE_DEFAULT)
    ap.add_argument("--dry-run", action="store_true",
                    help="say what would happen; write nothing at all")
    ap.add_argument("--force", action="store_true", help="re-copy unchanged sources")
    ap.add_argument("--verify-hashes", action="store_true",
                    help="re-read the share and check it against the manifests")
    ap.add_argument("--no-production-files", action="store_true",
                    help="omit the production_files/ symlinks")
    ap.add_argument("--jobs", type=int, default=1,
                    help="kept at 1 on purpose: this export is the bottleneck, "
                         "and two readers each get half of it")
    a = ap.parse_args(argv)
    dry = a.dry_run

    sources = [os.path.abspath(s) for s in a.sources]
    if a.all or not sources:
        roots = a.root or default_roots()
        if not roots and not sources:
            print("\n  No roots for --all. Add your folders under `roots:` in "
                  f"{os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sources.yaml')},\n"
                  "  or pass --root DIR, or name MD_RESULTS directories.\n")
            return 2
        sources += discover(roots)
    if not sources:
        print("\n  No MD_RESULTS directories found. Name one, or use --all.\n")
        return 2

    try:
        sources = [S.check_source(s) for s in sorted(set(sources))]
        dest = S.check_destination(a.dest, sources)
    except S.Refused as exc:
        print(f"\n  Refused: {exc}\n")
        return 2

    folders = [os.path.dirname(s) for s in sources]
    resolved, clashes = P.partition(folders)
    by_folder = {f: (p, b) for f, p, b in resolved}
    skipped = set()
    for exc in clashes:
        print(f"\n  Skipped: {exc}")
        skipped.update((exc.first, exc.second))
    skipped -= set(by_folder)
    sources = [s for s in sources if os.path.dirname(s) in by_folder]

    print(f"\n  destination : {dest}{'   (dry run)' if dry else ''}")
    print(f"  sources     : {len(sources)}"
          + (f"  ({len(skipped)} skipped: prefix clash)" if skipped else "") + "\n")

    # Planned UNDER the lock when anything will be written: a plan made before
    # it could be stale by the time another collector let go.
    held = contextlib.ExitStack()
    if not dry:
        L.ensure(dest)
        L.ensure(L.internal(dest))
        try:
            held.enter_context(S.Lock(L.lock(dest)))
        except S.Refused as exc:
            print(f"\n  Refused: {exc}\n")
            return 2
        # Under the lock no other collector is running, so anything left in
        # tmp/ is a crashed run's: a `.partial` (named by a pid that will never
        # come back) or a displaced `.old` copy.
        tmpdir = L.tmp(dest)
        if os.path.isdir(tmpdir):
            for name in os.listdir(tmpdir):
                if name.endswith((".partial", ".old")):
                    shutil.rmtree(os.path.join(tmpdir, name), ignore_errors=True)

    with held:
        return _collect(a, dest, sources, by_folder, dry, skipped)


def _collect(a, dest, sources, by_folder, dry, skipped):
    plans = []
    for s in sources:
        pfx, basis = by_folder[os.path.dirname(s)]
        plans.append(plan_one(dest, s, pfx, basis, a.force,
                              a.no_production_files))

    for p in plans:
        n_new = len(p["added"])
        n_rep = len(p["replaced"])
        n_del = len(p["swept"])
        if p["unchanged"]:
            state = "unchanged"
        elif p["old"] is None:
            state = f"new    {len(p['files'])} files, {len(p['links'])} links, " \
                    f"{p['total_bytes'] / 1e6:.1f} MB"
        else:
            state = f"update +{n_new} added, ~{n_rep} replaced, -{n_del} gone"
        print(f"    {p['prefix']:<28} {state}")
        if p["basis"] == "derived":
            print(f"    {'':<28} prefix derived from the path; declare "
                  f"`collection_prefix:` to name it")

    if dry:
        print("\n  NOTHING WAS WRITTEN (--dry-run).\n")
        return 1 if skipped else 0

    rows, failed = [], []
    try:
        for p in plans:
            try:
                doc = do_one(dest, p, dry)
            except Exception as exc:                       # noqa: BLE001
                failed.append((p["prefix"], exc))
                print(f"    {p['prefix']:<28} FAILED: {exc}")
                # Its half-built tree is ours and worthless: remove it now
                # rather than leave it for the next run to find.
                shutil.rmtree(L.partial(dest, p["prefix"]), ignore_errors=True)
                continue
            doc = doc or M.read(L.manifest(dest, p["prefix"]))
            if doc:
                rows.append(row_for(dest, p["prefix"], doc))
        # EVERY collection in the share, not only the ones this run
        # touched. Building the index from `rows` alone dropped every
        # untouched collection from INDEX.md and index.json, so collecting
        # one folder made the other nine disappear from the index while
        # their copies sat right there on disk.
        rows = []
        for pfx in sorted(L.collections_in(dest)):
            doc = M.read(L.manifest(dest, pfx))
            if doc:
                rows.append(row_for(dest, pfx, doc))
            else:
                rows.append({"prefix": pfx, "prefix_basis": "unknown",
                             "source_folder": "", "collected": "",
                             "records": 0, "total_ns": 0, "identifiers": [],
                             "files": 0, "links": 0, "bytes": 0,
                             "unresolved_links": 0, "orphan": True})
        changed = IDX.write_all(dest, rows, dry)
    except S.Refused as exc:
        print(f"\n  Refused: {exc}\n")
        return 2

    if a.verify_hashes:
        print()
        for r in rows:
            doc = M.read(L.manifest(dest, r["prefix"]))
            bad = C.verify(L.collected(dest, r["prefix"]), (doc or {}).get("files"))
            print(f"    {r['prefix']:<28} "
                  f"{'ok' if not bad else str(len(bad)) + ' file(s) DIFFER'}")

    total = sum(r["records"] for r in rows)
    print(f"\n  {len(rows)} collection(s) in {dest}")
    print(f"  {total} record(s) · index: "
          f"{', '.join(changed) if changed else 'unchanged'}")
    if failed:
        print(f"  {len(failed)} FAILED: {', '.join(p for p, _ in failed)}")
    if skipped:
        print(f"  {len(skipped)} SKIPPED (prefix clash): {', '.join(sorted(skipped))}")
    print()
    return 1 if (failed or skipped) else 0


if __name__ == "__main__":
    raise SystemExit(main())
