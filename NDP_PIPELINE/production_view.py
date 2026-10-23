#!/usr/bin/env python3
"""
production_view.py -- one folder holding exactly what each record publishes.

    python3 production_view.py <folder>            # build it
    python3 production_view.py <folder> --remove   # take it away again

    <folder>/MD_RESULTS/production_files/
      README.md
      <record>/<replica>/topology.psf    -> the topology the record names
      <record>/<replica>/01_pr.dcd       -> the production segments, IN ORDER
      <record>/<replica>/PR.conf         -> the config those segments ran under
      <record>/<replica>/PR.log          -> its log, when one sits beside it

WHY LINKS AND NOT A COPY. The ask behind this (2026-09-22) was a folder holding
"only the production MD files with the topology and the .conf", to screen what
goes out. Moving the files would mean moving ~450 GB, breaking the restart
chain on the cluster (`run.sh` resumes from `pr.part*.restart.xsc` NEXT TO the
trajectories), and losing the one thing a record needs most -- where the data
actually lives. Every entry here is a symlink, so the view costs nothing, can
be rebuilt at any time, and `--remove` leaves no trace.

Nothing is copied, moved or deleted outside `production_files/`, and every link
is read-only by intention: never write through one.

CAUTION. `rsync -L`, `cp -L` and `tar -h` FOLLOW symlinks. Sync this directory
without those flags, or you will copy the whole dataset by accident.
"""

import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import glossary, layout, ndp                        # noqa: E402

VIEW = "production_files"


def _records(base):
    """Every record this folder publishes: (file stem, the record document)."""
    import json
    out = []
    for d in (layout.results(base), layout.records_dir(base, 99)):
        if not os.path.isdir(d):
            continue
        for name in sorted(os.listdir(d)):
            if not name.endswith(".json") or name == "catalog.json":
                continue
            path = os.path.join(d, name)
            try:
                with open(path) as fh:
                    doc = json.load(fh)
            except Exception:                                 # noqa: BLE001
                continue
            if isinstance(doc, dict) and doc.get("lab_glossary"):
                out.append((name[:-5], doc))
    return out


def _link(target, link_path, report):
    """One symlink, replacing only a symlink we made before."""
    if os.path.islink(link_path):
        if os.path.realpath(link_path) == os.path.realpath(target):
            return False
        os.unlink(link_path)
    elif os.path.exists(link_path):
        report.append(f"kept a real file, no link made: {link_path}")
        return False
    os.symlink(target, link_path)
    return True


def _config_and_log(row):
    """The production config, and the log beside it if its stem matches."""
    cfg = str(row.get("config_path") or "").strip()
    if not cfg or cfg.lower() == "nan":
        return []
    root = row.get("source_root")
    path = cfg if os.path.isabs(cfg) else os.path.join(root, cfg)
    if not os.path.exists(path):
        return []
    out = [path]
    stem = os.path.splitext(path)[0]
    for ext in (".log", ".out"):
        for cand in (stem + ext, stem.lower() + ext, stem.upper() + ext):
            if os.path.exists(cand) and cand not in out:
                out.append(cand)
                break
    return out


def build(base):
    g = glossary.load(base)
    rows = {str(r.md_key): r for _, r in g.iterrows()}
    root = os.path.join(layout.results(base), VIEW)
    os.makedirs(root, exist_ok=True)
    report, wanted, links, gb = [], {os.path.join(root, "README.md")}, 0, 0.0

    for stem, doc in _records(base):
        for key in doc["lab_glossary"]["md_keys"]:
            row = rows.get(str(key))
            if row is None:
                report.append(f"{stem}: {key} is not in the glossary")
                continue
            here = os.path.join(root, stem, str(row.md_id))
            os.makedirs(here, exist_ok=True)
            wanted.add(here)
            top, segs = glossary.input_paths(row)
            pairs = []
            if top and os.path.exists(top):
                pairs.append((top, "topology" + os.path.splitext(top)[1]))
            width = max(2, len(str(len(segs))))
            for i, s in enumerate(segs, 1):
                if os.path.exists(s):
                    pairs.append((s, f"{i:0{width}d}_{os.path.basename(s)}"))
                    gb += os.path.getsize(s) / 1e9
                else:
                    report.append(f"{stem}/{row.md_id}: missing {s}")
            for p in _config_and_log(row):
                pairs.append((p, os.path.basename(p)))
            for target, name in pairs:
                dest = os.path.join(here, name)
                wanted.add(dest)
                links += 1
                _link(target, dest, report)

    # Anything left from an earlier run is ours and regenerable -- but only
    # ever a link or an empty directory of links.
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        for name in filenames:
            p = os.path.join(dirpath, name)
            if p not in wanted and os.path.islink(p):
                os.unlink(p)
        if dirpath != root and dirpath not in wanted and not os.listdir(dirpath):
            os.rmdir(dirpath)

    with open(os.path.join(root, "README.md"), "w") as fh:
        fh.write(_readme(base, links, gb))
    return links, gb, report


def _readme(base, links, gb):
    return f"""# Production files, and nothing else

Every entry below is a **symlink** into this folder's simulations: the topology
each record names, its production segments in time order, and the config (and
log) they ran under. Minimisation, heating and NPT equilibration are not here,
and neither is anything else the records do not publish.

- {links} link(s), {gb:,.1f} GB of trajectory referenced, 0 bytes copied.
- Rebuild: `python3 {os.path.join(HERE, 'production_view.py')} {base}`
- Remove : the same command with `--remove` (it deletes links, never data).

`rsync -L`, `cp -L` and `tar -h` follow symlinks. Sync this directory WITHOUT
those flags, or you will copy every trajectory by accident.

Read through these links; never write through them. The files they point at are
the live simulations, and on the cluster their restart chain still resumes from the
`pr.part*.restart.xsc` that sit beside them.
"""


def remove(base):
    root = os.path.join(layout.results(base), VIEW)
    if not os.path.isdir(root):
        return 0
    real = [os.path.join(dp, f) for dp, _, fs in os.walk(root) for f in fs
            if not os.path.islink(os.path.join(dp, f)) and f != "README.md"]
    if real:
        print("\n  Not removing: these are real files, not links:")
        for p in real[:10]:
            print(f"    {p}")
        return 1
    shutil.rmtree(root)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=os.getcwd())
    ap.add_argument("--remove", action="store_true",
                    help="delete the view (links only) and stop")
    a = ap.parse_args(argv)
    base = layout.base(a.folder)

    if a.remove:
        code = remove(base)
        if code == 0:
            print(f"\n  Removed {os.path.join(layout.results(base), VIEW)}/ "
                  f"-- links only; no simulation was touched.\n")
        return code

    try:
        links, gb, report = build(base)
    except FileNotFoundError:
        print(f"\n  No MD glossary under {base}")
        print(f"  Run {os.path.join(HERE, 'ndp_pipeline.sh')} {base} first.\n")
        return 2

    where = os.path.join(layout.results(base), VIEW)
    print(f"\n  {links} link(s) under {os.path.relpath(where, base)}/ "
          f"-- {gb:,.1f} GB referenced, nothing copied")
    for line in report[:15]:
        print(f"    {line}")
    print(f"  Read through them; never write through them. rsync/cp/tar must "
          f"NOT be given -L or -h here.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
