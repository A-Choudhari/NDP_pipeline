#!/usr/bin/env python3
"""
mddb_build.py -- the MDDB deposition files, from the glossary the NDP pipeline
already built.

    python3 mddb_build.py <folder> [--datasets id,id]

    .md_pipeline/mddb/exports/<dataset>/inputs.yaml      mwf input, pristine
    .md_pipeline/mddb/exports/<dataset>/provenance.yaml  where each input lives
    .md_pipeline/mddb/work/<dataset>/<md_id>/            link-only mdir for mwf

WHAT CHANGED. This used to scan the data itself and build the MD glossary on
the way past, so a metadata record could only be produced by also assembling
deposition files. The glossary is now the NDP pipeline's product
(`NDP_PIPELINE/glossary_build.py`), and this reads it: the same MDs, the same
segments in the same order, the same topologies -- one measurement, two uses.

Run the NDP pipeline first; this never writes settings.yaml or the glossary.
The links under work/ point at your data read-only; nothing here writes through
them, and mwf must only ever be run from that working area.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import audit, glossary, layout, mwf, settings, workdir   # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=os.getcwd())
    ap.add_argument("--datasets", help="comma-separated ids; default all")
    a = ap.parse_args(argv)

    base = layout.base(a.folder)
    if not settings.exists(base):
        print(f"\n  No {settings.path(base)}")
        print(f"  Run the NDP pipeline on {base} first.\n")
        return 2
    try:
        g = glossary.load(base)
    except FileNotFoundError:
        print(f"\n  No MD glossary under {base}")
        print("  Run the NDP pipeline first; this describes what it found.\n")
        return 2

    doc = settings.load(base)
    root = settings.config_of(base, doc)["root"]
    specs = {str(s.get("id") or s.get("name")): s
             for s in settings.datasets_of(base, doc)}
    if a.datasets:
        keep = {s.strip() for s in a.datasets.split(",")}
        specs = {k: v for k, v in specs.items() if k in keep}
        g = g[g.dataset_id.astype(str).isin(keep)]
    if not len(g):
        print("\n  The glossary holds no MDs for that selection.\n")
        return 1

    mddb = layout.mddb_dir(base)
    print(f"\n  folder   : {base}")
    print(f"  glossary : {len(g)} MD(s) in {g.dataset_id.nunique()} dataset(s)")
    print(f"  writing  : {os.path.relpath(mddb, base)}/\n")

    written = 0
    for did, rows in g.groupby("dataset_id", sort=True):
        spec = specs.get(str(did)) or {"id": did}
        blocked = rows[rows.state == "blocked"]
        usable = rows[rows.state != "blocked"]
        if len(blocked):
            print(f"        {did}: {len(blocked)} blocked MD(s) left out "
                  f"(mwf would reject them): "
                  f"{', '.join(blocked.md_id.astype(str)[:3])}")
        if not len(usable):
            print(f"        {did}: nothing depositable")
            continue
        doc_out, cover = mwf.build_dataset(spec, usable, root)
        try:
            runnable, prov = workdir.materialise(mddb, str(did), doc_out, root)
        except workdir.UnsafeWorkdir as exc:
            print(f"        {did}: CANNOT PREPARE A SAFE WORKING DIRECTORY\n"
                  f"        {exc}\n\n  Nothing further was written.\n")
            return 1
        out = os.path.join(layout.exports(base), str(did), "inputs.yaml")
        mwf.dump(runnable, out, dataset_id=str(did))
        mwf.write_provenance(
            os.path.join(layout.exports(base), str(did), "provenance.yaml"),
            str(did), prov)
        staged = workdir.stage_inputs(mddb, str(did), out)
        cover += audit.declaration(spec, len(doc_out["mds"]))
        for field, status, note in cover:
            if status == "WARNING":
                print(f"        {did}: WARNING `{field}` {note}")
        filled = sum(1 for _, st, _ in cover if st == "filled")
        absent = sum(1 for _, st, _ in cover if st == "absent")
        print(f"        {did}: {len(doc_out['mds'])} MD(s), {filled} fields "
              f"filled, {absent} absent -> "
              f"{os.path.relpath(out, base)} (work copy {staged})")
        written += 1

    if not a.datasets:
        # A dataset no longer declared must not sit in exports/ looking
        # depositable; these are our own regenerable output.
        import shutil
        declared = set(g.dataset_id.astype(str))
        exports = layout.exports(base)
        if os.path.isdir(exports):
            for name in sorted(os.listdir(exports)):
                if name not in declared and os.path.isdir(os.path.join(exports, name)):
                    shutil.rmtree(os.path.join(exports, name))
                    print(f"        removed stale export '{name}'")
        removed, kept = workdir.prune(mddb, {str(d): set(rows.md_id.astype(str))
                                             for d, rows in g.groupby("dataset_id")})
        for name in removed:
            print(f"        removed stale working area '{name}' (links only)")
        for name in kept:
            print(f"        KEPT stale working area '{name}': it holds more "
                  f"than links (workflow output?). Inspect before deleting.")

    problems = workdir.verify(mddb)
    if problems:
        print(f"\n  WORKING AREA CHECK FAILED ({len(problems)}):")
        for p in problems[:20]:
            print(f"    {p}")
        return 1

    print(f"\n  Done. {written} dataset(s) under "
          f"{os.path.relpath(layout.exports(base), base)}/")
    print(f"  Run mwf ONLY from the working area, never on your data:")
    print(f"      cd {os.path.join(mddb, 'work', '<dataset>')} && "
          f"mwf run -dir . -inp inputs.yaml\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
