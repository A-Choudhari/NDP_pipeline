#!/usr/bin/env python3
"""
glossary_build.py -- one stable, named entry per MD, from the index.

    python3 glossary_build.py <folder> [--datasets id,id]

Reads the trajectory index (index_build.py) and settings.yaml, scans each run's
MD config and log for what no trajectory header carries, and writes:

    MD_RESULTS/glossary.md           every MD, for reading
    MD_RESULTS/glossary.parquet      the same, for scripts
    .md_pipeline/index/mdconfigs.parquet, engines.parquet   what each said, and where
    .md_pipeline/coverage.md         what could not be established, and why

This is the step that used to sit inside the MDDB exporter, so a metadata
record could only be built by also assembling deposition files. It needs
neither: the MDs come from the index, the files' own headers and the configs.

Writes only inside <folder>/MD_RESULTS/ and <folder>/.md_pipeline/.
"""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import config as cfglib                            # noqa: E402
from core import (atomic, audit, engine, glossary, identity, layout,   # noqa: E402
                  mdconfig, mds, ndp, settings, stages)


def select_runs(spec, runs, production_only=True):
    """The indexed runs one declared dataset selects.

    With `production_only` (the default) this is what can be described: runs the
    index credits with at least one production file. The gap between the two is
    the honest answer to "where did the others go", and coverage.md states it.
    """
    prod = runs
    if production_only and "n_production" in runs.columns:
        prod = runs[runs.n_production > 0]
    if spec.get("system"):
        sel = prod[prod.system.astype(str) == str(spec["system"])]
        if spec.get("condition"):
            sel = sel[sel.condition.astype(str) == str(spec["condition"])]
        return sel
    if spec.get("paths"):
        wanted = {str(p).strip("/") for p in spec["paths"]}
        return prod[prod.rel_path.astype(str).isin(wanted)]
    return prod.iloc[0:0]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=os.getcwd())
    ap.add_argument("--datasets", help="comma-separated ids; default all")
    a = ap.parse_args(argv)

    base = layout.base(a.folder)
    try:   # one build per folder at a time: concurrent runs clobber and sweep
        atomic.hold(layout.lock(base, "build"), what=base)
    except atomic.Locked as exc:
        print(f"\n  Refused: {exc}\n")
        return 2
    runs_p = os.path.join(layout.registry(base), "runs.parquet")
    if not os.path.exists(runs_p):
        print(f"\n  No index at {runs_p}")
        print(f"  Run:  python3 {os.path.join(HERE, 'index_build.py')} {base}\n")
        return 2
    if not settings.exists(base):
        print(f"\n  No {settings.path(base)}\n")
        return 2
    if not cfglib.check_dependencies():
        return 2

    import pandas as pd
    doc = settings.load(base)
    cfg = settings.config_of(base, doc)
    root = cfg["root"]
    classify = stages.classifier(cfg)
    try:
        replica_rules = cfglib.compile_replica_rules(cfg)
    except ValueError as exc:
        print(f"\n  {exc}\n")
        return 1
    runs = pd.read_parquet(runs_p)
    specs = settings.datasets_of(base, doc)
    if a.datasets:
        keep = {s.strip() for s in a.datasets.split(",")}
        specs = [s for s in specs if str(s.get("id")) in keep]
    if not specs:
        print(f"\n  No datasets declared in {settings.path(base)}\n")
        return 1

    # Dataset ids name the records and prefix every md_key, so two datasets
    # sharing one would make `<dataset>/<md_id>` ambiguous.
    try:
        identity.check_unique(
            [(str(s.get("id") or s.get("name") or "dataset"),
              f"settings.yaml datasets entry {i}") for i, s in enumerate(specs)],
            "dataset id")
    except identity.IdentityCollision as exc:
        print(f"\n  {exc}\n")
        return 1

    print(f"\n  folder   : {base}")
    print(f"  index    : {len(runs)} run director"
          f"{'ies' if len(runs) != 1 else 'y'}")
    print(f"  datasets : {len(specs)}\n")

    needed = sorted({str(p) for s in specs
                     for p in select_runs(s, runs).rel_path.astype(str)})
    if not needed:
        print("  No production runs in any declared dataset.\n")
        return 1

    print(f"  [1/3] reading MD configs for {len(needed)} run dir(s) ...")
    t0 = time.time()
    hdr_dt = {k: v for k, v in zip(runs.rel_path.astype(str),
                                   runs.dt_ps_production) if v == v}
    cfgs = mdconfig.scan(root, needed, verbose=True, header_dt_ps=hdr_dt)
    found = int(cfgs.config_path.notna().sum())
    inherited = int((cfgs.config_scope == "inherited").sum())
    print(f"        {found}/{len(cfgs)} have a usable config "
          f"({found - inherited} beside the trajectory, {inherited} inherited) "
          f"({time.time() - t0:.0f}s)")
    reg = layout.registry(base)
    os.makedirs(reg, exist_ok=True)
    atomic.to_parquet(cfgs, os.path.join(reg, "mdconfigs.parquet"),
                      compression="zstd", index=False)

    print("\n  [2/3] reading program and version from run logs ...")
    t0 = time.time()
    cfg_engine = dict(zip(cfgs.rel_path.astype(str), cfgs.config_engine)) \
        if "config_engine" in cfgs.columns else {}
    engs = engine.scan(root, needed, cfg_engine, verbose=True)
    by_log = int((engs.source == "log").sum())
    print(f"        {by_log} from a log banner (with version), "
          f"{len(engs) - by_log} without ({time.time() - t0:.0f}s)")
    atomic.to_parquet(engs, os.path.join(reg, "engines.parquet"),
                      compression="zstd", index=False)
    eng_by_run = {str(r.rel_path): {"program": r.program, "version": r.version}
                  for _, r in engs.iterrows()}

    print("\n  [3/3] building the glossary ...")
    rows, cover_all = [], []
    for spec in specs:
        did = str(spec.get("id") or spec.get("name") or "dataset")
        sel = select_runs(spec, runs)
        matched = select_runs(spec, runs, production_only=False)
        if not len(sel):
            print(f"        {did}: NO RUNS MATCHED -- check its selector")
            cover_all.append((did, [("mds", "absent", "selector matched no runs")]))
            continue
        md_root = identity.dataset_root(spec, cfg.get("systems") or {})
        ids = identity.assign_md_ids([str(p) for p in sel.rel_path.astype(str)],
                                     md_root)
        recs, notes = mds.md_records(root, sel, classify, cfglib.TRAJ_EXT, ids,
                                     basis=getattr(classify, "basis", None),
                                     replicas=replica_rules)
        # Which topology CONTENT each MD ran from: same content, same system.
        groups = mds.topology_groups(
            [os.path.join(root, r["topology"]) for r in recs if r["topology"]])
        by_rel = {r["topology"]: groups.get(os.path.join(root, r["topology"]))
                  for r in recs if r["topology"]}
        n_groups = len(set(by_rel.values()))
        sub = cfgs[cfgs.rel_path.isin(sel.rel_path.astype(str))]
        # The frame spacing a person declared for this dataset's systems, for
        # the MDs whose headers carry none. The glossary is where the analysis
        # pipeline reads the spacing from, so a declaration that stopped at the
        # record left the two disagreeing about the same MD.
        declared = ndp.declared_intervals(settings.ndp_of(base, doc), did)
        rows += glossary.build_rows(
            did, spec, recs, sel, sub, root, engines=eng_by_run,
            lengths=spec.get("production_lengths_ns")
            or settings.production_lengths(base, doc)
            or glossary.PRODUCTION_LENGTHS_NS,
            groups=by_rel, heterogeneous=n_groups > 1, declared_ps=declared)
        cover = list(notes) + audit.declaration(spec, len(recs), matched, sel)
        if n_groups > 1:
            cover.append(("mds", "note",
                          f"{n_groups} distinct topologies by content: this "
                          f"dataset holds more than one molecular system, and "
                          f"gets one record per system."))
        cover_all.append((did, cover))
        print(f"        {did}: {len(recs)} MD(s), {n_groups} system(s)")

    try:
        identity.check_unique(
            [(str(r.get("md_key")), str(r.get("source_dir") or r.get("md_key")))
             for r in rows], "glossary md_key")
    except identity.IdentityCollision as exc:
        print(f"\n  CANNOT WRITE THE GLOSSARY\n\n  {exc}\n")
        return 1

    g = glossary.write(layout.mddb_dir(base), rows,
                       page=layout.glossary_md(base))
    # The machine-readable twin, beside the page a person reads: an analysis
    # should not have to know where the pipeline keeps its registry.
    atomic.to_parquet(g, os.path.join(layout.results(base), "glossary.parquet"),
                      compression="zstd", index=False)
    _coverage(base, cover_all, found, len(cfgs))
    # Which settings this glossary was built from, so ndp_build.py can tell a
    # person that an edit to them has not reached the records yet.
    atomic.write_json(os.path.join(layout.mddb_registry(base), "glossary_inputs.json"),
                      {"digest": ndp.glossary_inputs(doc),
                       "keys": list(ndp.GLOSSARY_INPUTS)}, indent=1)

    states = g.state.value_counts().to_dict()
    print(f"\n  Glossary: {len(g)} MD(s) -- {states.get('ok', 0)} ok, "
          f"{states.get('review', 0)} review, {states.get('blocked', 0)} blocked")
    print(f"      {layout.glossary_md(base)}")
    print(f"      {os.path.join(layout.results(base), 'glossary.parquet')}\n")
    return 0


def _coverage(base, cover_all, cfg_found, cfg_total):
    """What could not be established, and why. Hidden: it is a working note."""
    L = ["# What the glossary could not establish", "",
         f"Generated {time.strftime('%Y-%m-%dT%H:%M:%S')}. "
         f"{cfg_found} of {cfg_total} run directories had a readable MD config.",
         ""]
    for did, cover in cover_all:
        L += [f"## {did}", "", "| field | status | why |", "|---|---|---|"]
        for field, status, note in cover:
            L.append(f"| `{field}` | {status} | {note} |")
        L.append("")
    out = os.path.join(layout.internal(base), "coverage.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    atomic.write_text(out, "\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
