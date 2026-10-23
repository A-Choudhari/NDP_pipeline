#!/usr/bin/env python3
"""
glossary_query.py -- look up MDs in the glossary and hand their inputs to analyses.

    python3 glossary_query.py <folder> summary
    python3 glossary_query.py <folder> list  [filters] [--columns a,b,c]
    python3 glossary_query.py <folder> show  <md_key>
    python3 glossary_query.py <folder> paths [md_key ...] [filters] [--json]

`<folder>` is your data folder, or the lab-wide glossary directory
(`<folder>`), whose keys are
`<folder>/<md_key>` and can be given to `show` and `paths` as such.

Filters (combine freely; a comma means "any of"):

    --dataset gp120-open,gp120-occluded-open
    --system  gp120-gp41-glyco
    --condition HD4
    --state   ok            (ok | review | blocked)
    --issue   ordering_unverified
    --where   "temp_K == 310 and simulated_ns > 400"

`paths` prints, per MD, the topology and the trajectory segments in time order
as absolute paths -- what an analysis script needs:

    python3 glossary_query.py . paths --system gp120-occluded-open --state ok \\
        --links --json > selection.json

Without --links the paths are the original files under their original names.
With --links they go through the read-only links in MDDB_substrate/work/, which
reach the same bytes but carry each file's TRUE extension: two LC3B trajectories
are NetCDF files named .mdcrd, and tools that pick a reader by extension fail on
the original names. Prefer --links for analysis.

Those files are for reading. Nothing here writes anywhere; the glossary itself is
built by glossary_build.py. By default `paths` refuses MDs whose state is `blocked`,
because their inputs demonstrably do not fit together; add --include-blocked to
override, and it will say so on stderr.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import glossary                      # noqa: E402

DEFAULT_LIST = ["md_key", "state", "system", "condition", "natom",
                "n_segments", "simulated_ns", "framestep_ns", "temp_K",
                "ensemble", "issues"]


def _split(v):
    return [x.strip() for x in v.split(",") if x.strip()] if v else None


def _selection(g, a, keys=None):
    return glossary.select(g, dataset=_split(a.dataset), system=_split(a.system),
                           condition=_split(a.condition), state=_split(a.state),
                           issue=_split(a.issue), where=a.where, keys=keys)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", help="your data folder (the one holding MDDB_substrate/)")
    ap.add_argument("command", choices=("summary", "list", "show", "paths"))
    ap.add_argument("keys", nargs="*", help="md_key(s), e.g. gp120-open/CDR_1-1F4_0-MD_1-R1")
    for f in ("dataset", "system", "condition", "state", "issue"):
        ap.add_argument(f"--{f}")
    ap.add_argument("--where", help="pandas query over glossary columns")
    ap.add_argument("--columns", help="comma-separated columns for `list`")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--include-blocked", action="store_true",
                    help="let `paths` return MDs whose inputs do not fit together")
    ap.add_argument("--links", action="store_true",
                    help="`paths` through the read-only work links, whose "
                         "extensions match each file's real format")
    a = ap.parse_args(argv)

    try:
        g = glossary.load(a.folder)
    except FileNotFoundError as exc:
        print(f"\n  {exc}\n", file=sys.stderr)
        return 2

    import pandas as pd
    pd.set_option("display.width", 220)
    pd.set_option("display.max_colwidth", 70)
    pd.set_option("display.max_rows", 500)

    if a.command == "summary":
        if "member" in g.columns and not a.json:
            # The lab-wide glossary: one line per folder first.
            bym = g.groupby("member").agg(
                owner=("owner", "first"),
                MDs=("md_key", "size"),
                ok=("state", lambda s: int((s == "ok").sum())),
                review=("state", lambda s: int((s == "review").sum())),
                blocked=("state", lambda s: int((s == "blocked").sum())),
                run=("member_status", "first"),
                built=("member_built_at", "first"))
            print(f"\n{len(g)} MDs from {g.member.nunique()} folders\n")
            print(bym.to_string())
        by = g.groupby("dataset_id").agg(
            MDs=("md_key", "size"),
            ok=("state", lambda s: int((s == "ok").sum())),
            review=("state", lambda s: int((s == "review").sum())),
            blocked=("state", lambda s: int((s == "blocked").sum())),
            simulated_us=("simulated_ns", lambda s: round(s.sum() / 1000, 2)),
            trajectory_gb=("trajectory_gb", lambda s: round(s.sum(), 1)),
            topology_groups=("topology_group", "nunique"))
        if a.json:
            print(by.reset_index().to_json(orient="records", indent=2))
        else:
            print(f"\n{len(g)} MDs in {g.dataset_id.nunique()} datasets\n")
            print(by.to_string())
            counts = {}
            for s in g.issues.fillna(""):
                for i in [x.strip() for x in s.split(",") if x.strip()]:
                    counts[i] = counts.get(i, 0) + 1
            if counts:
                print("\nissues:")
                for k, n in sorted(counts.items(), key=lambda kv: -kv[1]):
                    sev, text = glossary.ISSUES.get(k, ("?", ""))
                    print(f"  {n:4d}  {k:22s} [{sev}] {text}")
            print()
        return 0

    sel = _selection(g, a, keys=a.keys or None)
    if a.keys:
        known = set(sel.md_key)
        if "lab_key" in sel.columns:
            known |= set(sel.lab_key)
        missing = set(a.keys) - known
        if missing:
            print(f"  unknown md_key(s): {sorted(missing)}", file=sys.stderr)
            return 2

    if a.command == "list":
        cols = _split(a.columns) or DEFAULT_LIST
        bad = [c for c in cols if c not in g.columns]
        if bad:
            print(f"  unknown column(s) {bad}. Available: {list(g.columns)}",
                  file=sys.stderr)
            return 2
        if a.json:
            print(sel[cols].to_json(orient="records", indent=2))
        else:
            print(sel[cols].to_string(index=False) if len(sel) else "  no MDs match")
            print(f"\n  {len(sel)} of {len(g)} MDs")
        return 0

    if a.command == "show":
        if len(sel) != 1:
            print("  `show` needs exactly one md_key", file=sys.stderr)
            return 2
        row = sel.iloc[0]
        if a.json:
            print(row.to_json(indent=2))
            return 0
        for col in g.columns:
            print(f"  {col:15s} {row[col]}")
        for i in [x.strip() for x in str(row.issues or "").split(",") if x.strip()]:
            sev, text = glossary.ISSUES.get(i, ("?", ""))
            print(f"  -> {i}: {text} [{sev}]")
        return 0

    # paths
    blocked = sel[sel.state == "blocked"]
    if len(blocked) and not a.include_blocked:
        print(f"  skipping {len(blocked)} blocked MD(s) whose inputs do not fit "
              f"together (use --include-blocked to override): "
              f"{list(blocked.md_key[:5])}", file=sys.stderr)
        sel = sel[sel.state != "blocked"]
    elif len(blocked):
        print(f"  WARNING: including {len(blocked)} blocked MD(s)", file=sys.stderr)
    review = sel[sel.state == "review"]
    if len(review):
        print(f"  note: {len(review)} MD(s) are in `review`; check their issues "
              f"before trusting results", file=sys.stderr)

    relabelled = sel[sel.issues.fillna("").str.contains("relabelled_format")]
    if len(relabelled) and not a.links:
        print(f"  note: {len(relabelled)} MD(s) have a segment whose extension "
              f"misstates its format; add --links to get paths tools read "
              f"correctly: {list(relabelled.md_key[:3])}", file=sys.stderr)

    out = []
    for _, md in sel.iterrows():
        top, trajs = glossary.input_paths(md, via="links" if a.links else "source")
        out.append({"md_key": md.md_key, "state": md.state, "issues": md.issues,
                    "topology": top, "trajectories": trajs,
                    "formats": [f for f in str(md.formats or "").split("; ") if f]})
    if a.json:
        print(json.dumps(out, indent=2))
    else:
        for rec in out:
            print(f"# {rec['md_key']}  [{rec['state']}]")
            print(f"topology\t{rec['topology']}")
            for t in rec["trajectories"]:
                print(f"trajectory\t{t}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
