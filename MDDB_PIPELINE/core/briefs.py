#!/usr/bin/env python3
"""
briefs.py -- the one page a person reads: MD_RESULTS/README.md.

What was found, what has to be checked, and what to do next -- every number
read back out of the built Parquet, so the page cannot drift from the data.

Design rules, learned from what actually goes wrong here:
  * Lead with what NOT to trust. A number without its caveat gets misused.
  * State how each number was obtained, not just its value.
  * Name the gaps explicitly. "unknown" and null are information.
  * Give copy-pasteable queries, not prose descriptions of queries.
"""

import datetime as dt
import os

from . import config

import pandas as pd


def _load(reg, name):
    p = os.path.join(reg, f"{name}.parquet")
    if not os.path.exists(p):
        return None
    try:
        return pd.read_parquet(p)
    except Exception:                                        # noqa: BLE001
        return None


def _fmt_gb(x):
    return f"{x/1000:.2f} TB" if x >= 1000 else f"{x:.1f} GB"


def gaps_of(subst, cfg, runs, docs=None, tabs=None):     # noqa: ARG001
    """What a person has to act on, worst first.

    Read from the index, so the page always states what the data says.
    """
    gaps = []
    if runs is not None and len(runs):
        u = int((runs.system.astype(str) == "unknown").sum())
        if u:
            # Name the exact key. This line is written FOR AN AGENT TO ACT ON,
            # and "add regex rules to substrate.yaml" sent a reader looking for
            # a setting whose name it never gave.
            gaps.append(f"{u} of {len(runs)} runs have no declared `system` — if the "
                        f"paths encode one, add a `systems:` block to "
                        f"`settings.yaml` under `folder:` mapping a label to a "
                        f"regex, then rebuild")
        m = int((runs.n_frames_method == "unavailable").sum())
        if m:
            gaps.append(f"{m} runs have no frame count (format has no cheap header)")
        # A 0-byte trajectory is a job that died before writing a frame. That is
        # a failed RUN, not an indexing problem, and it is worth seeing without
        # having to run the acceptance suite to find it.
        if "header_errors" in runs.columns:
            _e = runs[runs.header_errors.astype(str).str.contains(
                "empty file", na=False)]
            if len(_e):
                gaps.append(
                    f"**{len(_e)} run(s) contain a 0-byte trajectory** — a job "
                    f"that produced nothing: "
                    + ", ".join(f"`{p}`" for p in list(_e.rel_path)[:4]))
    try:
        unc = [n for n, st, _ in config.coverage(cfg) if st == "UNCOVERED"]
    except Exception:                                            # noqa: BLE001
        unc = []
    if unc:
        gaps.insert(0, f"**{len(unc)} top-level director"
                       f"{'y is' if len(unc) == 1 else 'ies are'} neither indexed "
                       f"nor pruned — nothing here describes them: "
                       + ", ".join(f"`{u}`" for u in unc[:8])
                       + ". Add them to `corpora:` or `prune:` under `folder:` "
                         "in settings.yaml.")
    return gaps or ["none detected"]


def build_readme(folder, cfg, summary, out_path=None):
    """MD_RESULTS/README.md -- the one page a person opens.

    What was found, what to do about it, and where the five files are. The long
    long diagnostics stay in .md_pipeline/ (coverage.md): this page names them
    but does not repeat them, because a front page nobody finishes reading is
    the same as no front page.
    """
    from . import layout                                    # noqa: PLC0415
    base = layout.base(folder)
    subst = layout.index_dir(base)
    reg = os.path.join(subst, "registry")
    runs = _load(reg, "runs")
    gl = _load(os.path.join(layout.mddb_dir(base), "registry"), "mds")
    L = []
    A = L.append

    A(f"# {os.path.basename(base)} — MD results")
    A("")
    A(f"Everything this pipeline found under `{base}`, and what is left to do.")
    A(f"Owner: {cfg.get('owner')}. Generated "
      f"{dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}.")
    A("")

    A("## What is here")
    A("")
    if runs is not None and len(runs):
        known = runs[runs.prod_ns_header.notna()]
        A(f"- **{len(runs)} run directories** indexed, "
          f"{_fmt_gb(runs.traj_gb.sum())} of trajectory"
          + (f", **{known.prod_ns_header.sum()/1000:.2f} µs** of production sampling"
             if len(known) else ""))
    if gl is not None and len(gl):
        states = gl.state.value_counts().to_dict()
        A(f"- **{len(gl)} MDs** in the glossary — {states.get('ok', 0)} ok, "
          f"{states.get('review', 0)} to review, {states.get('blocked', 0)} blocked")
        if "nominal_ns" in gl.columns and gl.nominal_ns.notna().any():
            for length, n in sorted(gl.nominal_ns.dropna().value_counts().items()):
                A(f"    - {n} of {_ns(length)}")
    A("")

    A("## Your files")
    A("")
    A("| file | what it is |")
    A("|---|---|")
    A("| `settings.yaml` | **the one file to edit**: what is indexed, the science "
      "behind each dataset, and the catalogue record |")
    A("| `glossary.md` | every MD, with how each value was established |")
    A("| `glossary.parquet` | the same, for scripts |")
    A("| `catalog.json` | every record of this folder, for the lab catalogue |")
    A("| `NDP_registration.md` | what to paste into the NDP form, and what is missing |")
    A("| `*.json` | the record to host and link from the NDP catalogue |")
    A("| `standard_formats/` | each record again as schema.org/DCAT JSON-LD and "
      "as an RO-Crate, and the folder as one DCAT catalogue |")
    A("| `README.md` | this page |")
    A("")
    A(f"The machinery — the index, the caches and the logs — is in "
      f"`{layout.INTERNAL}/`. Nothing in it needs opening, and it can be "
      f"deleted and rebuilt.")
    A("")

    records = _records(base)
    if records:
        A("## The records")
        A("")
        A("| record | identifier | replicas | sampling | volume |")
        A("|---|---|---|---|---|")
        for r in records:
            done = r.get("complete_replicas")
            reps = (f"{r.get('replicas')}"
                    + (f" ({done} complete)" if done is not None else ""))
            A(f"| `{r.get('record')}` | {r.get('identifier') or '—'} | {reps} | "
              f"{_ns(r['total_ns']) if r.get('total_ns') else '—'} | "
              f"{r.get('size_gb') or '—'} GB |")
        A("")
        A("`catalog.json` holds the same rows, for the lab catalogue.")
        A("")

    A("## What to check")
    A("")
    for g in gaps_of(subst, cfg, runs):
        A(f"- {g}")
    if gl is not None and len(gl) and "issues" in gl.columns:
        counts = (gl.issues.fillna("").str.split(", ").explode()
                  .replace("", None).dropna().value_counts())
        for issue, n in counts.items():
            meaning = _issue_meaning(issue)
            A(f"- {n} MD(s) flagged `{issue}`"
              + (f" — {meaning}" if meaning else "")
              + " (see `glossary.md`)")
    A("")

    A("## Next")
    A("")
    A("1. Fill in what only you know in `settings.yaml` — the paper, the "
      "construct, the glycan profile, the ions, and under `ndp: defaults:` who "
      "may get the data (`access_level`) and whether it may train AI models "
      "(`ai_training`).")
    A("2. Re-run the pipeline; every measured value is refreshed and your text "
      "is left alone.")
    A("3. Host the `.json` record and register it: `NDP_registration.md` has the "
      "text for every field of the form.")
    A("")
    A("## For scripts and agents")
    A("")
    A("```python")
    A("import pandas as pd")
    A(f'runs = pd.read_parquet("{os.path.join(reg, "runs.parquet")}")')
    A(f'mds  = pd.read_parquet("{os.path.join(layout.mddb_dir(base), "registry", "mds.parquet")}")')
    A("```")
    A("")
    A(f"- `{os.path.join(layout.results(base), 'glossary.parquet')}` — the "
      f"glossary as a table: one row per MD, every measured value")
    A(f"- `{os.path.join(layout.internal(base), 'coverage.md')}` — what could "
      f"not be established, and why")
    A("")

    out = out_path or layout.readme(base)
    return layout.write_page(out, "\n".join(L) + "\n")


def _records(base):
    """The records this folder has, from catalog.json -- the file that already
    holds exactly what the catalogue shows, so the page cannot disagree."""
    from . import layout                                    # noqa: PLC0415
    import json                                             # noqa: PLC0415
    try:
        with open(os.path.join(layout.results(base), "catalog.json")) as fh:
            return json.load(fh).get("datasets") or []
    except Exception:                                       # noqa: BLE001
        return []


def _ns(value):
    """1000.0 -> '1 µs', 500.0 -> '500 ns'."""
    v = float(value)
    return f"{v/1000:g} µs" if v >= 1000 else f"{v:g} ns"


def _issue_meaning(issue):
    """One line for a glossary issue, from the glossary's own vocabulary.

    The vocabulary lives with the glossary; a page built before one exists
    costs the explanation, not the page.
    """
    try:
        from .glossary import ISSUES                          # noqa: PLC0415
        return ISSUES.get(issue, ("", ""))[1]
    except Exception:                                        # noqa: BLE001
        return ""
