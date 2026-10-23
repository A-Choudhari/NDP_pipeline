#!/usr/bin/env python3
"""
ndp_build.py -- one NDP metadata record (JSON) per simulated system.

    python3 ndp_build.py /path/to/<folder> [--datasets id,id] [--check]

The National Data Platform catalogue holds metadata only. The lab's registration
guide links each catalogue entry to one small JSON file per simulated system,
with the same field names in every file. This script writes those files from
the MD glossary, so nothing in them is typed twice or from memory.

Reads  the MD glossary (glossary_build.py) and MD_RESULTS/settings.yaml: the
       dataset's contact, licence and citation, and the `ndp:` section.

Writes MD_RESULTS/settings.yaml      first run: one empty block per system under
                                     `ndp: systems:`. Later runs: new systems
                                     appended, nothing else touched. FILL IT IN,
                                     then re-run.
       MD_RESULTS/<name>.json        the record to host and link (more than
                                     three: MD_RESULTS/NDP_records/)
       MD_RESULTS/NDP_registration.md  per record, the text to paste into the
                                     NDP form, and everything still missing
       MD_RESULTS/standard_formats/  every record again as schema.org/DCAT
                                     JSON-LD and as an RO-Crate, plus the folder
                                     as one DCAT catalogue (core/standard.py)

A system is the MDs of one dataset that share a topology; they are its replicas.
Blocked MDs are not counted (--include-blocked to override). Values the MDs of
a system disagree on are left null and reported, never settled by majority.

Writes only inside <folder>/MDDB_substrate/. It reads the MD config named in
the glossary (for the target pressure) and nothing else of your data. It
uploads nothing: registration is done by hand on nationaldataplatform.org, and
records stay private until the paper they support is published.

--check writes nothing and exits 1 if any record is incomplete or has warnings.
"""

import argparse
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import (atomic, catalog, glossary, lab, layout, ndp, protocol,  # noqa: E402
                  settings, standard)


def _read(path):
    with open(path) as fh:
        return fh.read()


def _is_record(path):
    """True for a file this pipeline wrote as a record, not any old JSON.

    Belt and braces beside the top-level rule: a record always carries its
    `record_id` in `lab_glossary`, so anything else in the folder -- a note, a
    settings export, someone's scratch file -- is left where it is.
    """
    try:
        doc = json.loads(_read(path))
    except (ValueError, OSError):
        return False
    return isinstance(doc, dict) and isinstance(doc.get("lab_glossary"), dict) \
        and "record_id" in doc["lab_glossary"]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=os.getcwd())
    ap.add_argument("--datasets", help="comma-separated ids; default all")
    ap.add_argument("--include-blocked", action="store_true",
                    help="count MDs whose topology and trajectory do not fit together")
    ap.add_argument("--check", action="store_true",
                    help="write nothing; exit 1 unless every record is ready")
    a = ap.parse_args(argv)

    base = layout.base(a.folder)
    try:   # one build per folder at a time: concurrent runs clobber and sweep
        atomic.hold(layout.lock(base, "build"), what=base)
    except atomic.Locked as exc:
        print(f"\n  Refused: {exc}\n")
        return 2
    try:
        g = glossary.load(base)
    except FileNotFoundError:
        print(f"\n  No MD glossary under {base}")
        print(f"  Run:  python3 {os.path.join(HERE, 'glossary_build.py')} {base}\n")
        return 2
    if not settings.exists(base):
        print(f"\n  No settings.yaml at {settings.path(base)}\n")
        return 2
    doc_all = settings.load(base)
    specs = {str(d.get("id") or d.get("name")): d
             for d in settings.datasets_of(base, doc_all)}

    if a.datasets:
        keep = {s.strip() for s in a.datasets.split(",")}
        g = g[g.dataset_id.astype(str).isin(keep)]
    if not len(g):
        print("\n  The glossary holds no MDs for that selection.\n")
        return 1

    print(f"\n  folder         : {base}")
    print(f"  glossary       : {len(g)} MDs in {g.dataset_id.nunique()} dataset(s)")
    stale = _glossary_stale(base, doc_all)
    if stale:
        print(f"\n  NOTE: {', '.join(ndp.GLOSSARY_INPUTS)} changed in settings.yaml "
              f"since the glossary\n  was built. The records below still use the "
              f"old values. Rerun the whole pipeline:\n      "
              f"{os.path.join(HERE, 'ndp_pipeline.sh')} {base} --yes\n")

    # Which MDs are production MDs, by length. Read from settings.yaml BEFORE
    # grouping, because it decides what a system's replicas are.
    cur_p = settings.path(base)
    pre = settings.ndp_of(base, doc_all)
    mode = str(ndp.merged_defaults(pre.get("defaults")).get("lengths")
               or "all").strip().lower()
    if mode not in ndp.LENGTH_CHOICES:
        print(f"  `lengths: {mode}` in settings.yaml is not one of "
              f"{', '.join(ndp.LENGTH_CHOICES)}; using `all`.")
        mode = "all"
    short = g.iloc[0:0]
    if mode == "nominal":
        if "nominal_ns" not in g.columns:
            print("  lengths        : this glossary predates nominal_ns -- rebuild it "
                  "with glossary_build.py; every MD is kept for now")
        else:
            is_short = g.nominal_ns.isna() & g.simulated_ns.notna()
            short, g = g[is_short], g[~is_short]
            print(f"  lengths        : {len(g)} MD(s) reach a declared production "
                  f"length; {len(short)} shorter one(s) left out")
    if not len(g):
        print("\n  No MD reaches a declared production length. Set `lengths: all` "
              "under `ndp: defaults:`,\n  or correct `production_lengths_ns:`, in "
              "settings.yaml, and rebuild the glossary.\n")
        return 1

    # --- 1. group into systems and measure each ---------------------------
    all_systems = ndp.systems_of(g, include_blocked=a.include_blocked)
    skipped = [(rid, nb) for rid, rows, nb in all_systems if not len(rows)]
    systems = [(rid, rows, nb) for rid, rows, nb in all_systems if len(rows)]
    measured, warnings = {}, {}
    for rid, rows, nb in systems:
        measured[rid], warnings[rid] = ndp.measure(rows)
        if nb:
            warnings[rid].append(f"{nb} blocked MD(s) of this system are not counted.")
    left_out = {}
    for rid, rows, _ in ndp.systems_of(short, include_blocked=True):
        left_out[rid] = [f"{k} ({ns:g} ns)" for k, ns in zip(rows.md_key, rows.simulated_ns)]
    print(f"  systems        : {len(systems)} (MDs sharing one topology = replicas)")

    # --- 2. the curation file ---------------------------------------------
    if a.check:
        doc = pre
        added, gone = [], []
    else:
        try:
            doc, added, gone = ndp.sync_curation(base, systems, measured)
        except ndp.CurationLayout as exc:
            print(f"\n  {exc}\n")
            return 1
    defaults = ndp.merged_defaults(doc.get("defaults"))

    # --- 3. the records ----------------------------------------------------
    rec_dir = layout.records_dir(base, len(systems))
    index_p = layout.ndp_index(base)
    old = {}
    if os.path.exists(index_p):
        try:
            old = {r["record_id"]: r for r in json.loads(_read(index_p))["records"]}
        except (ValueError, KeyError, TypeError):
            old = {}
    today = time.strftime("%Y-%m-%d")
    entries, sheets, seen_files = [], [], {}
    # rid -> (entry, record), so an incumbent renamed by a later collision is
    # corrected in BOTH places. `taken_ids` is what makes the identifier
    # de-duplication actually run: it used to be fed the entry dicts, which
    # have never carried an `identifier` key, so the list was all None and the
    # `-002` bump never fired. 18 identifiers are duplicated in the lab-wide
    # catalogue because of it.
    by_rid, taken_ids = {}, []
    suggested = {}                   # rid -> {field: value the pipeline guessed}
    notes = {}                       # rid -> [note]: said, but never a warning
    for rid, rows, _ in systems:
        cur = ndp.curation_for(rid, doc)
        spec = specs.get(rid.split("/", 1)[0]) or {}
        warn = warnings[rid]
        notes[rid] = ndp.declared_notes(rows, cur, warn)
        # What this record's own MDs have in common below the folder
        # (`1-1F4`): what tells one system of a folder from another, in its
        # name, its identifier and its tags.
        where = ndp.common_subpath(list(rows.source_dir))
        record = ndp.build_record(rid, rows, measured[rid], cur, spec, defaults, warn,
                                  stages=_protocol(rows))
        # What the lab catalogue files this dataset under. Suggested from the
        # folder's own names on the first run, then whatever the person wrote.
        # A suggestion fills the record exactly as before, but is remembered,
        # so the registration sheet can ask a person to confirm it: an
        # identifier guessed from `TO_SHARE_FOR_PUBLICATION` read `<PREFIX>-COV2-TO-`
        # and counted as filled in (SPIKE_CLEAVED, 2026-09-30).
        guessed = {}
        if not record.get("topic"):
            record["topic"] = catalog.topic_of(base, doc_all.get("folder"), rid)
            if record["topic"]:
                guessed["topic"] = record["topic"]
        if not record.get("identifier"):
            record["identifier"] = catalog.suggest_identifier(
                record["topic"], base, rid, what=where, taken=taken_ids)
            if record["identifier"]:
                guessed["identifier"] = record["identifier"]
        if record.get("identifier"):
            taken_ids.append(record["identifier"])
        if not record.get("tags"):
            record["tags"] = catalog.suggest_tags(record, cur, what=where)
            if record["tags"]:
                guessed["tags"] = ", ".join(map(str, record["tags"]))
        if guessed:
            suggested[rid] = guessed
        keywords, rejected = ndp.keywords_of(cur, defaults)
        warn += ndp.check_text(record)
        if rejected:
            warn.append("keyword(s) NDP will refuse (letters, numbers, spaces, "
                        "- _ . only): " + ", ".join(rejected))
        # The record's own sampling, not the raw measurement: a length
        # recovered from a declared frame interval names the file too.
        rel = ndp.file_name(rid, cur, folder=base, measured=record["sampling"],
                            paths=list(rows.source_dir),
                            existing=(old.get(rid) or {}).get("file"))
        kept = (old.get(rid) or {}).get("file")
        if kept and rel == kept and ndp.slug_would_rename(rid, cur, record["sampling"],
                                                          kept):
            notes[rid].append(f"file name `{kept}` kept although the sampling "
                              f"changed; `rename: true` in its settings.yaml "
                              f"block renames it.")
        # Two systems of one folder with the same sampling: keep both, told
        # apart by the topology id they already carry in the glossary.
        rel, incumbent = _claim_name(rel, rid, seen_files)
        if incumbent and incumbent in by_rid:
            old_entry, old_record = by_rid[incumbent]
            moved, _ = _claim_name(old_entry["file"], incumbent, seen_files)
            old_entry["file"] = moved
            # The record's OWN copy of its name, which catalog.json reads.
            # Missing this is what put dead pointers in the catalogue.
            old_record["record"] = moved
        record["record"] = rel
        text = ndp.dumps(record)
        sha = ndp.digest(text)
        prev = old.get(rid) or {}
        missing = ndp.missing_fields(record)
        entries.append({
            "record_id": rid,
            "file": rel,
            "title": record["title"],
            "status": "incomplete" if missing else "ready",
            "missing": missing,
            "warnings": warn,
            "keywords": keywords,
            # MDs of this system shorter than every declared production length.
            "left_out_short": left_out.get(rid, []),
            "sha256_16": sha,
            "first_written": prev.get("first_written") or today,
            "last_changed": (prev.get("last_changed") or today)
            if prev.get("sha256_16") == sha else today,
        })
        by_rid[rid] = (entries[-1], record)
        sheets.append([entries[-1], record, measured[rid], text])

    # A record renamed by a LATER collision was dumped before it was renamed,
    # so the text and its digest are rebuilt once every name is settled.
    for sheet in sheets:
        entry, record, _, text = sheet
        if record.get("record") != entry["file"]:
            record["record"] = entry["file"]
        fresh = ndp.dumps(record)
        if fresh != text:
            sheet[3] = fresh
            sha = ndp.digest(fresh)
            prev = old.get(entry["record_id"]) or {}
            entry["sha256_16"] = sha
            entry["last_changed"] = ((prev.get("last_changed") or today)
                                     if prev.get("sha256_16") == sha else today)

    # Every field still missing gets an empty line in its own system's block,
    # so a person always finds the key to answer it on.
    lines_added, lines_skipped = ({}, {}) if a.check else ndp.ensure_lines(
        base, {e["record_id"]: e["missing"] for e in entries if e["missing"]})

    ready = sum(1 for e in entries if e["status"] == "ready")
    n_warn = sum(len(e["warnings"]) for e in entries)
    for rid in sorted(notes):
        for n in notes[rid]:
            print(f"  NOTE: {rid}: {n}")
    for n in ndp.swallowed_comments(doc_all):
        print(f"  NOTE: settings.yaml: {n}")
    if suggested:
        print(f"\n  {len(suggested)} record(s) carry a SUGGESTED identifier/topic/tags "
              f"-- confirm them in settings.yaml\n  before registering (listed in "
              f"{os.path.relpath(layout.registration_md(base), base)}).")
    if a.check:
        print(f"\n  {ready}/{len(entries)} record(s) ready, {n_warn} warning(s). "
              f"Nothing was written (--check).\n")
        return 0 if ready == len(entries) and not n_warn else 1

    for entry, _, _, text in sheets:
        path = os.path.join(rec_dir, entry["file"])
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not (os.path.exists(path) and _read(path) == text):
            atomic.write_text(path, text)

    # A record whose system is gone, or whose slug changed, must not sit in
    # records/ looking current. These are our own regenerable output, so they
    # are removed -- loudly, and only on a full build.
    if not a.datasets:
        wanted = {os.path.join(rec_dir, e["file"]) for e in entries}
        # catalog.json is the folder's index, not a record; it lives beside
        # them and is rewritten by catalog.write.
        wanted.add(os.path.join(layout.results(base), "catalog.json"))
        # ONLY the top level, and only files that ARE records. This walked
        # the whole tree, and `rec_dir` is MD_RESULTS itself -- so a full
        # build would have deleted every .json under `analyses/` that it had
        # not written: 39 files on 2026-09-24, including each system's
        # `frames.json` (the time axis every plot's frame provenance is read
        # from) and 35 stage manifests. Nothing under `analyses/` belongs to
        # this pipeline, and a sweep that removes what it did not write is
        # not a sweep.
        # BOTH locations, because the one in use CHANGES. Records live beside
        # the other files until there are more than three, then they move into
        # NDP_records/. A folder whose count crosses that line leaves its old
        # copies behind in the location it stopped using, and nothing looked
        # there: 1-HD5, 3-HD3 and 4-HD2 each carried a full set of superseded
        # records in NDP_records/ after their record count fell, invisible to
        # the catalogue and counted by anything that globbed the tree.
        other = (os.path.join(layout.results(base), "NDP_records")
                 if rec_dir == layout.results(base) else layout.results(base))
        for d in (rec_dir, other):
            if not os.path.isdir(d):
                continue
            for name in sorted(os.listdir(d)):
                path = os.path.join(d, name)
                if not name.endswith(".json") or path in wanted:
                    continue
                if not os.path.isfile(path) or not _is_record(path):
                    continue
                os.remove(path)
                where = "" if d == rec_dir else " (from the location no longer in use)"
                print(f"        removed stale record '{name}'{where}")
        # An emptied NDP_records/ is our own directory and should not linger.
        empty = os.path.join(layout.results(base), "NDP_records")
        if os.path.isdir(empty) and not os.listdir(empty):
            os.rmdir(empty)
    else:
        # A partial build replaces its own datasets' entries and keeps the rest.
        done = {e["record_id"] for e in entries}
        entries = sorted([r for rid, r in old.items() if rid not in done] + entries,
                         key=lambda e: e["record_id"])

    catalog.write(os.path.join(layout.results(base), "catalog.json"), base,
                  [rec for _, rec, _, _ in sheets],
                  base_url=defaults.get("base_url"))
    # The same records in the formats other catalogues read natively.
    for name in standard.write_all(layout.results(base), base,
                                   [(e, rec, text) for e, rec, _, text in sheets],
                                   base_url=defaults.get("base_url"),
                                   lab_info=lab.info(), sweep=not a.datasets):
        print(f"        removed stale standard_formats/{name}/")
    # The front page is written by the index, before any of this existed; now
    # that the glossary and the records do, write it again so it says so.
    try:
        from core import briefs
        briefs.build_readme(base, settings.config_of(base, doc_all), {})
    except Exception as exc:                                 # noqa: BLE001
        print(f"  (README not refreshed: {type(exc).__name__}: {exc})")
    os.makedirs(os.path.dirname(index_p), exist_ok=True)
    atomic.write_json(index_p, {"template": "NDP metadata registration guide (2026-09-21)",
                                "base_url": defaults.get("base_url"),
                                "records": entries}, indent=2, ensure_ascii=False)
    only_short = {rid: keys for rid, keys in left_out.items()
                  if rid not in {e["record_id"] for e in entries}}
    _registration_report(base, sheets, defaults, skipped, gone, cur_p, only_short,
                         suggested=suggested, notes=notes)

    print()
    for rid, keys in sorted(lines_added.items()):
        print(f"  {rid}: empty line(s) added to its block in settings.yaml for "
              f"{', '.join(keys)}")
    for rid, why in sorted(lines_skipped.items()):
        print(f"  NOTE: {rid}: no empty lines added to its block in settings.yaml "
              f"({why}); type the missing keys there by hand.")
    if added:
        print(f"  {len(added)} new system(s) added to "
              f"{os.path.relpath(cur_p, base)} -- fill them in and re-run.")
    for rid in gone:
        print(f"  NOTE: {rid} is listed in settings.yaml but no longer in the glossary.")
    for rid, nb in skipped:
        print(f"  NOTE: {rid} has only blocked MDs ({nb}); no record written.")
    for rid, keys in only_short.items():
        print(f"  NOTE: {rid}: every MD is shorter than the production lengths "
              f"({', '.join(keys[:3])}{' ...' if len(keys) > 3 else ''}); no record.")
    for e in entries:
        if e["left_out_short"]:
            print(f"  {e['record_id']}: left out as shorter than the production "
                  f"lengths: {', '.join(e['left_out_short'])}")
    if n_warn:
        print(f"\n  {n_warn} WARNING(S) -- fix before registering:")
        for e in entries:
            for w in e["warnings"]:
                print(f"    {e['record_id']}: {w}")
    print(f"\n  Done. {len(sheets)} record(s) under "
          f"{os.path.relpath(rec_dir, base)}/ -- {ready} ready, "
          f"{len(sheets) - ready} incomplete.")
    print(f"  For the lab catalogue    : "
          f"{os.path.relpath(os.path.join(layout.results(base), 'catalog.json'), base)}")
    print(f"  schema.org/DCAT, RO-Crate: "
          f"{os.path.relpath(os.path.join(layout.results(base), standard.DIR), base)}/")
    print(f"  What to paste into the NDP form, and what is missing: "
          f"{os.path.relpath(layout.registration_md(base), base)}")
    print(f"  Nothing was uploaded. Host the record(s) (GitHub or Zenodo), set "
          f"`base_url` under `ndp: defaults:` in settings.yaml,\n  and keep "
          f"entries Unlisted until the paper they support is published.\n")
    return 0


def _glossary_stale(base, doc):
    """True when the settings the glossary was built from have changed since.

    False when it cannot tell -- a glossary built before the stamp existed --
    so a folder built earlier says nothing new until its next full run."""
    p = os.path.join(layout.mddb_registry(base), "glossary_inputs.json")
    try:
        built = json.loads(_read(p)).get("digest")
    except (OSError, ValueError, AttributeError):
        return False
    return bool(built) and built != ndp.glossary_inputs(doc)


def _protocol(rows):
    """The merged protocol of one system's replicas, keyed by replica label."""
    labels = [str(r) for r in rows.replica_index] \
        if "replica_index" in rows.columns else []
    if len(set(labels)) != len(rows) or "None" in labels:
        labels = [str(m) for m in rows.md_id]
    per, whole = {}, {}
    for lab, r in zip(labels, rows.itertuples()):
        root, rel = getattr(r, "source_root", None), getattr(r, "source_dir", None)
        if root and rel:
            per[lab] = protocol.stages_of(os.path.join(root, rel))
        whole[lab] = _production_trajectory(r)
    stages = protocol.merge({k: v for k, v in per.items() if v})
    return [_whole_production(s, whole) for s in stages]


def _production_trajectory(r):
    """What a replica's production really is: every segment, not the first.

    A stage describes the file its config's `outputName` names -- `pr.part1`
    for a run that restarted eight times, which reads as if production were
    244 frames of 2.25 GB. Now that the record publishes the production stage
    and nothing else (2026-09-22), that single segment IS the reader's whole
    picture, so it is replaced by the glossary's count for the run.
    """
    segs = [s for s in str(getattr(r, "segments", "") or "").split("; ") if s]
    out = {"segments": int(getattr(r, "n_segments", 0) or len(segs)) or None,
           "files": (segs[0] if len(segs) == 1 else
                     f"{segs[0]} ... {segs[-1]}" if segs else None),
           "frames": _int(getattr(r, "frames", None)),
           "atoms": _int(getattr(r, "natom", None)),
           "frame_interval_ps": _int(round(float(r.framestep_ns) * 1000))
           if getattr(r, "framestep_ns", None) == getattr(r, "framestep_ns", None)
           and getattr(r, "framestep_ns", None) else None,
           "size_gb": round(float(r.trajectory_gb), 2)
           if getattr(r, "trajectory_gb", None) == getattr(r, "trajectory_gb", None)
           and getattr(r, "trajectory_gb", None) else None,
           "length_ns": round(float(r.simulated_ns), 1)
           if getattr(r, "simulated_ns", None) == getattr(r, "simulated_ns", None)
           and getattr(r, "simulated_ns", None) else None}
    return {k: v for k, v in out.items() if v is not None}


def _int(v):
    try:
        return None if v is None or v != v else int(v)
    except (TypeError, ValueError):
        return None


def _whole_production(stage, whole):
    """Give a production stage the run's own trajectory numbers."""
    if stage.get("stage") != "production" or not whole:
        return stage
    kept = {k: v for k, v in whole.items() if v}
    if not kept:
        return stage
    stage = dict(stage)
    # `outputName pr.part1` names the window, not the stage a reader means.
    stage["name"] = re.sub(r"\.part\d+$", "", str(stage.get("name") or "")) or "pr"
    stage["trajectory"] = (list(kept.values())[0] if len(kept) == 1
                           else {"per_replica": kept})
    return stage


def _claim_name(rel, rid, claimed):
    """A record file name nothing else holds, and the incumbent to rename.

    `claimed` maps every name ever handed out to the record that took it, and
    NOTHING IS EVER REMOVED FROM IT. That is the fix. The previous version
    popped the bare name when two records collided, so a THIRD record wanting
    it found it free and took it -- while the record that had been renamed went
    on carrying the old name inside its own JSON, because only its index entry
    was updated. `catalog.json` reads the name out of the record, so the
    catalogue pointed at a file that was never written: on 2026-09-26 one stem
    appeared four times in the lab-wide catalogue and all four pointers were
    dead.

    The bare name is burned on the first collision -- neither the incumbent nor
    anyone later keeps it -- so every record of a colliding group is told apart
    by the topology id it already carries, rather than one of them silently
    keeping the ambiguous name.
    """
    if rel not in claimed:
        claimed[rel] = rid
        return rel, None
    incumbent = claimed[rel]
    stem, gid = rel[:-5], rid.split("/", 1)[1]
    name, n = f"{stem}_{gid}.json", 2
    while name in claimed:
        name, n = f"{stem}_{gid}_{n}.json", n + 1
    claimed[name] = rid
    claimed[rel] = None          # burned: nobody holds the bare name now
    return name, incumbent


def _publisher_line():
    org = lab.info()["ndp_organization"]
    if org:
        return f"Publisher: **{org}**."
    who = lab.short_name()
    return ("Publisher: **Individual** until NDP lists "
            + (f"an organization for the {who}." if who else "your lab's organization.")
            + " (Set `ndp_organization:` in lab.yaml once it does.)")


def _owner_text(owner):
    owner = owner or {}
    if not owner.get("name"):
        return None
    extra = ", ".join(x for x in (owner.get("affiliation"), owner.get("email"),
                                  owner.get("orcid")) if x)
    return owner["name"] + (f" ({extra})" if extra else "")


def _location_text(loc):
    loc = loc or {}
    text = ": ".join(f"`{x}`" if x == loc.get("path") else x
                     for x in (loc.get("storage"), loc.get("path")) if x)
    if loc.get("uri"):
        text += f" — {loc['uri']}" if text else loc["uri"]
    return text or None


def _registration_report(base, sheets, defaults, skipped, gone, cur_p,
                         only_short=None, suggested=None, notes=None):
    """Per record: the NDP form, field by field, and what still blocks it."""
    base_url = str(defaults.get("base_url") or "").rstrip("/")
    ready = [s for s in sheets if s[0]["status"] == "ready"]
    L = ["# NDP registration sheet", "",
         f"Generated {time.strftime('%Y-%m-%dT%H:%M:%S')} by ndp_build.py. "
         f"{len(sheets)} record(s): {len(ready)} ready, "
         f"{len(sheets) - len(ready)} incomplete.", "",
         "One section per simulated system, in the order of the form at "
         "nationaldataplatform.org (Catalog, Add to Catalog). Visibility: "
         "**Unlisted** until the supporting paper is published. Asset type: "
         "**Dataset**. " + _publisher_line(), "",
         f"Missing values are filled in `{cur_p}`: the `ndp:` section for the "
         "record, the `datasets:` section for contact, licence and citation. "
         "Data owner, access level, AI training and the storage name fall back "
         f"to the lab's own values in `{lab.path()}`. "
         "Then re-run ndp_build.py.", ""]
    if not base_url:
        L += ["**`base_url` is not set** under `ndp: defaults:` in settings.yaml, "
              "so the item URLs below are relative. Host the records first, then "
              "set it.", ""]
    L += ["| record | status | replicas | total ns | missing |", "|---|---|---|---|---|"]
    for e, rec, m, _ in sheets:
        L.append(f"| `{e['record_id']}` | {e['status']} | {m['replicas']} | "
                 f"{m['total_ns'] if m['total_ns'] is not None else '—'} | "
                 f"{', '.join(e['missing']) or '—'} |")
    L.append("")
    for rid, nb in skipped:
        L.append(f"- `{rid}`: only blocked MDs ({nb}); no record written.")
    for rid in gone:
        L.append(f"- `{rid}`: listed in settings.yaml but no longer in the glossary.")
    for rid, keys in (only_short or {}).items():
        L.append(f"- `{rid}`: every MD is shorter than the production lengths "
                 f"({', '.join(keys)}); no record written.")
    if skipped or gone or only_short:
        L.append("")
    if suggested:
        # Status is not changed by these; the lab catalogue is the registry of
        # identifiers, so a person copies each value (or a better one) onto the
        # key already in settings.yaml, and it stops being listed here.
        L += ["## Suggested by the pipeline: confirm before registering", "",
              "These values were guessed from folder names, not written by a "
              "person. Type the value you want on the matching key of that "
              f"record's block under `ndp: systems:` in `{cur_p}`, then re-run "
              "ndp_build.py.", "",
              "| record | field | suggested value |", "|---|---|---|"]
        for rid in sorted(suggested):
            for field, value in suggested[rid].items():
                L.append(f"| `{rid}` | `{field}` | {value} |")
        L.append("")

    for e, rec, m, _ in sheets:
        # `file` is a bare name, or `NDP_records/<name>` once a folder holds
        # more than three records. Splitting unconditionally raised IndexError
        # on every folder with three or fewer, so setting `base_url` -- the one
        # thing you do just before registering -- killed the final step.
        url = (f"{base_url}/{e['file'].split('/', 1)[-1]}" if base_url
               else e["file"])
        L += [f"## {e['record_id']} — {e['status']}", ""]
        if e["missing"]:
            L += ["Missing: " + ", ".join(f"`{f}`" for f in e["missing"]), ""]
        if e["left_out_short"]:
            L += ["Left out, shorter than the production lengths: "
                  + ", ".join(f"`{k}`" for k in e["left_out_short"]), ""]
        for w in e["warnings"]:
            L.append(f"- WARNING: {w}")
        for n in (notes or {}).get(e["record_id"], []):
            L.append(f"- Note: {n}")
        if e["warnings"] or (notes or {}).get(e["record_id"]):
            L.append("")
        L += ["| form field | value |", "|---|---|",
              f"| Title | {ndp.catalog_title(rec, m) or '**missing: title**'} |",
              f"| Description | {ndp.catalog_description(rec, m)} |",
              f"| Keywords | {', '.join(e['keywords'])} |",
              f"| Contact | {rec['contact'] or '**missing**'} |",
              "| Item: Title | MD metadata record (JSON) |",
              f"| Item: Description | Machine-readable metadata for this simulation. "
              f"{str(rec['data_availability'] or '').rstrip('.')}. |",
              "| Item: Format / Access | JSON / Direct Download |",
              f"| Item: URL | {url} |",
              f"| Item: Issue Date / Last Updated | {e['first_written']} / {e['last_changed']} |",
              f"| License | {rec['license'] or '**missing**'} |",
              f"| Data owner (PI) | {_owner_text(rec.get('data_owner')) or '**missing**'} |",
              f"| Access level | {rec.get('access_level') or '**missing**'} |",
              f"| AI training | {rec.get('ai_training') or '**missing**'} |",
              f"| Raw data location | {_location_text(rec.get('data_location')) or '**missing**'} |",
              f"| MDs | {', '.join(f'`{k}`' for k in rec['lab_glossary']['md_keys'])} |",
              ""]
    layout.write_page(layout.registration_md(base), "\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
