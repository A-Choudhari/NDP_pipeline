#!/usr/bin/env python3
"""
glossary.py -- one stable, named entry per MD, to run analyses against.

WHAT THE GLOSSARY IS FOR. ~100 TB of simulation lives under directory names
that were chosen for the person running the job, not for anyone reading them
later: `pr`, `rep2`, `REP1`, `Top1/5/MODEL35`. Before anything can be analysed
consistently -- or cited, or deposited -- each MD needs an identity that does
not move, and the facts about it need to sit in one place with their evidence.

    MD_RESULTS/glossary.md         the glossary, one row per MD, for reading
    MD_RESULTS/glossary.parquet    the same, for scripts

The key is `md_key = <dataset_id>/<md_id>`. It is derived from the run's own
path against its dataset's DECLARED root (identity.py), so it survives a rebuild
and does not change when someone adds a run beside it.

HOW TO USE IT FOR ANALYSIS. Select the MDs you mean, then ask for their inputs:

    from core import glossary
    g = glossary.load("/path/to/<folder>")
    sel = glossary.select(g, system="gp120-occluded-open", state="ok")
    for _, md in sel.iterrows():
        top, trajs = glossary.input_paths(md)
        # e.g. MDAnalysis.Universe(top, *trajs) -- segments already in time order

or from the shell with glossary_query.py (list / show / paths / summary).

STATE. Every MD carries one of three states and the issues behind it:

    ok       every property this toolkit can check was checked and holds.
    review   usable, but something must be confirmed before trusting a result
             (segment order was guessed, the topology could not be checked,
             the dataset mixes different molecules, ...).
    blocked  the MDDB workflow will reject it as it stands (topology atom count
             differs from the trajectory, trajectory formats are mixed).

A state is a statement about evidence, not about science. `ok` does not mean
the simulation is correct, only that its bookkeeping is.
"""

import os

#: Issue -> (severity, one-line meaning). Severity decides the state.
ISSUES = {
    "topology_mismatch": ("blocked",
                          "topology declares a different atom count from the trajectory"),
    "mixed_formats": ("blocked",
                      "segments are in different trajectory formats; mwf requires one"),
    "no_topology": ("blocked", "no topology found for this MD"),
    "topology_unchecked": ("review",
                           "topology format declares no atom count, or the trajectory's is unknown"),
    "ordering_unverified": ("review",
                            "segment order could not be read from headers; natural filename order used"),
    # A NOTE, not a review item, and deliberately so. Whether a dataset mixes
    # molecules is a decision about the DATASET -- it decides what gets
    # deposited together -- not a defect of the MD. At `review` severity it
    # marked 270 of 273 MDs, which made the state column useless for choosing
    # analysis targets while saying nothing about any individual simulation.
    # The dataset-level view (overview table, `topology_groups`) is where it is
    # acted on; this note keeps it visible on every affected entry.
    "grouping_unreviewed": ("note",
                            "dataset holds more than one molecular system; confirm it belongs together"),
    "mixed_spacing": ("review",
                      "production segments were written at different frame spacings"),
    "no_config": ("review", "no MD config found; timestep, temperature and ensemble unknown"),
    "no_frame_count": ("review", "frame count not readable cheaply for this format"),
    "stage_assumed": ("review",
                      "a trajectory is production only by default (a .dcd whose name "
                      "says nothing, with no other production file beside it); "
                      "confirm it, or say what it is, under `folder: stages:` in settings.yaml"),
    "split_directory": ("note",
                        "one of several replicas stored side by side in one directory, "
                        "told apart by the replica label in its file name"),
    "short_of_nominal": ("review",
                         "shorter than every production length declared for its "
                         "dataset (production_lengths_ns): a run that stopped, a "
                         "test, or a stage -- not a full production MD"),
    "incomplete": ("note",
                   "production stopped short of the length its config asks for "
                   "(see completion); still usable, but not the planned run"),
    "segment_overlap": ("note",
                        "a restart segment re-writes step(s) the previous segment "
                        "already holds; simulated_ns counts each step once"),
    "partial_segment": ("review",
                        "a segment holds FEWER frames than its own header "
                        "claims: a copy still in flight, or a write that was "
                        "cut off. Counts here use what is on disk -- re-run "
                        "once the transfer has finished"),
    "segment_gap": ("review",
                    "steps are MISSING between segments: the production is not "
                    "one continuous trajectory. A segment was never copied, or "
                    "the run restarted from an earlier checkpoint. gap_ns says "
                    "how much is absent, and simulated_ns does not include it"),
    "inherited_config": ("note",
                         "config taken from a template directory, cross-checked against the header"),
    "chamber_topology": ("note",
                         "topology is a CHARMM-in-AMBER (chamber) prmtop; MDAnalysis cannot read it"),
    "relabelled_format": ("note",
                          "a segment's content contradicts its extension (e.g. NetCDF named .mdcrd); "
                          "its work link carries the true extension -- read via --links"),
}

#: Column order for mds.parquet and the CLI.
COLUMNS = [
    "md_key", "dataset_id", "md_id", "state", "issues",
    "system", "condition", "collections", "replica_index", "stage_dir",
    "engine", "program", "version", "natom", "topology_group",
    "n_segments", "frames", "simulated_ns", "gap_ns", "nominal_ns", "target_ns",
    "completion",
    "framestep_ns", "framestep_basis",
    "timestep_fs", "temp_K", "ensemble", "ff", "wat", "boxtype",
    "config_scope", "config_path",
    "ordering", "trajectory_gb", "topology", "topology_natom",
    "source_dir", "segments", "formats", "source_root",
    "work_dir", "links", "topology_link",
]


def _nan(v):
    return v is None or (isinstance(v, float) and v != v)


def _clean(v):
    """Registry values -> plain Python, with NaN as None."""
    if _nan(v):
        return None
    if hasattr(v, "item"):
        try:
            return v.item()
        except (AttributeError, ValueError):
            pass
    if isinstance(v, (list, tuple)) or (hasattr(v, "tolist") and not isinstance(v, str)):
        return ", ".join(str(x) for x in list(v))
    return v


def state_of(issues):
    """The worst severity among an MD's issues."""
    sev = {ISSUES[i][0] for i in issues if i in ISSUES}
    if "blocked" in sev:
        return "blocked"
    if "review" in sev:
        return "review"
    return "ok"


def segment_coverage(paths):
    """What the segments of one MD really cover: duplicates, and holes.

    Returns (duplicate_frames, missing_frames, blocks, partial_segments).
    `partial_segments` counts files holding fewer frames than their own header
    claims -- a transfer still in flight. `missing_frames` counts
    frame slots inside the span that NO segment holds -- CDR 2-JL427/3R on
    2026-09-22 held pr.part6 then pr.part10, 108.4 ns apart, because parts 7-9
    had not been copied from the cluster yet. A record that called that one 395 ns
    trajectory would be wrong in a way nobody could see from the outside.

    (0, 0, 1, 0) when the segments cannot be read or their step counters reset.
    """
    try:
        from core.dcdheader import read_dcd_header
    except ImportError:                                    # pragma: no cover
        return 0, 0, 1, 0
    heads = []
    for p in paths:
        if not str(p).lower().endswith(".dcd"):
            return 0, 0, 1, 0
        h = read_dcd_header(p)
        if "error" in h or not h.get("nsavc") or not h.get("n_frames"):
            return 0, 0, 1, 0
        heads.append(h)
    if len({h["nsavc"] for h in heads}) > 1:
        return 0, 0, 1, 0
    # Only a CONTINUING step counter says anything. Segments restarted with
    # firsttimestep 0 all begin at the same step (pr.dcd, pr2.dcd ... in the
    # HD4 trimers), which is a reset, not repeated or missing frames.
    if len({h["istart"] for h in heads}) == 1:
        return 0, 0, 1, 0
    total = sum(h["n_frames"] for h in heads)
    if total > 5_000_000:      # a set of that many steps is not worth the RAM
        return 0, 0, 1, 0
    nsavc = heads[0]["nsavc"]
    steps = set()
    for h in heads:
        steps.update(range(h["istart"],
                           h["istart"] + h["n_frames"] * nsavc, nsavc))
    ordered = sorted(steps)
    blocks = 1 + sum(1 for a, b in zip(ordered, ordered[1:]) if b - a > nsavc)
    span = (ordered[-1] - ordered[0]) // nsavc + 1
    return (total - len(steps), span - len(steps), blocks,
            sum(1 for h in heads if h.get("truncated")))


def overlap_frames(paths):
    """Frames at the joins of consecutive DCD segments that repeat a step.

    A restart from a checkpoint written before the previous segment's last
    frame writes that step again: 2-REPLICA/1R's pr.dcd ends at step
    264,500,000 and pr_r1.dcd (firsttimestep 264,499,998) begins there. Read
    from the headers (first step, save interval, frame count), so it is exact;
    0 when any segment is not a readable DCD or the intervals differ.

    Counted over every segment at once, as the set of step numbers they hold.
    Pairwise comparison in file order got CDR wrong on 2026-09-22: a window
    that died after one frame (`pr.part10.dcd`, 1 frame) shares its first step
    with the segment that resumed there, and a pair whose step counter does not
    ADVANCE used to abandon the whole MD's overlap accounting -- so YZ23/3R was
    reported 1.8 ns long and the apo 1R 1.4 ns long.
    """
    return segment_coverage(paths)[0]


#: The lab's production lengths, in ns, when datasets.yaml declares none.
#: Stated by the user on 2026-09-21: "our MDs could be 1 microsecond, 500
#: nanoseconds, 1.5 microseconds". Declare others (100 for a screening set, 250
#: for LC3B) under `production_lengths_ns:` at the top of datasets.yaml, or on
#: one dataset.
PRODUCTION_LENGTHS_NS = (500, 1000, 1500)

#: How far short of a declared length an MD may fall and still be that length:
#: 2 %, so 1008 ns is 1 us and 500.0 ns of a planned 504 is 500 ns.
NOMINAL_TOLERANCE = 0.02


def nominal_of(simulated_ns, lengths=PRODUCTION_LENGTHS_NS):
    """The largest declared length this MD reaches, or None.

    "Reaches" means at least that length, less the tolerance. An MD longer than
    it was declared (a 585 ns replica of a 500 ns campaign) is still that
    length's MD; one shorter than the smallest declared length is not a
    production MD of this dataset at all.
    """
    if simulated_ns is None or _nan(simulated_ns):
        return None
    hits = [float(L) for L in lengths or () if L and
            float(simulated_ns) >= float(L) * (1 - NOMINAL_TOLERANCE)]
    if not hits:
        return None
    best = max(hits)
    return int(best) if best == int(best) else best


def build_rows(dataset_id, spec, records, runs, configs, source_root,
               engines=None, lengths=PRODUCTION_LENGTHS_NS, groups=None,
               heterogeneous=False, work=None, declared_ps=None):
    """Glossary rows for one exported dataset.

    `records` is what mds.md_records produced: one per MD, with its files in
    time order, its topology and whether a directory was split into replicas.
    `groups` maps a topology path to its content group, `work` (optional) is
    the MDDB working area, added by the deposition pipeline when it runs.

    `declared_ps` maps a topology group to the frame spacing a person wrote in
    settings.yaml (`ndp.declared_intervals`). It is used ONLY where the headers
    state none -- a declaration never replaces a measurement, for the same
    reason a curated field never overwrites a measured one -- and every row it
    fills says so in `framestep_basis`.
    """
    from .mds import topology_atoms
    try:
        from core.trajscan import frame_info, sniff_format
    except ImportError:                                  # pragma: no cover
        frame_info = sniff_format = None
    project = work

    run_by = {str(r.rel_path): r for _, r in runs.iterrows()}
    cfg_by = {}
    if configs is not None and len(configs):
        cfg_by = {str(r.rel_path): r for _, r in configs.iterrows()}
    groups = groups or {}
    declared_ps = declared_ps or {}
    heterogeneous = bool(heterogeneous) \
        and str(spec.get("heterogeneous", "")).lower() != "reviewed"

    rows = []
    for rec in records:
        rel = rec["source_dir"]
        r = run_by.get(rel)
        c = cfg_by.get(rel)
        issues = []

        natom = None if r is None else _clean(r.natom)
        top = rec.get("topology")
        top_natom = topology_atoms(os.path.join(source_root, top)) if top else None
        if not top:
            issues.append("no_topology")
        elif top_natom is None or natom is None:
            issues.append("topology_unchecked")
        elif int(top_natom) != int(natom):
            issues.append("topology_mismatch")
        if top and str(top).lower().endswith((".prmtop", ".parm7")):
            try:
                from core.trajscan import is_chamber_prmtop
                if is_chamber_prmtop(os.path.join(source_root, top)):
                    issues.append("chamber_topology")
            except ImportError:
                pass

        segs = rec.get("trajectories") or []
        # Formats by CONTENT where it can be read, by extension otherwise. The
        # extension alone would call `pr.dcd` + `MD1.mdcrd` mixed and miss a
        # NetCDF file named .mdcrd entirely.
        formats = []
        for s in segs:
            f = sniff_format(os.path.join(source_root, rel, s)) if sniff_format else None
            formats.append(f or os.path.splitext(s)[1].lower().lstrip("."))
        if len(set(formats)) > 1:
            issues.append("mixed_formats")
        if rec.get("relabelled"):
            issues.append("relabelled_format")
        if rec.get("ordering") == "natural":
            issues.append("ordering_unverified")
        if heterogeneous:
            issues.append("grouping_unreviewed")
        if rec.get("stage_assumed"):
            issues.append("stage_assumed")
        replica = rec.get("replica")
        if replica is not None:
            issues.append("split_directory")

        dt_mixed = bool(r.dt_mixed) if r is not None and not _nan(r.dt_mixed) else False
        if dt_mixed:
            issues.append("mixed_spacing")
        frames = None if r is None else _clean(r.n_frames_production)
        dt_ps = None if r is None else _clean(r.dt_ps_production)
        simulated_ns = None if r is None else _clean(r.prod_ns_header)

        # Where the spacing came from, before anything is derived from it.
        # A declaration fills only the silence: if the headers stated a
        # spacing it stands, and if they CONTRADICTED each other
        # (`dt_mixed`) the contradiction stands too -- one number written in
        # settings.yaml is a person supplying what was lost, not a licence to
        # paper over segments that disagree.
        framestep_basis = None if dt_ps is None else "headers"
        if dt_ps is None and not dt_mixed:
            declared = declared_ps.get(groups.get(top))
            if declared:
                dt_ps, framestep_basis = float(declared), "declared"

        if replica is not None:
            # The index counts per DIRECTORY; this MD is one file of several.
            # Its own frames come from its own header, and its length from
            # those frames at the directory's spacing -- one spacing, since
            # dt_mixed would have withheld it.
            frames, simulated_ns = None, None
            if frame_info is not None:
                counts = [frame_info(os.path.join(source_root, rel, s))[0]
                          for s in segs]
                if counts and all(c is not None for c in counts):
                    frames = int(sum(counts))
                    if dt_ps and not dt_mixed:
                        simulated_ns = round(frames * dt_ps / 1000.0, 2) or None
        if frames is None:
            issues.append("no_frame_count")
        # A frame count and a declared spacing are a length. Only ever filled
        # in, never recomputed over a header's own figure.
        if simulated_ns is None and frames and framestep_basis == "declared":
            simulated_ns = round(frames * dt_ps / 1000.0, 2) or None
        dup, missing, _blocks, partial = segment_coverage(
            [os.path.join(source_root, rel, s) for s in segs]) \
            if len(segs) > 1 else (0, 0, 1, 0)
        if dup:
            issues.append("segment_overlap")
            if simulated_ns and dt_ps and not dt_mixed:
                simulated_ns = round(simulated_ns - dup * dt_ps / 1000.0, 2)
        if partial:
            issues.append("partial_segment")
        gap_ns = None
        if missing:
            issues.append("segment_gap")
            if dt_ps and not dt_mixed:
                gap_ns = round(missing * dt_ps / 1000.0, 2)

        # How long this MD was MEANT to be, from its production configs, and
        # how much of that the trajectories hold. A 1 us campaign whose replica
        # holds 70 ns is a run that stopped, not a 70 ns simulation.
        from .protocol import planned_production_ns
        target_ns = planned_production_ns(os.path.join(source_root, rel))
        if target_ns and simulated_ns and simulated_ns > target_ns * 1.01:
            # Holding MORE than the config asks for means that config describes
            # one restart segment, not the whole run (HD4 trimers restarted from
            # step 0: 585 ns under a 504 ns config). Not a target, then.
            target_ns = None
        completion = (round(simulated_ns / target_ns, 3)
                      if target_ns and simulated_ns else None)
        # One frame of tolerance: the last frame lands on the final step.
        if completion is not None and dt_ps and \
                simulated_ns + dt_ps / 1000.0 < target_ns * 0.999:
            issues.append("incomplete")
        nominal_ns = nominal_of(simulated_ns, lengths)
        if simulated_ns is not None and nominal_ns is None and lengths:
            issues.append("short_of_nominal")

        cfg_path = None if c is None else _clean(c.config_path)
        scope = None if c is None else _clean(c.config_scope)
        if not cfg_path:
            issues.append("no_config")
        elif scope == "inherited":
            issues.append("inherited_config")

        size_gb = 0.0
        for s in segs:
            try:
                size_gb += os.path.getsize(os.path.join(source_root, rel, s)) / 1e9
            except OSError:
                pass

        rows.append({
            "md_key": f"{dataset_id}/{rec['md_id']}",
            "dataset_id": dataset_id,
            "md_id": rec["md_id"],
            "state": state_of(issues),
            "issues": ", ".join(issues),
            "system": None if r is None else _clean(r.system),
            "condition": None if r is None else _clean(r.condition),
            "collections": _clean(spec.get("collections")),
            "replica_index": f"{replica}R" if replica is not None
            else (None if r is None or "replica_index" not in r
                  else _clean(r.replica_index)),
            "stage_dir": None if r is None or "stage_dir" not in r
            else _clean(r.stage_dir),
            # `engine` is the trajectory FORMAT FAMILY the index recorded
            # (`NAMD/CHARMM` names two programs); `program` and `version` are
            # the software that actually ran, from its own start-up banner.
            "engine": None if r is None else _clean(r.engine),
            "program": (engines or {}).get(rel, {}).get("program"),
            "version": (engines or {}).get(rel, {}).get("version"),
            "natom": natom,
            "topology_group": groups.get(top),
            "n_segments": len(segs),
            "frames": frames,
            "simulated_ns": simulated_ns,
            "gap_ns": gap_ns,
            "nominal_ns": nominal_ns,
            "target_ns": target_ns,
            "completion": completion,
            "framestep_ns": None if (dt_ps is None or dt_mixed) else round(dt_ps / 1000.0, 6),
            "framestep_basis": None if (dt_ps is None or dt_mixed) else framestep_basis,
            "timestep_fs": None if c is None else _clean(c.timestep),
            "temp_K": None if c is None else _clean(c.temp),
            "ensemble": None if c is None else _clean(c.ensemble),
            "ff": None if c is None else _clean(c.ff),
            "wat": None if c is None else _clean(c.wat),
            "boxtype": None if c is None else _clean(c.boxtype),
            "config_scope": scope,
            "config_path": cfg_path,
            "ordering": rec.get("ordering"),
            "trajectory_gb": round(size_gb, 3),
            "topology": top,
            "topology_natom": top_natom,
            "source_dir": rel,
            "segments": "; ".join(segs),
            "formats": "; ".join(formats),
            "source_root": os.path.realpath(source_root),
            # Where the deposition workflow runs for this MD, when it has been
            # prepared: link-only directories, added by the MDDB pipeline.
            "work_dir": project,
            "links": "; ".join(os.path.join(rec["md_id"], name)
                               for name in rec.get("links") or []),
            "topology_link": rec.get("topology_link"),
        })
    return rows


#: Columns coerced to numbers on the way out, so a value that arrived as a
#: string from a registry does not make the whole column object-typed.
NUMERIC = ("natom", "frames", "simulated_ns", "gap_ns", "nominal_ns", "target_ns",
           "completion",
           "framestep_ns", "timestep_fs",
           "temp_K", "trajectory_gb", "topology_natom")


def _atomic(path, writer):
    """Write to a temporary name beside `path`, then rename into place.

    A reader -- the nightly combine, a colleague's analysis script -- never
    sees a half-written file, and a build that dies after this point leaves
    the previous file intact instead of a truncated one.
    """
    from .atomic import atomic                               # noqa: PLC0415
    atomic(path, writer)


def write(mddb, rows, columns=None, page=None, **render_kw):
    """Write the glossary: registry/mds.parquet under `mddb`, and its page.

    `rows` is a list of dicts or a DataFrame. `columns` defaults to COLUMNS;
    the lab-wide glossary passes its own list. `page` is where the Markdown
    goes -- MD_RESULTS/glossary.md for a data folder, and `<mddb>/glossary.md`
    for the lab-wide copy, which is a directory of its own. `render_kw` is
    handed to render_markdown (group=, title=, intro=). Both are written
    atomically.
    """
    import pandas as pd
    columns = list(columns or COLUMNS)
    if isinstance(rows, pd.DataFrame):
        df = rows.reindex(columns=columns)
    else:
        df = pd.DataFrame(rows, columns=columns)
    for col in NUMERIC:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    reg = os.path.join(mddb, "registry")
    out_page = page or os.path.join(mddb, "glossary.md")
    os.makedirs(reg, exist_ok=True)
    os.makedirs(os.path.dirname(out_page), exist_ok=True)
    _atomic(os.path.join(reg, "mds.parquet"),
            lambda tmp: df.to_parquet(tmp, index=False))
    text = render_markdown(df, **render_kw)

    def _write_text(tmp):
        with open(tmp, "w") as fh:
            fh.write(text)
    from core import layout                                  # noqa: PLC0415
    if not layout.same_page(out_page, text):
        _atomic(out_page, _write_text)
    return df


def _fmt(v, nd=None):
    if _nan(v):
        return "—"
    if nd is not None:
        return f"{float(v):,.{nd}f}"
    if isinstance(v, float) and v.is_integer():
        return f"{int(v):,}"
    return str(v)


def render_markdown(df, group="dataset_id", title="MD glossary", intro=None):
    """The glossary as Markdown.

    `group` is the column, or list of columns, that makes one section:
    `dataset_id` for one folder's glossary; `["member", "dataset_id"]` for the
    lab-wide one, where two members may both have a dataset called
    `gp120-open` and must not be merged. `intro` is extra Markdown (a string
    or a list of lines) placed before the overview; the lab glossary puts its
    table of folders there.
    """
    import re
    import time
    keys = [group] if isinstance(group, str) else list(group)
    by = keys[0] if len(keys) == 1 else keys

    def _label(v):
        v = v if isinstance(v, tuple) else (v,)
        return "/".join(str(x) for x in v)

    def _slug(s):
        return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")

    L = [f"# {title}", "",
         f"Generated {time.strftime('%Y-%m-%dT%H:%M:%S')}. One entry per MD, "
         f"{len(df)} in total. Machine-readable copy: "
         "`glossary.parquet` beside this report. Query it with "
         "`glossary_query.py`.", "",
         "## Terms", "",
         "| term | meaning |", "|---|---|",
         "| **MD** | one continuous simulation of one system: a topology plus its "
         "production trajectory segments, in time order |",
         "| **md_key** | `<dataset_id>/<md_id>`, the stable name of an MD. Derived "
         "from the run's path against its dataset's declared root, so it does not "
         "change on rebuild or when runs are added |",
         "| **dataset** | a group of MDs declared in `datasets.yaml`; one MDDB "
         "project |",
         "| **segment** | one trajectory file of an MD; restarts produce several |",
         "| **ordering** | how segment order was established: `header` (first-frame "
         "step read from the file), `single` (one segment), `natural` (filename "
         "order, a guess) |",
         "| **replica_index** | the repeat number read from the path (`1R`, "
         "`rep2`, `3-R3`), or empty |",
         "| **stage_dir** | a directory name that denotes where a run lives (`pr`, "
         "`rep`, `md`), not which repeat it is |",
         "| **topology_group** | first 8 characters of the topology file's content "
         "hash. MDs sharing a group use the same topology file content; different "
         "groups are different molecular systems or builds |",
         "| **simulated_ns** | production time, summed per segment from each "
         "file's own header; where the headers carry no timestep, the frame "
         "count at the spacing declared in settings.yaml -- see "
         "**framestep_basis** |",
         "| **nominal_ns** | the production length this MD reaches, among those "
         "declared for its dataset (`production_lengths_ns`, default 500, 1000 "
         "and 1500 ns; 2 % tolerance): 1008 ns is 1000. Empty when the MD is "
         "shorter than all of them, which is the issue `short_of_nominal` |",
         "| **target_ns** | how long production was meant to run, from the "
         "production configs (first step to the furthest `run` any of them asks "
         "for); empty for engines or layouts it cannot read |",
         "| **completion** | simulated_ns / target_ns. Below 1 (by more than one "
         "frame) the MD carries the note `incomplete` |",
         "| **framestep_ns** | time between saved frames; empty when segments "
         "disagree |",
         "| **framestep_basis** | where **framestep_ns** came from: `headers` "
         "if the trajectories stated it, `declared` if a person wrote "
         "`frame_interval_ps` under that system in settings.yaml because a "
         "stripping or concatenation tool had thrown the original away. Anything "
         "derived from a `declared` spacing -- simulated_ns, and the record's "
         "lengths and total -- rests on that declaration, not on a measurement |",
         "| **work_dir** | where the MDDB workflow runs for this MD; holds only "
         "read-only links to the source |", "",
         "## States", "",
         "| state | meaning |", "|---|---|",
         "| `ok` | every property this toolkit can check was checked and holds. "
         "Says nothing about whether the science is right |",
         "| `review` | usable, but confirm the listed issue before trusting a "
         "result |",
         "| `blocked` | the MDDB workflow will reject this MD as it stands |", "",
         "| issue | state | meaning |", "|---|---|---|"]
    for k, (sev, text) in ISSUES.items():
        L.append(f"| `{k}` | {sev} | {text} |")
    if intro:
        L += [""] + ([intro] if isinstance(intro, str) else list(intro))
    L += ["", "## Overview", "",
          "`topology groups` counts distinct topology files by content. More than "
          "one means the dataset holds different molecular systems or builds, "
          "which is a grouping decision to make in `datasets.yaml` before "
          "deposition: split it, or mark it `heterogeneous: reviewed`.", "",
          f"| {' / '.join(keys)} | MDs | ok | review | blocked | topology groups "
          "| grouping | simulated | trajectory |",
          "|---|---:|---:|---:|---:|---:|---|---:|---:|"]
    for lab, g in df.groupby(by, sort=True):
        did = _label(lab)
        ns = g.simulated_ns.sum(skipna=True)
        n_groups = g.topology_group.nunique(dropna=True)
        flagged = g.issues.fillna("").str.contains("grouping_unreviewed").any()
        grouping = ("**unreviewed**" if flagged
                    else ("reviewed" if n_groups > 1 else "single system"))
        L.append(f"| [{did}](#{_slug(did)}) | {len(g)} | "
                 f"{int((g.state == 'ok').sum())} | "
                 f"{int((g.state == 'review').sum())} | "
                 f"{int((g.state == 'blocked').sum())} | {n_groups} | {grouping} | "
                 f"{ns / 1000:.2f} µs | {g.trajectory_gb.sum():,.1f} GB |")
    counts = {}
    for s in df.issues.fillna(""):
        for i in [x.strip() for x in s.split(",") if x.strip()]:
            counts[i] = counts.get(i, 0) + 1
    if counts:
        L += ["", "| issue | MDs |", "|---|---:|"]
        for k, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            L.append(f"| `{k}` | {n} |")
    for lab, g in df.groupby(by, sort=True):
        did = _label(lab)
        L += ["", f"## {did}", ""]
        n_groups = g.topology_group.nunique(dropna=True)
        L.append(f"{len(g)} MDs across {n_groups} topology group(s). "
                 f"System `{', '.join(sorted({str(x) for x in g.system.dropna()}))}`.")
        L += ["", "| md_id | state | atoms | topo group | segments | ns | "
              "framestep | T / ensemble | issues | source |",
              "|---|---|---:|---|---:|---:|---:|---|---|---|"]
        for _, r in g.sort_values("md_id").iterrows():
            te = " / ".join(x for x in (
                None if _nan(r.temp_K) else f"{int(r.temp_K)} K",
                None if _nan(r.ensemble) else str(r.ensemble)) if x) or "—"
            L.append(
                f"| `{r.md_id}` | {r.state} | {_fmt(r.natom)} | "
                f"{r.topology_group or '—'} | {r.n_segments} | "
                f"{_fmt(r.simulated_ns, 1)} | {_fmt(r.framestep_ns)} | {te} | "
                f"{r.issues or '—'} | `{r.source_dir}` |")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------
# reading it back
# --------------------------------------------------------------------------
def glossary_path(folder):
    """Where the glossary of `folder` is, or None.

    Accepts a data folder in the current layout
    (`.md_pipeline/mddb/registry/mds.parquet`) or the one before it
    (`MDDB_substrate/registry/...`, still on colleagues' folders), the lab-wide
    copy (a bare `registry/mds.parquet`), or a `.parquet` file directly.
    """
    from core import layout
    folder = os.path.abspath(str(folder).rstrip("/"))
    if folder.endswith(".parquet") and os.path.isfile(folder):
        return folder
    return layout.glossary_parquet(folder)


def load(folder):
    """The glossary for a data folder (or the lab-wide one), as a DataFrame."""
    import pandas as pd
    path = glossary_path(folder)
    if path is None:
        raise FileNotFoundError(
            f"no glossary under {folder} -- run mddb_build.py on it first")
    return pd.read_parquet(path)


def select(df, dataset=None, system=None, condition=None, state=None,
           issue=None, where=None, keys=None):
    """Filter the glossary. Every argument is optional; lists mean "any of".

    `where` is a pandas query string for anything else, e.g.
    "temp_K == 310 and simulated_ns > 400".
    """
    def _any(col, val):
        vals = [val] if isinstance(val, str) else list(val)
        return df[col].astype(str).isin([str(v) for v in vals])

    mask = df.md_key.notna()
    if keys:
        hit = _any("md_key", keys)
        if "lab_key" in df.columns:          # the lab-wide glossary
            hit |= _any("lab_key", keys)
        mask &= hit
    if dataset:
        mask &= _any("dataset_id", dataset)
    if system:
        mask &= _any("system", system)
    if condition:
        mask &= _any("condition", condition)
    if state:
        mask &= _any("state", state)
    if issue:
        wanted = [issue] if isinstance(issue, str) else list(issue)
        mask &= df.issues.fillna("").apply(
            lambda s: any(w in [x.strip() for x in s.split(",")] for w in wanted))
    out = df[mask]
    if where:
        out = out.query(where)
    return out


def input_paths(md, via="source"):
    """(topology, [trajectory segments in time order]) as absolute paths.

    via="source"  the original files, under their original names.
    via="links"   the read-only links under work_dir. Same bytes, but a file
                  whose content contradicts its extension (a NetCDF named
                  .mdcrd) is reached through a link with the TRUE extension,
                  so tools that pick a reader from the name open it correctly.
                  Prefer this for analysis unless you need the source names.

    Either way these are for reading: never write through them. Pass a glossary
    row (a pandas Series or a dict).
    """
    get = md.get if hasattr(md, "get") else (lambda k: md[k])
    if via == "links":
        work = get("work_dir")
        links = [s for s in str(get("links") or "").split("; ") if s]
        tl = get("topology_link")
        return (os.path.join(work, tl) if tl else None,
                [os.path.join(work, s) for s in links])
    root = get("source_root")
    rel = get("source_dir")
    top = get("topology")
    segs = [s for s in str(get("segments") or "").split("; ") if s]
    return (os.path.join(root, top) if top else None,
            [os.path.join(root, rel, s) for s in segs])
