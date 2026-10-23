#!/usr/bin/env python3
"""
ndp.py -- one NDP metadata record (JSON) per simulated system, from the glossary.

WHAT NDP WANTS. The National Data Platform catalogue is metadata-only: an entry
is a title, a description, keywords and links, and the lab's registration guide
fills the link with one small JSON file per simulated system, always with the
same field names (title, lab, contact, system, setup, sampling, ...).

WHAT A "SIMULATED SYSTEM" IS HERE. The glossary has one row per MD. The guide
wants one entry per system, with its replicas counted inside it ("3 x 500 ns").
MDs that share a `topology_group` were run from byte-identical topologies, so
they are replicas of one system; that group, inside its dataset, is the unit:

    record_id = <dataset_id>/<topology_group>      e.g. gp120-open/d2830a2a

It is as stable as the md_key it sits beside: it changes only if the topology
file itself does.

WHERE EACH VALUE COMES FROM. Three sources, never mixed up:

    measured   the glossary (atoms, force field, program, ensemble, temperature,
               timestep, replicas, lengths, frame spacing) and the MD config
               (pressure). A value the MDs of one system DISAGREE on is left
               null and reported; it is never settled by majority.
    declared   settings.yaml `datasets:` (contact, licence, citation, pdb_ids, wat).
    curated    settings.yaml `ndp:` -- what only a person knows: the title, the construct,
               the starting structure, the glycan profile, the ions.

Every record carries every field of the template. A value nobody established is
null and makes the record `incomplete`; it is never guessed. (The mwf inputs
OMIT such keys instead, because mwf reads a null as a statement. The NDP
template asks for the same field names in every file, so here the key stays.)
"""

import hashlib
import json
import os
import re

from . import lab

#: Human-only fields of one system, in the order ndp.yaml lists them.
CURATED = ["identifier", "topic", "tags", "title", "description", "construct",
           "starting_structure", "glycan_profile", "ions",
           "related_publication", "keywords", "slug"]

#: What the lab catalogue (lab.yaml `catalog_url`) shows for every dataset:
#: a stable identifier, the research
#: topic it is filed under, free tags, a title, a description, the volume and
#: the file count. The first three are a person's to decide, so the record is
#: `incomplete` until they are filled; the rest are measured.
WEBSITE = ["identifier", "topic", "tags", "title", "description",
           "size_gb", "files"]

#: Measured fields. ndp.yaml asks for one ONLY where the measurement came back
#: empty (a run with no config beside it), and a curated value never replaces a
#: measured one: a wrong measurement is repaired where it was made.
FALLBACK = ["force_field", "water_model", "software", "ensemble",
            "temperature_K", "pressure_bar", "timestep_fs", "frame_interval_ps"]

#: Who may get the raw files. `licensable`: on terms, through a licence.
ACCESS_CHOICES = ("public", "private", "licensable")

#: Whether the data may be used to train AI/ML models -- a permission of its
#: own, because the citation licence (CC-BY) says nothing about it.
AI_CHOICES = ("allowed", "prohibited", "by-agreement")

#: How people write those answers. YAML reads a bare yes/no as true/false.
AI_ALIASES = {"true": "allowed", "yes": "allowed", "permitted": "allowed",
              "false": "prohibited", "no": "prohibited", "not-allowed": "prohibited",
              "not allowed": "prohibited", "forbidden": "prohibited",
              "on-request": "by-agreement", "on request": "by-agreement",
              "by agreement": "by-agreement", "ask": "by-agreement",
              "licensable": "by-agreement"}

#: Fields that are part of a record but whose absence never makes it
#: `incomplete`: an ORCID nobody has, a download URL for data that is not
#: online.
OPTIONAL = {"data_owner.email", "data_owner.orcid", "data_owner.affiliation",
            "data_location.uri"}

#: Defaults written into a new settings.yaml's `ndp: defaults:`. The lab's own
#: (name, keywords, data owner, access, AI training, storage) come from
#: lab.yaml; nothing here names a lab. `base_url` is where records/ will be
#: hosted (a GitHub raw URL or a Zenodo record); it feeds the item URL on the
#: registration sheet and the identifiers of the standard formats.
_LAB = lab.record_defaults()
DEFAULTS = {
    "lab": lab.name(),
    "data_availability": _LAB.get("data_availability"),
    "keywords": list(_LAB.get("keywords") or []),
    "data_owner": _LAB.get("data_owner"),
    "access_level": _LAB.get("access_level"),
    "ai_training": _LAB.get("ai_training"),
    "data_location": _LAB.get("data_location"),
    "base_url": None,
    # Which stages of the preparation protocol the record lists:
    #   production  production stages only (pr, pr_r1, ...) -- the default
    #   all         minimisation, heating, equilibration and production too
    #   none        no protocol section content
    # Production only, because that is what is published: "we want only the
    # information regarding the production phase" (2026-09-22). `all` puts the
    # preparation history back, for a reviewer who asks how it was built.
    # Sampling (replicas, lengths, total) is production-only whatever this says.
    "protocol": "production",
    # Which MDs count as the record's replicas:
    #   all      every MD of the system, whatever its length, each one marked
    #            complete or not in `sampling.replicas_detail`. The catalogue
    #            lists what exists, so this is the default.
    #   nominal  only those that reach a declared production length (the
    #            glossary's nominal_ns). Shorter ones are listed as left out.
    # An MD whose length could not be measured is kept either way, and flagged.
    "lengths": "all",
}



def merged_defaults(section):
    """DEFAULTS under a folder's own `ndp: defaults:`. A key the folder lists
    with no value (`access_level:`) is unanswered there, so the lab's value
    stands; `keywords: []` is an answer, and empties the list."""
    def unanswered(v):
        return v is None or (isinstance(v, dict) and all(_blank(x) for x in v.values()))
    return {**DEFAULTS, **{k: v for k, v in (section or {}).items() if not unanswered(v)}}


#: Values `lengths:` may take in ndp.yaml.
LENGTH_CHOICES = ("nominal", "all")

#: Values `protocol:` may take in ndp.yaml.
PROTOCOL_CHOICES = ("all", "production", "none")

#: datasets.yaml spells licences the way MDDB does; the template wants SPDX.
SPDX = {
    "creative commons attribution 4.0 international": "CC-BY-4.0",
    "creative commons attribution-sharealike 4.0 international": "CC-BY-SA-4.0",
    "creative commons attribution-noncommercial 4.0 international": "CC-BY-NC-4.0",
    "creative commons zero v1.0 universal": "CC0-1.0",
}

#: NDP's keyword box accepts nothing else (registration guide, step 3).
KEYWORD_OK = re.compile(r"^[A-Za-z0-9 _.\-]+$")

#: Replica lengths within this fraction of one another read as "n x L ns".
UNIFORM = 0.02


def _nan(v):
    return v is None or (isinstance(v, float) and v != v)


def _blank(v):
    return _nan(v) or (isinstance(v, str) and not v.strip()) or v == []


def _num(v):
    """A float that is a whole number becomes an int, so 310.0 K reads 310."""
    if _nan(v):
        return None
    v = float(v)
    return int(v) if v == int(v) else round(v, 5)


def _float(v):
    return None if _nan(v) else round(float(v), 5)


def _abs_dir(root, rel):
    """`root`/`rel` as an absolute, normalised path; None without a root."""
    if _blank(root) or _blank(rel):
        return None
    return os.path.normpath(os.path.join(str(root), str(rel)))


def raw_location(detail):
    """The directory every replica's files sit under, or None."""
    dirs = [d["location"] for d in detail or [] if d.get("location")]
    if not dirs:
        return None
    try:
        common = os.path.commonpath(dirs)
    except ValueError:                 # relative and absolute mixed
        return None
    # Replicas on different filesystems share only `/`: that is not where
    # anything lives, and stating it would read as an answer.
    if not [p for p in common.split(os.sep) if p]:
        return None
    return common


def _distinct(values):
    """Distinct established values, in first-seen order. None is "this MD does
    not say", not a competing value."""
    out = []
    for v in values:
        if _blank(v):
            continue
        if hasattr(v, "item"):
            v = v.item()
        if v not in out:
            out.append(v)
    return out


# --------------------------------------------------------------------------
# grouping
# --------------------------------------------------------------------------
def systems_of(glossary_df, include_blocked=False):
    """[(record_id, rows, n_blocked)] -- one per dataset and topology group.

    Blocked MDs are left out by default: their topology and trajectory
    demonstrably do not fit together, so counting them as sampling of this
    system would overstate the record. A system whose MDs are ALL blocked
    yields no record at all, and the caller reports it.
    """
    out = []
    df = glossary_df.sort_values(["dataset_id", "md_id"])
    for did, ds in df.groupby("dataset_id", sort=True):
        group = ds.topology_group.where(ds.topology_group.notna(),
                                        "md-" + ds.md_id.astype(str))
        for gid, rows in ds.groupby(group, sort=True):
            blocked = rows[rows.state == "blocked"]
            kept = rows if include_blocked else rows[rows.state != "blocked"]
            out.append((f"{did}/{gid}", kept, len(blocked)))
    return out


# --------------------------------------------------------------------------
# measured
# --------------------------------------------------------------------------
def pressure_bar(source_root, config_path, program):
    """Target pressure as the MD config states it, in bar, or None.

    Only a stated value counts: AMBER's undeclared default of 1 bar is a fact
    about AMBER, not something this config says.
    """
    if _blank(source_root) or _blank(config_path):
        return None
    from . import mdconfig
    try:
        with open(os.path.join(source_root, config_path), errors="replace") as fh:
            raw = fh.read()
    except OSError:
        return None
    try:
        if str(program or "").upper() == "NAMD":
            text, _ = mdconfig.resolve_tcl(mdconfig.strip_disabled_blocks(raw))
            val = mdconfig._directive(text, "langevinPistonTarget")
        else:
            val = mdconfig._amber_nml(raw, "pres0")
        return round(float(val), 5) if val is not None else None
    except ValueError:
        return None


def measure(rows, read_pressure=pressure_bar):
    """({template field: value}, [warnings]) for the MDs of one system."""
    warn = []

    def one(column, label, cast=lambda v: v, values=None):
        vals = _distinct(rows[column] if values is None else values)
        if len(vals) > 1:
            warn.append(f"`{label}`: the MDs disagree ({', '.join(map(str, vals))}); "
                        f"left null. Split the system or correct the MD that is wrong.")
            return None
        return cast(vals[0]) if vals else None

    if "issues" in rows.columns:
        assumed = [str(k) for k, i in zip(rows.md_key, rows.issues.fillna(""))
                   if "stage_assumed" in str(i)]
        if assumed:
            warn.append(f"{len(assumed)} MD(s) are production only by default "
                        f"(stage_assumed): confirm the files are the simulation, "
                        f"under `folder: stages:` in settings.yaml, before registering.")

    ensemble = one("ensemble", "ensemble")
    pressure = None
    if ensemble == "NPT":
        pressure = one(None, "pressure_bar", _float, values=[
            read_pressure(r.source_root, r.config_path, r.program)
            for r in rows.drop_duplicates("config_path").itertuples()])

    software = []
    for r in rows.itertuples():
        if _blank(r.program):
            continue
        name = str(r.program) if _blank(r.version) else f"{r.program} {r.version}"
        if name not in software:
            software.append(name)
    # "NAMD" alone adds nothing beside "NAMD 3.0.2" of the same system.
    software = [s for s in software
                if not any(o != s and o.startswith(s + " ") for o in software)]

    nominal = one("nominal_ns", "nominal_length_ns", _num) \
        if "nominal_ns" in rows.columns else None

    spacing = one(None, "frame_interval_ps", _num, values=[
        None if _nan(v) else round(float(v) * 1000.0, 3) for v in rows.framestep_ns])

    detail, complete = [], 0
    for r in rows.itertuples():
        ns = None if _nan(r.simulated_ns) else round(float(r.simulated_ns), 1)
        target = None if _nan(getattr(r, "target_ns", None)) else \
            round(float(r.target_ns), 1)
        # "Complete" means: it reached the length its own config asked for, or,
        # where no config says, one of the lab's production lengths.
        done = (None if ns is None else
                (ns + 0.5 >= target if target else
                 bool(not _nan(getattr(r, "nominal_ns", None)))))
        entry = {"id": str(r.md_id), "length_ns": ns,
                 "planned_ns": target, "complete": done,
                 "frames": None if _nan(r.frames) else int(r.frames),
                 "path": str(r.source_dir)}
        # Where this replica's files ARE, not just where they sit in the
        # folder: what the federation needs to point at the raw data.
        where = _abs_dir(getattr(r, "source_root", None), r.source_dir)
        if where:
            entry["location"] = where
        # Time that is missing INSIDE the run: segments that never arrived, or
        # a restart from an earlier checkpoint. `length_ns` is what is there,
        # so without this a reader would take it for one continuous run.
        gap = getattr(r, "gap_ns", None)
        if gap is not None and not _nan(gap) and float(gap) > 0:
            entry["gap_ns"] = round(float(gap), 1)
            entry["continuous"] = False
        detail.append(entry)
        complete += bool(done)

    lengths = [round(float(v), 1) for v in rows.simulated_ns if not _nan(v)]
    if len(lengths) < len(rows):
        warn.append(f"{len(rows) - len(lengths)} MD(s) have no readable length and "
                    f"are counted as replicas but not as sampling.")
    if lengths and min(lengths) < 0.5 * max(lengths):
        warn.append(f"replica lengths are very uneven ({min(lengths)} to "
                    f"{max(lengths)} ns): check these are all replicas of one "
                    f"simulation and not a test run or a restart stage.")

    return {
        "total_atoms": one("natom", "total_atoms", _num),
        "force_field": one("ff", "force_field", str),
        "water_model": one("wat", "water_model", str),
        "software": "; ".join(software) or None,
        "ensemble": ensemble,
        "temperature_K": one("temp_K", "temperature_K", _num),
        "pressure_bar": pressure,
        "timestep_fs": one("timestep_fs", "timestep_fs", _float),
        "replicas": int(len(rows)),
        "complete_replicas": complete,
        "replicas_detail": detail,
        "planned_length_ns": one("target_ns", "planned_length_ns", _num)
        if "target_ns" in rows.columns else None,
        "size_gb": round(float(sum(v for v in rows.trajectory_gb if not _nan(v))), 2)
        if "trajectory_gb" in rows.columns else None,
        "files": int(sum(len([x for x in str(s).split("; ") if x])
                         for s in rows.segments)) if "segments" in rows.columns else None,
        "nominal_ns": nominal,
        "length_per_replica_ns": round(sum(lengths) / len(lengths), 1) if lengths else None,
        "replica_lengths_ns": lengths,
        "total_ns": round(sum(lengths), 1) if lengths else None,
        "frame_interval_ps": spacing,
        # WHERE THE LENGTHS REST. The glossary says whether its frame spacing
        # was read from the trajectory headers or supplied by a person in
        # settings.yaml, and a length computed from a declared spacing is a
        # declared length. Carrying that through kept the record honest once
        # the declaration reached the glossary: before it did,
        # `lengths_from_spacing` filled the lengths here and marked them
        # declared; afterwards the lengths arrived already filled and the
        # mark was silently lost, so the record stated 5 x 1050 ns with
        # nothing to say the 100 ps behind it was somebody's word.
        "length_basis": ("declared frame interval"
                         if "framestep_basis" in rows.columns
                         and "declared" in _distinct(rows.framestep_basis)
                         else None),
    }, warn


# --------------------------------------------------------------------------
# curated: ndp.yaml
# --------------------------------------------------------------------------
SECTION_HEAD = """
# --- the catalogue record ----------------------------------------------------
# What only you know about each simulated system. One entry per system: the MDs
# of one dataset that share a topology are its replicas. Atoms, force field,
# program, ensemble, temperature, pressure, timestep, replicas, lengths and
# frame spacing are MEASURED and never asked for here; contact, license and
# citation come from the dataset above.
#
# This section is yours. Re-running never rewrites it: entries for systems that
# have appeared since are appended at the end, and nothing else is touched.
# A field left empty stays null in the record and keeps it `incomplete`.
# Write `none` where a field does not apply (a glycan profile for LC3B).
#
#   title               names the construct and the key setup. No length or
#                       replica count: those are measured and added for you.
#   starting_structure  PDB id or model source. Empty = the dataset's pdb_ids.
#   force_field, water_model, software, ensemble, temperature_K, pressure_bar,
#   timestep_fs, frame_interval_ps
#                       listed ONLY where the glossary or the MD config could
#                       not establish them. A measured value always wins.
#   related_publication DOI or 'in preparation'. Empty = the dataset's citation.
#   keywords            added to the defaults. Letters, numbers, spaces, - _ .
#   slug                optional readable file name instead of the topology id.
#                       Do not change it once the URL is registered on NDP.
#
# Any other template field (lab, contact, license, data_availability) can be
# set under `defaults:` for every record, or inside one system to override it.
# So can the four the federation needs, whose lab-wide values are in lab.yaml:
#   data_owner          the PI who owns the data (name, email, orcid, affiliation)
#   access_level        public | private | licensable
#   ai_training         allowed | prohibited | by-agreement -- whether the data
#                       may train AI/ML models, separate from the licence
#   data_location       storage (which system) and uri (where to fetch it);
#                       the path is measured
# A key left empty here falls back to lab.yaml.
#
# `lengths:` chooses which MDs are a record's replicas: `nominal` (default) keeps
# only MDs that reach a length in `production_lengths_ns` above, `all` keeps
# every MD. The ones left out are listed on the registration sheet.
#
# `protocol:` chooses which stages the record lists: `production` (the default
# -- minimisation, heating and NPT equilibration are not published), `all` (the
# preparation history too) or `none`. It never changes sampling, which counts
# production only either way.
#
# `datasets:` holds what every system of one dataset shares (ions, water model,
# glycan profile, publication, keywords), so it is typed once. A system's own
# value wins; keywords add up. A dataset that appears later is not added there
# for you -- add its block by hand, above `systems:`.
#
# `systems:` must stay the LAST key of this section: new entries are appended
# to the end of the file.

ndp:
"""


#: What the systems of one dataset usually share.
SHARED = ["glycan_profile", "ions", "water_model", "related_publication", "keywords"]


#: `ndp: defaults:` keys a new settings.yaml leaves empty, so they keep
#: following lab.yaml; write a value there only where this folder differs.
LAB_LIVE = ("data_owner", "access_level", "ai_training", "data_location")


def _lab_value_text(v):
    """A short reminder of what lab.yaml answers, for the comment."""
    if isinstance(v, dict):
        v = v.get("name") or v.get("storage")
    return str(v).replace("#", "") if not _blank(v) else None


def scaffold_defaults(dataset_ids=(), indent="  "):
    """The `ndp:` section's own keys, ready to sit inside settings.yaml."""
    L = ["defaults:"]
    hints = {
        "base_url": "    # where records/ is hosted, e.g. https://raw.githubusercontent.com"
                    "/<org>/<repo>/main/records",
        "access_level": f"    # {' | '.join(ACCESS_CHOICES)}: who may get the raw files",
        "ai_training": f"    # {' | '.join(AI_CHOICES)}: AI/ML training, apart "
                       f"from the licence",
        "data_owner": "    # the PI who owns the data; only `name` is required",
        "data_location": "    # where the raw files live; `path` is measured",
    }
    for key, val in DEFAULTS.items():
        if key in LAB_LIVE:
            # Left empty, so lab.yaml stays the one place these are set: a
            # value copied here would no longer follow an edit to lab.yaml.
            said = _lab_value_text(val)
            L.append(f"  {key}:{hints.get(key, '')}"
                     + (f" -- empty: lab.yaml ({said})" if said else " -- empty: unanswered"))
        elif isinstance(val, list):
            L.append(f"  {key}: {json.dumps([str(v) for v in val], ensure_ascii=False)}")
        elif isinstance(val, dict):
            L.append(f"  {key}:{hints.get(key, '')}")
            L += [f"    {k}:" + (f" {json.dumps(v)}" if v is not None else "")
                  for k, v in val.items()]
        elif val is None:
            L.append(f"  {key}:{hints.get(key, '')}")
        elif key == "protocol":
            L.append(f"  {key}: {val}    # production (published) | all "
                     f"(+ min/heating/NPT) | none")
        elif key == "lengths":
            L.append(f"  {key}: {val}    # all (every replica, marked) | nominal (only full-length ones)")
        else:
            L.append(f"  {key}: {json.dumps(val)}")
    L += ["", "datasets:"]
    for did in dataset_ids:
        L.append(f"  {did}:")
        L += [f"    {k}: []" if k == "keywords" else f"    {k}:" for k in SHARED]
    L += ["", "systems:"]
    return "".join(f"{indent}{line}\n" if line else "\n" for line in L)


def curation_for(record_id, doc):
    """One system's curated values over its dataset's shared ones."""
    did = record_id.split("/", 1)[0]
    shared = dict((doc.get("datasets") or {}).get(did) or {})
    own = next((s for s in (doc.get("systems") or [])
                if isinstance(s, dict) and str(s.get("id")) == record_id), {})
    # `_said`, not `_blank`: an empty block in the system (`data_location:`
    # with an empty `path:` under it) is unanswered there, and must not hide
    # what the dataset says.
    out = {**shared, **{k: v for k, v in own.items() if _said(v) is not None}}
    # Where the files are is answered key by key -- a dataset's storage and
    # a system's own URL are both kept. Who owns them never is: one person's
    # name with another's e-mail would be worse than either alone.
    loc = [v for v in (shared.get("data_location"), own.get("data_location"))
           if isinstance(v, dict)]
    if len(loc) == 2:
        out["data_location"] = {**loc[0], **{k: v for k, v in loc[1].items()
                                             if _said(v) is not None}}
    out["keywords"] = list(shared.get("keywords") or []) + list(own.get("keywords") or [])
    return out


def declared_intervals(doc, dataset_id):
    """`{topology_group: frame_interval_ps}` a person declared for one dataset.

    The glossary is what this pipeline and MD_ANALYSES agree an MD *is*, and
    its `framestep_ns` is where the other side reads the frame spacing from.
    A spacing that reached only the record left the two disagreeing: on
    2026-09-24 `1-FULL-GP120`'s record said 100 ps while the glossary said
    nothing, and the analysis toolkit fell back to the post-processing
    placeholder in the DCD header -- 4.9e-05 ns -- and reported a 1050 ns
    replica as 0.513 ns. So the declaration is handed to the glossary too.

    Keyed by topology group because a record is `<dataset>/<topology_group>`
    and the glossary's rows carry that group; the dataset-level block has no
    `frame_interval_ps`, so only `systems:` is read.
    """
    out = {}
    for s in (doc.get("systems") or []):
        if not isinstance(s, dict):
            continue
        rid = str(s.get("id") or "")
        if "/" not in rid:
            continue
        did, group = rid.split("/", 1)
        if did != str(dataset_id) or _blank(s.get("frame_interval_ps")):
            continue
        try:
            ps = float(s["frame_interval_ps"])
        except (TypeError, ValueError):
            continue
        if ps > 0:
            out[group] = ps
    return out


#: settings.yaml keys that reach a record ONLY through the glossary, so an edit
#: to them needs glossary_build.py (ndp_pipeline.sh) before ndp_build.py sees it.
GLOSSARY_INPUTS = ("folder: stages:", "folder: replicas:", "production_lengths_ns:",
                   "frame_interval_ps")


def glossary_inputs(doc):
    """Digest of the settings the glossary is built from (GLOSSARY_INPUTS).

    glossary_build.py stores it beside the glossary; ndp_build.py recomputes it
    and says so when they differ -- running ndp_build.py alone after editing
    `stages:` left the old values in every record while reporting `ready`
    (SPIKE_CLEAVED, 2026-09-30)."""
    folder = doc.get("folder") or {}
    nd = dict(doc.get("ndp") or {})
    per = {}
    for d in doc.get("datasets") or []:
        did = str(d.get("id") or d.get("name") or "dataset")
        per[did] = [d.get("production_lengths_ns"), declared_intervals(nd, did)]
    payload = {"stages": folder.get("stages") or {},
               "replicas": folder.get("replicas") or [],
               "production_lengths_ns": doc.get("production_lengths_ns"),
               "datasets": per}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str)
                          .encode()).hexdigest()[:16]


def scaffold_entry(record_id, rows, measured, indent="  "):
    """One system's block: what was measured as a comment, the fields empty."""
    atoms = measured.get("total_atoms")
    head = (f"{len(rows)} MD(s), "
            + (f"{atoms:,} atoms, " if atoms else "")
            + (f"{measured['total_ns']} ns, " if measured.get("total_ns")
               else "length unknown, ")
            + (measured.get("software") or "program unknown"))
    L = ["", f"# {'-' * 70}", f"# {record_id}: {head}"]
    L += [f"#   {r.md_id}  ({r.source_dir})" for r in rows.itertuples()]
    # A length nothing could measure is usually one number away: these
    # trajectories were written by a tool that stated no timestep, and the
    # frame counts are solid. Say so where the field to fill actually is.
    if not measured.get("total_ns") and measured.get("frame_interval_ps") is None \
            and any(d.get("frames") for d in measured.get("replicas_detail") or []):
        frames = [d.get("frames") for d in measured["replicas_detail"] if d.get("frames")]
        L += ["#", f"#   No timestep in these trajectory headers, but the frames "
                   f"were counted ({', '.join(f'{f:,}' for f in frames[:5])}).",
              "#   Write `frame_interval_ps:` below and the lengths follow from it."]
    L += [f"# {'-' * 70}", f"- id: {record_id}"]
    asked = [k for k in FALLBACK if measured.get(k) is None
             and not (k == "pressure_bar" and measured.get("ensemble") in ("NVT", "NVE"))]
    for key in CURATED[:5] + asked + CURATED[5:]:
        L.append(f"  {key}: []" if key == "keywords" else f"  {key}:")
    return "".join(f"{indent}{line}\n" if line else "\n" for line in L)


class CurationLayout(Exception):
    """The `ndp:` section is laid out so appended entries would not join it."""


def sync_curation(folder, systems, measured_by_id):
    """Add a block to settings.yaml for every system it does not list yet.

    Returns (the `ndp:` section, added ids, ids listed but no longer in the
    glossary). Existing text is never rewritten -- comments and edits survive
    because the section is only ever appended to.
    """
    import yaml
    from core import settings
    doc = settings.ndp_of(folder)
    listed = {str(s.get("id")) for s in (doc.get("systems") or [])
              if isinstance(s, dict)}
    current = [rid for rid, _, _ in systems]
    missing = [rid for rid in current if rid not in listed]
    if missing:
        rows_by_id = {rid: rows for rid, rows, _ in systems}
        text = ""
        if not doc:
            text += scaffold_defaults(sorted({rid.split("/", 1)[0]
                                              for rid in current}))
        elif "systems" not in doc:
            text += "  systems:\n"
        text += "".join(scaffold_entry(rid, rows_by_id[rid], measured_by_id[rid])
                        for rid in missing)
        # Appending only lands inside `ndp: systems:` while that is the last key
        # of the last section. Check on the text BEFORE touching the file.
        before = settings.read(folder)
        try:
            after = yaml.safe_load(before + text) or {}
            now = {str(s.get("id")) for s in ((after.get("ndp") or {}).get("systems") or [])
                   if isinstance(s, dict)}
        except yaml.YAMLError:
            now = set()
        if not set(missing) <= now:
            raise CurationLayout(
                f"{settings.path(folder)}: cannot append {len(missing)} new "
                f"system(s). Keep `systems:` as the LAST key of the `ndp:` "
                f"section, and `ndp:` as the last section. Nothing was written.")
        settings.write(folder, before + text)
        doc = settings.ndp_of(folder)
    gone = sorted(listed - set(current))
    return doc, missing, gone


#: Record fields a person can answer in a system's own block. Anything else a
#: record can be missing (atoms, lengths) is measured, and no line in
#: settings.yaml could supply it.
FILLABLE = (set(CURATED) | set(FALLBACK)
            | {"authors", "contact", "license", "lab", "data_availability",
               "data_owner", "access_level", "ai_training", "data_location"})


def _key_of(field):
    """`system.construct` -> ("construct", None); `data_location.storage` ->
    ("data_location", "storage"): the settings.yaml key a missing field is
    answered on, and the sub-key when it is a block."""
    head, _, tail = field.partition(".")
    if head in ("data_owner", "data_location"):
        return head, tail or None
    key = field.rsplit(".", 1)[-1]
    return (key, None) if key in FILLABLE else (None, None)


def ensure_lines(folder, missing_by_id):
    """Give every field a record still misses an empty line in that system's
    own block of settings.yaml, so the place to answer it is always there.

    A block is written once, when its system first appears; a field that goes
    missing later (a frame spacing a re-index could no longer read, a field
    added to the template since) had no line, and a person had to know the
    key to type it. Existing text is never rewritten: the lines go at the end
    of the block, one block at a time, and each insertion is kept only if the
    file then parses to exactly the same values plus the new, empty keys.
    Line endings are kept as they are.

    Returns ({record_id: [keys added]}, {record_id: why it was left alone}).
    """
    import yaml
    from core import settings
    with open(settings.path(folder), newline="") as fh:
        raw = fh.read()
    eol = "\r\n" if "\r\n" in raw else "\n"
    lines = raw.split(eol)

    def parse(text_lines):
        try:
            return json.loads(json.dumps(yaml.safe_load(eol.join(text_lines)) or {},
                                         default=str))
        except yaml.YAMLError:
            return None

    current = parse(lines)
    if current is None:
        return {}, {}
    added, skipped = {}, {}
    for rid, missing in missing_by_id.items():
        blocks = [b for b in ((current.get("ndp") or {}).get("systems") or [])
                  if isinstance(b, dict) and str(b.get("id")) == rid]
        if len(blocks) != 1:
            continue
        block = blocks[0]
        want = {}
        for field in missing:
            key, sub = _key_of(field)
            if key is None or key in block:
                # Present, even empty, is where the person answers; a second
                # key would silently replace it.
                continue
            want.setdefault(key, [])
            if sub and sub not in want[key]:
                want[key].append(sub)
        if not want:
            continue
        at = [i for i, l in enumerate(lines)
              if re.match(rf"^(\s*)- id:\s*['\"]?{re.escape(rid)}['\"]?\s*(#.*)?$", l)]
        if len(at) != 1:
            skipped[rid] = "its `- id:` line was not found exactly once"
            continue
        at = at[0]
        pad = len(lines[at]) - len(lines[at].lstrip()) + 2

        def inside(line):
            return line.strip() and len(line) - len(line.lstrip()) >= pad

        end, last = at + 1, at
        while end < len(lines):
            if inside(lines[end]):
                last = end
            elif lines[end].strip():
                break
            else:
                # A blank line ends the block only if what follows is not
                # part of it (a `|` text block can hold blank lines).
                nxt = next((l for l in lines[end:] if l.strip()), "")
                if not inside(nxt):
                    break
            end += 1
        new = []
        for key, subs in want.items():
            if subs:
                new.append(" " * pad + f"{key}:")
                new += [" " * (pad + 2) + f"{sub}:" for sub in subs]
            else:
                new.append(" " * pad + f"{key}:" + " " * max(1, 22 - len(key))
                           + "# missing: fill in for this system")
        trial = lines[:last + 1] + new + lines[last + 1:]
        after = parse(trial)
        expect = json.loads(json.dumps(current))
        if after is not None:
            for b in ((after.get("ndp") or {}).get("systems") or []):
                if isinstance(b, dict) and str(b.get("id")) == rid:
                    for key in want:
                        b.pop(key, None)
        if after != expect:
            skipped[rid] = "the block's layout could not take new lines safely"
            continue
        lines, current = trial, parse(trial)
        added[rid] = list(want)
    if added:
        settings.write(folder, eol.join(lines))
    return added, skipped


# --------------------------------------------------------------------------
# the record
# --------------------------------------------------------------------------
def _said(v):
    """A curated value as written; 'none'/'n/a' is an answer, empty is not.
    A block whose every key is empty (`data_owner:` with nothing under it) is
    empty too, so the level below it still answers."""
    if _blank(v) or (isinstance(v, dict) and all(_blank(x) for x in v.values())):
        return None
    return v.strip() if isinstance(v, str) else v


def _authors(v):
    """The author line a reader sees. A YAML list becomes one string.

    `datasets:` holds authors as a list, because that is how a person types
    them and how a citation manager wants them. Every consumer of the record
    -- the catalogue, the analysis page -- shows one line, so the joining is
    done once, here, rather than differently in each of them.
    """
    if _blank(v):
        return None
    if isinstance(v, str):
        return v.strip() or None
    names = [str(x).strip() for x in v if not _blank(x)]
    return ", ".join(names) or None


def _choice(key, v, choices, aliases, warn):
    """One of `choices`, from what a person wrote; None (and a warning) for
    anything else. None stays None: unanswered, not refused."""
    if v is None:
        return None
    said = str(v).strip().lower()
    said = aliases.get(said, said)
    if said in choices:
        return said
    if warn is not None:
        warn.append(f"`{key}: {v}` is not one of {', '.join(choices)}; left null.")
    return None


def _owner(v):
    """The data owner as the record carries it: always the same four keys.
    A plain string is a name."""
    keys = ("name", "email", "orcid", "affiliation")
    if isinstance(v, dict):
        return {k: _said(v.get(k)) for k in keys}
    return {"name": _said(v) if not isinstance(v, (list, dict)) else None,
            **{k: None for k in keys[1:]}}


def _location(levels, detail):
    """Where the raw files live: `storage` and `uri` as declared, each from
    the most specific level that says it (a dataset's own URL over the lab's
    storage name), and `path` as MEASURED from the replicas -- declared only
    where no replica could say where it is."""
    got = {"storage": None, "uri": None, "path": None}
    for v in reversed(levels):              # lowest precedence first
        if isinstance(v, str) and _said(v):
            v = {"storage": v}
        if isinstance(v, dict):
            got.update({k: _said(v.get(k)) for k in got if _said(v.get(k)) is not None})
    return {"storage": got["storage"], "path": raw_location(detail) or got["path"],
            "uri": got["uri"]}


def _pdb_text(spec):
    ids = [str(p) for p in (spec.get("pdb_ids") or []) if not _blank(p)]
    return ("PDB " + ", ".join(ids)) if ids else None


def describe(record_id, measured, curated, spec, defaults):
    """The record's own description, composed from what was measured.

    The catalogue shows a paragraph under every title, and a person should not
    have to write it before the dataset can be listed. `description:` in
    settings.yaml replaces it whenever they would rather say it themselves.
    """
    stub = {
        "title": _said(curated.get("title")),
        "system": {"construct": _said(curated.get("construct")),
                   "starting_structure": _said(curated.get("starting_structure"))
                   or _pdb_text(spec),
                   "glycan_profile": _said(curated.get("glycan_profile")),
                   "total_atoms": measured.get("total_atoms")},
        "setup": {"force_field": measured.get("force_field"),
                  "water_model": measured.get("water_model")
                  or _said(curated.get("water_model")) or spec.get("wat"),
                  "ions": _said(curated.get("ions")),
                  "software": measured.get("software"),
                  "ensemble": measured.get("ensemble"),
                  "temperature_K": measured.get("temperature_K"),
                  "pressure_bar": measured.get("pressure_bar"),
                  "timestep_fs": measured.get("timestep_fs")},
        "sampling": {"frame_interval_ps": measured.get("frame_interval_ps")},
        "data_availability": _said(curated.get("data_availability"))
        or defaults.get("data_availability"),
        "related_publication": _said(curated.get("related_publication"))
        or spec.get("citation"),
        "contact": _said(curated.get("contact")) or spec.get("contact"),
        "protocol": [],
    }
    return catalog_description(stub, measured)


def lengths_from_spacing(sampling):
    """Fill the lengths a declared frame interval makes knowable.

    Trajectories written by a post-processing tool carry no timestep:
    `1-FULL-GP120` holds five 10,500-frame DCDs whose headers state no `delta`,
    so every length in its record was null while the frame counts were solid.
    One number in settings.yaml -- `frame_interval_ps:` under this system's
    `ndp:` block -- turns those counts into the lengths the catalogue asks for.
    Each value says it was derived, because it rests on what a person declared.
    """
    ps = sampling.get("frame_interval_ps")
    if not ps or sampling.get("total_ns") is not None:
        return sampling
    lengths = []
    for d in sampling.get("replicas_detail") or []:
        if d.get("length_ns") is None and d.get("frames"):
            d["length_ns"] = round(d["frames"] * float(ps) / 1000.0, 1)
            d["length_basis"] = "declared frame interval"
        if d.get("length_ns") is not None:
            lengths.append(d["length_ns"])
    if not lengths:
        return sampling
    sampling["replica_lengths_ns"] = sampling.get("replica_lengths_ns") or lengths
    sampling["total_ns"] = round(sum(lengths), 1)
    if sampling.get("length_per_replica_ns") is None and len(set(lengths)) == 1:
        sampling["length_per_replica_ns"] = lengths[0]
    sampling["length_basis"] = "declared frame interval"
    return sampling


def build_record(record_id, rows, measured, curated, spec, defaults, warn=None,
                 stages=None):
    """The JSON document for one system, in the template's field order.

    `warn`, when given, collects curated values that were ignored because the
    measurement says otherwise.
    """
    def pick(key, *fallbacks):
        for v in (curated.get(key), *fallbacks, defaults.get(key)):
            if _said(v) is not None:
                return _said(v)
        return None

    def est(key, *fallbacks):
        """Measured first; settings.yaml only where nothing was measured."""
        got, said = measured[key], _said(curated.get(key))
        if got is None:
            return pick(key, *fallbacks)
        if said is not None and str(said) != str(got) and warn is not None:
            warn.append(f"`{key}: {said}` in settings.yaml is ignored: the MDs say "
                        f"{got}. If they are wrong, fix it where it was read.")
        return got

    lic = pick("license", spec.get("license"))
    keep = str(pick("protocol") or DEFAULTS["protocol"]).strip().lower()
    if keep not in PROTOCOL_CHOICES:
        if warn is not None:
            warn.append(f"`protocol: {keep}` is not one of "
                        f"{', '.join(PROTOCOL_CHOICES)}; all stages are listed.")
        keep = "all"
    stages = [] if keep == "none" else [
        s for s in (stages or []) if keep == "all" or s.get("stage") == "production"]
    record = {
        # --- what the lab catalogue lists ---------------------------------
        "identifier": pick("identifier"),
        "topic": pick("topic"),
        "tags": list(curated.get("tags") or defaults.get("tags") or []),
        "title": pick("title"),
        "description": pick("description") or describe(record_id, measured,
                                                        curated, spec, defaults),
        "size_gb": measured.get("size_gb"),
        "files": measured.get("files"),
        # --- the NDP template ---------------------------------------------
        "lab": pick("lab"),
        "contact": pick("contact", spec.get("contact")),
        "system": {
            "construct": pick("construct"),
            "starting_structure": pick("starting_structure", _pdb_text(spec)),
            "glycan_profile": pick("glycan_profile"),
            "total_atoms": measured["total_atoms"],
        },
        "setup": {
            "force_field": est("force_field"),
            "water_model": est("water_model", spec.get("wat")),
            "ions": pick("ions"),
            "software": est("software"),
            "ensemble": est("ensemble"),
            "temperature_K": est("temperature_K"),
            "pressure_bar": est("pressure_bar"),
            "timestep_fs": est("timestep_fs"),
        },
        "sampling": {
            "replicas": measured["replicas"],
            "complete_replicas": measured.get("complete_replicas"),
            "planned_length_ns": measured.get("planned_length_ns"),
            "length_per_replica_ns": measured["length_per_replica_ns"],
            "frame_interval_ps": est("frame_interval_ps"),
            # The declared production length every replica reaches (glossary
            # nominal_ns): 1000 for 3 x 1008 ns. Null when unknown or mixed.
            "nominal_length_ns": measured.get("nominal_ns"),
            "replica_lengths_ns": measured["replica_lengths_ns"],
            "total_ns": measured["total_ns"],
            # "declared frame interval" when the lengths rest on a spacing a
            # person supplied rather than one the headers stated. Null means
            # measured, which is the ordinary case.
            "length_basis": measured.get("length_basis"),
            # Every replica, named, with how far it got: a run that stopped is
            # part of the record, not something the catalogue hides.
            "replicas_detail": measured.get("replicas_detail") or [],
        },
        "data_availability": pick("data_availability"),
        "related_publication": pick("related_publication", spec.get("citation")),
        "license": SPDX.get(str(lic).lower(), lic) if lic else None,
        # --- federation and licensing ---------------------------------------
        # Who owns the data (the PI), who may get the raw files, whether they
        # may train an AI model -- a permission the citation licence does not
        # grant or refuse -- and where the files actually are. Lab-wide values
        # come from lab.yaml; any level of settings.yaml overrides them.
        "data_owner": _owner(pick("data_owner", spec.get("data_owner"))),
        "access_level": _choice("access_level",
                                pick("access_level", spec.get("access_level")),
                                ACCESS_CHOICES, {}, warn),
        "ai_training": _choice("ai_training",
                               pick("ai_training", spec.get("ai_training")),
                               AI_CHOICES, AI_ALIASES, warn),
        "data_location": _location(
            [curated.get("data_location"), spec.get("data_location"),
             defaults.get("data_location")],
            measured.get("replicas_detail")),
        # Who to cite. NOT one of the NDP form's fields -- that form asks for a
        # contact and a publisher -- so it sits after the template rather than
        # inside it. The lab catalogue and the analysis page both show an author
        # line, and both were showing "to be added" for a value settings.yaml
        # had held all along: `authors:` under the dataset. One source, joined
        # by _authors(), so the two can never print a different list.
        "authors": _authors(pick("authors", spec.get("authors"))),
        # Whether the paragraph in `description` is a PERSON'S or this tool's. The
        # catalogue shows one either way, so a composed one is not a gap --
        # but a published page should print somebody's own words and not a
        # sentence assembled from the atom count, and it cannot tell them
        # apart from the text alone.
        "description_basis": ("declared" if _said(pick("description"))
                              else "composed"),
        # Every stage from minimisation to the end of production, read from
        # each stage's config, log and trajectory (protocol.py). [] when none
        # could be read. Not in the guide's template: it is what "how was this
        # prepared" means for a reader, and it costs nothing to carry.
        "protocol": stages,
        "lab_glossary": {
            "record_id": record_id,
            "md_keys": [str(k) for k in rows.md_key],
        },
    }
    lengths_from_spacing(record["sampling"])
    return record


def missing_fields(record):
    """Dotted names of template fields still null. NVT has no pressure to give."""
    out = []
    for key, val in record.items():
        if key in ("lab_glossary", "protocol", "tags", "size_gb", "files"):
            continue
        if isinstance(val, dict):
            # Derived or measured values are never a person's to fill in: an
            # absent planned length means the config did not state one.
            # `length_basis` joins these because NULL IS ITS NORMAL VALUE:
            # it is filled only when the lengths rest on a declared frame
            # interval, so an empty one means "measured", not "unanswered",
            # and listing it as missing made every ordinary record
            # incomplete for a question nobody had asked.
            derived = ("nominal_length_ns", "planned_length_ns",
                       "complete_replicas", "replicas_detail", "length_basis")
            out += [f"{key}.{k}" for k, v in val.items()
                    if _blank(v) and k not in derived]
        elif _blank(val):
            out.append(key)
    if record["setup"]["ensemble"] in ("NVT", "NVE") and "setup.pressure_bar" in out:
        out.remove("setup.pressure_bar")
    return [f for f in out if f not in OPTIONAL]


def keywords_of(curated, defaults):
    """(keywords, rejected) -- defaults first, then the system's own."""
    good, bad = [], []
    for kw in list(defaults.get("keywords") or []) + list(curated.get("keywords") or []):
        kw = str(kw).strip()
        if not kw or kw.lower() in (g.lower() for g in good):
            continue
        (good if KEYWORD_OK.match(kw) else bad).append(kw)
    return good, bad


def check_text(record):
    """Warnings the NDP form would not give: stray quotes, a contact that is
    not an address."""
    warn = []
    title = record.get("title")
    if title and (title[0] in "\"'“”‘’" or title[-1] in "\"'“”‘’"):
        warn.append(f"`title` starts or ends with a quote character: {title!r}")
    contact = record.get("contact")
    if contact and not re.match(r"^[^@\s]+@[^@\s]+\.[^@\s]+$", str(contact)):
        warn.append(f"`contact` is not an e-mail address: {contact!r}")
    return warn


#: Path components that name a place, not a system: they say nothing in a file
#: name, and every folder has them.
GENERIC = {"md", "mds", "run", "runs", "replica", "replicas", "rep", "reps",
           "data", "sim", "sims", "simulation", "simulations", "trajectories",
           "trajectory", "production", "prod", "analysis", "analyses"}


def _slug(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(text)).strip("-. ")


def _is_generic(name):
    """`MD`, `0-MD`, `REPLICA`, `runs2` -- a place, not a system.

    The lab numbers its directories (`0-MD`, `1-GLYC`, `3-YZ23`), so the index
    prefix is stripped before the word is judged: `0-MD` is as generic as `MD`,
    while `3-YZ23` is not.
    """
    word = re.sub(r"^\d+[-_]|[-_]?\d+$", "", str(name)).strip("-_ ").lower()
    return word in GENERIC or not word


def folder_tail(folder, parts=2):
    """The last `parts` telling components of a folder's path.

    `.../ABZYME/3-OCCLUDED-OPEN/MDs/2-REPLICA` -> `3-OCCLUDED-OPEN_2-REPLICA`.
    A person named these directories after what they hold, so they are the
    readable part of any name we can build without being told one.
    """
    names = [n for n in str(folder).strip("/").split("/")
             if n and not _is_generic(n)]
    return "_".join(_slug(n) for n in names[-parts:] if _slug(n))


def sampling_tag(measured):
    """`3x1us`, `5x500ns`, `5-replicas` -- replicas and how long each ran.

    Takes either the measurement or a finished record's `sampling` block, so a
    length recovered from a declared frame interval names the file too.
    """
    n = measured.get("replicas") or 0
    length = (measured.get("nominal_ns") or measured.get("nominal_length_ns")
              or measured.get("length_per_replica_ns"))
    if not length:
        # No length anywhere: say so in words rather than in a field that
        # looks like a measurement ("5x-unknown" read as a length of unknown).
        return f"{n}-replicas" if n else "" 
    length = float(length)
    return (f"{n}x{length/1000:g}us" if length >= 1000 else f"{n}x{length:g}ns")


def common_subpath(paths, parts=2):
    """What the MDs of one record have in common, below the folder.

    `1-1F4/0-MD/{1-R1,2-R2,3-R3}` -> `1-1F4`. Four systems in one folder would
    otherwise all be named after the folder and collide; what tells them apart
    is exactly the directory a person put each system in.
    """
    if not paths:
        return ""
    split = [str(p).strip("/").split("/") for p in paths]
    if len(split) == 1 and len(split[0]) > 1:
        # One replica: its own directory IS the whole path, and the last
        # component is the replica (`1R`), not the system.
        split = [split[0][:-1]]
    common = []
    for i in range(min(len(p) for p in split)):
        names = {p[i] for p in split}
        if len(names) != 1:
            break
        common.append(names.pop())
    telling = [c for c in common if not _is_generic(c)]
    # Everything generic (`MD/REPLICA`) is still the only thing distinguishing
    # this record, so keep its first component rather than nothing.
    chosen = (telling or common[:1])[-parts:]
    return "_".join(_slug(c) for c in chosen if _slug(c))


def file_name(record_id, curated, folder=None, measured=None, paths=(),
              existing=None):
    """The record's file name: where it lives and how much of it there is.

    `1a918746.json` -- the topology hash -- was unreadable, which is the whole
    complaint it answers. In order: an explicit `slug:`, then the name this
    record ALREADY has, then the folder plus what this record's own MDs have in
    common (`CDR_1-1F4`); each followed by the sampling, so two records of one
    system are still told apart. The topology id only ever breaks a remaining
    tie.

    A NAME, ONCE GIVEN, STAYS. The title used to come second, which meant that
    filling in `title:` renamed every record in the folder -- on 2026-09-24
    `CDR_1-1F4_3x500ns.json` became
    `HIV---Env-Open-State---Tribound-1F4-antibody_3x500ns.json`, orphaning the
    analysis bundle and the production files keyed on the old stem, and
    turning " - " into "---" on the way. A file name is part of a published
    record's identity: a title is a sentence for people to read, not an
    address. Only an explicit `slug:` renames a record now, which is a thing
    someone did on purpose.
    """
    _, gid = record_id.split("/", 1)
    curated = curated or {}
    stem = _slug(_said(curated.get("slug")) or "")
    if not stem:
        if existing:
            return existing
        else:
            stem = "_".join(x for x in (folder_tail(folder or "", parts=1),
                                        common_subpath(paths)) if x)
    # The same rule for a slugged name: once `<slug>_6x1us.json` exists, a
    # changed replica count or frame interval no longer renames it -- the
    # hosted URL is the record's address (SPIKE_CLEAVED went through three
    # names in two days, 2026-09-30). A new slug, or `rename: true`, renames.
    elif existing and not _said(curated.get("rename")) \
            and _slugged_by(existing, stem):
        return existing
    tag = sampling_tag(measured or {})
    stem = "_".join(x for x in (stem, tag) if x)
    return f"{stem or gid}.json"


def _slugged_by(existing, stem):
    name = os.path.basename(existing)
    return name == f"{stem}.json" or name.startswith(f"{stem}_")


def slug_would_rename(record_id, curated, measured, existing):
    """True when a kept slugged name differs from the one the sampling gives now."""
    curated = curated or {}
    stem = _slug(_said(curated.get("slug")) or "")
    if not stem or not existing or not _slugged_by(existing, stem):
        return False
    tag = sampling_tag(measured or {})
    fresh = "_".join(x for x in (stem, tag) if x) + ".json"
    return os.path.basename(existing) != fresh


def declared_notes(rows, curated, warn):
    """Notes for what a person declared in the `ndp:` block that the MDs could
    not show, and the two warnings such a declaration answers, taken out of
    `warn` in place. Nothing here changes a record.

      * `nominal_lengths_ns: [1000, 400]` -- a system whose replicas were run
        to more than one length on purpose (conventional + adaptive sampling).
        When every replica reaches one of them, and each is reached, the
        "MDs disagree" and "very uneven" warnings are expected and dropped.
      * `config_declared` -- MDs with no MD config beside them (`no_config`)
        whose ensemble, temperature and timestep the block states by hand.
        The glossary keeps them at `review`; the note says why that is fine.
    """
    curated = curated or {}
    notes = []
    declared = curated.get("nominal_lengths_ns")
    if isinstance(declared, (int, float)):
        declared = [declared]
    if isinstance(declared, list) and declared:
        try:
            targets = sorted({float(v) for v in declared if float(v) > 0})
        except (TypeError, ValueError):
            targets = []
        lengths = [float(v) for v in rows.simulated_ns if not _nan(v)]
        hit = {t: 0 for t in targets}
        fits = bool(targets) and len(lengths) == len(rows)
        for v in lengths:
            reach = [t for t in targets if v >= (1 - UNIFORM) * t]
            if not reach:
                fits = False
                break
            hit[max(reach)] += 1
        if fits and all(hit.values()):
            before = len(warn)
            warn[:] = [w for w in warn
                       if not w.startswith("`nominal_length_ns`: the MDs disagree")
                       and not w.startswith("replica lengths are very uneven")]
            notes.append("replicas run to the declared nominal lengths "
                         + " + ".join(f"{hit[t]} x {t:g} ns" for t in reversed(targets))
                         + (f"; {before - len(warn)} length warning(s) expected and "
                            f"not counted" if before != len(warn) else ""))
        else:
            notes.append(f"`nominal_lengths_ns: {declared}` does not match the "
                         f"replica lengths; the length warnings stand")
    if "issues" in rows.columns:
        n_cfg = sum("no_config" in str(i) for i in rows.issues.fillna(""))
        stated = [k for k in ("ensemble", "temperature_K", "timestep_fs")
                  if _said(curated.get(k)) is not None]
        if n_cfg and len(stated) == 3:
            notes.append(f"config_declared: {n_cfg} MD(s) have no MD config "
                         f"(no_config); ensemble, temperature and timestep are "
                         f"taken from settings.yaml")
    return notes


def swallowed_comments(doc):
    """Values in settings.yaml that swallowed a comment: `a@b.edu# note` is
    one string to YAML, because a comment needs a space before its `#`."""
    out = []

    def walk(node, where):
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{where}{k}:" if not where else f"{where} {k}:")
        elif isinstance(node, list):
            for v in node:
                walk(v, where)
        elif isinstance(node, str) and re.search(r"\S#\s", node):
            out.append(f"`{where} {node}` looks like a comment was read as part "
                       f"of the value -- put a space before the `#`, or quote it")

    for section in ("datasets", "ndp"):
        walk(doc.get(section), f"{section}:")
    return out


def dumps(record):
    """Deterministic text: a re-run that changes nothing changes no file."""
    return json.dumps(record, indent=2, ensure_ascii=False) + "\n"


def digest(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


# --------------------------------------------------------------------------
# the registration sheet: what to paste into the NDP form
# --------------------------------------------------------------------------
def sampling_phrase(measured):
    n, lengths = measured["replicas"], measured["replica_lengths_ns"]
    if not lengths:
        return f"{n} replica(s)"
    if len(lengths) == n and min(lengths) >= (1 - UNIFORM) * max(lengths):
        return f"{n} × {round(measured['length_per_replica_ns']):g} ns"
    return f"{n} replicas, {round(measured['total_ns']):g} ns in total"


def _sampling_of(record, measured):
    """The sampling the prose should describe: the record's, once it has one.

    `measured` is what the trajectory headers said. The record's `sampling`
    is that PLUS the lengths a declared frame interval makes knowable
    (`lengths_from_spacing`), and the prose must not be written from the
    poorer of the two: on 2026-09-24 `1-FULL-GP120`'s own JSON said
    `total_ns: 5250.0` while its description said "Sampling: 5 replica(s)"
    and the registration sheet printed `total ns` as a dash -- one document
    contradicting itself, in the very box a person pastes into the NDP form.
    """
    got = record.get("sampling") or {}
    # Only when it is the REAL sampling block. `describe()` composes the
    # auto-description from a stub record whose `sampling` holds just the
    # frame interval, and reading the replica count out of that one raised
    # KeyError on every build.
    if all(k in got for k in ("replicas", "replica_lengths_ns")):
        return got
    return measured


def catalog_title(record, measured):
    return (f"{record['title']} — all-atom MD, "
            f"{sampling_phrase(_sampling_of(record, measured))}"
            if record.get("title") else None)


def catalog_description(record, measured):
    """The Description box, composed only from what the record establishes."""
    s, u = record["system"], record["setup"]
    P = []
    lead = f"Molecular dynamics simulations of {s['construct'] or record['title'] or 'this system'}"
    if s["total_atoms"]:
        lead += f" ({s['total_atoms']:,} atoms)"
    P.append(lead + ".")
    if s["starting_structure"]:
        P.append(f"Starting structure: {s['starting_structure']}.")
    if s["glycan_profile"] and str(s["glycan_profile"]).lower() not in ("none", "n/a"):
        P.append(f"Glycan profile: {s['glycan_profile']}.")
    setup = [f"force field {u['force_field']}" if u["force_field"] else None,
             f"{u['water_model']} water" if u["water_model"] else None,
             u["ions"], u["software"]]
    if any(setup):
        P.append("Setup: " + ", ".join(x for x in setup if x) + ".")
    cond = [u["ensemble"],
            f"{u['temperature_K']} K" if u["temperature_K"] is not None else None,
            f"{u['pressure_bar']} bar" if u["pressure_bar"] is not None else None,
            f"{u['timestep_fs']} fs timestep" if u["timestep_fs"] is not None else None]
    if any(cond):
        P.append("Conditions: " + ", ".join(x for x in cond if x) + ".")
    samp = f"Sampling: {sampling_phrase(_sampling_of(record, measured))}"
    if record["sampling"]["frame_interval_ps"] is not None:
        samp += f", one frame every {record['sampling']['frame_interval_ps']} ps"
    P.append(samp + ".")
    if record.get("protocol"):
        from .protocol import summary
        P.append(f"Protocol: {summary(record['protocol'])}.")
    if record["related_publication"]:
        P.append(f"Related publication: {record['related_publication']}.")
    avail = str(record["data_availability"] or "").strip().rstrip(".")
    P.append("This entry is metadata only."
             + (f" {avail}." if avail else "")
             + (f" Contact: {record['contact']}." if record["contact"] else ""))
    return " ".join(P)
