#!/usr/bin/env python3
"""
mwf.py -- assemble one MDDB-workflow input file from the glossary.

The workflow (`mwf`) reads a YAML with the project's metadata and a list of
`mds`, each naming a directory and its trajectory files. Spec:
https://mddb-workflow.readthedocs.io/en/latest/input.html

WHAT THIS DOES AND DOES NOT DECIDE. Every value here was either written by a
person in `settings.yaml` (`datasets:`) or measured by the glossary. Nothing is
inferred, and a field with no evidence is OMITTED rather than written as null:
mwf reads a null as "measured and found to be nothing", which is a different
claim from "not established".

WHERE A VALUE IS PLACED is decided by the evidence, not assumed. One value
shared by every MD of a dataset is a property of the dataset and is written
once; two or more are properties of the individual MDs and are written on their
own `mds` entries, which the spec allows and the workflow uploads as that MD's
own metadata. Dropping them, as an earlier version did, threw away ninety
measured frame spacings to avoid asserting a single wrong one.
"""

import os

from . import atomic, identity

#: Fields a person writes in settings.yaml, passed through as they are.
CURATED_KEYS = (
    "name", "description", "authors", "contact", "type", "method", "license",
    "linkcense", "citation", "thanks", "collections", "groups", "pdb_ids",
    "forced_references", "ligands", "wat", "multimeric", "interactions",
    "pbc_selection", "accession", "framestep",
)

#: Measured fields, and the glossary column each comes from. `mds` entries may
#: carry any of them; the project carries the ones every MD agrees on.
MEASURED = {
    "program": "program", "version": "version", "temp": "temp_K",
    "timestep": "timestep_fs", "ensemble": "ensemble", "ff": "ff",
    "wat": "wat", "boxtype": "boxtype", "framestep": "framestep_ns",
}

#: What the spec allows on an individual `mds` entry.
PER_MD_KEYS = frozenset({"program", "version", "temp", "timestep", "ensemble",
                         "ff", "wat", "boxtype", "framestep"})


def _clean(value):
    if value is None or (isinstance(value, float) and value != value):
        return None
    if hasattr(value, "item"):
        try:
            value = value.item()
        except (AttributeError, ValueError):
            pass
    return value


def _coerce(key, value):
    """The type the MDDB spec gives this field.

    `temp` and `timestep` are numbers and `ff` is a string or a list of them,
    so a temperature read as 310.0 is written 310 and a single force field is
    written as a one-item list rather than a bare string.
    """
    value = _clean(value)
    if value is None or value == "":
        return None
    if key == "temp":
        return int(round(float(value)))
    if key == "timestep":
        return float(value)
    if key == "ff":
        if isinstance(value, str):
            return [v.strip() for v in value.split(",") if v.strip()]
        listed = [str(v) for v in list(value)]
        return listed or None
    return value


def _compare_key(value):
    """A hashable stand-in for a measured value, so lists compare by content."""
    if isinstance(value, (list, tuple)):
        return tuple(str(v) for v in value)
    return str(value)


def _same(values):
    """The one value they all carry, or None when they disagree or are empty."""
    seen = [v for v in values if v is not None]
    if not seen or len(seen) != len(values):
        return None
    first = seen[0]
    key = (lambda v: tuple(v) if isinstance(v, (list, tuple)) else v)
    return first if all(key(v) == key(first) for v in seen) else None


def build_dataset(spec, rows, root):
    """(document, coverage) for one dataset, from its glossary rows.

    `rows` is the glossary's own view of the MDs: ids, source directories,
    segments in time order, topologies and every measured value. The document
    still points at the SOURCE; workdir.materialise rewrites it to the
    link-only working area before anything runs.
    """
    doc, cover = {}, []

    def record(field, value, source, note=""):
        if value is None or (isinstance(value, (list, dict)) and not value):
            cover.append((field, "absent", note or "no evidence"))
            return
        doc[field] = value
        cover.append((field, "filled", source))

    for key in CURATED_KEYS:
        if key in spec and spec[key] not in (None, "", [], {}):
            record(key, spec[key], "settings.yaml (datasets:)")
        elif key in ("name", "authors", "license", "citation"):
            cover.append((key, "absent", "add it to settings.yaml"))

    mds, per_md = [], {}
    for r in rows.sort_values("md_id").itertuples():
        segs = [s for s in str(r.segments).split("; ") if s]
        if not segs:
            continue
        entry = {"name": str(r.md_id), "mdir": str(r.source_dir),
                 "input_trajectory_filepaths": segs}
        top = _clean(getattr(r, "topology", None))
        if top:
            entry["input_topology_filepath"] = top
        mds.append(entry)
        per_md[str(r.md_id)] = {k: _coerce(k, getattr(r, col, None))
                                for k, col in MEASURED.items()}

    if not mds:
        cover.append(("mds", "absent", "no MD in this dataset holds segments"))
        return {**doc, "mds": [], "mdref": 0}, cover

    # Fail BEFORE anything is written: mwf rejects the whole project for a
    # duplicate MD name, and only after it has rewritten inputs.yaml in place.
    identity.check_unique([(m["name"], m["mdir"]) for m in mds], "MD name")

    # One topology for the project when every MD ran from the same file;
    # otherwise each entry keeps its own, which is what tells two molecular
    # systems apart inside one dataset.
    tops = [m.get("input_topology_filepath") for m in mds]
    shared_top = _same(tops)
    if shared_top:
        record("input_topology_filepath", shared_top, "the glossary: every MD "
               "runs from this topology")
        for m in mds:
            m.pop("input_topology_filepath", None)
    elif any(tops):
        cover.append(("input_topology_filepath", "per-MD",
                      f"{len({t for t in tops if t})} different topologies; "
                      f"written on each `mds` entry"))

    for key in sorted(MEASURED):
        values = [per_md[m["name"]].get(key) for m in mds]
        shared = _same(values)
        if shared is not None:
            record(key, shared, "the glossary: every MD agrees")
        elif any(v is not None for v in values) and key in PER_MD_KEYS:
            for m in mds:
                v = per_md[m["name"]].get(key)
                if v is not None:
                    m[key] = v
            cover.append((key, "per-MD", "the MDs differ, or only some were "
                                         "measured; written on each entry"))
        else:
            cover.append((key, "absent", "not established for any MD"))

    doc["mds"] = mds
    doc["mdref"] = 0
    cover.append(("mds", "filled", f"{len(mds)} MD(s) from the glossary"))
    return doc, cover


def dump(doc, path, dataset_id=None):
    """Write the YAML with `mds` last, so a person reading the file sees the
    descriptive metadata before a long list of directories."""
    import yaml
    ordered = {k: v for k, v in doc.items()
               if k not in ("mds", "mdref") and not k.startswith("_")}
    ordered["mdref"] = doc.get("mdref", 0)
    ordered["mds"] = doc.get("mds", [])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    did = dataset_id or os.path.basename(os.path.dirname(path))

    def _w(tmp):
        with open(tmp, "w") as fh:
            fh.write("# mwf input file -- generated from the MD glossary.\n"
                     "# Spec: https://mddb-workflow.readthedocs.io/en/latest/input.html\n"
                     "# Fields absent here had no evidence on disk.\n"
                     "#\n"
                     "# DO NOT RUN mwf ON THIS FILE IN PLACE. Every mdir and input\n"
                     "# below names a link inside .md_pipeline/mddb/work/" + did + "/,\n"
                     "# which points read-only at the real data. Run from there:\n"
                     "#   cd .md_pipeline/mddb/work/" + did + " && mwf run -dir . -inp inputs.yaml\n"
                     "# Where each input really lives: provenance.yaml beside this file.\n\n")
            yaml.safe_dump(ordered, fh, sort_keys=False, default_flow_style=False,
                           allow_unicode=True)
    return atomic.atomic(path, _w)


def write_provenance(path, dataset_id, prov):
    """Where every input of an export really lives. Not an mwf input file."""
    import yaml
    def _w(tmp):
        with open(tmp, "w") as fh:
            fh.write(f"# provenance.yaml -- where the inputs of {dataset_id} live.\n"
                     "# Not an mwf input. inputs.yaml beside this file refers to links\n"
                     "# under work_dir; each link points at source_root/source_dir/<file>.\n"
                     "# To run on another host, change source_root and rebuild the links.\n\n")
            yaml.safe_dump(prov, fh, sort_keys=False, default_flow_style=False)
    return atomic.atomic(path, _w)
