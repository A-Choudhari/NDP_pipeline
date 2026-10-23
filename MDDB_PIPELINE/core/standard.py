#!/usr/bin/env python3
"""
standard.py -- every record again, in formats other catalogues read natively.

Our JSON is ours: its field names are this pipeline's. A catalogue that holds
every lab's records has to read them without learning each lab's keys, so each
record is also written in three standard shapes, built from the record alone:

    MD_RESULTS/standard_formats/
      catalog.dcat.jsonld            the folder as one DCAT catalogue
      <record>/dataset.jsonld        schema.org Dataset + DCAT Dataset, one node
      <record>/ro-crate-metadata.json  RO-Crate 1.1, Process Run Crate profile
      <record>/record.json           the record itself, so the crate is whole

WHY THESE. The National Data Platform runs on CKAN, whose DCAT extension reads
DCAT (and writes schema.org); Google Dataset Search reads schema.org. One
JSON-LD node typed as both carries both vocabularies, so either reader finds
its own terms. The DCAT catalogue is what a harvester points at: one file per
folder, every dataset in it. RO-Crate packages the record with its provenance
-- which engine produced which replica, and which tool read them -- as
CreateActions, which is what the Process Run Crate profile asks for.

NOTHING NEW IS SAID HERE. Every value is the record's; a value the record has
as null is left out, never filled. The vocabulary for "may this train an AI
model" is not settled anywhere, so it is written twice: as an ODRL policy
(W3C; `purpose` = AI model training) and as plain text in `usageInfo`.
"""

import json
import os
import re

#: Terms used. schema.org is the default vocabulary; DCAT, Dublin Core, FOAF,
#: vCard and ODRL are prefixed.
CONTEXT = {
    "@vocab": "https://schema.org/",
    "dcat": "http://www.w3.org/ns/dcat#",
    "dct": "http://purl.org/dc/terms/",
    "foaf": "http://xmlns.com/foaf/0.1/",
    "vcard": "http://www.w3.org/2006/vcard/ns#",
    "odrl": "http://www.w3.org/ns/odrl/2/",
    "xsd": "http://www.w3.org/2001/XMLSchema#",
}

#: access_level -> the EU access-right vocabulary DCAT-AP uses.
ACCESS_RIGHTS = {
    "public": "http://publications.europa.eu/resource/authority/access-right/PUBLIC",
    "licensable": "http://publications.europa.eu/resource/authority/access-right/RESTRICTED",
    "private": "http://publications.europa.eu/resource/authority/access-right/NON_PUBLIC",
}

ACCESS_TEXT = {
    "public": "Raw trajectories are public.",
    "licensable": "Raw trajectories are available under a licence agreement.",
    "private": "Raw trajectories are not public.",
}

AI_TEXT = {
    "allowed": "Use of this data to train AI/ML models is allowed.",
    "prohibited": "Use of this data to train AI/ML models is not allowed.",
    "by-agreement": "Use of this data to train AI/ML models requires the data owner's agreement.",
}

AI_PURPOSE = "AI model training"

RO_CRATE = "https://w3id.org/ro/crate/1.1"
PROCESS_RUN = "https://w3id.org/ro/wfrun/process/0.5"

DIR = "standard_formats"
CATALOG = "catalog.dcat.jsonld"
DATASET = "dataset.jsonld"
CRATE = "ro-crate-metadata.json"
RECORD = "record.json"
#: Everything this module writes into a record's directory. A directory
#: holding anything else is not ours, and is never removed.
OURS = {DATASET, CRATE, RECORD}


# --------------------------------------------------------------------------
# small pieces
# --------------------------------------------------------------------------
def _clean(node):
    """Drop None, empty strings and empty containers, all the way down."""
    if isinstance(node, dict):
        out = {k: _clean(v) for k, v in node.items()}
        return {k: v for k, v in out.items() if v not in (None, "", [], {})}
    if isinstance(node, list):
        out = [_clean(v) for v in node]
        return [v for v in out if v not in (None, "", [], {})]
    return node


def _slug(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(text)).strip("-.") or "x"


def _date(d):
    return {"@value": d, "@type": "xsd:date"} if d else None


def license_url(lic):
    """An SPDX id becomes its spdx.org URL; a URL stays; free text has none."""
    if not lic:
        return None
    lic = str(lic).strip()
    if re.match(r"^https?://", lic):
        return lic
    if re.match(r"^[A-Za-z0-9.+-]+$", lic):
        return f"https://spdx.org/licenses/{lic}"
    return None


def orcid_url(v):
    if not v:
        return None
    v = str(v).strip()
    m = re.search(r"(\d{4}-\d{4}-\d{4}-\d{3}[\dX])$", v)
    return f"https://orcid.org/{m.group(1)}" if m else None


def authors_of(record):
    """The record's author line back into names."""
    return [a.strip() for a in str(record.get("authors") or "").split(",") if a.strip()]


def software_of(record):
    """`NAMD 2.14 (Frontera, TACC); AMBER 22` -> [(name, version, note)]."""
    out = []
    for s in str((record.get("setup") or {}).get("software") or "").split(";"):
        m = re.match(r"^\s*(\S+)(?:\s+(\d[\w.\-]*))?\s*(.*?)\s*$", s)
        if m and m.group(1):
            out.append((m.group(1), m.group(2), m.group(3).strip("() ") or None))
    return out


def _one(refs):
    """RO-Crate writes a single value, not a list of one."""
    refs = [r for r in refs if r]
    return refs[0] if len(refs) == 1 else refs


def record_url(entry, base_url):
    if not base_url or not entry.get("file"):
        return None
    return f"{str(base_url).rstrip('/')}/{os.path.basename(entry['file'])}"


def dataset_id(record, entry, base_url):
    """The dataset's IRI: the hosted record's URL once `base_url` is set,
    otherwise a fragment named for the record's file -- unique in a folder,
    which a person's identifier is not guaranteed to be."""
    url = record_url(entry, base_url)
    if url:
        return url + "#dataset"
    name = (os.path.basename(entry["file"])[:-5] if entry.get("file") else
            (record.get("lab_glossary") or {}).get("record_id") or "dataset")
    return "#" + _slug(name)


EMAIL = re.compile(r"^[^@\s,;]+@[^@\s,;]+\.[^@\s,;]+$")


def mailto(addr):
    """`mailto:` for one plain address; None for anything else (a list of
    several, a name with an address), which no IRI can hold."""
    addr = str(addr or "").strip()
    return f"mailto:{addr}" if EMAIL.match(addr) else None


def ai_policy(record, target):
    """The AI-training answer as an ODRL policy; None when unanswered."""
    ai = record.get("ai_training")
    if ai not in AI_TEXT:
        return None
    rule = {"odrl:action": {"@id": "odrl:use"},
            "odrl:target": {"@id": target},
            "odrl:constraint": {"odrl:leftOperand": {"@id": "odrl:purpose"},
                                "odrl:operator": {"@id": "odrl:eq"},
                                "odrl:rightOperand": AI_PURPOSE}}
    if ai == "prohibited":
        body = {"odrl:prohibition": [rule]}
    else:
        if ai == "by-agreement":
            rule["odrl:duty"] = {"odrl:action": {"@id": "odrl:obtainConsent"}}
        body = {"odrl:permission": [rule]}
    return {"@id": f"{target}-ai-training-policy", "@type": "odrl:Set",
            "odrl:uid": {"@id": f"{target}-ai-training-policy"}, **body}


def _owner_person(owner):
    owner = owner or {}
    if not owner.get("name"):
        return None
    return {"@type": "Person", "name": owner["name"],
            "identifier": orcid_url(owner.get("orcid")),
            "email": owner.get("email"),
            "affiliation": ({"@type": "Organization", "name": owner["affiliation"]}
                            if owner.get("affiliation") else None)}


def _raw_description(record):
    loc = record.get("data_location") or {}
    where = ": ".join(x for x in (loc.get("storage"), loc.get("path")) if x)
    parts = [f"Raw trajectory files ({record['files']} files)" if record.get("files")
             else "Raw trajectory files",
             f"on {where}" if where else None]
    text = " ".join(x for x in parts if x) + "."
    avail = str(record.get("data_availability") or "").strip().rstrip(".")
    return text + (f" {avail}." if avail else "")


# --------------------------------------------------------------------------
# schema.org + DCAT
# --------------------------------------------------------------------------
def dataset(record, entry=None, base_url=None, lab_info=None):
    """One node typed schema:Dataset and dcat:Dataset, without @context."""
    entry, lab_info = entry or {}, lab_info or {}
    did = dataset_id(record, entry, base_url)
    url = record_url(entry, base_url)
    lic = license_url(record.get("license"))
    owner = record.get("data_owner") or {}
    access = record.get("access_level")
    keywords = list(entry.get("keywords") or []) or list(record.get("tags") or [])
    authors = authors_of(record)
    loc = record.get("data_location") or {}
    size_gb = record.get("size_gb")
    publisher = record.get("lab") or lab_info.get("name")
    contact = record.get("contact")

    # Each distribution is one node, listed under both vocabularies' keys.
    record_dist = {
        "@id": f"{did}-record",
        "@type": ["DataDownload", "dcat:Distribution"],
        "name": "MD metadata record (JSON)", "dct:title": "MD metadata record (JSON)",
        "encodingFormat": "application/json", "dcat:mediaType": "application/json",
        "contentUrl": url, "dcat:downloadURL": {"@id": url} if url else None,
        "dcat:accessURL": {"@id": url} if url else None,
    }
    raw_dist = {
        "@id": f"{did}-raw",
        "@type": ["DataDownload", "dcat:Distribution"],
        "name": "Raw trajectories", "dct:title": "Raw trajectories",
        "description": _raw_description(record),
        "dct:description": _raw_description(record),
        "contentSize": f"{size_gb} GB" if size_gb is not None else None,
        "dcat:byteSize": ({"@value": str(int(round(size_gb * 1e9))),
                           "@type": "xsd:nonNegativeInteger"}
                          if size_gb is not None else None),
        "contentUrl": loc.get("uri"),
        "dcat:accessURL": {"@id": loc["uri"]} if loc.get("uri") else None,
        "dct:accessRights": {"@id": ACCESS_RIGHTS[access]} if access in ACCESS_RIGHTS else None,
    }
    node = {
        "@id": did,
        "@type": ["Dataset", "dcat:Dataset"],
        # --- schema.org --------------------------------------------------
        "identifier": record.get("identifier"),
        "name": record.get("title"),
        "description": record.get("description"),
        "keywords": keywords,
        "about": ({"@type": "DefinedTerm", "name": record["topic"]}
                  if record.get("topic") else None),
        "measurementTechnique": "Molecular dynamics simulation",
        "url": url,
        "creator": [{"@type": "Person", "name": a} for a in authors],
        "publisher": {"@type": "Organization", "name": publisher} if publisher else None,
        "accountablePerson": _owner_person(owner),
        "contactPoint": ({"@type": "ContactPoint", "email": contact,
                          "contactType": "data requests"} if contact else None),
        "license": lic or record.get("license"),
        "citation": record.get("related_publication"),
        "isAccessibleForFree": (access == "public") if access in ACCESS_RIGHTS else None,
        "conditionsOfAccess": ACCESS_TEXT.get(access),
        "usageInfo": ({"@type": "CreativeWork", "name": "AI/ML training",
                       "text": AI_TEXT[record["ai_training"]]}
                      if record.get("ai_training") in AI_TEXT else None),
        "contentSize": f"{size_gb} GB" if size_gb is not None else None,
        "dateCreated": entry.get("first_written"),
        "dateModified": entry.get("last_changed"),
        "isPartOf": ({"@type": "DataCatalog", "name": lab_info.get("catalog_name"),
                      "url": lab_info.get("catalog_url")}
                     if lab_info.get("catalog_name") else None),
        "distribution": [record_dist, raw_dist],
        # --- DCAT / Dublin Core --------------------------------------------
        "dct:identifier": record.get("identifier"),
        "dct:title": record.get("title"),
        "dct:description": record.get("description"),
        "dcat:keyword": keywords,
        "dct:creator": [{"@type": "foaf:Person", "foaf:name": a} for a in authors],
        "dct:publisher": ({"@type": "foaf:Organization", "foaf:name": publisher}
                          if publisher else None),
        "dct:rightsHolder": ({"@type": "foaf:Person", "foaf:name": owner["name"],
                              "foaf:mbox": ({"@id": mailto(owner.get("email"))}
                                            if mailto(owner.get("email")) else None)}
                             if owner.get("name") else None),
        "dcat:contactPoint": ({"@type": "vcard:Kind",
                               "vcard:fn": owner.get("name") or publisher,
                               "vcard:hasEmail": ({"@id": mailto(contact)}
                                                  if mailto(contact) else contact)}
                              if contact else None),
        "dct:license": {"@id": lic} if lic else None,
        "dct:accessRights": {"@id": ACCESS_RIGHTS[access]} if access in ACCESS_RIGHTS else None,
        "dct:bibliographicCitation": record.get("related_publication"),
        "dct:issued": _date(entry.get("first_written")),
        "dct:modified": _date(entry.get("last_changed")),
        "dcat:landingPage": {"@id": url} if url else None,
        "dcat:distribution": [{"@id": record_dist["@id"]}, {"@id": raw_dist["@id"]}],
        "odrl:hasPolicy": ai_policy(record, did),
    }
    return _clean(node)


def dataset_doc(record, entry=None, base_url=None, lab_info=None):
    return {"@context": CONTEXT, **dataset(record, entry, base_url, lab_info)}


def catalog_doc(folder, nodes, lab_info=None, base_url=None):
    """The folder's datasets as one dcat:Catalog (and schema:DataCatalog)."""
    lab_info = lab_info or {}
    name = os.path.basename(str(folder).rstrip("/"))
    title = " — ".join(x for x in (lab_info.get("catalog_name") or lab_info.get("name"), name) if x)
    cid = f"#catalog-{_slug(name)}"
    return {"@context": CONTEXT, **_clean({
        "@id": cid,
        "@type": ["DataCatalog", "dcat:Catalog"],
        "name": title, "dct:title": title,
        "url": lab_info.get("catalog_url"),
        "dcat:landingPage": ({"@id": lab_info["catalog_url"]}
                             if lab_info.get("catalog_url") else None),
        "publisher": ({"@type": "Organization", "name": lab_info["name"]}
                      if lab_info.get("name") else None),
        "dct:publisher": ({"@type": "foaf:Organization", "foaf:name": lab_info["name"]}
                          if lab_info.get("name") else None),
        "dataset": [{"@id": n["@id"]} for n in nodes],
        "dcat:dataset": nodes,
    })}


# --------------------------------------------------------------------------
# RO-Crate
# --------------------------------------------------------------------------
def crate(record, entry=None, base_url=None, lab_info=None, record_bytes=None):
    """ro-crate-metadata.json for one record's directory."""
    entry, lab_info = entry or {}, lab_info or {}
    did = dataset_id(record, entry, base_url)
    lic = license_url(record.get("license"))
    owner = record.get("data_owner") or {}
    loc = record.get("data_location") or {}
    access = record.get("access_level")
    sampling = record.get("sampling") or {}
    graph = []

    def add(node):
        node = _clean(node)
        if node and not any(n["@id"] == node["@id"] for n in graph):
            graph.append(node)
        return {"@id": node["@id"]} if node else None

    people = []
    for name in authors_of(record):
        people.append(add({"@id": f"#person-{_slug(name)}", "@type": "Person", "name": name}))
    owner_ref = None
    if owner.get("name"):
        owner_ref = add({"@id": orcid_url(owner.get("orcid")) or f"#owner-{_slug(owner['name'])}",
                         "@type": "Person", "name": owner["name"],
                         "email": owner.get("email"),
                         "affiliation": owner.get("affiliation")})
    publisher = record.get("lab") or lab_info.get("name")
    org_ref = add({"@id": lab_info.get("url") or f"#org-{_slug(publisher)}",
                   "@type": "Organization", "name": publisher,
                   "url": lab_info.get("url")}) if publisher else None
    contact_ref = add({"@id": mailto(record["contact"]) or "#contact", "@type": "ContactPoint",
                       "email": record["contact"], "contactType": "data requests"}) \
        if record.get("contact") else None
    lic_ref = add({"@id": lic, "@type": "CreativeWork", "name": record.get("license")}) \
        if lic else None
    profile_ref = add({"@id": PROCESS_RUN, "@type": "CreativeWork",
                       "name": "Process Run Crate", "version": "0.5"})

    engines = []
    for name, version, note in software_of(record):
        engines.append(add({"@id": f"#software-{_slug(name)}-{_slug(version or 'unknown')}",
                            "@type": "SoftwareApplication", "name": name,
                            "softwareVersion": version, "description": note}))
    tool_ref = add({"@id": "#ndp-pipeline", "@type": "SoftwareApplication",
                    "name": "NDP_PIPELINE (MD publishing pipeline)",
                    "description": "Reads trajectory headers, MD configs and run logs "
                                   "and writes one metadata record per simulated system."})

    # One data entity and one CreateAction per replica: the engine produced
    # the trajectory files; where they are is the record's measured path.
    raw, runs, used = [], [], set()
    for d in sampling.get("replicas_detail") or []:
        # Unique per replica: two ids that slug alike (`r 1`, `r/1`) must
        # not become one entity, which `add` would silently keep once.
        rid, n = _slug(d.get("id")), 2
        while rid in used:
            rid, n = f"{_slug(d.get('id'))}-{n}", n + 1
        used.add(rid)
        where = d.get("location") or d.get("path")
        url = None
        if loc.get("uri") and d.get("location") and loc.get("path"):
            rel = os.path.relpath(d["location"], loc["path"])
            url = str(loc["uri"]).rstrip("/") + ("" if rel == "." else "/" + rel)
        raw_ref = add({"@id": url or f"#raw-{rid}", "@type": "Dataset",
                       "name": f"Raw trajectory files, replica {d.get('id')}",
                       "description": ": ".join(x for x in (loc.get("storage"), where) if x)
                       or None,
                       "contentUrl": url,
                       "conditionsOfAccess": ACCESS_TEXT.get(access)})
        raw.append(raw_ref)
        desc = [f"{d['length_ns']} ns" if d.get("length_ns") is not None else None,
                f"{d['frames']} frames" if d.get("frames") else None,
                f"planned {d['planned_ns']} ns" if d.get("planned_ns") else None,
                f"{d['gap_ns']} ns missing inside the run" if d.get("gap_ns") else None]
        runs.append(add({"@id": f"#run-{rid}", "@type": "CreateAction",
                         "name": f"Production MD, replica {d.get('id')}",
                         "description": ", ".join(x for x in desc if x) or None,
                         "instrument": _one(engines), "result": raw_ref,
                         "actionStatus": ({"@id": "http://schema.org/CompletedActionStatus"}
                                          if d.get("complete") else None)}))
    files = [add({"@id": RECORD, "@type": "File", "name": "MD metadata record",
                  "encodingFormat": "application/json",
                  "contentSize": str(record_bytes) if record_bytes else None}),
             add({"@id": DATASET, "@type": "File",
                  "name": "schema.org / DCAT description",
                  "encodingFormat": "application/ld+json"})]
    runs.append(add({"@id": "#metadata-extraction", "@type": "CreateAction",
                     "name": "Metadata read from the trajectories, configs and logs",
                     "description": "Every value in the record was read from the "
                                    "replicas' own files or declared in settings.yaml; "
                                    "the record says which.",
                     "instrument": tool_ref, "object": _one(raw), "result": files[0],
                     "endTime": entry.get("last_changed")}))

    root = _clean({
        "@id": "./", "@type": "Dataset",
        "identifier": record.get("identifier") or did,
        "name": record.get("title") or (record.get("lab_glossary") or {}).get("record_id"),
        "description": record.get("description"),
        "keywords": list(entry.get("keywords") or []) or list(record.get("tags") or []),
        "datePublished": entry.get("first_written"),
        "dateModified": entry.get("last_changed"),
        "license": lic_ref,
        "author": people, "publisher": org_ref, "accountablePerson": owner_ref,
        "contactPoint": contact_ref,
        "citation": record.get("related_publication"),
        "isAccessibleForFree": (access == "public") if access in ACCESS_RIGHTS else None,
        "conditionsOfAccess": ACCESS_TEXT.get(access),
        "usageInfo": AI_TEXT.get(record.get("ai_training")),
        "conformsTo": profile_ref,
        "hasPart": files,
        "mentions": runs,
        "sameAs": {"@id": did} if did.startswith("http") else None,
    })
    descriptor = {"@id": CRATE, "@type": "CreativeWork",
                  "conformsTo": {"@id": RO_CRATE}, "about": {"@id": "./"}}
    return {"@context": "https://w3id.org/ro/crate/1.1/context",
            "@graph": [descriptor, root] + graph}


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def dumps(doc):
    return json.dumps(doc, indent=2, ensure_ascii=False) + "\n"


def _write(path, text):
    """Only when it says something new: an unchanged build changes no file."""
    try:
        with open(path) as fh:
            if fh.read() == text:
                return False
    except OSError:
        pass
    from .atomic import write_text                           # noqa: PLC0415
    write_text(path, text)
    return True


def stem(entry):
    return os.path.basename(entry["file"])[:-5]


def write_all(results_dir, folder, pairs, base_url=None, lab_info=None, sweep=False):
    """Write every record's standard files and the folder's catalogue.

    `pairs` is [(index entry, record, record text)]. With `sweep`, a record
    directory under standard_formats/ that no record owns any more is removed
    -- only when it holds nothing but what this module writes. Returns the
    names of removed directories.
    """
    top = os.path.join(results_dir, DIR)
    wanted = set()
    for entry, record, text in pairs:
        d = os.path.join(top, stem(entry))
        wanted.add(stem(entry))
        _write(os.path.join(d, RECORD), text)
        _write(os.path.join(d, DATASET),
               dumps(dataset_doc(record, entry, base_url, lab_info)))
        _write(os.path.join(d, CRATE),
               dumps(crate(record, entry, base_url, lab_info,
                           record_bytes=len(text.encode("utf-8")))))
    removed = []
    if sweep and os.path.isdir(top):
        for name in sorted(os.listdir(top)):
            d = os.path.join(top, name)
            if name in wanted or not os.path.isdir(d) or os.path.islink(d):
                continue
            inside = set(os.listdir(d))
            if not inside <= OURS:
                continue
            for f in inside:
                os.remove(os.path.join(d, f))
            os.rmdir(d)
            removed.append(name)
    # The catalogue is every record directory present, so a partial build
    # (--datasets) still leaves the whole folder in it.
    nodes = []
    if os.path.isdir(top):
        for name in sorted(os.listdir(top)):
            p = os.path.join(top, name, DATASET)
            try:
                with open(p) as fh:
                    doc = json.load(fh)
            except (OSError, ValueError):
                continue
            doc.pop("@context", None)
            nodes.append(doc)
    _write(os.path.join(top, CATALOG),
           dumps(catalog_doc(folder, nodes, lab_info, base_url)))
    return removed
