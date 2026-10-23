#!/usr/bin/env python3
"""
catalog.py -- the lab catalogue's own view of a folder's records.

The lab's catalogue (lab.yaml `catalog_url`) lists every dataset under a
research topic, with a title, a description, tags, its volume, its file count
and a stable identifier (`SMITHLAB-HIV-ENV-001`). A record carries all of those, and
`catalog.json` gathers them so the site can ingest a folder in one read instead
of opening every record.

The identifier and the topic are a PERSON's to decide -- they have to be unique
across the lab, and the topic decides where the dataset is filed. What this
module does is propose them from the folder's own names, so the first run leaves
something to correct rather than a blank.
"""

import json
import os
import re
import time

from . import lab

#: The catalogue's topics, read from lab.yaml (`topics:`):
#:   TOPICS          folder or topic word -> topic slug, as the site uses them
#:   TOPIC_CODE      topic slug -> the identifier's topic segment
#:   TOPIC_MOLECULE  what a title names after the topic, where the topic is a
#:                   PATHOGEN and not the molecule that was simulated. `p53` or
#:                   `apobec3b` already name it, and repeating the word there
#:                   would only make the title longer.
TOPICS, TOPIC_CODE, TOPIC_MOLECULE = lab.topics()

#: First word of every suggested identifier: lab.yaml `identifiers: prefix:`,
#: or NDP_ID_PREFIX.
ID_PREFIX = lab.id_prefix()

#: Conformational states a directory name can carry, spelled as a title spells
#: them. What distinguishes two records of one construct is usually the state,
#: and it is usually in the path (`ABZYME/4-OPEN/CDR`).
STATES = {"open": "Open State", "closed": "Closed State",
          "occluded": "Occluded State", "prefusion": "Prefusion State",
          "postfusion": "Postfusion State"}


def _words(text):
    return [w for w in re.split(r"[^A-Za-z0-9]+", str(text)) if w]


def _human_ns(ns):
    """A duration a person reads: ns below a microsecond, us above it.

    The same rule as `human_ns` in MD_ANALYSES `core/io.py` -- copied, not
    imported, because the two trees are separate and stay that way; they must
    agree, because a record and the page describing it must never state
    different lengths for one run. Tags are ASCII, so this copy writes `us`
    where the page writes the µ sign.
    """
    if abs(ns) >= 1000:
        return f"{ns / 1000:.1f}".rstrip("0").rstrip(".") + " us"
    return f"{ns:.1f}".rstrip("0").rstrip(".") + " ns"


def topic_of(folder, cfg=None, record_id=""):
    """The catalogue topic this folder most likely belongs to, or None.

    Read from the folder's own path and the topics the index proposed. A guess
    a person can see and correct in settings.yaml; never written into a record
    without passing through there.
    """
    hay = " ".join([str(folder), " ".join((cfg or {}).get("topics") or []),
                    str(record_id)]).lower()
    for word in _words(hay):
        hit = TOPICS.get(word.lower())
        if hit:
            return hit
    return None


def suggest_identifier(topic, folder, record_id, taken=(), what=""):
    """`SMITHLAB-HIV-CDR-1F4-001` -- topic, where it lives, what it is, a number.

    `what` is the record's own sub-path (`1-1F4`), which is what tells one
    system in a folder from another; the topology hash never appears in a name
    a person has to read. Unique within what this run knows; the lab's
    catalogue is the real registry, so a person confirms it before listing.
    """
    code = TOPIC_CODE.get(topic or "", (topic or "LAB").upper()[:6])
    parts = [w.upper() for w in _words(os.path.basename(str(folder).rstrip("/")))]
    where = parts[0][:12] if parts else "MD"
    words = [w.upper() for w in _words(what) if not w.isdigit()]
    stem = "-".join(x for x in (ID_PREFIX, code, where, (words or [""])[0][:12]) if x)
    n = 1
    while f"{stem}-{n:03d}" in set(taken):
        n += 1
    return f"{stem}-{n:03d}"


def suggest_title(record, folder, curated=None, what=""):
    """`HIV - Env Open State - Tribound 1F4 antibody` -- what, in what state,
    with what bound to it.

    The owner's rule, 2026-09-24. Three clauses in a fixed order, because the
    catalogue lists a whole project's datasets one under the other: read down
    the column and the words that repeat are the ones the systems share, and
    the only word that changes is the one the science changes. The last clause
    is dropped entirely when nothing is bound -- an apo run that advertised a
    "Tribound  antibody" with a hole in it would be worse than no clause at
    all, and it is the apo run every comparison is measured against. A
    proposal, like the identifier: a person confirms it in settings.yaml.
    """
    curated = curated or {}
    topic = record.get("topic") or ""
    construct = str(curated.get("construct")
                    or (record.get("system") or {}).get("construct")
                    or TOPIC_MOLECULE.get(topic) or "")
    state = ""
    for part in reversed(str(folder).rstrip("/").split("/")):
        hit = [STATES[w.lower()] for w in _words(part) if w.lower() in STATES]
        if hit:
            state = "" if hit[0].split()[0].lower() in construct.lower() else hit[0]
            break
    # A one- or two-letter word is a folder's bookkeeping (`MD`, `1`), never
    # the name of what was bound -- which is how the apo system stays apo.
    ligand = next((w for w in _words(what) if len(w) > 2), "")
    return " - ".join(x for x in (
        TOPIC_CODE.get(topic, topic.upper()),
        " ".join(y for y in (construct, state) if y),
        f"Tribound {ligand} antibody" if ligand else "") if x)


def suggest_tags(record, curated=None, what=""):
    """A few honest tags: where it lives, the construct's words, the scale."""
    tags = list((curated or {}).get("tags") or [])
    system = (record.get("system") or {})
    for source in (what, system.get("construct"), record.get("title")):
        for w in _words(source or "")[:3]:
            if len(w) > 2 and w.lower() not in [t.lower() for t in tags]:
                tags.append(w)
    total = (record.get("sampling") or {}).get("total_ns")
    if total and total >= 1000:
        tags.append(_human_ns(total))
    return tags[:6]


def entry(record):
    """One row of catalog.json: exactly what the catalogue page shows."""
    sampling = record.get("sampling") or {}
    return {
        "identifier": record.get("identifier"),
        "topic": record.get("topic"),
        "tags": record.get("tags") or [],
        "title": record.get("title"),
        "description": record.get("description"),
        "size_gb": record.get("size_gb"),
        "files": record.get("files"),
        "replicas": sampling.get("replicas"),
        "complete_replicas": sampling.get("complete_replicas"),
        "total_ns": sampling.get("total_ns"),
        "record": record.get("record"),
        "data_availability": record.get("data_availability"),
        "license": record.get("license"),
        # Federation and licensing: who owns it, who may get it, and whether
        # it may train a model -- the last separate from the citation licence.
        "data_owner": (record.get("data_owner") or {}).get("name"),
        "access_level": record.get("access_level"),
        "ai_training": record.get("ai_training"),
    }


def write(path, folder, records, base_url=None):
    """catalog.json: every record of this folder, for the site to ingest."""
    doc = {
        "catalog": lab.info()["catalog_name"] or "MD data catalog",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "folder": str(folder),
        "base_url": base_url,
        "datasets": [entry(r) for r in records],
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    text = json.dumps(doc, indent=2, ensure_ascii=False) + "\n"
    if os.path.exists(path):
        with open(path) as fh:
            # Only `generated` differs on a run that changed nothing; rewriting
            # the file then would make every run look like a change.
            old = fh.read()
        if _without_time(old) == _without_time(text):
            return path
    from .atomic import write_text                           # noqa: PLC0415
    return write_text(path, text)


def _without_time(text):
    return re.sub(r'"generated": "[^"]*"', '"generated": ""', text)
