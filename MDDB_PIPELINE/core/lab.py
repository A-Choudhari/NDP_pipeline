#!/usr/bin/env python3
"""
lab.py -- the lab this kit is set up for, read from ONE file: <kit>/lab.yaml.

Everything that names a lab used to sit in the code: the lab's name and
keywords in ndp.py, the identifier prefix and the topic vocabulary in
catalog.py, the directory-name topics in topics.py, one person's authors and
licence in datasets_init.py. Another lab had to find and edit all of them. Now
they are read from lab.yaml at the top of the kit, and nothing else.

    MD_PIPELINE_LAB=/path/to/other_lab.yaml     use that file instead

WITHOUT THE FILE, NOTHING IS INVENTED. Every lab-specific value is empty: no
lab name, no topics, the identifier prefix `LAB`. A record then lists what is
missing rather than carrying another lab's name.
"""

import copy
import os

#: Environment variable naming a lab file other than the kit's own.
ENV = "MD_PIPELINE_LAB"

#: NDP_PIPELINE/core/lab.py (or MDDB_PIPELINE/core/lab.py) -> the kit.
KIT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

#: The shape of lab.yaml, with every lab-specific value empty.
NEUTRAL = {
    "lab": {"name": None, "short_name": None, "url": None,
            "catalog_name": None, "catalog_url": None, "ndp_organization": None},
    "identifiers": {"prefix": "LAB"},
    "record_defaults": {
        "keywords": ["molecular dynamics"],
        "data_availability": None,
        "data_owner": {"name": None, "email": None, "orcid": None, "affiliation": None},
        "access_level": None,
        "ai_training": None,
        "data_location": {"storage": None, "uri": None},
    },
    "people": {},
    "topics": {},
    "index_topics": {},
    "storage": {"data_roots": []},
}

_cache = {}


def path():
    """The lab file in use: $MD_PIPELINE_LAB, else <kit>/lab.yaml."""
    return os.path.expanduser(os.environ.get(ENV) or os.path.join(KIT, "lab.yaml"))


#: Sections that are lists of the lab's own entries, replaced whole, never
#: merged key by key into the neutral (empty) ones.
_WHOLE = ("people", "topics", "index_topics")


def _merge(base, over):
    """`over` on top of `base`, one level of mapping at a time. A key the lab
    file leaves out keeps its neutral value; a key it sets empty is empty."""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(out.get(k), dict) and v is None:
            # `lab:` with nothing (or only comments) under it: an empty
            # section, not a section replaced by nothing.
            continue
        if isinstance(out.get(k), dict) and isinstance(v, dict) and k not in _WHOLE:
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def _list(v):
    """A YAML list as written; one plain value is a list of one, never a
    string to be iterated letter by letter."""
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [x for x in v if x is not None]
    return [v]


def load():
    """The whole lab file over NEUTRAL. Cached per path; a broken file is an
    error, because silently falling back would publish records under no lab."""
    p = path()
    if p not in _cache:
        doc = {}
        if os.path.exists(p):
            import yaml
            with open(p) as fh:
                doc = yaml.safe_load(fh) or {}
            if not isinstance(doc, dict):
                raise ValueError(f"{p}: expected a mapping at the top level")
        elif os.environ.get(ENV):
            raise FileNotFoundError(f"{ENV}={p}: no such file")
        got = _merge(NEUTRAL, doc)
        for section in ("lab", "identifiers", "record_defaults", "storage"):
            if not isinstance(got.get(section), dict):
                raise ValueError(f"{p}: `{section}:` must be a mapping")
        got["record_defaults"]["keywords"] = [
            str(k) for k in _list(got["record_defaults"].get("keywords"))]
        got["storage"]["data_roots"] = [str(r) for r in _list(got["storage"].get("data_roots"))]
        _cache[p] = got
    return copy.deepcopy(_cache[p])


def reload():
    _cache.clear()
    return load()


def _text(v):
    return str(v).strip() if v is not None and str(v).strip() else None


# --------------------------------------------------------------------------
# what the rest of the kit asks for
# --------------------------------------------------------------------------
def name():
    return _text(load()["lab"].get("name"))


def short_name():
    lab = load()["lab"]
    return _text(lab.get("short_name")) or _text(lab.get("name"))


def info():
    """The `lab:` block: name, short_name, url, catalog_name, catalog_url,
    ndp_organization (each None when unset)."""
    return {k: _text(v) for k, v in load()["lab"].items()}


def id_prefix():
    """First word of suggested identifiers. NDP_ID_PREFIX still wins."""
    got = os.environ.get("NDP_ID_PREFIX") or load()["identifiers"].get("prefix") or "LAB"
    return str(got).strip().upper() or "LAB"


def record_defaults():
    return load()["record_defaults"]


def people():
    """{account: {authors, contact, license, license_url}}."""
    return {str(k): (v if isinstance(v, dict) else {})
            for k, v in (load().get("people") or {}).items()}


def topics():
    """(word -> topic, topic -> code, topic -> molecule) for the catalogue."""
    words, codes, molecules = {}, {}, {}
    for slug, spec in (load().get("topics") or {}).items():
        slug, spec = str(slug), (spec or {})
        if not isinstance(spec, dict):
            spec = {}
        for w in [slug] + _list(spec.get("words")):
            words.setdefault(str(w).lower(), slug)
        if spec.get("code"):
            codes[slug] = str(spec["code"])
        if spec.get("molecule"):
            molecules[slug] = str(spec["molecule"])
    return words, codes, molecules


def index_topics():
    """{research area: [regex]} that init.py proposes from directory names."""
    return {str(k): [str(p) for p in _list(v)]
            for k, v in (load().get("index_topics") or {}).items()}


def data_roots():
    return [str(r) for r in (load().get("storage") or {}).get("data_roots") or []]


if __name__ == "__main__":
    import json
    print(f"# {path()}")
    print(json.dumps(load(), indent=2, ensure_ascii=False))
