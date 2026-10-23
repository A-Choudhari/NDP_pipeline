#!/usr/bin/env python3
"""
layout.py -- where everything a run writes goes. The only module that knows.

WHAT A RUN LEAVES BEHIND. Two directories, and a person only ever opens one:

    <folder>/MD_RESULTS/                 five files, all of them readable
        README.md                        what was found, what is missing, next steps
        settings.yaml                    THE file you edit
        glossary.md                      every MD, one row each
        NDP_registration.md              what to paste into the NDP form
        <name>.json                      the record to host and link
        NDP_records/<name>.json          ... when there are more than three
        standard_formats/                each record as schema.org/DCAT and
                                         RO-Crate, the folder as one DCAT
                                         catalogue (core/standard.py)

    <folder>/.md_pipeline/               machinery; hidden, and skipped by the scan
        index/     registry/*.parquet, tables/, LEGEND.md, notes.yaml
        mddb/      registry/*.parquet, exports/<ds>/, work/<ds>/, cache/, ndp_index.json
        logs/      pipeline-<stamp>.log, cron.log, status.json
        backups/   settings.before-<stamp>.yaml
        migrated-<date>/   the old SUBSTRATE/ and MDDB_substrate/ config files

WHY A HIDDEN DIRECTORY. The scan skips dot-directories (config.walk), so the
machinery cannot index itself, and `ls` in a simulation folder shows one new
name instead of two trees and forty-seven files.

THE OLD LAYOUT (SUBSTRATE/ and MDDB_substrate/) is still recognised, because a
colleague's folder may hold one until they re-run: `legacy()` locates it and
`migrate.py` moves it. Nothing here deletes anything.
"""

import os
import re

#: The two directory names, in one place.
RESULTS = "MD_RESULTS"
INTERNAL = ".md_pipeline"

#: What the previous layout called them.
LEGACY_INDEX = "SUBSTRATE"
LEGACY_MDDB = "MDDB_substrate"

#: More than this many records and they get their own directory, so the five
#: files a person reads stay visible.
RECORDS_INLINE = 6

_OWN = (RESULTS, INTERNAL, LEGACY_INDEX, LEGACY_MDDB)


#: The two directories inside .md_pipeline/ that scripts are handed directly.
_INSIDE = ("index", "mddb")


def base(folder):
    """The data folder, given any of its own output directories.

    Handed MD_RESULTS/, .md_pipeline/, .md_pipeline/mddb/, an old SUBSTRATE/ or
    the folder itself, every one of them answers the same folder -- which is
    what lets a script take `folder` from a person without asking which of its
    directories they meant.
    """
    p = os.path.abspath(str(folder).rstrip("/"))
    while True:
        name, parent = os.path.basename(p), os.path.dirname(p)
        if name in _OWN:
            p = parent
        elif name in _INSIDE and os.path.basename(parent) == INTERNAL:
            p = os.path.dirname(parent)
        else:
            return p


def results(folder):
    return os.path.join(base(folder), RESULTS)


def internal(folder):
    return os.path.join(base(folder), INTERNAL)


# --- what a person opens ---------------------------------------------------
def settings(folder):
    return os.path.join(results(folder), "settings.yaml")


def readme(folder):
    return os.path.join(results(folder), "README.md")


def glossary_md(folder):
    return os.path.join(results(folder), "glossary.md")


def registration_md(folder):
    return os.path.join(results(folder), "NDP_registration.md")


def records_dir(folder, n_records=1):
    """Where the NDP records go: beside the other files, or in their own
    directory once there are more than three."""
    return (results(folder) if n_records <= RECORDS_INLINE
            else os.path.join(results(folder), "NDP_records"))


# --- machinery -------------------------------------------------------------
def index_dir(folder):
    """The index: what the substrate layer writes (was <folder>/SUBSTRATE)."""
    return os.path.join(internal(folder), "index")


def mddb_dir(folder):
    """The deposition layer's own directory (was <folder>/MDDB_substrate)."""
    return os.path.join(internal(folder), "mddb")


def registry(folder):
    return os.path.join(index_dir(folder), "registry")


def tables(folder):
    return os.path.join(index_dir(folder), "tables")


def legend(folder):
    return os.path.join(index_dir(folder), "LEGEND.md")


def notes(folder):
    """Facts, caveats and document status a person records (was curated.yaml)."""
    return os.path.join(index_dir(folder), "notes.yaml")


def mddb_registry(folder):
    return os.path.join(mddb_dir(folder), "registry")


def exports(folder):
    return os.path.join(mddb_dir(folder), "exports")


def cache(folder):
    return os.path.join(mddb_dir(folder), "cache")


def ndp_index(folder):
    return os.path.join(mddb_dir(folder), "ndp_index.json")


def logs(folder):
    return os.path.join(internal(folder), "logs")


def status(folder):
    return os.path.join(logs(folder), "status.json")


def backups(folder):
    return os.path.join(internal(folder), "backups")


def lock(folder, name):
    return os.path.join(internal(folder), f".{name}.lock")


# --- the layout that came before -------------------------------------------
def legacy(folder):
    """(old index dir, old mddb dir) -- whether or not they exist."""
    b = base(folder)
    return os.path.join(b, LEGACY_INDEX), os.path.join(b, LEGACY_MDDB)


def legacy_present(folder):
    """True when a folder still holds the old layout and not the new one."""
    old_index, old_mddb = legacy(folder)
    if os.path.exists(settings(folder)):
        return False
    return (os.path.exists(os.path.join(old_index, "substrate.yaml"))
            or os.path.exists(os.path.join(old_mddb, "datasets.yaml")))


def excluded(folder):
    """Every directory a scan of this folder must not index: ours, old and new.

    Returned for `config.load`'s `exclude_paths`, so the index never describes
    its own output -- and so `MD_RESULTS/` is reported as *output* rather than
    as a directory nobody has covered.
    """
    b = base(folder)
    return [os.path.join(b, name) for name in _OWN]


_STAMP = re.compile(r"Generated \d{4}-\d\d-\d\dT[\d:.+\-]+")


def same_page(path, text):
    """True when `path` already says this, give or take the time it was made."""
    if not os.path.exists(path):
        return False
    with open(path) as fh:
        old = fh.read()
    return _STAMP.sub("Generated", old) == _STAMP.sub("Generated", text)


def write_page(path, text):
    """Write a generated page, but only when it really says something new.

    Every page carries the time it was made, so a comparison that counted the
    stamp would rewrite all of them on every run and make an unchanged folder
    look like a changed one to whoever syncs it (or reads `ls -lt`). The same
    rule the catalogue already applies to its `generated` field.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if same_page(path, text):
        return path
    from .atomic import write_text                           # noqa: PLC0415
    return write_text(path, text)


def glossary_parquet(folder):
    """The MD glossary, wherever this folder keeps it.

    Three layouts answer: the current one, the old MDDB_substrate/ one, and the
    lab-wide copy, which is a bare registry/ + reports/ pair.
    """
    b = base(folder)
    for path in (os.path.join(mddb_registry(b), "mds.parquet"),
                 os.path.join(b, LEGACY_MDDB, "registry", "mds.parquet"),
                 os.path.join(b, "registry", "mds.parquet")):
        if os.path.exists(path):
            return path
    return None
