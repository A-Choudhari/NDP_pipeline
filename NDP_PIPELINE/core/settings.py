#!/usr/bin/env python3
"""
settings.py -- MD_RESULTS/settings.yaml: the one file a person edits.

It replaces three: SUBSTRATE/substrate.yaml (what to index), MDDB_substrate/
datasets.yaml (the science behind each dataset) and MDDB_substrate/ndp.yaml
(what only a person knows about each simulated system). One file, three
sections, and the same rule everywhere: the tools fill in what they MEASURE
and never touch a value a person wrote.

    folder:                  what to index -- corpora, labels, stages
    production_lengths_ns:   how long a full production MD is here
    datasets:                one block per dataset: authors, paper, licence
    ndp:                     defaults, per-dataset values, per-system curation

Readers call `load(folder)`; `config_of(folder)` hands the `folder:` section to
the index layer in the shape lib/config.py has always used. Writers go through
`add_block` (append, for a dataset or system that has just appeared) and
`fill(...)` (fill empty fields only) so comments and edits survive; every write
takes a timestamped copy into .md_pipeline/backups/ first.
"""

import os
import re
import shutil
import time

from . import config, layout, yamledit

SECTIONS = ("folder", "production_lengths_ns", "datasets", "ndp")

HEADER = """# settings.yaml -- the one file to edit for this folder.
#
# Written by the pipeline, edited by you. Anything the tools can MEASURE they
# fill in; anything only you know stays yours and is never overwritten. Re-run
# the pipeline after editing.
#
#   folder:                 what is indexed: corpora, system labels, and what a
#                           trajectory file really is (`stages:`)
#   production_lengths_ns:  how long a full production MD is here. An MD that
#                           reaches one of these is a production MD of that
#                           length; a shorter one is flagged and left out of the
#                           NDP records.
#   datasets:               one block per dataset -- the science: authors, the
#                           paper, licence, PDB entries, water model
#   ndp:                    the catalogue record: title, construct, glycans,
#                           ions, and which stages and MDs it lists
#
# Everything the run produces is in MD_RESULTS/ beside this file; the machinery
# lives in .md_pipeline/ and can be deleted and rebuilt at any time.

"""


def path(folder):
    return layout.settings(folder)


def exists(folder):
    return os.path.exists(path(folder))


def read(folder):
    with open(path(folder)) as fh:
        return fh.read()


def load(folder):
    """The whole document, with every section present (empty when unset)."""
    import yaml
    if not exists(folder):
        raise FileNotFoundError(
            f"no settings.yaml at {path(folder)} -- run the pipeline on "
            f"{layout.base(folder)} first")
    doc = yaml.safe_load(read(folder)) or {}
    doc.setdefault("folder", {})
    doc.setdefault("datasets", [])
    doc.setdefault("ndp", {})
    return doc


class RootMoved(ValueError):
    """`folder: root:` no longer names the folder the settings file sits in."""


def config_of(folder, doc=None):
    """The `folder:` section as the index layer's cfg dict.

    `lib.config` has always taken a validated mapping with `root`, `corpora`,
    `prune`, the label rules and `exclude_paths`. It still does; only where the
    text comes from has changed.

    THE ROOT IS CHECKED, NOT TRUSTED. `root` was filled in only when missing
    and never compared to the folder again, so a settings.yaml that had been
    copied, or a folder that had been renamed or moved, made the build INDEX
    THE OLD TREE while writing records into the new folder -- silently, and
    with every path in the glossary pointing somewhere else. It becomes a live
    hazard the moment MD_RESULTS folders are being copied around, which is
    exactly what the collector does. Copying a settings.yaml as a starting
    point for a new folder is a reasonable thing to do; inheriting its root is
    never what was meant.
    """
    doc = doc if doc is not None else load(folder)
    cfg = dict(doc.get("folder") or {})
    here = layout.base(folder)
    cfg.setdefault("root", here)
    declared = os.path.realpath(os.path.expanduser(str(cfg["root"])))
    if declared != os.path.realpath(here):
        raise RootMoved(
            f"{path(folder)}: `folder: root:` says\n"
            f"      {declared}\n"
            f"  but this settings file is in\n"
            f"      {os.path.realpath(here)}\n"
            f"  The build would read the first and write into the second. If "
            f"this folder was moved or renamed, or these settings were copied "
            f"from another folder, re-run with --rebase-root to point it here. "
            f"If the root is right and the location is wrong, move the folder "
            f"back.")
    return config.validate(cfg, path(folder), folder=here)


def rebase_root(folder):
    """Point `folder: root:` at the folder this settings file is actually in.

    Rewrites one value and leaves every comment and every other key alone,
    through the same comment-preserving editor the rest of the file uses. The
    previous value is printed, because a root that was wrong may mean the
    records in this folder describe someone else's data.
    """
    doc = load(folder)
    here = os.path.realpath(layout.base(folder))
    was = ((doc.get("folder") or {}).get("root"))
    if was and os.path.realpath(os.path.expanduser(str(was))) == here:
        return was, here, False
    text = open(path(folder)).read()
    fixed = re.sub(r"^(\s*root:).*$", lambda m: f"{m.group(1)} {here}",
                   text, count=1, flags=re.M)
    if fixed == text:
        raise ValueError(f"{path(folder)}: no `root:` line to rebase")
    write(folder, fixed)          # backs up what was there, writes atomically
    return was, here, True


def datasets_of(folder, doc=None):
    doc = doc if doc is not None else load(folder)
    return list(doc.get("datasets") or [])


def ndp_of(folder, doc=None):
    doc = doc if doc is not None else load(folder)
    return dict(doc.get("ndp") or {})


def production_lengths(folder, doc=None):
    doc = doc if doc is not None else load(folder)
    lengths = doc.get("production_lengths_ns")
    return list(lengths) if lengths else None


# --- writing ---------------------------------------------------------------
def backup(folder):
    """A timestamped copy of settings.yaml, kept out of sight."""
    src = path(folder)
    if not os.path.exists(src):
        return None
    out = os.path.join(layout.backups(folder),
                       f"settings.before-{time.strftime('%Y%m%dT%H%M%S')}"
                       f"-{time.time_ns() % 10**9:09d}.yaml")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    shutil.copy2(src, out)
    return out


def write(folder, text):
    """Replace settings.yaml atomically, after backing up what was there."""
    backup(folder)
    target = path(folder)
    os.makedirs(os.path.dirname(target), exist_ok=True)
    from .atomic import write_text                           # noqa: PLC0415
    return write_text(target, text, newline="")  # line endings as given


LENGTHS_SECTION = """
# --- how long a full production MD is here -----------------------------------
# An MD that reaches one of these lengths (2 % tolerance) is a production MD of
# that length; one shorter than all of them is flagged `short_of_nominal` and
# left out of the NDP records. A dataset can declare its own list.

production_lengths_ns: [500, 1000, 1500]
"""


DATASETS_SECTION = """
# --- the science behind each dataset ----------------------------------------
# Proposed from the index: the selector and the counts were measured, every
# other field is yours. A field left empty is OMITTED from the deposition files
# rather than written as null -- an absent key means "not established", a null
# means "measured and found to be nothing". Only say what you know.
#
# Spec: https://mddb-workflow.readthedocs.io/en/latest/input.html
#
# Selector: `system:` reuses a label declared under `folder: systems:`, so new
# replicas are picked up as they are produced. `paths:` lists run directories
# when a dataset cannot be expressed as a label.

datasets:
"""


def section_head(name):
    """A section's explanatory comment and key.

    Both live in core/ so every pipeline that shares core/ writes the same
    file. They used to be imported from `mddb_init` and `mddb_lib.ndp` -- module
    names that no longer exist -- behind a bare except, so every new
    settings.yaml silently lost both explanations.
    """
    if name == "datasets":
        return DATASETS_SECTION
    if name == "ndp":
        from .ndp import SECTION_HEAD                        # noqa: PLC0415
        return SECTION_HEAD
    raise KeyError(name)


def create(folder, folder_section="", lengths_section=None, datasets_section=None,
           ndp_section=None):
    """Write a settings.yaml that has all four sections, in order."""
    from . import layout                                     # noqa: PLC0415
    if datasets_section is None:
        datasets_section = section_head("datasets")
    if ndp_section is None:
        ndp_section = section_head("ndp")
    text = (HEADER
            + (folder_section or "folder:\n")
            + (lengths_section if lengths_section is not None else LENGTHS_SECTION)
            + datasets_section
            + ndp_section)
    os.makedirs(layout.results(folder), exist_ok=True)
    return write(folder, text)


def add_block(folder, section, block):
    """Append a dataset or system block to its section, keeping the rest."""
    text = yamledit.insert_under(read(folder), section, block)
    return write(folder, text)


def fill(folder, updates, path_in_doc=("datasets",), key="id"):
    """Fill EMPTY fields of the named entries. Never replaces a written value."""
    text = yamledit.fill_empty(read(folder), updates, path_in_doc, key)
    return write(folder, text)
