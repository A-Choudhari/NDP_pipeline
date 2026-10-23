#!/usr/bin/env python3
"""
migrate.py -- move a folder from the old layout to MD_RESULTS/ + .md_pipeline/.

The layout before 2026-09-22 put two visible trees in every data folder,
`SUBSTRATE/` and `MDDB_substrate/`, and asked a person to edit three separate
YAML files. This moves such a folder, once, keeping every value verbatim:

    SUBSTRATE/substrate.yaml      -> MD_RESULTS/settings.yaml, `folder:`
    MDDB_substrate/datasets.yaml  -> MD_RESULTS/settings.yaml, `datasets:`
    MDDB_substrate/ndp.yaml       -> MD_RESULTS/settings.yaml, `ndp:`
    SUBSTRATE/curated.yaml        -> .md_pipeline/index/notes.yaml
    SUBSTRATE/registry, tables    -> .md_pipeline/index/
    MDDB_substrate/registry,
      exports, work, cache        -> .md_pipeline/mddb/
    MDDB_substrate/reports/
      pipeline_status.json        -> .md_pipeline/logs/status.json

RULES. Nothing is deleted and nothing is rewritten: the three configs are
copied into `.md_pipeline/migrated-<date>/` before anything else, the
directories are MOVED (`os.rename` -- `work/` is a link tree and `tables/` can
be gigabytes), and the old directories are left in place holding a `MOVED.md`
that says where everything went. Running it twice does nothing the second time.
"""

import os
import time

from . import layout, settings

MOVED_NAME = "MOVED.md"


def needed(folder):
    """True when this folder still holds the old layout."""
    return layout.legacy_present(folder)


def _indent(text, prefix="  "):
    return "".join(f"{prefix}{line}" if line.strip() else line
                   for line in text.splitlines(keepends=True))


def _strip_leading_comments(text):
    """Drop a file's own header comment block, keep everything it said after."""
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith("#")):
        i += 1
    return "".join(lines[i:])


def _section_text(name, default):
    """A section's explanatory comment, from core/settings.py."""
    try:
        return settings.section_head(name)
    except KeyError:
        return default


def build_settings_text(substrate_yaml=None, datasets_yaml=None, ndp_yaml=None):
    """The new settings.yaml, from whichever old files exist.

    Values are carried across as TEXT, so a person's comments, ordering and
    spacing survive the move; only the two-space indent of a section changes.
    """
    import yaml
    folder_section = "folder:\n"
    if substrate_yaml:
        body = _strip_leading_comments(substrate_yaml)
        folder_section = "folder:\n" + _indent(body)
    lengths = settings.LENGTHS_SECTION
    datasets_section = _section_text("datasets", "\ndatasets:\n")
    if datasets_yaml:
        doc = yaml.safe_load(datasets_yaml) or {}
        if doc.get("production_lengths_ns"):
            lengths = ("\n# --- how long a full production MD is here "
                       "-----------------------------------\n\n"
                       f"production_lengths_ns: {doc['production_lengths_ns']}\n")
        body = _strip_leading_comments(datasets_yaml)
        # The old file began with `datasets:`; that key now heads the section.
        body = body.split("datasets:", 1)[1].lstrip("\n") if "datasets:" in body else body
        datasets_section = _section_text("datasets", "\ndatasets:\n") + body
        if not datasets_section.endswith("\n"):
            datasets_section += "\n"
    ndp_section = _section_text("ndp", "\nndp:\n")
    if ndp_yaml:
        body = _strip_leading_comments(ndp_yaml)
        ndp_section = _section_text("ndp", "\nndp:\n") + _indent(body)
    return settings.HEADER + folder_section + lengths + datasets_section + ndp_section


def _move(src, dst, moved):
    """Move a file or directory, keeping what is already at the destination."""
    if not os.path.exists(src) or os.path.exists(dst):
        return False
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    os.rename(src, dst)
    moved.append((src, dst))
    return True


def migrate(folder, verbose=True):
    """Move `folder` to the current layout. Returns what was moved."""
    base = layout.base(folder)
    old_index, old_mddb = layout.legacy(base)
    result = {"folder": base, "moved": [], "settings": None, "kept": []}
    if not needed(base):
        return result

    def read(path):
        try:
            with open(path) as fh:
                return fh.read()
        except OSError:
            return None

    substrate_yaml = read(os.path.join(old_index, "substrate.yaml"))
    datasets_yaml = read(os.path.join(old_mddb, "datasets.yaml"))
    ndp_yaml = read(os.path.join(old_mddb, "ndp.yaml"))

    # 1. Keep the originals, verbatim, before touching anything.
    keep = os.path.join(layout.internal(base),
                        f"migrated-{time.strftime('%Y-%m-%d')}")
    os.makedirs(keep, exist_ok=True)
    for name, text in (("substrate.yaml", substrate_yaml),
                       ("datasets.yaml", datasets_yaml),
                       ("ndp.yaml", ndp_yaml)):
        if text is not None:
            with open(os.path.join(keep, name), "w") as fh:
                fh.write(text)
            result["kept"].append(os.path.join(keep, name))

    # 2. One settings.yaml with everything those three said.
    os.makedirs(layout.results(base), exist_ok=True)
    result["settings"] = settings.write(
        base, build_settings_text(substrate_yaml, datasets_yaml, ndp_yaml))

    # 3. Move the machinery. Renames, so a 6 TB tables/ and a work/ link tree
    #    cost nothing and no link is rewritten.
    moved = result["moved"]
    for src, dst in (
            (os.path.join(old_index, "registry"), layout.registry(base)),
            (os.path.join(old_index, "tables"), layout.tables(base)),
            (os.path.join(old_index, "curated.yaml"), layout.notes(base)),
            (os.path.join(old_index, "LEGEND.md"), layout.legend(base)),
            (os.path.join(old_index, "logs"), layout.logs(base)),
            (os.path.join(old_mddb, "registry"), layout.mddb_registry(base)),
            (os.path.join(old_mddb, "exports"), layout.exports(base)),
            (os.path.join(old_mddb, "work"),
             os.path.join(layout.mddb_dir(base), "work")),
            (os.path.join(old_mddb, "cache"), layout.cache(base)),
            (os.path.join(old_mddb, "backups"), layout.backups(base)),
            (os.path.join(old_mddb, "reports", "pipeline_status.json"),
             layout.status(base)),
            (os.path.join(old_mddb, "ndp", "index.json"), layout.ndp_index(base)),
            (os.path.join(old_mddb, "logs"),
             os.path.join(layout.logs(base), "mddb")),
    ):
        _move(src, dst, moved)

    # 4. Say where everything went, in the directories a person will still open.
    for old in (old_index, old_mddb):
        if os.path.isdir(old):
            with open(os.path.join(old, MOVED_NAME), "w") as fh:
                fh.write(_moved_note(base, old))

    if verbose:
        print(f"\n  Moved {os.path.basename(base)} to the current layout:")
        print(f"    settings          -> {result['settings']}")
        for src, dst in moved:
            print(f"    {os.path.relpath(src, base):<34} -> "
                  f"{os.path.relpath(dst, base)}")
        print(f"    the three old config files are kept, unchanged, in "
              f"{os.path.relpath(keep, base)}/")
        print(f"    {layout.LEGACY_INDEX}/ and {layout.LEGACY_MDDB}/ are left "
              f"in place, holding {MOVED_NAME}; delete them when you are happy.\n")
    return result


def _moved_note(base, old):
    return f"""# This directory has moved

The pipeline now writes one visible folder and one hidden one:

    {layout.RESULTS}/       README.md, settings.yaml, glossary.md,
                     NDP_registration.md, and the NDP record(s)
    {layout.INTERNAL}/   the index, the caches, the workflow's link tree, logs

What was in `{os.path.basename(old)}/` is in `{layout.INTERNAL}/`, and the three
YAML files that were edited by hand are now the sections of
`{layout.RESULTS}/settings.yaml`. Copies of them, exactly as they were, are in
`{layout.INTERNAL}/migrated-*/`.

Nothing was deleted. This directory can be removed once you are happy with
`{layout.RESULTS}/`.
"""
