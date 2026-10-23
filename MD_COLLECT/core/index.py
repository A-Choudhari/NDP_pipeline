#!/usr/bin/env python3
"""
index.py -- what the share says to someone who opens it cold.

Three files at the top of MD_NDP_DATA, and the first thing all three say is
what this is NOT: it is not the trajectories. Somebody will otherwise take the
2.4 GB, find `production_files/`, and conclude the simulations are in here.

The per-collection rows are read out of each copy's own `catalog.json`, not out
of our manifest, so the index cannot disagree with the thing it indexes. The
manifest supplies only what the catalogue cannot know: where it came from, when
it was collected, and how it was named.

An unchanged run does not touch these files. Rewriting a byte-identical INDEX.md
would bump its mtime and make every folder look changed to whoever syncs the
share, which is the sort of thing that turns a 30-second no-op into a nightly
2 GB transfer.
"""

import json
import os
import re
import time

from . import layout


def _same(path, text):
    try:
        with open(path) as fh:
            return fh.read() == text
    except OSError:
        return False


def write_page(path, text):
    """Write only if the content differs. Returns True when it changed."""
    if _same(path, text):
        return False
    # Atomically: index.json is what a reader of the share opens first, and a
    # half-written one reads as an empty or broken collection.
    from .manifest import _atomic
    _atomic(path, text)
    return True


def catalogue_of(collected_dir):
    """(records, total_ns, identifiers) from the copy's own catalog.json."""
    p = os.path.join(collected_dir, "catalog.json")
    try:
        with open(p) as fh:
            doc = json.load(fh)
    except (OSError, ValueError):
        return 0, 0.0, []
    ds = doc.get("datasets") or []
    return (len(ds),
            sum(d.get("total_ns") or 0 for d in ds),
            [d.get("identifier") for d in ds if d.get("identifier")])


def _uncomment(text):
    """A YAML line without its comment: a `#` after a space, outside quotes."""
    quote = None
    for i, ch in enumerate(text):
        if quote:
            quote = None if ch == quote else quote
        elif ch in "'\"":
            quote = ch
        elif ch == "#" and (i == 0 or text[i - 1] in " \t"):
            return text[:i].rstrip()
    return text.rstrip()


def _scalar(val):
    val = val.strip()
    if len(val) >= 2 and val[0] == val[-1] and val[0] in "'\"":
        return val[1:-1]
    return val


def lab_name(path=None):
    """The lab's short name from the kit's lab.yaml (or $MD_PIPELINE_LAB), or
    None. PyYAML when it is there; otherwise the `short_name:` (else
    `name:`) of the top-level `lab:` block is read by hand, so the collector
    keeps needing nothing but the standard library."""
    kit = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.expanduser(path or os.environ.get("MD_PIPELINE_LAB")
                              or os.path.join(kit, "lab.yaml"))
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError:
        return None
    try:
        import yaml
        lab = (yaml.safe_load(text) or {}).get("lab") or {}
        got = lab.get("short_name") or lab.get("name") if isinstance(lab, dict) else None
        return str(got).strip() or None if got else None
    except Exception:                                        # noqa: BLE001
        pass
    found, inside, lines = {}, False, text.splitlines()
    for i, line in enumerate(lines):
        body = _uncomment(line)
        if not body.strip():
            continue
        if not line.startswith((" ", "\t", "-")):
            key, _, rest = body.partition(":")
            inside = key.strip() == "lab"
            rest = rest.strip()
            if inside and rest.startswith("{") and rest.endswith("}"):
                # lab: {name: ..., short_name: ...}
                for part in re.split(r",(?=(?:[^'\"]|'[^']*'|\"[^\"]*\")*$)", rest[1:-1]):
                    k, sep, v = part.partition(":")
                    if sep and k.strip() in ("name", "short_name"):
                        found[k.strip()] = _scalar(v)
                inside = False
            continue
        key, sep, val = body.strip().partition(":")
        if not (inside and sep and key in ("name", "short_name")):
            continue
        val = val.strip()
        if val[:1] in (">", "|"):
            # a folded or literal block: the indented lines that follow
            more = []
            for nxt in lines[i + 1:]:
                if nxt.strip() and len(nxt) - len(nxt.lstrip()) <= len(line) - len(line.lstrip()):
                    break
                if nxt.strip():
                    more.append(nxt.strip())
            val = " ".join(more)
        found[key] = _scalar(val)
    return found.get("short_name") or found.get("name") or None


def human_ns(ns):
    if not ns:
        return "-"
    return f"{ns / 1000:.2f} µs" if ns >= 1000 else f"{ns:.0f} ns"


README_TEXT = """# MD_NDP_DATA

Catalogue records and analysis bundles for the lab's finished MD, gathered in
one place. Each folder is a real, self-contained copy of one project's
`MD_RESULTS/`.

**These are not the trajectories.** They are the metadata, the tables and the
pages built from the simulations. The simulations themselves stay where they
were run. `production_files/` inside a collection holds **symlinks** pointing
at them.

Those symlinks are why one rule matters here: **`rsync -L`, `cp -L` and
`tar -h` follow symlinks.** Copying this share with any of them will try to
pull in every trajectory they point at. Use `rsync -a` (no `-L`), `cp -a` or
`tar` without `-h`.

## What is here

- `INDEX.md` — one row per collection: where it came from, when, how much.
- `index.json` — the same, for scripts.
- `<PREFIX>_MD_RESULTS/` — one project's records, glossary, registration
  sheet, settings and analysis bundles, exactly as they are at the source.

Nothing of the collector's is inside a collection: a folder taken from here
and handed to someone is byte for byte what the source folder held. The
provenance for all of them lives in `index.json`.

## Refreshing it

```bash
~/MD_COLLECT/collect_pipeline.sh --all --dry-run   # what would change
~/MD_COLLECT/collect_pipeline.sh --all
```

A collection whose source has not changed is skipped without reading a byte.
A source that has changed is re-copied and the previous copy replaced only
once the new one is complete. `--verify-hashes` re-reads the share and checks
it against what was recorded when it was collected.

## Naming

A prefix comes from `collection_prefix:` at the top of that project's
`MD_RESULTS/settings.yaml`. Where nobody declared one it is derived from the
path, and `INDEX.md` says `derived` in its row. Declaring one is a single line.
"""


def render_index(root, rows):
    """INDEX.md: what someone scanning the share wants in one screen."""
    L = ["# Collections", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')}. "
         f"{len(rows)} collection(s).", "",
         "| collection | records | sampling | files | links | collected | prefix | source |",
         "|---|---:|---:|---:|---:|---|---|---|"]
    t_rec = t_ns = t_files = t_links = 0
    for r in rows:
        L.append(f"| `{r['prefix']}` | {r['records']} | {human_ns(r['total_ns'])} "
                 f"| {r['files']} | {r['links']} | {r['collected'][:10]} "
                 f"| {r['prefix_basis']} | `{r['source_folder']}` |")
        t_rec += r["records"]
        t_ns += r["total_ns"]
        t_files += r["files"]
        t_links += r["links"]
    L += ["", f"**Total:** {t_rec} records, {human_ns(t_ns)}, "
              f"{t_files} files, {t_links} symlinks.", ""]
    derived = [r["prefix"] for r in rows if r["prefix_basis"] == "derived"]
    if derived:
        L += ["> " + ", ".join(f"`{d}`" for d in derived)
              + " took a prefix derived from the path because nobody declared"
                " one. Add `collection_prefix:` at the top of that project's"
                " `MD_RESULTS/settings.yaml`, above `folder:`, to name it"
                " yourself.", ""]
    unresolved = sum(r.get("unresolved_links", 0) for r in rows)
    if unresolved:
        L += [f"> {unresolved} symlink(s) in `production_files/` do not resolve "
              f"from here. The data they point at has moved or is not readable "
              f"by this account.", ""]
    return "\n".join(L)


def write_all(root, rows, dry=False):
    """README.md, INDEX.md and index.json. Returns what changed."""
    if dry:
        return []
    changed = []
    if write_page(os.path.join(root, layout.README), README_TEXT):
        changed.append(layout.README)
    if write_page(os.path.join(root, layout.INDEX_MD), render_index(root, rows)):
        changed.append(layout.INDEX_MD)
    # `generated` is the one field that differs on a run that changed nothing,
    # and rewriting the file for it alone bumps its mtime and makes the share
    # look changed to whoever syncs it. So the comparison ignores it, and the
    # stamp is only refreshed when something else actually moved.
    path = os.path.join(root, layout.INDEX_JSON)
    doc = {"collection": " ".join(x for x in (lab_name(), "MD metadata") if x),
           "root": root,
           "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "collections": rows}
    was = None
    try:
        with open(path) as fh:
            was = json.load(fh)
    except (OSError, ValueError):
        was = None
    if was is not None:
        a = {k: v for k, v in was.items() if k != "generated"}
        b = {k: v for k, v in doc.items() if k != "generated"}
        if a == b:
            return changed
    if write_page(path, json.dumps(doc, indent=1, sort_keys=True) + "\n"):
        changed.append(layout.INDEX_JSON)
    return changed
