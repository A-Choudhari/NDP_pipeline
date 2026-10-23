#!/usr/bin/env python3
"""
yamledit.py -- change a YAML file without losing what a person wrote around it.

settings.yaml is one file that both a person and the tools write: the tools
fill in what they measure, the person writes the science, and neither may
clobber the other. Round-tripping through PyYAML would drop every comment and
reflow every block, so both operations here work on the TEXT.

    fill_empty(text, updates, path=["datasets"])   replace only EMPTY values,
                                                   and only the ones named
    insert_under(text, "datasets", block)          add a block to a section,
                                                   before the next top-level key

`fill_empty` came from enrich_pdb_metadata.py, which has filled datasets.yaml
this way since 2026-09; it is generalised here to any section of any document
and shared by both writers. It refuses aliases, merge keys and duplicate keys
rather than risk editing the wrong entry, and it re-parses the result and
compares the whole document before returning it.
"""

import json
from copy import deepcopy

import yaml


def empty(value):
    """True when a field says nothing: absent, null, blank, or an empty list."""
    return value is None or value == "" or value == [] or value == {}


def _section(document, node, path):
    """(list of entries, list of their nodes) for `path` in the document."""
    for name in path:
        if not isinstance(node, yaml.MappingNode):
            raise ValueError(f"{'.'.join(path)}: not a mapping")
        found = [v for k, v in node.value if k.value == name]
        if len(found) != 1:
            raise ValueError(f"the document must contain exactly one `{name}:`")
        node = found[0]
        document = (document or {}).get(name)
    if not isinstance(node, yaml.SequenceNode):
        raise ValueError(f"`{'.'.join(path)}` must be a list")
    return document or [], node.value


def fill_empty(text, updates, path=("datasets",), key="id"):
    """Fill empty fields of the entries named in `updates`; return the new text.

    `updates` is {entry id: {field: value}}. A field that already has a value is
    never replaced -- that is a person's writing, and this raises rather than
    overwrite it.
    """
    document = yaml.safe_load(text)
    root = yaml.compose(text)
    if not isinstance(root, yaml.MappingNode):
        raise ValueError("the document must be a mapping")
    entries, nodes = _section(document, root, list(path))
    expected = deepcopy(document)
    target_entries, _ = _section(expected, root, list(path))
    edits = []
    for spec, mapping, target in zip(entries, nodes, target_entries):
        changes = updates.get(str(spec.get(key)), {})
        if not changes:
            continue
        if not isinstance(mapping, yaml.MappingNode) or mapping.flow_style:
            raise ValueError("use block-style entries so comments survive")
        pairs = {k.value: (k, v) for k, v in mapping.value}
        if len(pairs) != len(mapping.value) or "<<" in pairs:
            raise ValueError("duplicate keys and YAML merges cannot be updated")
        additions = []
        for field, value in changes.items():
            if not empty(spec.get(field)):
                raise ValueError(f"refusing to replace nonempty "
                                 f"{spec.get(key)}.{field}")
            target[field] = value
            encoded = json.dumps(value, ensure_ascii=False)
            if field in pairs:
                k, v = pairs[field]
                if v.start_mark.index < k.end_mark.index:
                    raise ValueError("YAML aliases cannot be updated safely")
                suffix = "\n" if v.end_mark.line > v.start_mark.line else ""
                prefix = " " if text[v.start_mark.index - 1:v.start_mark.index] == ":" else ""
                edits.append((v.start_mark.index, v.end_mark.index,
                              prefix + encoded + suffix))
            else:
                additions.append((field, encoded))
        if additions:
            id_key, id_value = pairs[key]
            if id_value.start_mark.line != id_value.end_mark.line:
                raise ValueError(f"`{key}` must be a single-line scalar")
            lines = text.splitlines(keepends=True)
            end = sum(len(line) for line in lines[:id_value.end_mark.line + 1])
            indent = " " * id_key.start_mark.column
            block = ("" if text[:end].endswith("\n") else "\n")
            block += "".join(f"{indent}{f}: {v}\n" for f, v in additions)
            edits.append((end, end, block))
    patched = text
    for start, end, block in sorted(edits, reverse=True):
        patched = patched[:start] + block + patched[end:]
    if yaml.safe_load(patched) != expected:
        raise ValueError("YAML preservation check failed; nothing was changed")
    return patched


def insert_under(text, section, block):
    """Add `block` at the end of a top-level `section:`, before the next one.

    Appending to the end of the file would put a new dataset under whatever
    section happens to come last. Top-level keys start in column 0, which is
    all this needs to find where the section ends.
    """
    lines = text.splitlines(keepends=True)
    start = next((i for i, line in enumerate(lines)
                  if line.startswith(f"{section}:")), None)
    if start is None:
        raise ValueError(f"no `{section}:` section in this document")
    end = len(lines)
    for i in range(start + 1, len(lines)):
        stripped = lines[i].rstrip("\n")
        if stripped and not stripped[0].isspace() and not stripped.startswith("#"):
            end = i
            break
    # Keep the blank line that separates sections below the inserted block.
    while end > start + 1 and not lines[end - 1].strip():
        end -= 1
    if not block.endswith("\n"):
        block += "\n"
    return "".join(lines[:end]) + block + "".join(lines[end:])
