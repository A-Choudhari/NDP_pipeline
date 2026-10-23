#!/usr/bin/env python3
"""Fill missing MDDB descriptions and structure citations from local PDB evidence.

    python3 enrich_pdb_metadata.py /path/to/data --apply
    python3 enrich_pdb_metadata.py /path/to/data --output-dir ./metadata-preview

Without --apply, writes a proposed settings.yaml beside the reports only. Existing
values and YAML comments are preserved. Public RCSB responses are cached; local
paths and simulation data are never sent to the API. No LLM/API key is needed.
"""

import argparse
from collections import defaultdict
from copy import deepcopy
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import urllib.error
import urllib.request

import pandas as pd
import yaml

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from glossary_build import select_runs  # noqa: E402 -- the same selectors
from datasets_init import resolve  # noqa: E402

PDB_ID = re.compile(r"[1-9][A-Za-z0-9]{3}\Z")
API = "https://data.rcsb.org/rest/v1/core/entry/"
FIELDS = ("description", "pdb_ids", "citation")
MAX_HEADER = 65536


def empty(value):
    return value is None or value == "" or value == [] or value == {}


def inside(path, root):
    return path == root or root in path.parents


def structure_ids(path):
    """Read bounded PDB/mmCIF headers; filenames alone are only candidates."""
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if not name.endswith((".pdb", ".ent", ".cif", ".mmcif")):
        return []
    stem = name.rsplit(".", 1)[0]
    found = []
    m = re.fullmatch(r"(?:pdb)?([1-9][a-z0-9]{3})", stem)
    if m:
        found.append((m[1].upper(), "filename"))
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        head = fh.read(MAX_HEADER)
    for line in head.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            break
        if line.startswith("HEADER"):
            code = line[62:66].strip().upper()
            if PDB_ID.fullmatch(code):
                found.append((code, "PDB HEADER"))
    for match in re.finditer(r"(?m)^_entry\.id\s+['\"]?([1-9][A-Za-z0-9]{3})(?=['\"\s]|$)", head):
        found.append((match[1].upper(), "mmCIF _entry.id"))
    return sorted(set(found))


def scan_sources(specs, selections, root, max_up, max_files, source_map):
    """Never inherit another dataset's ancestor or leave the requested root.

    All ancestors participate in ownership, even those beyond the scan depth.
    Otherwise a shallow run could inherit a shared project folder incorrectly.
    Nearby files are suggestions. Only headers beside selected trajectories or
    explicitly curated metadata_sources establish automatic PDB associations.
    """
    owners = defaultdict(set)
    run_dirs = {}
    for spec in specs:
        did = spec["id"]
        dirs = {root / str(p) for p in selections[did].rel_path}
        run_dirs[did] = dirs
        for directory in dirs:
            while inside(directory, root):
                owners[directory].add(did)
                if directory == root:
                    break
                directory = directory.parent
    header_cache = {}
    result = {}
    for spec in specs:
        did = spec["id"]
        sources = {}
        warnings = []
        for directory in sorted(run_dirs[did]):
            for level in range(max_up + 1):
                if not inside(directory, root) or directory == root:
                    break
                if owners[directory] == {did}:
                    sources[directory] = sources.get(directory, False) or level == 0
                directory = directory.parent
        explicit = spec.get("metadata_sources") or []
        extra = source_map.get(did) or []
        if not isinstance(explicit, list) or not isinstance(extra, list):
            raise ValueError(f"{did}: metadata_sources and --sources entries must be lists")
        files = {}
        for directory, direct in sorted(sources.items()):
            try:
                for path in directory.iterdir():
                    files[path] = (direct, False)
            except OSError as exc:
                warnings.append(f"Cannot list {directory.relative_to(root)}: {exc.strerror}")
        for item in explicit + extra:
            path = (root / str(item)).resolve()
            if not inside(path, root):
                raise ValueError(f"{did}: metadata source outside data root: {item}")
            if not path.is_file():
                raise ValueError(f"{did}: metadata source is not a file: {item}")
            files[path] = (True, True)
        evidence = defaultdict(list)
        scanned = 0
        # Explicit sources take priority if the automatic discovery limit is hit.
        for path, (direct, curated) in sorted(files.items(), key=lambda kv: (not kv[1][1], str(kv[0]))):
            if not re.search(r"\.(?:pdb|ent|cif|mmcif)(?:\.gz)?$", path.name, re.I):
                continue
            if not inside(path.resolve(), root):
                warnings.append(f"Skipped structure symlink outside data root: {path.relative_to(root)}")
                continue
            if scanned >= max_files and not curated:
                warnings.append(f"Stopped after {max_files} structure files; discovery is incomplete")
                break
            scanned += 1
            try:
                if path not in header_cache:
                    header_cache[path] = structure_ids(path)
                hits = header_cache[path]
            except (OSError, EOFError) as exc:
                warnings.append(f"Cannot read {path.relative_to(root)}: {exc}")
                continue
            headers = {code for code, kind in hits if kind != "filename"}
            for code, kind in hits:
                confirmed = (curated or (direct and kind != "filename"))
                if headers and (code not in headers or len(headers) > 1):
                    confirmed = False
                evidence[code].append({"path": str(path.relative_to(root)), "kind": kind,
                                       "confirmed": confirmed, "curated_source": curated})
        result[did] = {"evidence": dict(evidence), "warnings": warnings, "files_scanned": scanned}
    return result


def atomic_write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = path.stat().st_mode & 0o777 if path.exists() else 0o644
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class RCSB:
    def __init__(self, cache, offline=False, timeout=15):
        self.cache, self.offline, self.timeout = cache, offline, timeout
        self.memo = {}

    def get(self, code):
        if code in self.memo:
            return self.memo[code]
        path = self.cache / f"{code}.json"
        error = None
        data = None
        try:
            if path.exists():
                data = json.loads(path.read_text())
                if not isinstance(data, dict) or data.get("rcsb_id", "").upper() != code:
                    data = None
            if data is None and self.offline:
                error = "No cached RCSB entry (offline)"
            elif data is None:
                request = urllib.request.Request(API + code, headers={"User-Agent": "SUBSTRATE-PDB-metadata/1.0"})
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    data = json.load(response)
                if not isinstance(data, dict) or data.get("rcsb_id", "").upper() != code:
                    raise ValueError("RCSB response does not match the requested PDB ID")
                atomic_write(path, json.dumps(data, indent=2) + "\n")
        except (OSError, ValueError, urllib.error.URLError) as exc:
            data, error = None, str(exc)
        self.memo[code] = (data, error)
        return data, error


def citation_text(entries):
    """Cite the structure publication, never claim it reports these MD runs."""
    groups = {}
    for code, entry in sorted(entries.items()):
        citation = entry.get("rcsb_primary_citation") or {}
        title = citation.get("title")
        doi = citation.get("pdbx_database_id_PubMed")
        if not title or title.lower().strip().rstrip(".") == "to be published":
            continue
        doi = citation.get("pdbx_database_id_DOI") or (f"PMID:{doi}" if doi else "")
        key = doi.lower() or title
        groups.setdefault(key, {"codes": [], "citation": citation})["codes"].append(code)
    parts = []
    for group in groups.values():
        c = group["citation"]
        authors = c.get("rcsb_authors") or []
        author = (authors[0] + (" et al." if len(authors) > 1 else "")) if authors else ""
        bits = [author, str(c.get("year") or ""), c["title"], c.get("journal_abbrev") or ""]
        reference = ". ".join(str(b).strip().rstrip(".") for b in bits if b) + "."
        if c.get("pdbx_database_id_DOI"):
            reference += " https://doi.org/" + c["pdbx_database_id_DOI"]
        elif c.get("pdbx_database_id_PubMed"):
            reference += f" https://pubmed.ncbi.nlm.nih.gov/{c['pdbx_database_id_PubMed']}/"
        parts.append(f"PDB {', '.join(group['codes'])}: {reference}")
    if not parts:
        return None
    return "Structure reference(s); these do not establish a publication for the simulations: " + " (br) ".join(parts)


def propose(spec, selected, discovery, client):
    evidence = discovery["evidence"]
    declared = spec.get("pdb_ids") or []
    if isinstance(declared, str):
        declared = [declared]
    if not isinstance(declared, list) or any(not PDB_ID.fullmatch(str(c)) for c in declared):
        raise ValueError(f"{spec['id']}: pdb_ids must contain four-character PDB accessions")
    declared = {str(code).upper() for code in declared}
    confirmed = {code for code, rows in evidence.items() if any(r["confirmed"] for r in rows)}
    # A nonempty curated PDB list is authoritative and is never extended.
    eligible = declared or confirmed
    entries, lookups = {}, {}
    for code in sorted(set(evidence) | declared):
        entry, error = client.get(code)
        lookups[code] = {"url": API + code, "error": error,
                         "title": (entry.get("struct") or {}).get("title") if entry else None,
                         "eligible": code in eligible}
        if entry and code in eligible:
            entries[code] = entry
    updates = {}
    if empty(spec.get("pdb_ids")) and entries:
        updates["pdb_ids"] = sorted(entries)
    if empty(spec.get("description")) and len(selected):
        name = spec.get("name") or spec["id"]
        # No run count, and no engine label.
        #
        # This text is written once into datasets.yaml, is human-owned from then
        # on, and is deposited verbatim. So it must state only what will not go
        # stale underneath it.
        #
        # A count does. It was measured here and frozen; when the stage
        # classifier was repaired, LC3B's export went from 157 MDs to 2 while
        # this sentence went on claiming 157, and CD4's from 4 to 3.
        # mddb_build.py now WARNS when a description quotes a count its export
        # contradicts, but the fix at the source is not to put a volatile number
        # into prose: how many MDs a dataset holds is `mds`, regenerated every
        # build.
        #
        # The engine label goes stale differently: `selected.engine` is a
        # trajectory FORMAT FAMILY, and its `NAMD/CHARMM` names no real program.
        # The record now carries `program` and `version` read from the engine's
        # own start-up banner, so a description repeating the family label would
        # contradict the fields beside it.
        description = f"Molecular dynamics dataset '{name}'"
        if spec.get("condition"):
            description += f", condition '{spec['condition']}'"
        description += "."
        if entries:
            refs = [f"PDB {code} ({(entry.get('struct') or {}).get('title', 'structure entry')})" for code, entry in sorted(entries.items())]
            description += " Referenced structures: " + "; ".join(refs) + "."
        updates["description"] = description
    if empty(spec.get("citation")) and entries:
        citation = citation_text(entries)
        if citation:
            updates["citation"] = citation
    discovery.update({"lookups": lookups, "updates": updates,
                      "unresolved": [field for field in FIELDS if empty(updates.get(field, spec.get(field)))]})
    if not len(selected):
        discovery["warnings"].append("No production runs match this dataset selector")
    return updates


def patch_yaml(original, updates):
    """Replace just empty YAML value spans; retain comments and other edits.

    Reject aliases/merges/duplicates rather than accidentally changing another
    dataset. Parse the result and compare the complete document before writing.
    """
    document = yaml.safe_load(original)
    node = yaml.compose(original)
    if not isinstance(node, yaml.MappingNode):
        raise ValueError("datasets.yaml must be a mapping")
    dataset_nodes = [v for k, v in node.value if k.value == "datasets"]
    if len(dataset_nodes) != 1 or not isinstance(dataset_nodes[0], yaml.SequenceNode):
        raise ValueError("datasets.yaml must contain one datasets list")
    expected = deepcopy(document)
    edits = []
    for spec, mapping, target in zip(document["datasets"], dataset_nodes[0].value, expected["datasets"]):
        changes = updates.get(spec["id"], {})
        if not changes:
            continue
        if not isinstance(mapping, yaml.MappingNode) or mapping.flow_style:
            raise ValueError("Use block-style dataset mappings to preserve YAML comments")
        pairs = {k.value: (k, v) for k, v in mapping.value}
        if len(pairs) != len(mapping.value) or "<<" in pairs:
            raise ValueError("Duplicate keys and YAML merges are not supported for updates")
        additions = []
        for key, value in changes.items():
            if not empty(spec.get(key)):
                raise ValueError(f"Refusing to replace nonempty {spec['id']}.{key}")
            target[key] = value
            encoded = json.dumps(value, ensure_ascii=False)
            if key in pairs:
                k, v = pairs[key]
                if v.start_mark.index < k.end_mark.index:
                    raise ValueError("YAML aliases cannot be updated safely")
                # For a block scalar/sequence, retain a newline before the next key.
                suffix = "\n" if v.end_mark.line > v.start_mark.line else ""
                prefix = " " if original[v.start_mark.index - 1:v.start_mark.index] == ":" else ""
                edits.append((v.start_mark.index, v.end_mark.index, prefix + encoded + suffix))
            else:
                additions.append((key, encoded))
        if additions:
            id_key, id_value = pairs["id"]
            if id_value.start_mark.line != id_value.end_mark.line:
                raise ValueError("Dataset ids must be single-line scalars")
            lines = original.splitlines(keepends=True)
            # Insert after the id line, unambiguously inside this dataset.
            end = sum(len(line) for line in lines[:id_value.end_mark.line + 1])
            indent = " " * id_key.start_mark.column
            text = ("" if original[:end].endswith("\n") else "\n")
            text += "".join(f"{indent}{key}: {value}\n" for key, value in additions)
            edits.append((end, end, text))
    patched = original
    for start, end, text in sorted(edits, reverse=True):
        patched = patched[:start] + text + patched[end:]
    if yaml.safe_load(patched) != expected:
        raise ValueError("YAML preservation check failed; no changes written")
    return patched


def report_markdown(report):
    lines = ["# PDB metadata enrichment", "", f"Generated {report['generated_at']}", "",
             "Only empty fields are filled. Structure citations do not identify the simulation paper.",
             "Nearby files are candidates, not established starting structures; verify them before adding pdb_ids.",
             "Descriptions summarize the indexed runs, not scientific conclusions. Counts reflect this snapshot.", ""]
    for did, row in report["datasets"].items():
        lines += [f"## {did}", "", "Fields filled: " + (", ".join(row["updates"]) or "none"),
                  "Still empty: " + (", ".join(row["unresolved"]) or "none"), ""]
        for code, lookup in row["lookups"].items():
            status = "accepted association" if lookup["eligible"] else "candidate; needs review"
            lines.append(f"- PDB {code}: {status}. {lookup['title'] or lookup['error'] or 'No title'}")
            for source in row["evidence"].get(code, [])[:5]:
                lines.append(f"  - `{source['path']}` ({source['kind']})")
        for warning in row["warnings"]:
            lines.append(f"- {warning}")
        lines.append("")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--apply", action="store_true", help="fill empty fields in settings.yaml, with a timestamped backup")
    ap.add_argument("--output-dir", type=Path, help="write preview, reports and cache here instead of .md_pipeline/mddb")
    ap.add_argument("--sources", type=Path, help="YAML mapping of dataset id to curated structure-file paths, relative to data root")
    ap.add_argument("--offline", action="store_true", help="use cached RCSB responses only")
    ap.add_argument("--max-up", type=int, default=4, help="ancestor levels for candidate discovery (default 4)")
    ap.add_argument("--max-files", type=int, default=2000, help="structure files per dataset (default 2000)")
    ap.add_argument("--timeout", type=float, default=15, help="seconds per RCSB request")
    a = ap.parse_args(argv)
    if a.apply and a.output_dir:
        ap.error("--apply and --output-dir are mutually exclusive")
    if a.max_up < 0 or a.max_files < 1 or a.timeout <= 0:
        ap.error("max-up must be nonnegative; max-files and timeout must be positive")
    base, subst, mddb = map(Path, resolve(a.folder))
    from core import layout as _layout, settings as _settings
    declaration = Path(_settings.path(str(base)))
    original = declaration.read_text()
    data = yaml.safe_load(original)
    specs = data.get("datasets") if isinstance(data, dict) else None
    if not isinstance(specs, list) or not specs or any(not isinstance(s, dict) or not s.get("id") for s in specs):
        raise ValueError("settings.yaml must contain a nonempty `datasets:` list with ids")
    if len({s['id'] for s in specs}) != len(specs):
        raise ValueError("Duplicate dataset ids")
    root = Path((data.get("folder") or {}).get("root") or base).resolve()
    if root != base.resolve():
        raise ValueError("settings.yaml `folder: root:` differs from the requested folder; rebuild the index first")
    runs = pd.read_parquet(subst / "registry" / "runs.parquet")
    selections = {s["id"]: select_runs(s, runs) for s in specs}
    source_map = yaml.safe_load(a.sources.read_text()) if a.sources else {}
    if not isinstance(source_map, dict) or set(source_map) - set(selections):
        raise ValueError("--sources must map declared dataset ids to lists of structure paths")
    out = a.output_dir.resolve() if a.output_dir else Path(_layout.internal(str(base)))
    client = RCSB(out / "cache" / "rcsb", a.offline, a.timeout)
    reports = out
    print(f"Scanning local structure metadata for {len(specs)} datasets ...", flush=True)
    discoveries = scan_sources(specs, selections, root, a.max_up, a.max_files, source_map)
    changes = {}
    for spec in specs:
        did = spec["id"]
        changes[did] = propose(spec, selections[did], discoveries[did], client)
        print(f"  {did}: {', '.join(changes[did]) or 'no empty fields resolved'}", flush=True)
    patched = patch_yaml(original, changes)
    report = {"generated_at": datetime.now(timezone.utc).isoformat(), "source": str(declaration),
              "datasets": discoveries}
    atomic_write(reports / "pdb_metadata.json", json.dumps(report, indent=2) + "\n")
    atomic_write(reports / "pdb_metadata.md", report_markdown(report))
    proposal = reports / "datasets.enriched.yaml"
    atomic_write(proposal, patched)
    if a.apply and patched != original:
        if declaration.read_text() != original:
            raise ValueError("settings.yaml changed during discovery; rerun to preserve those edits")
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        backup = Path(_layout.backups(str(base))) / f"settings.before-pdb-{stamp}.yaml"
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(declaration, backup)
        atomic_write(declaration, patched)
        print(f"Updated {declaration}; backup: {backup}")
    else:
        print(f"{'No changes needed.' if a.apply else 'Preview only.'} Proposal: {proposal}")
    errors = {code: error for code, (_, error) in client.memo.items() if error}
    if errors:
        print(f"RCSB metadata unavailable for {len(errors)} ID(s); unresolved values remain empty. See report.")
    print(f"Evidence report: {reports / 'pdb_metadata.md'}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, yaml.YAMLError) as exc:
        print(f"Metadata enrichment failed: {exc}", file=sys.stderr)
        sys.exit(2)
