#!/usr/bin/env python3
"""
io.py -- one row per measurement, and a manifest that says where it came from.

TIDY, LONG FORMAT, because that is what this lab already writes
(`ANALYSIS/agg/ALL.tsv`, 596,602 rows) and what every plotting and aggregation
step can consume without knowing the analysis. One row = one number, carrying
enough context to be read on its own.

THREE COLUMNS THAT ARE NEVER OPTIONAL:

    ns_per_frame   because the projects differ -- OPEN is 0.2 ns/frame,
                   OCCLUDED-OPEN and CLOSED are 0.1 -- and a table that pools
                   them without saying so has already produced a wrong figure
                   in this workspace.
    stride         so a pilot run at stride 100 can never be mistaken for a
                   production run at stride 1.
    resolved       whether the sampling actually supports the observable. A
                   figure greys out an unresolved value; a caption refuses to
                   quote it. Better than a plausible number with no warning.

DETERMINISM. The data file carries no timestamp and no absolute path, and rows
are written in a fixed column order with fixed float formatting, so running the
same analysis twice produces byte-identical output. Only the manifest records
when it ran -- and the metadata record deliberately does not copy that field,
or every rebuild would mark every record changed.
"""

import hashlib
import json
import os
import tempfile

#: The leading block every row carries, in this order.
COLUMNS = [
    "system", "replica", "scope", "scope_id", "stage", "observable",
    "subset", "value", "unit",
    "n_frames", "ns_sampled", "ns_per_frame", "stride", "resolved",
    "floor_note", "param_name", "param_value",
    "segid", "resid", "resname", "label",
    # An AGGREGATE row: `value` is the mean, `sd` the spread it was taken
    # over, `n` how many independent things went into it. Written by any stage
    # that pools protomers and replicas, and read by the page to draw a band.
    "sd", "sem", "n",
]

#: What `scope` means, so every stage agrees.
#:
#:   frame      one value at one time point of one replica
#:   residue    one value per residue, averaged over the run
#:   run        one value for one replica
#:   aggregate  one value pooled over protomers AND replicas: `value` is the
#:              mean, `sd` the standard deviation across the pooled units, `n`
#:              their count. The three protomers of one trimer share a box and
#:              a thermostat -- they are pseudoreplicates, not independent
#:              samples -- so `sem` is reported over REPLICAS only, and a row
#:              that pools protomers says so in `param_name`.
SCOPES = ("frame", "residue", "run", "aggregate", "bin", "pair", "site")

#: What a missing value means, kept as the same sentinel the existing tidy
#: tables use, so the two can be concatenated.
NA = "NOT_APPLICABLE"


def human_ns(ns, unit=True):
    """A duration a person reads: ns below a microsecond, µs above it.

    One definition, used by the pages and by the metadata pipeline's tags, so
    a record and the page describing it can never disagree about how long a
    run was. 1512 ns -> "1.5 µs"; 504 ns -> "504 ns"; 2000 -> "2 µs".
    """
    if ns is None:
        return ""
    try:
        ns = float(ns)
    except (TypeError, ValueError):
        return ""
    if abs(ns) >= 1000:
        text = f"{ns / 1000:.1f}".rstrip("0").rstrip(".")
        return f"{text} \u00b5s" if unit else text
    text = f"{ns:.1f}".rstrip("0").rstrip(".")
    return f"{text} ns" if unit else text


def absent(v):
    """True when a column holds no value, whichever reader produced it.

    `core.io.read` maps the sentinel to None; a reader that parses the TSV
    directly sees the string "NOT_APPLICABLE". A filter written against one
    of those matched nothing in the other -- which is how every structure
    figure came out empty while the pages were fine.
    """
    return v is None or str(v).strip() in ("", "None", NA)


def _fmt(v):
    if v is None:
        return NA
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        if v != v:                                   # NaN
            return NA
        return f"{v:.6g}"
    return str(v)


class Rows:
    """A growing tidy table. Extra keys are kept, appended after the block."""

    def __init__(self, **defaults):
        self.defaults = defaults
        self.rows = []
        self.extra = []

    def add(self, **kw):
        row = dict(self.defaults)
        row.update(kw)
        for k in row:
            if k not in COLUMNS and k not in self.extra:
                self.extra.append(k)
        self.rows.append(row)
        return self

    def __len__(self):
        return len(self.rows)

    @property
    def columns(self):
        return COLUMNS + self.extra

    def to_tsv(self):
        cols = self.columns
        out = ["\t".join(cols)]
        for r in self.rows:
            out.append("\t".join(_fmt(r.get(c)) for c in cols))
        return "\n".join(out) + "\n"

    def to_frame(self):                                    # pragma: no cover
        import pandas as pd
        return pd.DataFrame(self.rows, columns=self.columns)


def _atomic(path, writer):
    """Call `writer(tmp)` on a temporary name beside `path`, then rename."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    os.close(fd)
    try:
        writer(tmp)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return path


def _atomic_write(path, text):
    def _w(tmp):
        with open(tmp, "w") as fh:
            fh.write(text)
    return _atomic(path, _w)


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def write(out_dir, stem, rows, manifest=None, parquet=True):
    """Write one analysis's rows, and the manifest beside them.

    Returns the list of files written, relative to `out_dir`.
    """
    written = []
    tsv = os.path.join(out_dir, f"{stem}.tsv")
    _atomic_write(tsv, rows.to_tsv())
    written.append(os.path.basename(tsv))

    if parquet and len(rows):
        # A column that holds both floats and strings -- `param_value` carries
        # a s.e.m. on one row and "fit_on: receptor" on the next -- becomes
        # pandas `object`, and `to_parquet` raises on it. That used to be
        # swallowed here, so `structure` wrote its TSV and silently NO parquet
        # for weeks. Mixed columns are written as text, and a failure is said
        # out loud rather than costing a file nobody notices is missing.
        try:                                               # pragma: no cover
            import pandas as pd
            df = rows.to_frame()
            for col in df.columns:
                if df[col].dtype == object:
                    kinds = {type(v) for v in df[col] if v is not None}
                    if len(kinds) > 1:
                        df[col] = df[col].map(lambda v: None if v is None else str(v))
            pq = os.path.join(out_dir, f"{stem}.parquet")
            _atomic(pq, lambda tmp: df.to_parquet(tmp, index=False))
            written.append(os.path.basename(pq))
        except Exception as exc:                           # noqa: BLE001
            # The previous run's parquet must not sit beside this run's TSV
            # looking current: a reader that prefers parquet would get the
            # old numbers.
            stale = os.path.join(out_dir, f"{stem}.parquet")
            if os.path.exists(stale):
                os.unlink(stale)
            print(f"        note: no parquet for {stem} "
                  f"({type(exc).__name__}: {str(exc)[:80]}); the TSV is complete")

    if manifest is not None:
        man = dict(manifest)
        man["outputs"] = {f: digest(os.path.join(out_dir, f)) for f in written}
        man["rows"] = len(rows)
        _atomic_write(os.path.join(out_dir, f"{stem}.manifest.json"),
                      json.dumps(man, indent=1, sort_keys=True) + "\n")
        written.append(f"{stem}.manifest.json")
    return written


def read(out_dir, stem=None):
    """Read tidy rows back from disk, so a page can be rebuilt without a
    trajectory. Values that were written as numbers come back as numbers; the
    sentinel comes back as None, not as the string "NOT_APPLICABLE"."""
    rows = []
    for name in sorted(os.listdir(out_dir)):
        if not name.endswith(".tsv") or (stem and not name.startswith(stem)):
            continue
        with open(os.path.join(out_dir, name)) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                vals = line.rstrip("\n").split("\t")
                row = {}
                for k, v in zip(cols, vals):
                    if v == NA or v == "":
                        row[k] = None
                        continue
                    try:
                        row[k] = int(v) if v.lstrip("-").isdigit() else float(v)
                    except ValueError:
                        row[k] = v
                rows.append(row)
    return rows


def manifest(stage, spec, inputs, stride=None, extra=None):
    """What produced a result: the spec it followed and the files it read.

    Inputs are recorded as [size, mtime], keyed by their path relative to the
    directory all of them share. Keying by file name made `R1/merged.dcd` and
    `R2/merged.dcd` one entry -- the last replica overwrote the others, so a
    change to any other replica went unrecorded. Relative, so the manifest
    still names no absolute path and survives the folder moving.
    """
    inputs = list(inputs)
    dirs = [os.path.dirname(os.path.abspath(p)) for p in inputs]
    try:
        common = os.path.commonpath(dirs) if dirs else ""
    except ValueError:                                     # different drives
        common = ""
    ins = {}
    for p in inputs:
        key = (os.path.relpath(os.path.abspath(p), common) if common
               else os.path.abspath(p)).replace(os.sep, "/")
        try:
            st = os.stat(p)
            ins[key] = [st.st_size, int(st.st_mtime)]
        except OSError:
            ins[key] = None
    import time
    return {"stage": stage, "system": spec.id, "bundle": spec.bundle,
            "spec_hash": spec_hash(spec), "stride": stride,
            "selections": {k: v for k, v in sorted(spec.selections.items())},
            "inputs": ins, "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
            **(extra or {})}


def spec_hash(spec):
    """A hash over everything that decides what an analysis computes.

    Selections, components, strides and the MD list -- not labels, not paths
    outside the inputs, and never the time. Change what is measured and every
    manifest built from it is invalidated; rename the system and nothing is.
    """
    d = spec.as_dict()
    keep = {"components": d["components"], "bundle": d["bundle"],
            "selections": d["selections"], "stride": d["stride"],
            "mds": [{"id": m["id"], "framestep_ns": m["framestep_ns"],
                     "trajectories": [os.path.basename(t)
                                      for t in m["trajectories"]]}
                    for m in d["mds"]]}
    return hashlib.sha1(
        json.dumps(keep, sort_keys=True).encode()).hexdigest()[:12]


def needs_run(out_dir, stem, man):
    """True when the result is missing or its inputs/spec have changed."""
    path = os.path.join(out_dir, f"{stem}.manifest.json")
    if not os.path.exists(path):
        return True
    try:
        with open(path) as fh:
            old = json.load(fh)
    except Exception:                                      # noqa: BLE001
        return True
    for key in ("spec_hash", "stride", "inputs", "stage"):
        if old.get(key) != man.get(key):
            return True
    for name, dig in (old.get("outputs") or {}).items():
        p = os.path.join(out_dir, name)
        if not os.path.exists(p) or digest(p) != dig:
            return True
    return False
