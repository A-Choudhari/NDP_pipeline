# MD_COLLECT — finished records, gathered where the lab can see them

The two pipelines each write into the folder they were run on. That is right:
the records belong beside the data they describe. It also means that after
twenty projects there are twenty `MD_RESULTS/` folders scattered eight
directories deep, and no way to see what the lab has.

This kit copies them into one place.

```bash
export SUBSTRATE_PYTHON=$(which python3)
$KIT/MD_COLLECT/collect_pipeline.sh --all --dry-run     # what would change
$KIT/MD_COLLECT/collect_pipeline.sh --all
```

Each project's `MD_RESULTS/` becomes
`$MD_COLLECT_DEST/<PREFIX>_MD_RESULTS/` (default `~/MD_NDP_DATA`).

**It is a separate pipeline on purpose.** It builds no records and analyses
nothing. It reads finished output and moves it. `NDP_PIPELINE` and
`MD_ANALYSES` know nothing about it, and it never writes into a source.

## Naming

Put one line at the **top** of that project's `MD_RESULTS/settings.yaml`,
above `folder:`:

```yaml
collection_prefix: SARS-CoV-2-Spike
```

That project then collects as `SARS-CoV-2-Spike_MD_RESULTS`. The metadata pipeline
never rewrites this file wholesale, so the line survives every rebuild — but
it has to be above `folder:`, because the editors that fill in the sections
append at the end of the section they are given and stop at the next
column-zero key.

Where nobody declared one, a prefix is derived from the path and `INDEX.md`
says `derived` in that row, every run, until someone names it. Two projects
reaching the same prefix stops the whole run before a byte moves, naming both
— overwriting one project's collection with another's is not a thing to
resolve automatically.

## What it copies, and what it refuses

A **real copy**, byte for byte, including `data/*.tsv`. Those TSVs are about
29% of the bytes and hold the same rows as the parquet beside them; they are
kept anyway, because a rule that drops files ends the byte-for-byte property
that makes a collection safe to hand to someone.

**Symlinks stay symlinks.** `production_files/` points at the raw trajectories
with absolute paths. Following them would turn a 2.4 GB share into terabytes
and would replace pointers-at-the-run with frozen duplicates. So:

> **`rsync -L`, `cp -L` and `tar -h` follow symlinks.** Sync this share with
> `rsync -a` (no `-L`), `cp -a`, or `tar` without `-h`. Otherwise you will
> pull in every trajectory the share points at.

It will not follow a symlinked directory, write anything outside the
destination, copy a source that has not changed, or remove anything it did not
put there.

## Re-running

A source whose shape has not changed is skipped without reading a byte: the
comparison is `lstat` only — size and nanosecond mtime per file, link targets,
directory list — against the manifest written when it was last collected.
Hashes are computed while bytes are copied, which is the only place they are
free.

**What that cannot catch:** a file rewritten to the same length within the same
second. `--verify-hashes` re-reads the share and checks it against the recorded
digests, and is the only thing that finds such a rewrite. It is not run by
default because it reads everything.

## When a source changes

The new tree is built under `.collect/tmp/`, verified, and only then renamed
into place. The previous copy is displaced first and removed only once the new
one is in position, and the manifest is written last — so "finished" means
*manifest present and its fingerprint matches the source*. A crash leaves
either the old complete collection or a `.partial` nobody claims.

Files that a person added inside a collection are **carried forward** into the
new copy and reported as kept. The manifest is the only authority for what the
collector put there, and it removes nothing else. The metadata pipeline once
came one build away from deleting 39 files by sweeping a tree it had not
written; a sweep that removes what it did not write is not a sweep.

## What lands in the share

```
MD_NDP_DATA/
  README.md  INDEX.md  index.json    what is here, from where, when, how much
  <PREFIX>_MD_RESULTS/               a pristine copy; nothing of ours inside
  .collect/                          manifests, temp trees, logs, the lock
```

Nothing of the collector's goes inside a collection, so a folder taken out of
the share is byte for byte what the source held. The cost is that such a folder
alone no longer says where it came from — which is what `index.json` is for.

## Why `--jobs` is 1

The kit's own measurements on this export: about 110 MB/s, and two hosts
reading at once each got half. Reading is the whole budget, so parallelism
buys nothing here and costs forty people their seek time.

## Tests

```bash
$SUBSTRATE_PYTHON -m pytest $KIT/MD_COLLECT/tests -q
```

Every one of them exists because something produced a result that still looked
right: a copy whose symlinks were followed, a half-copied tree named like a
finished one, a sweep that ate a colleague's note, a dry run that created the
destination, an index that quietly forgot the collections it had not touched.
