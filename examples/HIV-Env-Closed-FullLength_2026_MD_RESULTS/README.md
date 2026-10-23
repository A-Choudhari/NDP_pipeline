# MDs — MD results

Everything this pipeline found under `/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs`, and what is left to do.
Owner: maalbani. Generated 2026-10-01T21:37:02+00:00.

## What is here

- **3 run directories** indexed, 200.2 GB of trajectory
- **5 MDs** in the glossary — 0 ok, 5 to review, 0 blocked
    - 5 of 1 µs

## Your files

| file | what it is |
|---|---|
| `settings.yaml` | **the one file to edit**: what is indexed, the science behind each dataset, and the catalogue record |
| `glossary.md` | every MD, with how each value was established |
| `glossary.parquet` | the same, for scripts |
| `catalog.json` | every record of this folder, for the lab catalogue |
| `NDP_registration.md` | what to paste into the NDP form, and what is missing |
| `*.json` | the record to host and link from the NDP catalogue |
| `standard_formats/` | each record again as schema.org/DCAT JSON-LD and as an RO-Crate, and the folder as one DCAT catalogue |
| `README.md` | this page |

The machinery — the index, the caches and the logs — is in `.md_pipeline/`. Nothing in it needs opening, and it can be deleted and rebuilt.

## The records

| record | identifier | replicas | sampling | volume |
|---|---|---|---|---|
| `1-FULL-GP120_5-replicas.json` | AMARO-HIV-MDS-001 | 5 (5 complete) | 5.25 µs | 161.62 GB |

`catalog.json` holds the same rows, for the lab catalogue.

## What to check

- **1 top-level directory is neither indexed nor pruned — nothing here describes them: `PICTURES`. Add them to `corpora:` or `prune:` under `folder:` in settings.yaml.
- 3 of 3 runs have no declared `system` — if the paths encode one, add a `systems:` block to `settings.yaml` under `folder:` mapping a label to a regex, then rebuild
- 5 MD(s) flagged `split_directory` — one of several replicas stored side by side in one directory, told apart by the replica label in its file name (see `glossary.md`)
- 5 MD(s) flagged `no_config` — no MD config found; timestep, temperature and ensemble unknown (see `glossary.md`)

## Next

1. Fill in what only you know in `settings.yaml` — the paper, the construct, the glycan profile, the ions, and under `ndp: defaults:` who may get the data (`access_level`) and whether it may train AI models (`ai_training`).
2. Re-run the pipeline; every measured value is refreshed and your text is left alone.
3. Host the `.json` record and register it: `NDP_registration.md` has the text for every field of the form.

## For scripts and agents

```python
import pandas as pd
runs = pd.read_parquet("/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs/.md_pipeline/index/registry/runs.parquet")
mds  = pd.read_parquet("/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs/.md_pipeline/mddb/registry/mds.parquet")
```

- `/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs/MD_RESULTS/glossary.parquet` — the glossary as a table: one row per MD, every measured value
- `/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs/.md_pipeline/coverage.md` — what could not be established, and why

