# TO_SHARE_FOR_PUBLICATION — MD results

Everything this pipeline found under `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION`, and what is left to do.
Owner: lcasalino. Generated 2026-09-30T22:08:40+00:00.

## What is here

- **3 run directories** indexed, 88.7 GB of trajectory
- **15 MDs** in the glossary — 0 ok, 15 to review, 0 blocked
    - 6 of 400 ns
    - 3 of 500 ns
    - 6 of 1 µs

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
| `spike-fulllength-closed_3x500ns.json` | AMARO-SARSCOV2-SPIKE-CLOSED | 3 (3 complete) | 1.7314 µs | 15.12 GB |
| `spike-fulllength-open-N165A-N234A_6x715.2ns.json` | AMARO-SARSCOV2-SPIKE-N165A_N234A_MUTANT | 6 (6 complete) | 4.2912 µs | 36.72 GB |
| `spike-fulllength-open_6x703.7ns.json` | AMARO-SARSCOV2-SPIKE-OPEN | 6 (6 complete) | 4.2222 µs | 36.87 GB |

`catalog.json` holds the same rows, for the lab catalogue.

## What to check

- 3 of 3 runs have no declared `system` — if the paths encode one, add a `systems:` block to `settings.yaml` under `folder:` mapping a label to a regex, then rebuild
- 15 MD(s) flagged `split_directory` — one of several replicas stored side by side in one directory, told apart by the replica label in its file name (see `glossary.md`)
- 15 MD(s) flagged `no_config` — no MD config found; timestep, temperature and ensemble unknown (see `glossary.md`)

## Next

1. Fill in what only you know in `settings.yaml` — the paper, the construct, the glycan profile, the ions, and under `ndp: defaults:` who may get the data (`access_level`) and whether it may train AI models (`ai_training`).
2. Re-run the pipeline; every measured value is refreshed and your text is left alone.
3. Host the `.json` record and register it: `NDP_registration.md` has the text for every field of the form.

## For scripts and agents

```python
import pandas as pd
runs = pd.read_parquet("/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/.md_pipeline/index/registry/runs.parquet")
mds  = pd.read_parquet("/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/.md_pipeline/mddb/registry/mds.parquet")
```

- `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/MD_RESULTS/glossary.parquet` — the glossary as a table: one row per MD, every measured value
- `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/.md_pipeline/coverage.md` — what could not be established, and why

