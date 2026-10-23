# MD glossary

Generated 2026-09-30T15:08:39. One entry per MD, 15 in total. Machine-readable copy: `glossary.parquet` beside this report. Query it with `glossary_query.py`.

## Terms

| term | meaning |
|---|---|
| **MD** | one continuous simulation of one system: a topology plus its production trajectory segments, in time order |
| **md_key** | `<dataset_id>/<md_id>`, the stable name of an MD. Derived from the run's path against its dataset's declared root, so it does not change on rebuild or when runs are added |
| **dataset** | a group of MDs declared in `datasets.yaml`; one MDDB project |
| **segment** | one trajectory file of an MD; restarts produce several |
| **ordering** | how segment order was established: `header` (first-frame step read from the file), `single` (one segment), `natural` (filename order, a guess) |
| **replica_index** | the repeat number read from the path (`1R`, `rep2`, `3-R3`), or empty |
| **stage_dir** | a directory name that denotes where a run lives (`pr`, `rep`, `md`), not which repeat it is |
| **topology_group** | first 8 characters of the topology file's content hash. MDs sharing a group use the same topology file content; different groups are different molecular systems or builds |
| **simulated_ns** | production time, summed per segment from each file's own header; where the headers carry no timestep, the frame count at the spacing declared in settings.yaml -- see **framestep_basis** |
| **nominal_ns** | the production length this MD reaches, among those declared for its dataset (`production_lengths_ns`, default 500, 1000 and 1500 ns; 2 % tolerance): 1008 ns is 1000. Empty when the MD is shorter than all of them, which is the issue `short_of_nominal` |
| **target_ns** | how long production was meant to run, from the production configs (first step to the furthest `run` any of them asks for); empty for engines or layouts it cannot read |
| **completion** | simulated_ns / target_ns. Below 1 (by more than one frame) the MD carries the note `incomplete` |
| **framestep_ns** | time between saved frames; empty when segments disagree |
| **framestep_basis** | where **framestep_ns** came from: `headers` if the trajectories stated it, `declared` if a person wrote `frame_interval_ps` under that system in settings.yaml because a stripping or concatenation tool had thrown the original away. Anything derived from a `declared` spacing -- simulated_ns, and the record's lengths and total -- rests on that declaration, not on a measurement |
| **work_dir** | where the MDDB workflow runs for this MD; holds only read-only links to the source |

## States

| state | meaning |
|---|---|
| `ok` | every property this toolkit can check was checked and holds. Says nothing about whether the science is right |
| `review` | usable, but confirm the listed issue before trusting a result |
| `blocked` | the MDDB workflow will reject this MD as it stands |

| issue | state | meaning |
|---|---|---|
| `topology_mismatch` | blocked | topology declares a different atom count from the trajectory |
| `mixed_formats` | blocked | segments are in different trajectory formats; mwf requires one |
| `no_topology` | blocked | no topology found for this MD |
| `topology_unchecked` | review | topology format declares no atom count, or the trajectory's is unknown |
| `ordering_unverified` | review | segment order could not be read from headers; natural filename order used |
| `grouping_unreviewed` | note | dataset holds more than one molecular system; confirm it belongs together |
| `mixed_spacing` | review | production segments were written at different frame spacings |
| `no_config` | review | no MD config found; timestep, temperature and ensemble unknown |
| `no_frame_count` | review | frame count not readable cheaply for this format |
| `stage_assumed` | review | a trajectory is production only by default (a .dcd whose name says nothing, with no other production file beside it); confirm it, or say what it is, under `folder: stages:` in settings.yaml |
| `split_directory` | note | one of several replicas stored side by side in one directory, told apart by the replica label in its file name |
| `short_of_nominal` | review | shorter than every production length declared for its dataset (production_lengths_ns): a run that stopped, a test, or a stage -- not a full production MD |
| `incomplete` | note | production stopped short of the length its config asks for (see completion); still usable, but not the planned run |
| `segment_overlap` | note | a restart segment re-writes step(s) the previous segment already holds; simulated_ns counts each step once |
| `partial_segment` | review | a segment holds FEWER frames than its own header claims: a copy still in flight, or a write that was cut off. Counts here use what is on disk -- re-run once the transfer has finished |
| `segment_gap` | review | steps are MISSING between segments: the production is not one continuous trajectory. A segment was never copied, or the run restarted from an earlier checkpoint. gap_ns says how much is absent, and simulated_ns does not include it |
| `inherited_config` | note | config taken from a template directory, cross-checked against the header |
| `chamber_topology` | note | topology is a CHARMM-in-AMBER (chamber) prmtop; MDAnalysis cannot read it |
| `relabelled_format` | note | a segment's content contradicts its extension (e.g. NetCDF named .mdcrd); its work link carries the true extension -- read via --links |

## Overview

`topology groups` counts distinct topology files by content. More than one means the dataset holds different molecular systems or builds, which is a grouping decision to make in `datasets.yaml` before deposition: split it, or mark it `heterogeneous: reviewed`.

| dataset_id | MDs | ok | review | blocked | topology groups | grouping | simulated | trajectory |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| [CLOSED](#closed) | 3 | 0 | 3 | 0 | 1 | single system | 1.73 µs | 15.1 GB |
| [MUTANT](#mutant) | 6 | 0 | 6 | 0 | 1 | single system | 4.29 µs | 36.7 GB |
| [OPEN](#open) | 6 | 0 | 6 | 0 | 1 | single system | 4.22 µs | 36.9 GB |

| issue | MDs |
|---|---:|
| `split_directory` | 15 |
| `no_config` | 15 |

## CLOSED

3 MDs across 1 topology group(s). System `unknown`.

| md_id | state | atoms | topo group | segments | ns | framestep | T / ensemble | issues | source |
|---|---|---:|---|---:|---:|---:|---|---|---|
| `CLOSED_R1` | review | 72759 | d8dbb86a | 1 | 543.5 | 0.1 | — | split_directory, no_config | `CLOSED` |
| `CLOSED_R2` | review | 72759 | d8dbb86a | 1 | 573.8 | 0.1 | — | split_directory, no_config | `CLOSED` |
| `CLOSED_R3` | review | 72759 | d8dbb86a | 1 | 614.1 | 0.1 | — | split_directory, no_config | `CLOSED` |

## MUTANT

6 MDs across 1 topology group(s). System `unknown`.

| md_id | state | atoms | topo group | segments | ns | framestep | T / ensemble | issues | source |
|---|---|---:|---|---:|---:|---:|---|---|---|
| `MUTANT_R1` | review | 71294 | 783dd3bb | 1 | 1,001.3 | 0.1 | — | split_directory, no_config | `MUTANT` |
| `MUTANT_R2` | review | 71294 | 783dd3bb | 1 | 1,036.2 | 0.1 | — | split_directory, no_config | `MUTANT` |
| `MUTANT_R3` | review | 71294 | 783dd3bb | 1 | 1,018.4 | 0.1 | — | split_directory, no_config | `MUTANT` |
| `MUTANT_R4` | review | 71294 | 783dd3bb | 1 | 411.7 | 0.1 | — | split_directory, no_config | `MUTANT` |
| `MUTANT_R5` | review | 71294 | 783dd3bb | 1 | 407.0 | 0.1 | — | split_directory, no_config | `MUTANT` |
| `MUTANT_R6` | review | 71294 | 783dd3bb | 1 | 416.6 | 0.1 | — | split_directory, no_config | `MUTANT` |

## OPEN

6 MDs across 1 topology group(s). System `unknown`.

| md_id | state | atoms | topo group | segments | ns | framestep | T / ensemble | issues | source |
|---|---|---:|---|---:|---:|---:|---|---|---|
| `OPEN_R1` | review | 72759 | e0812b9a | 1 | 1,000.5 | 0.1 | — | split_directory, no_config | `OPEN` |
| `OPEN_R2` | review | 72759 | e0812b9a | 1 | 1,000.3 | 0.1 | — | split_directory, no_config | `OPEN` |
| `OPEN_R3` | review | 72759 | e0812b9a | 1 | 1,006.3 | 0.1 | — | split_directory, no_config | `OPEN` |
| `OPEN_R4` | review | 72759 | e0812b9a | 1 | 404.3 | 0.1 | — | split_directory, no_config | `OPEN` |
| `OPEN_R5` | review | 72759 | e0812b9a | 1 | 404.3 | 0.1 | — | split_directory, no_config | `OPEN` |
| `OPEN_R6` | review | 72759 | e0812b9a | 1 | 406.5 | 0.1 | — | split_directory, no_config | `OPEN` |
