# NDP_PIPELINE — one folder in, one JSON per simulated system out

Point it at a folder of finished simulations. It reads every trajectory header,
every MD config and every run log once, and leaves a record per simulated
system: enough to register the dataset on the
[National Data Platform](https://nationaldataplatform.org) and to list it in the
lab's data catalogue (`catalog_url` in `../lab.yaml`), plus the same record in
standard formats (schema.org/DCAT, RO-Crate) for any other catalogue.

```bash
export SUBSTRATE_PYTHON=$(which python3)   # a python3 with pandas, pyarrow, pyyaml
$KIT/NDP_PIPELINE/ndp_pipeline.sh /path/to/<folder> --yes
```

It reads your data and never writes to it. Everything it produces is inside the
folder you name, and it reads nothing outside that folder.

## What you get

```
<folder>/MD_RESULTS/
  README.md                  what was found, what to check, what to do next
  settings.yaml              THE file you edit
  glossary.md                every MD: atoms, sampling, force field, evidence
  glossary.parquet           the same, for scripts
  catalog.json               every record of this folder, for the lab catalogue
  NDP_registration.md        what to paste into each field of the NDP form
  spike-fulllength-open_6x703.7ns.json   one record per simulated system
  ...                        (more than six: they move into NDP_records/)
  standard_formats/          every record as schema.org/DCAT and RO-Crate:
    catalog.dcat.jsonld        the folder as one DCAT catalogue
    <record>/dataset.jsonld    schema.org Dataset + dcat:Dataset
    <record>/ro-crate-metadata.json   RO-Crate 1.1, Process Run Crate
    <record>/record.json       the record, so the crate is self-contained

<folder>/.md_pipeline/       the index, caches and logs; never needs opening
```

A record's name says where it came from and how much of it there is:
`<folder>_<system>_<replicas>x<length>.json`. Write a `title:` and the title
names it instead; write a `slug:` and that wins.

## The record

Each JSON carries what the NDP registration guide asks for **and** what the lab
catalogue shows:

```json
{ "identifier": "SMITHLAB-SARSCOV2-SPIKE-OPEN",   "topic": "covid-19",
  "tags": ["spike", "open-state", "glycans", "SARS-CoV-2"],
  "title": "SARS-CoV-2 full-length spike glycoprotein, open state: all-atom MD",
  "description": "All-atom MD simulations of the full-length SARS-CoV-2 spike ...",
  "size_gb": 36.87, "files": 6,
  "system":   {"construct": ..., "starting_structure": "PDB 6VSB (chain A RBD 'up')",
               "glycan_profile": ..., "total_atoms": 72759},
  "setup":    {"force_field": "CHARMM36 all-atom additive", "water_model": "TIP3P",
               "ions": ..., "software": "NAMD 2.14 (Frontera, TACC)", "ensemble": "NPT",
               "temperature_K": 310, "pressure_bar": 1.01325, "timestep_fs": 2},
  "sampling": {"replicas": 6, "complete_replicas": 6, "frame_interval_ps": 100,
               "replica_lengths_ns": [1000.5, 1000.3, 1006.3, 404.3, 404.3, 406.5],
               "total_ns": 4222.2,
               "replicas_detail": [{"id": "OPEN_R1", "length_ns": 1000.5,
                                    "complete": true, "frames": 10005, "path": "OPEN"}, ...]},
  "protocol": [ ... the production stage (all stages with `protocol: all`) ... ],
  "data_availability": ..., "related_publication": "https://doi.org/10.1021/acscentsci.0c01056",
  "license": "CC-BY-4.0",
  "data_owner":    {"name": "Pat Investigator", "email": ..., "orcid": null,
                    "affiliation": "Example University"},
  "access_level":  null,     # a person's answer: public | private | licensable
  "ai_training":   null,     # a person's answer: allowed | prohibited | by-agreement
  "data_location": {"storage": "Smith Lab cluster storage",
                    "path": "/data/.../TO_SHARE_FOR_PUBLICATION/OPEN", "uri": null},
  "lab_glossary": {"record_id": ..., "md_keys": [...]} }
```

The full file is in `../examples/SARS-CoV-2-Spike_Casalino2020_MD_RESULTS/`.

Everything measured is filled in for you. `identifier`, `topic` and `title` are
yours: the first two are **suggested** in `settings.yaml` from the folder's own
names, and the record counts as `incomplete` until you confirm them.

`access_level` and `ai_training` are yours too, and nothing guesses them: the
record is `incomplete` until someone answers, lab-wide in `../lab.yaml` (the
kit ships `ai_training: by-agreement` as the lab policy and `access_level`
unanswered) or under
`ndp: defaults:`, a dataset's block, or one system. `ai_training` is a
permission of its own: a CC-BY licence says how to cite the data, not whether
it may train a model. `data_owner` and `data_location.storage` come from
`lab.yaml` unless a folder says otherwise; `data_location.path` is measured.

## settings.yaml — the one file you edit

Four sections. The tools fill in what they measure and never overwrite a value
you wrote; re-run after editing.

Every field a record still misses has an empty line in that system's block
under `ndp: systems:` (the build adds one, marked `# missing: fill in for this
system`, when a field goes missing after the block was written). Answer it
there, per system.

| section | what it holds |
|---|---|
| `folder:` | what is indexed: corpora, system labels, and `stages:` — what a trajectory file really is |
| `production_lengths_ns:` | how long a full production MD is here (default `[500, 1000, 1500]`) |
| `datasets:` | the science per dataset: authors, contact, licence, paper, PDB entries, water model |
| `ndp:` | per system: identifier, topic, tags, title, construct, glycans, ions, publication; under `defaults:` the lab-wide values, including `data_owner`, `access_level`, `ai_training`, `data_location` |

## What is published: the production phase, and only that

A record describes the production simulations. Minimisation, heating and NPT
equilibration are read (they are how the system was prepared, and the glossary
keeps them), but they are **not** in the JSON. The production stage reports the
whole run — every segment of a restarted job, its frames, its size, its length
— not the one window the config's `outputName` happens to name.

Put the preparation history back, for a reviewer who asks how it was built:

```yaml
ndp:
  defaults:
    protocol: all      # production (default) | all | none
```

### Seeing exactly what goes out

```bash
python3 $K/production_view.py <folder>            # build it
python3 $K/production_view.py <folder> --remove   # take it away
```

`MD_RESULTS/production_files/<record>/<replica>/` then holds **symlinks** to
just those files: `topology.psf`, the production segments numbered in time
order, and the `.conf` (and `.log`) they ran under. Nothing is copied or moved
— moving ~450 GB would break the restart chain on the cluster, which resumes from the
`pr.part*.restart.xsc` beside each trajectory, and would cost the record the
one thing it needs most: where the data really lives.

`rsync -L`, `cp -L` and `tar -h` follow symlinks. Sync that directory without
those flags.

## What counts as production

- The stage is read from the file's name and its directory: `min*`, `heating`,
  `npt*`/`eq*`, `*stripped*`, `*merged*`, `-old`, `frames`/`stride` are not the
  simulation; `pr*`, `prod*`, `md*` are.
- A `.dcd` whose name says nothing is production **by default**, unless another
  file in its directory is already named as production, or the directory is
  named for something else (`ANALYSIS/`, `SCRIPTS/`, `07_EQUILIBRATION/`). Such
  an MD is `review` with the issue `stage_assumed` until you confirm it.
- A file with fewer atoms than the production beside it is a solute-only
  extract, not a second simulation.
- Correct any of it under `folder: stages:` — a stage, then a regex (single
  quotes), matched on the file's path relative to the folder:

```yaml
folder:
  stages:
    derived:  ['/RMSF/all\.dcd$', '_movie\.dcd$']
    stripped: '/traj/.*_prot\.dcd$'
    production: 'base_model\..*\.R[0-9]+\.dcd$'
```

## One system, several replicas

MDs that ran from a **byte-identical topology** are replicas of one system, and
one record. Every replica is listed with how long it ran and whether it reached
the length its config asked for — a run that stopped is part of the record, not
something the catalogue hides. `lengths: nominal` under `ndp: defaults:` keeps
only the full-length ones instead.

Replicas side by side in ONE directory are told apart by name: `R1`, `1R`,
`rep1`, `replica_2`, or a bare `_1 … _6` closing the name. The bare number counts
only when every file's header shows it starting at the same step, so restart
segments stay one MD. Anything else is declared under `folder:`:
`replicas: ['_run_(\d+)\.dcd$']` (the capture group is the replica number).
A system run to more than one length on purpose (3 × 1 µs + 3 × 400 ns adaptive)
declares `nominal_lengths_ns: [1000, 400]` in its `ndp: systems:` block, which
turns the two length warnings into a note.

`stages:`, `replicas:`, `production_lengths_ns:` and `frame_interval_ps` reach a
record only through the glossary: after editing them, run `ndp_pipeline.sh`, not
`ndp_build.py` alone (which says so when it notices). Identifier, topic and tags
left blank are guessed from folder names; the guesses are listed in
`NDP_registration.md` for a person to confirm. A record's file name, once
written, stays even if its sampling changes; `rename: true` in its block lets it
move.

## When a trajectory cannot say how long it is

Files written by a stripping or concatenation tool carry no timestep, so their
length cannot be measured — the record then says `5-replicas` rather than
inventing one. Write the spacing you used under that system in `settings.yaml`:

```yaml
ndp:
  systems:
    - id: MDs/70c12cdf
      frame_interval_ps: 100     # 10,500 frames -> 1050 ns per replica
```

The lengths, the total and the file name follow from it, and every value it
produced is marked `"length_basis": "declared frame interval"` so a reader can
tell it from a measurement.

## The steps, if you prefer to run them yourself

```bash
K=$KIT/NDP_PIPELINE
python3 $K/init.py           <folder> --yes    # settings.yaml
python3 $K/index_build.py    <folder>          # read every trajectory header
python3 $K/datasets_init.py  <folder> --yes    # propose datasets
python3 $K/enrich_pdb.py     <folder> --apply  # fill EMPTY descriptions, PDB ids
python3 $K/glossary_build.py <folder>          # the MD glossary
python3 $K/ndp_build.py      <folder>          # the records + catalog.json
python3 $K/production_view.py <folder>         # links to just what is published
python3 $K/glossary_query.py <folder> summary  # read the glossary back
```

`ndp_build.py --check` writes nothing and exits non-zero unless every record is
ready — useful before you register anything.

## Registering, and the lab catalogue

Nothing is uploaded. Host the records (a lab GitHub repository, or Zenodo for a
DOI), set `base_url` under `ndp: defaults:` so the sheet shows full item URLs,
and register by hand on nationaldataplatform.org — `NDP_registration.md` has the
text for every field of the form. `catalog.json` is the same set of records in
the shape the lab catalogue lists them. Keep entries *Unlisted* until the paper
they support is published.

## The other pipeline

`../MDDB_PIPELINE/` prepares the MDDB deposition files (`inputs.yaml` and the
link-only working area for `mwf`) from the glossary this one builds. Run this
pipeline first; that one never writes `settings.yaml` or the glossary.

## Tests

```bash
for t in $K/tests/*.py; do python3 "$t"; done
```

They use synthetic fixtures only: nothing reads your data, nothing writes
outside a temporary directory.
