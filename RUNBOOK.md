# RUNBOOK.md: the runbook for the MD → NDP publishing pipeline

**What this file is.** This is the one file that runs the whole publishing pipeline on a
folder of **finished** MD simulations. Here "the pipeline" means describe → curate → analyse →
share → register with the National Data Platform (NDP).

You can use it in two ways:

- **With Claude Code** (recommended). Open Claude Code anywhere and type:

  ```
  claude --add-dir /path/to/MD_NDP_PIPELINE
  > Read /path/to/MD_NDP_PIPELINE/RUNBOOK.md and run it on /path/to/<project-folder>
  ```

  Claude then runs the checks and steps below in order, and stops to ask you whenever it
  needs something only you know.
- **By hand.** Copy the commands below one block at a time.

The plain-language guide, with explanations for every level, is `NDP_instructions.md`
(or `.pdf`) in this same folder. Read it first if the words below are new to you.

Authors: Marco Albani and Rommie E. Amaro, UC San Diego. Questions: Marco Albani, maalbani@ucsd.edu

---

**Worked examples.** `$K/examples/` holds the finished `MD_RESULTS/` of two published systems:
the SARS-CoV-2 spike (Casalino et al. 2020) and the full-length HIV-1 Env (Nat. Commun. 2026).
Open them to see what each step produces.

---

## 0. The kit map (paths are relative to where you unzipped the archive, `$K`)

| Kit | Path | Job | Python it needs |
|---|---|---|---|
| NDP_PIPELINE | `$K/NDP_PIPELINE` | describes the data: glossary + one NDP JSON record per system | `$SUBSTRATE_PYTHON` |
| MD_ANALYSES | `$K/MD_ANALYSES` | analyses + web pages per record | `$P` (conda env `mdviz`) |
| MD_COLLECT | `$K/MD_COLLECT` | copies finished `MD_RESULTS/` into the shared lab collection | `$SUBSTRATE_PYTHON` (read by the shell script itself) |
| your collection | `$MD_COLLECT_DEST` (default `~/MD_NDP_DATA`) | the shared collection MD_COLLECT writes (data only) | none |
| MDDB_PIPELINE | `$K/MDDB_PIPELINE` | *optional*: MDDB deposition inputs for `mwf` | `$SUBSTRATE_PYTHON`, then `mwf_env` |

**Order:** NDP_PIPELINE → (edit settings.yaml) → NDP_PIPELINE `ndp_build.py` → MD_ANALYSES →
MD_COLLECT → NDP registration. MDDB comes after NDP_PIPELINE, and only if you want to deposit to MDDB.

**What goes where.** Everything is written inside your project folder:

- `<folder>/MD_RESULTS/` holds what you read and share.
- `<folder>/.md_pipeline/` holds the machinery. You can delete it; it is rebuilt.

The pipeline **never writes to your trajectories or inputs.**

---

## 1. Setup (once per shell)

```bash
export K=/path/to/MD_NDP_PIPELINE                         # where you unzipped the archive
export SUBSTRATE_PYTHON=$(which python3)                  # needs pandas + pyarrow + pyyaml (requirements-ndp.txt)
export P=$(conda run -n mdviz which python)               # analysis env: conda env create -f $K/MD_ANALYSES/env/environment.yml
export F=/path/to/<project-folder>                        # ONE project folder of finished MDs
```

Optional, when your defaults differ (the README in `$K` has the full table):

```bash
export CPPTRAJ=/path/to/cpptraj            # if cpptraj is not on your PATH (AmberTools)
export MD_COLLECT_DEST=/path/to/collection # where MD_COLLECT copies results (default ~/MD_NDP_DATA)
export MD_PIPELINE_LAB=/path/to/lab.yaml   # only to use a lab file other than $K/lab.yaml
```

**Your lab, once per installation:** fill in `$K/lab.yaml` (lab name, identifier prefix, topics,
people, data owner, storage) before the first run. It is the only lab-specific file;
`examples/lab.yaml` is a filled-in model and `lab.template.yaml` a clean copy.
`python3 $K/NDP_PIPELINE/core/lab.py` prints what is in force. A folder's own `settings.yaml` can
override any of it under `ndp: defaults:`.

**Two answers no tool can give**, needed before any record is `ready`: `access_level`
(`public` | `private` | `licensable`) and `ai_training` (`allowed` | `prohibited` | `by-agreement`).
Ask the data owner; write the answer under `ndp: defaults:` (the whole folder), a
dataset, or one system in `settings.yaml`, or once in `lab.yaml` if it holds lab-wide. Never guess.

---

## 2. Preflight checks (run all of them; stop on any FAIL)

```bash
# a) the folder exists and is readable
test -d "$F" && test -r "$F" && echo "OK folder" || echo "FAIL: $F not readable"

# b) ONE folder, not a parent of other built folders
find "$F" -mindepth 2 -maxdepth 6 -type d -name MD_RESULTS -not -path '*/.md_pipeline/*' | head
#    -> must print nothing. If it lists folders, point F at one of THEM instead.

# c) the Pythons work
$SUBSTRATE_PYTHON -c "import pandas, pyarrow, yaml; print('OK substrate python')"
$P -c "import MDAnalysis, mdtraj, plotly; print('OK mdviz')"

# d) the simulations are FINISHED (no trajectory changed in the last 60 minutes)
find "$F" \( -name '*.dcd' -o -name '*.xtc' -o -name '*.nc' \) -mmin -60 | head
#    -> must print nothing. Anything listed is still being written: wait and rerun later.

# e) the lab file is filled in
$SUBSTRATE_PYTHON -c "import sys; sys.path.insert(0, '$K/NDP_PIPELINE'); from core import lab; \
print('OK lab:', lab.name()) if lab.name() else print('FAIL: fill in', lab.path())"
```

---

## 3. Step 1: describe the data (NDP_PIPELINE)

**Folder already built?** If `$F/MD_RESULTS/settings.yaml` already exists, rerunning Step 1 is safe:
it keeps every value in `settings.yaml` and only picks up new or changed simulations. Rerun it if
simulations were added since the last build. Otherwise go straight to Step 2.

*Terms:* an **MD** is one replica run. A **system** (= one **record**, one JSON) groups the MDs
that share the same topology. For example, 4 systems × 3 replicas = 12 MDs.

```bash
$K/NDP_PIPELINE/ndp_pipeline.sh "$F" --yes
```

Useful flags:

| Flag | What it does |
|---|---|
| `--production-view` | links the production files into `MD_RESULTS/production_files/` |
| `--pdb-offline` | uses cached RCSB data only (no internet) |
| `--skip-pdb-metadata` | skips the PDB lookup step |
| `--deep` | counts frames in `.xtc`/`.mdcrd` (slow) |
| `--rebase-root` | the folder was moved or renamed |
| `--include-nested` | builds anyway when the folder holds other built folders; only if you mean it |

It runs six steps: set up → index trajectory headers → propose datasets → PDB metadata →
glossary → records.

**Success looks like:**

- `$F/.md_pipeline/logs/status.json` says `"status": "finished"`, `"exit": 0`
- `$F/MD_RESULTS/` contains:
  - `README.md`, `settings.yaml`
  - `glossary.md`/`.parquet`
  - `catalog.json`, `NDP_registration.md`
  - one `*.json` per system, in `NDP_records/` when there are more than six
  - `standard_formats/`: per record `dataset.jsonld` (schema.org + DCAT), `ro-crate-metadata.json`
    (RO-Crate) and `record.json`; plus `catalog.dcat.jsonld` for the folder

**Stop and report** if:

- it says `Refused: this folder contains N folder(s) that have already been built`. Pick one of those folders.
- it finds no production runs (`datasets_init` exit 1). The folder holds no finished production MD, or the `stages:` rules are needed (Step 2).

---

## 4. Step 2: curate (fill in what only a person knows)

Open `$F/MD_RESULTS/settings.yaml`. **It is the only file anyone edits.** Its sections:

| Section | What goes there |
|---|---|
| `collection_prefix:` | at the very **top**, above `folder:`. The name in the lab collection, e.g. `IgG-1F4` |
| `folder:` | `corpora:` (index these dirs) / `prune:` (ignore these); `stages:` rules saying what a trajectory really is; `replicas:` rules saying which files in one directory are replicas (see below) |
| `production_lengths_ns:` | the length of a full production run here, e.g. `[500, 1000]`. Shorter runs are left out |
| `datasets:` | authors, contact, licence, paper/citation, PDB entries, water model |
| `ndp: defaults:` | for every record of the folder: `access_level`, `ai_training`, `data_owner`, `data_location` (storage, uri), `base_url`, `protocol`, `lengths`. Empty = `lab.yaml`'s value |
| `ndp: systems:` | per system: **`identifier`** (`<PREFIX>-<TOPIC>-<DESCRIPTOR>-###`, from your lab's catalogue), **`topic`** (hiv, influenza, covid-19, …), `title`, `slug`, `construct`, `starting_structure`, `glycan_profile`, `frame_interval_ps`; optional `nominal_lengths_ns: [1000, 400]` (replicas run to more than one length on purpose), `rename: true` (let a slugged file name change) |

Then rebuild the records and check them. **Which command depends on what you edited:**

- You edited `folder: stages:`, `folder: replicas:`, `production_lengths_ns:` or `frame_interval_ps`?
  These reach the records only through the glossary, so rerun the **whole pipeline**
  (`$K/NDP_PIPELINE/ndp_pipeline.sh "$F" --yes`). `ndp_build.py` prints a NOTE when you forget.
- Anything else (titles, identifiers, authors, …): `ndp_build.py` is enough:

```bash
$SUBSTRATE_PYTHON $K/NDP_PIPELINE/ndp_build.py "$F"
$SUBSTRATE_PYTHON $K/NDP_PIPELINE/ndp_build.py "$F" --check     # writes nothing; read the "N/M record(s) ready" line
```

`--check` **exits 1 while any record is incomplete.** That means "not ready yet"; it is not a crash.
It prints only the count. **What is missing** is listed per record in `$F/MD_RESULTS/NDP_registration.md`.

**Suggested values.** When `identifier`, `topic` or `tags` is blank, the pipeline fills in a
guess taken from folder names, e.g. `SMITHLAB-COV2-TO-OPEN-001` from `TO_SHARE_FOR_PUBLICATION`. A
guess does **not** make a record incomplete. Instead, `NDP_registration.md` lists every guess
under *"Suggested by the pipeline: confirm before registering"*. Copy the value you want onto
that key in `settings.yaml`, and it stops being listed.

**Replicas side by side in one directory** are recognised when the names say so: `R1`, `1R`,
`rep1`, `replica_2`, or a bare `_1 … _6` closing the name. The bare number counts only when every
file's header shows it starting at the same step, so restart segments stay one run. For any
other naming, declare it:

```yaml
folder:
  replicas:
    - '_run_(\d+)\.dcd$'          # the part in ( ) is the replica number
```

Four rules to remember:

- **Set a value on the key that already exists** in the block. YAML keeps the *last* duplicate key,
  so a second, empty `slug:` further down silently wins.
- An unnamed `.dcd` counts as **production** until a `stages:` rule says otherwise. The stage names are:
  `production`, `equilibration`, `heating`, `minimisation`, `merged`, `stripped`, `superseded`, `derived`.
- Only the **production phase** goes into the NDP JSON. `ndp: defaults: protocol: all` restores
  the other phases for a reviewer.
- A record stays `incomplete` until construct, starting structure, glycan profile, publication and
  `access_level` are filled in. This is deliberate. **Every missing field has an empty line in that
  system's own block** (the build adds one, marked `# missing: fill in for this system`, if a field
  goes missing later): answer it there, per system. Values that differ between systems (construct,
  structure, glycans, ions, water, force field, software) belong in each system's block. Identifier, topic and tags are guessed; confirm them (above).
- Put a space before any `#` after a value (`contact: a@b.edu  # note`). Otherwise YAML reads the
  comment as part of the value. `ndp_build.py` prints a NOTE when it sees one.

---

## 5. Step 3: analyse (MD_ANALYSES)

```bash
$P $K/MD_ANALYSES/analyse.py "$F" --detect-only                        # writes nothing: what is each system?
$P $K/MD_ANALYSES/analyse.py "$F" --only composition,structure          # the standard first pass
```

Options:

- **More stages** (comma-separated in `--only`): `interface_contacts`, `binding_energy`, `ions`, `water_shell`,
  `water_bridges`, `water_residence`, `density_grids`, `rdf`.
- **One system:** `--system <record-name>`
- **Fewer frames:** `--stride N`
- **Redraw pages only:** `--pages-only [--page story|plotly|paper|both]`

**Output:** `$F/MD_RESULTS/analyses/<record>/`, which holds:

- `page_story.html`, `page_plotly.html`
- `frames.json`
- `data/*.tsv|*.parquet`
- `figures_user/`: drop your own VMD renders here

Optional: `$P $K/MD_ANALYSES/import_slides.py "$F" deck.pptx --dry-run` and
`$P $K/MD_ANALYSES/import_paper.py "$F" paper.pdf --dry-run`.

Long jobs: run inside `tmux`/`screen`, or with `nohup … &`.

---

## 6. Step 4: share with the lab (MD_COLLECT, optional)

MD_COLLECT copies each project's finished `MD_RESULTS/` into one shared collection, named by
`collection_prefix:` at the top of its `settings.yaml`. Whoever owns the collection runs it:

```bash
$K/MD_COLLECT/collect_pipeline.sh "$F/MD_RESULTS" --dry-run         # what would be copied
$K/MD_COLLECT/collect_pipeline.sh "$F/MD_RESULTS"                   # copy into $MD_COLLECT_DEST
```

To collect many projects at once, list their parent folders under `roots:` in
`$K/MD_COLLECT/sources.yaml` (or pass `--root DIR`) and use `--all`.

**Never** copy `production_files/` with `rsync -L`, `cp -L` or `tar -h`. Those follow the symlinks
and pull in every trajectory. Use `rsync -a`, `cp -a` or plain `tar`.

---

## 7. Step 5: register on the NDP (by hand, by a person)

1. Open `$F/MD_RESULTS/NDP_registration.md`. It lists, per record, the text to paste into the NDP
   form, in the form's order, and anything still missing.
2. Register on the NDP (nationaldataplatform.org):
   - **Prepare**: one JSON per system, hosted at a URL. Set `ndp: defaults: base_url:`. Host
     `standard_formats/` beside it: `catalog.dcat.jsonld` is what a DCAT harvester (CKAN, which NDP
     runs on) reads.
   - **Register**: Catalog → Add to Catalog.
   - **Approval**
   - **Attach**: your lab's workspace and catalogue, if it has them.
   - **Publish**
3. **Keep everything private (`private: true`) while the work is unpublished.** Check with the PI
   before anything goes public.

---

## 8. Optional: MDDB deposition (MDDB_PIPELINE)

```bash
$SUBSTRATE_PYTHON $K/MDDB_PIPELINE/mddb_build.py "$F"            # needs Step 1 done first
cd "$F/.md_pipeline/mddb/work/<dataset>"
conda activate mwf_env                                           # one-time install: $K/MDDB_PIPELINE/install_mwf.sh (~10 GB)
mwf run -dir . -inp inputs.yaml                                  # ONLY from this work directory
```

---

## 9. Rules for Claude (or anyone) running this file

1. Run the preflight (§2) first. Stop on any FAIL and tell the person why.
2. Run the steps **in order**. Never skip Step 1, because everything else reads its glossary.
3. **Never invent** `identifier`, `topic`, `title`, `construct`, authors, citations, licences,
   `access_level` or `ai_training`. Ask the person.
   Measured values are filled in by the tools; leave them alone.
4. Before editing `settings.yaml`, show the person the exact lines you will change. The tools
   keep timestamped backups in `.md_pipeline/backups/`.
5. **Never** delete, move or rename anything in the data folder. Never write outside
   `MD_RESULTS/` and `.md_pipeline/`.
6. **Never** register, publish or upload to NDP/MDDB/Zenodo/GitHub. Registration is a person's job (§7).
7. Never run MD_COLLECT without `--dry-run` unless the person owns the collection and asked for it.
8. Before saying you are done, run `ndp_build.py "$F" --check`, then report:
   - records ready / incomplete
   - what is still missing
   - where the pages are
   - the next human step

---

## 10. If something goes wrong

| It says | It means / do this |
|---|---|
| `Refused: this folder contains N folder(s) that have already been built` | you are one directory too high; point at one of the folders it lists |
| `Refused: folder: root: says ... but this settings file is in ...` | folder moved or settings copied; rerun with `--rebase-root` |
| `... cannot run this. Name one that can: export SUBSTRATE_PYTHON=...` | your Python lacks pandas/pyarrow/yaml; see §1 |
| trust-check lines, then exit 2 | the kit's files (or a folder above them) are writable by others; `chmod -R go-w $K`, then rerun |
| `N top-level directories are neither indexed nor pruned` | add each to `corpora:` (index) or `prune:` (ignore) |
| a record appeared for something that is not a simulation | **delete it from `corpora:`**; `prune:` alone does not remove a declared corpus |
| `N MD(s) flagged short_of_nominal` | set `production_lengths_ns:` to the length a full run reaches here |
| `N MD(s) flagged stage_assumed` | add a `stages:` rule saying what those files are |
| `SKIPPED: n replica(s) have no frame spacing` | add `frame_interval_ps:` for that system |
| `N MD(s) flagged no_topology` | trajectory atom count matches no topology; a file there is probably not from this run |
| record is `incomplete` | a human-only field is empty in `settings.yaml` |
| `NOTE: ... changed in settings.yaml since the glossary was built` | rerun the whole pipeline, `ndp_pipeline.sh "$F" --yes` |
| `N record(s) carry a SUGGESTED identifier/topic/tags` | confirm them in `settings.yaml`; they are listed in `NDP_registration.md` |
| `the MDs disagree (1000, 400)` + `replica lengths are very uneven` | if the replicas were run to two lengths on purpose, add `nominal_lengths_ns: [1000, 400]` to that system |
| `NOTE: ... config_declared` | no MD config beside the trajectories; the values you typed are used. Fine |
| six files `name_1 … name_6` became one long MD | their starts differ or can't be read; declare `folder: replicas:` |
| you edited `settings.yaml` and nothing changed | a later empty duplicate key overrides yours; edit the existing key |
| `nothing to rebuild from: run the stages first` | `--pages-only` before any analysis ran |
| `WARNING: ... holds N frames, the segments hold M` | the stripped (`STRIP*/merged.dcd`) copy and the raw segments disagree, e.g. restart frames kept twice or a partial strip. The path in the warning names the replica directory. Ask the person who made the strip; do not re-strip on their behalf |

More detail:

- `NDP_instructions.md` (this folder)
- `$K/README.md` and `$K/examples/README.md`
- the README of each kit
- `$F/MD_RESULTS/README.md`
