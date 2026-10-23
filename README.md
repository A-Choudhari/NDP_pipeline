# MD publishing pipeline: from finished simulations to NDP records

These tools take a folder of **finished molecular dynamics simulations** and turn it into two things:

1. **One metadata record (JSON) per simulated system.** It is ready to register on the
   [National Data Platform (NDP)](https://nationaldataplatform.org) or to list in a lab catalogue,
   and is written again in standard formats (schema.org/DCAT JSON-LD and RO-Crate) so any
   catalogue can read it.
2. **Analyses and interactive web pages** for each system: RMSD, RMSF, Rg, contacts, binding
   energy, water and ions. On a curve you can click any point and get the exact frame.

The tools **never write to your data**. Everything they produce goes into a new `MD_RESULTS/`
folder beside your trajectories, plus a hidden `.md_pipeline/` folder.

Any lab can use it: everything lab-specific is in one file, `lab.yaml`.

Authors: Marco Albani and Rommie E. Amaro, University of California San Diego (2026).
Contact: Marco Albani, maalbani@ucsd.edu.

---

## What is in this archive

| Folder | What it does | Needs |
|---|---|---|
| `NDP_PIPELINE/` | **Step 1–2.** Indexes the trajectories, builds a glossary with one entry per MD, writes one NDP record per system | Python ≥ 3.9 with `pandas`, `pyarrow`, `pyyaml` (`scipy` for AMBER NetCDF) |
| `MD_ANALYSES/` | **Step 3.** Analyses and web pages, read from the glossary | conda env from `MD_ANALYSES/env/environment.yml`; `cpptraj` (AmberTools) for some stages |
| `MD_COLLECT/` | **Step 4 (optional).** Copies finished `MD_RESULTS/` folders into one shared collection | Python 3 standard library |
| `MDDB_PIPELINE/` | **Optional.** Deposition inputs for the [MDDB](https://mddbr.eu) `mwf` workflow | `mwf` (`MDDB_PIPELINE/install_mwf.sh`) |
| `examples/` | Two complete, **published** outputs to look at before running anything | none |
| `lab.yaml` | **Everything lab-specific**: lab name, identifier prefix, topics, people, data owner, storage. **Fill it in first**; it is the only file a lab edits. `lab.template.yaml` is a clean copy; `examples/lab.yaml` a filled-in one | |
| `RUNBOOK.md` | **The runbook**: every step, check and command. Give it to Claude Code or follow it by hand | |
| `NDP_instructions.md` / `.pdf` | **The guide**: the same workflow explained for every level, undergraduate to PI | |
| `CHANGELOG.md` | What changed in each version | |
| `docs/COMPARISON.md` | How this compares with other MD metadata tools and schemas | |
| `requirements-ndp.txt` | Python packages for NDP_PIPELINE and MD_COLLECT | |

## Quick start

```bash
unzip MD_NDP_PIPELINE_v1.1_2026-10-01.zip && cd MD_NDP_PIPELINE
export KIT=$PWD
#    edit $KIT/lab.yaml: your lab's name, identifier prefix, topics, data owner, storage

# Python for the metadata side
python3 -m pip install -r requirements-ndp.txt          # or use any env that has them
export SUBSTRATE_PYTHON=$(which python3)

# Python for the analysis side
conda env create -f MD_ANALYSES/env/environment.yml     # creates `mdviz`
export P=$(conda run -n mdviz which python)

# one project folder of FINISHED simulations
export F=/path/to/your/project
$KIT/NDP_PIPELINE/ndp_pipeline.sh "$F" --yes                         # 1. describe
#    edit $F/MD_RESULTS/settings.yaml                                 # 2. curate
$SUBSTRATE_PYTHON $KIT/NDP_PIPELINE/ndp_build.py "$F"                 #    rebuild records
$P $KIT/MD_ANALYSES/analyse.py "$F" --detect-only                    # 3. what is it?
$P $KIT/MD_ANALYSES/analyse.py "$F" --only composition,structure     #    analyse
```

Then open `$F/MD_RESULTS/analyses/<record>/page_story.html` in a browser.
`RUNBOOK.md` gives every step in full, with checks and stop conditions.

## Configuring it for your lab

**Everything that belongs to one lab is in `lab.yaml`**, at the top of the kit. Nothing in the
code names a lab: a lab fills in that one file, once, before the first run, and edits nothing
else. It ships empty (`lab.template.yaml` is the same clean copy, to start over from);
`examples/lab.yaml` is the filled-in file the two examples were built with.

| `lab.yaml` section | What it controls |
|---|---|
| `lab:` | the lab's name (`lab:` in every record, the publisher in the standard formats), short name, catalogue name and URL, NDP organization |
| `identifiers: prefix:` | first word of suggested identifiers (`SMITHLAB-HIV-ENV-001`) |
| `record_defaults:` | what every record falls back to: keywords, data availability, the data owner (PI), access level, AI-training permission, the storage the raw files are on |
| `people:` | per Unix account: authors, contact and licence pre-filled into that person's new dataset blocks |
| `topics:` | the catalogue's research topics: slug, identifier code, the words that file a folder under it |
| `index_topics:` | research areas proposed from directory names when a folder is first indexed |
| `storage: data_roots:` | the lab's data tree, which `install_mwf.sh` refuses to install into |

`MD_PIPELINE_LAB=/path/to/other_lab.yaml` uses another file (one kit, several labs). Without a
lab file nothing is invented: no lab name, no topics, identifiers start with `LAB`, and records
list what is missing. `python3 NDP_PIPELINE/core/lab.py` prints what is in force.

Per machine, set these where your defaults differ:

| Setting | What it controls | Default |
|---|---|---|
| `SUBSTRATE_PYTHON` | Python used by the shell entry points | `python3` |
| `CPPTRAJ` | serial `cpptraj` binary | `cpptraj` on your PATH |
| `CPPTRAJ_OMP`, `AMBER_SH` (or `AMBERHOME`) | optional OpenMP build of cpptraj, and the script to source first | `cpptraj.OMP` on PATH; none |
| `NDP_ID_PREFIX` | overrides `identifiers: prefix:` for one shell | unset |
| `MD_COLLECT_DEST` | where `MD_COLLECT` writes the shared collection | `~/MD_NDP_DATA` |
| `MD_COLLECT/sources.yaml` → `roots:` (or `MD_COLLECT_ROOTS`) | where `collect_pipeline.sh --all` looks | none: pass folders or `--root` |

Each folder's own `settings.yaml` overrides any `record_defaults:` value under `ndp: defaults:`,
per dataset (`datasets:`), or per system. A key left empty falls back to `lab.yaml`.

## What every record says about rights and storage

Besides the science, each record carries what a federated catalogue and a licensing decision need:

| Field | Values | Where it comes from |
|---|---|---|
| `data_owner` | `name`, `email`, `orcid`, `affiliation` of the PI who owns the data | `lab.yaml`, overridable per dataset or system |
| `access_level` | `public` · `private` · `licensable` | **a person decides**: lab-wide in `lab.yaml`, or per folder, dataset or system |
| `ai_training` | `allowed` · `prohibited` · `by-agreement` — separate from the citation licence | **a person decides**, the same way |
| `data_location` | `storage` (which system), `path` (**measured**: the directory every replica sits under), `uri` (where to fetch it, if online) | storage and URI declared; path measured. Each replica's own directory is in `sampling.replicas_detail[].location` |

A record stays `incomplete` until `access_level` and `ai_training` are answered, at whichever
level suits: once in `lab.yaml` for a lab-wide policy, or per folder, dataset or system.
Every field a record still misses gets an empty line in that system's own block of
`settings.yaml` on each build (`# missing: fill in for this system`), so there is always a
place to type the answer for that system.
Note that `data_location.path` is an absolute path on your storage, and it is published with the record.

## Standard formats

Every build also writes `MD_RESULTS/standard_formats/`:

| File | Format | Read by |
|---|---|---|
| `<record>/dataset.jsonld` | one JSON-LD node typed **schema.org `Dataset`** and **DCAT `dcat:Dataset`**, with licence, access rights (EU access-right vocabulary), rights holder, and AI training as an ODRL policy | CKAN/NDP (ckanext-dcat), Google Dataset Search, any DCAT harvester |
| `catalog.dcat.jsonld` | the folder as one `dcat:Catalog` holding every dataset | a harvester pointed at a lab or folder |
| `<record>/ro-crate-metadata.json` | **RO-Crate 1.1**, Process Run Crate 0.5 profile: each replica's production run as a `CreateAction` (engine → trajectory files), and the metadata extraction as another | RO-Crate tools, WorkflowHub, provenance-aware repositories |
| `<record>/record.json` | the record itself, so the crate is self-contained | |

Every value comes from the record; nothing null is written. Set `base_url` under
`ndp: defaults:` and the identifiers become the hosted record's URL. The crates pass the
`rocrate-validator` REQUIRED checks for both `ro-crate-1.1` and `process-run-crate-0.5`.

**Shared installations.** `trust_check.sh` refuses to run if the kit's files, or the folders above
them, can be changed by anyone other than the owner. If several people share one copy, keep it
owned by one account and not group-writable (`chmod -R go-w`).

## Tests

```bash
for t in NDP_PIPELINE/tests/*.py MDDB_PIPELINE/tests/*.py; do $SUBSTRATE_PYTHON "$t"; done
$SUBSTRATE_PYTHON -m pytest MD_COLLECT/tests -q
for t in MD_ANALYSES/tests/*.py; do $P "$t"; done     # cpptraj tests skip unless cpptraj is found
```

The fixtures are synthetic. The tests read no real data and write only to temporary directories,
and run as the fictional lab in `NDP_PIPELINE/tests/lab_fixture.yaml`, so they pass whatever your
`lab.yaml` says. All suites passed on the machine this archive was built on: 135 NDP, 87 MDDB,
25 MD_COLLECT and 259 MD_ANALYSES tests.

## Examples (both published)

- `examples/SARS-CoV-2-Spike_Casalino2020_MD_RESULTS/`: full-length SARS-CoV-2 spike: open,
  closed and the N165A/N234A mutant. 15 replicas, 10.2 µs in total. Casalino et al., *ACS Cent.
  Sci.* 2020, 6, 1722. https://doi.org/10.1021/acscentsci.0c01056
- `examples/HIV-Env-Closed-FullLength_2026_MD_RESULTS/`: full-length, membrane-embedded HIV-1
  Env (BG505) in the closed state. 5 × 1050 ns. *N-Glycans modulate tilting of HIV-1 envelope
  glycoprotein*, Nat. Commun. 2026. https://doi.org/10.1038/s41467-026-71874-2

`examples/README.md` explains what each file is. These folders hold **metadata and analyses only**;
the trajectories are not included.

## Licence and citation

MIT License (see `LICENSE`). The examples describe published work; please cite those papers when
you use them. If these tools help your work, please cite or acknowledge the authors
(Marco Albani and Rommie E. Amaro, UC San Diego).
