# MDDB_PIPELINE — deposition files for the MDDB workflow

Turns the MD glossary the [NDP pipeline](../NDP_PIPELINE/README.md) built into
what `mwf` (the MDDB workflow) consumes: one `inputs.yaml` per dataset, and a
working area made of read-only links so the workflow never writes into your
data.

```bash
export SUBSTRATE_PYTHON=$(which python3)
$KIT/NDP_PIPELINE/ndp_pipeline.sh <folder> --yes   # first
python3 $KIT/MDDB_PIPELINE/mddb_build.py <folder>  # then
```

## What it writes

```
<folder>/.md_pipeline/mddb/
  exports/<dataset>/inputs.yaml      the mwf input, pristine -- never run in place
  exports/<dataset>/provenance.yaml  where every input really lives
  work/<dataset>/inputs.yaml         the copy mwf is allowed to rewrite
  work/<dataset>/<md_id>/src.*       read-only links to your trajectories
```

Nothing visible changes in `MD_RESULTS/`: the glossary, the records and
`settings.yaml` belong to the other pipeline, and this one only reads them.

## Why the working area exists

`mwf` treats each MD's `mdir` as scratch: it creates it, writes
`structure.pdb`, `trajectory.xtc`, `metadata.json`, an `inpro/` directory and
every `mda.*.json` into it, and at the end deletes `time.txt`, `#*` and
`f*.trajectory.xtc` from it. Pointing it at a run directory would mean writing
into, and deleting from, raw simulation data. So `mdir` is always a directory of
symlinks named `src.<original name>`, and the safety checks refuse anything that
resolves outside `work/` or inside your data.

**Run mwf only from there:**

```bash
cd <folder>/.md_pipeline/mddb/work/<dataset>
mwf run -dir . -inp inputs.yaml
```

`install_mwf.sh` installs the workflow and its conda environment in your home
(no sudo, about 10 GB).

## What goes into inputs.yaml

- From `settings.yaml` (`datasets:`): name, description, authors, contact,
  licence, citation, PDB entries, references, ligands, interactions. The
  NDP-only keys in that block (`data_owner`, `access_level`, `ai_training`) are
  not part of the mwf spec and are not written.
- From the glossary: program and version, temperature, timestep, ensemble,
  force field, water model, box type, frame spacing, the topology, and each
  MD's trajectory segments **in time order**.

A value every MD agrees on is written once, for the project; where the MDs
differ, it is written on each `mds` entry, which the spec allows. A field with
no evidence is left out rather than written as null — mwf reads a null as
"measured and found to be nothing".

MDs the glossary marks `blocked` (their topology and trajectory do not fit
together) are left out and named, because mwf would reject the project for
them.

## Sharing code with the NDP pipeline

`core/` is a copy of `NDP_PIPELINE/core/` plus `mwf.py` and `workdir.py`, which
only this pipeline needs. Keep it in step after changing the other one:

```bash
./sync_core.sh          # shows the diff, then copies with --apply
```

## Tests

```bash
python3 tests/test_repairs.py            # 72 regressions, one per defect that
python3 tests/test_default_production.py # ever reached an export
```
