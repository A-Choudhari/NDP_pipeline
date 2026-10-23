# Changes

## Unreleased: correctness and safety fixes (phase 0 of the scale-up plan)

- **Datasets with more than 200 runs kept them all.** `datasets_init.py` wrote at most 200
  entries into a dataset's `paths:` selector, so the remaining runs never reached the glossary
  or the records. Folders built with v1.1 that have such a dataset: check its `paths:` list.
- **New `settings.yaml` files keep their section explanations.** They were imported from module
  names that no longer existed, behind a bare `except`, so `datasets:` and `ndp:` were written
  bare. The text now lives in `core/settings.py`.
- **DCD frame counts.** Files with fixed atoms (CHARMM/NAMD `NAMNF`) were sized as if every frame
  held every atom, read as truncated and undercounted; `NSET=0` with frames on disk gave 0
  frames; a corrupt title-block size is now an error, not a frame count; 4D DCDs are sized right.
- **No more partial files.** Records, `catalog.json`, `ndp_index.json`, the standard formats, the
  generated pages, `runs.parquet`, `glossary.parquet`, `settings.yaml` and the mwf inputs are
  written to a unique temporary name and renamed into place (`core/atomic.py`).
- **One build per folder at a time.** `index_build`, `glossary_build` and `ndp_build` take a
  lock in `.md_pipeline/` and refuse, naming the holder, while another run holds it. A lock left
  by a process that died on the same host is taken over.
- Log ranking no longer treats `protein*.log` as a production log (`startswith("pr")`).
- NetCDF4 (HDF5) trajectories are read through `h5py` when installed, and otherwise say what is
  missing instead of `error:TypeError`. The directory-depth limit now applies on Windows.
  `enrich_pdb.py` prints the report path it actually wrote.
- **Analyses:** manifest inputs are keyed by path, not file name, so `R1/merged.dcd` and
  `R2/merged.dcd` are both recorded. A failed parquet write no longer leaves the previous run's
  parquet beside the new TSV. The rdf, grid and ion stages accept `OW` (GROMACS/TIP4P/SPC) as
  well as `OH2`; a bare `O` is still refused. `CPPTRAJ_TIMEOUT_S` sets a ceiling on every cpptraj
  call; rdf uses the OpenMP build when present, with threads from this machine's core count
  (`CPPTRAJ_OMP_THREADS` overrides). A stage failure now says which earlier stages wrote files.
- **MD_COLLECT:** a prefix clash skips that pair of folders and collects the rest (exit 1, both
  named) instead of refusing the whole run. Planning happens under the lock. `.partial` trees
  from crashed runs are cleared. The lock is an O_EXCL file rather than `flock` (unreliable on
  NFS). Copies are fsynced, and the manifest keeps the full SHA-256; older 16-character digests
  still verify. `INDEX.md` and `index.json` are written atomically.
- **Zero-padded segids no longer break a whole system's analysis.** `detect` wrote glycan
  segids `G001…G020` as the range `segid G1 to G20` (the padding was lost through `int()`), which
  matched no atom, so every stage of that system was skipped. `detect._compress` keeps each
  segid's own text and `selections._expand_range` keeps the width.
- **Tests that never ran now run.** Four test files called `unittest.main()` above their last
  classes, so 26 tests were never executed by `python3 tests/test_x.py`; they all pass. 28 new
  regression tests, each checked to fail on v1.1 and pass now (two exceptions: the walk-depth bug
  shows only on Windows, and one test guards behaviour that must not change). Totals: 172 NDP,
  87 MDDB, 30 MD_COLLECT, 271 MD_ANALYSES = 560 (was 506); 12 skip without cpptraj.

## v1.1 (2026-10-01)

- **For any lab.** The kit is now `MD_NDP_PIPELINE/` and the runbook `RUNBOOK.md` (was
  `Amaro_NDP_data.md`). The guide, the runbook and the READMEs describe no particular lab. The kit
  ships `lab.yaml` empty, to fill in once; `lab.template.yaml` is a clean copy and
  `examples/lab.yaml` the filled-in file the examples were built with. A test checks that no
  value of that lab file appears in the code.

- **`lab.yaml`: one file for everything lab-specific.** Lab name, identifier prefix, catalogue
  topics, directory-name topics, per-person authors and licence, record defaults and the lab's data
  tree used to be spread through the code. They are now read from `lab.yaml` (or the file named by
  `MD_PIPELINE_LAB`). Another lab edits that file and nothing else. `NDP_ID_PREFIX` still works.
- **Four new fields in every record:** `data_owner` (the PI), `access_level`
  (`public` | `private` | `licensable`), `ai_training` (`allowed` | `prohibited` | `by-agreement`,
  separate from the citation licence) and `data_location` (storage, measured path, optional URL;
  each replica's own directory is in `sampling.replicas_detail[].location`). A record is
  `incomplete` until `access_level` and `ai_training` are answered. They also appear in
  `catalog.json` and on the registration sheet.
- **Standard formats:** every build writes `MD_RESULTS/standard_formats/`: schema.org + DCAT
  JSON-LD per record, a DCAT catalogue for the folder, and an RO-Crate 1.1 (Process Run Crate)
  per record.
- **A line for every missing field.** Each build gives every field a record still misses an empty
  key in that system's own block of `settings.yaml` (`# missing: fill in for this system`), so the
  place to answer it is always there. Existing text is never rewritten; the insertion is kept only
  if the file reads back to the same values plus the new empty keys. `--check` adds nothing.
- New `settings.yaml` files leave `data_owner`, `access_level`, `ai_training` and `data_location`
  empty under `ndp: defaults:` (the comment shows lab.yaml's value), so they keep following
  `lab.yaml`. An empty key or empty block at any level falls through to the level below.
- Fixes from a review before release: an empty block added to a system no longer hides its
  dataset's values; `data_location` merges key by key across levels; a `lab.yaml` section left
  empty, or a single value where a list is expected, no longer crashes; keywords with numbers or
  `:` are written safely; standard-format ids are unique per record and per replica; Windows line
  endings and multi-line text in `settings.yaml` are kept; `install_mwf.sh` still refuses the data
  tree when Python cannot read `lab.yaml`; `MD_COLLECT` reads the lab name from any YAML form.
- `MDDB_PIPELINE/core/` brought back in step with `NDP_PIPELINE/core/` (`sync_core.sh --apply`).
- Examples rebuilt with this version.
- Tests: 135 NDP (33 new), 87 MDDB, 25 MD_COLLECT (2 new), 259 MD_ANALYSES. The NDP tests run
  as a fictional lab (`NDP_PIPELINE/tests/lab_fixture.yaml`).

**Upgrading a folder built with v1.0:** run `ndp_build.py` on it (no re-index needed). Records
gain the new fields and `standard_formats/`, and stay `incomplete` until `access_level` is set.

## v1.0 (2026-09-30)

First release.
