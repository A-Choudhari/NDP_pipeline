# Examples: two finished, published MD_RESULTS folders

Each folder is what the pipeline leaves beside a project's simulations, copied as it was on
30 Sep 2026 and rebuilt with this version on 1 Oct 2026 (new record fields, `standard_formats/`). Trajectories are **not** included, and neither is `production_files/` (symlinks to
them). Large `.tsv` tables were dropped where a `.parquet` file holds the same rows. Paths inside
the files still name the storage where the data lives. `lab.yaml` here is the lab file they
were built with: a filled-in model for your own.

## SARS-CoV-2-Spike_Casalino2020_MD_RESULTS/

Casalino L. et al., *Beyond Shielding: The Roles of Glycans in the SARS-CoV-2 Spike Protein*,
ACS Cent. Sci. 2020, 6 (10), 1722–1734. https://doi.org/10.1021/acscentsci.0c01056

| system | replicas | lengths (ns) | total |
|---|---:|---|---:|
| open | 6 | 1000.5, 1000.3, 1006.3, 404.3, 404.3, 406.5 | 4.22 µs |
| N165A/N234A mutant | 6 | 1001.3, 1036.2, 1018.4, 411.7, 407.0, 416.6 | 4.29 µs |
| closed | 3 | 543.5, 573.8, 614.1 | 1.73 µs |

The simulations are protein and glycans only (water, ions and membrane were stripped before
sharing), run in NAMD 2.14 with CHARMM36 at 310 K and saved one frame every 100 ps.

What it shows:
- **Replicas side by side in one directory:** `…_1.dcd … _6.dcd` are read as six
  replicas because every file starts at step 0.
- **Replicas run to two lengths on purpose:** 3 × ~1 µs conventional MD plus 3 × ~400 ns adaptive
  sampling. With `nominal_lengths_ns: [1000, 400]` in a system's block, the length warnings
  become a note.
- **No NAMD configs archived:** setup values are typed in `settings.yaml` (`config_declared`).
- **A fully curated `settings.yaml`:** identifiers, titles, authors, licence and paper are all filled in.

Analysis pages: `analyses/<record>/page_story.html`. The mutant's analysis was still running when
this copy was taken, so its folder has data but no page yet. The `*.manifest.json` files are not
included.

## HIV-Env-Closed-FullLength_2026_MD_RESULTS/

*N-Glycans modulate tilting of HIV-1 envelope glycoprotein*, Nat. Commun. 2026.
https://doi.org/10.1038/s41467-026-71874-2

One system: full-length BG505 Env trimer (gp160), glycosylated and palmitoylated, in an asymmetric
lipid bilayer. 256,522 atoms, 5 × 1050 ns, one frame every 100 ps.

What it shows:
- **A record with every field written out:** `1-FULL-GP120_5-replicas.json`.
- **Three pages:**
  - `page_story.html`, for readers;
  - `page_plotly.html`, the working view;
  - `page_paper_story.html`, which places the article's own figures (`figures_paper/`, CC-BY-4.0,
    made with `MD_ANALYSES/import_paper.py`) beside the measurements.
- **Your own renders on the page:** `figures_user/` holds VMD images, and the page shows them.
- **A guessed identifier:** the identifier `…-HIV-MDS-001` was suggested from the folder name `MDs` and never
  confirmed. The pipeline now lists such guesses in `NDP_registration.md` under *"Suggested by the
  pipeline: confirm before registering"*.

## Files in every MD_RESULTS folder

| file | what it is |
|---|---|
| `README.md` | a summary written about that project's data |
| `settings.yaml` | the one file a person edits |
| `glossary.md` / `glossary.parquet` | one entry per MD: frames, frame spacing, topology, engine, and how each was established |
| `<record>.json` | the NDP metadata record, one per simulated system |
| `catalog.json` | the records in one file, for a lab website catalogue |
| `NDP_registration.md` | per record, the text to paste into the NDP form, and anything still missing |
| `standard_formats/` | each record as schema.org/DCAT JSON-LD (`dataset.jsonld`) and as an RO-Crate (`ro-crate-metadata.json`), plus `catalog.dcat.jsonld` for the folder |
| `analyses/<record>/` | the pages, `data/` tables, `frames.json` (the time axis), and figures |

## Rights and storage in these records

Both examples carry the four federation fields. `data_owner`, `ai_training: by-agreement` and the
storage name come from `lab.yaml`; `data_location.path` was measured. **`access_level` is
deliberately unanswered**, so every record here reads `incomplete` with `access_level` as the only
missing field: that is what a record looks like before the data owner has decided who may get the
raw files.

