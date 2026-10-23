# MD_ANALYSES — analyses that know which frame they came from

Point it at a folder the metadata pipeline has already indexed. It reads the
glossary, works out what each system *is*, runs the analyses that suit it, and
leaves a bundle beside the record: the numbers as tidy tables, and a page where
clicking a point on a curve tells you the frame it came from.

```bash
P=$(conda run -n mdviz which python)
$P $KIT/MD_ANALYSES/analyse.py <folder> --detect-only   # what is this?
$P $KIT/MD_ANALYSES/analyse.py <folder> --system spike-fulllength-open_6x703.7ns
```

It reads your data and never writes to it. Everything it produces is inside
`<folder>/MD_RESULTS/analyses/<record>/`.

## What you get

```
MD_RESULTS/analyses/spike-fulllength-open_6x703.7ns/
  page_story.html       the page for the website: one column, read top to bottom
  page_paper_story.html the same page with the published paper's figures in it
  page_plotly.html      the working page: two dense columns, everything at once
  frames.json           every replica's time axis and where each frame lives
  data/*.tsv *.parquet  one row per measurement, with its stride and spacing
  data/*.manifest.json  what produced it: the spec, the inputs, the hashes
  figures_user/         YOUR pictures and movies — both pages show them
  figures_paper/        the article's figures, if a paper is attached
  paper.yaml            the citation, the licence, and a caption per figure
```

## Changing the words on a page

Three kinds of text, three different places, and only one of them is a file you
edit.

| what | where it comes from |
|---|---|
| the title, and the paragraph under it | `description:` and `title:` for that system in `MD_RESULTS/settings.yaml` |
| everything in the **Record** table at the foot | the same file: identifier, topic, contact, licence, publication, authors, construct, starting structure, glycan profile |
| captions for your own renders | `credits.yaml` beside them in `figures_user/` |
| the section headings and their one-line ledes | `pages/story.py`, the same on every page |
| the sentence under each plot | written from the rows in `data/`, on purpose |

Write a `description:` and it appears under the title. Leave it empty and the
metadata pipeline composes one from the atom count and the setup for the
catalogue listing, and the page shows nothing — a composed sentence at the top
of a published page restates the record table three screens down and reads as
though somebody had chosen to say it. `description_basis` on the record is what
tells the two apart. `**bold**` works, and the text is escaped before the
markup is applied, so a residue name or a path in it cannot open a tag.

**The plot captions are deliberately not yours to edit.** Each one is written
from the same `scope="aggregate"` rows the curve beside it is drawn from, so it
carries the measured number and cannot drift from the plot. That is the whole
point of them. If a caption says something wrong, the row it read is wrong, and
that is worth knowing.

## Two pages, one set of numbers

Both are built from the same rows in `data/` and drawn by the same code, so
they cannot disagree about a value — they differ only in what a reader is
assumed to know.

```bash
$P analyse.py <folder> --pages-only                  # story page only (default)
$P analyse.py <folder> --pages-only --page plotly    # the working page
$P analyse.py <folder> --pages-only --page both
```

`page_story.html` is the published one: the supplied renders first, then one
analysis at a time, each under a sentence that carries its own measured number
— *"the trimer settles at <mean ± sd> Å after about <t> ns and stays there"* —
read straight out of the same `scope="aggregate"` rows the curve is drawn
from, so a caption cannot drift from its plot. A section whose rows do not
exist is not drawn at all: the apo trimer simply has no interface section.

`page_plotly.html` is the working view and is **only** written when you ask
for it by name, so rebuilding the story page can never replace it.

## The paper, on the same page as the measurements

A dataset that has been published has more to show than this kit measures. The
article's figures go on the page beside ours, each carrying its number, its
caption as published, and a link to the article:

```bash
$P import_paper.py <folder> <paper.pdf> --dry-run   # what it found
$P import_paper.py <folder> <paper.pdf>
$P analyse.py <folder> --system <record> --pages-only --page paper
```

It is a **separate file**, `page_paper_story.html`, built by the same builder
from the same rows — so the two pages cannot disagree about a value, and
publishing the combined one never touches the plain one.

**The licence decides whether a figure is copied.** `import_paper.py` reads the
article's own licence statement. Creative Commons or CC0 and the figures are
written into `figures_paper/` and shown; anything else and the captions and the
citation are still imported, but each figure is a link to the article instead of
a copy of it. There is no flag to override that.

**Which section a figure lands in** is `section:` in `paper.yaml` — `system`,
`moves`, `touch`, `solvent`, or `paper` for the article's own section at the
end. The importer proposes one from the caption's own words; you correct it by
editing one line. A figure routed at a section this page does not draw falls
back to the paper's rather than vanishing.

**Two things on the page, told apart.** A figure from the article is tagged with
its number and its source; a section that holds both kinds says `measured here`
on ours. A page that blends a published result into a measurement without
saying so is claiming someone else's work as its output.

**The figures are re-rendered, not lifted.** A PDF may store an image at a
different horizontal and vertical resolution and let the page geometry undo it:
five of the six figures of the first paper attached here are stored at 300 × 150
dpi, so the pixel array is squashed to half height and only the on-page
rectangle says so. The figure's rectangle is re-rendered at one dpi in both
directions, which also keeps the vector labels drawn *on top* of the bitmap —
the panel letters and the leader lines are not in the stored array at all.
`--dpi 200` is the default and is sharp at the page's full width; `--dpi 300`
doubles the page weight for a print-quality one.

**What extraction gets wrong, and says so.** A caption ends at the journal's
source-data sentence or at a blank line that follows a full stop — a column
break lands mid-sentence, which is how the two are told apart. Where neither
fires, `caption_basis` reads `run-on` and that caption is the one to check.
A PDF printed from a browser hides `fi` and `fl` in the Unicode private-use
area, so every extractor drops them and "significant" becomes "signicant"; the
mapping is inferred by scoring each candidate against the system word list, and
`ligature_basis` says whether that succeeded.

## Pictures from a deck

```bash
$P import_slides.py <folder> <deck.pptx> --dry-run
```

Takes the renders out of a PowerPoint and files each one under the record its
slide is about, reading the system from the slide title. Captions are matched
by **position** — a label belongs to the picture whose left edge it shares,
which is how a deck says which is which — and labels stacked at one x are
recognised as a legend for both pictures rather than a caption for each.
`credits.yaml` beside the images records the slide, title and caption every
one came from, is rewritten rather than appended on a second run, and retires
the files a previous run imported.

## Clicking a point

The provenance panel answers the question a curve cannot:

```
frame 2016 · t = 403.4 ns · replica 1-R1
step 201,700,000 · lives in merged.dcd, frame 2019 of 2525
RMSD <value> Å · Rg <value> Å               [copy VMD command]
```

It names **the file you actually open** — the merged trajectory the analysis
ran on — and the frame's index inside it. The time and the step come from the
original `pr*.dcd` headers (`istart`, `nsavc`, `delta`), because cpptraj writes
`merged.dcd` with a placeholder timestep and a timeline built from it would
report a 504 ns run as 0.003 ns. Two sources, one answer, and the `pr.part*`
segment behind it stays out of the way.

## What it decides for you, and what you confirm

`--detect-only` reads the topology and says what the system is made of:

```
spike-fulllength-open_6x703.7ns            (examples/: SARS-CoV-2 spike, Casalino et al. 2020)
  detected : protein, glycans  -> glycoprotein
  system   : 72,759 solute atoms, 0 waters, 0 ions, 595 glycan residues
  receptor : segid AS1 AS2 BS1 BS2 CS1 CS2
  glycans  : segid G1 to G70
  runs     : rmsd, rmsf, rg
```

Components are read from residue names, protein chains are grouped into
entities by their sequence, and — the part a size rule cannot get right — the
PSF's bond list says which chains carry N-glycans. gp41 is 153 residues, the
size of an antibody domain pair, and would have been called a ligand; it
carries two glycans per protomer, and the Fab chains carry none. The proposal
is yours to correct: components select a **bundle**, and a bundle is a list of
stages, so a new kind of system is a table entry in `core/spec.py`, not new code.

## Which files it reads

A replica holds a lot of trajectories and only some of them are the run. The
pairing is made on **atom counts**, never on names:

- the solute trajectory is any `.dcd` under a `STRIP*/` directory — whether it
  is called `merged.dcd` or `<name>_stripped.dcd`;
- its topology is the file whose `!NATOM` matches that trajectory's header,
  looked for in the replica, its parent and grandparent — one solute topology
  often serves three replicas from `MDs/1-INPUTs/merged.psf`;
- the candidate whose frame count matches the run wins, so a half-finished
  strip from an earlier attempt cannot be analysed as if it were the whole
  thing, and the frame count compared against is the run's **unique** frames:
  a restart re-simulates its checkpoint frame, so segments that sum to 10,081
  are a 10,080-frame run and a strip of 10,080 is complete, not short.

`solvated_ions.psf` sits right beside the stripped trajectory in every
replica and has fifteen times its atoms. Nothing here can pick it by accident.

## Selections

Written in VMD/CHARMM syntax, because that is what everyone here types, and
translated with the differences that actually bite:

| you write | it runs | why |
|---|---|---|
| `noh` | `not name H*` | heavy atoms by NAME — these topologies have broken masses |
| `within 5 of X` | `(around 5 (X)) or (X)` | VMD's `within` **includes** X; MDAnalysis's `around` does not |
| `segid GA1 to GA23 GB1 to GB23` | all 46 | every range is expanded, not just the first |
| `water`, `ions`, `glycan` | explicit resname lists | the same atoms a VMD macro would pick |
| `chain A` | refused-ish (warned) | chainID collapses the three antibody copies and every glycan |

Every selection is resolved against the topology **before** a trajectory is
opened, and an empty one is an error — the lab's own TCL scripts carry the
warning that the wrong line "silently gives zero contacts".

## Stages

| stage | reads | gives |
|---|---|---|
| `composition` | topology + NAMD `.xst` | TIP3P waters, ions per species, box edges and volume, measured salt concentration, water density |
| `structure` | the stripped trajectory | RMSD and Rg per frame, RMSF per residue, per subset (receptor / ligand / glycans), fitted on the receptor |
| `interface_contacts` | the stripped trajectory | residue–residue contacts across the interface, and their occupancy |
| `binding_energy` | results already on disk | an end-state estimate per component; no trajectory, so it costs nothing |
| `ions` | the solvated trajectory | ion counts and sites |
| `water_shell`, `water_bridges`, `water_residence` | the solvated trajectory | one walk, three observables: the shell, what it bridges, how long it stays |
| `density_grids` | the solvated trajectory | `.dx` volumes per species, in multiples of bulk |
| `rdf` | the solvated trajectory | radial distribution functions |

`core/spec.py` holds the whole map: `BUNDLES` is what each kind of system
should get, `IMPLEMENTED` is what exists as code, and `--detect-only` prints
both — `runs:` and `planned:` — so the difference is on screen rather than
implied.

Restart seams are dropped: a window that resumed from an earlier checkpoint
re-emits a step with *different* coordinates, and those few frames per
replica would otherwise be a 0.2 ns spike nobody could explain.

**The antibody that flies away.** NAMD runs these with `wrapAll on`, so every
molecule is wrapped into the cell independently and a bound Fab can be written
a box length from the trimer — which ruins a render and would put a ~200 Å step
in its RMSD. Before anything is measured, each group is moved back into the
receptor's periodic image, whole and by segment (never atom by atom, which is
what tears a chain in half), using explicit minimum-image arithmetic — nothing
is ever `wrap()`ped, because `cellOrigin` is not the origin. Each run reports
`segments_reimaged`, so a trajectory that needs a lot of it is visible rather
than silently corrected. For pictures, re-image with
`cpptraj: autoimage anchor <receptor mask>` — which is what the
`STRIP_2026-09-23/strip.in` files here do.

Still to come, in this order: the surface set (`bsa`, `epitope_sasa`,
`glycan_shielding`), the membrane set (`leaflet`, `insertion_depth`), the two
core stages nothing writes yet (`dssp`, `pca`), `glycan_hydration` and
`solvent_coupling`, and the `analyses` block in the record itself.

## Environment

```bash
conda env create -f env/environment.yml && conda activate mdviz
```

One env because the machine's two stacks are split: `mwf_env` has MDAnalysis
2.10 but no plotly or pyarrow; base anaconda has the viz stack but MDAnalysis
2.7 and no netCDF4, and `kaleido` is missing everywhere. cpptraj and VMD stay
external (`cpptraj` from AmberTools -- set `CPPTRAJ` if it is not on your PATH -- and VMD 1.9.4
headless with Tachyon — OptiX is broken on this driver).

## Cost

Reading is the budget: this NFS mount gives ~94 MB/s cold. A stripped 504 ns
replica is 1.7 GB (~18 s); the solvated one is 23 GB (~4 min). So every stage
that can share a walk of the trajectory does, and per-frame results are
persisted rather than recomputed. Random access beats a sequential walk above
about stride 3 (~215 ms per seek against 92 ms per frame), which is why the
solvent stages will read strided rather than walking.
