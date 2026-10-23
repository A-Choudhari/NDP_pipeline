---
title: "Where this sits"
subtitle: "What else exists for publishing molecular dynamics, and how this compares"
date: "2026-09-27"
---

## The short version

Nothing found does all three of what this does — measure metadata off the
files, analyse to a reader-facing page, and collect across a lab — over
private storage, across engines, for simulations nobody instrumented in
advance. But **each stage separately has real prior art**, and one thing we
have described as novel is not.

Read this before the meeting. If someone asks "why didn't you just use X", X
is almost certainly in here.

---

## The closest neighbours

**MDDB — Molecular Dynamics Data Bank.** EU project, the nearest thing
architecturally: normalise raw MD, run quality control and analyses, load to a
database, browse in a viewer with plots linked to structure.
<https://mddbr.eu/> · <https://github.com/mmb-irb/MDDB-workflow>

The important difference is the opposite of what you would guess. **Their
metadata is author-declared.** In the workflow's own inputs template, force
field, water model, ensemble, temperature, timestep and frame step are
optional YAML fields a person types, marked for display. Only the topology,
the trajectory and the list of MDs are mandatory. Automatic derivation is
limited to structural things.

> **Carry this caveat.** The code here refers to an "mwf" format. I could only
> verify `mwf` as the *command name* of MDDB-workflow. I found no published
> specification of a "molecular workflow format". Treat any claim that we emit
> a standard format as unverified until someone produces the spec.

**GROMACS MetaDump.** The single closest published match to our first stage:
automated annotation of GROMACS runs, reading the binary run input, the
topology and the structure, emitting force field, water model, temperature,
pressure, timestep and composition as JSON or YAML. GROMACS only, and it
defines its own schema. *J Cheminform* 2025.
<https://doi.org/10.1186/s13321-025-01082-5>

**This is the one to cite and compare against.** Our difference is breadth:
multiple engines, multiple system kinds, in one record.

**MDverse.** Indexed roughly 250,000 files across Zenodo, Figshare and OSF and
**parsed metadata out of the files** — integrator, timestep, length,
thermostat, temperature from run parameters; composition from structures;
frame counts from trajectories. Same instinct as ours, aimed at public
repositories rather than lab storage. *eLife* 2024.
<https://elifesciences.org/articles/90061>

They state our problem as unsolved, which is worth quoting: it is often
impossible to find even the identity of the molecular system, the temperature,
or the length.

**AmberMeta.** Parses AMBER inputs and outputs, reconstructs the simulation as
phases and steps, validates atom-count agreement and continuity between steps,
and flags gaps as findings. Uncomfortably close to our first stage, for AMBER.
<https://github.com/MicheleBonus/ambermeta>
*Unverified:* no paper or DOI on the repository, and it carries a
business-source licence rather than an open-source one.

**MDRepo.** Deposition with a **hand-written** metadata file; the automatic
part is post-processing. Overlaps our collector and analysis stages, not our
metadata stage. *NAR* 2025 · <https://mdrepo.org>

**BioExcel-COVID-19.** MDDB's prototype. Metadata by deposition form with
mandatory items, plus analyses and streaming with plots linked to snapshots.
*NAR* 2024.

---

## Schema and infrastructure

- **NOMAD / FAIRmat** — a real, published metadata schema with engine parsers
  and workflow provenance, at scale. Materials-science framing. This is the
  strongest "why didn't you emit this?" candidate.
  <https://nomad-lab.eu/>
- **AiiDA / Materials Cloud** — automatic provenance graphs, but captured *at
  run time*. Ours is forensic: it reads simulations that finished years ago and
  were never instrumented. That is a real difference, and it is the one to say
  out loud.
- **MOLSIM ontology** — a biomolecular-simulation vocabulary covering force
  fields, water models, ensembles, thermostats. Alpha, no paper found, but it
  is the vocabulary our JSON keys should map to.
  <https://github.com/CPCLab/molsim-ontology>
- **MoSDeF / TRUE** and **H5MD** — build-time provenance and a self-describing
  trajectory format. Adjacent, not competing.
- *Unverified:* I could not confirm a MolSSI project named "MDInterop". MolSSI's
  MDI is a runtime code-to-code API and MMIC is component I/O schema; neither
  is a deposition metadata schema.

---

## The viewer side, and one claim we should drop

**MDsrv** has offered click-navigable time traces since 2022: clicking a point
on an RMSD or distance trace loads the matching 3D frame. MDposit does the
same, and VIAMD does it on the desktop. *NAR* 2022.
<https://nglviewer.org/mdsrv/>

**So "click a point, get the frame" is not novel, and we should stop implying
it is.** What is unusual is the deployment: a static, dependency-free HTML file
that needs no server, works from a filesystem or an email attachment, and hands
back a VMD command rather than rendering the structure itself. That is a real
difference and an honest one, but it is a difference in *how it is delivered*,
not in *what it can do*.

---

## What actually looks unusual

Tested sceptically, in descending order of confidence.

1. **Per-value provenance.** Every value carries how it was established — read
   from a header, declared by a person, inherited from a template config — and
   the glossary records that per field, per simulation. I found nothing
   MD-specific that does this. The nearest published analogue is general
   guidance on recording raw metadata before structuring it (*Sci Data* 2025,
   <https://arxiv.org/abs/2408.17309>). **Lead with this.**
2. **Null and `incomplete` rather than a guess**, enforced rather than advised.
   MDverse names exactly this as the failure mode nobody has solved. Expect a
   reviewer to point out that it is discipline rather than technology, which is
   fair and does not make it less useful.
3. **Refusing a number the sampling cannot support** — residence times at one
   frame per 20 ns are not approximate, they are meaningless, and the rows say
   so rather than carrying a plausible figure.
4. **Multi-engine, multi-system-kind, post-hoc.** Not novel in any single
   piece, but the combination over storage nobody prepared is where the other
   tools do not go.

---

## What a reviewer would say to reuse

MDAnalysis and MDTraj for parsing. MetaDump for the GROMACS branch. MDDB's
analyses and viewer rather than a new page. MDsrv or nglview for viewing.
RO-Crate and schema.org for packaging. MOLSIM or NOMAD terms instead of our own
JSON keys. CKAN's DCAT extension for the collector index.

The defensible answer is narrow and true: none of them run over finished
simulations that were never instrumented, and none record per-field provenance.

---

## Three things worth emitting alongside our JSON

Cheap, and they make our records interoperable rather than ours alone.

> **Update 2026-10-01:** 1 and 2 are now written for every record, in
> `MD_RESULTS/standard_formats/` (see the README). 3 is still to do.

1. **schema.org `Dataset` / DCAT JSON-LD.** The National Data Platform runs on
   CKAN, which speaks this natively, and it buys indexing in Google Dataset
   Search for nothing.
2. **RO-Crate, Workflow Run profile.** One JSON-LD manifest per record folder.
   Our provenance glossary maps onto PROV-O terms almost directly. **Highest
   value of the three.** <https://doi.org/10.1371/journal.pone.0309210>
3. **A term-mapping table** from our keys to MOLSIM, and to NOMAD's schema
   where one exists.

---

## The context to state plainly

No community metadata standard for molecular dynamics has actually been
adopted. That is said repeatedly across this literature — Abraham et al.,
*Sharing Data from Molecular Simulations*, *JCIM* 2019
(<https://doi.org/10.1021/acs.jcim.9b00665>), the ELIXIR data-management
guidance, and the BioExcel and MDDB community work.

That gap is where this stands. It is also why interoperating with one of the
existing schemas, rather than waiting for a standard, is the sensible next
move.

---

*Compiled 2026-09-27 from a web survey. Items marked unverified could not be
confirmed from a primary source; everything else carries its link.*
