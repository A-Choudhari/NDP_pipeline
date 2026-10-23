# NDP registration sheet

Generated 2026-10-01T14:37:02 by ndp_build.py. 1 record(s): 0 ready, 1 incomplete.

One section per simulated system, in the order of the form at nationaldataplatform.org (Catalog, Add to Catalog). Visibility: **Unlisted** until the supporting paper is published. Asset type: **Dataset**. Publisher: **Individual** until NDP lists an organization for the Amaro Lab. (Set `ndp_organization:` in lab.yaml once it does.)

Missing values are filled in `/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs/MD_RESULTS/settings.yaml`: the `ndp:` section for the record, the `datasets:` section for contact, licence and citation. Data owner, access level, AI training and the storage name fall back to the lab's own values in `/net/gpfs-amarolab/maalbani/AMARO_MD_PIPELINE_release/AMARO_MD_PIPELINE/lab.yaml`. Then re-run ndp_build.py.

**`base_url` is not set** under `ndp: defaults:` in settings.yaml, so the item URLs below are relative. Host the records first, then set it.

| record | status | replicas | total ns | missing |
|---|---|---|---|---|
| `MDs/70c12cdf` | incomplete | 5 | 5250.0 | access_level |

## Suggested by the pipeline: confirm before registering

These values were guessed from folder names, not written by a person. Type the value you want on the matching key of that record's block under `ndp: systems:` in `/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs/MD_RESULTS/settings.yaml`, then re-run ndp_build.py.

| record | field | suggested value |
|---|---|---|
| `MDs/70c12cdf` | `identifier` | AMARO-HIV-MDS-001 |
| `MDs/70c12cdf` | `tags` | Full, length, HIV, Env, Closed, 5.2 us |

## MDs/70c12cdf — incomplete

Missing: `access_level`

- Note: config_declared: 5 MD(s) have no MD config (no_config); ensemble, temperature and timestep are taken from settings.yaml

| form field | value |
|---|---|
| Title | HIV - Env Closed State - Full-Length Membrane-Embedded - 5.25 µs — all-atom MD, 5 × 1050 ns |
| Description | Molecular dynamics simulations of Full-length HIV-1 clade A BG505 Env trimer (gp160), three protomers, glycosylated and palmitoylated, embedded in an asymmetric lipid bilayer mimicking the composition of the viral budding site. Ectodomain HXB2 residues 31-662, with the SOSIP mutations and the furin cleavage site reverted to the native clade A sequence (UniProt Q2N0S6); MPER, transmembrane domain and cytoplasmic tail included; the two cytoplasmic-tail cysteines C764 and C837 palmitoylated. (256,522 atoms). Starting structure: Composite model. Ectodomain built on the BG505 SOSIP.664 crystal structure (PDB 5T3Z, 3.5 A), with the missing gaps modelled as disordered loops in Modeller 10.3. MPER-TMD-CT built by i-TASSER homology modelling on the clade D NMR structure PDB 7LOI, carrying the clade A sequence. The two halves were joined through the HR2 helix (residues 637-668) re-modelled with AlphaFold2, and the junction relaxed with autoIMD in VMD. The finished MPER was checked against the crystal structure PDB 4G6F.. Glycan profile: Heterogeneous, site-specific profile derived from glycoanalytic data for clade A BG505 viral Env: 23 N-linked sites glycosylated per monomer, not uniform across protomers. N185e, N185h, N618 and N625 are left unoccupied, as they have no consistent evidence of occupancy in clade A BG505. Glycans were built with Glycan Reader & Modeler in CHARMM-GUI, preserving the GlcNAc / GlcNAc2 / ManGlcNAc2 cores resolved in the cryo-EM structures as a scaffold where feasible. Per-chain detail: Supplementary Tables 1-3 of the related publication.. Setup: force field CHARMM36m, TIP3P water, 150 mM NaCl, NAMD3. Conditions: NPT, 310 K, 1.01325 bar, 2 fs timestep. Sampling: 5 × 1050 ns, one frame every 100 ps. Related publication: https://doi.org/10.1038/s41467-026-71874-2. This entry is metadata only. Trajectories stored on lab storage; available on request. Contact: maalbani@ucsd.edu. |
| Keywords | molecular dynamics, Amaro Lab |
| Contact | maalbani@ucsd.edu |
| Item: Title | MD metadata record (JSON) |
| Item: Description | Machine-readable metadata for this simulation. Trajectories stored on lab storage; available on request. |
| Item: Format / Access | JSON / Direct Download |
| Item: URL | 1-FULL-GP120_5-replicas.json |
| Item: Issue Date / Last Updated | 2026-09-22 / 2026-10-01 |
| License | CC-BY-4.0 |
| Data owner (PI) | Rommie E. Amaro (UC San Diego, ramaro@ucsd.edu) |
| Access level | **missing** |
| AI training | by-agreement |
| Raw data location | Amaro Lab GPFS storage, UC San Diego: `/net/gpfs-amarolab/maalbani/ABZYME/1-FULL-GP120/MDs` |
| MDs | `MDs/md_R1`, `MDs/md_R2`, `MDs/md_R3`, `MDs/md_R4`, `MDs/md_R5` |

