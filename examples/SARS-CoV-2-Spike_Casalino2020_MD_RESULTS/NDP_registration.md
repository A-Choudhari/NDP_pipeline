# NDP registration sheet

Generated 2026-10-01T14:40:19 by ndp_build.py. 3 record(s): 0 ready, 3 incomplete.

One section per simulated system, in the order of the form at nationaldataplatform.org (Catalog, Add to Catalog). Visibility: **Unlisted** until the supporting paper is published. Asset type: **Dataset**. Publisher: **Individual** until NDP lists an organization for the Amaro Lab. (Set `ndp_organization:` in lab.yaml once it does.)

Missing values are filled in `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/MD_RESULTS/settings.yaml`: the `ndp:` section for the record, the `datasets:` section for contact, licence and citation. Data owner, access level, AI training and the storage name fall back to the lab's own values in `/net/gpfs-amarolab/maalbani/AMARO_MD_PIPELINE_release/AMARO_MD_PIPELINE/lab.yaml`. Then re-run ndp_build.py.

**`base_url` is not set** under `ndp: defaults:` in settings.yaml, so the item URLs below are relative. Host the records first, then set it.

| record | status | replicas | total ns | missing |
|---|---|---|---|---|
| `CLOSED/d8dbb86a` | incomplete | 3 | 1731.4 | access_level |
| `MUTANT/783dd3bb` | incomplete | 6 | 4291.2 | access_level |
| `OPEN/e0812b9a` | incomplete | 6 | 4222.2 | access_level |

## CLOSED/d8dbb86a — incomplete

Missing: `access_level`

- Note: config_declared: 3 MD(s) have no MD config (no_config); ensemble, temperature and timestep are taken from settings.yaml

| form field | value |
|---|---|
| Title | SARS-CoV-2 full-length spike glycoprotein, closed state: all-atom MD — all-atom MD, 3 replicas, 1731 ns in total |
| Description | Molecular dynamics simulations of Full-length SARS-CoV-2 spike glycoprotein homotrimer, residues 16-1273 per chain (chains A, B, C; segments S1 16-685, S2 686-1273), in the closed state (all RBDs 'down'). Non-canonical residues: HSD/HSE neutral histidines, CYP palmitoylated cysteine. Shared trajectories hold protein + glycans only and are not aligned; water, ions and the POPC/POPS/POPI/POPE/CHL bilayer were stripped. (72,759 atoms). Starting structure: PDB 6VXX (all RBDs 'down'). Glycan profile: 66 N-glycans (22 per chain: G1-G22 chain A, G25-G46 chain B, G48-G69 chain C) and 4 O-glycans (G23, G24 chain A; G47 chain B; G70 chain C), all on glycan chain G. Residues used: BGLN, AGAN, BGAL, BMAN, AMAN, AFUC, ANE5 (complex, hybrid and oligomannose N-glycans).. Setup: force field CHARMM36 all-atom additive, TIP3P water, 150 mM NaCl (Na+/Cl-), neutralising, NAMD 2.14 (Frontera, TACC). Conditions: NPT, 310 K, 1.01325 bar, 2 fs timestep. Sampling: 3 replicas, 1731 ns in total, one frame every 100 ps. Related publication: https://doi.org/10.1021/acscentsci.0c01056. This entry is metadata only. Trajectories stored on lab storage; available on request. Contact: ramaro@ucsd.edu. |
| Keywords | molecular dynamics, Amaro Lab, SARS-CoV-2, spike glycoprotein, glycan shield, NAMD, CHARMM36 |
| Contact | ramaro@ucsd.edu |
| Item: Title | MD metadata record (JSON) |
| Item: Description | Machine-readable metadata for this simulation. Trajectories stored on lab storage; available on request. |
| Item: Format / Access | JSON / Direct Download |
| Item: URL | spike-fulllength-closed_3x500ns.json |
| Item: Issue Date / Last Updated | 2026-09-29 / 2026-10-01 |
| License | CC-BY-4.0 |
| Data owner (PI) | Rommie E. Amaro (UC San Diego, ramaro@ucsd.edu) |
| Access level | **missing** |
| AI training | by-agreement |
| Raw data location | Amaro Lab GPFS storage, UC San Diego: `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/CLOSED` |
| MDs | `CLOSED/CLOSED_R1`, `CLOSED/CLOSED_R2`, `CLOSED/CLOSED_R3` |

## MUTANT/783dd3bb — incomplete

Missing: `access_level`

- WARNING: `nominal_length_ns`: the MDs disagree (1000, 400); left null. Split the system or correct the MD that is wrong.
- WARNING: replica lengths are very uneven (407.0 to 1036.2 ns): check these are all replicas of one simulation and not a test run or a restart stage.
- Note: config_declared: 6 MD(s) have no MD config (no_config); ensemble, temperature and timestep are taken from settings.yaml

| form field | value |
|---|---|
| Title | SARS-CoV-2 full-length spike glycoprotein, open state, N165A/N234A: all-atom MD — all-atom MD, 6 replicas, 4291 ns in total |
| Description | Molecular dynamics simulations of Full-length SARS-CoV-2 spike glycoprotein homotrimer, residues 16-1273 per chain (chains A, B, C; segments S1 16-685, S2 686-1273), in the open state (chain A RBD 'up') bearing the N165A and N234A substitutions. Non-canonical residues: HSD/HSE neutral histidines, CYP palmitoylated cysteine. Shared trajectories hold protein + glycans only and are not aligned; water, ions and the POPC/POPS/POPI/POPE/CHL bilayer were stripped. (71,294 atoms). Starting structure: PDB 6VSB (chain A RBD 'up'), with N165A and N234A substitutions. Glycan profile: 60 N-glycans (20 per chain: G6/G7, G30/G31 and G53/G54 removed with the N165A and N234A substitutions) and 4 O-glycans (G23, G24 chain A; G47 chain B; G70 chain C), all on glycan chain G. Residues used: BGLN, AGAN, BGAL, BMAN, AMAN, AFUC, ANE5.. Setup: force field CHARMM36 all-atom additive, TIP3P water, 150 mM NaCl (Na+/Cl-), neutralising, NAMD 2.14 (Frontera, TACC). Conditions: NPT, 310 K, 1.01325 bar, 2 fs timestep. Sampling: 6 replicas, 4291 ns in total, one frame every 100 ps. Related publication: https://doi.org/10.1021/acscentsci.0c01056. This entry is metadata only. Trajectories stored on lab storage; available on request. Contact: ramaro@ucsd.edu. |
| Keywords | molecular dynamics, Amaro Lab, SARS-CoV-2, spike glycoprotein, glycan shield, NAMD, CHARMM36 |
| Contact | ramaro@ucsd.edu |
| Item: Title | MD metadata record (JSON) |
| Item: Description | Machine-readable metadata for this simulation. Trajectories stored on lab storage; available on request. |
| Item: Format / Access | JSON / Direct Download |
| Item: URL | spike-fulllength-open-N165A-N234A_6x715.2ns.json |
| Item: Issue Date / Last Updated | 2026-09-29 / 2026-10-01 |
| License | CC-BY-4.0 |
| Data owner (PI) | Rommie E. Amaro (UC San Diego, ramaro@ucsd.edu) |
| Access level | **missing** |
| AI training | by-agreement |
| Raw data location | Amaro Lab GPFS storage, UC San Diego: `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/MUTANT` |
| MDs | `MUTANT/MUTANT_R1`, `MUTANT/MUTANT_R2`, `MUTANT/MUTANT_R3`, `MUTANT/MUTANT_R4`, `MUTANT/MUTANT_R5`, `MUTANT/MUTANT_R6` |

## OPEN/e0812b9a — incomplete

Missing: `access_level`

- WARNING: `nominal_length_ns`: the MDs disagree (1000, 400); left null. Split the system or correct the MD that is wrong.
- WARNING: replica lengths are very uneven (404.3 to 1006.3 ns): check these are all replicas of one simulation and not a test run or a restart stage.
- Note: config_declared: 6 MD(s) have no MD config (no_config); ensemble, temperature and timestep are taken from settings.yaml

| form field | value |
|---|---|
| Title | SARS-CoV-2 full-length spike glycoprotein, open state: all-atom MD — all-atom MD, 6 replicas, 4222 ns in total |
| Description | Molecular dynamics simulations of Full-length SARS-CoV-2 spike glycoprotein homotrimer, residues 16-1273 per chain (chains A, B, C; segments S1 16-685, S2 686-1273), in the open state (chain A RBD 'up'). Non-canonical residues: HSD/HSE neutral histidines, CYP palmitoylated cysteine. Shared trajectories hold protein + glycans only and are not aligned; water, ions and the POPC/POPS/POPI/POPE/CHL bilayer were stripped. (72,759 atoms). Starting structure: PDB 6VSB (chain A RBD 'up'). Glycan profile: 66 N-glycans (22 per chain: G1-G22 chain A, G25-G46 chain B, G48-G69 chain C) and 4 O-glycans (G23, G24 chain A; G47 chain B; G70 chain C), all on glycan chain G. Residues used: BGLN, AGAN, BGAL, BMAN, AMAN, AFUC, ANE5 (complex, hybrid and oligomannose N-glycans).. Setup: force field CHARMM36 all-atom additive, TIP3P water, 150 mM NaCl (Na+/Cl-), neutralising, NAMD 2.14 (Frontera, TACC). Conditions: NPT, 310 K, 1.01325 bar, 2 fs timestep. Sampling: 6 replicas, 4222 ns in total, one frame every 100 ps. Related publication: https://doi.org/10.1021/acscentsci.0c01056. This entry is metadata only. Trajectories stored on lab storage; available on request. Contact: ramaro@ucsd.edu. |
| Keywords | molecular dynamics, Amaro Lab, SARS-CoV-2, spike glycoprotein, glycan shield, NAMD, CHARMM36 |
| Contact | ramaro@ucsd.edu |
| Item: Title | MD metadata record (JSON) |
| Item: Description | Machine-readable metadata for this simulation. Trajectories stored on lab storage; available on request. |
| Item: Format / Access | JSON / Direct Download |
| Item: URL | spike-fulllength-open_6x703.7ns.json |
| Item: Issue Date / Last Updated | 2026-09-29 / 2026-10-01 |
| License | CC-BY-4.0 |
| Data owner (PI) | Rommie E. Amaro (UC San Diego, ramaro@ucsd.edu) |
| Access level | **missing** |
| AI training | by-agreement |
| Raw data location | Amaro Lab GPFS storage, UC San Diego: `/net/gpfs-amarolab/lcasalino/CORONAVIRUS/SPIKE_CLEAVED/TO_SHARE_FOR_PUBLICATION/OPEN` |
| MDs | `OPEN/OPEN_R1`, `OPEN/OPEN_R2`, `OPEN/OPEN_R3`, `OPEN/OPEN_R4`, `OPEN/OPEN_R5`, `OPEN/OPEN_R6` |

