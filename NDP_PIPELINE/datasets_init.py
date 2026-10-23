#!/usr/bin/env python3
"""
datasets_init.py -- propose this folder's datasets, into settings.yaml.

    python3 datasets_init.py <folder> [--yes] [--force]

Reads that folder's SUBSTRATE registry and groups production runs into
candidate datasets -- one per declared `system` label, or per common parent
directory where nothing is labelled. It PROPOSES; it never assigns. The output
is a scaffold in which every judgement field is present but empty, because the
fields that matter for deposition (authors, licence, citation, PDB ids,
ligands, interactions) are exactly the ones a machine must not invent.

Writes only the `datasets:` section of MD_RESULTS/settings.yaml, and only
blocks for datasets it does not already list. Nothing else is touched.

EXIT CODES. The wrapper has to tell "this folder legitimately holds no
production runs" apart from "this failed", and both used to return 1. A caller
that cannot distinguish them reports a crash as a clean, empty success.

    0  dataset blocks written
    1  no eligible production runs -- a legitimate, empty result
    2  precondition failure: no substrate, or an unreadable registry
    3  declined: the file exists without --force, or the user said no
"""

#: See the EXIT CODES block above. Named so callers never test bare integers.
EXIT_OK = 0
EXIT_NO_RUNS = 1
EXIT_PRECONDITION = 2
EXIT_DECLINED = 3

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def resolve(folder):
    """(data folder, index directory, deposition directory).

    The names live in lib/layout.py; this stays because every script here calls
    it. Given any of a folder's own output directories -- MD_RESULTS/,
    .md_pipeline/, or the older SUBSTRATE/ and MDDB_substrate/ -- it answers for
    the data folder they belong to.
    """
    from core import layout
    base = layout.base(folder)
    return base, layout.index_dir(base), layout.mddb_dir(base)


from core.settings import DATASETS_SECTION as SECTION_HEAD   # noqa: E402,F401


# Attribution and licence are decisions a PERSON takes for THEIR data, so they
# are pre-filled only for an owner who has taken them: the accounts listed under
# `people:` in lab.yaml. Everyone else gets the keys empty, which exports as
# "not established". Pre-filling one list for every folder put one person's
# name and e-mail on every colleague's datasets.
BLANK_ATTRIBUTION = [
    "authors: []            # who should be cited",
    "contact:               # an e-mail that will still work in 5 years",
]
BLANK_LICENCE = [
    "license:               # e.g. Creative Commons Attribution 4.0 International",
    "linkcense:             # the licence URL",
]


def _q(v):
    """A YAML scalar that stays one, whatever punctuation it holds."""
    import json
    return json.dumps(str(v), ensure_ascii=False)


def prefill(owner):
    """(attribution lines, licence lines) for `owner`, from lab.yaml `people:`;
    the blank scaffold for anyone not listed there."""
    from core import lab
    who = lab.people().get(str(owner or ""))
    if not who:
        return BLANK_ATTRIBUTION, BLANK_LICENCE
    authors = who.get("authors") or []
    authors = [str(a) for a in ([authors] if isinstance(authors, str) else authors)]
    attribution = (["authors:"] + [f"  - {_q(a)}" for a in authors]
                   if authors else [BLANK_ATTRIBUTION[0]])
    attribution.append(f"contact: {who['contact']}" if who.get("contact")
                       else BLANK_ATTRIBUTION[1])
    licence = [f"license: {who['license']}" if who.get("license") else BLANK_LICENCE[0],
               f"linkcense: {who['license_url']}" if who.get("license_url")
               else BLANK_LICENCE[1]]
    return attribution, licence


def folder_owner(subst, base):
    """The owner declared in settings.yaml; the folder's owner on disk otherwise."""
    from core import settings
    try:
        data = settings.load(base).get("folder") or {}
        # `unknown` is what an init under cron or an `env -i` shell used to
        # record (no $USER); it is an absence, not a name.
        owner = str(data.get("owner") or "").strip()
        if owner and owner.lower() != "unknown":
            return owner
    except Exception:                                    # noqa: BLE001
        pass
    try:
        import pwd
        return pwd.getpwuid(os.stat(base).st_uid).pw_name
    except (KeyError, OSError):
        return None


ENTRY_COMMENT = """
  # ---------------------------------------------------------------------
  # {label}: {nruns} production runs, {gb:.1f} GB, {us} replicas {reps}
  # ---------------------------------------------------------------------"""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default=os.getcwd())
    ap.add_argument("--yes", action="store_true",
                    help="accept the proposal without prompting")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing datasets.yaml")
    a = ap.parse_args(argv)

    from core.precheck import require  # noqa: E402  (local helper below)
    from core import lab
    base, subst, mddb = resolve(a.folder)
    runs_p = os.path.join(subst, "registry", "runs.parquet")
    if not require(subst, runs_p):
        return EXIT_PRECONDITION
    owner = folder_owner(subst, base)
    attribution, licence = prefill(owner)

    import pandas as pd
    runs = pd.read_parquet(runs_p)
    prod = runs[runs.n_production > 0] if "n_production" in runs.columns else runs
    if not len(prod):
        print(f"\n  No production runs in {subst}.")
        print("  Nothing to deposit. This folder may hold analysis only.\n")
        return EXIT_NO_RUNS

    labelled = prod[prod.system.astype(str) != "unknown"]
    groups, selector = {}, {}
    if len(labelled):
        for label, g in labelled.groupby(labelled.system.astype(str)):
            conds = sorted({str(c) for c in g.condition.astype(str)} - {"unknown"})
            # A "dataset" is one scientific comparison unit, not a campaign. A
            # system spanning several declared conditions -- different docking
            # poses, ligand chemistries, truncations -- is several datasets that
            # happen to share a protein, and one `mds` list of 236 replicas
            # mixing them would present them as replicates of each other.
            # Splitting is proposed automatically; merge them back by hand if
            # they really are one thing.
            if len(conds) > 1:
                for cond, sub in g.groupby(g.condition.astype(str)):
                    if str(cond) == "unknown":
                        key = f"{label}-unspecified"
                    else:
                        key = f"{label}-{cond}"
                    groups[key] = sub
                    selector[key] = ("system_condition", (label, str(cond)))
            else:
                groups[label] = g
                selector[label] = ("system", label)
    unlabelled = prod[prod.system.astype(str) == "unknown"]
    from core.trajscan import REPLICA_RE
    from core.identity import sanitize

    def top_of(rel):
        # A run lying loose in the folder itself has the path "."; and in a
        # folder that holds only 1R/ 2R/ 3R/, the first component is a REPLICA,
        # which would make one dataset per replica. Both belong to the folder.
        first = str(rel).split("/")[0]
        if first in (".", "") or REPLICA_RE.match(first):
            return sanitize(os.path.basename(base)) or "root"
        return first

    for parent, g in unlabelled.groupby(unlabelled.rel_path.astype(str).map(top_of)):
        # Named for the folder, not `<folder>-unlabelled`: the suffix said
        # something true about the settings and nothing about the data, and it
        # was in every md_key and every record name a person had to read.
        key = parent
        groups[key] = g
        selector[key] = ("paths", sorted(g.rel_path.astype(str)))

    print(f"\n  {len(prod)} production runs in {base}")
    print(f"  proposing {len(groups)} dataset(s):\n")
    for label, g in sorted(groups.items(), key=lambda kv: -kv[1].traj_gb.sum()):
        reps = sorted({str(x) for x in g.replica.astype(str)} - {"unknown", "nan"})
        us = g.prod_ns_header.sum() / 1000.0
        print(f"    {label:28s} {len(g):3d} runs  {g.traj_gb.sum():8.1f} GB  "
              + (f"{us:5.2f} us  " if us else "   --    ")
              + f"replicas={reps[:5] or 'unnamed'}")
    if len(unlabelled):
        print(f"\n  {len(unlabelled)} run(s) carry no `system` label. Declaring one"
              f"\n  under `folder: systems:` in settings.yaml gives better names.")
    if attribution is BLANK_ATTRIBUTION:
        print("\n  authors, contact and license are left empty for you to fill"
              "\n  (or list this account under `people:` in lab.yaml).")

    from core import settings   
    out = settings.path(base)
    if not settings.exists(base):
        print(f"\n  No settings.yaml at {out}. Run substrate_init.py first.\n")
        return EXIT_DECLINED
    declared = {str(d.get("id")) for d in settings.datasets_of(base)}
    if declared and not a.force:
        # Never rewrite a block a person has filled in; only add what is new.
        groups = {k: v for k, v in groups.items() if k not in declared}
        if not groups:
            print(f"\n  every dataset is already declared in {out}\n")
            return EXIT_OK
        print(f"\n  {len(declared)} dataset(s) already declared; adding "
              f"{len(groups)}.")
    if not a.yes:
        try:
            if input("\n  add these to settings.yaml? [y/N]: ").strip().lower() \
                    not in ("y", "yes"):
                print("  nothing written.\n")
                return EXIT_DECLINED
        except EOFError:
            print("  not a terminal; re-run with --yes\n")
            return EXIT_DECLINED

    lines = []
    for label, g in sorted(groups.items(), key=lambda kv: -kv[1].traj_gb.sum()):
        reps = sorted({str(x) for x in g.replica.astype(str)} - {"unknown", "nan"})
        us = g.prod_ns_header.sum() / 1000.0
        lines.append(ENTRY_COMMENT.format(
            label=label, nruns=len(g), gb=g.traj_gb.sum(),
            us=(f"{us:.2f} us," if us else "sampling not readable from headers,"),
            reps=reps[:5] or "unnamed"))
        kind, val = selector[label]
        lines.append(f"  - id: {label}")
        if kind == "system":
            lines.append(f"    system: {label}          # label from `folder: systems:`")
        elif kind == "system_condition":
            sysname, cond = val
            lines.append(f"    system: {sysname}")
            lines.append(f"    condition: {cond}       # one comparison unit")
            lines.append(f"    collections: {sysname}  # keeps the set linked")
        else:
            lines.append("    paths:")
            # Every path: the selector IS the dataset. A cap here (it was 200)
            # silently dropped the remaining runs from the glossary and records.
            for p in val:
                lines.append(f"      - {p}")
        lines.append(f"    name: {label}")
        lines.append("    description: ")
        # Pre-filled only for an owner who has decided them (lab.yaml `people:`);
        # otherwise blank, and so absent from the export until filled.
        for line in attribution:
            lines.append("    " + line)
        lines.append("    type: trajectory        # trajectory | ensemble")
        lines.append("    method: Classical MD")
        for line in licence:
            lines.append("    " + line)
        lines.append("    citation: ")
        lines.append("    thanks: ")
        if kind != "system_condition":
            # Split datasets already carry a collections id above. Emitting the
            # key twice made YAML take the LAST one -- the empty one -- so the
            # link between related datasets was silently dropped.
            lines.append("    collections: ")
        lines.append("    groups:                 # author group, e.g. "
                     + (lab.short_name() or "Smith Lab"))
        lines.append("    pdb_ids: []")
        lines.append("    forced_references: []   # UniProt; 'noref' for constructs")
        lines.append("    ligands: []             # {pubchem: 1986} etc.")
        # AMBER declares the water model in the tleap script beside the
        # trajectory, so the config scan fills it. NAMD never declares it --
        # the only clue is a parameter FILENAME, and reading TIP3P out of
        # `toppar_water_ions_namd.update.str` is an inference about the
        # science. For NAMD this stays empty unless a person fills it.
        lines.append("    wat:                    # water model - AMBER only; "
                     "set by hand for NAMD")
        lines.append("    multimeric:             # monomer | dimer | trimer ...")
        lines.append("    interactions: []        # {name, agent_1, selection_1, ...}")
        lines.append("    pbc_selection: auto     # verify for membrane systems")
        # The federation's questions, per dataset when they differ from the
        # lab's (lab.yaml). Empty falls through to `ndp: defaults:`, then lab.yaml.
        lines.append("    data_owner:             # the PI who owns it; empty = lab.yaml")
        lines.append("    access_level:           # public | private | licensable")
        lines.append("    ai_training:            # allowed | prohibited | by-agreement "
                     "(apart from the licence)")
        lines.append("")
    os.makedirs(mddb, exist_ok=True)
    settings.add_block(base, "datasets", "\n".join(lines))
    print(f"\n  wrote {len(groups)} dataset block(s) into {out}")
    print(f"\n  Next: fill in the empty fields, then\n"
          f"      python3 {os.path.join(HERE, 'glossary_build.py')} {base}\n")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
