#!/usr/bin/env python3
"""
engine.py -- read the program that ran each MD, and its version, from its logs.

WHY THIS EXISTS. The MDDB input spec has two separate fields:

    program: NAMD        the software that carried the trajectory
    version: 3.0.2       its version

and warns, for `program`, to "check already existing values for this field to
avoid making duplicates" -- because the value becomes a facet everyone else's
records are grouped by.

The substrate index does not record either. Its `engine` column is a TRAJECTORY
FORMAT FAMILY, not a program: `NAMD/CHARMM` means "something that writes DCD and
reads PSF". Depositing that string would create a program nobody can search for,
matching neither NAMD nor CHARMM, in 11 of 13 datasets here. And `version` was
simply never emitted, though every NAMD and AMBER run writes its own version
into the first lines of its log.

WHAT IT READS. Only the banner every engine prints when it starts:

    NAMD     Info: NAMD 3.0.2 for Linux-x86_64-multicore-CUDA
    AMBER              Amber 22 SANDER                   2022
    GROMACS  GROMACS version: 2024.2
    CHARMM   Chemistry at HARvard Macromolecular Mechanics ... Version 47b1

This is evidence of what actually ran, which is the same standard the rest of
this toolkit holds to: the trajectory header beats the config file, and the log
beats both for identifying the binary. Where no log carries a banner, both
fields are OMITTED -- an absent key means "not established".

WHY NOT THE CONFIG. A config file says what was intended and is routinely copied
between runs; `mdconfigs.parquet:config_engine` is a good fallback for the
program NAME (a NAMD config can only have been run by NAMD) but it can never
give a version. So: log banner first, config engine second, nothing third.
"""

import os
import re

#: Bytes read from the head of each candidate log. Every engine prints its
#: banner in the first lines; reading further would be paying for nothing.
HEAD_BYTES = 64 * 1024

#: At most this many candidate files per run directory.
MAX_FILES = 8

#: Extensions worth opening at all.
LOG_EXTS = (".log", ".out", ".mdout")

#: Basenames that never carry an engine banner, however they are spelled. These
#: are analysis tools that also write `.log`; opening them is wasted work.
SKIP_STEMS = ("cpptraj", "ambpdb", "antechamber", "parmed", "mmpbsa",
              "vmd", "catdcd", "result")

#: (program, regex) in priority order. The first banner found wins. Each regex
#: must capture the version in group 1.
BANNERS = (
    # `Info: NAMD 3.0.2 for Linux-x86_64-multicore-CUDA`
    ("NAMD", re.compile(r"^\s*Info:\s*NAMD\s+(\S+)\s+for\s", re.M)),
    # `          Amber 22 SANDER                              2022`
    ("AMBER", re.compile(r"^\s*Amber\s+(\d+)\s+[A-Z]", re.M)),
    # `GROMACS version:    2024.2`
    ("GROMACS", re.compile(r"^\s*GROMACS version:\s*(\S+)", re.M)),
    # CHARMM prints its version in the banner block
    ("CHARMM", re.compile(r"Chemistry at HARvard Macromolecular Mechanics"
                          r"(?:.|\n){0,400}?Version\s+(\S+)", re.I)),
)

#: What a substrate `engine` label may be trusted to say about the PROGRAM.
#: `NAMD/CHARMM` is a format family and deliberately maps to nothing: it is the
#: exact ambiguity this module exists to resolve.
ENGINE_PROGRAM = {"AMBER": "AMBER", "GROMACS": "GROMACS", "NAMD": "NAMD",
                  "CHARMM": "CHARMM"}


def _candidates(run_dir):
    """Log files in a run directory, production-looking ones first.

    A production log is the one that describes the trajectory being deposited.
    An equilibration log usually names the same binary, so it is an acceptable
    fallback, but it is read second so that a run whose stages used different
    builds reports the one that produced the data.
    """
    try:
        names = os.listdir(run_dir)
    except OSError:
        return []
    files = []
    for n in names:
        low = n.lower()
        if not low.endswith(LOG_EXTS):
            continue
        if any(low.startswith(s) for s in SKIP_STEMS):
            continue
        path = os.path.join(run_dir, n)
        if not os.path.isfile(path):
            continue
        # production first, then anything else, alphabetically within each rank.
        # The trajectory classifier decides, by whole words: a bare
        # startswith("pr") ranked `protein_min.log` as production -- the bug
        # trajscan fixed for trajectories, back again for logs.
        from .trajscan import stage_of                       # noqa: PLC0415
        rank = 0 if (stage_of(n) == "production" or low.startswith("driver")) else 1
        files.append((rank, n, path))
    files.sort()
    return [p for _, _, p in files[:MAX_FILES]]


def read_banner(path):
    """(program, version) from one log file, or (None, None)."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(HEAD_BYTES).decode("utf-8", "ignore")
    except OSError:
        return None, None
    for program, pattern in BANNERS:
        m = pattern.search(head)
        if m:
            return program, m.group(1).strip()
    return None, None


def of_run(root, rel_path, fallback_engine=None):
    """Identify the program and version behind one run directory.

    Returns a dict with `program`, `version`, `evidence` and `source`:

        source = 'log'      both fields read from a banner; `evidence` names
                            the file, relative to `root`.
        source = 'config'   only the program is known, from the config the
                            scan attributed to this run. No version.
        source = None       nothing established; both fields are None.
    """
    run_dir = os.path.join(root, rel_path)
    for path in _candidates(run_dir):
        program, version = read_banner(path)
        if program:
            return {"program": program, "version": version,
                    "evidence": os.path.relpath(path, root), "source": "log"}
    program = ENGINE_PROGRAM.get(str(fallback_engine or ""))
    if program:
        return {"program": program, "version": None,
                "evidence": "config file engine", "source": "config"}
    return {"program": None, "version": None, "evidence": None, "source": None}


def scan(root, rel_paths, config_engines=None, verbose=False):
    """Identify program and version for many runs. Returns a DataFrame.

    `config_engines` maps rel_path -> the engine the config scan attributed to
    it, used only when no log carries a banner.
    """
    import pandas as pd
    config_engines = config_engines or {}
    rows = []
    for i, rel in enumerate(rel_paths, 1):
        if verbose and i % 100 == 0:
            print(f"        ... {i}/{len(rel_paths)} run dirs", flush=True)
        rec = of_run(root, rel, config_engines.get(rel))
        rec["rel_path"] = rel
        rows.append(rec)
    return pd.DataFrame(rows, columns=["rel_path", "program", "version",
                                       "evidence", "source"])
