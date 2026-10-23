#!/usr/bin/env python3
"""
index_build.py -- read every trajectory once, and write down what it says.

    python3 index_build.py <folder> [--deep]

One row per run directory: how many files, which stage each is, atom and frame
counts, frame spacing, the topology that matches, the replica label. Nothing is
interpreted here -- that is the glossary's job. Reads your data and writes only
.md_pipeline/index/registry/runs.parquet inside the folder.

    --deep     scan .xtc/ASCII .mdcrd for frame counts (SLOW, off by default)
"""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import atomic, config, layout, settings              # noqa: E402


def find_folder(arg):
    """The data folder, however it was named on the command line.

    A path to the folder itself, or to any of its own output directories --
    MD_RESULTS/, .md_pipeline/, or the SUBSTRATE/ of the layout before them,
    which is what every older note and script passes.
    """
    return layout.base(os.path.abspath(arg) if arg else os.getcwd())


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("folder", nargs="?")
    p.add_argument("--deep", action="store_true")
    p.add_argument("--quiet-period", type=int, default=0, metavar="SECONDS",
                   help="skip sources modified within SECONDS (unattended runs)")
    a = p.parse_args(argv)

    if not config.check_dependencies():
        return 2

    folder = find_folder(a.folder)
    try:   # one build per folder at a time: concurrent runs clobber and sweep
        atomic.hold(layout.lock(folder, "build"), what=folder)
    except atomic.Locked as exc:
        print(f"\n  Refused: {exc}\n")
        return 2
    cfg_path = settings.path(folder)
    if not os.path.exists(cfg_path):
        print(f"\n  No settings.yaml at {cfg_path}")
        print(f"  Run first:  python3 {os.path.join(HERE,'substrate_init.py')} {folder}\n")
        return 2
    # substrate_init refuses a folder you cannot write with a clear message;
    # this used to be the one step without that check, so a colleague pointing
    # at a folder that already had an index got a PermissionError traceback.
    if not config.is_writable(folder):
        print(f"\n  {folder} is not writable by you.")
        print("  The index is written inside the folder it describes; run this "
              "on a folder you own.\n")
        return 2

    from core import briefs, trajscan

    cfg = settings.config_of(folder)
    if a.deep:
        cfg["deep_scan"] = True
    if a.quiet_period:
        cfg["quiet_period_s"] = a.quiet_period
    registry = layout.registry(folder)
    os.makedirs(registry, exist_ok=True)
    os.makedirs(layout.results(folder), exist_ok=True)

    print(f"\n  results   : {layout.results(folder)}")
    print(f"  root      : {cfg['root']}")
    print(f"  topics    : {', '.join(cfg.get('topics') or ['(none)'])}")
    print(f"  corpora   : {len(config.corpus_roots(cfg))}\n")

    summary = {}

    print("  [1/2] scanning trajectories ...")
    t0 = time.time()
    runs = trajscan.scan(cfg, verbose=True)
    # ALWAYS write it, even with zero rows. Everything downstream starts from
    # runs.parquet; a missing file reads as "the build broke", while an empty
    # one truthfully says "no trajectories here".
    atomic.to_parquet(runs, os.path.join(registry, "runs.parquet"),
                      compression="zstd", index=False)
    _count_runs(summary, runs)
    print(f"        {len(runs)} runs, {summary['traj_gb']} GB "
          f"({time.time()-t0:.0f}s)\n")

    print("  [2/2] writing the page a person reads ...")
    briefs.build_readme(folder, cfg, summary)

    print("\n  Done.\n")
    print(f"      {layout.readme(folder)}   <- start here")
    print(f"      {settings.path(folder)}   <- the file to edit")
    print(f"      {os.path.join(registry, 'runs.parquet')}   (the index)\n")
    return 0


def _count_runs(summary, runs):
    summary["runs"] = len(runs)
    summary["traj_gb"] = round(float(runs.traj_gb.sum()), 1) if len(runs) else 0


if __name__ == "__main__":
    sys.exit(main())
