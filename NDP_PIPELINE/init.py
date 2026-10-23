#!/usr/bin/env python3
"""
init.py -- set up a SUBSTRATE in your own folder.

    python3 init.py                 # interactive, in the current dir
    python3 init.py /path/to/mine   # somewhere else
    python3 init.py --yes           # accept all proposals, no prompts

Creates <root>/SUBSTRATE/substrate.yaml plus a curated.yaml scaffold. Writes
NOTHING outside <root>/SUBSTRATE/, and reads nothing outside <root>.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from core import config, layout, settings, topics                                # noqa: E402


# Files that mark a directory as a software checkout rather than science data.
SOFTWARE_MARKERS = {
    "setup.py", "pyproject.toml", "CMakeLists.txt", "configure", "configure.ac",
    "Makefile.am", "package.json", "Cargo.toml", "setup.cfg", "meson.build",
    "environment.yml", "requirements.txt", "LICENSE", "CONTRIBUTING.md",
}
ARCHIVE_EXT = (".tar", ".gz", ".tgz", ".zip", ".bz2", ".xz", ".zst", ".whl")
# Force-field / parameter distributions. CHARMM36m's toppar files proposed
# Glycans-Mucins, Lipids-Aerosols and Nucleic-acids for a purely HIV project.
FF_EXT = (".rtf", ".prm", ".str", ".inp", ".par", ".top", ".itp", ".frcmod", ".lib")
# Names that are always toolchain, never a science corpus.
ALWAYS_SKIP = {"SUBSTRATE", "MDDB_substrate", "MD_RESULTS", "NDP_PIPELINE",
               "MDDB_PIPELINE", "NDP", "DOWNLOAD", "Downloads",
               "software", "Software", "bin", "tools", "toppar"}


#: Below this many loose files, a "most of them look like X" test is not
#: a measurement. Three is the smallest count where a majority means
#: anything at all.
MIN_FILES_TO_JUDGE = 3


def looks_like_software(path):
    """True when a directory is a code checkout, force field or download cache.

    Indexing these adds thousands of vendor files and, worse, drags in false
    topics -- alphafold3's own test data proposed 'SARS-CoV-2' for a lab member
    who has never worked on it.
    """
    try:
        entries = os.listdir(path)
    except OSError:
        return True
    names = set(entries)
    if names & SOFTWARE_MARKERS:
        return True
    files = [e for e in entries if os.path.isfile(os.path.join(path, e))]
    # A MAJORITY OF ONE IS NOT EVIDENCE. These tests ask whether most of the
    # loose files look like vendor material, which says nothing at all when
    # there are one or two of them. `B40T77/newversion/` holds a single zipped
    # postscript file beside the directory with all twelve of its production
    # runs in it -- 338 GB, 5.2 µs -- and "1 of 1 files is an archive" was
    # enough to classify the whole subtree as a software checkout. The build
    # then reported "No production runs. This folder may hold analysis only",
    # which is the worst kind of wrong: confident, quiet, and about absence.
    if len(files) >= MIN_FILES_TO_JUDGE:
        if sum(1 for f in files if f.endswith(ARCHIVE_EXT)) > len(files) / 2:
            return True
        if sum(1 for f in files
               if f.endswith(FF_EXT) or f.startswith("toppar")) > len(files) / 2:
            return True
    if len(entries) >= MIN_FILES_TO_JUDGE and \
            sum(1 for e in entries if e.startswith("toppar")) > len(entries) / 2:
        return True
    return False


def pick_corpora(root, prune, limit=40):
    """Top-level directories worth scanning, largest fan-out first."""
    out, skipped = [], []
    for d in sorted(os.listdir(root)):
        full = os.path.join(root, d)
        if not os.path.isdir(full) or d.startswith(".") or d in set(prune):
            continue
        if d in ALWAYS_SKIP:
            skipped.append((d, "toolchain/download directory"))
            continue
        if looks_like_software(full):
            skipped.append((d, "looks like a software checkout"))
            continue
        try:
            n = len(os.listdir(full))
        except OSError:
            continue
        out.append((d, n))
    out.sort(key=lambda t: -t[1])
    return [d for d, _ in out[:limit]], skipped


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("root", nargs="?", default=os.getcwd())
    p.add_argument("--yes", action="store_true", help="accept proposals without prompting")
    p.add_argument("--force", action="store_true", help="overwrite an existing substrate.yaml")
    a = p.parse_args(argv)

    if not config.check_dependencies():
        return 2

    root = os.path.abspath(os.path.expanduser(a.root))
    if not os.path.isdir(root):
        print(f"  not a directory: {root}")
        return 2
    if not config.is_writable(root):
        print(f"\n  {root} is not writable by you.")
        print("  Run this inside your OWN folder -- a substrate must be built by its owner.")
        return 2

    subst = layout.index_dir(root)
    cfg_path = settings.path(root)
    if os.path.exists(cfg_path) and not a.force:
        print(f"\n  {cfg_path} already exists. Use --force to overwrite,")
        print("  or just run: python3 index_build.py")
        return 0

    print(f"\n  Setting up {layout.RESULTS}/ for: {root}\n")

    cfg = config.default_config(root)
    corpora, skipped = pick_corpora(root, cfg["prune"])
    print(f"  Found {len(corpora)} candidate top-level directories to index:")
    print("    " + ", ".join(corpora))
    if skipped:
        print(f"\n  Excluded {len(skipped)} (add them under `folder: corpora:` in "
              f"settings.yaml if you disagree):")
        for d, why in skipped:
            print(f"    {d:<24} {why}")

    print("\n  Scanning names for research topics (this reads names only, not data)...")
    names = []
    for c in corpora:
        names.append(c)
        names += topics.collect_names(os.path.join(root, c), max_depth=2)
    proposals = topics.propose(names)
    print()
    print(topics.format_proposal(proposals))

    chosen = sorted(proposals, key=lambda t: -len(proposals[t]))
    if not a.yes and sys.stdin.isatty():
        print("\n  These are PROPOSALS from directory and file names. Nothing is assigned")
        print("  until you confirm. Enter a comma-separated list to accept or correct,")
        print("  or press Enter to accept all of the above.")
        reply = input(f"\n  topics [{', '.join(chosen) or 'none'}]: ").strip()
        if reply:
            chosen = [t.strip() for t in reply.split(",") if t.strip()]
    else:
        print("\n  --yes: accepting the proposals above.")

    cfg["topics"] = chosen
    cfg["corpora"] = [{"name": c, "path": c, "engine": "auto"} for c in corpora]

    os.makedirs(subst, exist_ok=True)
    settings.create(root, folder_section=config.folder_section(cfg))

    print(f"\n  wrote {cfg_path}")
    print(f"\n  Topics recorded: {', '.join(chosen) or '(none)'}")
    print(f"  Corpora to index: {len(corpora)}")
    print("\n  Next:\n")
    print(f"      python3 {os.path.join(HERE, 'index_build.py')} {root}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
