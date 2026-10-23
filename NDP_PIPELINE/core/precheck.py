#!/usr/bin/env python3
"""
precheck.py -- one precondition, stated the same way everywhere.

Every step after the index reads what the index already found, so "you have not
built one yet" is the most likely thing to go wrong on a new folder. It gets an
explicit, actionable message rather than a FileNotFoundError, because the fix is
one command and the traceback does not say so.
"""

import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def require(index_dir, runs_parquet):
    """True when the index is usable; otherwise explain and return False."""
    from . import layout
    folder = layout.base(index_dir)
    if not os.path.isdir(index_dir):
        print(f"\n  No index at {index_dir}")
        print("\n  Every step here describes what the index already found,")
        print("  so build that first:\n")
        print(f"      {os.path.join(HERE, 'ndp_pipeline.sh')} {folder}\n")
        return False
    if not os.path.exists(runs_parquet):
        print(f"\n  {index_dir} exists but has no registry/runs.parquet")
        print("  The index did not finish. Re-run:\n")
        print(f"      python3 {os.path.join(HERE, 'index_build.py')} {folder}\n")
        return False
    return True
