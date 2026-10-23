#!/usr/bin/env python3
"""
audit.py -- check what a record CLAIMS against what it actually contains.

Two claims drift away from the data without anything failing, and both end up
in what gets published:

1. **A description quoting a count.** These were generated with the run count of
   the day and then frozen in settings.yaml, where a person may since have
   edited them. When the stage classifier was repaired, LC3B went from 157
   exported MDs to 2 while its description went on saying 157. The text belongs
   to whoever wrote it, so it is never rewritten here -- but a record that
   overstates itself by a factor of 78 must not leave quietly.

2. **The gap between the runs a selector matches and the runs that hold
   production.** Everything else in a coverage report describes the runs that
   were used; without this row nothing anywhere says that 155 of LC3B's 157
   were dropped upstream as equilibration.
"""

import re

#: A number followed by a word that counts simulations. Deliberately narrow:
#: `100 ns` and `CHARMM36m` must not match.
COUNT_CLAIM = re.compile(
    r"\b(\d[\d,]*)\s+(?:indexed\s+)?(?:production\s+)?"
    r"(run directories|runs|replicas|simulations|trajectories|MDs)\b",
    re.I)


def declaration(spec, n_mds, matched=None, selected=None):
    """[(field, status, note)] -- what this dataset's own text gets wrong."""
    rows = []
    for field in ("description", "name"):
        text = spec.get(field)
        if not isinstance(text, str):
            continue
        for claim, noun in COUNT_CLAIM.findall(text):
            stated = int(claim.replace(",", ""))
            if stated == n_mds:
                continue
            rows.append((field, "WARNING",
                         f"says '{stated} {noun}' but this dataset contains "
                         f"{n_mds} MD(s). settings.yaml is where this text "
                         f"lives and it is never rewritten for you; correct it "
                         f"before publishing, or the record overstates itself."))
    if matched is not None and selected is not None:
        dropped = len(matched) - len(selected)
        if dropped > 0:
            rows.append((
                "mds", "note",
                f"the selector matches {len(matched)} indexed run director"
                f"{'ies' if len(matched) != 1 else 'y'}, of which "
                f"{len(selected)} hold production data; the other {dropped} "
                f"hold only minimisation, heating, equilibration or derived "
                f"files and are not part of this dataset."))
    return rows
