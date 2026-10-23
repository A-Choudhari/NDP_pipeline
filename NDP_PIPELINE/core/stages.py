#!/usr/bin/env python3
"""
stages.py -- which trajectory files count as the simulation itself.

The classifier lives in `trajscan` and is re-exported here rather than copied:
one vocabulary, used by the index and by everything downstream. What belongs
HERE is the POLICY -- which stages may be published as the simulation. That is
a different question from what a file is.

For the record, since this is the layer that had the bug: `stage_of` used to
match filename prefixes, so `protein_sol_nacl_eq.mdcrd` read as production via
*pr*otein and 468 equilibration files were listed as production data across the
catalogue. See trajscan.py for the full account.
"""

import os

from .trajscan import classify_dir, frame_info, stage_of   # noqa: F401

#: Stages that may be deposited as the simulation itself.
#:
#: Deliberately just `production`. `merged` is built FROM the production files,
#: so listing it would present the same sampling twice; `stripped` and `derived`
#: are re-encodings that are not a bit-for-bit stand-in for the raw run; and
#: `unresolved` means nobody has established what the file is yet, which is a
#: reason to look rather than a reason to publish.
DEPOSITABLE = {"production"}


def is_depositable(filename):
    """True when this file may be listed in an mwf inputs.yaml."""
    return stage_of(filename) in DEPOSITABLE


def classify(filenames):
    """{stage: [filename, ...]} for one run directory, for reporting."""
    out = {}
    for name in filenames:
        out.setdefault(stage_of(name), []).append(name)
    return out


class Classifier:
    """stage_of for export, with what the index knows and a name does not.

    Called exactly like stage_of(name, parent), so everything that takes a
    classifier takes this. The difference is that it answers for the whole
    directory the file sits in (trajscan.classify_dir): the author's `stages:`
    rules from substrate.yaml first, then the name, then the .dcd default --
    which depends on the siblings, and is why a one-file function cannot give it.
    Without this the export would disagree with the index about which files a
    run's production is, and the index is the authority.

    `basis(name, parent)` says how a stage was established: author, name,
    atoms or default. Only `default` is an assumption, and the glossary flags it.
    """

    def __init__(self, stage_rules=(), root=None, traj_ext=None):
        from . import config
        self.rules = list(stage_rules)
        self.root = os.path.abspath(root) if root else None
        self.traj_ext = tuple(traj_ext or config.TRAJ_EXT)
        self._skip = config.skip_symlink
        self._dirs = {}

    def _dir(self, parent):
        key = os.path.abspath(str(parent))
        if key not in self._dirs:
            try:
                # The same files the index counted: trajectories that are not
                # links it was told to skip. A different listing would let the
                # default apply here and not there.
                names = [n for n in os.listdir(key) if n.endswith(self.traj_ext)
                         and not self._skip(os.path.join(key, n))]
            except OSError:
                names = []
            rel = os.path.relpath(key, self.root) if self.root else None
            if rel in (".", None):
                rel = ""
            # Atom counts from the headers, as the index read them: a file
            # holding fewer atoms than its neighbours is a subset, not the run.
            natoms = {n: frame_info(os.path.join(key, n))[1] for n in names}
            self._dirs[key] = classify_dir(key, names, self.rules, rel, natoms)
        return self._dirs[key]

    def __call__(self, name, parent=None):
        if parent is None:
            return stage_of(name)
        return self._dir(parent).get(name, (stage_of(name, parent=parent), "name"))[0]

    def basis(self, name, parent):
        return self._dir(parent).get(name, (None, "name"))[1]


def classifier(cfg):
    """This folder's classifier, from its settings (`folder:` section)."""
    from . import config
    return Classifier(config.compile_stage_rules(cfg), cfg.get("root"))
