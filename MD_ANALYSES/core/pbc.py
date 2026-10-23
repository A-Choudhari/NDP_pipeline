#!/usr/bin/env python3
"""
pbc.py -- the periodic box, handled the one way that is safe here.

THE PROBLEM, in this workspace's own words. NAMD runs these systems with
`wrapAll on`, so every molecule is wrapped into the primary cell
INDEPENDENTLY: a Fab bound to the trimer can be written on the far side of the
box from it. The complex has not moved -- its image has. Anything that measures
a distance, an RMSD or a contact without correcting for that gets a number that
is wrong by a box length, and anything that RENDERS it produces a picture with
the antibody flying off into space, which is what prompted this module.

WHY NOT JUST WRAP. Because `cellOrigin` is not the origin -- for 1F4 it is
(169.7, 169.0, 158.7) -- so MDAnalysis's `AtomGroup.wrap()` and
`pack_into_box()` move atoms relative to the wrong reference and silently tear
the solute apart. The rule inherited from `70_water_shell.py` is absolute:
NOTHING IS EVER WRAPPED. Only explicit minimum-image arithmetic, applied to a
COPY:

    d = x - c;  d -= L * rint(d / L)

MOLECULES MOVE WHOLE. A group is shifted by ONE box vector chosen from its
centre, never atom by atom: imaging atoms individually is exactly what breaks a
molecule across the boundary.
"""

from . import threads                                      # noqa: F401

import numpy as np


class NonOrthorhombic(ValueError):
    """The cell is not 90/90/90 and this module's arithmetic does not apply."""


def box_lengths(dimensions, tol=1e-3):
    """(Lx, Ly, Lz) from an MDAnalysis `ts.dimensions`, orthorhombic only."""
    if dimensions is None or len(dimensions) < 6:
        return None
    L = np.asarray(dimensions[:3], dtype=np.float64)
    angles = np.asarray(dimensions[3:6], dtype=np.float64)
    if np.any(np.abs(angles - 90.0) > tol * 90.0):
        raise NonOrthorhombic(
            f"cell angles {angles} are not 90/90/90; minimum-image arithmetic "
            f"here assumes an orthorhombic box (these runs set "
            f"`useFlexibleCell no`). Use a triclinic-aware distance routine.")
    if np.any(L <= 0):
        return None
    return L


def minimum_image(d, L):
    """Displacements folded into [-L/2, L/2). Operates on a copy."""
    d = np.asarray(d, dtype=np.float64).copy()
    d -= L * np.rint(d / L)
    return d


def shift_into_image(X, reference_point, L, groups=None):
    """Move X, whole, into the periodic image nearest `reference_point`.

    `groups` is an array of labels (one per atom) naming the units that must
    move together -- residues for a solvent molecule, segments for a protein
    chain. Without it the whole of X is treated as one unit, which is right for
    a single chain and wrong for a set of them.

    Returns (imaged coordinates, number of units that had to be moved).
    """
    X = np.asarray(X, dtype=np.float64)
    if L is None or not len(X):
        return X, 0
    ref = np.asarray(reference_point, dtype=np.float64)
    out = X.copy()
    moved = 0
    if groups is None:
        c = X.mean(axis=0)
        n = np.rint((c - ref) / L)
        if np.any(n):
            out -= L * n
            moved = 1
        return out, moved
    starts, counts = groups if isinstance(groups, tuple) else runs_of(groups)
    # One vectorised reduction for every unit, not a boolean mask per unit per
    # frame: the glycans of one trimer are 75 segments, and masking each of
    # them on every one of 2,520 frames turned a two-minute analysis into
    # twenty. `reduceat` sums each contiguous run in a single pass.
    sums = np.add.reduceat(X, starts, axis=0)
    centres = sums / counts[:, None]
    n = np.rint((centres - ref) / L)
    hit = np.any(n != 0, axis=1)
    if np.any(hit):
        out -= np.repeat(n * L, counts, axis=0)
        moved = int(hit.sum())
    return out, moved


def runs_of(labels):
    """(start index, length) of each contiguous run of equal labels.

    Atoms come out of MDAnalysis in index order and a segment's atoms are
    contiguous in these topologies, so a run IS a segment. A label that
    reappears after an interruption starts a new run, which is the safe
    reading: it moves that stretch on its own rather than assuming the two
    stretches belong together.
    """
    labels = np.asarray(labels)
    if not len(labels):
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    change = np.flatnonzero(labels[1:] != labels[:-1]) + 1
    starts = np.r_[0, change]
    counts = np.diff(np.r_[starts, len(labels)])
    return starts, counts


def needs_imaging(X, reference_point, L, factor=0.5):
    """A cheap test that answers "is anything out of place?" in two reductions.

    The guard has to be free when nothing is wrong, because on a trajectory
    that cpptraj has already imaged -- which is most of them -- the answer is
    always no, and paying for a per-segment check on every frame anyway cost
    14 minutes of CPU on a run that otherwise takes two.

    True when the group either spans more than half the box (split across the
    boundary) or sits more than half a box from the reference (wrong image).
    Both are necessary conditions for there being any work to do.
    """
    if L is None or not len(X):
        return False
    lo = X.min(axis=0)
    hi = X.max(axis=0)
    if np.any((hi - lo) > factor * L):
        return True
    centre = 0.5 * (lo + hi)
    return bool(np.any(np.abs(centre - reference_point) > factor * L))


def assembly_center(X):
    """The reference everything else is imaged about: the receptor's centre."""
    return np.asarray(X, dtype=np.float64).mean(axis=0)


def looks_split(X, L, factor=0.5):
    """True when a group spans more than half the box -- i.e. it is split.

    A cheap, honest test: no molecule in these systems is 100 A across in every
    direction, so a span above half the box means the group is living in two
    images at once and any centre of mass computed from it is meaningless.
    """
    if L is None or len(X) < 2:
        return False
    span = X.max(axis=0) - X.min(axis=0)
    return bool(np.any(span > factor * L))
