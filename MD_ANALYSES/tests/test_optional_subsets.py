#!/usr/bin/env python3
"""
test_optional_subsets.py -- a subset a system does not have is not a mistake.

    python3 tests/test_optional_subsets.py

The structure stage resolves receptor, ligand, glycans, epitope, solute and
protein before it opens a trajectory, and an empty selection is an error --
deliberately, because the lab's own TCL scripts carry the warning that the
wrong selection line silently gives zero contacts.

That rule was applied to the GENERIC defaults too, so the first folder with no
sugars in it aborted the whole stage:

    core.selections.EmptySelection: glycans: 'glycan' matched no atoms in this
    topology (translated to '(resname AFU AFUC AGAL ... NAG)')

Ten replicas of a peptide system, every one of them readable, and nothing was
measured. The distinction that fixes it is not "is this subset important" but
WHO WROTE THE SELECTION: `glycans: glycan` is this kit's default and matching
nothing means the system has none; `glycans: segid GA1 to GA23` is somebody's
sentence and matching nothing means they got it wrong.
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np                                          # noqa: E402

from core import selections as S                            # noqa: E402
from core import spec as SP                                 # noqa: E402
from stages.structure import basic                          # noqa: E402


def protein_only_universe():
    """Two chains of alanine. No sugars, no ligand resnames, no epitope."""
    import MDAnalysis as mda
    atoms = [("PA1", 1, "ALA", "CA"), ("PA1", 2, "ALA", "CA"),
             ("PA2", 1, "ALA", "CA"), ("PA2", 2, "ALA", "CA")]
    u = mda.Universe.empty(len(atoms), n_residues=len(atoms),
                           n_segments=2,
                           atom_resindex=np.arange(len(atoms)),
                           residue_segindex=np.array([0, 0, 1, 1]),
                           trajectory=True)
    u.add_TopologyAttr("names", [a[3] for a in atoms])
    u.add_TopologyAttr("resnames", [a[2] for a in atoms])
    u.add_TopologyAttr("resids", [a[1] for a in atoms])
    u.add_TopologyAttr("segids", ["PA1", "PA2"])
    u.atoms.positions = np.zeros((len(atoms), 3), dtype=np.float32)
    return u


class DeclaredSelections(unittest.TestCase):
    def test_a_spec_remembers_which_selections_a_person_wrote(self):
        s = SP.Spec(id="d/a", selections={"receptor": "segid PA1"})
        self.assertEqual(s.declared, {"receptor"})
        # The generic defaults are present but were not declared.
        self.assertEqual(s.selection("glycans"), "glycan")
        self.assertNotIn("glycans", s.declared)

    def test_a_blank_selection_is_not_a_declaration(self):
        s = SP.Spec(id="d/a", selections={"receptor": "segid PA1",
                                          "epitope": "   ", "ligand": None})
        self.assertEqual(s.declared, {"receptor"})

    def test_a_round_trip_does_not_turn_defaults_into_declarations(self):
        """as_dict() hands back EVERY selection, defaults included. Rebuilding
        from that would mark them all as declared and make an absent subset
        fatal again -- one folder later, and only on a rebuild."""
        s = SP.Spec(id="d/a", selections={"receptor": "segid PA1"})
        again = SP.Spec.from_dict(s.as_dict())
        self.assertEqual(again.declared, {"receptor"})


class AbsentSubsets(unittest.TestCase):
    def test_a_system_with_no_glycans_still_gets_measured(self):
        u = protein_only_universe()
        s = SP.Spec(id="d/a", selections={"receptor": "segid PA1",
                                          "ligand": "segid PA2"})
        groups = basic._groups(s, u)
        self.assertIn("receptor", groups)
        self.assertIn("ligand", groups)
        self.assertNotIn("glycans", groups)   # dropped, not invented
        self.assertNotIn("epitope", groups)

    def test_a_selection_a_person_wrote_is_still_fatal_when_it_matches_nothing(self):
        """The whole point of the distinction. Someone typing a segid range
        that does not exist must hear about it, exactly as before."""
        u = protein_only_universe()
        s = SP.Spec(id="d/a", selections={"receptor": "segid PA1",
                                          "glycans": "segid GA1 to GA23"})
        with self.assertRaises(S.EmptySelection) as caught:
            basic._groups(s, u)
        self.assertIn("glycans", str(caught.exception))

    def test_an_empty_receptor_is_always_fatal(self):
        """`receptor` is not optional under any circumstances: it is what
        everything else is fitted and measured against."""
        u = protein_only_universe()
        s = SP.Spec(id="d/a", selections={"receptor": "segid NOPE"})
        with self.assertRaises(S.EmptySelection):
            basic._groups(s, u)


if __name__ == "__main__":
    unittest.main(verbosity=2)
