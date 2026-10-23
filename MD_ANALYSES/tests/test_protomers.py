#!/usr/bin/env python3
"""
test_protomers.py -- one test per way a protomer split can be quietly wrong.

    python3 tests/test_protomers.py

Every segid here is deliberately mis-spelt for the job: the copies are called
Q7/Q3/Q9 and H9/H1/H5, never A/B/C. A rule that reads the first letter, or that
sorts the segids alphabetically, passes on the real 1F4 topology and fails here
-- which is the point, because it also fails on the real JL427 topology, where
the three Fab heavy chains are JLK, JLX and JLY and the one built first is the
one bound to the protomer built LAST.

Synthetic fixtures only; nothing here reads the lab's data.
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from core import protomers as P                             # noqa: E402
from core import topology as T                              # noqa: E402
from test_core import write_psf                             # noqa: E402


def _chain(segid, resname, n):
    return [(segid, r, resname, "CA", 12.0) for r in range(1, n + 1)]


#: A C3 assembly whose names say nothing. Two receptor entities in three copies
#: (the gp120/gp41 pattern), a two-chain ligand in three copies (the Fab
#: pattern), and one stray chain that belongs to no protomer.
#: The FILE ORDER is the build order, and it is the only evidence there is.
ASSEMBLY = (
    _chain("Q7", "ALA", 5) + _chain("Z2", "GLY", 3) +        # protomer 0
    _chain("Q3", "ALA", 5) + _chain("Z8", "GLY", 3) +        # protomer 1
    _chain("Q9", "ALA", 5) + _chain("Z5", "GLY", 3) +        # protomer 2
    _chain("H9", "VAL", 4) + _chain("L4", "TRP", 2) +
    _chain("H1", "VAL", 4) + _chain("L6", "TRP", 2) +
    _chain("H5", "VAL", 4) + _chain("L2", "TRP", 2)
)

RECEPTOR = ["Q7", "Q3", "Q9", "Z2", "Z8", "Z5"]
LIGAND = ["H9", "H1", "H5", "L4", "L6", "L2"]


def entities_of(atoms):
    with tempfile.TemporaryDirectory() as d:
        p = write_psf(os.path.join(d, "s.psf"), atoms)
        return T.entities(T.segments(T.read_psf(p)))


class TestGrouping(unittest.TestCase):
    def test_a_protomer_is_the_i_th_copy_of_every_entity(self):
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND)
        self.assertEqual(s.n, 3)
        self.assertEqual([p["id"] for p in s], ["A", "B", "C"])
        self.assertEqual(s[0]["receptor"], "segid Q7 Z2")
        self.assertEqual(s[1]["receptor"], "segid Q3 Z8")
        self.assertEqual(s[2]["receptor"], "segid Q9 Z5")
        self.assertEqual(s[0]["ligand"], "segid H9 L4")
        self.assertEqual(s[1]["ligand"], "segid H1 L6")
        self.assertEqual(s[2]["ligand"], "segid H5 L2")

    def test_the_answer_is_not_the_alphabetical_one(self):
        """Sorting the segids would put Q3 first and pair it with Z2 -- a
        protomer assembled out of two different ones, with no error."""
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND)
        self.assertNotEqual(s[0]["receptor"], "segid Q3 Z2")

    def test_every_segid_lands_in_exactly_one_protomer(self):
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND)
        seen = [sg for p in s
                for sg in p["receptor_segids"] + p["ligand_segids"]]
        self.assertEqual(sorted(seen), sorted(RECEPTOR + LIGAND))
        self.assertEqual(len(seen), len(set(seen)))

    def test_equivalent_residues_share_an_entity_label(self):
        """Residue 3 of Q7, Q3 and Q9 is the same residue of three copies, and
        the aggregate is indexed by (entity, resid) because of it."""
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND)
        self.assertEqual(s.entity_of["Q7"], s.entity_of["Q3"])
        self.assertEqual(s.entity_of["Q7"], s.entity_of["Q9"])
        self.assertNotEqual(s.entity_of["Q7"], s.entity_of["Z2"])


class TestMonomer(unittest.TestCase):
    def test_one_copy_of_everything_is_one_protomer(self):
        """A monomer is a trimer with n=1: no caller needs a special case."""
        atoms = _chain("Q7", "ALA", 5) + _chain("Z2", "GLY", 3) + \
            _chain("H9", "VAL", 4) + _chain("L4", "TRP", 2)
        s = P.split(entities_of(atoms), ["Q7", "Z2"], ["H9", "L4"])
        self.assertEqual(s.n, 1)
        self.assertEqual(len(s), 1)
        self.assertEqual(s[0]["id"], "A")
        self.assertEqual(sorted(s[0]["receptor_segids"]), ["Q7", "Z2"])
        self.assertEqual(sorted(s[0]["ligand_segids"]), ["H9", "L4"])

    def test_an_apo_assembly_has_no_ligand_rather_than_an_empty_one(self):
        """The apo trimer must produce NO ligand rows -- not empty ones, and
        certainly not zeros, which would plot as a perfectly rigid antibody."""
        rec = _chain("Q7", "ALA", 5) + _chain("Z2", "GLY", 3) + \
            _chain("Q3", "ALA", 5) + _chain("Z8", "GLY", 3) + \
            _chain("Q9", "ALA", 5) + _chain("Z5", "GLY", 3)
        s = P.split(entities_of(rec), RECEPTOR, [])
        self.assertEqual(s.n, 3)
        for p in s:
            self.assertIsNone(p["ligand"])
            self.assertEqual(p["ligand_segids"], [])


class TestExcluded(unittest.TestCase):
    def test_a_single_chain_belongs_to_no_protomer_and_is_reported(self):
        """One CD4 on a trimer is not a third of anything. Splitting it three
        ways would give three protomers of which two are missing a piece."""
        atoms = ASSEMBLY + _chain("X1", "SER", 6)
        s = P.split(entities_of(atoms), RECEPTOR + ["X1"], LIGAND)
        self.assertEqual(s.n, 3)
        for p in s:
            self.assertNotIn("X1", p["receptor_segids"])
            self.assertNotIn("X1", p["ligand_segids"])
        self.assertEqual([e["segids"] for e in s.excluded], [["X1"]])
        self.assertEqual(s.excluded[0]["copies"], 1)
        self.assertEqual(s.excluded[0]["role"], "receptor")
        self.assertTrue(any("X1" in n for n in s.notes))

    def test_the_stray_chains_do_not_decide_the_copy_count(self):
        """Three copies of one subunit plus two different stray chains is
        (3, 1, 1) by entity: on a show of hands the mode is ONE, and the split
        would call this trimer a monomer and report nothing per protomer. The
        mode is taken over atoms, so the subunit that is most of the mass
        decides."""
        atoms = (_chain("Q7", "ALA", 5) + _chain("Q3", "ALA", 5) +
                 _chain("Q9", "ALA", 5) + _chain("X1", "SER", 6) +
                 _chain("X2", "THR", 7))
        s = P.split(entities_of(atoms), ["Q7", "Q3", "Q9", "X1", "X2"], [])
        self.assertEqual(s.n, 3)
        self.assertEqual(sorted(e["segids"][0] for e in s.excluded),
                         ["X1", "X2"])


class TestContactPairing(unittest.TestCase):
    """The JL427 case: the ligand copies are NOT in the receptor's build
    order, and only the coordinates say so."""

    def test_contacts_overrule_the_copy_index_and_say_so(self):
        counts = {"H9": [0, 0, 150], "H1": [280, 0, 0], "H5": [0, 190, 0],
                  "L4": [0, 0, 40], "L6": [95, 0, 0], "L2": [0, 60, 0]}
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND, contacts=counts)
        self.assertEqual(s.rule, "contact")
        self.assertEqual(s[0]["ligand"], "segid H1 L6")
        self.assertEqual(s[1]["ligand"], "segid H5 L2")
        self.assertEqual(s[2]["ligand"], "segid H9 L4")
        # and the receptor is untouched: contacts never redefine the units
        self.assertEqual(s[0]["receptor"], "segid Q7 Z2")
        self.assertTrue(any("disagrees with the copy index" in n
                            for n in s.notes))

    def test_agreeing_contacts_leave_the_answer_alone_and_are_quiet(self):
        counts = {"H9": [150, 0, 0], "H1": [0, 280, 0], "H5": [0, 0, 190],
                  "L4": [40, 0, 0], "L6": [0, 95, 0], "L2": [0, 0, 60]}
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND, contacts=counts)
        self.assertEqual(s.rule, "contact")
        self.assertEqual(s[0]["ligand"], "segid H9 L4")
        self.assertFalse(any("disagrees" in n for n in s.notes))

    def test_two_copies_peaking_on_one_protomer_are_still_placed(self):
        """The 1F4 replica 2-R2 case: a light chain that barely touches the
        receptor peaks on the SAME protomer as its sister copy, so one argmax
        per copy is not even a permutation. The whole assignment reads the same
        numbers and gets the only consistent answer -- and it has to, because a
        protomer that means different chains in different replicas would pool
        unlike things into one mean."""
        counts = {"H9": [150, 0, 0], "H1": [280, 0, 0], "H5": [0, 0, 190],
                  "L4": [40, 0, 0], "L6": [0, 95, 0], "L2": [0, 0, 60]}
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND, contacts=counts)
        self.assertEqual(s.rule, "contact")
        # H1 outbids H9 for protomer A, so H9 takes what is left.
        self.assertEqual(s[0]["ligand"], "segid H1 L4")
        self.assertEqual(s[1]["ligand"], "segid H9 L6")
        self.assertEqual(s[2]["ligand"], "segid H5 L2")

    def test_a_copy_that_touches_nothing_is_placed_by_elimination(self):
        """IFQ makes no heavy-atom contact with any gp120 on 1F4 frame 0 -- it
        sits on its protomer's glycans. One unplaced copy and one free
        protomer is not a guess, and the note says it happened."""
        counts = {"H9": [150, 0, 0], "H1": [0, 280, 0], "H5": [0, 0, 0],
                  "L4": [40, 0, 0], "L6": [0, 95, 0], "L2": [0, 0, 60]}
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND, contacts=counts)
        self.assertEqual(s.rule, "contact")
        self.assertEqual(s[2]["ligand"], "segid H5 L2")
        self.assertTrue(any("H5 touches no receptor protomer" in n
                            for n in s.notes))

    def test_an_undecidable_pairing_is_refused_not_invented(self):
        """Two copies that touch nothing can be swapped without changing a
        single contact. Picking one anyway would build a protomer out of two
        chains that are not on it, and every number after that is fiction."""
        counts = {"H9": [0, 0, 0], "H1": [0, 0, 0], "H5": [0, 0, 190],
                  "L4": [40, 0, 0], "L6": [0, 95, 0], "L2": [0, 0, 60]}
        s = P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND, contacts=counts)
        self.assertEqual(s.rule, "copy_index")
        self.assertEqual(s[0]["ligand"], "segid H9 L4")
        self.assertTrue(any("equally good placements" in n for n in s.notes))

    def test_a_callable_is_given_the_receptor_protomers(self):
        seen = {}

        def contacts(groups):
            seen["groups"] = groups
            return {"H9": [1, 0, 0], "H1": [0, 1, 0], "H5": [0, 0, 1],
                    "L4": [1, 0, 0], "L6": [0, 1, 0], "L2": [0, 0, 1]}

        P.split(entities_of(ASSEMBLY), RECEPTOR, LIGAND, contacts=contacts)
        self.assertEqual(seen["groups"],
                         [["Q7", "Z2"], ["Q3", "Z8"], ["Q9", "Z5"]])


class TestContactCounts(unittest.TestCase):
    def test_contacts_are_counted_between_the_right_groups(self):
        import numpy as np
        X = np.array([[0.0, 0, 0], [1.0, 0, 0],        # protomer 0
                      [50.0, 0, 0], [51.0, 0, 0],      # protomer 1
                      [2.0, 0, 0],                     # ligand near protomer 0
                      [52.0, 0, 0]])                   # ligand near protomer 1
        c = P.contact_counts(X, [[0, 1], [2, 3]], [[4], [5]], cutoff=3.0)
        self.assertGreater(c[0][0], c[0][1])
        self.assertGreater(c[1][1], c[1][0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
