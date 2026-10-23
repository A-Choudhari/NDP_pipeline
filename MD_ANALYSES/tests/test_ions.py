#!/usr/bin/env python3
"""
test_ions.py -- one test per way the ion stage could be quietly wrong.

    python3 tests/test_ions.py

Synthetic fixtures only: a Universe built in memory with a geometry whose
answer is known by construction, so "the first shell holds six waters" is a
fact about the test and not about TIP3P. Nothing here opens the lab's data.

Every failure below is silent in the sense that matters -- it produces a
number, the number is plausible, and nothing raises. A coordination number of
4.9 instead of 5.6 looks like a different force field; a contact map built
from hydrogens looks like a contact map.
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np                                         # noqa: E402

from core import geom, spec as SP                          # noqa: E402
from core import selections as S                           # noqa: E402
from stages.solvent import ions                            # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def universe(atoms, positions, box):
    """A Universe in memory. `atoms` are (segid, resid, resname, name)."""
    import MDAnalysis as mda
    residues, res_of_atom = [], []
    seen = {}
    for segid, resid, resname, _ in atoms:
        key = (segid, resid)
        if key not in seen:
            seen[key] = len(residues)
            residues.append((segid, resid, resname))
        res_of_atom.append(seen[key])
    segids = []
    seg_of_res = []
    for segid, _, _ in residues:
        if segid not in segids:
            segids.append(segid)
        seg_of_res.append(segids.index(segid))
    u = mda.Universe.empty(len(atoms), n_residues=len(residues),
                           n_segments=len(segids),
                           atom_resindex=np.asarray(res_of_atom),
                           residue_segindex=np.asarray(seg_of_res),
                           trajectory=True)
    u.add_TopologyAttr("names", [a[3] for a in atoms])
    u.add_TopologyAttr("resnames", [r[2] for r in residues])
    u.add_TopologyAttr("resids", [r[1] for r in residues])
    u.add_TopologyAttr("segids", segids)
    u.atoms.positions = np.asarray(positions, dtype=np.float32)
    u.dimensions = np.asarray(box, dtype=np.float32)
    return u


#: The six faces of an octahedron, the directions a first shell is built along.
FACES = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0],
                  [0, -1., 0], [0, 0, 1.], [0, 0, -1.]])

#: Eight body diagonals, for a second shell that cannot collide with the first.
CORNERS = np.array([[x, y, z] for x in (1., -1) for y in (1., -1)
                    for z in (1., -1)]) / np.sqrt(3.0)


def hydrated_ion(box=(60., 60., 60.), first=2.5, second=4.5, far=9.0):
    """One SOD with 6 waters in its first shell, 8 in its second, 2 outside.

    Plus a protein residue whose BACKBONE CARBONYL sits 2.0 A from the ion --
    closer than any water -- so a stage that selected water by `name O` would
    count it and come out with seven.
    """
    centre = np.array([30., 30., 30.])
    atoms = [("ION", 1, "SOD", "SOD")]
    xyz = [centre]
    n = 0
    for shell, dirs in ((first, FACES), (second, CORNERS),
                        (far, FACES[:2])):
        for d in dirs:
            n += 1
            atoms += [("WT1", n, "TIP3", "OH2"), ("WT1", n, "TIP3", "H1")]
            xyz += [centre + shell * d, centre + shell * d + [0.4, 0.5, 0.6]]
    atoms += [("A120", 7, "ASP", "N"), ("A120", 7, "ASP", "CA"),
              ("A120", 7, "ASP", "C"), ("A120", 7, "ASP", "O"),
              ("A120", 7, "ASP", "HB1")]
    xyz += [centre + [8.0, 0, 0], centre + [7.0, 0, 0], centre + [6.0, 0, 0],
            centre + [2.0, 0, 0], centre + [1.0, 0, 0]]
    return universe(atoms, xyz, list(box) + [90., 90., 90.])


def test_spec(**selections):
    sel = {"ions": "resname SOD CLA", "solute": "not water and not ions"}
    sel.update(selections)
    return SP.Spec(id="synthetic/1", label="synthetic",
                   components=["protein", "water", "ions"], selections=sel)


# --------------------------------------------------------------------------
class TestWaterSelection(unittest.TestCase):
    """`name O` returns backbone carbonyls, and nothing complains."""

    def test_water_is_selected_by_OH2_not_by_O(self):
        u = hydrated_ion()
        sel = ions._groups(test_spec(), u)
        ag = u.select_atoms(sel["water_oxygen"].mdanalysis)
        self.assertEqual(len(ag), 16)                 # 6 + 8 + 2 waters
        self.assertEqual(set(ag.names), {"OH2"})
        # what the mistake would have picked up: the ASP carbonyl, at 2.0 A
        self.assertIn("ASP", set(u.select_atoms("name O").resnames))

    def test_a_topology_without_OH2_is_refused_by_name(self):
        """A water named O is not a water this stage can count: it would be
        indistinguishable from a carbonyl, so it is an error and not a
        silently smaller shell."""
        u = universe([("WT1", 1, "TIP3", "O"), ("ION", 1, "SOD", "SOD"),
                      ("A120", 1, "ASP", "CB")],
                     [[1., 1, 1], [3., 1, 1], [9., 1, 1]],
                     [40., 40, 40, 90, 90, 90])
        with self.assertRaises(S.EmptySelection) as cm:
            ions._groups(test_spec(), u)
        self.assertIn("OH2", str(cm.exception))


    def test_gromacs_water_named_OW_is_counted(self):
        """OW can only be a water oxygen; the ion stage accepted OH2 alone, so
        a GROMACS/TIP4P/SPC system failed here after passing the water stage."""
        u = hydrated_ion()
        u.atoms.names = ["OW" if n == "OH2" else n for n in u.atoms.names]
        sel = ions._groups(test_spec(), u)
        ag = u.select_atoms(sel["water_oxygen"].mdanalysis)
        self.assertEqual((len(ag), set(ag.names)), (16, {"OW"}))


class TestCoordination(unittest.TestCase):
    def test_the_shells_hold_exactly_what_was_put_in_them(self):
        u = hydrated_ion()
        sel = ions._groups(test_spec(), u)
        ion = u.select_atoms(sel["ions"].mdanalysis)
        wat = u.select_atoms(sel["water_oxygen"].mdanalysis)
        n1, n2 = ions.coordination(ion.positions, wat.positions, u.dimensions)
        self.assertEqual(int(n1[0]), 6)               # 2.5 A < 3.2 A
        self.assertEqual(int(n2[0]), 14)              # + 8 at 4.5 A < 5.6 A
        # and the 2 waters at 9.0 A are in neither

    def test_the_first_shell_is_a_subset_of_the_second(self):
        """One search at the outer radius, thresholded: if the two ever came
        from different searches they could disagree, and a coordination number
        larger in the first shell than the second would be reported without
        comment."""
        u = hydrated_ion()
        sel = ions._groups(test_spec(), u)
        n1, n2 = ions.coordination(
            u.select_atoms(sel["ions"].mdanalysis).positions,
            u.select_atoms(sel["water_oxygen"].mdanalysis).positions,
            u.dimensions)
        self.assertTrue(bool(np.all(n1 <= n2)))

    def test_a_water_across_the_boundary_still_counts(self):
        """The minimum image, which is the whole reason nothing here is
        wrapped: by raw coordinates this water is 18.5 A away, by the
        convention it is 1.5 A and in the first shell."""
        u = universe([("ION", 1, "SOD", "SOD"),
                      ("WT1", 1, "TIP3", "OH2")],
                     [[1., 1, 1], [19.5, 1, 1]], [20., 20, 20, 90, 90, 90])
        ion = u.select_atoms("resname SOD").positions
        wat = u.select_atoms("name OH2").positions
        n1, _ = ions.coordination(ion, wat, u.dimensions)
        self.assertEqual(int(n1[0]), 1)
        # the same geometry with no box is the wrong answer, which is what
        # makes passing `box=` a convention rather than a detail
        n1_nobox, _ = ions.coordination(ion, wat, None)
        self.assertEqual(int(n1_nobox[0]), 0)


class TestPartition(unittest.TestCase):
    """The four classes must partition the ions exactly, or the fractions do
    not sum to one and `bulk` becomes whatever was left over."""

    def setUp(self):
        # glycan at x=0, protein at x=20, in a 60 A box.
        atoms = [("GA1", 1, "AMAN", "C1"), ("A120", 1, "ASP", "CB"),
                 ("A120", 1, "ASP", "HB1")]
        xyz = [[0., 30, 30], [20., 30, 30], [10., 30, 30]]
        for i, x in enumerate((2.0, 11.0, 18.0, 45.0), start=1):
            atoms.append(("ION", i, "SOD", "SOD"))
            xyz.append([x, 30., 30.])
        # Bulk water, far from everything: this stage refuses a system with no
        # water at all, and rightly -- it would have nothing to hydrate.
        for i in (1, 2):
            atoms.append(("WT1", i, "TIP3", "OH2"))
            xyz.append([30. + i, 50., 50.])
        self.u = universe(atoms, xyz, [60., 60, 60, 90, 90, 90])
        self.spec = test_spec(glycans="segid GA1")

    def test_every_ion_gets_exactly_one_class(self):
        sel = ions._groups(self.spec, self.u)
        solute = self.u.select_atoms(sel["solute"].mdanalysis + " and "
                                     + S.HEAVY)
        is_gly = np.isin(solute.indices,
                         self.u.select_atoms(sel["glycans"].mdanalysis + " and "
                                             + S.HEAVY).indices)
        ion = self.u.select_atoms(sel["ions"].mdanalysis)
        code, _, _ = ions.shell_classes(ion.positions, solute.positions,
                                        is_gly, self.u.dimensions)
        got = [ions.PARTITION_CLASSES[c] for c in code]
        self.assertEqual(got, ["glycan_shell", "bulk", "protein_shell",
                               "bulk"])
        self.assertEqual(len(code), len(ion))

    def test_an_ion_touching_both_is_both_not_counted_twice(self):
        u = self.u
        # move the second ion between the two, within 4 A of each
        pos = u.atoms.positions
        pos[4] = [3.5, 30., 30.]                      # 3.5 A from the glycan
        pos[0] = [0., 30, 30]
        pos[1] = [7.0, 30, 30]                        # 3.5 A from the protein
        u.atoms.positions = pos
        sel = ions._groups(self.spec, u)
        solute = u.select_atoms(sel["solute"].mdanalysis + " and " + S.HEAVY)
        is_gly = np.isin(solute.indices,
                         u.select_atoms(sel["glycans"].mdanalysis + " and "
                                        + S.HEAVY).indices)
        ion = u.select_atoms(sel["ions"].mdanalysis)
        code, _, _ = ions.shell_classes(ion.positions, solute.positions,
                                        is_gly, u.dimensions)
        self.assertEqual(ions.PARTITION_CLASSES[code[1]], "both")

    def test_a_hydrogen_is_not_a_contact(self):
        """Heavy atoms by NAME. The ASP hydrogen sits 1 A from an ion whose
        nearest heavy atom is 9 A away; counting it would put that residue in
        a contact map it has no business being in."""
        u = self.u
        pos = u.atoms.positions
        pos[5] = [11.0, 30., 30.]                     # 1 A from HB1 at x=10
        u.atoms.positions = pos
        sel = ions._groups(self.spec, u)
        heavy = u.select_atoms(sel["solute"].mdanalysis + " and " + S.HEAVY)
        self.assertNotIn("HB1", set(heavy.names))
        is_gly = np.isin(heavy.indices,
                         u.select_atoms(sel["glycans"].mdanalysis + " and "
                                        + S.HEAVY).indices)
        ion = u.select_atoms(sel["ions"].mdanalysis)
        code, _, _ = ions.shell_classes(ion.positions, heavy.positions,
                                        is_gly, u.dimensions)
        self.assertEqual(ions.PARTITION_CLASSES[code[1]], "bulk")


class TestEnrichment(unittest.TestCase):
    def test_the_box_volume_cancels(self):
        """Two systems, same local excess, different boxes: raw counts differ
        and the enrichment does not. That is the only reason this number is
        written at all."""
        small = ions.enrichment(10, 100, 50, 10000)
        big = ions.enrichment(20, 200, 100, 20000)
        self.assertAlmostEqual(small, 20.0, 9)
        self.assertAlmostEqual(small, big, 9)

    def test_an_empty_probe_is_none_not_infinity(self):
        self.assertIsNone(ions.enrichment(5, 100, 0, 10000))
        self.assertIsNone(ions.enrichment(0, 0, 10, 10000))

    def test_no_excess_reads_as_one(self):
        self.assertAlmostEqual(ions.enrichment(37, 1000, 370, 10000), 1.0, 9)


class TestDwell(unittest.TestCase):
    def test_episodes_come_from_the_shared_run_length_kernel(self):
        """A contact broken for one frame is two episodes, not one long one.
        geom.run_lengths is what every lifetime in this toolkit uses, so ion
        dwell and water residence cannot disagree about that."""
        mask = np.array([1, 1, 1, 0, 1, 1, 0, 0, 1], dtype=bool)
        got = [int(n) for _, n in geom.run_lengths(mask)]
        self.assertEqual(got, [3, 2, 1])
        st = ions.dwell_stats(got)
        self.assertEqual(st["n"], 3)
        self.assertAlmostEqual(st["median"], 2.0)

    def test_a_distribution_at_the_sampling_floor_is_not_resolved(self):
        """Three of four episodes over before the next frame: the mean is the
        frame spacing, and a table that reports it as a residence time is
        worse than one that reports nothing."""
        st = ions.dwell_stats([1, 1, 1, 5])
        ok, note = ions.dwell_floor(st, 10.0)
        self.assertFalse(ok)
        self.assertIn("single frame", note)
        self.assertIn("10 ns", note)

    def test_a_resolved_distribution_carries_no_note(self):
        st = ions.dwell_stats([4, 5, 6, 1])
        ok, note = ions.dwell_floor(st, 0.2)
        self.assertTrue(ok)
        self.assertIsNone(note)

    def test_the_spacing_verdict_survives_into_the_note(self):
        """check_spacing's message is carried, not raised: the three time
        averages in this stage are unharmed by a coarse stride, and killing
        them to protect the dwell time would lose the cheap answers for the
        expensive one."""
        try:
            SP.check_spacing("ion_dwell", 0.2, 100, ions.MIN_SPACING_NS)
            self.fail("20 ns spacing should not pass a 0.5 ns requirement")
        except SP.SpacingTooCoarse as exc:
            prior = str(exc)
        ok, note = ions.dwell_floor(ions.dwell_stats([4, 5, 6]), 20.0, prior)
        self.assertFalse(ok)
        self.assertIn("0.5 ns or finer", note)


class TestInputs(unittest.TestCase):
    def test_this_stage_refuses_a_stripped_trajectory(self):
        """There are no ions in the stripped file. Reporting zeros from it
        would be a measurement of the strip mask."""
        md = SP.MD(id="1R", topology="merged.psf", trajectories=["merged.dcd"])
        with self.assertRaises(ValueError):
            md.inputs(solvated=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
