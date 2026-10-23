#!/usr/bin/env python3
"""
test_contacts.py -- one test per way the interface stage could be quietly wrong.

    python3 tests/test_contacts.py

Synthetic fixtures only: two small chains at a separation chosen so the contact
set is exactly computable by hand, a periodic case where the true partner is
across the box, and a twelve-chain toy assembly written out as a PSF so the seam
with `core.protomers` is exercised for real. Nothing here opens the lab's data,
so "this pair is in contact and that one is not" is a fact about the test and
not about 1F4.

EVERY FAILURE BELOW IS SILENT IN THE SENSE THAT MATTERS. A contact map built
from hydrogens looks like a contact map. A pair counted once per close atom
looks like a contact map in which big residues matter more. A minimum image
that is not applied looks like an interface that fell apart. An aggregate that
averages only over the protomers where something happened looks like a
universal interaction. A ligand chain placed on the wrong protomer looks like a
contact map with an unexplained asymmetry. None of them raise, all of them make
a figure, and the figure is wrong.
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np                                          # noqa: E402

from core import geom, spec as SP                           # noqa: E402
from core import selections as S                            # noqa: E402
from stages.contacts import interface                       # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def universe(atoms, positions, box):
    """A Universe in memory. `atoms` are (segid, resid, resname, name)."""
    import MDAnalysis as mda
    residues, res_of_atom, seen = [], [], {}
    for segid, resid, resname, _ in atoms:
        key = (segid, resid)
        if key not in seen:
            seen[key] = len(residues)
            residues.append((segid, resid, resname))
        res_of_atom.append(seen[key])
    segids, seg_of_res = [], []
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


def write_psf(path, atoms):
    """The smallest CHARMM PSF core.topology will read: just the !NATOM block.

    The stage hands `core.protomers` a topology PATH, so the seam can only be
    tested with a real file. Atom ORDER is what protomers takes the copy index
    from, so it is the same order the Universe above is built in.
    """
    with open(path, "w") as fh:
        fh.write("PSF\n\n       1 !NTITLE\n* synthetic\n\n")
        fh.write(f"{len(atoms):8d} !NATOM\n")
        for i, (segid, resid, resname, name) in enumerate(atoms, start=1):
            fh.write(f"{i:8d} {segid:<8s} {resid:<8d} {resname:<8s} "
                     f"{name:<8s} CT1  0.000000  12.011000           0\n")
        fh.write("\n       0 !NBOND: bonds\n")
    return path


def two_chains(gap=4.0, box=(60., 60., 60.)):
    """One receptor chain and one ligand chain, contacts known by construction.

    Three receptor residues 20 A apart along y, three ligand residues facing
    them across a gap in x. With gap = 4.0 and a 4.5 A cutoff the contact set is
    EXACTLY {(R 1, L 101), (R 2, L 102)}: the third ligand residue is placed
    10 A out, and no residue can reach its neighbour's partner because 20 A
    along y dwarfs the cutoff. Each chain's hydrogens point at the other chain
    and are closer than any heavy atom, so a stage that selected by mass would
    report contacts the geometry does not contain.
    """
    atoms, xyz = [], []
    for i, resid in enumerate((1, 2, 3)):
        y = 20.0 * i
        atoms += [("RA", resid, "ALA", "CA"), ("RA", resid, "ALA", "CB"),
                  ("RA", resid, "ALA", "HB1")]
        xyz += [[0.0, y, 0.0], [0.0, y, 1.0], [1.0, y, 0.0]]
    for i, resid in enumerate((101, 102, 103)):
        y = 20.0 * i
        x = gap if resid != 103 else 10.0
        atoms += [("LX", resid, "GLY", "CA"), ("LX", resid, "GLY", "N"),
                  ("LX", resid, "GLY", "HA1")]
        xyz += [[x, y, 0.0], [x, y, 1.0], [x - 1.0, y, 0.0]]
    return universe(atoms, xyz, list(box) + [90., 90., 90.])


def pairs_of(u, receptor="segid RA", ligand="segid LX", box=True,
             cutoff=interface.CONTACT_A):
    """The stage's residue-pair answer for one frame, as readable labels."""
    rec = u.select_atoms(S.translate(receptor) + " and " + S.HEAVY)
    lig = u.select_atoms(S.translate(ligand) + " and " + S.HEAVY)
    r_of, r_seg, r_resid, _ = interface._residue_table(rec)
    l_of, l_seg, l_resid, _ = interface._residue_table(lig)
    r_idx, l_idx = interface.residue_pairs(
        rec.positions, lig.positions, r_of, l_of, len(l_seg),
        u.dimensions if box else None, cutoff=cutoff)
    return sorted((f"{r_seg[a]}:{r_resid[a]}", f"{l_seg[b]}:{l_resid[b]}")
                  for a, b in zip(r_idx.tolist(), l_idx.tolist()))


def pairs_of_including_hydrogens(u):
    """What the map would hold if heavy atoms were not filtered -- the mistake
    the test below exists to catch, spelled out so it is visible."""
    rec = u.select_atoms("segid RA")
    lig = u.select_atoms("segid LX")
    r_of, _r_seg, _r_resid, _ = interface._residue_table(rec)
    l_of, l_seg, _l_resid, _ = interface._residue_table(lig)
    a, b = interface.residue_pairs(rec.positions, lig.positions, r_of, l_of,
                                   len(l_seg), u.dimensions)
    return list(zip(a.tolist(), b.tolist()))


def test_spec(**selections):
    sel = {"receptor": "segid RA", "ligand": "segid LX"}
    sel.update(selections)
    return SP.Spec(id="synthetic/1", label="synthetic",
                   components=["protein", "antibody"], selections=sel)


# --------------------------------------------------------------------------
class TestTheContactSet(unittest.TestCase):
    """Two chains at a known separation: arithmetic, not biology."""

    def test_exactly_the_pairs_that_are_inside_the_cutoff(self):
        got = pairs_of(two_chains(gap=4.0))
        self.assertEqual(got, [("RA:1", "LX:101"), ("RA:2", "LX:102")])

    def test_one_angstrom_further_out_and_there_are_no_contacts(self):
        """4.6 A is outside a 4.5 A cutoff. If this ever passes, the cutoff is
        not being applied where the docstring says it is."""
        self.assertEqual(pairs_of(two_chains(gap=4.6)), [])

    def test_a_hydrogen_is_not_a_contact(self):
        """Heavy atoms by NAME. At a 6 A gap the heavy atoms of the two chains
        are 6.0 A apart and there is no interface, while their hydrogens -- each
        pointing at the other chain -- are 4.0 A apart and inside the cutoff. A
        selection by mass, or `not type H` on a topology whose masses are
        guessed (MDAnalysis reads the atom name CA as carbon), puts two contacts
        in the map that the geometry does not contain."""
        u = two_chains(gap=6.0)
        self.assertEqual(pairs_of(u), [])
        self.assertEqual(len(pairs_of_including_hydrogens(u)), 2,
                         "the fixture must have hydrogens inside the cutoff "
                         "for this test to mean anything")

    def test_a_pair_counts_once_however_many_atoms_are_close(self):
        """The observable is "are these two residues touching", not "how many
        of their atoms are". Counting atom pairs would make the map a map of
        residue size: a TRP against a TYR would outweigh a salt bridge."""
        atoms, xyz = [], []
        for k in range(6):                    # one fat receptor residue
            atoms.append(("RA", 1, "TRP", f"C{k}"))
            xyz.append([0.0, 0.5 * k, 0.0])
        for k in range(5):                    # one fat ligand residue
            atoms.append(("LX", 101, "TYR", f"C{k}"))
            xyz.append([3.0, 0.5 * k, 0.0])
        u = universe(atoms, xyz, [60., 60, 60, 90, 90, 90])
        from MDAnalysis.lib.distances import capped_distance
        n_atom_pairs = len(capped_distance(
            u.select_atoms("segid RA").positions.astype(np.float32),
            u.select_atoms("segid LX").positions.astype(np.float32),
            interface.CONTACT_A, box=u.dimensions, return_distances=False))
        self.assertGreater(n_atom_pairs, 10)
        self.assertEqual(pairs_of(u), [("RA:1", "LX:101")])


class TestAcrossTheBoundary(unittest.TestCase):
    """The minimum image, which is the whole reason nothing here is wrapped."""

    def test_the_true_partner_is_across_the_box(self):
        """By raw coordinates these two residues are 18.5 A apart and there is
        no interface. By the minimum image they are 1.5 A apart and in contact.
        `wrapAll on` puts a bound antibody on the far side of the cell as a
        matter of routine, so this is the normal case, not the exotic one."""
        u = universe([("RA", 1, "ALA", "CA"), ("LX", 101, "GLY", "CA")],
                     [[1.0, 1.0, 1.0], [19.5, 1.0, 1.0]],
                     [20., 20, 20, 90, 90, 90])
        self.assertEqual(pairs_of(u), [("RA:1", "LX:101")])
        # The same geometry with no box is the wrong answer, which is what makes
        # passing `box=` a convention here rather than an implementation detail.
        self.assertEqual(pairs_of(u, box=False), [])

    def test_nothing_was_moved_to_get_that_answer(self):
        """capped_distance applies the convention internally. If any coordinate
        changed, something wrapped -- and `cellOrigin` for these systems is not
        the origin, so wrapping tears the solute apart."""
        u = universe([("RA", 1, "ALA", "CA"), ("LX", 101, "GLY", "CA")],
                     [[1.0, 1.0, 1.0], [19.5, 1.0, 1.0]],
                     [20., 20, 20, 90, 90, 90])
        before = u.atoms.positions.copy()
        pairs_of(u)
        self.assertTrue(np.array_equal(before, u.atoms.positions))

    def test_a_partner_across_the_box_is_not_invented(self):
        """The other half of the same guard: 10 A apart by minimum image is 10 A
        apart, and a stage that folded distances wrongly could produce contacts
        out of nothing at all."""
        u = universe([("RA", 1, "ALA", "CA"), ("LX", 101, "GLY", "CA")],
                     [[1.0, 1.0, 1.0], [11.0, 1.0, 1.0]],
                     [20., 20, 20, 90, 90, 90])
        self.assertEqual(pairs_of(u), [])

    def test_the_split_evidence_is_a_minimum_image_count_too(self):
        """The reason this stage does not call core.protomers.contact_counts,
        which uses plain Euclidean distance and says it needs coordinates
        already brought into one image: a chain that `wrapAll on` wrote a box
        length from the receptor it is bound to would be counted as touching
        nothing, and the split would then place it by elimination."""
        u = universe([("RA", 1, "ALA", "CA"), ("LX", 101, "GLY", "CA")],
                     [[1.0, 1.0, 1.0], [19.5, 1.0, 1.0]],
                     [20., 20, 20, 90, 90, 90])
        rec = u.select_atoms("segid RA")
        lig = u.select_atoms("segid LX")
        with_box = interface.neighbour_counts(lig.positions, [0], 1,
                                              rec.positions, [0], 1,
                                              u.dimensions, cutoff=3.0)
        without = interface.neighbour_counts(lig.positions, [0], 1,
                                             rec.positions, [0], 1, None,
                                             cutoff=3.0)
        self.assertEqual(int(with_box[0, 0]), 1)
        self.assertEqual(int(without[0, 0]), 0)


class TestSelections(unittest.TestCase):
    def test_a_system_with_no_ligand_is_refused(self):
        """An apo system has no interface. Returning an empty table would make
        "nothing bound" and "you forgot to say what the ligand is" the same
        result."""
        u = two_chains()
        spec = SP.Spec(id="apo/1", label="apo",
                       selections={"receptor": "segid RA"})
        with self.assertRaises(S.EmptySelection) as cm:
            interface._groups(spec, u)
        self.assertIn("ligand", str(cm.exception))

    def test_a_selection_that_matches_nothing_is_an_error(self):
        u = two_chains()
        with self.assertRaises(S.EmptySelection):
            interface._groups(test_spec(ligand="segid NOPE"), u)

    def test_the_spec_text_goes_through_the_translator(self):
        """This stage never calls select_atoms on the spec's raw text: it goes
        through core/selections.py, which is what makes `noh`, `segid A to C`
        and the chainID warning behave the same here as everywhere else. The
        resolved selection reports its atom count, so a person can see that the
        ligand really is 5,376 atoms and not 0."""
        resolved = interface._groups(test_spec(ligand="segid LX and noh"),
                                     two_chains())
        self.assertIn("not name H*", resolved["ligand"].mdanalysis)
        self.assertEqual(resolved["ligand"].n_atoms, 6)
        self.assertEqual(resolved["receptor"].n_atoms, 9)


# --------------------------------------------------------------------------
# the seam with core.protomers
# --------------------------------------------------------------------------
#: Three receptor protomers 100 A apart, and three two-chain ligands placed so
#: that the COPY INDEX IS WRONG: ligand copy 1 (H1+L1) sits on receptor copy 3,
#: copy 2 on copy 1, copy 3 on copy 2. This is the JL427 case in miniature, and
#: it is the only thing that tells the two rules apart -- a fixture where the
#: copy index happens to be right would pass whether or not the contact
#: evidence ever reached `split`.
LIGAND_HOME = {"H1": 2, "L1": 2, "H2": 0, "L2": 0, "H3": 1, "L3": 1}


def mini_assembly(tmpdir, box=(400., 60., 60.)):
    """A twelve-chain toy trimer: a PSF on disk and the Universe beside it.

    Two receptor entities (R of 5 residues, G of 3) and two ligand entities
    (heavy of 4, light of 2), three copies each, written in copy order so the
    topology carries the build order `core.protomers` takes its copy index from.
    """
    seqs = {"R": ["GLU"] * 5, "G": ["ASP"] * 3,
            "H": ["ALA", "GLY", "SER", "THR"], "L": ["VAL", "LEU"]}
    atoms, xyz = [], []

    def chain(segid, kind, origin):
        for i, rn in enumerate(seqs[kind]):
            for k, nm in enumerate(("CA", "CB")):
                atoms.append((segid, i + 1, rn, nm))
                xyz.append([origin[0], origin[1] + 3.0 * i, origin[2] + k])

    for family, dx in (("R", -6.0), ("G", -6.0)):
        for c in range(3):
            chain(f"{family}{'ABC'[c]}", family,
                  (100.0 * c + dx, 0.0 if family == "R" else 15.0, 0.0))
    for family, dx in (("H", 0.0), ("L", 3.0)):
        for c in range(3):
            home = LIGAND_HOME[f"{family}{c + 1}"]
            chain(f"{family}{c + 1}", family, (100.0 * home + dx, 0.0, 0.0))
    psf = write_psf(os.path.join(tmpdir, "mini.psf"), atoms)
    return psf, universe(atoms, xyz, list(box) + [90., 90., 90.])


class TestTheProtomerSplit(unittest.TestCase):
    """What this stage owes core.protomers, and what it must not decide."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.psf, self.u = mini_assembly(self.tmp.name)
        self.rec = self.u.select_atoms("segid RA RB RC GA GB GC")
        self.lig = self.u.select_atoms("segid H1 H2 H3 L1 L2 L3")

    def tearDown(self):
        self.tmp.cleanup()

    def _split(self, cutoff=interface.PAIRING_A):
        return interface.split_assembly(
            self.psf,
            np.asarray([str(s) for s in self.rec.segids]), self.rec.positions,
            np.asarray([str(s) for s in self.lig.segids]), self.lig.positions,
            self.u.dimensions, cutoff=cutoff)

    def test_the_receptor_protomers_come_from_the_copy_index(self):
        """Never from proximity, and never from the segid. On the real Env
        trimers gp120 of one protomer packs against gp41 of the NEXT, so a
        proximity rule builds a unit that does not exist -- which is exactly why
        this stage asks core.protomers instead of deciding for itself."""
        split, _ = self._split()
        self.assertEqual(split.n, 3)
        self.assertEqual([p["receptor_segids"] for p in split],
                         [["RA", "GA"], ["RB", "GB"], ["RC", "GC"]])

    def test_the_contact_evidence_reaches_split_and_overrules_the_index(self):
        """The whole point of supplying `contacts`. If the evidence never
        arrived, or arrived in the wrong shape, `split` would fall back to the
        copy index without failing -- and this fixture's ligands would each be
        bolted onto a protomer 100 A away."""
        split, evidence = self._split()
        self.assertEqual(split.rule, "contact")
        self.assertEqual({p["id"]: p["ligand_segids"] for p in split},
                         {"A": ["H2", "L2"], "B": ["H3", "L3"],
                          "C": ["H1", "L1"]})
        self.assertTrue(any("disagrees with the copy index" in n
                            for n in split.notes),
                        "core.protomers must be given the chance to say that "
                        "the two rules differ; a silent agreement is not the "
                        "same as agreement")
        # and the evidence really is per protomer, in protomer order
        self.assertEqual(len(evidence["H1"]), 3)
        self.assertEqual(int(np.argmax(evidence["H1"])), LIGAND_HOME["H1"])

    def test_the_entity_labels_are_what_makes_resid_equivalence_work(self):
        """Residue 384 of A120 and of B120 are the same residue only because
        the two segments hold the same SEQUENCE. `split.entity_of` is where that
        comes from, and the aggregate keys are built on it rather than on any
        property of the name."""
        split, _ = self._split()
        ent = split.entity_of
        self.assertEqual(ent["RA"], ent["RB"])
        self.assertEqual(ent["RB"], ent["RC"])
        self.assertNotEqual(ent["RA"], ent["GA"])
        self.assertEqual(ent["H1"], ent["H3"])
        self.assertNotEqual(ent["H1"], ent["L1"])


class TestProtomerOffsets(unittest.TestCase):
    """The aggregation key: a difference round the ring, not a boolean."""

    def setUp(self):
        self.units = ["A", "B", "C"]
        self.index_of = {"A": 0, "B": 1, "C": 2}

    def test_the_ligands_own_protomer_is_offset_zero(self):
        off = interface.protomer_offsets(["A", "B", "C"], self.index_of,
                                         self.units, 3)
        self.assertEqual([int(off[i, i]) for i in range(3)], [0, 0, 0])

    def test_the_same_neighbour_relation_gets_the_same_offset(self):
        """Under three-fold symmetry, A's ligand touching B and B's ligand
        touching C are THE SAME measurement. Pooling them is what turns a
        neighbour contact into a number with a spread instead of three
        unrelated rows; a cognate/cross flag would pool the A-to-C contact in
        with them, which is different geometry."""
        off = interface.protomer_offsets(["A", "B", "C"], self.index_of,
                                         self.units, 3)
        self.assertEqual(int(off[0, 1]), int(off[1, 2]))     # A->B == B->C
        self.assertEqual(int(off[0, 1]), 1)
        self.assertEqual(int(off[0, 2]), 2)                  # A->C is not that
        self.assertNotEqual(int(off[0, 1]), int(off[0, 2]))

    def test_a_chain_in_no_protomer_has_no_place_in_the_ring(self):
        """One CD4 on a trimer is not a third of anything. core.protomers
        excludes it, and a -1 here keeps it out of every aggregate rather than
        being pooled as if it were a repeating unit."""
        off = interface.protomer_offsets(["A", None, "C"], self.index_of,
                                         self.units, 3)
        self.assertTrue(bool(np.all(off[:, 1] == -1)))
        self.assertTrue(bool(np.all(off[:, 0] >= 0)))

    def test_a_monomer_is_a_trimer_with_n_equal_one(self):
        off = interface.protomer_offsets(["A", "A"], {"A": 0}, ["A"], 1)
        self.assertEqual(off.shape, (1, 2))
        self.assertTrue(bool(np.all(off == 0)))


# --------------------------------------------------------------------------
# pooling
# --------------------------------------------------------------------------
class TestPooling(unittest.TestCase):
    def test_sd_is_over_the_units_and_sem_is_over_the_replicas(self):
        """Nine values, three replicas of three protomers. The sd describes what
        the nine did; the s.e.m. describes what an independent repeat would do,
        and it is taken over the three replica means ONLY because protomers in
        one box are pseudoreplicates. The two numbers differ by a factor of
        nearly two here, and quoting the wrong one is the whole point of
        core/io.py's SCOPES note."""
        values = [1., 2, 3, 4, 5, 6, 7, 8, 9]
        replicas = ["r1"] * 3 + ["r2"] * 3 + ["r3"] * 3
        st = interface.pool(values, replicas)
        self.assertAlmostEqual(st["value"], 5.0, 9)
        self.assertEqual(st["n"], 9)
        self.assertEqual(st["n_replicas"], 3)
        self.assertAlmostEqual(st["sd"], float(np.sqrt(7.5)), 9)
        # replica means 2, 5, 8 -> sd 3.0 -> sem 3/sqrt(3)
        self.assertAlmostEqual(st["sem"], 3.0 / np.sqrt(3.0), 9)
        # and it is NOT sd/sqrt(9), which would claim three times the precision
        self.assertGreater(st["sem"], st["sd"] / 3.0)

    def test_one_replica_has_no_standard_error(self):
        st = interface.pool([1., 2, 3], ["r1"] * 3)
        self.assertIsNone(st["sem"])
        self.assertEqual(st["n_replicas"], 1)
        self.assertIsNotNone(st["sd"])

    def test_a_unit_that_never_saw_the_pair_contributes_its_zero(self):
        """THE aggregation bug. A contact present in one protomer of nine and
        absent from the rest is an occupancy of about 0.1, not of 1.0. Averaging
        over only the units that happened to show it is how a figure comes to
        claim a universal interaction from a single-protomer excursion."""
        agg = interface.Pooled()
        for replica in ("r1", "r2", "r3"):
            for protomer in ("A", "B", "C"):
                agg.register("pair", replica, protomer)
        agg.add("pair", ("E", 384, "H", 100, 0), "r1", "A", 0.9,
                observable="contact_occupancy")
        (key, meta, values, replicas), = list(agg.items("pair"))
        self.assertEqual(len(values), 9)
        self.assertEqual(sorted(values), [0.0] * 8 + [0.9])
        self.assertAlmostEqual(interface.pool(values, replicas)["value"],
                               0.1, 9)

    def test_groups_do_not_borrow_each_others_units(self):
        """A receptor residue is pooled over the copies of ITS entity, a pair
        over the ligand copies. Filling zeros from the wrong list would divide
        by nine something that only ever had three chances to happen."""
        agg = interface.Pooled()
        agg.add("residue:receptor:GP120:residue_occupancy", ("GP120", 384),
                "r1", "A", 0.6, observable="residue_occupancy")
        for pid in ("A", "B", "C"):
            agg.register("residue:receptor:GP120:residue_occupancy", "r1", pid)
        agg.add("residue:receptor:CD4:residue_occupancy", ("CD4", 20),
                "r1", "A", 0.3, observable="residue_occupancy")
        gp120 = list(agg.items("residue:receptor:GP120:residue_occupancy"))
        cd4 = list(agg.items("residue:receptor:CD4:residue_occupancy"))
        self.assertEqual(len(gp120[0][2]), 3)
        self.assertEqual(len(cd4[0][2]), 1)

    def test_every_group_is_registered_before_anything_is_added(self):
        """The registration that makes the zero-fill possible. A group created
        only when its first value arrives has only the units that had values,
        and the fill would then be a no-op exactly when it matters."""
        groups = interface.pooled_groups(
            ["A120", "B120", "A41"], ["IFH", "IFL"],
            {"A120": "A120", "B120": "A120", "A41": "A41", "IFH": "IFH",
             "IFL": "IFL"})
        self.assertIn("pair", groups)
        self.assertIn("run", groups)
        self.assertIn("residue:receptor:A120:residue_occupancy", groups)
        self.assertIn("residue:receptor:A41:partners_per_frame", groups)
        self.assertIn("residue:ligand:IFL:residue_occupancy", groups)
        # one group per (side, entity, observable), and gp120's two copies share
        # one group rather than getting one each
        self.assertEqual(len(groups), 2 + 2 * 2 + 2 * 2)

    def test_keys_come_back_in_a_fixed_order(self):
        """Determinism: same input, byte-identical TSV. A dict iteration order
        that drifted would reorder the rows and nothing would look wrong."""
        agg = interface.Pooled()
        for key in (("E", 400, "H", 9, 0), ("E", 384, "H", 100, 0),
                    ("E", 384, "H", 9, 0)):
            agg.add("pair", key, "r1", "A", 0.5, observable="x")
        self.assertEqual([k for k, _, _, _ in agg.items("pair")],
                         [("E", 384, "H", 9, 0), ("E", 384, "H", 100, 0),
                          ("E", 400, "H", 9, 0)])


# --------------------------------------------------------------------------
# lifetimes and the sampling floor
# --------------------------------------------------------------------------
class TestLifetimes(unittest.TestCase):
    def test_episodes_come_from_the_shared_run_length_kernel(self):
        """A contact broken for one frame is two contacts, not one long one.
        geom.run_lengths is what every lifetime in this toolkit uses, so contact
        lifetime, ion dwell and water residence cannot disagree about that."""
        mask = np.array([1, 1, 1, 0, 1, 1, 0, 0, 1], dtype=bool)
        got = [int(n) for _, n in geom.run_lengths(mask)]
        self.assertEqual(got, [3, 2, 1])
        st = interface.lifetime_stats(got)
        self.assertEqual(st["n"], 3)
        self.assertAlmostEqual(st["median"], 2.0)

    def test_a_distribution_at_the_sampling_floor_is_not_resolved(self):
        """Three of four episodes over before the next sampled frame: the mean
        is the frame spacing, and a table that reports it as a lifetime is worse
        than one that reports nothing."""
        st = interface.lifetime_stats([1, 1, 1, 5])
        ok, note = interface.lifetime_floor(st, 1.0)
        self.assertFalse(ok)
        self.assertIn("single sampled frame", note)
        self.assertIn("1 ns", note)

    def test_a_resolved_distribution_carries_no_note(self):
        ok, note = interface.lifetime_floor(
            interface.lifetime_stats([4, 5, 6, 1]), 0.2)
        self.assertTrue(ok)
        self.assertIsNone(note)

    def test_occupancy_has_no_floor_and_the_lifetime_does(self):
        """The rule this stage is built on: a time average needs no particular
        spacing, so MIN_SPACING_NS is None and a coarse stride never kills the
        map -- but anything with dynamic content declares one."""
        self.assertIsNone(interface.MIN_SPACING_NS)
        self.assertTrue(interface.LIFETIME_MIN_SPACING_NS)
        # one frame per 100 ns still passes for the map, and must
        self.assertEqual(SP.check_spacing("interface_contacts", 0.2, 500,
                                          interface.MIN_SPACING_NS), 100.0)
        # the same spacing is refused for anything with dynamic content
        with self.assertRaises(SP.SpacingTooCoarse):
            SP.check_spacing("contact_lifetime", 0.2, 500,
                             interface.LIFETIME_MIN_SPACING_NS)

    def test_the_spacing_verdict_survives_into_the_note(self):
        """At DEFAULT_STRIDE, which is chosen for the map, the lifetimes are NOT
        resolved -- and they say so in their own rows rather than being dropped,
        because dropping them would cost the occupancies the same walk paid
        for."""
        with self.assertRaises(SP.SpacingTooCoarse) as cm:
            SP.check_spacing("contact_lifetime", 0.2, interface.DEFAULT_STRIDE,
                             interface.LIFETIME_MIN_SPACING_NS)
        ok, note = interface.lifetime_floor(
            interface.lifetime_stats([4, 5, 6]), 1.0, str(cm.exception))
        self.assertFalse(ok)
        self.assertIn("0.5 ns or finer", note)


class TestInputs(unittest.TestCase):
    def test_this_stage_reads_the_stripped_trajectory(self):
        """An interface is made of solute. Reading the solvated file would cost
        thirteen times the I/O for the same geometry, and `inputs(solvated=...)`
        is how a stage says which it needs instead of guessing."""
        md = SP.MD(id="1R", topology="merged.psf", trajectories=["merged.dcd"],
                   solvated_topology="solv.psf",
                   solvated_trajectories=["solv.dcd"])
        top, trajs = md.inputs(solvated=False)
        self.assertEqual((top, trajs), ("merged.psf", ["merged.dcd"]))

    def test_the_extra_columns_are_declared_so_the_header_is_stable(self):
        """Two tables whose columns depend on what a system happened to hold
        cannot be concatenated without a realignment nobody would notice."""
        from core import io
        rows = io.Rows(system="x", stage=interface.STAGE,
                       **interface.EXTRA_COLUMNS)
        rows.add(observable="contact_occupancy", value=0.5)
        header = rows.to_tsv().splitlines()[0].split("\t")
        for col in ("partner_segid", "partner_resid", "partner_resname",
                    "entity", "partner_entity", "protomer", "protomer_offset",
                    "cutoff_A", "sd", "sem", "n"):
            self.assertIn(col, header)


if __name__ == "__main__":
    unittest.main(verbosity=2)
