#!/usr/bin/env python3
"""
test_core.py -- one test per defect that has already cost something.

    python3 tests/test_core.py

Synthetic fixtures only: nothing here reads the lab's data, and every file it
writes lives in a temporary directory. The point is not coverage, it is that
each of these failures was silent -- a wrong answer with no error -- and a
silent failure only stays fixed if a test fails when it returns.
"""

import os
import struct
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from core import frames as F                               # noqa: E402
from core import pbc, spec, topology as T                  # noqa: E402
from core import selections as S                           # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def write_psf(path, atoms):
    """A minimal CHARMM PSF: header, !NATOM block, and the bonds given."""
    lines = ["PSF", "", "       1 !NTITLE", " REMARKS synthetic", "",
             f"{len(atoms):>8} !NATOM"]
    for i, (segid, resid, resname, name, mass) in enumerate(atoms, 1):
        lines.append(f"{i:>8} {segid:<5} {resid:<5} {resname:<5} {name:<5} "
                     f"{'X':<5} {0.0:>10.6f} {mass:>13.4f}{0:>12}")
    lines.append("")
    lines.append(f"{0:>8} !NBOND: bonds")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def write_psf_with_bonds(path, atoms, bonds):
    write_psf(path, atoms)
    text = open(path).read().replace(f"{0:>8} !NBOND: bonds",
                                     f"{len(bonds):>8} !NBOND: bonds")
    flat = []
    for a, b in bonds:
        flat += [f"{a:>8}", f"{b:>8}"]
    rows = ["".join(flat[i:i + 8]) for i in range(0, len(flat), 8)]
    with open(path, "w") as fh:
        fh.write(text + "\n".join(rows) + "\n")
    return path


def write_dcd(path, n_frames=10, natom=6, istart=0, nsavc=500, delta=0.04091):
    """A CHARMM-flavoured DCD with real (zeroed) frames, unit cell included."""
    icntrl = [n_frames, istart, nsavc, n_frames * nsavc, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 24]
    head = struct.pack("<i4s", 84, b"CORD") + struct.pack("<9i", *icntrl[:9])
    head += struct.pack("<f", delta) + struct.pack("<10i", *icntrl[10:])
    head += struct.pack("<i", 84)
    with open(path, "wb") as fh:
        fh.write(head)
        fh.write(struct.pack("<ii", 84, 1) + b"synthetic".ljust(80)
                 + struct.pack("<i", 84))
        fh.write(struct.pack("<iii", 4, natom, 4))
        frame = struct.pack("<i", 48) + b"\0" * 48 + struct.pack("<i", 48)
        frame += (struct.pack("<i", natom * 4) + b"\0" * (natom * 4)
                  + struct.pack("<i", natom * 4)) * 3
        for _ in range(n_frames):
            fh.write(frame)
    return path


TRIMER = (
    [("A120", r, "ALA", "CA", 12.0) for r in range(1, 6)] +
    [("A41", r, "GLY", "CA", 12.0) for r in range(1, 4)] +
    [("GA1", 1, "BGLC", "C1", 12.0), ("GA1", 2, "AMAN", "C1", 12.0)] +
    [("IFH", r, "VAL", "CA", 12.0) for r in range(1, 4)] +
    [("WT1", r, "TIP3", "OH2", 16.0) for r in range(1, 21)] +
    [("ION", 1, "SOD", "SOD", 23.0), ("ION", 2, "CLA", "CLA", 35.5)]
)


# --------------------------------------------------------------------------
class TestSelections(unittest.TestCase):
    """The translator, whose failures are all quiet ones."""

    def test_every_range_is_expanded_not_just_the_first(self):
        """`segid GA1 to GA23 GB1 to GB23` matched 33 of 75 segments on the
        real 1F4 topology: the first range expanded, the rest were read as
        literal segids called `to`. Not empty, so nothing complained."""
        out = S.translate("segid GA1 to GA3 GB1 to GB3 GAE")
        self.assertEqual(out, "segid GA1 GA2 GA3 GB1 GB2 GB3 GAE")

    def test_an_unexpanded_range_is_an_error(self):
        with self.assertRaises(ValueError):
            S.translate("segid GA1 to GB3")      # different stems, no meaning

    def test_within_keeps_the_reference_group(self):
        """VMD's `within` includes the reference set; MDAnalysis's `around`
        excludes it. Translating one to the other loses every atom of it."""
        out = S.translate("water and within 3.5 of (segid A120)")
        self.assertIn("around 3.5 (segid A120)", out)
        self.assertIn("or (segid A120)", out)

    def test_exwithin_is_the_exclusive_form(self):
        out = S.translate("ions and exwithin 4 of glycan")
        self.assertIn("around 4 (", out)
        self.assertNotIn("or (resname", out.split("around 4")[1][:40])

    def test_a_condition_after_within_is_not_swallowed(self):
        out = S.translate("water and within 3.5 of (segid A120) and name OH2")
        self.assertTrue(out.rstrip().endswith("and name OH2"))

    def test_noh_is_by_name_never_by_mass(self):
        self.assertIn("not name H*", S.translate("protein and noh"))

    def test_macros_mean_the_same_atoms_as_in_vmd(self):
        self.assertIn("TIP3", S.translate("water"))
        self.assertIn("SOD", S.translate("ions"))
        self.assertIn("AMAN", S.translate("glycan"))

    def test_an_uppercase_residue_name_is_not_a_macro(self):
        """`ION`, `MAN` and `CL` are residue names; expanding them as macros
        would turn one selection into a recursive mess."""
        self.assertEqual(S.translate("resname ION MAN CL"), "resname ION MAN CL")


class TestTopology(unittest.TestCase):
    def test_composition_counts_molecules_not_atoms(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_psf(os.path.join(d, "s.psf"), TRIMER)
            atoms = T.read_psf(p)
            comp = T.composition(atoms)
            self.assertEqual(comp["water_molecules"], 20)
            self.assertEqual(comp["ions_by_name"], {"CLA": 1, "SOD": 1})
            self.assertEqual(comp["glycan_residues"], 2)
            self.assertEqual(comp["solute_atoms"], len(TRIMER) - 22)

    def test_copies_of_one_chain_are_one_entity(self):
        atoms = [("A", r, "ALA", "CA", 12.0) for r in range(1, 6)] + \
                [("B", r, "ALA", "CA", 12.0) for r in range(1, 6)] + \
                [("H", r, "VAL", "CA", 12.0) for r in range(1, 4)]
        with tempfile.TemporaryDirectory() as d:
            p = write_psf(os.path.join(d, "s.psf"), atoms)
            ents = T.entities(T.segments(T.read_psf(p)))
            self.assertEqual(ents[0].copies, 2)
            self.assertEqual(sorted(ents[0].segids), ["A", "B"])
            self.assertEqual(ents[1].segids, ["H"])

    def test_glycan_links_say_which_chain_carries_them(self):
        """The evidence that keeps gp41 in the receptor: it carries glycans,
        and a Fab of the same size does not."""
        with tempfile.TemporaryDirectory() as d:
            p = write_psf_with_bonds(os.path.join(d, "s.psf"), TRIMER,
                                     [(1, 9), (6, 10)])   # A120-GA1, A41-GA1
            atoms = T.read_psf(p)
            segs = T.segments(atoms)
            hosts = T.glycan_hosts(T.bonded_segments(p, atoms), segs)
            self.assertEqual(set(hosts), {"A120", "A41"})
            self.assertNotIn("IFH", hosts)


class TestFrames(unittest.TestCase):
    def test_a_restart_seam_is_counted_once(self):
        with tempfile.TemporaryDirectory() as d:
            a = write_dcd(os.path.join(d, "pr.dcd"), n_frames=5,
                          istart=100000, nsavc=100000)
            b = write_dcd(os.path.join(d, "pr.part1.dcd"), n_frames=4,
                          istart=500000, nsavc=100000)
            t = F.build("1R", [a, b])
            self.assertEqual(len(t), 8)          # 9 raw frames, one repeated
            self.assertEqual(t.seams, [5])
            self.assertEqual(len(t) + len(t.seams), 9)

    def test_a_frame_knows_which_file_and_index_it_is(self):
        with tempfile.TemporaryDirectory() as d:
            a = write_dcd(os.path.join(d, "pr.dcd"), n_frames=5,
                          istart=100000, nsavc=100000)
            b = write_dcd(os.path.join(d, "pr.part1.dcd"), n_frames=4,
                          istart=600000, nsavc=100000)
            t = F.build("1R", [a, b])
            w = t.where(6)
            self.assertEqual(w["file"], "pr.part1.dcd")
            self.assertEqual(w["local_frame"], 1)
            self.assertEqual(w["raw_frame"], 6)

    def test_order_comes_from_the_headers_not_the_names(self):
        with tempfile.TemporaryDirectory() as d:
            late = write_dcd(os.path.join(d, "pr.part2.dcd"), n_frames=2,
                             istart=900000, nsavc=100000)
            early = write_dcd(os.path.join(d, "pr.part10.dcd"), n_frames=2,
                              istart=100000, nsavc=100000)
            t = F.build("1R", [late, early])
            self.assertEqual(t.segments[0].name, "pr.part10.dcd")

    def test_a_hole_in_the_timeline_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            a = write_dcd(os.path.join(d, "a.dcd"), n_frames=3,
                          istart=100000, nsavc=100000)
            b = write_dcd(os.path.join(d, "b.dcd"), n_frames=3,
                          istart=900000, nsavc=100000)
            g = F.gaps(F.build("1R", [a, b]))
            self.assertEqual(len(g), 1)
            self.assertEqual(g[0][0], 300000)


class TestFrameSpacing(unittest.TestCase):
    """Where "how long is a frame" is allowed to come from.

    Everything the pages show is plotted against time and every caption
    quotes one, so a wrong frame spacing is not a wrong axis label -- it is
    every number on the page. The failure below happened: the FULL-GP120
    production files are post-processed copies, cpptraj wrote them with
    `nsavc = 1` and a default delta, and a 1,050 ns replica came out as
    0.513 ns -- wrong by a factor of two thousand, with nothing raised, while
    the stage tables beside it wrote `time_ns = 0.0` for every frame. The two
    halves of the same bundle disagreed and neither complained.
    """

    def placeholder(self, d):
        """What cpptraj leaves behind: no real step counter, no real delta."""
        return F.build("md_R1", [write_dcd(os.path.join(d, "post.dcd"),
                                           n_frames=10, istart=0, nsavc=1)])

    def test_a_placeholder_is_refused_rather_than_used(self):
        with tempfile.TemporaryDirectory() as d:
            t = self.placeholder(d)
            self.assertIsNone(t.framestep_ns)
            self.assertIsNone(t.framestep_basis)
            self.assertIsNone(t.as_dict()["framestep_ns"])

    def test_real_headers_still_state_the_spacing_themselves(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_dcd(os.path.join(d, "pr.dcd"), n_frames=10,
                          istart=100000, nsavc=100000)
            t = F.build("1R", [p])
            self.assertEqual(t.framestep_basis, "headers")
            self.assertGreater(t.framestep_ns, 0.001)

    def test_a_declared_interval_times_frames_by_position(self):
        # The step counters are placeholders too, so they cannot be mixed in:
        # frame i is at i x the declared interval, and nowhere else.
        with tempfile.TemporaryDirectory() as d:
            t = F.build("md_R1", [write_dcd(os.path.join(d, "post.dcd"),
                                            n_frames=10, nsavc=1)],
                        framestep_ns=0.1)
            self.assertEqual(t.framestep_basis, "declared")
            self.assertAlmostEqual(t.time_ns(0), 0.0)
            self.assertAlmostEqual(t.time_ns(7), 0.7)
            self.assertAlmostEqual(t.duration_ns(), 1.0)

    def test_the_document_says_which_clock_was_used(self):
        # A reader of frames.json -- or of the provenance panel -- must be
        # able to tell a measured interval from one somebody typed.
        with tempfile.TemporaryDirectory() as d:
            t = F.build("md_R1", [write_dcd(os.path.join(d, "post.dcd"),
                                            n_frames=10, nsavc=1)],
                        framestep_ns=0.1)
            self.assertEqual(t.as_dict()["framestep_basis"], "declared")
            self.assertEqual(t.as_dict()["duration_ns"], 1.0)

    def test_a_real_run_is_timed_by_its_step_counter_not_its_position(self):
        """A restart must not reset the clock.

        Frame 0 of a segment that resumed at step 1,000,000 is not at time
        zero. This is why position-timing is only used where there is no
        real step counter at all.
        """
        with tempfile.TemporaryDirectory() as d:
            a = write_dcd(os.path.join(d, "a.dcd"), n_frames=5,
                          istart=100000, nsavc=100000)
            b = write_dcd(os.path.join(d, "b.dcd"), n_frames=5,
                          istart=600000, nsavc=100000)
            t = F.build("1R", [a, b])
            self.assertEqual(t.framestep_basis, "headers")
            self.assertGreater(t.time_ns(5), t.time_ns(4))
            self.assertGreater(t.time_ns(0), 0)      # not frame 0 of the run


class TestAnalysedIndexing(unittest.TestCase):
    """Which frame number names a frame inside the file you open.

    The provenance panel exists to answer exactly this, and it is the one
    thing on the page that cannot be checked by looking at it: "frame 2019 of
    2525" is either right or quietly off by one.

    A restart re-simulates its checkpoint frame. A strip either keeps that
    duplicate -- the CDR strips did, so `merged.dcd` is the RAW sequence --
    or drops it with `trajin ... 2 last`, as the occluded-open strips did, so
    their file is the DE-DUPLICATED sequence and is one frame shorter per
    seam. Assuming the first put every frame after a seam one off in the
    second.
    """

    def timeline(self, d):
        """Two segments that overlap by one frame: a restart seam."""
        a = write_dcd(os.path.join(d, "a.dcd"), n_frames=5, istart=0, nsavc=100)
        b = write_dcd(os.path.join(d, "b.dcd"), n_frames=5, istart=400,
                      nsavc=100)
        return F.build("1R", [a, b], framestep_ns=0.1)

    def test_a_strip_that_kept_the_duplicate_is_the_raw_sequence(self):
        with tempfile.TemporaryDirectory() as d:
            t = self.timeline(d)
            t.n_analysed = len(t.frames) + len(t.seams)
            self.assertEqual(t.analysed_indexing, "raw")

    def test_a_strip_that_dropped_it_is_the_unique_sequence(self):
        with tempfile.TemporaryDirectory() as d:
            t = self.timeline(d)
            t.n_analysed = len(t.frames)
            self.assertEqual(t.analysed_indexing, "unique")

    def test_a_seam_is_what_makes_the_two_differ(self):
        with tempfile.TemporaryDirectory() as d:
            t = self.timeline(d)
            # if this fixture stops overlapping, the test above proves nothing
            self.assertGreater(len(t.seams), 0)

    def test_a_length_matching_neither_is_reported_as_unknown(self):
        # Better to say the frame cannot be located than to name a wrong one:
        # the panel hides the VMD command rather than hand over a command
        # that opens the wrong picosecond.
        with tempfile.TemporaryDirectory() as d:
            t = self.timeline(d)
            t.n_analysed = len(t.frames) + 17
            self.assertIsNone(t.analysed_indexing)
            t.n_analysed = None
            self.assertIsNone(t.analysed_indexing)

    def test_the_document_carries_both_numbers(self):
        with tempfile.TemporaryDirectory() as d:
            t = self.timeline(d)
            t.analysed = ["merged.dcd"]
            t.analysed_path = os.path.join(d, "merged.dcd")
            t.n_analysed = len(t.frames)
            doc = t.as_dict()
            self.assertEqual(doc["n_analysed"], len(t.frames))
            self.assertEqual(doc["analysed_indexing"], "unique")
            self.assertNotEqual(doc["n_raw"], doc["n_analysed"])


class TestPBC(unittest.TestCase):
    def test_a_group_across_the_boundary_is_brought_back_whole(self):
        """The Fab that appears to fly away: `wrapAll on` writes it a box
        length from the trimer it is bound to."""
        import numpy as np
        L = np.array([200.0, 210.0, 180.0])
        X = np.array([[5.0, 5, 5], [6.0, 5, 5]])
        ref = np.array([195.0, 5, 5])
        Y, moved = pbc.shift_into_image(X, ref, L)
        self.assertEqual(moved, 1)
        self.assertAlmostEqual(float(np.linalg.norm(Y.mean(0) - ref)), 10.5, 3)
        # and the group is still whole: the bond length did not change
        self.assertAlmostEqual(float(np.linalg.norm(Y[1] - Y[0])), 1.0, 6)

    def test_segments_move_independently(self):
        import numpy as np
        L = np.array([100.0, 100, 100])
        X = np.array([[1.0, 1, 1], [2.0, 1, 1], [99.0, 1, 1]])
        g = np.array(["H", "H", "L"])
        Y, moved = pbc.shift_into_image(X, np.array([2.0, 1, 1]), L, groups=g)
        self.assertEqual(moved, 1)
        self.assertAlmostEqual(Y[2][0], -1.0, 6)
        self.assertAlmostEqual(Y[0][0], 1.0, 6)

    def test_a_triclinic_cell_is_refused_not_guessed(self):
        with self.assertRaises(pbc.NonOrthorhombic):
            pbc.box_lengths([100, 100, 100, 90, 60, 90])

    def test_a_split_group_is_detected(self):
        import numpy as np
        L = np.array([200.0, 200, 200])
        self.assertTrue(pbc.looks_split(np.array([[0.0, 0, 0], [150.0, 0, 0]]), L))
        self.assertFalse(pbc.looks_split(np.array([[0.0, 0, 0], [40.0, 0, 0]]), L))


class TestSpec(unittest.TestCase):
    def test_components_choose_the_bundle(self):
        self.assertEqual(spec.bundle_for(["protein"]), "apo-protein")
        self.assertEqual(spec.bundle_for(["protein", "glycans"]), "glycoprotein")
        self.assertEqual(spec.bundle_for(["protein", "glycans", "antibody"]),
                         "antibody-complex")

    def test_a_solvated_system_gets_the_solvent_stages(self):
        s = spec.Spec("x/1", components=["protein", "glycans", "water", "ions"])
        self.assertIn("composition", s.stages())
        self.assertIn("water_residence", s.stages())
        self.assertIn("rmsd", s.stages())

    def test_an_observable_the_sampling_cannot_support_is_refused(self):
        """At 0.2 ns/frame, stride 100 is one frame per 20 ns: fine for
        occupancy, meaningless for exchange dynamics."""
        with self.assertRaises(spec.SpacingTooCoarse):
            spec.check_spacing("water_residence", 0.2, 100, 0.1)
        self.assertEqual(spec.check_spacing("composition", 0.2, 100, None), 20.0)

    def test_a_solvent_stage_refuses_a_stripped_trajectory(self):
        md = spec.MD(id="1R", topology="m.psf", trajectories=["merged.dcd"])
        with self.assertRaises(ValueError):
            md.inputs(solvated=True)

    def test_the_time_axis_comes_from_the_original_segments(self):
        """cpptraj writes merged.dcd with a placeholder timestep, so a
        timeline built from it reports 504 ns as 0.003 ns."""
        md = spec.MD(id="1R", topology="m.psf", trajectories=["merged.dcd"],
                     solvated_topology="s.psf",
                     solvated_trajectories=["pr.dcd", "pr.part1.dcd"])
        self.assertEqual(md.timeline_inputs(), ["pr.dcd", "pr.part1.dcd"])


class TestSegidRanges(unittest.TestCase):
    """detect writes `segid G001 to G020`; the selection translator must read
    back exactly those segids. Padding was lost both ways (int() then `G1`), so
    a glycoprotein with zero-padded glycan segids matched no atom and every
    analysis of it was skipped."""

    def _round_trip(self, segids):
        from core import detect
        text = S.translate("segid " + " ".join(detect._compress(segids)))
        return set(text.split()[1:])

    def test_padded_unpadded_gapped_and_mixed_families_round_trip(self):
        for segids in (
                [f"G{i:03d}" for i in range(1, 21)],          # G001..G020
                [f"GA{i}" for i in range(1, 24)],             # GA1..GA23
                ["G001", "G002", "G005"],                     # padded, gapped
                [f"G{i:03d}" for i in range(1, 5)] + ["G7", "G8", "G9", "GAE"]):
            with self.subTest(segids=segids[:3]):
                self.assertEqual(self._round_trip(segids), set(segids))


class TestPartialStages(unittest.TestCase):
    def test_a_later_stage_s_failure_still_reports_the_stages_that_wrote(self):
        """The message said "nothing was written" after earlier stages had."""
        import analyse
        from types import SimpleNamespace
        from stages.solvent import composition
        from stages.structure import basic
        saved = composition.run, basic.run
        composition.run = lambda spec, data: ["row"]

        def boom(*a, **k):
            raise RuntimeError("structure died")
        basic.run = boom
        try:
            done = {}
            with self.assertRaises(RuntimeError):
                analyse.run_stages(SimpleNamespace(has_solvent=True), "/nowhere",
                                   ["composition", "structure"], done=done)
            self.assertEqual(list(done), ["composition"])
        finally:
            composition.run, basic.run = saved


class TestManifestInputs(unittest.TestCase):
    def test_replicas_whose_files_share_a_name_are_each_recorded(self):
        """Keyed by basename, R1/merged.dcd and R2/merged.dcd were one entry."""
        from core import io
        from types import SimpleNamespace
        sp = SimpleNamespace(id="s", bundle="b", selections={},
                             as_dict=lambda: {"components": [], "bundle": "b",
                                              "selections": {}, "stride": {},
                                              "mds": []})
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for r in ("R1", "R2"):
                os.makedirs(os.path.join(tmp, r))
                p = os.path.join(tmp, r, "merged.dcd")
                open(p, "wb").write(b"x" * (len(paths) + 1))
                paths.append(p)
            ins = io.manifest("structure", sp, paths)["inputs"]
            self.assertEqual(sorted(ins), ["R1/merged.dcd", "R2/merged.dcd"])
            self.assertNotEqual(ins["R1/merged.dcd"], ins["R2/merged.dcd"])

    def test_a_failed_parquet_does_not_leave_the_old_one_beside_the_new_tsv(self):
        from core import io
        with tempfile.TemporaryDirectory() as tmp:
            rows = io.Rows(system="s")
            rows.add(value=1.0)
            io.write(tmp, "x", rows)
            self.assertTrue(os.path.exists(os.path.join(tmp, "x.parquet")))
            bad = io.Rows(system="s")
            bad.add(value=object())                  # cannot be written as parquet
            io.write(tmp, "x", bad)
            self.assertFalse(os.path.exists(os.path.join(tmp, "x.parquet")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
