#!/usr/bin/env python3
"""
test_default_production.py -- a silent .dcd is production until its author
says otherwise, and that default never overrules evidence.

    python3 tests/test_default_production.py     # or: python3 -m pytest tests/

Decided 2026-09-21. Before, a .dcd whose name said nothing was `unresolved` and
never counted, so ABZYME/1-FULL-GP120/MDs (five 32 GB replicas) produced an
empty glossary. Each test pins one edge of the new rule. Fixtures are synthetic:
tiny DCD headers written into a temp directory, nothing read from real data.
"""

import os
import struct
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from core import stages                              # noqa: E402
from core.mds import replica_split                   # noqa: E402
from core import config                                   # noqa: E402
from core.trajscan import classify_dir, frame_info        # noqa: E402


def write_dcd(path, n_frames=10, natom=100, istart=0, nsavc=500, delta=0.04091):
    """A CHARMM-flavoured DCD: header, then `n_frames` frames of zeroes.

    The frames used to be left out -- until the reader started counting what a
    file really holds (a transfer in flight announces frames it has not got
    yet), at which point every fixture looked truncated. They cost 1.3 kB each
    at the default 100 atoms."""
    icntrl = [n_frames, istart, nsavc, n_frames * nsavc, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 24]
    head = struct.pack("<i4s", 84, b"CORD") + struct.pack("<9i", *icntrl[:9])
    head += struct.pack("<f", delta) + struct.pack("<10i", *icntrl[10:])
    head += struct.pack("<i", 84)
    title = b"test".ljust(80)
    with open(path, "wb") as fh:
        fh.write(head)
        fh.write(struct.pack("<ii", 84, 1) + title + struct.pack("<i", 84))
        fh.write(struct.pack("<iii", 4, natom, 4))
        # unit cell (icntrl[10] is 1 above) + x, y, z, each a marked block
        frame = struct.pack("<i", 48) + b"\0" * 48 + struct.pack("<i", 48)
        frame += (struct.pack("<i", natom * 4) + b"\0" * (natom * 4)
                  + struct.pack("<i", natom * 4)) * 3
        for _ in range(max(0, n_frames)):
            fh.write(frame)


def rules(**stages_):
    return config.compile_stage_rules({"stages": stages_})


class TestDefault(unittest.TestCase):
    def test_a_silent_dcd_alone_is_production_by_default(self):
        names = [f"base_model.lipids_protein_glycans.R{i}.dcd" for i in range(1, 6)]
        out = classify_dir("/x/MDs", names)
        self.assertEqual({v for v in out.values()}, {("production", "default")})

    def test_names_that_say_something_keep_it(self):
        out = classify_dir("/x/run", ["heating.dcd", "eq1.dcd", "pr6900frames.dcd",
                                      "notes.dcd"])
        self.assertEqual(out["heating.dcd"], ("heating", "name"))
        self.assertEqual(out["eq1.dcd"], ("equilibration", "name"))
        self.assertEqual(out["pr6900frames.dcd"], ("derived", "name"))
        self.assertEqual(out["notes.dcd"], ("production", "default"))
        out = classify_dir("/x/SASA/stripped", ["1R.pr.dcd", "all.dcd"])
        self.assertEqual(out["1R.pr.dcd"], ("stripped", "name"))

    def test_never_beside_files_already_named_production(self):
        """pr.dcd ... pr7.dcd are the run; a silent 4092-5091.dcd beside them is a
        frame-range extract, and appending it would repeat frames."""
        out = classify_dir("/x/pr", ["pr.dcd", "pr2.dcd", "4092-5091.dcd"])
        self.assertEqual(out["4092-5091.dcd"], ("unresolved", "name"))

    def test_a_directory_named_for_another_stage_is_evidence(self):
        out = classify_dir("/x/PENTAMER_REBUILD/07_EQUILIBRATION_1", ["IF4_rebuilt_omm.dcd"])
        self.assertEqual(out["IF4_rebuilt_omm.dcd"], ("unresolved", "name"))
        out = classify_dir("/x/heating/", ["a.dcd"])
        self.assertEqual(out["a.dcd"], ("unresolved", "name"))

    def test_only_dcd(self):
        out = classify_dir("/x/amber", ["250ns.mdcrd", "traj.nc"])
        self.assertEqual({v for v in out.values()}, {("unresolved", "name")})

    def test_the_author_wins_over_name_and_default(self):
        r = rules(derived=r"/RMSF/all\.dcd$", production=r"4092-5091\.dcd$")
        out = classify_dir("/x/RMSF", ["all.dcd", "two.dcd"], r, "HD2/RMSF")
        self.assertEqual(out["all.dcd"], ("derived", "author"))
        self.assertEqual(out["two.dcd"], ("production", "default"))
        out = classify_dir("/x/pr", ["pr.dcd", "4092-5091.dcd"], r, "a/pr")
        self.assertEqual(out["4092-5091.dcd"], ("production", "author"))

    def test_author_production_is_production_for_the_default_too(self):
        r = rules(production=r"keep\.dcd$")
        out = classify_dir("/x/d", ["keep.dcd", "other.dcd"], r, "d")
        self.assertEqual(out["other.dcd"], ("unresolved", "name"))

    def test_an_unknown_stage_stops_the_build(self):
        with self.assertRaises(ValueError):
            rules(prodution=".*")

    def test_the_export_agrees_with_the_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = os.path.join(tmp, "a", "MDs")
            os.makedirs(run)
            for n in ("R1.dcd", "R2.dcd", "eq.dcd"):
                write_dcd(os.path.join(run, n))
            clf = stages.Classifier(rules(derived=r"R2\.dcd$"), tmp)
            self.assertEqual(clf("R1.dcd", run), "production")
            self.assertEqual(clf.basis("R1.dcd", run), "default")
            self.assertEqual(clf("R2.dcd", run), "derived")
            self.assertEqual(clf("eq.dcd", run), "equilibration")


class TestSubsetsAndAttempts(unittest.TestCase):
    def test_a_solute_only_subset_is_not_production(self):
        """2-REPLICA/1R: md.dcd (45,210 atoms) beside pr.dcd and pr_r1.dcd
        (711,846). It was ordered first and opened three published MDs."""
        names = ["md.dcd", "pr.dcd", "pr_r1.dcd", "md_merged.dcd", "npt2.dcd"]
        natoms = {"md.dcd": 45210, "pr.dcd": 711846, "pr_r1.dcd": 711846,
                  "md_merged.dcd": 45210, "npt2.dcd": 711846}
        out = classify_dir("/x/1R", names, natoms=natoms)
        self.assertEqual(out["md.dcd"], ("stripped", "atoms"))
        self.assertEqual(out["pr.dcd"], ("production", "name"))
        self.assertEqual(out["pr_r1.dcd"], ("production", "name"))
        self.assertEqual(out["md_merged.dcd"], ("merged", "name"))
        self.assertEqual(out["npt2.dcd"], ("equilibration", "name"))

    def test_the_author_can_keep_a_smaller_system(self):
        r = rules(production=r"small\.dcd$")
        out = classify_dir("/x/d", ["pr.dcd", "small.dcd"], r, "d",
                           natoms={"pr.dcd": 1000, "small.dcd": 10})
        self.assertEqual(out["small.dcd"], ("production", "author"))

    def test_a_named_subset_leaves_room_for_the_default(self):
        out = classify_dir("/x/d", ["md.dcd", "system.dcd"],
                           natoms={"md.dcd": 10, "system.dcd": 1000})
        self.assertEqual(out["md.dcd"], ("stripped", "atoms"))
        self.assertEqual(out["system.dcd"], ("production", "default"))

    def test_a_failed_attempt_is_superseded(self):
        out = classify_dir("/x/1R/failed_heating_nvt_2026-07-30", ["heating.dcd"])
        self.assertEqual(out["heating.dcd"], ("superseded", "name"))


class TestPlaceholderSpacing(unittest.TestCase):
    def test_a_tool_written_header_has_no_spacing(self):
        """DELTA 1.0 with NSAVC 1 is VMD's placeholder: 10,500 frames would read
        as 0.5 ns."""
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "R1.dcd")
            write_dcd(p, n_frames=10500, nsavc=1, delta=1.0)
            nf, na, dt, _ = frame_info(p)
            self.assertEqual((nf, na, dt), (10500, 100, None))
            write_dcd(p, n_frames=100, nsavc=25000, delta=0.0818193182349205)
            self.assertAlmostEqual(frame_info(p)[2], 100.0, places=3)


class TestReplicaSplit(unittest.TestCase):
    def make(self, tmp, specs):
        for name, istart in specs:
            write_dcd(os.path.join(tmp, name), istart=istart)
        return [n for n, _ in specs]

    def test_side_by_side_replicas_are_separate_mds(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = self.make(tmp, [(f"base.R{i}.dcd", 0) for i in (1, 2, 3, 4, 5)])
            groups = replica_split(tmp, "", files)
            self.assertEqual(groups, [(i, [f"base.R{i}.dcd"]) for i in (1, 2, 3, 4, 5)])

    def test_restart_segments_are_not_replicas(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = self.make(tmp, [("pr.dcd", 0), ("pr2.dcd", 500), ("pr3.dcd", 1000)])
            self.assertIsNone(replica_split(tmp, "", files))
            files = self.make(tmp, [("pr_r1.dcd", 500), ("pr_r2.dcd", 1000)])
            self.assertIsNone(replica_split(tmp, "", files),
                              "labels that continue one another are segments")
            files = self.make(tmp, [("prh01.dcd", 0), ("prh02.dcd", 0)])
            self.assertIsNone(replica_split(tmp, "", files))

    def test_amber_replicas_are_told_from_segments_by_their_time(self):
        """DESIGN-V3's md_rep2.nc and md_rep3.nc both run 105-600 ps: replicas,
        which had been one MD. Files that continue one another stay one MD."""
        from scipy.io import netcdf_file
        with tempfile.TemporaryDirectory() as tmp:
            def nc(name, t0):
                with netcdf_file(os.path.join(tmp, name), "w") as f:
                    f.createDimension("frame", None)
                    t = f.createVariable("time", "f4", ("frame",))
                    t[:] = [t0, t0 + 5.0]
            nc("md_rep2.nc", 105.0)
            nc("md_rep3.nc", 105.0)
            self.assertEqual(replica_split(tmp, "", ["md_rep2.nc", "md_rep3.nc"]),
                             [(2, ["md_rep2.nc"]), (3, ["md_rep3.nc"])])
            nc("md_rep3.nc", 610.0)
            self.assertIsNone(replica_split(tmp, "", ["md_rep2.nc", "md_rep3.nc"]))

    def test_anything_short_of_every_condition_is_one_md(self):
        with tempfile.TemporaryDirectory() as tmp:
            cases = [[("run.R1.dcd", 0), ("run.dcd", 0)],           # one unlabelled
                     [("a.R1.dcd", 0), ("b.R2.dcd", 0)],            # names differ
                     [("x.R1.dcd", 0), ("x.1R.dcd", 0)]]            # same replica twice
            for specs in cases:
                self.assertIsNone(replica_split(tmp, "", self.make(tmp, specs)), specs)


if __name__ == "__main__":
    unittest.main(verbosity=2)
