#!/usr/bin/env python3
"""
test_records.py -- the record a person reads: its name, and the catalogue row.

    python3 tests/test_records.py        # or: python3 -m pytest tests/

The complaint these answer, 2026-09-21: "the final .json file has a weird name,
which is not what I want" -- it was `1a918746.json`, a topology hash. And, the
next day, one folder holding four systems, each of which must get its own file
rather than four names that collide.
"""

import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
# The lab the tests run as, whatever the kit's own lab.yaml says.
os.environ["MD_PIPELINE_LAB"] = os.path.join(HERE, "lab_fixture.yaml")
os.environ.pop("NDP_ID_PREFIX", None)

import ndp_build                                             # noqa: E402
from core import catalog, layout, ndp                        # noqa: E402

CDR = "/data/lab/ABZYME/4-OPEN/CDR"
FAB = ["1-1F4/0-MD/1-R1", "1-1F4/0-MD/2-R2", "1-1F4/0-MD/3-R3"]
APO = ["MD/REPLICA/1R", "MD/REPLICA/2R", "MD/REPLICA/3R"]
SAMPLING = {"replicas": 3, "nominal_ns": 500}


class TestNames(unittest.TestCase):
    def test_four_systems_of_one_folder_get_four_names(self):
        names = {ndp.file_name(f"CDR/{gid}", {}, CDR, SAMPLING, paths)
                 for gid, paths in (
                     ("d2830a2a", FAB),
                     ("23d640e6", ["2-JL427/0-MD/REPLICA/1R",
                                   "2-JL427/0-MD/REPLICA/2R"]),
                     ("2ae05d56", ["3-YZ23/0-MD/REPLICA/1R",
                                   "3-YZ23/0-MD/REPLICA/2R"]),
                     ("e8232418", APO))}
        self.assertEqual(names, {"CDR_1-1F4_3x500ns.json",
                                 "CDR_2-JL427_3x500ns.json",
                                 "CDR_3-YZ23_3x500ns.json",
                                 "CDR_MD_3x500ns.json"})

    def test_the_name_says_the_folder_the_system_and_the_sampling(self):
        self.assertEqual(ndp.file_name("CDR/d2830a2a", {}, CDR, SAMPLING, FAB),
                         "CDR_1-1F4_3x500ns.json")
        self.assertEqual(
            ndp.file_name("CDR/x", {}, CDR, {"replicas": 3, "nominal_ns": 1000}, APO),
            "CDR_MD_3x1us.json")
        self.assertEqual(
            ndp.file_name("CDR/x", {}, CDR, {"replicas": 5}, APO),
            "CDR_MD_5-replicas.json")

    def test_a_title_does_not_rename_a_record_but_a_slug_does(self):
        """Writing a title must not move the file.

        It used to. On 2026-09-24 filling in the four CDR titles renamed every
        record in the folder -- `CDR_1-1F4_3x500ns.json` became
        `HIV---Env-Open-State---Tribound-1F4-antibody_3x500ns.json` (` - `
        slugs to `---`), which orphaned the analysis bundle and the production
        files keyed on the old stem and left four stale records beside the new
        ones. A title is a sentence for people to read; a file name is an
        address, and an address that changes when you describe the thing is
        not an address. Only `slug:` renames, because a person types it.
        """
        self.assertEqual(
            ndp.file_name("CDR/d2830a2a", {"title": "Env OPEN + 1F4 Fab"},
                          CDR, SAMPLING, FAB),
            "CDR_1-1F4_3x500ns.json")
        self.assertEqual(
            ndp.file_name("CDR/d2830a2a", {"slug": "hiv env 1f4"}, CDR,
                          SAMPLING, FAB),
            "hiv-env-1f4_3x500ns.json")

    def test_a_record_keeps_the_name_it_already_has(self):
        """Whatever it was called before, it is still called that.

        The name a record was published under is the one people and other
        tools hold; the builder passes the previous index entry's file name
        so a rule change cannot silently re-address an existing record. An
        explicit `slug:` still overrides it -- that is someone asking.
        """
        self.assertEqual(
            ndp.file_name("CDR/d2830a2a", {}, CDR, SAMPLING, FAB,
                          existing="something_older.json"),
            "something_older.json")
        self.assertEqual(
            ndp.file_name("CDR/d2830a2a", {"slug": "renamed"}, CDR, SAMPLING,
                          FAB, existing="something_older.json"),
            "renamed_3x500ns.json")

    def test_a_directory_that_only_says_where_still_names_the_record(self):
        """`MD/REPLICA` is generic all the way down; the record still needs a
        name, and `MD` is what tells it from `1-1F4`."""
        self.assertEqual(ndp.common_subpath(APO), "MD")
        self.assertEqual(ndp.common_subpath(FAB), "1-1F4")
        self.assertEqual(ndp.common_subpath(["3-YZ23/0-MD/REPLICA/1R"]), "3-YZ23")


class TestLengthFromDeclaredSpacing(unittest.TestCase):
    """A length the trajectory headers cannot state, from settings.yaml.

    2026-09-22, 1-FULL-GP120: five 10,500-frame DCDs written by a stripping
    tool, no `delta` in any header, so the record said `5x-unknown` and every
    length was null. Declaring the frame interval must recover them -- and must
    say that is where they came from.
    """

    def sampling(self, **kw):
        base = {"replicas": 5, "frame_interval_ps": 100, "total_ns": None,
                "length_per_replica_ns": None, "replica_lengths_ns": [],
                "replicas_detail": [{"id": f"md_R{i}", "length_ns": None,
                                     "frames": 10500} for i in range(1, 6)]}
        base.update(kw)
        return base

    def test_frames_times_the_declared_interval(self):
        s = ndp.lengths_from_spacing(self.sampling())
        self.assertEqual(s["total_ns"], 5250.0)
        self.assertEqual(s["length_per_replica_ns"], 1050.0)
        self.assertEqual(s["length_basis"], "declared frame interval")
        self.assertEqual(s["replicas_detail"][0]["length_ns"], 1050.0)
        self.assertEqual(s["replicas_detail"][0]["length_basis"],
                         "declared frame interval")

    def test_a_measured_length_is_never_replaced(self):
        s = ndp.lengths_from_spacing(self.sampling(total_ns=1395.0))
        self.assertEqual(s["total_ns"], 1395.0)
        self.assertNotIn("length_basis", s)

    def test_the_prose_says_what_the_record_says(self):
        """A record must not contradict its own description.

        `lengths_from_spacing` fills the RECORD's sampling, but the
        registration sheet and the auto-description were composed from the
        raw `measured` dict, which is still empty. On 2026-09-24
        `1-FULL-GP120`'s JSON said `total_ns: 5250.0` while its own
        description said "Sampling: 5 replica(s)" and the sheet printed
        `total ns` as a dash -- in the box a person pastes into the NDP form.
        """
        measured = {"replicas": 5, "replica_lengths_ns": [],
                    "length_per_replica_ns": None, "total_ns": None}
        record = {
            "title": "HIV - Env Closed State",
            "system": {"construct": None, "starting_structure": None,
                       "glycan_profile": None, "total_atoms": 256522},
            "setup": {"force_field": None, "water_model": None, "ions": None,
                      "software": None, "ensemble": None, "temperature_K": None,
                      "pressure_bar": None, "timestep_fs": None},
            "sampling": ndp.lengths_from_spacing({
                "replicas": 5, "frame_interval_ps": 100,
                "replica_lengths_ns": [], "length_per_replica_ns": None,
                "total_ns": None,
                "replicas_detail": [{"frames": 10500} for _ in range(5)]}),
            "protocol": [], "related_publication": None,
            "data_availability": None, "contact": None,
        }
        # the record knows its own length ...
        self.assertEqual(record["sampling"]["total_ns"], 5250.0)
        # ... and so, now, does every sentence written about it
        self.assertEqual(ndp.sampling_phrase(ndp._sampling_of(record, measured)),
                         "5 × 1050 ns")
        self.assertIn("5 × 1050 ns", ndp.catalog_title(record, measured))
        self.assertIn("one frame every 100 ps",
                      ndp.catalog_description(record, measured))
        self.assertNotIn("5 replica(s)",
                         ndp.catalog_description(record, measured))

    def test_a_partial_sampling_block_falls_back_to_the_measurement(self):
        """`describe()` composes from a stub whose sampling is not the real one.

        Reading a replica count out of that stub raised KeyError on every
        build, so the prose falls back to `measured` unless the block it is
        handed is the whole thing.
        """
        measured = {"replicas": 3, "replica_lengths_ns": [500.0] * 3,
                    "length_per_replica_ns": 500.0, "total_ns": 1500.0}
        stub = {"title": "x", "sampling": {"frame_interval_ps": 200}}
        self.assertEqual(ndp._sampling_of(stub, measured), measured)
        self.assertIn("3 × 500 ns", ndp.catalog_title(stub, measured))

    def test_no_interval_declared_leaves_it_unknown(self):
        s = ndp.lengths_from_spacing(self.sampling(frame_interval_ps=None))
        self.assertIsNone(s["total_ns"])

    def test_the_recovered_length_names_the_file(self):
        self.assertEqual(
            ndp.file_name("MDs/x", {}, "/data/1-FULL-GP120", self.sampling(), ["."]),
            "1-FULL-GP120_5-replicas.json")
        self.assertEqual(
            ndp.file_name("MDs/x", {}, "/data/1-FULL-GP120",
                          ndp.lengths_from_spacing(self.sampling()), ["."]),
            "1-FULL-GP120_5x1.05us.json")


class TestWhatIsPublished(unittest.TestCase):
    """Only the production phase goes out.

    2026-09-22: "min1.dcd, min2.dcd, heating.dcd etc are also recognized in the
    file, which is something that we don't want to be published". They reached
    the record through `protocol:`, whose default was `all`.
    """

    class Row:
        segments = "pr.dcd; pr.part1.dcd; pr.part2.dcd"
        n_segments, frames, natom = 3, 2491, 767697
        framestep_ns, trajectory_gb, simulated_ns = 0.2, 22.954, 497.23

    def test_the_default_record_lists_production_only(self):
        self.assertEqual(ndp.DEFAULTS["protocol"], "production")

    def test_a_production_stage_counts_every_segment_not_the_first(self):
        t = ndp_build._production_trajectory(self.Row())
        self.assertEqual((t["segments"], t["frames"], t["length_ns"]),
                         (3, 2491, 497.2))
        self.assertEqual(t["files"], "pr.dcd ... pr.part2.dcd")
        self.assertEqual((t["frame_interval_ps"], t["size_gb"]), (200, 22.95))

    def test_the_stage_is_named_for_the_stage_not_the_window(self):
        stage = ndp_build._whole_production(
            {"stage": "production", "name": "pr.part1", "trajectory": {"file": "x"}},
            {"1R": {"segments": 3}})
        self.assertEqual(stage["name"], "pr")
        self.assertEqual(stage["trajectory"], {"segments": 3})

    def test_an_equilibration_stage_is_left_alone(self):
        stage = {"stage": "equilibration", "name": "npt1", "trajectory": None}
        self.assertEqual(ndp_build._whole_production(stage, {"1R": {"segments": 3}}),
                         stage)


class TestGapsAreNotHidden(unittest.TestCase):
    """A hole in the middle of a run must show.

    CDR 2-JL427/3R, 2026-09-22: `pr.part6.dcd` then `pr.part10.dcd`, 108.4 ns
    apart, because parts 7-9 had not been copied from the cluster. The frames that
    ARE there add up to 370 ns, and nothing said they were not consecutive.
    """

    def test_missing_segments_are_counted_and_the_run_is_not_continuous(self):
        from core.glossary import segment_coverage
        import tempfile
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tests"))
        from test_protocol import write_dcd                  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            a, b = (os.path.join(tmp, f) for f in ("a.dcd", "b.dcd"))
            write_dcd(a, n_frames=5, istart=100000, nsavc=100000)   # 0.1-0.5 M
            write_dcd(b, n_frames=5, istart=900000, nsavc=100000)   # 0.9-1.3 M
            self.assertEqual(segment_coverage([a, b]), (0, 3, 2, 0))

    def test_a_continuous_run_has_no_gap(self):
        from core.glossary import segment_coverage
        import tempfile
        from test_protocol import write_dcd                  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            a, b = (os.path.join(tmp, f) for f in ("a.dcd", "b.dcd"))
            write_dcd(a, n_frames=5, istart=100000, nsavc=100000)
            write_dcd(b, n_frames=5, istart=600000, nsavc=100000)
            self.assertEqual(segment_coverage([a, b]), (0, 0, 1, 0))


class TestPartialFiles(unittest.TestCase):
    """A file still being copied must not be published at its intended length.

    CDR 2-JL427/3R, 2026-09-22: `pr.part8.dcd` announced 78 frames in its
    header while 3 were on disk -- an rsync from the cluster in flight. Frame counts
    now come from what the file actually holds.
    """

    def test_the_count_is_what_is_on_disk(self):
        import tempfile
        from core.dcdheader import read_dcd_header
        from test_protocol import write_dcd                  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "pr.part8.dcd")
            write_dcd(p, n_frames=78, istart=100000, nsavc=100000)
            full = read_dcd_header(p)
            self.assertEqual((full["n_frames"], full["truncated"]), (78, False))
            with open(p, "rb") as fh:
                head = fh.read(os.path.getsize(p) - full["n_frames_on_disk"] // 2
                               * ((56 if full["has_unitcell"] else 0)
                                  + 3 * (full["natom"] * 4 + 8)))
            with open(p, "wb") as fh:                 # half the frames arrive
                fh.write(head)
            cut = read_dcd_header(p)
            self.assertTrue(cut["truncated"])
            self.assertEqual(cut["n_frames_header"], 78)
            self.assertEqual(cut["n_frames"], cut["n_frames_on_disk"])
            self.assertLess(cut["n_frames"], 78)

    @staticmethod
    def _charmm_dcd(path, n_frames, natom=10, namnf=0, nset=None, titlesize=84):
        """A CHARMM DCD, fixed atoms written the way CHARMM/NAMD write them:
        the free-atom index list in the header, every atom in frame 1 and only
        the free ones after."""
        import struct
        icntrl = [n_frames if nset is None else nset, 0, 100, 100 * n_frames,
                  0, 0, 0, 0, namnf, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 24]
        rec = lambda b: struct.pack("<i", len(b)) + b + struct.pack("<i", len(b))
        out = (struct.pack("<i4s", 84, b"CORD") + struct.pack("<9i", *icntrl[:9])
               + struct.pack("<f", 0.04091) + struct.pack("<10i", *icntrl[10:])
               + struct.pack("<i", 84))
        out += (struct.pack("<i", titlesize) + struct.pack("<i", 1)
                + b" " * (titlesize - 4) + struct.pack("<i", 84))
        out += rec(struct.pack("<i", natom))
        nfree = natom - namnf
        if namnf:
            out += rec(struct.pack(f"<{nfree}i", *range(1, nfree + 1)))
        for k in range(n_frames):
            n = natom if (k == 0 or not namnf) else nfree
            out += rec(b"\0" * 48)                            # unit cell
            out += b"".join(rec(b"\0" * (4 * n)) for _ in range(3))
        with open(path, "wb") as fh:
            fh.write(out)

    def test_fixed_atoms_do_not_read_as_a_truncated_file(self):
        from core.dcdheader import read_dcd_header
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "fixed.dcd")
            self._charmm_dcd(p, n_frames=7, natom=10, namnf=4)
            h = read_dcd_header(p)
            self.assertEqual((h["n_frames"], h["n_frames_on_disk"], h["truncated"]),
                             (7, 7, False))
            self.assertEqual(h["n_fixed"], 4)

    def test_a_4d_dcd_is_sized_with_four_coordinate_blocks(self):
        import struct
        from core.dcdheader import read_dcd_header
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "four.dcd")
            self._charmm_dcd(p, n_frames=6)
            data = bytearray(open(p, "rb").read())
            data[8 + 11 * 4: 8 + 12 * 4] = struct.pack("<i", 1)   # icntrl[11]: 4D
            body_start = len(data) - 6 * (56 + 3 * (10 * 4 + 8))
            frames = b""
            for _ in range(6):                       # rewrite with a 4th block
                rec = lambda b: struct.pack("<i", len(b)) + b + struct.pack("<i", len(b))
                frames += rec(b"\0" * 48) + b"".join(rec(b"\0" * 40) for _ in range(4))
            open(p, "wb").write(bytes(data[:body_start]) + frames)
            h = read_dcd_header(p)
            self.assertEqual((h["n_frames"], h["truncated"], h["has_4d"]), (6, False, True))

    def test_an_xplor_dcd_s_timestep_is_not_read_as_a_unit_cell_flag(self):
        """In the X-PLOR form DELTA is a float64 over icntrl[9:11]; its upper
        half was read as the unit-cell flag, sizing every frame 56 bytes too
        big and undercounting."""
        import struct
        from core.dcdheader import read_dcd_header
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "xplor.dcd")
            n, natom = 9, 10
            icntrl = [n, 0, 100, 100 * n, 0, 0, 0, 0, 0]
            rec = lambda b: struct.pack("<i", len(b)) + b + struct.pack("<i", len(b))
            out = (struct.pack("<i4s", 84, b"CORD") + struct.pack("<9i", *icntrl)
                   + struct.pack("<d", 0.04091) + struct.pack("<9i", *[0] * 9)
                   + struct.pack("<i", 84))
            out += rec(struct.pack("<i", 1) + b" " * 80) + rec(struct.pack("<i", natom))
            for _ in range(n):
                out += b"".join(rec(b"\0" * (4 * natom)) for _ in range(3))
            open(p, "wb").write(out)
            h = read_dcd_header(p)
            self.assertEqual((h["n_frames"], h["truncated"], h["has_unitcell"]),
                             (9, False, False))

    def test_nset_left_at_zero_counts_the_frames_on_disk(self):
        from core.dcdheader import read_dcd_header
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "died.dcd")
            self._charmm_dcd(p, n_frames=5, nset=0)
            h = read_dcd_header(p)
            self.assertEqual((h["n_frames"], h["nset_unset"]), (5, True))

    def test_a_corrupt_title_size_is_an_error_not_a_frame_count(self):
        from core.dcdheader import read_dcd_header
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "bad.dcd")
            self._charmm_dcd(p, n_frames=3)
            with open(p, "r+b") as fh:
                fh.seek(92)
                fh.write((12345).to_bytes(4, "little"))
            self.assertIn("error", read_dcd_header(p))

    def test_a_partial_segment_is_reported(self):
        import tempfile
        from core.glossary import segment_coverage
        from test_protocol import write_dcd                  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            a, b = (os.path.join(tmp, f) for f in ("a.dcd", "b.dcd"))
            write_dcd(a, n_frames=5, istart=100000, nsavc=100000)
            write_dcd(b, n_frames=5, istart=600000, nsavc=100000)
            with open(b, "rb") as fh:
                keep = fh.read(os.path.getsize(b) - 1000)
            with open(b, "wb") as fh:
                fh.write(keep)
            self.assertEqual(segment_coverage([a, b])[3], 1)


class TestCatalogue(unittest.TestCase):
    def test_the_identifier_reads_like_the_lab_catalogue(self):
        first = catalog.suggest_identifier("hiv", CDR, "CDR/d2830a2a", what="1-1F4")
        self.assertEqual(first, "TESTLAB-HIV-CDR-1F4-001")
        # A second system of the same folder and name gets the next number.
        again = catalog.suggest_identifier("hiv", CDR, "CDR/other",
                                           taken=[first], what="1-1F4")
        self.assertEqual(again, "TESTLAB-HIV-CDR-1F4-002")

    def test_the_topic_comes_from_the_folder_and_is_one_the_site_uses(self):
        self.assertEqual(catalog.topic_of(CDR, {"topics": ["HIV-1"]}), "hiv")
        self.assertEqual(catalog.topic_of("/data/SPIKE_open", {}), "covid-19")
        self.assertIsNone(catalog.topic_of("/data/something", {}))

    def test_catalog_json_carries_what_the_page_shows(self):
        record = {"identifier": "TESTLAB-HIV-CDR-1F4-001", "topic": "hiv",
                  "tags": ["1F4"], "title": "Env + 1F4", "description": "...",
                  "size_gb": 64.36, "files": 29, "record": "CDR_1-1F4.json",
                  "license": "CC-BY-4.0", "data_availability": "on request",
                  "sampling": {"replicas": 3, "complete_replicas": 2,
                               "total_ns": 1395.0}}
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(layout.results(tmp), "catalog.json")
            catalog.write(path, tmp, [record])
            doc = json.loads(open(path).read())
            row = doc["datasets"][0]
            self.assertEqual(row["identifier"], "TESTLAB-HIV-CDR-1F4-001")
            self.assertEqual((row["size_gb"], row["files"]), (64.36, 29))
            self.assertEqual((row["replicas"], row["complete_replicas"]), (3, 2))
            self.assertEqual(row["record"], "CDR_1-1F4.json")
            # A run that changed nothing must not rewrite the file, or every
            # run would look like a change to whoever syncs it to the site.
            before = os.stat(path).st_mtime_ns
            catalog.write(path, tmp, [record])
            self.assertEqual(os.stat(path).st_mtime_ns, before)




class TestCollidingNames(unittest.TestCase):
    """Three or more records of one folder wanting the same file name.

    Two was handled. Three was not: the de-collision POPPED the bare name, so
    the third record found it free and took it, while the record that had been
    renamed went on carrying the old name inside its own JSON -- and
    `catalog.entry()` reads the name out of the record, not out of the index
    entry. The lab-wide catalogue ended up with one stem used four times and
    four dead pointers.
    """

    def claim(self, n):
        """n records all wanting `X.json`; returns [(name, incumbent), ...]."""
        claimed, out = {}, []
        for i in range(n):
            out.append(ndp_build._claim_name("X.json", f"d/gid{i}", claimed))
        return out, claimed

    def test_three_records_get_three_distinct_names(self):
        out, _ = self.claim(3)
        names = [name for name, _ in out]
        self.assertEqual(len(set(names)), 3, names)

    def test_nobody_keeps_the_bare_name_once_it_is_contested(self):
        """Not even the first. An ambiguous name held by one of a colliding
        group is the thing that made the catalogue unreadable."""
        out, claimed = self.claim(3)
        self.assertNotIn("X.json", [name for name, _ in out[1:]])
        self.assertIsNone(claimed["X.json"])

    def test_the_incumbent_is_reported_so_it_can_be_renamed(self):
        out, _ = self.claim(2)
        self.assertIsNone(out[0][1])
        self.assertEqual(out[1][1], "d/gid0")

    def test_a_fourth_still_gets_its_own_name(self):
        out, _ = self.claim(4)
        self.assertEqual(len(set(name for name, _ in out)), 4)

    def test_an_uncontested_name_is_left_exactly_as_it_was(self):
        claimed = {}
        name, incumbent = ndp_build._claim_name("CDR_1-1F4_3x500ns.json",
                                                "CDR/d2830a2a", claimed)
        self.assertEqual(name, "CDR_1-1F4_3x500ns.json")
        self.assertIsNone(incumbent)


class TestIdentifierDedupe(unittest.TestCase):
    """The `-002` bump never ran.

    `suggest_identifier(..., taken=...)` was fed the index entries, which have
    never carried an `identifier` key, so `taken` was a list of None and every
    record of a folder got `-001`. Eighteen identifiers are duplicated in the
    lab-wide catalogue because of it.
    """

    def test_a_taken_identifier_is_bumped(self):
        first = catalog.suggest_identifier("hiv", CDR, "CDR/a", what="1-1F4")
        second = catalog.suggest_identifier("hiv", CDR, "CDR/b", what="1-1F4",
                                            taken=[first])
        third = catalog.suggest_identifier("hiv", CDR, "CDR/c", what="1-1F4",
                                           taken=[first, second])
        self.assertEqual(len({first, second, third}), 3)
        self.assertTrue(second.endswith("-002"))
        self.assertTrue(third.endswith("-003"))

    def test_a_list_of_none_bumps_nothing(self):
        """What the old call actually passed. Pinned so the regression is
        visible as a difference rather than as an absence."""
        a = catalog.suggest_identifier("hiv", CDR, "CDR/a", what="1-1F4",
                                       taken=[None, None])
        b = catalog.suggest_identifier("hiv", CDR, "CDR/b", what="1-1F4",
                                       taken=[None, None])
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main(verbosity=2)
