#!/usr/bin/env python3
"""
test_ndp.py -- the NDP records say what the glossary says, and nothing more.

    python3 tests/test_ndp.py          # or: python3 -m pytest tests/

Fixtures are synthetic and self-contained: nothing reads production data, and
nothing writes outside a temp directory. Each test pins one way a record could
end up claiming something nobody established.
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

import pandas as pd                              # noqa: E402
import yaml                                      # noqa: E402

import ndp_build                                 # noqa: E402
from core import ndp                         # noqa: E402
from core import layout, settings                 # noqa: E402

TEMPLATE = {
    # The catalogue's own fields come first, then the NDP template in order.
    None: ["identifier", "topic", "tags", "title", "description", "size_gb",
           "files", "lab", "contact", "system", "setup", "sampling",
           "data_availability", "related_publication", "license"],
    "system": ["construct", "starting_structure", "glycan_profile", "total_atoms"],
    "setup": ["force_field", "water_model", "ions", "software", "ensemble",
              "temperature_K", "pressure_bar", "timestep_fs"],
    "sampling": ["replicas", "complete_replicas", "planned_length_ns",
                 "length_per_replica_ns", "frame_interval_ps"],
}


def md(dataset, md_id, group, ns=500.0, **over):
    row = {"md_key": f"{dataset}/{md_id}", "dataset_id": dataset, "md_id": md_id,
           "state": "ok", "issues": "", "topology_group": group,
           "program": "NAMD", "version": "3.0.2", "natom": 743875.0,
           "simulated_ns": ns, "framestep_ns": 0.1, "timestep_fs": 2.0,
           "temp_K": 310.0, "ensemble": "NPT", "ff": "CHARMM36m", "wat": None,
           "config_path": f"{md_id}/PR.conf", "source_root": "/nowhere",
           "source_dir": md_id,
           # what the catalogue counts: frames, volume and files per replica
           "frames": 5000.0, "trajectory_gb": 10.0, "segments": "pr.dcd",
           "target_ns": None, "nominal_ns": None}
    row.update(over)
    return row


def read(path):
    with open(path) as fh:
        return fh.read()


def frame(*rows):
    return pd.DataFrame(list(rows))


def no_pressure(*_):
    return None


class TestNominalLength(unittest.TestCase):
    def test_the_largest_declared_length_reached(self):
        from core.glossary import nominal_of
        self.assertEqual(nominal_of(1008.0), 1000)
        self.assertEqual(nominal_of(500.0), 500)       # 500 of a planned 504
        self.assertEqual(nominal_of(585.1), 500)       # ran longer: still 500
        self.assertEqual(nominal_of(1500.0), 1500)
        self.assertIsNone(nominal_of(489.0))           # beyond the 2 % tolerance
        self.assertIsNone(nominal_of(72.2))
        self.assertIsNone(nominal_of(None))
        self.assertEqual(nominal_of(100.8, [100]), 100)


class TestGrouping(unittest.TestCase):
    def test_replicas_of_one_topology_are_one_system(self):
        g = frame(md("open", "A_1R", "aaaa"), md("open", "A_2R", "aaaa"),
                  md("open", "B_1R", "bbbb"))
        systems = ndp.systems_of(g)
        self.assertEqual([rid for rid, _, _ in systems], ["open/aaaa", "open/bbbb"])
        self.assertEqual(len(systems[0][1]), 2)

    def test_same_topology_in_two_datasets_is_two_records(self):
        g = frame(md("models", "A", "aaaa"), md("aptasuite", "A", "aaaa"))
        self.assertEqual(len(ndp.systems_of(g)), 2)

    def test_blocked_mds_are_not_sampling(self):
        g = frame(md("open", "A_1R", "aaaa"),
                  md("open", "A_2R", "aaaa", state="blocked"))
        (_, rows, n_blocked), = ndp.systems_of(g)
        self.assertEqual((len(rows), n_blocked), (1, 1))
        (_, rows, _), = ndp.systems_of(g, include_blocked=True)
        self.assertEqual(len(rows), 2)


class TestMeasure(unittest.TestCase):
    def test_disagreement_is_null_and_reported_never_majority(self):
        g = frame(md("d", "1R", "a"), md("d", "2R", "a"),
                  md("d", "3R", "a", temp_K=300.0))
        m, warn = ndp.measure(g, no_pressure)
        self.assertIsNone(m["temperature_K"])
        self.assertTrue(any("temperature_K" in w for w in warn))

    def test_an_md_that_does_not_say_is_not_a_disagreement(self):
        g = frame(md("d", "1R", "a"), md("d", "2R", "a", ff=None, version=None))
        m, warn = ndp.measure(g, no_pressure)
        self.assertEqual(m["force_field"], "CHARMM36m")
        self.assertEqual(m["software"], "NAMD 3.0.2")
        self.assertEqual(warn, [])

    def test_sampling(self):
        g = frame(md("d", "1R", "a", ns=499.9), md("d", "2R", "a", ns=500.1),
                  md("d", "3R", "a", ns=499.9))
        m, _ = ndp.measure(g, no_pressure)
        self.assertEqual((m["replicas"], m["total_ns"], m["frame_interval_ps"]),
                         (3, 1499.9, 100))
        self.assertEqual(ndp.sampling_phrase(m), "3 × 500 ns")

    def test_uneven_replicas_are_not_called_n_times_l(self):
        g = frame(md("d", "1R", "a", ns=505.0), md("d", "2R", "a", ns=258.6))
        m, warn = ndp.measure(g, no_pressure)
        self.assertEqual(ndp.sampling_phrase(m), "2 replicas, 764 ns in total")
        g = frame(md("d", "1R", "a", ns=547.3), md("d", "pr", "a", ns=1.5))
        _, warn = ndp.measure(g, no_pressure)
        self.assertTrue(any("uneven" in w for w in warn))

    def test_pressure_is_read_from_the_config_not_assumed(self):
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, "PR.conf"), "w") as fh:
                fh.write("if {0} {\nlangevinPistonTarget 5.0\n}\n"
                         "langevinPiston on\nlangevinPistonTarget 1.01325\n")
            with open(os.path.join(tmp, "md.in"), "w") as fh:
                fh.write("&cntrl\n ntp=1, taup=2.0,\n/\n")
            self.assertEqual(ndp.pressure_bar(tmp, "PR.conf", "NAMD"), 1.01325)
            self.assertIsNone(ndp.pressure_bar(tmp, "md.in", "AMBER"))
            self.assertIsNone(ndp.pressure_bar(tmp, "absent.conf", "NAMD"))


class TestRecord(unittest.TestCase):
    def setUp(self):
        self.rows = frame(md("d", "1R", "a"), md("d", "2R", "a"))
        self.m, _ = ndp.measure(self.rows, lambda *_: 1.01325)
        self.spec = {"contact": "x@example.edu", "pdb_ids": ["5FYL"], "wat": "TIP3P",
                     "license": "Creative Commons Attribution 4.0 International"}

    def test_every_template_field_is_present_even_when_unknown(self):
        rec = ndp.build_record("d/a", self.rows, self.m, {}, {}, ndp.DEFAULTS)
        for block, keys in TEMPLATE.items():
            have = list(rec if block is None else rec[block])
            self.assertEqual(have[:len(keys)], keys)
        self.assertIsNone(rec["title"])
        self.assertIn("title", ndp.missing_fields(rec))

    def test_declared_values_fill_what_ndp_yaml_leaves_empty(self):
        rec = ndp.build_record("d/a", self.rows, self.m, {}, self.spec, ndp.DEFAULTS)
        self.assertEqual(rec["system"]["starting_structure"], "PDB 5FYL")
        self.assertEqual(rec["setup"]["water_model"], "TIP3P")
        self.assertEqual(rec["license"], "CC-BY-4.0")
        self.assertEqual(rec["setup"]["pressure_bar"], 1.01325)

    def test_a_curated_value_never_replaces_a_measured_one(self):
        warn = []
        rec = ndp.build_record("d/a", self.rows, self.m, {"temperature_K": 300},
                               self.spec, ndp.DEFAULTS, warn)
        self.assertEqual(rec["setup"]["temperature_K"], 310)
        self.assertTrue(any("ignored" in w for w in warn))

    def test_ready_only_when_nothing_is_missing(self):
        cur = {"title": "BG505 SOSIP.664 glycosylated trimer", "construct": "BG505",
               "glycan_profile": "none", "ions": "150 mM NaCl",
               "related_publication": "in preparation",
               # Who to cite. A record with no authors is not ready to
               # register, so this belongs in the set that makes one ready.
               "authors": ["A Researcher", "B Colleague"],
               # the lab catalogue needs a home and a stable name for it
               "identifier": "TESTLAB-HIV-TEST-001", "topic": "hiv",
               # the federation's questions: lab.yaml leaves both unanswered
               "access_level": "private", "ai_training": "prohibited"}
        rec = ndp.build_record("d/a", self.rows, self.m, cur, self.spec, ndp.DEFAULTS)
        self.assertEqual(ndp.missing_fields(rec), [])
        self.assertEqual(rec["authors"], "A Researcher, B Colleague")
        self.assertEqual(ndp.catalog_title(rec, self.m),
                         "BG505 SOSIP.664 glycosylated trimer — all-atom MD, 2 × 500 ns")
        self.assertIn("metadata only", ndp.catalog_description(rec, self.m))

    def test_authors_come_from_the_dataset_as_one_line(self):
        """A list in settings.yaml reaches the record as the line a page shows.

        The analysis page prints `doc["authors"]` straight into a table cell.
        While the record carried no such key at all it showed "to be added
        before registration" for a value the dataset block had held since the
        folder was first indexed -- a field can only be missing if something
        actually looked for it.
        """
        spec = dict(self.spec)
        spec["authors"] = ["Ada Researcher", "Ben Colleague"]
        rec = ndp.build_record("d/a", self.rows, self.m, {}, spec, ndp.DEFAULTS)
        self.assertEqual(rec["authors"], "Ada Researcher, Ben Colleague")
        self.assertNotIn("authors", ndp.missing_fields(rec))
        # An empty list is not an answer.
        spec["authors"] = []
        rec = ndp.build_record("d/a", self.rows, self.m, {}, spec, ndp.DEFAULTS)
        self.assertIsNone(rec["authors"])
        self.assertIn("authors", ndp.missing_fields(rec))

    def test_a_dataset_shares_values_and_a_system_overrides_them(self):
        doc = {"datasets": {"d": {"ions": "150 mM NaCl", "keywords": ["HIV-1 Env"]}},
               "systems": [{"id": "d/a", "ions": None, "keywords": ["gp120"]},
                           {"id": "d/b", "ions": "neutralising only"}]}
        self.assertEqual(ndp.curation_for("d/a", doc)["ions"], "150 mM NaCl")
        self.assertEqual(ndp.curation_for("d/a", doc)["keywords"], ["HIV-1 Env", "gp120"])
        self.assertEqual(ndp.curation_for("d/b", doc)["ions"], "neutralising only")
        self.assertEqual(ndp.curation_for("other/a", doc)["keywords"], [])

    def test_the_protocol_keyword_keeps_only_what_is_asked(self):
        stages = [{"name": "min1", "stage": "minimisation"},
                  {"name": "npt2", "stage": "equilibration"},
                  {"name": "pr", "stage": "production"},
                  {"name": "pr_r1", "stage": "production"}]
        def names(cur, defaults=ndp.DEFAULTS, warn=None):
            rec = ndp.build_record("d/a", self.rows, self.m, cur, {}, defaults,
                                   warn, stages)
            return [s["name"] for s in rec["protocol"]]
        # The default publishes the production phase and nothing else: the
        # owner does not want min/heating/NPT in what goes to the NDP.
        self.assertEqual(names({}), ["pr", "pr_r1"])
        self.assertEqual(names({"protocol": "production"}), ["pr", "pr_r1"])
        self.assertEqual(names({"protocol": "all"}),
                         ["min1", "npt2", "pr", "pr_r1"])
        self.assertEqual(names({}, {**ndp.DEFAULTS, "protocol": "none"}), [])
        warn = []
        # An unknown value falls back to listing everything, and says so.
        self.assertEqual(len(names({"protocol": "prod"}, warn=warn)), 4)
        self.assertTrue(warn)

    def test_nvt_is_not_missing_a_pressure(self):
        rows = frame(md("d", "1R", "a", ensemble="NVT"))
        m, _ = ndp.measure(rows, no_pressure)
        rec = ndp.build_record("d/a", rows, m, {}, {}, ndp.DEFAULTS)
        self.assertNotIn("setup.pressure_bar", ndp.missing_fields(rec))

    def test_keywords_ndp_would_refuse_are_caught(self):
        good, bad = ndp.keywords_of({"keywords": ["HIV-1 Env", "gp120/gp41", "test lab"]},
                                    ndp.DEFAULTS)
        self.assertEqual(good, ["molecular dynamics", "Test Lab", "HIV-1 Env"])
        self.assertEqual(bad, ["gp120/gp41"])


class TestBuild(unittest.TestCase):
    """End to end, on a folder holding nothing but a glossary."""

    def build(self, tmp, *rows):
        reg = layout.mddb_registry(tmp)
        os.makedirs(reg, exist_ok=True)
        frame(*rows).to_parquet(os.path.join(reg, "mds.parquet"), index=False)
        if not settings.exists(tmp):
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
        return ndp_build.main([tmp])

    def test_the_stale_sweep_does_not_touch_another_tool_s_files(self):
        """A full build must not delete what it did not write.

        The sweep walked the WHOLE of MD_RESULTS removing every `.json` it had
        not produced. Records sit at the top level, but `analyses/` sits under
        it -- so on 2026-09-24 a full build was one run away from deleting 39
        files belonging to the analysis kit, including each system's
        `frames.json`, which is the time axis every plot's frame provenance is
        read from. Nothing there is regenerable from this pipeline.
        """
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa")), 0)
            results = layout.results(tmp)
            bundle = os.path.join(results, "analyses", "SOMETHING", "data")
            os.makedirs(bundle, exist_ok=True)
            keep = [os.path.join(results, "analyses", "SOMETHING",
                                 "frames.json"),
                    os.path.join(bundle, "stage.manifest.json"),
                    os.path.join(results, "notes.json")]
            for f in keep:
                with open(f, "w") as fh:
                    fh.write('{"not": "a record"}')

            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa")), 0)
            for f in keep:
                self.assertTrue(os.path.exists(f),
                                f"a full build deleted {os.path.relpath(f, tmp)}")

    def test_a_title_written_later_does_not_rename_the_record(self):
        """Filling in `title:` must leave the record where it is.

        Every path anyone holds -- the analysis bundle beside it, the
        production files, a link already sent -- is keyed on the file name.
        """
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa")), 0)
            results = layout.results(tmp)
            before = sorted(f for f in os.listdir(results)
                            if f.endswith(".json") and f != "catalog.json")
            cur_p = settings.path(tmp)
            # Read it out BEFORE opening for write: `open(..., "w")` truncates,
            # so reading inside the block reads the empty file back.
            text = read(cur_p).replace("    title:\n",
                                       "    title: HIV - Env Open State\n")
            with open(cur_p, "w") as fh:
                fh.write(text)

            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa")), 0)
            after = sorted(f for f in os.listdir(results)
                           if f.endswith(".json") and f != "catalog.json")
            self.assertEqual(after, before)
            got = json.loads(read(os.path.join(results, after[0])))
            self.assertEqual(got["title"], "HIV - Env Open State")

    def test_rerun_keeps_edits_appends_new_systems_and_changes_no_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            cur_p = settings.path(tmp)
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa")), 0)
            text = read(cur_p).replace("    title:\n", "    title: My trimer  # mine\n")
            with open(cur_p, "w") as fh:
                fh.write(text)

            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa"),
                                        md("d", "X_1R", "bbbb")), 0)
            after = read(cur_p)
            self.assertTrue(after.startswith(text),
                            "text a person wrote in settings.yaml was rewritten")
            ids = [s["id"] for s in yaml.safe_load(after)["ndp"]["systems"]]
            self.assertEqual(ids, ["d/aaaa", "d/bbbb"])
            rec_p = os.path.join(layout.results(tmp), self._name(tmp, "d/aaaa"))
            self.assertEqual(json.loads(read(rec_p))["title"], "My trimer")

            stamp = os.stat(rec_p).st_mtime_ns
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa"),
                                        md("d", "X_1R", "bbbb")), 0)
            self.assertEqual(os.stat(rec_p).st_mtime_ns, stamp)

    def _name(self, tmp, record_id):
        """The file a record was written to, from the index this run wrote."""
        index = json.loads(read(layout.ndp_index(tmp)))
        return next(r["file"] for r in index["records"]
                    if r["record_id"] == record_id)

    def test_a_record_whose_system_is_gone_does_not_stay(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build(tmp, md("d", "1R", "aaaa"), md("d", "X_1R", "bbbb"))
            gone = os.path.join(layout.results(tmp), self._name(tmp, "d/bbbb"))
            self.assertTrue(os.path.exists(gone))
            self.build(tmp, md("d", "1R", "aaaa"))
            self.assertFalse(os.path.exists(gone))
            index = json.loads(read(layout.ndp_index(tmp)))
            self.assertEqual([r["record_id"] for r in index["records"]], ["d/aaaa"])

    def test_only_mds_of_a_production_length_are_replicas(self):
        """3 x ~1 us replicas and a 72 ns run that stopped.

        `lengths: all` is the default -- the catalogue lists what exists, each
        replica marked -- so this asks for `nominal` explicitly first.
        """
        from core.glossary import nominal_of
        with tempfile.TemporaryDirectory() as tmp:
            rows = [md("d", f"{i}R", "aaaa", ns=ns, nominal_ns=nominal_of(ns))
                    for i, ns in ((1, 1008.0), (2, 1008.0), (3, 1008.0), (4, 72.2))]
            self.assertEqual(self.build(tmp, *rows), 0)
            text = read(settings.path(tmp))
            with open(settings.path(tmp), "w") as fh:
                fh.write(text.replace("lengths: all", "lengths: nominal"))
            self.assertEqual(self.build(tmp, *rows), 0)
            rec_p = os.path.join(layout.results(tmp), self._name(tmp, "d/aaaa"))
            rec = json.loads(read(rec_p))
            self.assertEqual(rec["sampling"]["replicas"], 3)
            self.assertEqual(rec["sampling"]["nominal_length_ns"], 1000)
            index = json.loads(read(layout.ndp_index(tmp)))
            self.assertEqual(index["records"][0]["left_out_short"], ["d/4R (72.2 ns)"])
            text = read(settings.path(tmp))
            with open(settings.path(tmp), "w") as fh:
                fh.write(text.replace("lengths: nominal", "lengths: all"))
            self.assertEqual(self.build(tmp, *rows), 0)
            rec = json.loads(read(os.path.join(layout.results(tmp),
                                               self._name(tmp, "d/aaaa"))))
            self.assertEqual(rec["sampling"]["replicas"], 4)

    def test_check_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = layout.mddb_registry(tmp)
            os.makedirs(reg)
            frame(md("d", "1R", "aaaa")).to_parquet(os.path.join(reg, "mds.parquet"),
                                                    index=False)
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
            before = read(settings.path(tmp))
            self.assertEqual(ndp_build.main([tmp, "--check"]), 1)
            self.assertEqual(sorted(os.listdir(layout.mddb_dir(tmp))), ["registry"])
            self.assertEqual(read(settings.path(tmp)), before)
            self.assertEqual(os.listdir(layout.results(tmp)), ["settings.yaml"])


class TestGlossaryInputs(unittest.TestCase):
    """ndp_build.py says so when an edit needs the glossary rebuilt."""
    DOC = {"folder": {"root": "/x"}, "datasets": [{"id": "OPEN"}],
           "ndp": {"systems": [{"id": "OPEN/abc", "title": "t"}]}}

    def test_digest_moves_only_with_the_glossary_inputs(self):
        import copy
        d0 = ndp.glossary_inputs(self.DOC)
        other = copy.deepcopy(self.DOC)
        other["ndp"]["systems"][0]["title"] = "a better title"
        self.assertEqual(ndp.glossary_inputs(other), d0, "title reaches ndp_build directly")
        for change in (lambda d: d["folder"].update(stages={"derived": "x"}),
                       lambda d: d["folder"].update(replicas=["_(\\d+)$"]),
                       lambda d: d.update(production_lengths_ns=[400]),
                       lambda d: d["ndp"]["systems"][0].update(frame_interval_ps=100)):
            d = copy.deepcopy(self.DOC)
            change(d)
            self.assertNotEqual(ndp.glossary_inputs(d), d0)

    def test_stale_only_when_a_stamp_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(ndp_build._glossary_stale(tmp, self.DOC), "no stamp: silent")
            reg = layout.mddb_registry(tmp)
            os.makedirs(reg)
            with open(os.path.join(reg, "glossary_inputs.json"), "w") as fh:
                json.dump({"digest": ndp.glossary_inputs(self.DOC)}, fh)
            self.assertFalse(ndp_build._glossary_stale(tmp, self.DOC))
            edited = {**self.DOC, "production_lengths_ns": [400]}
            self.assertTrue(ndp_build._glossary_stale(tmp, edited))


class TestDeclaredNotes(unittest.TestCase):
    """Issues 4-7 of SPIKE_CLEAVED (2026-09-30): said, never silently changed."""
    OPEN = [1000.5, 1000.3, 1006.3, 404.3, 404.3, 406.5]
    WARN = ["`nominal_length_ns`: the MDs disagree (1000, 400); left null. Split "
            "the system or correct the MD that is wrong.",
            "replica lengths are very uneven (404.3 to 1006.3 ns): check ...",
            "something else"]

    def rows(self, lengths, issues=""):
        return pd.DataFrame({"simulated_ns": lengths, "issues": [issues] * len(lengths)})

    def test_declared_mixed_lengths_answer_the_two_length_warnings(self):
        warn = list(self.WARN)
        notes = ndp.declared_notes(self.rows(self.OPEN), {"nominal_lengths_ns": [1000, 400]}, warn)
        self.assertEqual(warn, ["something else"])
        self.assertIn("3 x 1000 ns + 3 x 400 ns", notes[0])

    def test_undeclared_or_wrong_lengths_change_nothing(self):
        warn = list(self.WARN)
        self.assertEqual(ndp.declared_notes(self.rows(self.OPEN), {}, warn), [])
        self.assertEqual(warn, self.WARN)
        notes = ndp.declared_notes(self.rows(self.OPEN), {"nominal_lengths_ns": [1000]}, warn)
        self.assertEqual(warn, self.WARN, "400 ns replicas reach no declared length")
        self.assertIn("does not match", notes[0])
        notes = ndp.declared_notes(self.rows(self.OPEN), {"nominal_lengths_ns": [1000, 400, 200]}, warn)
        self.assertEqual(warn, self.WARN, "a declared length nobody reaches")

    def test_config_declared_is_a_note_only_when_all_three_are_stated(self):
        rows = self.rows([1000.0, 1000.0], issues="no_config")
        full = {"ensemble": "NPT", "temperature_K": 310, "timestep_fs": 2}
        self.assertIn("config_declared", ndp.declared_notes(rows, full, [])[0])
        self.assertEqual(ndp.declared_notes(rows, {"ensemble": "NPT"}, []), [])
        self.assertEqual(ndp.declared_notes(self.rows([1000.0]), full, []), [])

    def test_a_slugged_name_stays_when_the_sampling_moves(self):
        m6 = {"replicas": 6, "replica_lengths_ns": self.OPEN,
              "length_per_replica_ns": 703.7, "total_ns": 4222.2}
        cur = {"slug": "spike-fulllength-open"}
        old = "spike-fulllength-open_6x1us.json"
        self.assertEqual(ndp.file_name("OPEN/a", cur, measured=m6, existing=old), old)
        self.assertTrue(ndp.slug_would_rename("OPEN/a", cur, m6, old))
        self.assertNotEqual(ndp.file_name("OPEN/a", {**cur, "rename": True},
                                          measured=m6, existing=old), old)
        self.assertNotEqual(ndp.file_name("OPEN/a", {"slug": "new-name"},
                                          measured=m6, existing=old), old,
                            "a new slug is a rename someone asked for")
        fresh = ndp.file_name("OPEN/a", cur, measured=m6)
        self.assertEqual(ndp.file_name("OPEN/a", cur, measured=m6, existing=fresh), fresh)
        self.assertFalse(ndp.slug_would_rename("OPEN/a", cur, m6, fresh))

    def test_a_swallowed_comment_is_named(self):
        doc = {"datasets": [{"contact": "pi@example.edu# an e-mail that will work"}],
               "ndp": {"systems": [{"id": "a/b", "title": "Spike #3 open, run 2"}]}}
        got = ndp.swallowed_comments(doc)
        self.assertEqual(len(got), 1)
        self.assertIn("pi@example.edu#", got[0])
        self.assertEqual(ndp.swallowed_comments({"datasets": [{"contact": "a@b.edu"}]}), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
