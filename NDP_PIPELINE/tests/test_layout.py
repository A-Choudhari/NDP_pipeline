#!/usr/bin/env python3
"""
test_layout.py -- one visible folder, one hidden one, and a move that loses
nothing.

    python3 tests/test_layout.py        # or: python3 -m pytest tests/

Decided 2026-09-22, after a run left 47 files in 24 directories across two
visible trees and asked a person to edit three separate YAML files. These pin
what a person now sees, and that a folder from the old layout arrives at it
with every value they wrote intact.
"""

import contextlib
import io
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
# The lab the tests run as, whatever the kit's own lab.yaml says.
os.environ["MD_PIPELINE_LAB"] = os.path.join(HERE, "lab_fixture.yaml")
os.environ.pop("NDP_ID_PREFIX", None)

import yaml                                          # noqa: E402

from core import config, layout, migrate, settings, yamledit   # noqa: E402
import init                                                  # noqa: E402


class TestWhereThingsGo(unittest.TestCase):
    def test_every_output_directory_answers_the_same_folder(self):
        for given in ("/f", "/f/MD_RESULTS", "/f/.md_pipeline",
                      "/f/.md_pipeline/mddb", "/f/.md_pipeline/index",
                      "/f/SUBSTRATE", "/f/MDDB_substrate", "/f/"):
            self.assertEqual(layout.base(given), "/f", given)

    def test_what_a_person_sees_is_five_names(self):
        visible = {os.path.basename(p) for p in (
            layout.settings("/f"), layout.readme("/f"), layout.glossary_md("/f"),
            layout.registration_md("/f"))}
        self.assertEqual(visible, {"settings.yaml", "README.md", "glossary.md",
                                   "NDP_registration.md"})
        for p in (layout.registry("/f"), layout.tables("/f"), layout.legend("/f"),
                  layout.notes("/f"), layout.mddb_registry("/f"),
                  layout.exports("/f"), layout.cache("/f"), layout.logs("/f"),
                  layout.status("/f"), layout.backups("/f"),
                  layout.ndp_index("/f")):
            self.assertTrue(p.startswith("/f/.md_pipeline/"), p)

    def test_records_move_out_of_the_way_once_there_are_many(self):
        self.assertEqual(layout.records_dir("/f", 1), layout.results("/f"))
        self.assertEqual(layout.records_dir("/f", layout.RECORDS_INLINE),
                         layout.results("/f"))
        self.assertEqual(layout.records_dir("/f", layout.RECORDS_INLINE + 1),
                         os.path.join(layout.results("/f"), "NDP_records"))

    def test_our_own_output_is_never_indexed(self):
        """Including the old layout's: a folder holding both must not have
        either walked into, or the index would describe itself."""
        with tempfile.TemporaryDirectory() as tmp:
            for d in (layout.RESULTS, layout.INTERNAL, layout.LEGACY_INDEX,
                      layout.LEGACY_MDDB, "DATA", "stray"):
                os.makedirs(os.path.join(tmp, d))
            cfg = config.validate({"owner": "t", "root": tmp,
                                   "corpora": [{"name": "DATA", "path": "DATA"}]},
                                  settings.path(tmp))
            status = {n: st for n, st, _ in config.coverage(cfg)}
            for d in (layout.RESULTS, layout.LEGACY_INDEX, layout.LEGACY_MDDB):
                self.assertEqual(status.get(d), "output", d)
            self.assertEqual(status.get("stray"), "UNCOVERED")
            walked = {os.path.basename(d) for d, _, _ in config.walk(
                tmp, cfg["prune"], exclude=cfg["exclude_paths"], max_depth=1)}
            self.assertNotIn(layout.RESULTS, walked)
            self.assertNotIn(layout.INTERNAL, walked)


class TestSettingsFile(unittest.TestCase):
    def test_one_file_with_every_section(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section=f"folder:\n  owner: t\n  root: {tmp}\n")
            doc = settings.load(tmp)
            self.assertEqual(sorted(doc), ["datasets", "folder",
                                           "ndp", "production_lengths_ns"])
            self.assertEqual(settings.config_of(tmp)["root"], tmp)
            self.assertEqual(settings.production_lengths(tmp), [500, 1000, 1500])

    def test_a_new_file_keeps_each_section_s_explanation(self):
        """The explanations were imported from modules that no longer existed,
        behind a bare except, so every new file got bare `datasets:`/`ndp:`."""
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
            text = settings.read(tmp)
            self.assertIn("the science behind each dataset", text)
            self.assertIn("the catalogue record", text)

    def test_a_dataset_of_many_runs_keeps_every_path(self):
        """`paths:` is the dataset's selector; it was cut at 200 entries."""
        import pandas as pd
        import datasets_init
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
            reg = os.path.join(layout.index_dir(tmp), "registry")
            os.makedirs(reg)
            n = 250
            pd.DataFrame({
                "rel_path": [f"big/run{i:03d}" for i in range(n)],
                "system": ["unknown"] * n, "condition": ["unknown"] * n,
                "n_production": [1] * n, "replica": ["unknown"] * n,
                "traj_gb": [1.0] * n, "prod_ns_header": [500.0] * n,
            }).to_parquet(os.path.join(reg, "runs.parquet"))
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(datasets_init.main([tmp, "--yes"]),
                                 datasets_init.EXIT_OK)
            (ds,) = settings.datasets_of(tmp)
            self.assertEqual(len(ds["paths"]), n)

    def test_a_written_value_is_never_overwritten_and_comments_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n",
                            datasets_section="\ndatasets:\n  - id: d\n    wat:  # mine\n"
                                             "    citation: my paper\n")
            settings.fill(tmp, {"d": {"wat": "TIP3P"}})
            text = settings.read(tmp)
            self.assertIn("# mine", text)
            self.assertEqual(settings.datasets_of(tmp)[0]["wat"], "TIP3P")
            with self.assertRaises(ValueError):
                settings.fill(tmp, {"d": {"citation": "someone else's paper"}})
            self.assertEqual(settings.datasets_of(tmp)[0]["citation"], "my paper")

    def test_a_new_block_joins_its_own_section(self):
        text = "folder:\n  a: 1\n\ndatasets:\n  - id: x\n\nndp:\n  q: 2\n"
        out = yamledit.insert_under(text, "datasets", "  - id: y\n")
        doc = yaml.safe_load(out)
        self.assertEqual([d["id"] for d in doc["datasets"]], ["x", "y"])
        self.assertEqual(doc["ndp"], {"q": 2})

    def test_every_write_leaves_a_copy_of_what_was_there(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section="folder:\n  root: /x\n")
            settings.write(tmp, "folder:\n  root: /y\n")
            kept = os.listdir(layout.backups(tmp))
            self.assertEqual(len(kept), 1, kept)
            self.assertIn("root: /x", open(os.path.join(layout.backups(tmp),
                                                        kept[0])).read())


class TestMigration(unittest.TestCase):
    def old_folder(self, tmp):
        """A folder as the pipeline left it before 2026-09-22."""
        old_index, old_mddb = layout.legacy(tmp)
        os.makedirs(os.path.join(old_index, "registry"))
        os.makedirs(os.path.join(old_mddb, "work", "ds", "1R"))
        os.makedirs(os.path.join(old_mddb, "reports"))
        with open(os.path.join(old_index, "substrate.yaml"), "w") as fh:
            fh.write(f"# mine\nowner: labmate\nroot: {tmp}\n"
                     "systems:\n  gp120: 'OPEN'   # my label\n")
        with open(os.path.join(old_index, "curated.yaml"), "w") as fh:
            fh.write("facts: []\n")
        with open(os.path.join(old_mddb, "datasets.yaml"), "w") as fh:
            fh.write("# generated\ndatasets:\n  - id: gp120\n"
                     "    citation: my paper   # do not touch\n")
        with open(os.path.join(old_mddb, "ndp.yaml"), "w") as fh:
            fh.write("# ndp\ndefaults:\n  lab: My Lab\nsystems:\n"
                     "  - id: gp120/aaaa\n    title: My trimer\n")
        with open(os.path.join(old_index, "registry", "runs.parquet"), "wb") as fh:
            fh.write(b"stub")
        with open(os.path.join(old_mddb, "work", "ds", "1R", "src.pr.dcd"), "wb") as fh:
            fh.write(b"stub")
        return old_index, old_mddb

    def test_everything_written_by_hand_survives_the_move(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_index, old_mddb = self.old_folder(tmp)
            self.assertTrue(migrate.needed(tmp))
            migrate.migrate(tmp, verbose=False)

            doc = settings.load(tmp)
            self.assertEqual(doc["folder"]["owner"], "labmate")
            self.assertEqual(doc["folder"]["systems"], {"gp120": "OPEN"})
            self.assertEqual(doc["datasets"][0]["citation"], "my paper")
            self.assertEqual(doc["ndp"]["defaults"]["lab"], "My Lab")
            self.assertEqual(doc["ndp"]["systems"][0]["title"], "My trimer")
            text = settings.read(tmp)
            self.assertIn("# my label", text)
            self.assertIn("# do not touch", text)

    def test_the_machinery_moves_and_the_originals_are_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_index, old_mddb = self.old_folder(tmp)
            result = migrate.migrate(tmp, verbose=False)
            self.assertTrue(os.path.exists(
                os.path.join(layout.registry(tmp), "runs.parquet")))
            self.assertTrue(os.path.exists(os.path.join(
                layout.mddb_dir(tmp), "work", "ds", "1R", "src.pr.dcd")))
            self.assertTrue(os.path.exists(layout.notes(tmp)))
            self.assertFalse(os.path.exists(os.path.join(old_index, "registry")))
            # Nothing is deleted: the old directories stay, and say where to look.
            self.assertTrue(os.path.exists(os.path.join(old_index, migrate.MOVED_NAME)))
            self.assertTrue(os.path.exists(os.path.join(old_mddb, migrate.MOVED_NAME)))
            self.assertEqual(sorted(os.path.basename(p) for p in result["kept"]),
                             ["datasets.yaml", "ndp.yaml", "substrate.yaml"])

    def test_running_it_twice_changes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.old_folder(tmp)
            migrate.migrate(tmp, verbose=False)
            before = settings.read(tmp)
            second = migrate.migrate(tmp, verbose=False)
            self.assertEqual(second["moved"], [])
            self.assertEqual(settings.read(tmp), before)
            self.assertFalse(migrate.needed(tmp))


class TestGeneratedPages(unittest.TestCase):
    """A run that found nothing new must leave nothing new behind.

    2026-09-22: re-running on CDR rewrote README.md and NDP_registration.md
    every time, because each carries the time it was generated -- so an
    unchanged folder looked changed to `ls -lt` and to whoever syncs it.
    """

    def test_only_the_timestamp_differing_does_not_rewrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            page = os.path.join(tmp, "MD_RESULTS", "README.md")
            layout.write_page(page, "Generated 2026-09-22T10:00:00.\n4 MDs.\n")
            before = os.stat(page).st_mtime_ns
            layout.write_page(page, "Generated 2026-09-22T11:30:05+00:00.\n4 MDs.\n")
            self.assertEqual(os.stat(page).st_mtime_ns, before)
            self.assertIn("T10:00:00", open(page).read())

    def test_a_real_change_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            page = os.path.join(tmp, "MD_RESULTS", "README.md")
            layout.write_page(page, "Generated 2026-09-22T10:00:00.\n4 MDs.\n")
            layout.write_page(page, "Generated 2026-09-22T10:00:00.\n5 MDs.\n")
            self.assertIn("5 MDs", open(page).read())




class TestRootIsChecked(unittest.TestCase):
    """`folder: root:` named a folder and was then trusted forever.

    A settings.yaml copied as a starting point, or a folder that was moved or
    renamed, made the build READ the old tree and WRITE into the new one, with
    nothing said. It is a live hazard as soon as MD_RESULTS folders are being
    copied around.
    """

    def moved_folder(self, tmp, root_value):
        d = os.path.join(tmp, "moved", "MD_RESULTS")
        os.makedirs(d)
        with open(os.path.join(d, "settings.yaml"), "w") as fh:
            fh.write("folder:\n"
                     f"  root: {root_value}\n"
                     "  owner: someone\n"
                     "  corpora: []\n"
                     "  prune: []\n")
        return os.path.dirname(d)

    def test_a_root_that_is_not_this_folder_stops_the_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            elsewhere = os.path.join(tmp, "elsewhere")
            os.makedirs(elsewhere)
            folder = self.moved_folder(tmp, elsewhere)
            with self.assertRaises(settings.RootMoved) as caught:
                settings.config_of(folder)
            # Both paths must be in the message, or it cannot be acted on.
            self.assertIn(os.path.realpath(elsewhere), str(caught.exception))
            self.assertIn(os.path.realpath(folder), str(caught.exception))
            self.assertIn("--rebase-root", str(caught.exception))

    def test_a_root_that_matches_is_left_alone(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = self.moved_folder(tmp, os.path.join(tmp, "moved"))
            cfg = settings.config_of(folder)
            self.assertEqual(os.path.realpath(cfg["root"]),
                             os.path.realpath(folder))

    def test_rebase_points_it_here_and_says_what_it_was(self):
        with tempfile.TemporaryDirectory() as tmp:
            elsewhere = os.path.join(tmp, "elsewhere")
            os.makedirs(elsewhere)
            folder = self.moved_folder(tmp, elsewhere)
            was, now, changed = settings.rebase_root(folder)
            self.assertTrue(changed)
            self.assertEqual(was, elsewhere)
            self.assertEqual(now, os.path.realpath(folder))
            settings.config_of(folder)          # now validates

    def test_rebasing_an_already_correct_root_changes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = self.moved_folder(tmp, os.path.join(tmp, "moved"))
            _, _, changed = settings.rebase_root(folder)
            self.assertFalse(changed)


class TestSoftwareHeuristic(unittest.TestCase):
    """A directory holding the science was classified as a software checkout.

    `B40T77/newversion/` holds one zipped postscript file beside the directory
    containing all twelve of its production runs. "1 of 1 loose files is an
    archive" satisfied the majority test, the whole subtree was skipped, and
    the build reported "No production runs. This folder may hold analysis
    only" -- confident, quiet, and about absence, which is the hardest kind of
    wrong to notice.
    """

    def make(self, tmp, files, dirs=()):
        d = os.path.join(tmp, "candidate")
        os.makedirs(d)
        for f in files:
            open(os.path.join(d, f), "w").close()
        for x in dirs:
            os.makedirs(os.path.join(d, x))
        return d

    def test_one_archive_beside_a_data_directory_is_not_software(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = self.make(tmp, ["24Nov14.ps.zip"], ["strutture"])
            self.assertFalse(init.looks_like_software(d))

    def test_a_real_archive_dump_still_is(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = self.make(tmp, ["a.tar.gz", "b.zip", "c.tgz", "readme.txt"])
            self.assertTrue(init.looks_like_software(d))

    def test_a_force_field_directory_still_is(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = self.make(tmp, ["par_all36.prm", "top_all36.rtf",
                                "toppar_water.str", "notes.txt"])
            self.assertTrue(init.looks_like_software(d))

    def test_an_explicit_marker_fires_however_few_the_files(self):
        """SOFTWARE_MARKERS is evidence, not a majority, so it is not gated."""
        with tempfile.TemporaryDirectory() as tmp:
            marker = sorted(init.SOFTWARE_MARKERS)[0]
            d = self.make(tmp, [marker])
            self.assertTrue(init.looks_like_software(d))


class TestRecordsDirMoves(unittest.TestCase):
    """Records change location when their count crosses three.

    `layout.records_dir` puts them beside the other files until there are more
    than three, then in `NDP_records/`. A folder whose count FALLS back leaves
    a full set of superseded records in the location it stopped using, and the
    sweep only ever looked at the location in use. Three aptamer folders each
    carried a stale set that way -- invisible to catalog.json, and counted by
    anything that globbed the tree.
    """

    def test_the_threshold_moves_the_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            inline = layout.records_dir(tmp, layout.RECORDS_INLINE)
            spilled = layout.records_dir(tmp, layout.RECORDS_INLINE + 1)
            self.assertEqual(inline, layout.results(tmp))
            self.assertEqual(spilled,
                             os.path.join(layout.results(tmp), "NDP_records"))

    def test_the_two_locations_are_the_pair_the_sweep_must_check(self):
        """Pinned as a pair: whichever is in use, the other is where the stale
        copies are. If a third location is ever added this fails, which is the
        point."""
        with tempfile.TemporaryDirectory() as tmp:
            seen = {layout.records_dir(tmp, n) for n in range(1, 12)}
            self.assertEqual(seen, {layout.results(tmp),
                                    os.path.join(layout.results(tmp),
                                                 "NDP_records")})


class TestPhase0Regressions(unittest.TestCase):
    def test_a_migrated_settings_file_keeps_the_explanations_too(self):
        text = migrate.build_settings_text()
        self.assertIn("the science behind each dataset", text)
        self.assertIn("the catalogue record", text)

    def test_two_backups_in_one_second_do_not_overwrite_each_other(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
            for v in ("a", "b", "c"):
                settings.write(tmp, f"folder:\n  root: {tmp}\n# {v}\n")
            self.assertEqual(len(os.listdir(layout.backups(tmp))), 3)
            self.assertFalse([n for n in os.listdir(layout.results(tmp))
                              if ".tmp" in n])

    def test_the_depth_limit_counts_directories_on_every_platform(self):
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "a", "b", "c", "d"))
            seen = {os.path.relpath(d, tmp)
                    for d, _, _ in config.walk(tmp, prune=(), max_depth=2)}
            self.assertEqual(seen, {".", "a", os.path.join("a", "b")})

    def test_a_build_refuses_a_folder_another_host_is_building(self):
        import json
        import ndp_build
        with tempfile.TemporaryDirectory() as tmp:
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
            lock = layout.lock(tmp, "build")
            os.makedirs(os.path.dirname(lock), exist_ok=True)
            with open(lock, "w") as fh:
                json.dump({"pid": 1, "host": "another-node"}, fh)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(ndp_build.main([tmp]), 2)
            self.assertIn("another-node", out.getvalue())
            self.assertTrue(os.path.exists(lock))      # never broken


class TestAtomicAndLock(unittest.TestCase):
    def test_a_failed_write_leaves_the_old_file_and_no_temporary(self):
        from core import atomic
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "catalog.json")
            atomic.write_text(p, "old")

            def boom(t):
                open(t, "w").write("half")
                raise RuntimeError("died mid-write")
            with self.assertRaises(RuntimeError):
                atomic.atomic(p, boom)
            self.assertEqual(open(p).read(), "old")
            self.assertEqual(os.listdir(tmp), ["catalog.json"])

    def test_a_second_build_of_the_same_folder_is_refused(self):
        import json
        import socket
        from core import atomic
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, ".build.lock")
            with atomic.lock(p):
                with self.assertRaises(atomic.Locked):
                    with atomic.lock(p):
                        pass
            self.assertFalse(os.path.exists(p))
            # A lock left by a dead process on this host is taken over.
            with open(p, "w") as fh:
                json.dump({"pid": 2 ** 22 + 12345, "host": socket.gethostname()}, fh)
            with atomic.lock(p):
                pass
            # One held from another host is never broken.
            with open(p, "w") as fh:
                json.dump({"pid": 1, "host": "elsewhere"}, fh)
            with self.assertRaises(atomic.Locked):
                with atomic.lock(p):
                    pass


if __name__ == "__main__":
    unittest.main(verbosity=2)
