#!/usr/bin/env python3
"""
test_collect.py -- the ways a collector can quietly ruin a shared folder.

    python3 tests/test_collect.py        # or: python3 -m pytest tests/

Every fixture is built in a temporary directory. Nothing here reads the lab's
data, so a failure is this code being wrong and never a folder having changed.

The failures worth pinning are the ones that leave something that still LOOKS
right: a copy whose symlinks were followed, a half-copied tree named like a
finished one, a sweep that ate a colleague's note, an index that quietly
forgot the collections it did not touch this run.
"""

import contextlib
import io
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import collect                                            # noqa: E402
from core import copy as C                                # noqa: E402
from core import index as IDX                             # noqa: E402
from core import layout as L                              # noqa: E402
from core import manifest as M                            # noqa: E402
from core import prefix as P                              # noqa: E402
from core import safety as S                              # noqa: E402


def make_source(root, name="proj", files=(("settings.yaml", "folder:\n"),),
                links=(), sub="MD_RESULTS"):
    folder = os.path.join(root, name)
    src = os.path.join(folder, sub)
    os.makedirs(src, exist_ok=True)
    for rel, text in files:
        p = os.path.join(src, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
    for rel, target in links:
        p = os.path.join(src, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        os.symlink(target, p)
    return folder, src


class Prefixes(unittest.TestCase):
    def test_declared_beats_derived(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder, src = make_source(tmp, "8-IgG", files=(
                ("settings.yaml", "collection_prefix: IgG-1F4\nfolder:\n  root: x\n"),))
            self.assertEqual(P.declared(folder), "IgG-1F4")
            got = P.resolve([folder])
            self.assertEqual(got, [(folder, "IgG-1F4", "declared")])

    def test_derived_skips_the_carrier_directories(self):
        """`MDs`, `2-REPLICA` and bare numbers say WHERE a run sits, not what
        it is. A prefix built from them names nothing."""
        self.assertEqual(P.derive("/a/ABZYME/1-FULL-GP120/MDs"),
                         "ABZYME-1-FULL-GP120")
        self.assertEqual(P.derive("/a/ABZYME/8-IgG/1-IF4"), "8-IgG-1-IF4")

    def test_a_collision_refuses_and_names_both(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, _ = make_source(tmp, "one", files=(
                ("settings.yaml", "collection_prefix: same\nfolder:\n"),))
            b, _ = make_source(tmp, "two", files=(
                ("settings.yaml", "collection_prefix: same\nfolder:\n"),))
            with self.assertRaises(P.Collision) as caught:
                P.resolve([a, b])
            self.assertIn(a, str(caught.exception))
            self.assertIn(b, str(caught.exception))

    def test_a_declared_prefix_wins_a_race_with_a_derived_one(self):
        """Resolved first on purpose: a name somebody asked for must not be
        taken by a name nobody chose."""
        with tempfile.TemporaryDirectory() as tmp:
            a, _ = make_source(tmp, "CDR", files=(("settings.yaml", "folder:\n"),))
            derived = P.derive(a)
            b, _ = make_source(tmp, "other", files=(
                ("settings.yaml", f"collection_prefix: {derived}\nfolder:\n"),))
            with self.assertRaises(P.Collision):
                P.resolve([a, b])


class Safety(unittest.TestCase):
    def test_a_source_must_be_an_md_results_a_build_wrote(self):
        with tempfile.TemporaryDirectory() as tmp:
            plain = os.path.join(tmp, "MD_RESULTS")
            os.makedirs(plain)
            with self.assertRaises(S.Refused):
                S.check_source(plain)            # no settings.yaml
            with self.assertRaises(S.Refused):
                S.check_source(os.path.join(tmp, "nope"))

    def test_destination_inside_a_source_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            with self.assertRaises(S.Refused):
                S.check_destination(os.path.join(src, "share"), [src])

    def test_a_source_inside_the_destination_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = os.path.join(tmp, "share")
            _, src = make_source(os.path.join(tmp, "share"), "proj")
            with self.assertRaises(S.Refused):
                S.check_destination(dest, [src])

    def test_the_real_lab_shape_is_accepted(self):
        """`/x/MD_RESULTS` as a source and `/x/MD_NDP_DATA` as the share are
        siblings. Testing containment against the source's FOLDER instead of
        the MD_RESULTS directory refuses every run this kit exists for."""
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "MD_RESULTS"))
            with open(os.path.join(tmp, "MD_RESULTS", "settings.yaml"), "w") as fh:
                fh.write("folder:\n")
            src = S.check_source(os.path.join(tmp, "MD_RESULTS"))
            S.check_destination(os.path.join(tmp, "MD_NDP_DATA"), [src])

    def test_a_link_target_cannot_cause_a_write_outside_the_share(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(S.Refused):
                S.check_write_target(tmp, os.path.join(tmp, "..", "escape"))


class Symlinks(unittest.TestCase):
    """`production_files/` points at the trajectories with absolute paths."""

    def test_a_link_is_copied_as_a_link_and_its_target_is_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            big = os.path.join(tmp, "huge.dcd")
            with open(big, "wb") as fh:
                fh.write(b"x" * 4096)
            _, src = make_source(tmp, "proj",
                                 links=(("production_files/a.dcd", big),))
            dirs, files, links = M.walk(src)
            self.assertIn("production_files/a.dcd", links)
            self.assertNotIn("production_files/a.dcd", files)
            dest = os.path.join(tmp, "out")
            C.build_tree(src, dest, dirs, files, links)
            made = os.path.join(dest, "production_files/a.dcd")
            self.assertTrue(os.path.islink(made))
            self.assertEqual(os.readlink(made), big)

    def test_the_target_of_a_link_is_never_read(self):
        """Following one would copy a 32 GB trajectory. The copy's size is the
        size of the real files and nothing more."""
        with tempfile.TemporaryDirectory() as tmp:
            big = os.path.join(tmp, "huge.dcd")
            with open(big, "wb") as fh:
                fh.write(b"x" * (1 << 20))
            _, src = make_source(tmp, "proj",
                                 links=(("production_files/a.dcd", big),))
            dirs, files, links = M.walk(src)
            dest = os.path.join(tmp, "out")
            C.build_tree(src, dest, dirs, files, links)
            total = sum(os.lstat(os.path.join(dp, f)).st_size
                        for dp, _, fs in os.walk(dest) for f in fs
                        if not os.path.islink(os.path.join(dp, f)))
            self.assertLess(total, 1 << 20)

    def test_a_broken_link_is_copied_not_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp, "proj",
                                 links=(("production_files/gone.dcd",
                                         "/no/such/file.dcd"),))
            dirs, files, links = M.walk(src)
            dest = os.path.join(tmp, "out")
            C.build_tree(src, dest, dirs, files, links)
            made = os.path.join(dest, "production_files/gone.dcd")
            self.assertTrue(os.path.islink(made))
            self.assertFalse(os.path.exists(made))


class Manifests(unittest.TestCase):
    def test_an_unchanged_source_has_an_unchanged_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            a = M.fingerprint(*M.walk(src))
            b = M.fingerprint(*M.walk(src))
            self.assertEqual(a, b)

    def test_one_changed_byte_changes_the_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            a = M.fingerprint(*M.walk(src))
            with open(os.path.join(src, "settings.yaml"), "w") as fh:
                fh.write("folder:\n  root: elsewhere\n")
            self.assertNotEqual(a, M.fingerprint(*M.walk(src)))

    def test_the_diff_names_what_was_added_replaced_and_swept(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp, files=(("settings.yaml", "a\n"),
                                             ("gone.json", "{}\n")))
            dirs, files, links = M.walk(src)
            old = M.build("p", "declared", src, tmp, dirs, files, links, {})
            os.unlink(os.path.join(src, "gone.json"))
            with open(os.path.join(src, "new.json"), "w") as fh:
                fh.write("{}\n")
            dirs2, files2, links2 = M.walk(src)
            added, replaced, _, swept = M.diff(old, dirs2, files2, links2)
            self.assertEqual(added, ["new.json"])
            self.assertEqual(swept, ["gone.json"])


class Sweeping(unittest.TestCase):
    def test_a_file_a_person_added_is_carried_forward_not_eaten(self):
        """The rule the metadata pipeline learned the hard way: a sweep that
        removes what it did not write is not a sweep."""
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            dirs, files, links = M.walk(src)
            previous = os.path.join(tmp, "previous")
            C.build_tree(src, previous, dirs, files, links)
            note = os.path.join(previous, "NOTE-FROM-A-COLLEAGUE.md")
            with open(note, "w") as fh:
                fh.write("do not delete me\n")
            old = M.build("p", "declared", src, tmp, dirs, files, links, {})
            carried = C.foreign(previous, old)
            self.assertEqual([rel for rel, _ in carried],
                             ["NOTE-FROM-A-COLLEAGUE.md"])
            fresh = os.path.join(tmp, "fresh")
            C.build_tree(src, fresh, dirs, files, links, carry=carried)
            self.assertTrue(os.path.exists(
                os.path.join(fresh, "NOTE-FROM-A-COLLEAGUE.md")))

    def test_what_we_did_write_is_not_reported_as_foreign(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            dirs, files, links = M.walk(src)
            previous = os.path.join(tmp, "previous")
            C.build_tree(src, previous, dirs, files, links)
            old = M.build("p", "declared", src, tmp, dirs, files, links, {})
            self.assertEqual(C.foreign(previous, old), [])


class Atomicity(unittest.TestCase):
    def test_a_failed_copy_leaves_the_previous_collection_intact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = os.path.join(tmp, "share")
            os.makedirs(root)
            final = L.collected(root, "p")
            os.makedirs(final)
            with open(os.path.join(final, "old.json"), "w") as fh:
                fh.write("previous\n")
            staged = os.path.join(tmp, "does-not-exist")
            with self.assertRaises(OSError):
                C.swap(root, "p", staged)
            self.assertTrue(os.path.exists(os.path.join(final, "old.json")),
                            "the previous copy was lost on a failed swap")

    def test_a_successful_swap_replaces_the_whole_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = os.path.join(tmp, "share")
            L.ensure(L.tmp(root))
            final = L.collected(root, "p")
            os.makedirs(final)
            with open(os.path.join(final, "old.json"), "w") as fh:
                fh.write("previous\n")
            staged = L.partial(root, "p")
            os.makedirs(staged)
            with open(os.path.join(staged, "new.json"), "w") as fh:
                fh.write("fresh\n")
            C.swap(root, "p", staged)
            self.assertTrue(os.path.exists(os.path.join(final, "new.json")))
            self.assertFalse(os.path.exists(os.path.join(final, "old.json")))


class DryRun(unittest.TestCase):
    def test_a_dry_run_creates_absolutely_nothing(self):
        """Not the share, not the lock, not a temp directory, not a log. The
        first version created MD_NDP_DATA because a path helper made the
        directory it was asked to name."""
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            dest = os.path.join(tmp, "share")
            rc = collect.main([src, "--dest", dest, "--dry-run"])
            self.assertEqual(rc, 0)
            self.assertFalse(os.path.exists(dest))


class Index(unittest.TestCase):
    def test_every_collection_in_the_share_is_listed(self):
        """Building the index from this run's rows alone dropped every
        collection the run did not touch."""
        with tempfile.TemporaryDirectory() as tmp:
            a_folder, a = make_source(tmp, "one")
            b_folder, b = make_source(tmp, "two")
            dest = os.path.join(tmp, "share")
            self.assertEqual(collect.main([a, "--dest", dest]), 0)
            self.assertEqual(collect.main([b, "--dest", dest]), 0)
            doc = json.load(open(os.path.join(dest, L.INDEX_JSON)))
            self.assertEqual(len(doc["collections"]), 2)

    def test_a_run_that_changes_nothing_touches_no_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            dest = os.path.join(tmp, "share")
            collect.main([src, "--dest", dest])
            before = {n: os.stat(os.path.join(dest, n)).st_mtime_ns
                      for n in (L.README, L.INDEX_MD, L.INDEX_JSON)}
            collect.main([src, "--dest", dest])
            after = {n: os.stat(os.path.join(dest, n)).st_mtime_ns
                     for n in (L.README, L.INDEX_MD, L.INDEX_JSON)}
            self.assertEqual(before, after)

    def test_nothing_of_ours_is_left_inside_a_collection(self):
        """A folder taken out of the share and handed to someone must be byte
        for byte what the source was."""
        with tempfile.TemporaryDirectory() as tmp:
            folder, src = make_source(tmp)
            dest = os.path.join(tmp, "share")
            collect.main([src, "--dest", dest])
            copied = L.collected(dest, P.derive(folder))
            names = set()
            for dp, _, fs in os.walk(copied):
                names.update(fs)
            self.assertEqual(names, {"settings.yaml"})


class AtScale(unittest.TestCase):
    """What breaks when a share holds thousands of collections."""

    def test_one_prefix_clash_skips_that_pair_and_collects_the_rest(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, a = make_source(tmp, "one", files=(
                ("settings.yaml", "collection_prefix: same\nfolder:\n"),))
            _, b = make_source(tmp, "two", files=(
                ("settings.yaml", "collection_prefix: same\nfolder:\n"),))
            _, c = make_source(tmp, "three")
            dest = os.path.join(tmp, "share")
            with contextlib.redirect_stdout(io.StringIO()):
                rc = collect.main([a, b, c, "--dest", dest])
            self.assertEqual(rc, 1)                       # said, not hidden
            self.assertEqual(sorted(L.collections_in(dest)), [P.derive(os.path.dirname(c))])

    def test_a_crashed_run_s_partial_tree_is_cleared(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            dest = os.path.join(tmp, "share")
            dead = os.path.join(L.tmp(dest), "proj.99999.partial")
            os.makedirs(dead)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(collect.main([src, "--dest", dest]), 0)
            self.assertFalse(os.path.exists(dead))

    def test_a_live_collector_s_lock_is_honoured(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, ".collect.lock")
            with S.Lock(path):
                with self.assertRaises(S.Refused):
                    with S.Lock(path):
                        pass
            self.assertFalse(os.path.exists(path))

    def test_a_collector_on_another_host_stops_this_one_before_it_plans(self):
        """Plans were made BEFORE the lock, so they could be stale by the time
        it was taken. Now nothing is planned, copied or indexed without it."""
        with tempfile.TemporaryDirectory() as tmp:
            _, src = make_source(tmp)
            dest = os.path.join(tmp, "share")
            os.makedirs(L.internal(dest))
            with open(L.lock(dest), "w") as fh:
                json.dump({"pid": 1, "host": "another-node"}, fh)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                self.assertEqual(collect.main([src, "--dest", dest]), 2)
            self.assertIn("another-node", out.getvalue())
            self.assertFalse(L.collections_in(dest))
            self.assertFalse(os.path.exists(os.path.join(dest, L.INDEX_JSON)))

    def test_the_full_digest_is_kept_and_old_short_ones_still_verify(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = os.path.join(tmp, "a"), os.path.join(tmp, "b")
            open(src, "wb").write(b"trajectory")
            digest, _ = C.copy_file(src, dst)
            self.assertEqual(len(digest), 64)
            self.assertEqual(C.verify(tmp, {"b": [10, 0, digest]}), [])
            self.assertEqual(C.verify(tmp, {"b": [10, 0, digest[:16]]}), [])


class TestLabName(unittest.TestCase):
    """The share is named for the lab in lab.yaml, read without PyYAML."""

    def test_short_name_from_the_lab_block_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "lab.yaml")
            with open(p, "w") as fh:
                fh.write('# header\nlab:\n  name: "Smith Lab, Elsewhere U"  # a note\n'
                         '  short_name: Smith Lab\nidentifiers:\n  short_name: NOT-THIS\n')
            self.assertEqual(IDX.lab_name(p), "Smith Lab")
            self.assertIsNone(IDX.lab_name(os.path.join(tmp, "absent.yaml")))

    def test_flow_style_and_hashes_inside_quotes(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "lab.yaml")
            with open(p, "w") as fh:
                fh.write('lab: {name: "Lab #5 Group", url: x}\n')
            self.assertEqual(IDX.lab_name(p), "Lab #5 Group")


if __name__ == "__main__":
    unittest.main(verbosity=2)
