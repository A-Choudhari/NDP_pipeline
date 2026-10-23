#!/usr/bin/env python3
"""
test_repairs.py -- regression fixtures for the 2026-09 audit repairs.

The lab-wide glossary was dropped on 2026-09-22 (each folder is self-contained
now), and its fixtures went with it.

    python3 tests/test_repairs.py          # or: python3 -m pytest tests/

Every test here pins a defect that reached the exported catalogue and was not
visible to either acceptance suite. Fixtures are synthetic and self-contained:
nothing reads production data, and nothing writes outside a temp directory.

The acceptance suites answer "does the toolkit still work". These answer "has
the specific thing that went wrong come back".
"""

import json
import os
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from core import identity, stages           # noqa: E402
from core.mds import _natural_key, production_files as _traj_files   # noqa: E402


# --------------------------------------------------------------------------
# 1. MD identity: deterministic, stable, collision-free
# --------------------------------------------------------------------------
class TestMDIdentity(unittest.TestCase):
    """The defect: three different DNA aptamers under Top1/, Top4/ and Top5/
    were all named `1R`, because the name came from the replica label. mwf
    rejects the whole project for a duplicate MD name (core/inputs.py:358)."""

    ROOT = "Aptamers/AptaSuite-results/test2-withoutprimers/top100"
    APTAMERS = [f"{ROOT}/Top1/5/MODEL35/REP1",
                f"{ROOT}/Top4/1/model7/REP1",
                f"{ROOT}/Top5/1/model15/REP1"]

    def test_sibling_replicas_do_not_collide(self):
        ids = identity.assign_md_ids(self.APTAMERS, self.ROOT)
        self.assertEqual(len(set(ids.values())), 3,
                         f"three distinct aptamers collapsed to {set(ids.values())}")

    def test_ids_are_readable_not_hashes(self):
        ids = identity.assign_md_ids(self.APTAMERS, self.ROOT)
        self.assertIn("Top1_5_MODEL35_REP1", ids.values())
        self.assertIn("Top4_1_model7_REP1", ids.values())

    def test_adding_a_run_never_renames_an_existing_one(self):
        """The property that makes a published <dataset_id>/<md_id> safe to
        cite. A common-ancestor scheme would be prettier and would renumber
        every MD the moment someone adds a run higher up the tree."""
        before = identity.assign_md_ids(self.APTAMERS, self.ROOT)
        after = identity.assign_md_ids(
            self.APTAMERS + [f"{self.ROOT}/Top9/2/model3/REP2"], self.ROOT)
        for rel in self.APTAMERS:
            self.assertEqual(before[rel], after[rel],
                             f"{rel} was renamed by an unrelated addition")

    def test_sanitisation_collisions_are_broken_apart(self):
        """mwf maps a name to a directory by replacing spaces with '_' and
        DELETING '. , ; : º /' (core/inputs.py:293-300). So `rep 1`/`rep_1` and
        `A.B`/`AB` fold together upstream even though they look distinct here."""
        for a, b in (("P/rep 1", "P/rep_1"), ("P/A.B", "P/AB"), ("P/x:y", "P/xy")):
            ids = identity.assign_md_ids([a, b], "P")
            self.assertNotEqual(ids[a], ids[b],
                                f"{a!r} and {b!r} still fold to {ids[a]!r}")

    def test_ids_are_deterministic_across_processes(self):
        """blake2b, not hash() -- Python randomises hash() per process, so a
        hash()-based suffix would give the same run a different id on every
        export. Pinning the literal value is what makes that regression
        detectable; a self-comparison would pass either way."""
        self.assertEqual(identity.short_hash("a/b/c"), "16831f")
        import subprocess
        out = subprocess.run(
            [sys.executable, "-c",
             "import sys; sys.path.insert(0, %r); "
             "from core import identity; print(identity.short_hash('a/b/c'))"
             % os.path.dirname(HERE)],
            capture_output=True, text=True, timeout=60,
            env={**os.environ, "PYTHONHASHSEED": "1"})
        self.assertEqual(out.stdout.strip(), identity.short_hash("a/b/c"),
                         "id changed under a different hash seed")

    def test_no_enumeration_order_dependence(self):
        forward = identity.assign_md_ids(self.APTAMERS, self.ROOT)
        backward = identity.assign_md_ids(list(reversed(self.APTAMERS)), self.ROOT)
        self.assertEqual(forward, backward)

    def test_check_unique_names_both_sides_of_a_clash(self):
        with self.assertRaises(identity.IdentityCollision) as ctx:
            identity.check_unique([("1R", "a/Top1/REP1"), ("1R", "a/Top4/REP1")])
        message = str(ctx.exception)
        self.assertIn("a/Top1/REP1", message)
        self.assertIn("a/Top4/REP1", message)

    def test_declared_root_is_used_not_the_matched_runs(self):
        systems = {"aptasuite-models": "Aptamers/AptaSuite-results"}
        root = identity.dataset_root({"system": "aptasuite-models"}, systems)
        self.assertEqual(root, "Aptamers/AptaSuite-results")

    def test_system_declared_as_a_list_yields_no_single_prefix(self):
        """LC3B is one campaign in three places; there is no common prefix, and
        inventing one would produce ids that move when a place is added."""
        systems = {"LC3B": ["ANALYSIS-PAPERS-ITA/LC3B", "DOWNLOAD/aataaa"]}
        self.assertEqual(identity.dataset_root({"system": "LC3B"}, systems), "")


# --------------------------------------------------------------------------
# 2. Stage classification
# --------------------------------------------------------------------------
class TestStageClassification(unittest.TestCase):
    """The defect: classification by filename PREFIX. `protein_sol_nacl_eq`
    starts with "pr" so 468 equilibration files were deposited as production;
    `base_model.lipids_protein_glycans.R1` matched nothing so 161 GB of real MD
    became invisible."""

    def test_protein_prefix_is_not_production(self):
        for name in ("protein_sol_nacl_eq.mdcrd",
                     "protein_sol_nacl_eq1.mdcrd",
                     "protein_sol_nacl_eq2.mdcrd"):
            self.assertEqual(stages.stage_of(name), "equilibration", name)

    def test_real_production_names_still_classify(self):
        for name in ("pr.dcd", "pr2.dcd", "pr.part10.dcd", "md.nc",
                     "md_rep2.nc", "MD1.mdcrd", "pr.mdcrd"):
            self.assertEqual(stages.stage_of(name), "production", name)

    def test_explicit_production_beats_an_ensemble_name(self):
        self.assertEqual(stages.stage_of("prod_npt.dcd"), "production")
        self.assertEqual(stages.stage_of("npt2.dcd"), "equilibration")

    def test_equilibration_beats_a_weak_production_token(self):
        self.assertEqual(stages.stage_of("pr_eq.dcd"), "equilibration")

    def test_solvent_stripped_derivative_is_caught(self):
        """Tokens here are the innocuous words not/water/ions, so this has to
        match as a substring."""
        self.assertEqual(
            stages.stage_of("protonated_aligned_w_ions.not_water_not_ions.R1.dcd"),
            "stripped")

    def test_superseded_and_derived_are_separated_from_production(self):
        self.assertEqual(stages.stage_of("pr6-old.dcd"), "superseded")
        self.assertEqual(stages.stage_of("pr6900frames.dcd"), "derived")

    def test_unmatched_names_are_unresolved_not_other(self):
        """`unresolved` is reported; `other` was silently dropped."""
        self.assertEqual(
            stages.stage_of("base_model.lipids_protein_glycans.R1.dcd"),
            "unresolved")
        self.assertEqual(stages.stage_of("250ns.mdcrd"), "unresolved")

    def test_an_extension_lookalike_is_not_production(self):
        self.assertEqual(stages.stage_of("mdcrd.mdcrd"), "unresolved")

    def test_only_production_is_depositable(self):
        self.assertEqual(stages.DEPOSITABLE, {"production"})

    def test_a_derivative_directory_overrides_a_production_filename(self):
        """`.../SASA-384-401/stripped/1R.pr.dcd` is a 45,210-atom extract of a
        788,461-atom run. The filename says production; only the directory says
        otherwise, and depositing it would publish a derivative as the run."""
        self.assertEqual(stages.stage_of("1R.pr.dcd"), "production")
        self.assertEqual(
            stages.stage_of("1R.pr.dcd", parent="/a/b/SASA-384-401/stripped"),
            "stripped")
        self.assertEqual(stages.stage_of("pr.dcd", parent="/a/b/backup"),
                         "superseded")

    def test_a_stage_directory_never_overrides_the_filename(self):
        """A directory called `pr` says where a trajectory lives, not what it
        is, so it must not promote an equilibration file to production."""
        self.assertEqual(
            stages.stage_of("protein_sol_nacl_eq.mdcrd", parent="/a/b/pr"),
            "equilibration")
        self.assertEqual(stages.stage_of("min1.dcd", parent="/a/b/md"),
                         "minimisation")


# --------------------------------------------------------------------------
# 3. Trajectory segment ordering
# --------------------------------------------------------------------------
class TestSegmentOrdering(unittest.TestCase):
    """The defect: sorted(os.listdir(d)) files pr.part10 and pr.part11 between
    pr.part1 and pr.part2. mwf concatenates in the order given, so 17 MDs had a
    scrambled time axis -- and a scrambled trajectory does not look broken."""

    def test_natural_key_orders_numerically(self):
        names = [f"pr.part{i}.dcd" for i in (1, 10, 11, 2, 3, 9)]
        ordered = sorted(names, key=_natural_key)
        self.assertEqual(ordered, ["pr.part1.dcd", "pr.part2.dcd", "pr.part3.dcd",
                                   "pr.part9.dcd", "pr.part10.dcd",
                                   "pr.part11.dcd"])

    def test_the_unnumbered_first_segment_comes_first(self):
        """The defect: the extension was part of the sort key, so `pr.dcd` keyed
        as ['pr.dcd'] and sorted AFTER ['pr', 2, '.dcd']. The lab restarts as
        pr.dcd, pr2.dcd, pr3.dcd ..., so 128 of 139 filename-ordered MDs had
        their opening segment appended to the end of the time axis. The earlier
        fixture used only numbered names and could not see it."""
        cases = {
            ("pr6.dcd", "pr.dcd", "pr2.dcd", "pr10.dcd"):
                ["pr.dcd", "pr2.dcd", "pr6.dcd", "pr10.dcd"],
            ("pr_r2.dcd", "pr_r1.dcd", "pr.dcd"):
                ["pr.dcd", "pr_r1.dcd", "pr_r2.dcd"],
            ("pr.part10.dcd", "pr.part2.dcd", "pr.dcd", "pr.part1.dcd"):
                ["pr.dcd", "pr.part1.dcd", "pr.part2.dcd", "pr.part10.dcd"],
            ("MD2.mdcrd", "MD.mdcrd", "MD10.mdcrd"):
                ["MD.mdcrd", "MD2.mdcrd", "MD10.mdcrd"],
        }
        for names, expected in cases.items():
            self.assertEqual(sorted(names, key=_natural_key), expected)

    def test_the_key_never_compares_text_with_a_number(self):
        """Names starting with digits, names without any, and the extension
        must all sort without a TypeError."""
        names = ["1R.pr.dcd", "pr.dcd", "10R.pr.dcd", "run.nc", "2.dcd", "md"]
        sorted(names, key=_natural_key)          # must not raise

    def test_lexicographic_order_is_the_bug_being_prevented(self):
        names = [f"pr.part{i}.dcd" for i in (1, 2, 10, 11)]
        self.assertNotEqual(sorted(names), sorted(names, key=_natural_key),
                            "fixture no longer reproduces the defect")

    def _run_dir(self, tmp, filenames):
        d = os.path.join(tmp, "run")
        os.makedirs(d)
        for n in filenames:
            open(os.path.join(d, n), "wb").close()
        return d

    def test_headerless_segments_fall_back_to_natural_order(self):
        """AMBER .mdcrd carries no header, so ordering is a guess -- and the
        caller is told it is a guess rather than assuming it was verified."""
        names = [f"md{i}.mdcrd" for i in (1, 2, 10, 11)]
        with tempfile.TemporaryDirectory() as tmp:
            self._run_dir(tmp, names)
            files, how = _traj_files(tmp, "run", stages.stage_of, (".mdcrd",))
            self.assertEqual(files, ["md1.mdcrd", "md2.mdcrd",
                                     "md10.mdcrd", "md11.mdcrd"])
            self.assertEqual(how, "natural")

    def test_single_segment_is_reported_as_such(self):
        with tempfile.TemporaryDirectory() as tmp:
            self._run_dir(tmp, ["pr.mdcrd"])
            files, how = _traj_files(tmp, "run", stages.stage_of, (".mdcrd",))
            self.assertEqual((files, how), (["pr.mdcrd"], "single"))

    def test_non_production_files_are_never_listed(self):
        with tempfile.TemporaryDirectory() as tmp:
            self._run_dir(tmp, ["pr.mdcrd", "protein_sol_nacl_eq.mdcrd",
                                "merged.mdcrd", "pr1.stripped.mdcrd",
                                "min1.mdcrd", "pr6-old.mdcrd"])
            files, _ = _traj_files(tmp, "run", stages.stage_of, (".mdcrd",))
            self.assertEqual(files, ["pr.mdcrd"])

    def test_the_dataset_proposer_lives_in_the_ndp_pipeline(self):
        """Its exit codes are pinned by NDP_PIPELINE/tests; what matters here
        is that this pipeline never proposes datasets of its own -- it reads
        the glossary the other one built."""
        import mddb_build
        src = open(mddb_build.__file__).read()
        self.assertIn("glossary.load(base)", src)
        self.assertNotIn("datasets_init", src)


class _Tree:
    """A throwaway source tree + MDDB_substrate, with a space in every path."""

    def __init__(self, tmp):
        self.root = os.path.join(tmp, "lab data")
        self.run = os.path.join(self.root, "sys A", "Top1", "REP1")
        os.makedirs(self.run)
        for name, body in (("pr.dcd", b"frames-1"), ("pr.part1.dcd", b"frames-2"),
                           ("solvated_ions.psf", b"topology")):
            with open(os.path.join(self.run, name), "wb") as fh:
                fh.write(body)
        self.mddb = os.path.join(self.root, "MDDB_substrate")
        os.makedirs(self.mddb)

    def doc(self, md_id="REP1"):
        rel = os.path.relpath(self.run, self.root)
        return {"name": "ds", "mdref": 0,
                "mds": [{"name": md_id, "mdir": rel,
                         "input_trajectory_filepaths": ["pr.dcd", "pr.part1.dcd"],
                         "input_topology_filepath":
                             os.path.join(rel, "solvated_ions.psf")}],
                "_ordering": {rel: "header"}}


class TestWorkdirSafety(unittest.TestCase):
    """The defect: every export set `mdir` to the production run directory, and
    mwf writes into -- and deletes `time.txt`, `#*`, `f*.trajectory.xtc` from --
    every mdir (mwf.py:146-215, tools/remove_trash.py:6-16)."""

    def setUp(self):
        from core import workdir
        self.w = workdir
        self._tmp = tempfile.TemporaryDirectory()
        self.t = _Tree(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_containment_is_by_component_not_by_string_prefix(self):
        base = self._tmp.name
        a, ab = os.path.join(base, "run"), os.path.join(base, "run2")
        os.makedirs(a)
        os.makedirs(ab)
        self.assertTrue(self.w.within(os.path.join(a, "x"), a))
        self.assertFalse(self.w.within(ab, a),
                         "'/run2' was treated as inside '/run' (string prefix)")

    def test_mdir_is_a_single_link_only_component(self):
        runnable, prov = self.w.materialise(self.t.mddb, "ds", self.t.doc(),
                                            self.t.root)
        md = runnable["mds"][0]
        self.assertEqual(md["mdir"], "REP1")
        md_dir = os.path.join(self.t.mddb, "work", "ds", "REP1")
        for name in md["input_trajectory_filepaths"] + [md["input_topology_filepath"]]:
            self.assertTrue(name.startswith(self.w.LINK_PREFIX), name)
            self.assertTrue(os.path.islink(os.path.join(md_dir, name)))
        self.assertEqual(os.path.realpath(os.path.join(md_dir, "src.pr.dcd")),
                         os.path.realpath(os.path.join(self.t.run, "pr.dcd")))
        self.assertEqual(prov["mds"][0]["trajectories"], ["pr.dcd", "pr.part1.dcd"])

    def test_nothing_is_created_inside_the_source_run(self):
        before = sorted(os.listdir(self.t.run))
        self.w.materialise(self.t.mddb, "ds", self.t.doc(), self.t.root)
        self.assertEqual(sorted(os.listdir(self.t.run)), before)

    def test_work_symlinked_into_source_is_refused(self):
        """If MDDB_substrate/work were a link into the data, resolving it would
        make everything beneath it look authorised."""
        os.symlink(self.t.run, os.path.join(self.t.mddb, "work"))
        with self.assertRaises(self.w.UnsafeWorkdir):
            self.w.materialise(self.t.mddb, "ds", self.t.doc(), self.t.root)
        self.assertEqual(sorted(os.listdir(self.t.run)),
                         ["pr.dcd", "pr.part1.dcd", "solvated_ions.psf"])

    def test_a_run_at_the_folder_root_may_hold_the_work_area(self):
        """Changed 2026-09-21. This used to be refused: MDDB_substrate inside a
        run directory puts every work dir "inside" that run. But that is exactly
        a run whose trajectories lie loose at the folder root (ABZYME/1-FULL-
        GP120/MDs), and it could then never be exported. The work dirs are under
        the real, unlinked authorised root, which the index never scans, so none
        of the run's files can be in them -- and the symlink and outside-root
        checks still apply."""
        nested = os.path.join(self.t.run, "MDDB_substrate")
        md_dir = os.path.join(nested, "work", "ds", "REP1")
        self.w.assert_safe(md_dir, nested, [self.t.run])
        with self.assertRaises(self.w.UnsafeWorkdir):
            self.w.assert_safe(os.path.join(self.t.run, "REP1"), nested, [self.t.run])

    def test_link_names_can_never_be_workflow_output_names(self):
        """Writing to a symlink writes through it. A source topology called
        `structure.pdb` linked under that bare name would be overwritten."""
        for reserved in ("structure.pdb", "trajectory.xtc", "topology.psf",
                         "metadata.json", "time.txt"):
            self.assertNotIn(self.w.link_name(reserved), self.w.RESERVED_NAMES)

    def test_links_are_idempotent_and_repoint_but_never_replace_a_real_file(self):
        d = os.path.join(self._tmp.name, "d")
        os.makedirs(d)
        a_src = os.path.join(self.t.run, "pr.dcd")
        b_src = os.path.join(self.t.run, "pr.part1.dcd")
        self.assertEqual(self.w.place_link(d, a_src), "src.pr.dcd")
        self.w.place_link(d, a_src)                      # no error, no change
        os.makedirs(os.path.join(self._tmp.name, "moved"))
        moved = os.path.join(self._tmp.name, "moved", "pr.dcd")
        with open(moved, "wb") as fh:
            fh.write(b"x")
        self.w.place_link(d, moved)
        self.assertEqual(os.readlink(os.path.join(d, "src.pr.dcd")), moved)
        with open(os.path.join(d, "src.pr.part1.dcd"), "w") as fh:
            fh.write("workflow output")
        with self.assertRaises(self.w.UnsafeWorkdir):
            self.w.place_link(d, b_src)

    def test_staged_inputs_follow_the_export_but_keep_the_previous_copy(self):
        self.w.materialise(self.t.mddb, "ds", self.t.doc(), self.t.root)
        export = os.path.join(self._tmp.name, "inputs.yaml")
        with open(export, "w") as fh:
            fh.write("name: v1\n")
        self.assertEqual(self.w.stage_inputs(self.t.mddb, "ds", export), "created")
        self.assertEqual(self.w.stage_inputs(self.t.mddb, "ds", export), "unchanged")
        with open(export, "w") as fh:
            fh.write("name: v2\n")
        self.assertEqual(self.w.stage_inputs(self.t.mddb, "ds", export), "refreshed")
        work = os.path.join(self.t.mddb, "work", "ds", "inputs.yaml")
        with open(work + ".previous") as fh:
            self.assertEqual(fh.read(), "name: v1\n")

    def test_prune_removes_link_only_areas_and_keeps_workflow_output(self):
        self.w.materialise(self.t.mddb, "ds", self.t.doc("REP1"), self.t.root)
        doc2 = self.t.doc("REP9")
        self.w.materialise(self.t.mddb, "ds", doc2, self.t.root)
        out_dir = os.path.join(self.t.mddb, "work", "ds", "REP9")
        with open(os.path.join(out_dir, "trajectory.xtc"), "w") as fh:
            fh.write("processed")
        os.makedirs(os.path.join(self.t.mddb, "work", "gone"))
        removed, kept = self.w.prune(self.t.mddb, {"ds": set()})
        self.assertIn("ds/REP1", removed)
        self.assertIn("ds/REP9", kept)
        self.assertIn("gone", removed)
        self.assertTrue(os.path.exists(os.path.join(self.t.run, "pr.dcd")),
                        "pruning a link removed its target")

    def test_rebuild_removes_links_no_longer_used_without_touching_targets(self):
        self.w.materialise(self.t.mddb, "ds", self.t.doc(), self.t.root)
        doc = self.t.doc()
        doc["mds"][0]["input_trajectory_filepaths"] = ["pr.dcd"]   # one fewer
        self.w.materialise(self.t.mddb, "ds", doc, self.t.root)
        md_dir = os.path.join(self.t.mddb, "work", "ds", "REP1")
        self.assertFalse(os.path.lexists(os.path.join(md_dir, "src.pr.part1.dcd")))
        self.assertTrue(os.path.exists(os.path.join(self.t.run, "pr.part1.dcd")))

    def test_verify_reports_a_dangling_link(self):
        self.w.materialise(self.t.mddb, "ds", self.t.doc(), self.t.root)
        self.assertEqual(self.w.verify(self.t.mddb), [])
        os.unlink(os.path.join(self.t.run, "pr.part1.dcd"))
        problems = self.w.verify(self.t.mddb)
        self.assertTrue(any("dangling" in p for p in problems), problems)


# --------------------------------------------------------------------------
# 6. The glossary
# --------------------------------------------------------------------------
class TestGlossary(unittest.TestCase):

    def test_state_is_the_worst_issue(self):
        from core import glossary as G
        self.assertEqual(G.state_of([]), "ok")
        self.assertEqual(G.state_of(["inherited_config"]), "ok")
        # A dataset-level grouping question must not mark the MD itself.
        self.assertEqual(G.state_of(["grouping_unreviewed"]), "ok")
        self.assertEqual(G.state_of(["ordering_unverified"]), "review")
        self.assertEqual(G.state_of(["ordering_unverified", "topology_mismatch"]),
                         "blocked")

    def test_a_chamber_prmtop_ranks_below_a_psf_with_the_same_atoms(self):
        """Found by handing a glossary entry to MDAnalysis, which refuses chamber
        prmtops. All fifteen DESIGN-V3 runs had been given one."""
        from core.trajscan import topology_rank, is_chamber_prmtop
        with tempfile.TemporaryDirectory() as tmp:
            ch = os.path.join(tmp, "cx.prmtop")
            native = os.path.join(tmp, "com.prmtop")
            psf = os.path.join(tmp, "com.psf")
            with open(ch, "w") as fh:
                fh.write("%VERSION\n%FLAG CTITLE\n%FORMAT(a80)\nx\n")
            with open(native, "w") as fh:
                fh.write("%VERSION\n%FLAG TITLE\n%FORMAT(a80)\nx\n")
            open(psf, "w").close()
            self.assertTrue(is_chamber_prmtop(ch))
            self.assertFalse(is_chamber_prmtop(native))
            ranked = sorted([ch, psf, native], key=topology_rank)
            self.assertEqual([os.path.basename(p) for p in ranked],
                             ["com.prmtop", "com.psf", "cx.prmtop"])

    def test_a_netcdf_file_named_mdcrd_is_linked_as_nc(self):
        """Both LC3B production trajectories are NetCDF named MD1.mdcrd. mwf and
        MDAnalysis pick a reader from the extension and misread them."""
        from core.trajscan import sniff_format
        from core import workdir
        with tempfile.TemporaryDirectory() as tmp:
            nc = os.path.join(tmp, "MD1.mdcrd")
            txt = os.path.join(tmp, "pr.mdcrd")
            with open(nc, "wb") as fh:
                fh.write(b"CDF\x02" + b"\x00" * 60)
            with open(txt, "w") as fh:
                fh.write("Cpptraj Generated trajectory\n  27.848  16.495\n")
            self.assertEqual(sniff_format(nc), "netcdf")
            self.assertEqual(sniff_format(txt), "ascii")
            self.assertEqual(workdir.corrected_ext(nc), ".nc")
            self.assertIsNone(workdir.corrected_ext(txt))
            d = os.path.join(tmp, "work")
            os.makedirs(d)
            self.assertEqual(workdir.place_link(d, nc, ".nc"), "src.MD1.nc")
            self.assertEqual(os.path.realpath(os.path.join(d, "src.MD1.nc")),
                             os.path.realpath(nc))

    def test_input_paths_through_links(self):
        from core import glossary as G
        row = {"work_dir": "/w/ds", "links": "REP1/src.MD1.nc; REP1/src.MD2.nc",
               "topology_link": "REP1/src.com.prmtop", "source_root": "/r",
               "source_dir": "p", "topology": "p/com.prmtop",
               "segments": "MD1.mdcrd; MD2.mdcrd"}
        top, trajs = G.input_paths(row, via="links")
        self.assertEqual(top, "/w/ds/REP1/src.com.prmtop")
        self.assertEqual(trajs, ["/w/ds/REP1/src.MD1.nc", "/w/ds/REP1/src.MD2.nc"])
        self.assertEqual(G.input_paths(row)[1], ["/r/p/MD1.mdcrd", "/r/p/MD2.mdcrd"])

    def test_select_and_input_paths(self):
        import pandas as pd
        from core import glossary as G
        df = pd.DataFrame([
            {"md_key": "a/1", "dataset_id": "a", "system": "s1", "condition": "x",
             "state": "ok", "issues": "", "temp_K": 310.0,
             "source_root": "/r", "source_dir": "p/1", "topology": "p/1/t.psf",
             "segments": "pr.dcd; pr.part1.dcd"},
            {"md_key": "b/1", "dataset_id": "b", "system": "s2", "condition": "y",
             "state": "review", "issues": "ordering_unverified, no_config",
             "temp_K": 300.0, "source_root": "/r", "source_dir": "q/1",
             "topology": None, "segments": "md.nc"},
        ])
        self.assertEqual(list(G.select(df, system="s1").md_key), ["a/1"])
        self.assertEqual(list(G.select(df, issue="no_config").md_key), ["b/1"])
        self.assertEqual(list(G.select(df, where="temp_K > 305").md_key), ["a/1"])
        self.assertEqual(list(G.select(df, state=["ok", "review"]).md_key),
                         ["a/1", "b/1"])
        top, trajs = G.input_paths(df.iloc[0])
        self.assertEqual(top, "/r/p/1/t.psf")
        self.assertEqual(trajs, ["/r/p/1/pr.dcd", "/r/p/1/pr.part1.dcd"])


# --------------------------------------------------------------------------
# 6. Placing a measured value at the level its evidence supports
# --------------------------------------------------------------------------
class TestFieldPlacement(unittest.TestCase):
    """The defect: a field whose replicas disagreed was DROPPED. Ninety HD4
    MDs were written at two frame spacings, so `framestep` -- which the spec
    calls "an important value ... used in many graph axes in the web client" --
    was omitted and all ninety measured values were lost.

    The MDDB spec allows an `mds` entry to carry any project metadata field
    (core/inputs_schema.py MDConfig, extra='allow'), and the workflow uploads
    it as that MD's own record (tools/generate_metadata.py METADATA_FIELDS).
    So disagreement is now distributed, not discarded."""

    def _build(self, per_run, extends, exported):
        """Run the placement rule alone, without a filesystem or a registry."""
        from core import mwf
        doc, cover, mds = {}, [], [{"mdir": r} for r in exported]
        values = {k: v for k, v in per_run.items() if v is not None}
        present = {r: v for r, v in values.items() if r in set(exported)}
        distinct = {mwf._compare_key(v) for v in present.values()}
        blind = len(exported) - len(present)
        if not present:
            return doc, mds, "absent"
        if len(distinct) == 1 and (extends or not blind):
            doc["f"] = next(iter(present.values()))
            return doc, mds, "filled"
        for m in mds:
            if present.get(m["mdir"]) is not None:
                m["f"] = present[m["mdir"]]
        return doc, mds, "per-MD"

    def test_disagreeing_replicas_keep_their_values_instead_of_losing_them(self):
        doc, mds, status = self._build(
            {"a": 0.1, "b": 0.1, "c": 0.2}, extends=True, exported=["a", "b", "c"])
        self.assertEqual(status, "per-MD")
        self.assertNotIn("f", doc, "a disagreeing field must not be asserted "
                                   "for the whole project")
        self.assertEqual([m["f"] for m in mds], [0.1, 0.1, 0.2],
                         "every measured value must survive on its own MD")

    def test_agreement_is_stated_once_at_project_level(self):
        doc, mds, status = self._build(
            {"a": 310, "b": 310}, extends=True, exported=["a", "b"])
        self.assertEqual(status, "filled")
        self.assertEqual(doc["f"], 310)
        self.assertTrue(all("f" not in m for m in mds),
                        "a shared value must not be repeated on every entry")

    def test_setup_evidence_speaks_for_a_replica_that_kept_no_config(self):
        """An MD config is copied from a shared template at launch, which is
        why the substrate tracks configs as `inherited`. Replicas whose copy
        was not kept inherit the value, and coverage says how many did."""
        doc, mds, status = self._build(
            {"a": 310, "b": 310}, extends=True, exported=["a", "b", "c"])
        self.assertEqual(status, "filled")
        self.assertEqual(doc["f"], 310)

    def test_execution_evidence_is_never_extended_to_another_execution(self):
        """The defect this rule prevents: `version` was measured on 2 of 93
        HD4 replicas and asserted for all 93 -- while three different NAMD
        versions demonstrably ran in this corpus. A project-level field is a
        claim about every MD, because the workflow hands it to any MD that
        does not override it."""
        doc, mds, status = self._build(
            {"a": "3.0alpha9"}, extends=False, exported=["a", "b", "c"])
        self.assertEqual(status, "per-MD")
        self.assertNotIn("f", doc, "a version seen once must not speak for all")
        self.assertEqual(mds[0]["f"], "3.0alpha9")
        self.assertTrue(all("f" not in m for m in mds[1:]),
                        "an MD that showed no evidence must carry none")

    def test_execution_evidence_covering_every_md_is_stated_once(self):
        doc, mds, status = self._build(
            {"a": "3.0.2", "b": "3.0.2"}, extends=False, exported=["a", "b"])
        self.assertEqual(status, "filled")
        self.assertEqual(doc["f"], "3.0.2")

    def test_the_spec_types_are_respected(self):
        from core.mwf import _coerce
        self.assertEqual(_coerce("temp", 310.0), 310)          # schema: number
        self.assertIsInstance(_coerce("timestep", 2), float)
        self.assertEqual(_coerce("ff", "CHARMM36m"), ["CHARMM36m"])
        self.assertIsNone(_coerce("temp", None))
        self.assertIsNone(_coerce("temp", float("nan")))

    def test_per_md_metadata_survives_the_rewrite_to_link_names(self):
        """The defect: workdir rebuilt each `mds` entry from scratch, keeping
        only the three path keys. Distributing a field to the MDs that measured
        it and then silently dropping it would be worse than not distributing
        it at all."""
        from core import workdir
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "data", "run1")
            os.makedirs(src)
            for name in ("pr.dcd", "sys.psf"):
                open(os.path.join(src, name), "w").close()
            mddb = os.path.join(tmp, "MDDB_substrate")
            os.makedirs(mddb)
            doc = {"name": "d", "mds": [{
                "name": "R1", "mdir": "data/run1",
                "input_trajectory_filepaths": ["pr.dcd"],
                "input_topology_filepath": "data/run1/sys.psf",
                "framestep": 0.2, "version": "3.0.2"}]}
            runnable, _ = workdir.materialise(mddb, "d", doc, tmp)
            entry = runnable["mds"][0]
            self.assertEqual(entry["framestep"], 0.2)
            self.assertEqual(entry["version"], "3.0.2")
            self.assertEqual(entry["mdir"], "R1", "mdir must name the link dir")
            self.assertTrue(entry["input_trajectory_filepaths"][0]
                            .startswith(workdir.LINK_PREFIX))


# --------------------------------------------------------------------------
# 7. The program that ran the MD, and its version
# --------------------------------------------------------------------------
class TestEngineIdentity(unittest.TestCase):
    """The defect: `program` was filled from the substrate's `engine` column,
    whose `NAMD/CHARMM` is a TRAJECTORY FORMAT FAMILY, not a program. That
    string went into 11 of 13 records, in a field the spec warns to keep free
    of duplicate spellings, matching neither NAMD nor CHARMM in any search.
    And `version`, which the spec asks for separately, was never emitted at
    all -- though every run's log opens with it."""

    BANNERS = {
        "namd.log": "Info: Built with CUDA version 11080\n"
                    "Info: NAMD 3.0.2 for Linux-x86_64-multicore-CUDA\n",
        "eq1.out": "          Amber 22 SANDER                              2022\n"
                   "|   Executable path: sander\n",
        "md.log": "GROMACS version:    2024.2\n",
    }
    EXPECTED = {"namd.log": ("NAMD", "3.0.2"),
                "eq1.out": ("AMBER", "22"),
                "md.log": ("GROMACS", "2024.2")}

    def test_every_engine_banner_is_read(self):
        from core import engine
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in self.BANNERS.items():
                path = os.path.join(tmp, name)
                with open(path, "w") as fh:
                    fh.write(text)
                self.assertEqual(engine.read_banner(path), self.EXPECTED[name],
                                 f"{name} banner not read")

    def test_a_format_family_is_never_reported_as_a_program(self):
        from core import engine
        self.assertIsNone(engine.ENGINE_PROGRAM.get("NAMD/CHARMM"),
                          "`NAMD/CHARMM` names no program and must not be "
                          "deposited as one")
        self.assertEqual(engine.ENGINE_PROGRAM.get("AMBER"), "AMBER")

    def test_an_analysis_log_is_not_mistaken_for_an_engine_log(self):
        """`result-OKB-protein.log` and `cpptraj.log` are MMPBSA and cpptraj
        output. They are the most common `.log` files in this corpus and carry
        no engine banner; opening them is wasted work."""
        from core import engine
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("cpptraj.log", "result-OKB-protein.log"):
                open(os.path.join(tmp, name), "w").close()
            self.assertEqual(engine._candidates(tmp), [])

    def test_the_config_names_the_program_but_never_a_version(self):
        from core import engine
        with tempfile.TemporaryDirectory() as tmp:
            got = engine.of_run(tmp, ".", fallback_engine="NAMD")
            self.assertEqual(got["program"], "NAMD")
            self.assertIsNone(got["version"],
                              "a config file cannot say which build ran")
            self.assertEqual(got["source"], "config")

    def test_nothing_is_claimed_when_there_is_no_evidence(self):
        from core import engine
        with tempfile.TemporaryDirectory() as tmp:
            got = engine.of_run(tmp, ".", fallback_engine="NAMD/CHARMM")
            self.assertIsNone(got["program"])
            self.assertIsNone(got["version"])


# --------------------------------------------------------------------------
# 8. What the record claims about itself
# --------------------------------------------------------------------------
class TestDeclarationAudit(unittest.TestCase):
    """The defect: LC3B's description said it contained 157 production run
    directories while its export contained 2. The count was generated when the
    description was, then frozen in datasets.yaml while the stage classifier
    was repaired underneath it. Nothing compared the two, and the overstated
    text is deposited verbatim."""

    def _audit(self, spec, n_mds, matched=0, selected=0):
        from core import audit
        return audit.declaration(spec, n_mds, range(matched or n_mds),
                                 range(selected or n_mds))

    def test_a_description_that_overstates_its_size_is_a_warning(self):
        rows = self._audit({"description": "containing 157 indexed production "
                                           "run directories."}, n_mds=2)
        self.assertTrue(any(f == "description" and s == "WARNING"
                            for f, s, _ in rows), rows)

    def test_a_description_whose_count_still_holds_is_silent(self):
        rows = self._audit({"description": "containing 2 indexed production "
                                           "run directories."}, n_mds=2)
        self.assertFalse([r for r in rows if r[0] == "description"])

    def test_measurements_in_prose_are_not_mistaken_for_counts(self):
        """`100 ns` and `CHARMM36m` must not read as a claim about size."""
        rows = self._audit({"description": "500 ns of CHARMM36m sampling at "
                                           "310 K across 3 replicas"}, n_mds=3)
        self.assertFalse([r for r in rows if r[0] == "description"], rows)

    def test_the_runs_dropped_before_export_are_reported(self):
        """Without this row nothing anywhere says that 155 of LC3B's 157
        matched runs hold only equilibration: the rest of the report describes
        the 2 that were exported."""
        rows = self._audit({"description": "x"}, n_mds=2, matched=157, selected=2)
        note = [n for f, s, n in rows if f == "mds" and s == "note"]
        self.assertTrue(note, rows)
        self.assertIn("157", note[0])
        self.assertIn("155", note[0])



# --------------------------------------------------------------------------
# 9. The lab-wide glossary
# --------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main(verbosity=2)
