#!/usr/bin/env python3
"""
test_protocol.py -- the preparation protocol is read from the stages' own files.

    python3 tests/test_protocol.py        # or: python3 -m pytest tests/

A synthetic NAMD run directory in the shape of ABZYME/3-OCCLUDED-OPEN/MDs/
2-REPLICA/1R: min1 -> heating -> npt1 -> pr -> pr_r1, plus a config prepared
and never run. Nothing reads production data; nothing writes outside a temp dir.
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
# The lab the tests run as, whatever the kit's own lab.yaml says.
os.environ["MD_PIPELINE_LAB"] = os.path.join(HERE, "lab_fixture.yaml")
os.environ.pop("NDP_ID_PREFIX", None)
sys.path.insert(0, HERE)

from core import protocol                        # noqa: E402
from core.glossary import overlap_frames          # noqa: E402
from test_default_production import write_dcd         # noqa: E402

HEAD = """set temp {temp}
set outname {out}
structure solvated_ions.psf
coordinates solvated_ions.pdb
outputName $outname
seed {seed}
firsttimestep {first}
timestep 2.0
langevin on
langevinTemp $temp
"""
BINCOORD = "if {{1}} {{\nbinCoordinates {inp}.coor\n}}\n"
PISTON = "if {1} {\nlangevinPiston on\nlangevinPistonTarget 1.01325\n}\n"
RESTRAIN = "if {1} {\nconstraints on\nconsref ref.pdb\nconskfile bb.cnst\nconskcol B\n}\n"
RAMP = ("for {set temp 10} {$temp <= 310} {incr temp 25} {\n  langevinTemp $temp\n"
        "  run 100\n}\n")


def config(d, name, out, body, temp=310, seed=1, first=0, inp=None):
    text = HEAD.format(temp=temp, out=out, seed=seed, first=first)
    if inp:
        text += BINCOORD.format(inp=inp)
    with open(os.path.join(d, name), "w") as fh:
        fh.write(text + body)


def log(d, out, first, last, completed=True):
    with open(os.path.join(d, f"{out}.log"), "w") as fh:
        if first:
            fh.write(f"Info: FIRST TIMESTEP         {first}\n")
        fh.write(f"ENERGY:  {last}   1.0   2.0\n")
        if completed:
            fh.write("[Partition 0][Node 0] End of program\n")


class TestProtocol(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        d = self.d = self._tmp.name
        config(d, "min1.conf", "min1", RESTRAIN + "minimize 5019\n", temp=10)
        log(d, "min1", 0, 5019)
        config(d, "heating.conf", "heating", PISTON + RESTRAIN + RAMP, temp=10, inp="min1")
        log(d, "heating", 0, 1300)
        config(d, "NPT1.conf", "npt1", PISTON + "run 252000\n", inp="heating")
        log(d, "npt1", 0, 252000)
        config(d, "PR.conf", "pr", PISTON + "run 500000\n", inp="npt1")
        log(d, "pr", 0, 250300, completed=False)
        write_dcd(os.path.join(d, "pr.dcd"), n_frames=5, istart=50000, nsavc=50000)
        config(d, "pr_r1.conf", "pr_r1", PISTON + "run 250000\n", inp="pr.restart",
               first=249998)
        log(d, "pr_r1", 249998, 500000)
        write_dcd(os.path.join(d, "pr_r1.dcd"), n_frames=6, istart=250000, nsavc=50000)
        config(d, "PR2.conf", "pr2", PISTON + "run 1\n", inp="pr")    # never run

    def tearDown(self):
        self._tmp.cleanup()

    def test_stages_follow_the_configs_chain_and_skip_what_never_ran(self):
        names = [s["name"] for s in protocol.stages_of(self.d)]
        self.assertEqual(names, ["min1", "heating", "npt1", "pr", "pr_r1"])

    def test_each_stage_says_what_it_did(self):
        st = {s["name"]: s for s in protocol.stages_of(self.d)}
        self.assertEqual((st["min1"]["stage"], st["min1"]["steps"],
                          st["min1"]["temperature_K"], st["min1"]["ensemble"]),
                         ("minimisation", 5019, None, None))
        self.assertEqual(st["min1"]["restraints"]["force_constants"], "bb.cnst")
        self.assertEqual(st["heating"]["temperature_ramp_K"],
                         {"from": 10, "to": 310, "step": 25})
        self.assertIsNone(st["heating"]["temperature_K"])
        self.assertEqual((st["npt1"]["ensemble"], st["npt1"]["pressure_bar"],
                          st["npt1"]["duration_ns"]), ("NPT", 1.01325, 0.504))
        self.assertIsNone(st["npt1"]["restraints"])

    def test_steps_and_status_come_from_the_log(self):
        st = {s["name"]: s for s in protocol.stages_of(self.d)}
        self.assertEqual(st["pr"]["status"], "interrupted at step 250300")
        self.assertEqual(st["pr_r1"]["steps"], 500000 - 249998)
        self.assertEqual(st["pr_r1"]["started_from"], "pr")

    def test_replicas_merge_and_differences_are_kept_per_replica(self):
        a = protocol.stages_of(self.d)
        b = [dict(s, seed=2) for s in a]
        merged = {s["name"]: s for s in protocol.merge({"1R": a, "2R": b})}
        self.assertEqual(merged["npt1"]["seed"], {"per_replica": {"1R": 1, "2R": 2}})
        self.assertEqual(merged["heating"]["temperature_ramp_K"],
                         {"from": 10, "to": 310, "step": 25})
        self.assertIn("10 to 310 K", protocol.summary(list(merged.values())))

    def test_a_restart_that_rewrites_a_step_is_counted_once(self):
        """pr ends at step 250,000; pr_r1 (firsttimestep 249,998) begins there."""
        pr, r1 = (os.path.join(self.d, f) for f in ("pr.dcd", "pr_r1.dcd"))
        self.assertEqual(overlap_frames([pr, r1]), 1)

    def test_a_window_that_died_after_one_frame_is_counted_once(self):
        """CDR 3-YZ23/3R: `pr.part10.dcd` holds a single frame at step
        159,100,000 and `pr.part11.dcd` resumes AT that step. The stub's frame
        is already in the next segment, and the segments around it still
        overlap by one frame each -- all of which used to be dropped."""
        a, b, c = (os.path.join(self.d, f) for f in ("s1.dcd", "s2.dcd", "s3.dcd"))
        write_dcd(a, n_frames=5, istart=100000, nsavc=100000)   # 100k..500k
        write_dcd(b, n_frames=1, istart=500000, nsavc=100000)   # the stub
        write_dcd(c, n_frames=4, istart=500000, nsavc=100000)   # 500k..800k
        self.assertEqual(overlap_frames([a, b, c]), 2)

    def test_order_does_not_change_the_count(self):
        a, b, c = (os.path.join(self.d, f) for f in ("s1.dcd", "s2.dcd", "s3.dcd"))
        write_dcd(a, n_frames=5, istart=100000, nsavc=100000)
        write_dcd(b, n_frames=1, istart=500000, nsavc=100000)
        write_dcd(c, n_frames=4, istart=500000, nsavc=100000)
        self.assertEqual(overlap_frames([c, a, b]), 2)

    def test_a_stub_window_does_not_make_the_order_unverifiable(self):
        """CDR 3-YZ23/3R: pr.part10 holds one frame at the step pr.part11
        resumes from. Both start at 159,100,000, and the run is still ordered:
        at one save interval the stub sits inside the longer segment."""
        import shutil
        from core import mds                                  # noqa: PLC0415
        d = os.path.join(self.d, "run")
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        write_dcd(os.path.join(d, "pr.dcd"), n_frames=5, istart=100000, nsavc=100000)
        write_dcd(os.path.join(d, "pr.part10.dcd"), n_frames=1,
                  istart=600000, nsavc=100000)
        write_dcd(os.path.join(d, "pr.part11.dcd"), n_frames=4,
                  istart=600000, nsavc=100000)
        files, how = mds.production_files(d, "", lambda n, p: "production", ".dcd")
        self.assertEqual(how, "header")
        self.assertEqual(files, ["pr.dcd", "pr.part10.dcd", "pr.part11.dcd"])

    def test_a_reset_step_counter_is_not_an_overlap(self):
        a, b = (os.path.join(self.d, f) for f in ("a.dcd", "b.dcd"))
        write_dcd(a, n_frames=5, istart=50000, nsavc=50000)
        write_dcd(b, n_frames=5, istart=50000, nsavc=50000)
        self.assertEqual(overlap_frames([a, b]), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
