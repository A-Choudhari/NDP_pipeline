#!/usr/bin/env python3
"""
test_binding.py -- one test per way an ingested free energy could be wrong.

    python3 tests/test_binding.py

Synthetic fixtures only: every `FINAL_MMGBSA.dat` here is written by this file
into a temporary directory, so what is tested is the READER, never the lab's
data. The one real number that appears is the DELTA TOTAL line of
`MMGBSA/A/3-R3/FINAL_MMGBSA.dat` (-5.0019, 11.1274, 0.7010), copied in as a
fixture: if the parser stops reproducing those three numbers off a file shaped
like that one, it has stopped reading MMPBSA.py output.

The failures below are the silent kind. An end-point free energy is a small
negative number in kcal/mol, and so is almost everything near it in the file --
the Complex block's TOTAL is five figures, but the Receptor-minus-Ligand
arithmetic, the per-frame series and the components all produce plausible
numbers. Nothing raises when the wrong one is taken.
"""

import math
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from core import spec as SP                                 # noqa: E402
from stages.contacts import binding as B                    # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
#: The Differences block of a real file, with the numbers left as parameters.
#: Component names with spaces (`1-4 VDW`, `DELTA G gas`) are kept because they
#: are what breaks a split-on-whitespace parser.
_DIFF = """\
Differences (Complex - Receptor - Ligand):
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
BOND                        -0.0000                0.0000              0.0000
1-4 VDW                      0.0000                0.0001              0.0000
VDWAALS                  {vdw:11.4f}            {vdwsd:11.4f}         {vdwsem:9.4f}
EEL                         34.4358               44.0192              2.7730
EGB                         45.7097               47.7924              3.0106
ESURF                      -10.9106                3.0544              0.1924

DELTA G gas                -39.8010               54.0393              3.4042
DELTA G solv                34.7991               45.7033              2.8790

DELTA TOTAL              {tot:11.4f}            {sd:11.4f}         {sem:9.4f}
"""

#: What sits ABOVE the differences block, and must never be read as one. The
#: Complex block has its own `TOTAL` line, five figures large; a parser that
#: takes the last `TOTAL` in the file, or the first, gets this.
_ABSOLUTE = """\
|Calculations performed using {n}.0 complex frames.
|
|All units are reported in kcal/mole.
-------------------------------------------------------------------------------

GENERALIZED BORN:

Complex:
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
VDWAALS                  -6772.4139               51.2982              3.2315
EEL                     -71019.9487              221.9662             13.9826

G gas                    26712.4745              204.6207             12.8899
G solv                  -13176.6624              165.8082             10.4449

TOTAL                    13535.8121              104.8069              6.6022


Receptor:
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
TOTAL                    15305.8733               92.9839              5.8574


Ligand:
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
TOTAL                    -1765.0593               46.0433              2.9005


"""


def write_final(path, total=-5.0019, sd=11.1274, sem=0.7010, n=252,
                vdw=-74.2368, vdwsd=19.9652, vdwsem=1.2577, with_delta=True):
    """A file shaped like MMPBSA.py's, with the DELTA TOTAL of our choosing."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    text = _ABSOLUTE.format(n=n)
    if with_delta:
        text += _DIFF.format(tot=total, sd=sd, sem=sem,
                             vdw=vdw, vdwsd=vdwsd, vdwsem=vdwsem)
    else:
        # A run killed part-way: the absolute blocks are written, the
        # differences are not. This is what a truncated file looks like.
        text += "\n"
    with open(path, "w") as fh:
        fh.write(text)
    return path


def write_per_frame(path, deltas):
    """A PER_FRAME.csv whose DELTA block holds the values given."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    head = "Frame #,VDWAALS,EEL,DELTA TOTAL"
    lines = ["GENERALIZED BORN:", "Complex Energy Terms", head]
    lines += [f"{i},0.0,0.0,{v + 1000:.4f}" for i, v in enumerate(deltas)]
    lines += ["", "DELTA Energy Terms", head]
    lines += [f"{i},-1.0,2.0,{v:.4f}" for i, v in enumerate(deltas)]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def write_in(path, igb=2, saltcon=0.150, idecomp=2):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("MM/GBSA, synthetic.\n"
                 "&general\n   startframe=1, interval=1, verbose=2,\n/\n"
                 f"&gb\n   igb={igb}, saltcon={saltcon:.3f},\n/\n"
                 f"&decomp\n   idecomp={idecomp}, dec_verbose=0,\n"
                 '   print_res="19,64,66-67"\n/\n')
    return path


def make_spec(tmp, md_ids=("1-R1", "2-R2", "3-R3"), label="SYN"):
    """A spec whose MDs sit at `<tmp>/0-MD/<id>/`, as the real ones do.

    Only the run DIRECTORY matters: this stage reads no trajectory, so the
    `.dcd` is an empty file and never opened.
    """
    mds = []
    for mid in md_ids:
        run = os.path.join(tmp, "0-MD", mid)
        os.makedirs(run, exist_ok=True)
        traj = os.path.join(run, "pr.dcd")
        open(traj, "a").close()
        mds.append(SP.MD(id=mid, topology=os.path.join(run, "m.psf"),
                         trajectories=[traj],
                         solvated_topology=os.path.join(run, "s.psf"),
                         solvated_trajectories=[traj], framestep_ns=0.2))
    return SP.Spec(id=f"SYN/{label}", label=label,
                   components=["protein", "antibody"], mds=mds)


def mmgbsa_tree(tmp, values, md_ids=("1-R1", "2-R2", "3-R3"),
                protomers=("A", "B", "C"), sds=None):
    """`<tmp>/0-MD/ANALYSIS/MMGBSA/<protomer>/<replica>/FINAL_MMGBSA.dat`.

    `values` is {(protomer, replica): delta_total}.
    """
    root = os.path.join(tmp, "0-MD", "ANALYSIS", "MMGBSA")
    write_in(os.path.join(root, "mmgbsa.in"))
    for (p, r), v in values.items():
        write_final(os.path.join(root, p, r, "FINAL_MMGBSA.dat"), total=v,
                    sd=(sds or {}).get((p, r), 11.0))
    return root


def rows_of(rows, **match):
    return [r for r in rows.rows
            if all(r.get(k) == v for k, v in match.items())]


# --------------------------------------------------------------------------
class TestParse(unittest.TestCase):
    """The reader, whose every failure produces a well-formed wrong number."""

    def test_the_summary_line_is_reproduced_exactly(self):
        """The whole stage rests on this: MMPBSA.py's own mean, standard
        deviation and standard error, taken verbatim and not recomputed."""
        with tempfile.TemporaryDirectory() as d:
            p = write_final(os.path.join(d, "FINAL_MMGBSA.dat"))
            got = B.parse_mmpbsa(p)["models"]["GB"]["DELTA_TOTAL"]
            self.assertEqual(got, (-5.0019, 11.1274, 0.7010))

    def test_the_differences_block_is_read_not_the_complex_block(self):
        """`TOTAL` appears four times in this file. Three of them are absolute
        energies of 13,535 / 15,305 / -1,765 kcal/mol, and taking any of them
        gives a number that is well formed, wrong by three orders of magnitude,
        and silent."""
        with tempfile.TemporaryDirectory() as d:
            p = write_final(os.path.join(d, "FINAL_MMGBSA.dat"))
            comps = B.parse_mmpbsa(p)["models"]["GB"]
            self.assertNotIn("TOTAL", comps)          # only DELTA TOTAL
            self.assertAlmostEqual(comps["DELTA_TOTAL"][0], -5.0019, 4)
            self.assertAlmostEqual(comps["VDWAALS"][0], -74.2368, 4)
            # The Complex block's VDWAALS is -6772; if it leaked in, the value
            # above would be it.
            self.assertGreater(comps["VDWAALS"][0], -1000)

    def test_a_component_name_with_spaces_survives(self):
        """`DELTA G gas` and `1-4 VDW` are one component each. A parser that
        splits on whitespace and takes field 0 as the name invents components
        called `DELTA` and `1-4`, and loses the ones it meant to keep."""
        with tempfile.TemporaryDirectory() as d:
            p = write_final(os.path.join(d, "FINAL_MMGBSA.dat"))
            comps = B.parse_mmpbsa(p)["models"]["GB"]
            self.assertIn("DELTA_G_GAS", comps)
            self.assertIn("DELTA_G_SOLV", comps)
            self.assertIn("1-4_VDW", comps)
            self.assertNotIn("DELTA", comps)

    def test_the_frame_count_is_read(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_final(os.path.join(d, "FINAL_MMGBSA.dat"), n=399)
            self.assertEqual(B.parse_mmpbsa(p)["n_frames"], 399)

    def test_a_malformed_file_yields_nothing_and_does_not_raise(self):
        """One unreadable file must not take a system's whole table with it."""
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "FINAL_MMGBSA.dat")
            with open(p, "w") as fh:
                fh.write("\x00\x01 not a data file at all\n"
                         "DELTA TOTAL  banana  apricot\n"
                         "Differences (Complex - Receptor - Ligand):\n"
                         "DELTA TOTAL  three  of  these  are  words\n")
            out = B.parse_mmpbsa(p)
            self.assertEqual(out["models"], {})
            self.assertIsNone(out["n_frames"])

    def test_a_file_that_does_not_exist_yields_nothing_and_does_not_raise(self):
        self.assertEqual(B.parse_mmpbsa("/nowhere/FINAL_MMGBSA.dat")["models"],
                         {})

    def test_the_per_frame_fallback_recomputes_and_admits_it(self):
        """Only when the summary block is missing, and the flag that says so
        is what keeps a recomputed s.e.m. -- sd/sqrt(n) over correlated frames,
        which understates the error -- from being read as the file's own."""
        vals = [-4.0, -6.0, -5.0, -5.0]
        with tempfile.TemporaryDirectory() as d:
            p = write_per_frame(os.path.join(d, "PER_FRAME.csv"), vals)
            out = B.parse_per_frame(p)
            self.assertTrue(out["recomputed"])
            mean, sd, sem = out["models"]["GB"]["DELTA_TOTAL"]
            self.assertAlmostEqual(mean, -5.0, 6)
            self.assertAlmostEqual(sd, math.sqrt(2 / 3), 6)
            self.assertAlmostEqual(sem, sd / 2.0, 6)
            # and the Complex block above it (values + 1000) is not what was
            # read: that block is not the difference.
            self.assertLess(mean, 0)


#: A CaFE result log, in the layout its `*_result.log` really uses: `name:
#: frames mean sd`, no standard error, and a Delta block that follows three
#: absolute blocks of the same shape.
_CAFE = """\
 ------------------------------------------------------------
   Title           Frames            Mean              SD
 ------------------------------------------------------------
 Complex:
   Elec:                4     -35986.1454         21.4894
   Total:               4     -47051.6993         26.6239
 ------------------------------------------------------------
 Receptor:
   Total:               4     -35171.1558         13.3743
 ------------------------------------------------------------
{delta} ------------------------------------------------------------
 * All energy values are in kcal/mol
"""

_CAFE_DELTA = """\
 Delta:
   Elec:                4         36.6817          5.1359
   Vdw:                 4        -33.2333         15.1660
   Total:               4        -12.4000          4.0000
"""


class TestCaFE(unittest.TestCase):
    def test_the_delta_block_is_read_and_the_sem_is_derived(self):
        """CaFE reports no standard error, so it is sd/sqrt(n) -- and that is
        the only place in this stage where a s.e.m. is derived rather than
        taken from the file."""
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "result_1.log")
            with open(p, "w") as fh:
                fh.write(_CAFE.format(delta=_CAFE_DELTA))
            out = B.parse_cafe(p)
            mean, sd, sem = out["models"]["PB"]["TOTAL"]
            self.assertAlmostEqual(mean, -12.4, 6)
            self.assertAlmostEqual(sd, 4.0, 6)
            self.assertAlmostEqual(sem, 2.0, 6)             # 4.0 / sqrt(4)
            self.assertEqual(out["n_frames"], 4)

    def test_a_log_with_no_delta_block_yields_nothing(self):
        """The real attempt in this workspace died in APBS: its logs hold the
        gas-phase MM term and no Delta block at all. Reading the Complex block
        instead would put -47,051 kcal/mol in a binding-energy column."""
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "result_1.log")
            with open(p, "w") as fh:
                fh.write(_CAFE.format(delta=""))
            self.assertEqual(B.parse_cafe(p), {})


class TestNamelists(unittest.TestCase):
    def test_the_method_parameters_are_read_from_the_input(self):
        with tempfile.TemporaryDirectory() as d:
            write_in(os.path.join(d, "mmgbsa.in"), igb=5, saltcon=0.1)
            params, conflicts, files = B.method_parameters(d)
            self.assertEqual(params["gb.igb"], "5")
            self.assertEqual(params["gb.saltcon"], "0.100")
            self.assertEqual(params["decomp.idecomp"], "2")
            self.assertEqual(conflicts, {})
            self.assertEqual(len(files), 1)
            # A per-protomer residue list is not a method setting, and two
            # protomers with different print_res are not two methods.
            self.assertNotIn("decomp.print_res", params)

    def test_inputs_that_disagree_about_the_solvent_model_are_flagged(self):
        """Two protomers run at different `igb` are two methods wearing one
        name. Pooling them would produce a mean of two incomparable numbers,
        and nothing in either file says so."""
        with tempfile.TemporaryDirectory() as d:
            write_in(os.path.join(d, "mmgbsa_A.in"), igb=2)
            write_in(os.path.join(d, "mmgbsa_B.in"), igb=5)
            _, conflicts, _ = B.method_parameters(d)
            self.assertIn("gb.igb", conflicts)
            self.assertEqual(conflicts["gb.igb"], ["2", "5"])


class TestLocate(unittest.TestCase):
    def test_the_replica_is_matched_against_the_spec_not_counted(self):
        p = "/x/MMGBSA/A/3-R3/FINAL_MMGBSA.dat"
        self.assertEqual(B.locate(p, ["1-R1", "2-R2", "3-R3"]), ("A", "3-R3"))

    def test_a_path_matching_no_md_id_is_reported_as_unplaced(self):
        """Better an unlabelled result than one labelled with the name of
        whatever directory it happened to sit in."""
        p = "/x/MMGBSA/scratch/old/FINAL_MMGBSA.dat"
        self.assertEqual(B.locate(p, ["1-R1", "2-R2"]), (None, None))


class TestDiscover(unittest.TestCase):
    def test_nine_result_files_are_one_finding_not_nine(self):
        """Grouping by directory would pool each file with itself and report
        nine aggregates of n=1, each with no spread and no s.e.m."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            mmgbsa_tree(d, {(p, r): -5.0 for p in "ABC"
                            for r in ("1-R1", "2-R2", "3-R3")})
            found, unreadable = B.discover(spec)
            self.assertEqual(len(found), 1)
            self.assertEqual(found[0]["method"], "mmpbsa.py")
            self.assertEqual(len(found[0]["results"]), 9)
            self.assertFalse(found[0]["attempt"])
            self.assertEqual(unreadable, [])

    def test_a_working_directory_with_no_result_is_an_attempt(self):
        """A CaFE run whose PB step died leaves the gas-phase MM logs and no
        Delta block. Ingesting the MM term would put a five-figure gas-phase
        energy in a column of binding free energies."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            work = os.path.join(d, "0-MD", "ANALYSIS", "CAFE")
            os.makedirs(work)
            for name in ("com_mm.log", "rec_mm.log", "lig_mm.log", "run.tcl"):
                with open(os.path.join(work, name), "w") as fh:
                    fh.write("# frame Elec Vdw\n")
            found, _ = B.discover(spec)
            self.assertEqual([f["method"] for f in found], ["cafe"])
            self.assertTrue(found[0]["attempt"])
            self.assertEqual(found[0]["results"], [])

            rows = B.run(spec, os.path.join(d, "out"), write=False)
            self.assertEqual(rows_of(rows, observable="binding_energy"), [])
            flag = rows_of(rows, observable="method_found")
            self.assertEqual(len(flag), 1)
            self.assertFalse(flag[0]["resolved"])
            self.assertIn("produced no free energy", flag[0]["floor_note"])


class TestRun(unittest.TestCase):
    def test_a_missing_method_yields_no_rows_and_no_exception(self):
        """Three of the four systems here have no end-point calculation. An
        absent method is not an error, and must not stop the run."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            out = os.path.join(d, "out")
            rows = B.run(spec, out, write=True)
            self.assertEqual(len(rows), 0)
            tsv = os.path.join(out, f"{spec.label}__{B.STAGE}.tsv")
            self.assertTrue(os.path.exists(tsv))
            with open(tsv) as fh:
                self.assertEqual(len(fh.read().splitlines()), 1)   # header only

    def test_every_run_row_carries_the_method_and_the_file_it_came_from(self):
        """`-7.6 kcal/mol` with no method is not a result, it is a rumour."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            mmgbsa_tree(d, {("A", "1-R1"): -5.0})
            rows = B.run(spec, os.path.join(d, "out"), write=False)
            run_rows = rows_of(rows, scope="run", subset="DELTA_TOTAL")
            self.assertEqual(len(run_rows), 1)
            r = run_rows[0]
            self.assertEqual(r["method"], "mmpbsa.py")
            self.assertEqual(r["source"], "A/1-R1/FINAL_MMGBSA.dat")
            self.assertEqual(r["unit"], "kcal/mol")
            self.assertEqual(r["protomer"], "A")
            self.assertEqual(r["replica"], "1-R1")
            self.assertEqual(r["label"], "GB")
            # The stage reads no trajectory, so it has no stride and no frame
            # spacing it can honestly claim.
            self.assertIsNone(r["stride"])
            self.assertIsNone(r["ns_per_frame"])

    def test_the_sem_is_over_replicas_and_the_sd_over_all_nine(self):
        """The three protomers of one trimer share a box and a thermostat, so
        n=9 is three independent samples wearing nine hats. Taking the s.e.m.
        as sd/sqrt(9) shrinks the error bar by 1.7x for nothing."""
        vals = {("A", "1-R1"): -18.6110, ("A", "2-R2"): -13.3519,
                ("A", "3-R3"): -5.0019,
                ("B", "1-R1"): -12.6304, ("B", "2-R2"): -14.5076,
                ("B", "3-R3"): -9.1760,
                ("C", "1-R1"): 4.2026, ("C", "2-R2"): -4.0141,
                ("C", "3-R3"): 4.4388}
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            mmgbsa_tree(d, vals)
            rows = B.run(spec, os.path.join(d, "out"), write=False)
            agg = rows_of(rows, scope="aggregate", subset="DELTA_TOTAL")
            self.assertEqual(len(agg), 1)
            a = agg[0]

            flat = list(vals.values())
            mean = sum(flat) / 9
            pooled_sd = math.sqrt(sum((v - mean) ** 2 for v in flat) / 8)
            rep_means = [sum(vals[(p, r)] for p in "ABC") / 3
                         for r in ("1-R1", "2-R2", "3-R3")]
            rm = sum(rep_means) / 3
            rep_sd = math.sqrt(sum((v - rm) ** 2 for v in rep_means) / 2)

            self.assertEqual(a["n"], 9)
            self.assertAlmostEqual(a["value"], mean, 6)
            self.assertAlmostEqual(a["sd"], pooled_sd, 6)
            self.assertAlmostEqual(a["sem"], rep_sd / math.sqrt(3), 6)
            # The whole point: NOT the naive pooled s.e.m.
            self.assertNotAlmostEqual(a["sem"], pooled_sd / 3.0, 3)
            self.assertIn("pseudoreplicate", a["param_value"])

    def test_a_single_replica_reports_no_replica_level_error(self):
        """Three protomers of one box give a spread, not a standard error."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d, md_ids=("1-R1",))
            mmgbsa_tree(d, {("A", "1-R1"): -5.0, ("B", "1-R1"): -9.0,
                            ("C", "1-R1"): -1.0}, md_ids=("1-R1",))
            rows = B.run(spec, os.path.join(d, "out"), write=False)
            a = rows_of(rows, scope="aggregate", subset="DELTA_TOTAL")[0]
            self.assertEqual(a["n"], 3)
            self.assertIsNone(a["sem"])
            self.assertIn("one replica only", a["floor_note"])

    def test_a_missing_delta_total_falls_back_and_the_row_says_so(self):
        """A truncated result is recoverable from PER_FRAME.csv -- but the row
        built that way must not look like one taken from the summary, because
        its s.e.m. treats consecutive frames as independent."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d, md_ids=("1-R1",))
            root = mmgbsa_tree(d, {}, md_ids=("1-R1",))
            run = os.path.join(root, "A", "1-R1")
            write_final(os.path.join(run, "FINAL_MMGBSA.dat"),
                        with_delta=False)
            write_per_frame(os.path.join(run, "PER_FRAME.csv"),
                            [-4.0, -6.0, -5.0, -5.0])
            rows = B.run(spec, os.path.join(d, "out"), write=False)
            r = rows_of(rows, scope="run", subset="DELTA_TOTAL")[0]
            self.assertAlmostEqual(r["value"], -5.0, 6)
            self.assertFalse(r["resolved"])
            self.assertIn("PER_FRAME", r["floor_note"])

    def test_a_truncated_file_with_no_fallback_yields_an_unresolved_row(self):
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d, md_ids=("1-R1",))
            root = mmgbsa_tree(d, {}, md_ids=("1-R1",))
            write_final(os.path.join(root, "A", "1-R1", "FINAL_MMGBSA.dat"),
                        with_delta=False)
            rows = B.run(spec, os.path.join(d, "out"), write=False)
            r = rows_of(rows, scope="run", subset="DELTA_TOTAL")[0]
            self.assertIsNone(r["value"])
            self.assertFalse(r["resolved"])
            self.assertEqual(rows_of(rows, scope="aggregate"), [])

    def test_a_result_that_belongs_to_no_replica_is_not_pooled(self):
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d, md_ids=("1-R1",))
            root = mmgbsa_tree(d, {("A", "1-R1"): -5.0}, md_ids=("1-R1",))
            write_final(os.path.join(root, "scrap", "old",
                                     "FINAL_MMGBSA.dat"), total=-999.0)
            rows = B.run(spec, os.path.join(d, "out"), write=False)
            a = rows_of(rows, scope="aggregate", subset="DELTA_TOTAL")[0]
            self.assertEqual(a["n"], 1)
            self.assertAlmostEqual(a["value"], -5.0, 6)     # -999 not pooled
            unplaced = rows_of(rows, observable="unplaced_result")
            self.assertEqual(len(unplaced), 1)
            self.assertFalse(unplaced[0]["resolved"])

    def test_the_caveats_are_in_the_data_not_only_in_the_docstring(self):
        """A reader who converts one of these into a Kd must have had to
        scroll past a row that says not to."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            mmgbsa_tree(d, {(p, r): -5.0 - i for i, (p, r) in enumerate(
                [(p, r) for p in "ABC" for r in ("1-R1", "2-R2", "3-R3")])})
            rows = B.run(spec, os.path.join(d, "out"), write=False)

            note = rows_of(rows, observable="interpretation")
            self.assertEqual(len(note), 1)
            self.assertFalse(note[0]["resolved"])
            for word in ("entropy", "Kd", "RANK", "igb"):
                self.assertIn(word, note[0]["floor_note"])

            # the method settings actually used, as their own rows
            params = {r["param_name"]: r["param_value"]
                      for r in rows_of(rows, observable="method_parameter")}
            self.assertEqual(params["gb.igb"], "2")
            self.assertEqual(params["gb.saltcon"], "0.150")

            # and what the within-run spread means for pooling
            spread = rows_of(rows, observable="spread_warning")
            self.assertEqual(len(spread), 1)
            self.assertFalse(spread[0]["resolved"])
            self.assertIn("LARGER THAN THE MEAN", spread[0]["floor_note"])
            self.assertIn("pseudoreplicated", spread[0]["floor_note"])

    def test_a_search_that_could_not_look_is_not_a_search_that_found_nothing(
            self):
        """An unreadable directory and an absent calculation are different
        facts, and only one of them means there is nothing to ingest."""
        with tempfile.TemporaryDirectory() as d:
            spec = make_spec(d)
            blocked = os.path.join(d, "0-MD", "ANALYSIS")
            os.makedirs(blocked)
            os.chmod(blocked, 0o000)
            try:
                rows = B.run(spec, os.path.join(d, "out"), write=False)
                bad = rows_of(rows, observable="search_incomplete")
                if os.geteuid() == 0:                   # root reads anything
                    self.skipTest("running as root: chmod 000 is not a block")
                self.assertEqual(len(bad), 1)
                self.assertFalse(bad[0]["resolved"])
            finally:
                os.chmod(blocked, 0o755)


if __name__ == "__main__":
    unittest.main(verbosity=2)
