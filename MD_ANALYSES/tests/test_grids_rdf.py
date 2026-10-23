#!/usr/bin/env python3
"""
test_grids_rdf.py -- one test per way these two stages can be quietly wrong.

    python3 tests/test_grids_rdf.py

Synthetic fixtures only: a 60 A box holding eight solute atoms, one sodium, one
chloride 2.92 A from it, and 64 waters on an 11 A lattice, written as a
three-model PDB that BOTH cpptraj and MDAnalysis read. Nothing here touches the
lab's data, and everything it writes lives in a temporary directory.

The point is not coverage. Every failure listed below is one that returns a
plausible number instead of an error:

  * cpptraj exits 0 after doing nothing, because a mask matched no atoms;
  * a grid is binned in the units of whatever `normdensity` happened to mean,
    so an isovalue of "5" is not five times bulk;
  * "the first peak" is the tallest peak, which on a protein is the carbon
    shell at 3.7 A and not the hydrogen bond at 2.8 A;
  * a `.dx` is read back with its axes transposed or truncated;
  * an atom-number mask is built from a selection that matched nothing.

The ion pair is placed at 2.92 A so the answer is arithmetic: with dr = 0.05 A
the bins are centred on 0.025 + 0.05k, so 2.92 falls in the bin centred at
2.925 and nowhere else, and with exactly one chloride per sodium the running
coordination number past that shell is exactly 1.
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np                                         # noqa: E402

from core import cpptraj as C                              # noqa: E402
from core import spec as SP                                # noqa: E402
from stages.solvent import grids as G                      # noqa: E402
from stages.solvent import rdf as R                        # noqa: E402

BOX = 60.0
#: The ion separation the g(r) peak has to reproduce, and the bin it must
#: land in at dr = 0.05 A.
ION_SEP = 2.92
ION_BIN = 2.925
N_WATER = 64
N_FRAMES = 3

HAVE_CPPTRAJ = os.path.exists(C.SERIAL_BIN)


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def _atom(serial, name, resname, resid, segid, x, y, z, element):
    """One ATOM record in the columns the PDB standard actually specifies.

    Written by hand rather than by a library because the column layout IS the
    fixture: MDAnalysis reads `segid` from columns 73-76 and this workspace
    selects by segid everywhere, so a record that is one column out produces a
    topology whose segments are named `RO` and `T1` -- which is not an error,
    just a selection that matches nothing.
    """
    return ("ATOM  "
            f"{serial:5d} {name:<4s} {resname:>3s} A{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}      "
            f"{segid:<4s}{element:>2s}\n")


def write_toy(path, box=BOX, n_frames=N_FRAMES):
    """A tiny solvated system with every number in it known in advance.

    Eight solute atoms symmetric about the box centre (so `autoimage` has
    nothing to move and `rms first` is the identity), one Na+/Cl- pair at a
    fixed separation, and 64 waters on a lattice, each far enough from the next
    that no two share a voxel at 1 A.

    Every frame is identical, deliberately: a fixture whose answer depends on
    thermal noise cannot assert an exact number, and the arithmetic is what is
    under test here.

    WHY THE HALF-ANGSTROM OFFSETS. The solute's bounding box is 26..34 in every
    axis, so with 30 A of padding and 1 A bins the grid origin is exactly -4.0
    and every voxel boundary falls on an integer. An atom placed AT an integer
    sits on a boundary, and the rounding that decides which side of it the atom
    lands on is not stable across frames -- the identity rotation a fit returns
    is only the identity to about 1e-7. Water and ions therefore sit on voxel
    CENTRES (x.5), which is what makes "one voxel per water" an exact count
    rather than a coin toss. The solute itself stays on integers, because it is
    what fixes the origin.
    """
    c = box / 2.0
    solute = []
    # Four CA and four N, symmetric in every axis and NOT coplanar, so the
    # least-squares fit has a well-defined rotation rather than a degenerate one.
    for k, (dx, dy, dz, nm, el) in enumerate([
            (-4.0, -2.0, -1.5, "CA", "C"), (4.0, 2.0, 1.5, "CA", "C"),
            (-2.0, 4.0, 1.5, "CA", "C"), (2.0, -4.0, -1.5, "CA", "C"),
            (-3.0, 0.0, 3.0, "N", "N"), (3.0, 0.0, -3.0, "N", "N"),
            (0.0, -3.0, 2.0, "N", "N"), (0.0, 3.0, -2.0, "N", "N")]):
        solute.append((nm, "ALA", k + 1, "PRO", c + dx, c + dy, c + dz, el))
    ions = [("SOD", "SOD", 101, "ION", 20.5, 20.5, 20.5, "NA"),
            ("CLA", "CLA", 201, "ION", 20.5 + ION_SEP, 20.5, 20.5, "CL")]
    water = []
    n = 0
    side = int(round(N_WATER ** (1.0 / 3.0)))
    for a in range(side):
        for b in range(side):
            for d in range(side):
                x, y, z = 6.5 + a * 12.0, 6.5 + b * 12.0, 6.5 + d * 12.0
                n += 1
                water += [("OH2", "WAT", 300 + n, "WT1", x, y, z, "O"),
                          ("H1", "WAT", 300 + n, "WT1", x + .76, y + .59, z, "H"),
                          ("H2", "WAT", 300 + n, "WT1", x - .76, y + .59, z, "H")]
    body = []
    for i, (nm, rn, ri, sg, x, y, z, el) in enumerate(solute + ions + water, 1):
        body.append(_atom(i, nm, rn, ri, sg, x, y, z, el))

    cryst = ("CRYST1" + f"{box:9.3f}{box:9.3f}{box:9.3f}"
             + f"{90.0:7.2f}{90.0:7.2f}{90.0:7.2f}" + " "
             + "P 1".ljust(11) + "   1\n")
    with open(path, "w") as fh:
        fh.write(cryst)
        for f in range(n_frames):
            fh.write(f"MODEL     {f + 1:4d}\n")
            fh.writelines(body)
            fh.write("ENDMDL\n")
        fh.write("END\n")
    return path, n


def toy_spec(pdb):
    md = SP.MD(id="TOY", topology=pdb, trajectories=[pdb], framestep_ns=0.1,
               solvated_topology=pdb, solvated_trajectories=[pdb])
    return SP.Spec(id="test/toy", label="toy",
                   components=["protein", "water", "ions"],
                   selections={"receptor": "segid PRO"}, mds=[md])


def write_dx(path, counts, origin, spacing, values):
    """An OpenDX grid in cpptraj's own layout, for testing the reader."""
    nx, ny, nz = counts
    with open(path, "w") as fh:
        fh.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
        fh.write("origin {} {} {}\n".format(*origin))
        fh.write(f"delta {spacing} 0 0\ndelta 0 {spacing} 0\n"
                 f"delta 0 0 {spacing}\n")
        fh.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
        fh.write(f"object 3 class array type double rank 0 items "
                 f"{len(values)} data follows\n")
        for i in range(0, len(values), 3):
            fh.write(" ".join(f"{v:g}" for v in values[i:i + 3]) + "\n")
        fh.write('\nobject "density [A^-3]" class field\n')
    return path


# --------------------------------------------------------------------------
class TestMasks(unittest.TestCase):
    """cpptraj has no `segid`, so every selection becomes atom numbers."""

    def test_contiguous_indices_become_one_range(self):
        self.assertEqual(C.mask_from_indices(np.arange(0, 10)), "@1-10")

    def test_gaps_become_separate_ranges_and_singletons(self):
        self.assertEqual(C.mask_from_indices([0, 1, 2, 6, 8, 9]),
                         "@1-3,7,9-10")

    def test_an_empty_selection_is_an_error_not_an_empty_mask(self):
        """The whole point: cpptraj would take `@` and do nothing, quietly."""
        with self.assertRaises(C.CpptrajError):
            C.mask_from_indices([])

    def test_a_mask_of_236838_ranges_is_refused(self):
        """Water oxygen by index is one range per molecule. The error has to
        name the alternative (`:WAT@OH2`), or the next person writes it again."""
        with self.assertRaises(C.CpptrajError) as cm:
            C.mask_from_indices(np.arange(0, 2000, 2), max_ranges=200)
        self.assertIn("NAME", str(cm.exception))

    def test_heavy_is_by_name_never_by_mass(self):
        """These topologies carry broken masses; `!@H*` is the only safe form."""
        self.assertEqual(C.heavy("@1-10"), "(@1-10)&(!@H*)")


class TestTrajin(unittest.TestCase):
    """The stride is applied by the reader. That is why these stages exist."""

    def test_every_segment_is_strided_at_read_time(self):
        out = C.trajin_lines(["a.dcd", "b.dcd"], stride=200)
        self.assertEqual(out, ["trajin a.dcd 1 last 200",
                               "trajin b.dcd 1 last 200"])

    def test_max_frames_caps_the_first_segment_only(self):
        """A cap spread over ten segments would change WHICH part of the run
        was sampled without saying so; a probe reads the head and admits it."""
        out = C.trajin_lines(["a.dcd", "b.dcd", "c.dcd"], stride=10,
                             max_frames=3)
        self.assertEqual(out, ["trajin a.dcd 1 21 10"])


class TestLogScanning(unittest.TestCase):
    """cpptraj says most of what went wrong in warnings, and exits 0 anyway."""

    BANNER = "CPPTRAJ: Trajectory Analysis. V6.24.0 (AmberTools)\n"
    OMP = ("CPPTRAJ: Trajectory Analysis. V6.4.4 (AmberTools) OpenMP\n"
           "| 24 OpenMP threads available.\n")

    def test_the_version_comes_from_the_banner_not_the_path(self):
        errs, warns, ver, flav, thr, *_ = C._scan(self.BANNER)
        self.assertEqual((ver, flav), ("V6.24.0", "serial"))
        self.assertIsNone(thr)

    def test_the_openmp_build_is_recognised_with_its_thread_count(self):
        _, _, ver, flav, thr, *_ = C._scan(self.OMP)
        self.assertEqual((ver, flav, thr), ("V6.4.4", "OpenMP", 24))

    def test_an_empty_mask_warning_is_promoted_to_an_error(self):
        """This is the failure mode the whole driver exists for: exit 0, no
        output file, one warning nobody reads."""
        errs, warns, *_ = C._scan(
            self.BANNER + "Warning: First mask has no atoms.\n"
                          "Warning: Setup incomplete for [radial]: Skipping\n")
        self.assertEqual(len(errs), 2)
        self.assertEqual(warns, [])

    def test_the_psf_and_dcd_boilerplate_is_not_surfaced(self):
        """Warnings that appear on EVERY run of these systems train people to
        ignore the list, which is how a real one gets missed."""
        errs, warns, *_ = C._scan(
            self.BANNER +
            "Warning: Determining bond length parameters from element types "
            "for 'solvated_ions.psf'.\n"
            "Warning: Trajectory box type is 'Orthorhombic' but topology box "
            "type is 'None'.\n")
        self.assertEqual((errs, warns), ([], []))

    def test_a_warning_that_is_neither_is_kept_for_the_manifest(self):
        _, warns, *_ = C._scan(self.BANNER + "Warning: something new\n")
        self.assertEqual(warns, ["Warning: something new"])

    def test_the_progress_bar_does_not_hide_the_rest_of_the_log(self):
        """Progress is written with carriage returns; splitting on '\\n' alone
        makes the whole run one line and every scan of it match nothing."""
        errs, _, _, _, _, nread, ndone = C._scan(
            self.BANNER + " 0%\r 50%\rRead 17 frames and processed 17 frames.\n")
        self.assertEqual((nread, ndone), (17, 17))
        self.assertEqual(errs, [])


class TestDXReader(unittest.TestCase):
    def test_a_grid_reads_back_with_its_geometry_and_values(self):
        with tempfile.TemporaryDirectory() as d:
            vals = np.arange(2 * 3 * 4, dtype=float)
            p = write_dx(os.path.join(d, "g.dx"), (2, 3, 4), (1.0, 2.0, 3.0),
                         0.5, vals)
            g = C.read_dx(p)
            self.assertEqual(g.counts, (2, 3, 4))
            self.assertEqual(g.origin, (1.0, 2.0, 3.0))
            self.assertEqual(g.spacing, (0.5, 0.5, 0.5))
            self.assertAlmostEqual(g.voxel_volume, 0.125)
            np.testing.assert_allclose(g.values, vals)

    def test_a_truncated_grid_is_refused_not_averaged(self):
        """A `.dx` cut short by a full disk averages to the wrong density and
        looks exactly like a correct one."""
        with tempfile.TemporaryDirectory() as d:
            p = write_dx(os.path.join(d, "g.dx"), (2, 3, 4), (0.0, 0.0, 0.0),
                         0.5, np.zeros(24))
            text = open(p).read().replace("items 24", "items 48")
            open(p, "w").write(text)
            with self.assertRaises(C.CpptrajError):
                C.read_dx(p)


class TestFirstShell(unittest.TestCase):
    """What "the first peak" means, on curves whose answer is known."""

    @staticmethod
    def _curve(peaks, r_max=15.0, dr=0.05):
        r = np.arange(dr / 2, r_max, dr)
        g = np.zeros_like(r)
        for centre, height, width in peaks:
            g += height * np.exp(-0.5 * ((r - centre) / width) ** 2)
        return r, g

    def test_the_first_peak_is_the_first_one_not_the_tallest(self):
        """The defect this replaces. Against all protein heavy atoms the
        tallest feature inside 6 A is the carbon shell at 3.7 A; the first
        hydration shell is the hydrogen bond at 2.8 A, and an `argmax` over the
        window reports the wrong one by a whole angstrom, with no warning."""
        r, g = self._curve([(2.8, 0.6, 0.12), (3.7, 1.4, 0.25)])
        out = R.first_shell(r, g, None)
        self.assertIsNotNone(out)
        self.assertAlmostEqual(out[0], 2.825, places=2)

    def test_a_curve_that_only_rises_has_no_first_shell(self):
        """Answering "6.0 A" because that is where the window ends would be
        inventing a shell that is not there."""
        r = np.arange(0.025, 15.0, 0.05)
        self.assertIsNone(R.first_shell(r, r / 15.0, None))

    def test_a_single_sharp_shell_still_reports_its_boundary(self):
        """One peak and nothing after it: the walk out to the first minimum
        must stop where the curve flattens, not run to r_max and give up."""
        r, g = self._curve([(2.925, 40000.0, 0.03)])
        nint = np.cumsum(g) / g.sum()
        out = R.first_shell(r, g, nint)
        self.assertIsNotNone(out)
        self.assertAlmostEqual(out[0], 2.925, places=2)
        self.assertGreater(out[2], out[0])

    def test_an_empty_curve_reports_nothing_rather_than_bin_zero(self):
        r = np.arange(0.025, 15.0, 0.05)
        self.assertIsNone(R.first_shell(r, np.zeros_like(r), None))


# --------------------------------------------------------------------------
@unittest.skipUnless(HAVE_CPPTRAJ, f"no cpptraj at {C.SERIAL_BIN}")
class TestOnAToySystem(unittest.TestCase):
    """End to end, on a box whose every answer is arithmetic."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.pdb, cls.n_water = write_toy(os.path.join(cls.tmp.name, "toy.pdb"))
        cls.spec = toy_spec(cls.pdb)
        cls.data = os.path.join(cls.tmp.name, "out", "data")
        os.makedirs(cls.data, exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @staticmethod
    def _value(rows, observable, subset=None, **match):
        for r in rows.rows:
            if r.get("observable") != observable:
                continue
            if subset is not None and r.get("subset") != subset:
                continue
            if all(r.get(k) == v for k, v in match.items()):
                return r
        return None

    # ------------------------------------------------------------ the driver
    def test_a_mask_that_matches_nothing_is_an_error_not_an_empty_result(self):
        """cpptraj answers this with two warnings, no output file and exit 0.
        Verified on the real 1F4 system, 2026-09-23: a 2,257-byte grid of
        zeros and a data file that was never written."""
        work = os.path.join(self.tmp.name, "empty")
        with self.assertRaises(C.CpptrajError) as cm:
            C.run([f"parm {self.pdb}", f"trajin {self.pdb}",
                   f"radial out {os.path.join(work, 'g.dat')} 0.05 15.0 "
                   f":NOSUCH@XX :SOD volume"],
                  work, expect=[os.path.join(work, "g.dat")])
        self.assertIn("no atoms", str(cm.exception).lower())

    def test_the_version_that_ran_is_recorded_from_the_banner(self):
        work = os.path.join(self.tmp.name, "banner")
        res = C.run([f"parm {self.pdb}", f"trajin {self.pdb}",
                     f"volume out {os.path.join(work, 'v.dat')}"],
                    work, expect=[os.path.join(work, "v.dat")])
        self.assertTrue(res.version.startswith("V"))
        self.assertEqual(res.frames_processed, N_FRAMES)
        self.assertIn("version", res.as_dict())

    # --------------------------------------------------------------- the rdf
    def test_the_ion_pair_peak_lands_in_the_bin_arithmetic_demands(self):
        rows = R.run(self.spec, self.data, stride=1, write=False)
        peak = self._value(rows, "first_peak_r", "CLA_around_SOD")
        self.assertIsNotNone(peak, "no cation-anion curve was produced")
        self.assertAlmostEqual(peak["value"], ION_BIN, places=3)

    def test_one_chloride_per_sodium_integrates_to_exactly_one(self):
        """The coordination number is the running integral per mask2 atom. One
        chloride at 2.92 A from the one sodium is N = 1, and nothing else."""
        rows = R.run(self.spec, self.data, stride=1, write=False)
        row = self._value(rows, "first_shell_coordination", "CLA_around_SOD")
        self.assertAlmostEqual(row["value"], 1.0, places=6)

    def test_the_measured_bulk_density_is_the_count_over_the_cell(self):
        """Not a textbook 0.0334 A^-3 and not what the setup script was asked
        for: atoms from the topology over the volume these frames had."""
        rows = R.run(self.spec, self.data, stride=1, write=False)
        row = self._value(rows, "bulk_density", "waterO")
        self.assertAlmostEqual(row["value"], self.n_water / BOX ** 3,
                               places=10)

    def test_every_bin_of_every_curve_is_written_once(self):
        rows = R.run(self.spec, self.data, stride=1, write=False)
        bins = [r for r in rows.rows if r.get("scope") == "bin"
                and r.get("observable") == "g_r"
                and r.get("subset") == "CLA_around_SOD"]
        self.assertEqual(len(bins), int(R.DEFAULT_RMAX / R.DEFAULT_DR))
        self.assertAlmostEqual(bins[0]["r_A"], R.DEFAULT_DR / 2, places=6)

    def test_a_stripped_trajectory_is_refused_by_both_stages(self):
        """A solvent stage pointed at the solute-only file would grid an empty
        box and report a g(r) of nothing, with no error anywhere."""
        dry = SP.Spec(id="test/dry", label="dry", components=["protein"],
                      selections={"receptor": "segid PRO"},
                      mds=[SP.MD(id="D", topology=self.pdb,
                                 trajectories=[self.pdb])])
        for stage in (R, G):
            with self.assertRaises(ValueError):
                stage.run(dry, self.data, stride=1, write=False)

    # ------------------------------------------------------------- the grids
    def test_the_grid_is_written_in_multiples_of_bulk(self):
        """The claim an isovalue rests on. One sodium, present in every frame
        and never moving, puts n_frames counts in one voxel; at that voxel the
        value must be exactly 1 / (voxel volume * bulk density), which for one
        ion in a 60 A box at 1 A spacing is 216,000."""
        out = os.path.join(self.tmp.name, "gridout", "data")
        os.makedirs(out, exist_ok=True)
        rows = G.run(self.spec, out, stride=1, write=False, spacing=1.0,
                     padding=30.0)
        rho = 1.0 / BOX ** 3
        peak = self._value(rows, "grid_peak_density", "SOD")
        self.assertAlmostEqual(peak["value"], 1.0 / (1.0 * rho), places=3)
        counts = self._value(rows, "grid_peak_counts", "SOD")
        self.assertAlmostEqual(counts["value"], float(N_FRAMES), places=6)

    def test_one_voxel_per_water_stands_above_bulk(self):
        """64 waters on a lattice, none sharing a voxel: the count above any
        threshold is 64, and a transposed or mis-scaled grid gives anything
        else."""
        out = os.path.join(self.tmp.name, "gridout2", "data")
        os.makedirs(out, exist_ok=True)
        rows = G.run(self.spec, out, stride=1, write=False, spacing=1.0,
                     padding=30.0)
        # `param_value` is written as TEXT throughout these stages: one column
        # holding both a float and "394x378x310" makes pandas call it `object`,
        # and core/io.py's parquet write then fails and is swallowed, leaving
        # only the TSV. The TSV is the contract, but losing the parquet without
        # a word is not worth a saved `str()`.
        for t in ("2", "5", "10"):
            row = self._value(rows, "grid_voxels_above", "water_OH2",
                              param_value=t)
            self.assertIsNotNone(row, f"no threshold row at {t}x bulk")
            self.assertEqual(row["value"], self.n_water, f"at {t}x bulk")

    def test_a_shot_noise_grid_is_marked_unresolved_with_what_it_needs(self):
        """Three frames of one ion is not a density, and a row that does not
        say so is the failure this workspace exists to prevent."""
        out = os.path.join(self.tmp.name, "gridout3", "data")
        os.makedirs(out, exist_ok=True)
        rows = G.run(self.spec, out, stride=1, write=False, spacing=1.0,
                     padding=30.0)
        row = self._value(rows, "grid_bulk_counts_expected", "SOD")
        self.assertFalse(row["resolved"])
        self.assertIn("frames", row["floor_note"])

    def test_the_dx_written_is_the_dx_read_back(self):
        out = os.path.join(self.tmp.name, "gridout4", "data")
        os.makedirs(out, exist_ok=True)
        rows = G.run(self.spec, out, stride=1, write=False, spacing=1.0,
                     padding=30.0)
        row = self._value(rows, "grid_file_mb", "SOD")
        # `grids/` is a sibling of `data/`, beside `cpptraj/`.
        dx = os.path.join(os.path.dirname(out), "grids", row["param_value"])
        self.assertTrue(os.path.exists(dx), dx)
        g = C.read_dx(dx)
        voxels = self._value(rows, "grid_voxels")
        self.assertEqual(g.n_voxels, voxels["value"])
        self.assertEqual(g.n_voxels, g.counts[0] * g.counts[1] * g.counts[2])

    def test_water_oxygen_must_be_named_OH2(self):
        """`O` in these topologies is a protein backbone carbonyl. A water mask
        that falls back to it selects the wrong molecule and still returns a
        smooth, plausible curve."""
        import MDAnalysis as mda
        u = mda.Universe(self.pdb)
        resnames, name = G.water_oxygen(u)
        self.assertEqual(name, "OH2")
        self.assertIn("WAT", resnames)
        self.assertEqual(R.water_oxygen_sel(u).mask, ":WAT@OH2")


@unittest.skipIf(os.name == "nt", "the stand-in cpptraj is a shell script")
class TestCpptrajTimeout(unittest.TestCase):
    """A hung read (NFS) blocked a cpptraj, and the whole run, forever."""

    def test_a_hung_cpptraj_is_stopped_and_reported(self):
        import stat
        with tempfile.TemporaryDirectory() as tmp:
            fake = os.path.join(tmp, "cpptraj")
            with open(fake, "w") as fh:
                fh.write("#!/bin/sh\nsleep 30\n")
            os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
            was, env = C.SERIAL_BIN, os.environ.get("CPPTRAJ_TIMEOUT_S")
            C.SERIAL_BIN, os.environ["CPPTRAJ_TIMEOUT_S"] = fake, "1"
            try:
                with self.assertRaises(C.CpptrajError) as cm:
                    C.run(["parm x"], os.path.join(tmp, "work"))
                self.assertIn("did not finish", str(cm.exception))
            finally:
                C.SERIAL_BIN = was
                if env is None:
                    os.environ.pop("CPPTRAJ_TIMEOUT_S", None)
                else:
                    os.environ["CPPTRAJ_TIMEOUT_S"] = env


class TestWaterOxygenNames(unittest.TestCase):
    """GROMACS, TIP4P and SPC files call the water oxygen OW. The water stage
    accepted it and the grid and rdf stages refused it, so the same system
    passed one solvent stage and failed the next."""

    @staticmethod
    def _universe(oxygen):
        import MDAnalysis as mda
        names = ["CA", "O"] + [oxygen, "HW1", "HW2"] * 3
        resnames = ["ALA", "SOL", "SOL", "SOL"]
        u = mda.Universe.empty(len(names), n_residues=4,
                               atom_resindex=[0, 0] + [1] * 3 + [2] * 3 + [3] * 3,
                               trajectory=True)
        u.add_TopologyAttr("name", names)
        u.add_TopologyAttr("resname", resnames)
        u.add_TopologyAttr("segid", ["A"])
        return u

    def test_OW_is_a_water_oxygen_for_grids_and_rdf(self):
        u = self._universe("OW")
        self.assertEqual(G.water_oxygen(u), (["SOL"], "OW"))
        self.assertEqual(R.water_oxygen_sel(u).mask, ":SOL@OW")

    def test_a_bare_O_is_still_refused(self):
        u = self._universe("O")
        from core import selections as S
        with self.assertRaises(S.EmptySelection):
            G.water_oxygen(u)
        with self.assertRaises(S.EmptySelection):
            R.water_oxygen_sel(u)


if __name__ == "__main__":
    unittest.main(verbosity=2)
