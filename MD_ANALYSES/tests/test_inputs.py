#!/usr/bin/env python3
"""
test_inputs.py -- choosing which files an analysis actually reads.

    python3 tests/test_inputs.py

Every stage downstream trusts this choice, and every way of getting it wrong
is silent: the analysis still runs, still fills its tables, and still draws a
page. The three that have actually happened here:

  * a topology found by NAME. `solvated_ions.psf` sits beside the stripped
    trajectory in every replica, and it has 712k atoms against the
    trajectory's 45k. The pairing is made on the ATOM COUNT, which the DCD
    header states, so it cannot be got wrong.

  * a header read too shallowly. psfgen writes one REMARKS line per patch and
    the CDR topologies carry 896 of them; a reader that gave up after 80 lines
    concluded three of the four systems had no solute topology and quietly
    analysed the SOLVATED trajectory instead -- four times the data, the wrong
    atoms, no error.

  * a trajectory found by NAME. The occluded-open replicas call theirs
    `STRIP/OCCOPEN_1R_stripped.dcd`, not `merged.dcd`, so a glob for one
    lab member's filename found nothing and fell back to the raw segments.

The fixtures are written by this file. No lab data is read.
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import analyse                                             # noqa: E402
from test_core import write_dcd                            # noqa: E402


def write_psf(path, natom, remarks=1):
    """A PSF header with `remarks` title lines before !NATOM, and no atoms.

    The atom rows are not needed: everything tested here reads the header.
    """
    lines = ["PSF EXT CMAP", "", f"{remarks:>10} !NTITLE"]
    lines += [f" REMARKS patch {i} applied to the molecule" for i in range(remarks)]
    lines += ["", f"{natom:>10} !NATOM"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


class NatomFromAVerboseHeader(unittest.TestCase):
    """The count must survive however chatty psfgen was."""

    def test_a_short_header(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_psf(os.path.join(d, "merged.psf"), 45210, remarks=1)
            self.assertEqual(analyse.psf_natoms(p), 45210)

    def test_eight_hundred_and_ninety_six_remarks(self):
        # The real CDR topologies. 80 lines was not enough, and being "not
        # enough" looked exactly like "this system has no solute topology".
        with tempfile.TemporaryDirectory() as d:
            p = write_psf(os.path.join(d, "merged.psf"), 55818, remarks=896)
            self.assertEqual(analyse.psf_natoms(p), 55818)

    def test_a_header_longer_than_the_budget_is_not_guessed_at(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_psf(os.path.join(d, "merged.psf"), 42, remarks=200)
            self.assertIsNone(analyse.psf_natoms(p, budget=64))

    def test_a_file_that_is_not_a_psf(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "notes.txt")
            with open(p, "w") as fh:
                fh.write("no NATOM record here\n")
            self.assertIsNone(analyse.psf_natoms(p))
            self.assertIsNone(analyse.psf_natoms(os.path.join(d, "absent.psf")))


class TopologyMatchedByAtomCount(unittest.TestCase):
    """The solute topology, wherever the person who built it put it."""

    def layout(self, d):
        """`MDs/2-REPLICA/1R` with the solute topology two levels up."""
        run = os.path.join(d, "MDs", "2-REPLICA", "1R")
        os.makedirs(os.path.join(run, "STRIP"))
        os.makedirs(os.path.join(d, "MDs", "1-INPUTs"))
        write_psf(os.path.join(run, "solvated_ions.psf"), 711846)
        write_psf(os.path.join(d, "MDs", "1-INPUTs", "merged.psf"), 45210)
        return run

    def test_a_topology_two_levels_up_is_found(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.layout(d)
            got = analyse.solute_topology(run, 45210)
            self.assertIsNotNone(got)
            self.assertEqual(os.path.basename(got), "merged.psf")
            self.assertIn("1-INPUTs", got)

    def test_the_solvated_topology_beside_it_is_not_taken(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.layout(d)
            # ask for an atom count nothing matches: better nothing than the
            # 712k-atom file sitting in the same directory
            self.assertIsNone(analyse.solute_topology(run, 99999))

    def test_the_replicas_own_copy_wins_when_there_is_one(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.layout(d)
            write_psf(os.path.join(run, "merged.psf"), 45210)
            self.assertEqual(analyse.solute_topology(run, 45210),
                             os.path.join(run, "merged.psf"))

    def test_it_does_not_reach_into_another_project(self):
        """A coincidental atom count in a sibling project is not a match.

        `ABZYME/` holds every project in the lab. A replica that sits directly
        in its own project folder is two levels below that, so an unbounded
        search would glob every other system's topologies -- and pairing a
        trajectory with a foreign topology of the same size raises nothing at
        all: every stage runs and every number is wrong.
        """
        with tempfile.TemporaryDirectory() as d:
            lab = os.path.join(d, "LAB")
            mine = os.path.join(lab, "MY-PROJECT", "MDs")
            other = os.path.join(lab, "SOMEONE-ELSES")
            os.makedirs(mine)
            os.makedirs(other)
            write_psf(os.path.join(other, "merged.psf"), 55632)
            # unbounded, the foreign topology is found
            self.assertIsNotNone(analyse.solute_topology(mine, 55632))
            # bounded by the folder being analysed, it is not
            self.assertIsNone(analyse.solute_topology(mine, 55632, root=mine))
            # and a topology inside the project still is
            write_psf(os.path.join(lab, "MY-PROJECT", "merged.psf"), 55632)
            self.assertIsNotNone(analyse.solute_topology(mine, 55632, root=mine))

    def test_the_bound_still_reaches_a_sibling_of_the_analysed_folder(self):
        """`MDs/1-INPUTs/merged.psf` for `MDs/2-REPLICA/1R` must still work.

        One solute topology serves three replicas and is kept beside them, not
        inside them -- the bound has to allow that or it breaks the system it
        was written for.
        """
        with tempfile.TemporaryDirectory() as d:
            run = os.path.join(d, "MDs", "2-REPLICA", "1R")
            os.makedirs(run)
            os.makedirs(os.path.join(d, "MDs", "1-INPUTs"))
            write_psf(os.path.join(d, "MDs", "1-INPUTs", "merged.psf"), 45210)
            root = os.path.join(d, "MDs", "2-REPLICA")
            self.assertIsNotNone(analyse.solute_topology(run, 45210, root=root))

    def test_no_atom_count_means_no_guess(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(analyse.solute_topology(self.layout(d), None))


class StrippedTrajectoryByContent(unittest.TestCase):
    """Any .dcd under STRIP*, chosen by how much of the run it covers."""

    def build(self, d, names_and_frames, natom=45210):
        run = os.path.join(d, "MDs", "2-REPLICA", "1R")
        os.makedirs(os.path.join(run, "STRIP"), exist_ok=True)
        os.makedirs(os.path.join(d, "MDs", "1-INPUTs"), exist_ok=True)
        write_psf(os.path.join(d, "MDs", "1-INPUTs", "merged.psf"), natom)
        for name, n in names_and_frames:
            write_dcd(os.path.join(run, "STRIP", name), n_frames=n, natom=natom)
        return run

    def test_a_file_not_called_merged_is_still_a_trajectory(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.build(d, [("OCCOPEN_1R_stripped.dcd", 10080)])
            top, dry = analyse.stripped_for(run, n_expected=10080)
            self.assertEqual(os.path.basename(dry), "OCCOPEN_1R_stripped.dcd")
            self.assertEqual(analyse.psf_natoms(top), 45210)

    def test_the_candidate_that_covers_the_run_wins_over_the_newer_one(self):
        # A replica can hold several generations of strip; the newest is not
        # the one that covers the whole run, and picking by date analysed a
        # third of a trajectory and said nothing about it.
        with tempfile.TemporaryDirectory() as d:
            run = self.build(d, [("merged.dcd", 2525)])
            os.makedirs(os.path.join(run, "STRIP_2026-09-23"))
            write_dcd(os.path.join(run, "STRIP_2026-09-23", "merged.dcd"),
                      n_frames=1293, natom=45210)
            _, dry = analyse.stripped_for(run, n_expected=2525)
            self.assertIn("STRIP/", dry.replace(os.sep, "/"))

    def test_no_stripped_pair_means_none_not_a_mismatched_pair(self):
        with tempfile.TemporaryDirectory() as d:
            run = os.path.join(d, "1R")
            os.makedirs(run)
            self.assertEqual(analyse.stripped_for(run), (None, None))

    def test_a_trajectory_whose_topology_is_missing_is_not_used(self):
        with tempfile.TemporaryDirectory() as d:
            run = os.path.join(d, "MDs", "2-REPLICA", "1R")
            os.makedirs(os.path.join(run, "STRIP"))
            write_dcd(os.path.join(run, "STRIP", "stripped.dcd"),
                      n_frames=100, natom=45210)
            # nothing on disk has 45,210 atoms: pairing it with whatever is
            # nearest would be worse than falling back to the raw segments
            self.assertEqual(analyse.stripped_for(run), (None, None))


if __name__ == "__main__":
    unittest.main(verbosity=2)
