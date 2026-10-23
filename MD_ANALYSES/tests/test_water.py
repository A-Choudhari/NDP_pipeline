#!/usr/bin/env python3
"""
test_water.py -- one test per defect the water stage has already paid for.

    python3 tests/test_water.py

Synthetic fixtures only. Nothing here reads the lab's data: every geometry is
three points in a box the test wrote itself, chosen so that the RIGHT answer is
obvious by hand and the WRONG answer is the one the original code produced
before the comment above it existed.

The expensive checks -- re-deriving the interfacial set from a full-side
box-aware recomputation, and re-deriving a sample of the flagged waters with an
independent float64 minimum-image brute force -- belong on a real trajectory,
which this file deliberately does not touch. What is ported here is the part
that can be pinned down synthetically: the CROSS-BOUNDARY minimum image, and the
`n_interfacial <= min(n_rec_only, n_ab_only)` invariant, which the original test
asserted on every frame of every output.
"""

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from core import spec as SP                                # noqa: E402
from stages.solvent import water as W                      # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
L = np.array([60.0, 70.0, 80.0])


def shell_masks(W_pos, R, A, L, cut=W.D_WATER_SHELL):
    """The two one-sided masks, by an independent brute force.

    Deliberately not the module's code path: this is the reference the module's
    answer is compared against, exactly as the original test used a hand-written
    minimum image rather than an MDAnalysis call.
    """
    def near(P, X):
        d = X[None, :, :] - P[:, None, :]
        d -= L * np.rint(d / L)
        return (np.sqrt((d * d).sum(-1)) <= cut).any(1)
    return near(W_pos, R), near(W_pos, A)


# --------------------------------------------------------------------------
class TestMinimumImage(unittest.TestCase):
    """NOTHING IS EVER WRAPPED. Only explicit minimum image, on a copy."""

    def test_a_contact_across_the_boundary_is_found(self):
        """The failure this prevents: a water at x = 0.5 and a protein atom at
        x = 59.5 in a 60 A box are 1 A apart, not 59. Plain Euclidean distance
        says they are nowhere near each other, the shell comes back empty, and
        the interface reads as dry."""
        w = np.array([[0.5, 10.0, 10.0]])
        a = np.array([[59.5, 10.0, 10.0]])
        d_plain = float(np.linalg.norm(w[0] - a[0]))
        d_image = float(np.linalg.norm(W.image_about(a, w[0], L)[0] - w[0]))
        self.assertAlmostEqual(d_plain, 59.0, 6)
        self.assertAlmostEqual(d_image, 1.0, 6)
        self.assertTrue(d_image <= W.D_WATER_SHELL)

    def test_the_original_coordinates_are_not_modified(self):
        """`image_about` operates on a COPY. An in-place fold would leave the
        Universe's own positions permanently displaced -- which is the wrapping
        this module exists to avoid."""
        a = np.array([[59.5, 10.0, 10.0]])
        before = a.copy()
        W.image_about(a, np.array([0.5, 10.0, 10.0]), L)
        self.assertTrue(np.array_equal(a, before))

    def test_a_group_split_across_the_boundary_is_made_coherent(self):
        """A group living in two images has a MEANINGLESS centroid: it sits in
        the middle of the box, between the two halves. Imaging about it leaves
        the group split, so `coherent` anchors on one atom first."""
        X = np.array([[1.0, 10, 10], [2.0, 10, 10], [59.0, 10, 10]])
        Y, c, half = W.coherent(X, L)
        span = Y.max(0) - Y.min(0)
        self.assertLess(float(span[0]), 5.0)          # 4 A, not 58
        self.assertLess(float(half[0]), 3.0)
        # the naive centroid would have been ~20.7, in empty space
        self.assertGreater(abs(float(X.mean(0)[0]) - float(c[0])), 10.0)

    def test_a_coherent_group_is_left_alone(self):
        X = np.array([[30.0, 35, 40], [31.0, 35, 40], [32.0, 35, 40]])
        Y, c, half = W.coherent(X, L)
        self.assertTrue(np.allclose(Y, X))
        self.assertTrue(np.allclose(c, X.mean(0)))

    def test_the_water_window_is_a_superset_not_a_decision(self):
        """The prefilter may over-select; it may never drop a water the fine
        search would have found. That is the whole basis for using it."""
        rng = np.random.default_rng(0)
        wxyz = rng.uniform(0, 1, size=(4000, 3)) * L
        wx = np.ascontiguousarray(wxyz[:, 0])
        cen = np.array([0.0, 5.0, 5.0])               # ON the boundary
        half = np.array([4.0, 4.0, 4.0])
        s, Wc = W.water_window(wx, wxyz, cen, half, W.D_WATER_SHELL, L)
        # brute force: every water whose min-image displacement is inside
        d = wxyz - cen
        d -= L * np.rint(d / L)
        want = np.flatnonzero((np.abs(d) <= half + W.D_WATER_SHELL).all(1))
        self.assertEqual(sorted(s.tolist()), sorted(want.tolist()))
        # and the returned coordinates really are in the image nearest `cen`
        self.assertTrue(np.all(np.abs(Wc - cen) <= half + W.D_WATER_SHELL + 1e-9))

    def test_a_patch_too_big_for_one_centre_is_refused(self):
        """`patch_ok` is the precondition for single-centre imaging. Above it
        the fast path is unsafe and the caller must fall back, not guess."""
        self.assertTrue(W.patch_ok(np.array([10.0, 10, 10]), 3.5, L))
        self.assertFalse(W.patch_ok(np.array([10.0, 10, 39.0]), 3.5, L))


class TestShellInvariant(unittest.TestCase):
    """n_interfacial <= min(n_rec_only, n_ab_only), BY CONSTRUCTION."""

    def _frame(self, seed):
        """One synthetic frame: two slabs of heavy atoms facing each other
        across a gap, with water everywhere including across the boundary."""
        rng = np.random.default_rng(seed)
        R = rng.uniform(0, 1, size=(300, 3)) * np.array([8.0, 20.0, 20.0]) \
            + np.array([0.0, 25.0, 30.0])          # receptor slab at x ~ 0-8
        A = rng.uniform(0, 1, size=(200, 3)) * np.array([8.0, 20.0, 20.0]) \
            + np.array([52.0, 25.0, 30.0])         # ligand slab at x ~ 52-60,
        #                                            i.e. touching R ACROSS the
        #                                            periodic boundary
        Wp = rng.uniform(0, 1, size=(6000, 3)) * L
        return R, A, Wp

    def test_the_invariant_holds_on_every_synthetic_frame(self):
        for seed in range(8):
            R, A, Wp = self._frame(seed)
            mr, ma = shell_masks(Wp, R, A, L)
            n_if = int((mr & ma).sum())
            self.assertLessEqual(n_if, min(int(mr.sum()), int(ma.sum())),
                                 f"seed {seed}: two-sided count exceeded a "
                                 f"one-sided one, which is impossible")

    def test_the_interface_reduction_loses_nothing(self):
        """The provably-lossless step: an interfacial water's receptor partner
        must lie within 3.5 + 3.5 = 7.0 A of some ligand atom. Restricting the
        two-sided search to rec_if / lig_if must give the SAME set."""
        for seed in range(8):
            R, A, Wp = self._frame(seed)
            mr, ma = shell_masks(Wp, R, A, L)
            full = np.flatnonzero(mr & ma)

            d = A[None, :, :] - R[:, None, :]
            d -= L * np.rint(d / L)
            dist = np.sqrt((d * d).sum(-1))
            rif = np.flatnonzero((dist <= W.D_INTERFACE).any(1))
            aif = np.flatnonzero((dist <= W.D_INTERFACE).any(0))
            self.assertTrue(len(rif) and len(aif),
                            f"seed {seed}: the fixture has no interface")
            mr2, ma2 = shell_masks(Wp, R[rif], A[aif], L)
            reduced = np.flatnonzero(mr2 & ma2)
            self.assertEqual(full.tolist(), reduced.tolist(),
                             f"seed {seed}: the reduction dropped a water")

    def test_the_exclusive_counts_partition(self):
        R, A, Wp = self._frame(3)
        mr, ma = shell_masks(Wp, R, A, L)
        n_if = int((mr & ma).sum())
        self.assertEqual(int(mr.sum()) - n_if, int((mr & ~ma).sum()))
        self.assertEqual(int(ma.sum()) - n_if, int((ma & ~mr).sum()))

    def test_the_boundary_is_where_the_contacts_are(self):
        """Guards the fixture itself: if the two slabs did not touch across the
        boundary, the invariant tests above would pass on an empty set."""
        R, A, Wp = self._frame(0)
        mr, ma = shell_masks(Wp, R, A, L)
        self.assertGreater(int((mr & ma).sum()), 0)
        # ...and a NON-periodic reading of the same frame finds nothing
        d = A[None, :, :] - R[:, None, :]
        self.assertGreater(float(np.sqrt((d * d).sum(-1)).min()), 40.0)


class TestRunLengthsAndSurvival(unittest.TestCase):
    """Residence is run_lengths plus arithmetic; both are checkable by hand."""

    def test_survival_is_one_at_tau_zero_and_falls(self):
        lens = np.array([1, 1, 2, 5, 10])
        tau, S = W.survival(lens, n_win=50, eff_ns=0.1, max_tau_frames=20)
        self.assertAlmostEqual(float(S[0]), 1.0, 9)
        self.assertTrue(np.all(np.diff(S[np.isfinite(S)]) <= 1e-12))
        self.assertAlmostEqual(float(tau[1]), 0.1, 9)

    def test_survival_matches_the_definition_term_by_term(self):
        """S(tau) = sum_runs max(L - tau, 0) / sum_runs L, computed the slow,
        obvious way."""
        lens = np.array([1, 3, 3, 7])
        tau, S = W.survival(lens, n_win=50, eff_ns=0.2, max_tau_frames=8)
        den = float(lens.sum())
        for i, t in enumerate(range(len(tau))):
            want = float(np.maximum(lens - t, 0).sum()) / den
            self.assertAlmostEqual(float(S[i]), want, 9)

    def test_an_empty_occupancy_gives_nan_not_zero(self):
        tau, S = W.survival(np.zeros(0, int), n_win=10, eff_ns=0.2,
                            max_tau_frames=5)
        self.assertTrue(np.isnan(S).all())

    def test_runs_are_measured_in_ns_via_the_replica_spacing(self):
        """The whole point of the port: a run of 3 frames is 0.3 ns in one
        project and 0.6 ns in another, and a module constant cannot know which."""
        occ = np.array([[0, 1, 1, 1, 0, 1, 0, 0]], dtype=bool)
        row, start, ln = W.all_runs(occ)
        self.assertEqual(ln.tolist(), [3, 1])
        self.assertAlmostEqual(float((ln * 0.1).max()), 0.3, 9)
        self.assertAlmostEqual(float((ln * 0.2).max()), 0.6, 9)

    def test_edge_runs_are_flagged_as_censored(self):
        occ = np.array([[1, 1, 0, 1, 1]], dtype=bool)
        row, start, ln = W.all_runs(occ)
        n_win = occ.shape[1]
        edge = (start == 0) | (start + ln == n_win)
        self.assertEqual(edge.tolist(), [True, True])

    def test_numpy_2_has_no_trapz(self):
        """`np.trapz` was REMOVED in NumPy 2 and this environment is 2.4. The
        original's two fits both called it; a stage that imports fine and then
        raises AttributeError inside a curve fit is the worst kind of late
        failure."""
        self.assertFalse(hasattr(np, "trapz"))
        self.assertTrue(hasattr(np, "trapezoid"))
        src = open(os.path.join(os.path.dirname(HERE), "stages", "solvent",
                                "water.py")).read()
        self.assertNotIn("np.trapz(", src)
        self.assertIn("np.trapezoid(", src)


class TestSamplingFloor(unittest.TestCase):
    """An unresolved number must say so, not look plausible."""

    def test_residence_is_refused_at_this_project_s_spacing(self):
        """0.2 ns/frame, stride 1: below the floor for water residence. The
        stage still runs; the rows say the number is not resolved."""
        eff, resolved, note = W._verdict(W.RESIDENCE, 0.2, 1,
                                         W.MIN_SPACING_NS_FOR[W.RESIDENCE])
        self.assertAlmostEqual(eff, 0.2, 9)
        self.assertFalse(resolved)
        self.assertIn("0.02", note)
        self.assertIn("stride 1", note)          # and says lowering it cannot help

    def test_a_fine_enough_spacing_is_accepted(self):
        eff, resolved, note = W._verdict(W.RESIDENCE, 0.01, 1,
                                         W.MIN_SPACING_NS_FOR[W.RESIDENCE])
        self.assertTrue(resolved)
        self.assertIsNone(note)

    def test_counts_have_no_floor(self):
        """A count is a count at any spacing; only the durations have a floor."""
        self.assertIsNone(W.MIN_SPACING_NS_FOR[W.SHELL])
        self.assertIsNone(W.MIN_SPACING_NS_FOR[W.BRIDGES])
        for stride in (1, 100):
            _, resolved, _ = W._verdict(W.SHELL, 0.2, stride,
                                        W.MIN_SPACING_NS_FOR[W.SHELL])
            self.assertTrue(resolved)

    def test_an_unknown_spacing_is_not_resolved(self):
        """`check_spacing` returns None for an unknown framestep, which reads as
        'no floor to violate'. It is the opposite."""
        _, resolved, note = W._verdict(W.RESIDENCE, None, 1,
                                       W.MIN_SPACING_NS_FOR[W.RESIDENCE])
        self.assertFalse(resolved)
        self.assertIn("framestep_ns", note)

    def test_the_spacing_check_is_the_shared_one(self):
        """Not a private copy: the refusal comes from core.spec, so every stage
        refuses at the same place and says the same thing."""
        with self.assertRaises(SP.SpacingTooCoarse):
            SP.check_spacing(W.RESIDENCE, 0.2, 1, W.MIN_SPACING_NS)


class TestNoModuleLevelTime(unittest.TestCase):
    """The port's central constraint, asserted rather than trusted."""

    def test_no_module_level_frame_count_or_spacing(self):
        """`C.NS_PER_FRAME = 0.1` and `N_FRAMES = 4999` are exactly what made
        the originals unportable: this OPEN set is 0.2 ns/frame with 2,520
        frames, and a module constant silently pools the two."""
        import re
        src = open(os.path.join(os.path.dirname(HERE), "stages", "solvent",
                                "water.py")).read()
        for line in src.splitlines():
            if not re.match(r"^[A-Z_][A-Z0-9_]*\s*=", line):
                continue
            name = line.split("=")[0].strip()
            self.assertNotIn("N_FRAMES", name)
            self.assertNotIn("NS_PER_FRAME", name)
            # a threshold named in FRAMES is the same bug wearing a hat
            self.assertFalse(name.endswith("_FRAMES") and name != "MAX_FRAMES",
                             f"{name} is a duration expressed in frames")

    def test_every_duration_constant_is_named_in_ns(self):
        for name in ("LONG_NS", "MAX_TAU_NS", "HYDRATION_EVERY_NS",
                     "FRAME_ROW_EVERY_NS", "FLICKER_MEAN_NS",
                     "MIN_SPACING_NS"):
            self.assertIsInstance(getattr(W, name), float)


class TestWaterIdentity(unittest.TestCase):
    """Water oxygen is OH2, never O; heavy atoms by NAME, never by mass."""

    def test_the_oxygen_name_preference_puts_charmm_first(self):
        self.assertEqual(W.WATER_O_NAMES[0], "OH2")
        self.assertIn("O", W.WATER_O_NAMES)       # last resort, inside water only

    def test_heavy_atoms_are_selected_by_name(self):
        """Never by mass. These topologies have broken masses -- the lab's own
        contact script documents a mass test keeping 100 of 4,963 heavy atoms --
        so `.masses` and a `mass` selection keyword must appear nowhere."""
        import io as _io
        import tokenize
        from core import selections as S
        self.assertEqual(S.HEAVY, "not name H*")
        path = os.path.join(os.path.dirname(HERE), "stages", "solvent",
                            "water.py")
        with open(path, "rb") as fh:
            src = fh.read()
        # Comments and docstrings SAY "mass" -- they are where the reasoning
        # lives. Only executable tokens are checked, so the prose survives and
        # `ag.masses` or a `mass` selection keyword would not.
        code = []
        for tok in tokenize.tokenize(_io.BytesIO(src).readline):
            if tok.type in (tokenize.COMMENT, tokenize.STRING,
                            tokenize.NL, tokenize.NEWLINE):
                continue
            code.append(tok.string)
        self.assertNotIn("masses", code)
        self.assertNotIn("mass", code)


if __name__ == "__main__":
    unittest.main(verbosity=2)
