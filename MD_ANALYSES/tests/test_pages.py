#!/usr/bin/env python3
"""
test_pages.py -- what a published page could say that the data does not.

    python3 tests/test_pages.py

The pages are the part of this kit that leaves the lab. A wrong number in a
TSV is found by whoever reads the TSV; a wrong number in a caption is read by
everyone and checked by nobody, so the failures worth testing here are the
ones that produce a page that still looks right:

  * a caption quoting a value the plot beside it does not draw -- the reason
    captions read `rows` and take no number from their caller;
  * a duration printed in the wrong unit, so a record and a page disagree
    about how long a run was;
  * a caption quoting a quantity its own stage marked unresolved;
  * an apo system showing an interface section with nothing in it;
  * markup in a caption escaping into the page as a tag.

Every fixture below is synthetic rows built in this file. Nothing reads the
lab's data, so a test failing here is this code being wrong, never a
trajectory having changed.
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from core import io                                        # noqa: E402
from pages import captions as C                            # noqa: E402
from pages import layout as L                              # noqa: E402
from pages import series as S                              # noqa: E402
from pages import story as ST                              # noqa: E402

NA = io.NA


def row(**kw):
    """A tidy row with the mandatory columns already filled."""
    base = {"scope": "aggregate", "subset": None, "value": None, "sd": None,
            "n": None, "ns_per_frame": 0.2, "stride": 1, "resolved": True,
            "protomer": NA, "replica": NA}
    base.update(kw)
    return base


def rmsd_rows(values, sd=0.5, n=9, subset="receptor"):
    """A pooled RMSD curve, one row per 0.2 ns frame."""
    return [row(stage="structure", observable="rmsd", subset=subset,
                time_ns=i * 0.2, value=v, sd=sd, n=n)
            for i, v in enumerate(values)]


class HumanNs(unittest.TestCase):
    """Nanoseconds below a microsecond, microseconds above it -- one rule.

    The record's length tag and the page's header line are written by two
    different trees. They must agree, so the rule is tested on the exact
    values the four CDR systems produce.
    """

    def test_below_a_microsecond_stays_ns(self):
        self.assertEqual(io.human_ns(504), "504 ns")
        self.assertEqual(io.human_ns(999), "999 ns")

    def test_at_and_above_a_microsecond(self):
        self.assertEqual(io.human_ns(1000), "1 µs")
        self.assertEqual(io.human_ns(1512), "1.5 µs")
        self.assertEqual(io.human_ns(2133), "2.1 µs")

    def test_the_sign_is_micro_not_the_letter_u(self):
        # U+00B5 MICRO SIGN, which the owner asked for by name. A latin "u"
        # renders as "1.5 us" and looks like a typo in a published header.
        self.assertIn("µ", io.human_ns(1512))
        self.assertNotIn("us", io.human_ns(1512))

    def test_unit_can_be_dropped_for_a_tag(self):
        self.assertEqual(io.human_ns(1512, unit=False), "1.5")


class CaptionNumbers(unittest.TestCase):
    """A caption's number must be the one the plot draws, or it is a lie."""

    def test_rmsd_caption_quotes_the_plateau_the_band_shows(self):
        # rises over 100 frames, then flat at 4.66 for 300 more
        vals = [1 + 3.66 * i / 100 for i in range(100)] + [4.66] * 300
        rows = rmsd_rows(vals)
        panel = S.panel(rows, "rmsd", "receptor", "RMSD", "Å",
                        "protein", "time (ns)")
        panel["short"] = "RMSD"
        text = C.for_panel(rows, panel)
        self.assertIn("4.66", text)
        # and the value it quotes is the one in the band, not a slice of it
        self.assertAlmostEqual(panel["band"]["mean"][-1], 4.66, places=2)
        self.assertIn("9 protomer", text)

    def test_a_curve_still_climbing_is_not_called_settled(self):
        rows = rmsd_rows([0.01 * i for i in range(400)])
        panel = S.panel(rows, "rmsd", "receptor", "RMSD", "Å",
                        "protein", "time (ns)")
        panel["short"] = "RMSD"
        text = C.for_panel(rows, panel)
        self.assertIn("still climbing", text)
        self.assertNotIn("stays there", text)

    def test_settling_time_is_rounded_because_about_means_about(self):
        vals = [1 + 3.0 * i / 100 for i in range(100)] + [4.0] * 300
        rows = rmsd_rows(vals)
        panel = S.panel(rows, "rmsd", "receptor", "RMSD", "Å",
                        "protein", "time (ns)")
        panel["short"] = "RMSD"
        text = C.for_panel(rows, panel)
        # the sampling is 0.2 ns; a settling time written to a tenth of a
        # nanosecond claims a precision the measurement does not have
        self.assertNotIn(".8 ns", text)
        self.assertNotIn(".2 ns", text)

    def test_binding_caption_quotes_the_pooled_total(self):
        rows = [row(stage="binding_energy", observable="binding_energy",
                    subset="DELTA_TOTAL", value=-7.63, sd=8.21, sem=2.74, n=9,
                    unit="kcal/mol"),
                row(stage="binding_energy", observable="binding_energy",
                    subset="VDWAALS", value=-60.0, sd=5.0, n=9),
                row(stage="binding_energy", scope="run", protomer="A",
                    replica="1-R1", observable="binding_energy",
                    subset="DELTA_TOTAL", value=-7.0),
                row(stage="binding_energy", param_name="gb.igb",
                    param_value="2", observable="binding_energy",
                    subset="DELTA_TOTAL", value=-7.63)]
        text = C.binding(rows)
        self.assertIn("-7.6", text)
        self.assertIn("8.2", text)
        # the method is named from its own settings, not left as a generic
        self.assertIn("MM/GBSA", text)

    def test_a_spread_that_rounds_to_zero_is_not_printed(self):
        # "0.0312 +/- 0.0000" says the replicas differ by something too small
        # to show, which is better said by not saying it
        self.assertNotIn("±", C._pm(0.0312, 1e-9, "Å⁻³", 4))
        self.assertIn("±", C._pm(0.0312, 0.002, "Å⁻³", 4))


class UnresolvedIsSaidInWords(unittest.TestCase):
    """A stage that could not resolve something must not be quoted anyway."""

    def test_contact_lifetime_is_described_not_quoted(self):
        rows = [
            row(stage="interface_contacts", observable="contact_pairs",
                value=1105),
            row(stage="interface_contacts", observable="contacts_per_frame",
                value=26.6, sd=9.8),
            row(stage="interface_contacts", observable="interface_residues",
                subset="receptor", value=59),
            row(stage="interface_contacts", observable="interface_residues",
                subset="ligand", value=52),
            row(stage="interface_contacts",
                observable="contact_lifetime_single_frame_fraction",
                value=0.52),
            row(stage="interface_contacts", observable="contact_lifetime_mean",
                value=0.21, resolved=False,
                floor_note="0.2 ns per frame cannot resolve a lifetime"),
        ]
        text = C.contact_map(rows)
        self.assertIn("not measurable", text)
        self.assertIn("52%", text)
        # the unresolvable number itself is nowhere in the sentence
        self.assertNotIn("0.21", text)

    def test_solvent_caption_refuses_to_report_residence(self):
        rows = [row(stage="rdf", scope="run", observable="bulk_density",
                    subset="waterO", value=0.0312),
                row(stage="rdf", scope="run", observable="first_peak_r",
                    subset="waterO_around_protein_polar", value=2.77)]
        text = C.solvent(rows)
        self.assertIn("2.77", text)
        # The refusal is what matters, not the words it is phrased in: with
        # no dwell measurement in these rows the sentence says the quantity
        # was not measured, and either way it quotes no residence time.
        self.assertIn("NOT reported", text)
        self.assertNotIn("ns residence", text)


class SamplingIsReadNotAssumed(unittest.TestCase):
    """How fine the sampling is comes from the rows, never from memory.

    The solvent caption said "at 0.2 ns per frame" on every page. True of the
    CDR campaign, wrong by a factor of two for the occluded-open runs, which
    write a frame every 100 ps -- and wrong in the one sentence whose job is
    to explain why a number is NOT being quoted.
    """

    def dwell(self, ns_per_frame, stride, fraction):
        return [row(stage="ions", scope="run", subset="SOD",
                    observable="ion_dwell_single_frame_fraction",
                    value=fraction, ns_per_frame=ns_per_frame, stride=stride),
                row(stage="rdf", scope="run", observable="bulk_density",
                    subset="waterO", value=0.0312,
                    ns_per_frame=ns_per_frame, stride=250)]

    def test_the_spacing_quoted_is_frame_interval_times_stride(self):
        text = C.solvent(self.dwell(0.2, 10, 0.98))
        self.assertIn("2 ns between sampled frames", text)
        text = C.solvent(self.dwell(0.1, 10, 0.96))
        self.assertIn("1 ns between sampled frames", text)

    def test_the_fraction_quoted_is_the_measured_one(self):
        self.assertIn("96%", C.solvent(self.dwell(0.1, 10, 0.96)))
        self.assertIn("98%", C.solvent(self.dwell(0.2, 10, 0.98)))

    def test_no_dwell_measurement_means_no_borrowed_evidence(self):
        # the apo systems skip the water-shell stage entirely: with nothing
        # measured, the sentence must not quote another system's figure
        rows = [row(stage="rdf", scope="run", observable="bulk_density",
                    subset="waterO", value=0.0312, ns_per_frame=0.1,
                    stride=250)]
        text = C.solvent(rows)
        self.assertIn("NOT reported", text)
        self.assertNotIn("%", text)

    def test_a_system_with_no_spacing_recorded_still_reads(self):
        rows = [row(stage="rdf", scope="run", observable="bulk_density",
                    subset="waterO", value=0.0312, ns_per_frame=None,
                    stride=None)]
        self.assertIn("NOT reported", C.solvent(rows))


class StoryPageShape(unittest.TestCase):
    """The layout decisions, made from the rows rather than from a flag."""

    def test_markup_is_bolded_but_html_cannot_escape(self):
        got = ST._rich("a **1.23 Å** value <script>x</script>")
        self.assertIn("<b>1.23 Å</b>", got)
        self.assertIn("&lt;script&gt;", got)
        self.assertNotIn("<script>", got)

    def test_an_unclosed_marker_does_not_leave_a_dangling_tag(self):
        got = ST._rich("**4.66")
        self.assertEqual(got.count("<b>"), got.count("</b>"))

    def test_an_empty_section_is_not_drawn_at_all(self):
        # the apo trimer has no interface; the section must be absent, not
        # present and empty -- an empty frame with a title reads as a failure
        self.assertEqual(ST._section("What it touches", "The interface",
                                     "lede", ""), "")
        self.assertIn("The interface",
                      ST._section("What it touches", "The interface", "lede",
                                  "<div></div>"))

    def test_apo_rows_report_no_ligand(self):
        apo = rmsd_rows([4.0] * 20)
        self.assertFalse(L.has_ligand(apo))
        self.assertTrue(L.has_ligand(apo + rmsd_rows([9.0] * 20,
                                                     subset="ligand")))


class TheOpeningSentence(unittest.TestCase):
    """What surrounds the solute is read, never assumed.

    The full-length gp120 trajectories are membrane-embedded and contain no
    water at all -- it was stripped before those files were written. The page
    opened by calling them "solvated and neutralised" anyway, because the
    phrase was hardcoded: an assertion nothing had measured, on a page whose
    entire claim is that it never makes one.
    """

    def water(self, n):
        return [row(stage="composition", scope="run",
                    observable="water_molecules", value=n, replica="R1")]

    def test_a_membrane_system_with_no_water_is_not_called_solvated(self):
        got = C.system([], {}, None, ("protein", "glycans", "membrane"))
        self.assertIn("embedded in a lipid bilayer", got)
        self.assertNotIn("solvated", got)

    def test_a_solvated_system_says_so_because_the_water_was_counted(self):
        got = C.system(self.water(236838), {}, None,
                       ("protein", "glycans", "water", "ions"))
        self.assertIn("solvated and neutralised", got)
        self.assertNotIn("lipid bilayer", got)

    def test_a_bound_ligand_is_named_before_its_surroundings(self):
        got = C.system(self.water(236838), {}, None,
                       ("protein", "glycans", "water", "ions", "antibody"))
        self.assertLess(got.index("antibody bound to it"),
                        got.index("solvated"))

    def test_a_membrane_system_that_also_kept_its_water_says_both(self):
        got = C.system(self.water(100), {}, None,
                       ("protein", "glycans", "water", "membrane"))
        self.assertIn("lipid bilayer", got)
        self.assertIn("solvated", got)

    def test_no_glycans_means_no_glycosylated(self):
        self.assertNotIn("glycosylated",
                         C.system([], {}, None, ("protein", "membrane")))


class CompositionTable(unittest.TestCase):
    """The one table where three replicas carry information -- and where not."""

    def test_counts_are_integers_and_identical_values_state_no_spread(self):
        rows = [row(stage="composition", scope="run", observable="ions",
                    subset="SOD", value=668.0, replica=f"R{i}")
                for i in range(3)]
        table = dict((lbl, (v, sd)) for lbl, v, sd, n in
                     L.composition_single_column(rows))
        value, spread = table["Na⁺ ions"]
        self.assertIsInstance(value, int)          # not "668.0" sodium ions
        self.assertIsNone(spread)                  # same box in all three

    def test_a_box_that_breathes_keeps_its_spread(self):
        rows = []
        for i, edge in enumerate((203.1, 204.2, 204.5)):
            rows += [row(stage="composition", scope="run",
                         observable="box_volume", value=8.4e6 + i * 1e4,
                         replica=f"R{i}")]
        table = dict((lbl, (v, sd)) for lbl, v, sd, n in
                     L.composition_single_column(rows))
        self.assertIsNotNone(table["Box volume (Å³)"][1])




class OwnersDescription(unittest.TestCase):
    """The paragraph under the title is the owner's, or there is none.

    Every record carries a `description`; where nobody wrote one the metadata
    pipeline composes it from the atom count and the setup. Printing that at
    the top of a published page would restate the record table further down
    and read as if somebody had chosen to say it.
    """

    def test_a_written_description_appears(self):
        from pages import story
        out = story._blurb({"description": "Three replicas of the bound Fab.",
                            "description_basis": "declared"})
        self.assertIn("Three replicas of the bound Fab.", out)
        self.assertIn('class="blurb"', out)

    def test_a_composed_one_does_not(self):
        from pages import story
        out = story._blurb({"description": "Molecular dynamics simulations of "
                                           "this system (20,862 atoms).",
                            "description_basis": "composed"})
        self.assertEqual(out, "")

    def test_an_absent_description_is_not_an_error(self):
        from pages import story
        self.assertEqual(story._blurb(None), "")
        self.assertEqual(story._blurb({}), "")
        self.assertEqual(story._blurb({"description_basis": "declared"}), "")

    def test_markup_in_a_description_cannot_open_a_tag(self):
        from pages import story
        out = story._blurb({"description": "<script>x</script> and **bold**",
                            "description_basis": "declared"})
        self.assertNotIn("<script>", out)
        self.assertIn("&lt;script&gt;", out)
        self.assertIn("<b>bold</b>", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
