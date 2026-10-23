#!/usr/bin/env python3
"""
test_paper.py -- the ways a paper could reach a page saying something untrue.

    python3 tests/test_paper.py

The combined page publishes somebody else's figures next to our measurements,
which makes a different class of mistake possible from the rest of this kit.
Every test below is one of those, and each exists because the obvious
implementation gets it wrong in a way that still renders perfectly:

  * a caption that runs on into the body text beside it -- the whole reason
    the boundary is two rules and not one;
  * a figure routed at a section this page does not draw, and so dropped;
  * a licence sentence broken across a line, read as "no licence", and a
    figure hosted that should only have been linked -- or the reverse;
  * a DOI taken from the reference list instead of the running header;
  * a figure matched to an image by order, so one stray logo shifts every
    figure by one;
  * a title that runs on through the by-line and into the abstract.

Nothing here opens the lab's data or the real article. The fixtures are
strings and a PDF built in memory.
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import import_paper as IP                                  # noqa: E402
from pages import paper as PA                              # noqa: E402


def lines_of(text):
    lines = text.split("\n")
    return lines, [0] * len(lines)


CAPTION_PAGE = """\
Fig. 1 | A title that ends here. a The first panel, described over
several lines as a column break falls where it falls, which is mid-

sentence and must not end the caption. b The second panel, and the
caption ends here.

body text resumes on the next block and must not be swept in, because it
sits in the column beside the figure and reads as if it belonged.
"""

SOURCE_DATA = """\
Fig. 2 | Another title. a A panel. b Another panel, with numbers in it.
Source data are provided as a Source Data file.

More body text that is not part of the caption at all.
"""


class Captions(unittest.TestCase):
    def test_a_blank_line_mid_sentence_does_not_end_a_caption(self):
        """A column break lands where it lands; a caption break lands after a
        full stop. Ending at the first blank line truncated every caption at
        its first column break, and the result read like a complete sentence."""
        caps = IP.captions(*lines_of(CAPTION_PAGE))
        self.assertEqual(len(caps), 1)
        text = caps[0]["caption"]
        self.assertIn("mid- sentence", text)
        self.assertTrue(text.rstrip().endswith("the caption ends here."))
        self.assertNotIn("body text resumes", text)
        self.assertIn("blank line", caps[0]["caption_basis"])

    def test_the_source_data_sentence_ends_a_caption(self):
        caps = IP.captions(*lines_of(SOURCE_DATA))
        self.assertTrue(caps[0]["caption"].endswith("a Source Data file."))
        self.assertNotIn("More body text", caps[0]["caption"])
        self.assertIn("source-data", caps[0]["caption_basis"])

    def test_the_title_is_the_first_sentence_and_the_number_is_kept(self):
        caps = IP.captions(*lines_of(CAPTION_PAGE))
        self.assertEqual(caps[0]["figure"], 1)
        self.assertEqual(caps[0]["title"], "A title that ends here.")

    def test_a_caption_with_no_end_marker_says_so(self):
        """It is left long and FLAGGED. A tool that silently guessed a
        boundary would be wrong in the one case nobody would think to check."""
        run_on = ("Fig. 9 | No marker anywhere. a A panel that simply keeps "
                  "going and going with no blank line and no source data "
                  "sentence to end it, for long enough to pass the minimum "
                  "length the rule needs before it will cut anything at all.")
        caps = IP.captions(*lines_of(run_on))
        self.assertIn("run-on", caps[0]["caption_basis"])


class Ligatures(unittest.TestCase):
    """A browser-printed PDF hides `fi` and `fl` in the private-use area.

    Every extractor drops them, so "significant" reads "signicant" on a
    published page. The mapping is a property of the font, so it is inferred
    rather than hardcoded -- and what it was inferred from is recorded.
    """

    VOCAB = {"significant", "field", "flexibility", "final", "file"}

    def test_the_mapping_is_scored_not_assumed(self):
        text = ("signicant " * 4 + "eld " * 4 + "nal " * 4
                + "exibility " * 6)
        ligs, basis = IP.infer_ligatures(text, vocabulary=self.VOCAB)
        self.assertEqual(ligs[""], "fi")
        self.assertEqual(ligs[""], "fl")
        self.assertIn("scored", basis)

    def test_repair_applies_the_map_and_drops_what_it_could_not_place(self):
        out = IP.repair("signicant ow ", {"": "fi",
                                                            "": "fl"})
        self.assertEqual(out, "significant flow ")

    def test_a_document_with_no_such_glyphs_says_so(self):
        ligs, basis = IP.infer_ligatures("ordinary text", vocabulary=self.VOCAB)
        self.assertEqual(ligs, {})
        self.assertIn("no ligature glyphs", basis)


class Metadata(unittest.TestCase):
    def test_a_licence_sentence_is_matched_across_a_line_break(self):
        """It wraps in the original. A pattern that does not flatten the
        whitespace first finds nothing, reads it as "no licence", and turns
        an embeddable figure into a link."""
        text = ("This article is licensed under a Creative Commons\n"
                "Attribution 4.0 International License, which permits use")
        spdx, may = IP.find_licence(text)
        self.assertEqual(spdx, "CC-BY-4.0")
        self.assertTrue(may)

    def test_no_licence_statement_means_no_embedding(self):
        spdx, may = IP.find_licence("All rights reserved.")
        self.assertIsNone(spdx)
        self.assertFalse(may)

    def test_the_doi_is_the_one_in_the_running_header_not_the_first_seen(self):
        """A reference list holds dozens, each once; the article's own is on
        every page. Taking the first match takes a cited paper's."""
        pages = ["Article https://doi.org/10.1000/theirs.1 body 10.1038/ours",
                 "10.1038/ours", "10.1038/ours", "10.1000/theirs.2"]
        doi, basis = IP.find_doi(pages)
        self.assertEqual(doi, "10.1038/ours")
        self.assertIn("most frequent", basis)

    def test_the_title_stops_at_the_by_line(self):
        """An author carries a superscript affiliation, so the by-line holds a
        letter against a digit. Without that the title ran through every
        author and into the abstract."""
        lines = ["Article", "https://doi.org/10.1038/x", "A Title Over",
                 "Two Lines", "Jane Doe1,3, John Roe 2", "Abstract begins"]
        self.assertEqual(IP.find_title(lines), "A Title Over Two Lines")

    def test_routing_falls_back_to_the_paper_when_nothing_matches(self):
        self.assertEqual(IP.route({"title": "Glycan shield of Env."}), "paper")
        self.assertEqual(IP.route({"title": "Env tilting motion over time."}),
                         "moves")


class Placement(unittest.TestCase):
    PAPER = {"url": "https://doi.org/10.1038/x", "embed": True,
             "licence": "CC-BY-4.0",
             "figures": [{"figure": 2, "file": "fig2.png", "section": "touch",
                          "title": "T.", "caption": "T. Body."},
                         {"figure": 1, "file": "fig1.png", "section": "moves",
                          "title": "M.", "caption": "M. Body."}]}

    def test_a_figure_routed_at_a_section_this_page_lacks_is_not_dropped(self):
        """An apo system draws no interface section. A figure sent there would
        vanish from the page with nothing to say it had."""
        live = {"moves", "paper"}
        self.assertEqual([f["figure"] for f in
                          PA.figures(self.PAPER, "touch", live)], [])
        self.assertEqual([f["figure"] for f in
                          PA.figures(self.PAPER, "paper", live)], [2])

    def test_figures_of_a_section_come_in_the_article_s_order(self):
        live = {"moves", "touch", "paper"}
        both = dict(self.PAPER)
        both["figures"] = [dict(f, section="moves") for f in self.PAPER["figures"]]
        self.assertEqual([f["figure"] for f in PA.figures(both, "moves", live)],
                         [1, 2])

    def test_an_embedded_figure_is_tagged_with_its_number_and_links_out(self):
        out = PA.figure_html(self.PAPER["figures"][1], self.PAPER)
        self.assertIn('src="figures_paper/fig1.png"', out)
        self.assertIn("Fig. 1, from the paper", out)
        self.assertIn("https://doi.org/10.1038/x", out)

    def test_without_a_licence_the_figure_is_a_link_and_not_a_copy(self):
        """The licence decides, not the caller. This is the one failure here
        that is somebody else's rights rather than our wrong number."""
        paper = dict(self.PAPER, embed=False, licence=None)
        out = PA.figure_html(paper["figures"][1], paper)
        self.assertNotIn("<img", out)
        self.assertIn("not reproduced here", out)
        self.assertIn("Fig. 1, from the paper", out)

    def test_a_caption_cannot_open_a_tag(self):
        evil = {"figure": 3, "file": "f.png", "section": "paper",
                "title": "", "caption": "<script>alert(1)</script>"}
        out = PA.figure_html(evil, self.PAPER)
        self.assertNotIn("<script>", out)
        self.assertIn("&lt;script&gt;", out)

    def test_no_paper_file_is_not_an_error(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(PA.read(tmp), {})


class Images(unittest.TestCase):
    """Matched by page, never by order of appearance."""

    def _pdf(self):
        import pymupdf
        doc = pymupdf.open()
        big = pymupdf.Pixmap(pymupdf.csRGB, pymupdf.IRect(0, 0, 700, 700), False)
        big.clear_with(200)
        small = pymupdf.Pixmap(pymupdf.csRGB, pymupdf.IRect(0, 0, 20, 20), False)
        small.clear_with(10)
        p0 = doc.new_page()
        p0.insert_image(pymupdf.Rect(10, 10, 300, 300), pixmap=small)
        p1 = doc.new_page()
        p1.insert_image(pymupdf.Rect(10, 10, 300, 300), pixmap=big)
        doc.new_page()
        return doc

    def test_a_logo_does_not_become_a_figure_or_steal_the_next_one(self):
        """Page one holds only a logo; page two holds a real figure and its
        own caption. The first caption must come away with nothing, rather
        than reaching over the break and taking the second figure's image."""
        doc = self._pdf()
        caps = [{"figure": 1, "page": 0}, {"figure": 2, "page": 1}]
        IP.extract_figures(doc, caps, "/nonexistent", write=False)
        self.assertIsNone(caps[0]["file"])
        self.assertIn("no unclaimed image", caps[0]["file_note"])
        self.assertEqual(caps[1]["file"], "fig2.png")

    def test_a_caption_that_flowed_over_the_page_break_finds_its_figure(self):
        """The figure sits at the foot of one page and its caption starts on
        the next. Looking only at the caption's own page loses it."""
        doc = self._pdf()
        caps = [{"figure": 1, "page": 2}]
        IP.extract_figures(doc, caps, "/nonexistent", write=False)
        self.assertEqual(caps[0]["file"], "fig1.png")
        self.assertIn("taken from page 2", caps[0]["file_note"])

    def test_the_figure_is_as_wide_as_the_page_shows_it_not_as_stored(self):
        """A PDF may store an image at 300 x 150 dpi and let the page geometry
        undo it. Five of the six figures of the article this was written for
        are stored that way, so saving the pixel array gives a figure twice as
        wide as it should be -- which looks like a CSS problem, and is not
        one. The rectangle on the page is the truth, and it is re-rendered."""
        import tempfile
        import pymupdf
        doc = pymupdf.open()
        square = pymupdf.Pixmap(pymupdf.csRGB, pymupdf.IRect(0, 0, 800, 800),
                                False)
        square.clear_with(180)
        page = doc.new_page()
        # A square bitmap laid into a rectangle twice as wide as it is tall:
        # the page says 2.0, the pixels say 1.0.
        page.insert_image(pymupdf.Rect(20, 20, 420, 220), pixmap=square)
        caps = [{"figure": 1, "page": 0}]
        with tempfile.TemporaryDirectory() as tmp:
            IP.extract_figures(doc, caps, tmp, write=True, dpi=100)
            out = pymupdf.open(os.path.join(tmp, "fig1.png"))[0].rect
        self.assertAlmostEqual(out.width / out.height, 2.0, delta=0.03)
        self.assertIn("re-rendered", caps[0]["render"])

    def test_one_image_is_never_given_to_two_figures(self):
        doc = self._pdf()
        caps = [{"figure": 1, "page": 1}, {"figure": 2, "page": 1}]
        IP.extract_figures(doc, caps, "/nonexistent", write=False)
        self.assertEqual(caps[0]["file"], "fig1.png")
        self.assertIsNone(caps[1]["file"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
