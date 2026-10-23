#!/usr/bin/env python3
"""
test_slides.py -- the ways a picture gets filed under the wrong caption.

    python3 tests/test_slides.py

Every deck here is written by this file into a temporary directory: a .pptx is
a zip of XML, so a fixture can carry the EXACT geometry of the real slides
without carrying the renders. The two geometries below are measured from
`Presentation_NDP_test1.pptx`, because both of its shapes broke an earlier,
reasonable-looking rule:

  slide 1   two renders that OVERLAP horizontally (x 33..674 and 599..1356)
            with "Top view" at x=95 and "Side view" at x=636 -- so "the label
            inside this picture's x-span" matched both labels to the first
            picture, and one render was published under the other's name.

  slide 2   a render starting at y = -2, ABOVE the slide's own title at y=16
            -- so "a heading is the text above every picture" found nothing,
            and the slide title "1F4 tri-bound/Occluded-Open" was used as the
            left render's caption.

Neither failure raises. Both produce a complete page with a mislabelled
figure, which is why they are tested rather than eyeballed.
"""

import os
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import import_slides as IS                                 # noqa: E402

EMU = IS.EMU


def sp(text, x, y, kind=None):
    """A text box at a position, optionally a placeholder of some kind."""
    ph = f'<p:ph type="{kind}"/>' if kind else ""
    return (f'<p:sp><p:nvSpPr><p:nvPr>{ph}</p:nvPr></p:nvSpPr>'
            f'<p:spPr><a:xfrm><a:off x="{x * EMU}" y="{y * EMU}"/></a:xfrm>'
            f'</p:spPr><p:txBody><a:p><a:r><a:t>{text}</a:t></a:r></a:p>'
            f'</p:txBody></p:sp>')


def pic(rid, x, y, w, h):
    return (f'<p:pic><p:blipFill><a:blip r:embed="{rid}"/></p:blipFill>'
            f'<p:spPr><a:xfrm><a:off x="{x * EMU}" y="{y * EMU}"/>'
            f'<a:ext cx="{w * EMU}" cy="{h * EMU}"/></a:xfrm></p:spPr></p:pic>')


def deck(slides, width=1280, height=720):
    """A .pptx on disk: {slide number: (shapes xml, {rId: media name})}."""
    path = os.path.join(tempfile.mkdtemp(), "deck.pptx")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("ppt/presentation.xml",
                   f'<p:presentation><p:sldSz cx="{width * EMU}" '
                   f'cy="{height * EMU}"/></p:presentation>')
        for n, (shapes, media) in slides.items():
            z.writestr(f"ppt/slides/slide{n}.xml",
                       f"<p:sld><p:cSld><p:spTree>{shapes}"
                       f"</p:spTree></p:cSld></p:sld>")
            rels = "".join(
                f'<Relationship Id="{rid}" Target="../media/{name}"/>'
                for rid, name in media.items())
            z.writestr(f"ppt/slides/_rels/slide{n}.xml.rels",
                       f"<Relationships>{rels}</Relationships>")
            for name in media.values():
                z.writestr(f"ppt/media/{name}", b"\x89PNG\r\n\x1a\n" + b"0" * 64)
    return path


#: Slide 1 of the real deck, to the pixel.
SLIDE_1 = (pic("rId2", 33, 98, 641, 583) + pic("rId3", 599, 106, 757, 583)
           + sp("Top view", 95, 110) + sp("Side view", 636, 110)
           + sp("OPEN STATE", 479, 35) + sp("416-433 Epitope", 8, 655)
           + sp("1", 992, 681, kind="sldNum"),
           {"rId2": "image1.png", "rId3": "image2.png"})

#: Slide 2 of the real deck: a picture above the title, and a stacked legend.
SLIDE_2 = (pic("rId3", -38, 107, 835, 593) + pic("rId4", 576, -2, 760, 720)
           + sp("1F4 tri-bound/Occluded-Open", 10, 16)
           + sp("Heavy chain", 546, 625) + sp("Light chain", 546, 664),
           {"rId3": "image3.png", "rId4": "image4.png"})


class OverlappingPictures(unittest.TestCase):
    """Two renders side by side, their x-spans overlapping."""

    def setUp(self):
        self.z = zipfile.ZipFile(deck({1: SLIDE_1}))
        self.parts = IS.slide_parts(self.z, 1)

    def test_each_label_lands_on_its_own_picture(self):
        got = {p["rid"]: p["caption"] for p in self.parts["pics"]}
        self.assertEqual(got["rId2"], "Top view")
        self.assertEqual(got["rId3"], "Side view")

    def test_the_filenames_follow_the_captions(self):
        names = [IS._label(p["caption"], i, len(self.parts["pics"]))
                 for i, p in enumerate(self.parts["pics"])]
        self.assertEqual(names, ["top-view", "side-view"])

    def test_the_heading_titles_the_slide_and_captions_nothing(self):
        # "OPEN STATE" sits in the top band. It names the slide; captioning a
        # render with it is how `open-state.png` came to exist.
        self.assertEqual(self.parts["title"], "OPEN STATE")
        self.assertNotIn("OPEN STATE",
                         [p["caption"] for p in self.parts["pics"]])

    def test_the_slide_number_is_not_a_title(self):
        # Every entry an earlier version wrote was titled "1".
        self.assertNotEqual(self.parts["title"], "1")


class TitleAboveAPicture(unittest.TestCase):
    """A render that starts above the title, and a legend shared by two."""

    def setUp(self):
        self.z = zipfile.ZipFile(deck({1: SLIDE_2}))
        self.parts = IS.slide_parts(self.z, 1)

    def test_the_title_is_read_from_the_top_band_not_from_the_pictures(self):
        self.assertEqual(self.parts["title"], "1F4 tri-bound/Occluded-Open")
        self.assertEqual([p["caption"] for p in self.parts["pics"]], ["", ""])

    def test_stacked_labels_are_one_legend_not_two_captions(self):
        self.assertEqual(self.parts["legend"], ["Heavy chain", "Light chain"])

    def test_uncaptioned_renders_are_named_by_where_they_sit(self):
        names = [IS._label(p["caption"], i, len(self.parts["pics"]))
                 for i, p in enumerate(self.parts["pics"])]
        self.assertEqual(names, ["view-left", "view-right"])

    def test_the_state_in_the_title_is_recognised(self):
        self.assertEqual(IS.state_of([self.parts["title"]]), "occluded-open")


class TitlePlaceholderWins(unittest.TestCase):
    """A real title placeholder beats the top-band guess."""

    def test_placeholder_is_preferred(self):
        shapes = (pic("rId2", 40, 200, 400, 300)
                  + sp("REAL TITLE", 60, 30, kind="title")
                  + sp("Decoration", 700, 10))
        z = zipfile.ZipFile(deck({1: (shapes, {"rId2": "a.png"})}))
        self.assertEqual(IS.slide_parts(z, 1)["title"], "REAL TITLE")


class SlideSize(unittest.TestCase):
    """The top band is a fraction of the slide, so the slide must be measured."""

    def test_size_is_read_from_the_deck(self):
        z = zipfile.ZipFile(deck({1: ("", {})}, width=1600, height=900))
        self.assertEqual(IS.slide_size(z), (1600, 900))

    def test_a_deck_without_a_size_falls_back_rather_than_failing(self):
        path = os.path.join(tempfile.mkdtemp(), "bare.pptx")
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("ppt/slides/slide1.xml", "<p:sld/>")
        self.assertEqual(IS.slide_size(zipfile.ZipFile(path)), (1280, 720))


class CreditsAreRewritten(unittest.TestCase):
    """Running twice must not leave two entries, or orphan the old files."""

    def test_only_this_decks_files_are_listed_for_retirement(self):
        d = tempfile.mkdtemp()
        path = os.path.join(d, "credits.yaml")
        with open(path, "w") as fh:
            fh.write("images:\n"
                     "  - file: old-name.png\n    from: deck.pptx\n"
                     "  - file: by-hand.png\n    from: a-person.png\n")
        got = IS._credited_files(path, "/somewhere/deck.pptx")
        # a picture someone put there by hand is not this importer's to delete
        self.assertEqual(got, ["old-name.png"])

    def test_a_missing_credits_file_is_not_an_error(self):
        self.assertEqual(IS._credited_files("/no/such/credits.yaml", "d.pptx"),
                         [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
