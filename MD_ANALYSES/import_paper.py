#!/usr/bin/env python3
"""
import_paper.py -- the published paper as page material, beside the measurements.

    python3 import_paper.py <folder> <paper.pdf> [--system <record>] [--dry-run]

Writes, into `<folder>/MD_RESULTS/analyses/<record>/`:

    figures_paper/fig1.png ...   the article's figures, at the resolution the
                                 PDF carries them
    paper.yaml                   the citation, the licence, the abstract, and
                                 one block per figure: its number, its title,
                                 its caption, and which section of the page it
                                 belongs in

`page_paper_story.html` is then built from this and from `data/` together
(`analyse.py --page paper`). Nothing here touches `page_story.html`.

WHY A FILE AND NOT A FLAG. The same reason `credits.yaml` exists beside the
imported renders: extraction gets most of a caption right and the boundary
wrong occasionally, and a person must be able to fix one line without re-running
anything. Every value carries how it was established -- `caption_basis`,
`ligature_basis` -- so what was read off the page and what was inferred stay
told apart, the same rule the metadata pipeline applies to a frame interval.

THE LICENCE IS CHECKED, NOT ASSUMED. Figures are written only when the article
states a licence that permits redistribution (Creative Commons BY, BY-SA, or
CC0). Otherwise the captions and the citation are still imported, `embed` is
set to false, and the page shows each figure as a link to the article instead
of a copy of it. A lab website is publication; a paper's figures are not ours
to re-host by default.

LIGATURES. A PDF printed from a browser can carry `fi` and `fl` as glyphs in
the Unicode private-use area, which every text extractor renders as nothing:
"significant" comes out "signicant" and "Source Data file" comes out "Source
Data le". The mapping is a property of the font, not a constant, so it is
INFERRED: each private-use character is expanded with every common ligature in
turn and scored against the system word list, and the expansion that turns the
most surrounding words into real words wins. `ligature_basis` says whether that
succeeded or a standard map was assumed.
"""

import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

RESULTS = "MD_RESULTS"
ANALYSES = "analyses"
FIGDIR = "figures_paper"

#: Tried in this order against the word list. Longest first: "ffi" must be
#: preferred to "ff" where both make a word, or "office" becomes "offce".
LIGATURES = ("ffi", "ffl", "fi", "fl", "ff", "ft", "st")

#: A licence that lets a lab website host a copy of the figure.
OPEN_LICENCES = [
    (re.compile(r"Creative Commons Attribution[- ]ShareAlike", re.I), "CC-BY-SA-4.0"),
    (re.compile(r"Creative Commons Attribution", re.I), "CC-BY-4.0"),
    (re.compile(r"CC0|public domain dedication", re.I), "CC0-1.0"),
]

#: Caption keyword -> the section of the page the figure argues for. Checked in
#: order; the first match wins. Everything unmatched stays with the paper, which
#: is the safe default: a figure in the wrong section is worse than a figure in
#: the paper's own.
ROUTES = [
    (re.compile(r"\b(tilt|tilting|hinge|rmsd|fluctuat|conformational landscape|"
                r"time evolution|dynamics)\b", re.I), "moves"),
    (re.compile(r"\b(workflow|assembly|modeling|modelling|embedded|"
                r"lipid bilayer|membrane patch|construct)\b", re.I), "system"),
    (re.compile(r"\b(contact|interface|epitope|paratope|binding)\b", re.I), "touch"),
    (re.compile(r"\b(water|hydration|solvent|ion density)\b", re.I), "solvent"),
]

SECTIONS = ("system", "moves", "touch", "solvent", "paper")


# --------------------------------------------------------------------------
# text
# --------------------------------------------------------------------------
def _words():
    for p in ("/usr/share/dict/words", "/usr/share/dict/linux.words",
              "/usr/dict/words"):
        if os.path.exists(p):
            try:
                with open(p, encoding="utf-8", errors="ignore") as fh:
                    return {w.strip().lower() for w in fh if len(w.strip()) > 2}
            except OSError:                                # pragma: no cover
                continue
    return set()


def infer_ligatures(text, vocabulary=None):
    """{private-use char: the letters it stands for}, and how that was decided.

    Scored rather than assumed. A PUA character is expanded with each candidate
    in turn inside the words it appears in, and the candidate that produces the
    most dictionary words wins -- provided it explains more than half of them,
    because a mapping that fixes two words in fifty has not been established.
    """
    vocab = _words() if vocabulary is None else vocabulary
    pua = sorted({c for c in text if 0xE000 <= ord(c) <= 0xF8FF})
    if not pua:
        return {}, "no ligature glyphs in the text"
    if not vocab:
        return ({c: lig for c, lig in zip(pua, ("fi", "fl", "ff", "ffi", "ffl"))},
                "assumed: no system word list to score against")
    out, scored = {}, True
    for ch in pua:
        samples = re.findall(r"\S*" + re.escape(ch) + r"\S*", text)[:400]
        samples = [re.sub(r"[^\w" + re.escape(ch) + "]", "", s).lower()
                   for s in samples]
        samples = [s for s in samples if s]
        best, best_hits = None, 0
        for lig in LIGATURES:
            hits = sum(1 for s in samples if s.replace(ch, lig) in vocab)
            if hits > best_hits:
                best, best_hits = lig, hits
        if best and samples and best_hits > len(samples) / 2:
            out[ch] = best
        else:
            scored = False
    if not out:
        return ({c: lig for c, lig in zip(pua, ("fi", "fl", "ff", "ffi", "ffl"))},
                "assumed: no candidate explained the words it appears in")
    basis = ("scored against the system word list" if scored else
             "scored against the system word list; some glyphs unresolved")
    return out, basis


def repair(text, ligatures):
    for ch, lig in (ligatures or {}).items():
        text = text.replace(ch, lig)
    # Anything still in the private-use area would render as a blank box.
    return re.sub(r"[-]", "", text)


# --------------------------------------------------------------------------
# the article
# --------------------------------------------------------------------------
FURNITURE = re.compile(r"^(Article|Nature Communications \|.*|https?://doi\.org/.*"
                       r"|[a-z]|\d{1,3}|1234567890.*)$")
FIGSTART = re.compile(r"^Fig\.?\s*(\d+)\s*\|\s*(.+)$")
DOI = re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")


def read_pdf(path):
    """(document, pages, lines, page per line, ligature map, how it was found).

    TWO READERS, ONE FOR EACH JOB. The text comes from `pdftotext`, whose
    reading order is what the caption rules below were written against: it
    follows a two-column page the way a person does, so a caption arrives as a
    run of lines rather than interleaved with the body text beside it. The
    images come from PyMuPDF, which pdftotext cannot give.

    Both see the same ligature glyphs, and neither knows what they mean: a PDF
    printed from a browser carries `fi` and `fl` in the private-use area, so
    "significant" arrives as "signi<U+E000>cant" and renders as "signicant" in
    anything that drops unknown glyphs. The mapping is inferred (see
    `infer_ligatures`) and applied here, once, to everything downstream.
    """
    import subprocess
    import pymupdf
    try:
        raw = subprocess.run(["pdftotext", path, "-"], check=True,
                             capture_output=True).stdout.decode("utf-8", "replace")
    except FileNotFoundError:                              # pragma: no cover
        raise SystemExit("pdftotext is not on PATH (poppler-utils). It reads a "
                         "two-column page in the order a person does, which is "
                         "what the caption boundaries depend on.")
    except subprocess.CalledProcessError as exc:           # pragma: no cover
        raise SystemExit(f"pdftotext failed on {path}: "
                         f"{exc.stderr.decode('utf-8', 'replace')[:300]}")
    ligs, basis = infer_ligatures(raw)
    # `\f` is pdftotext's page separator, and the only way to know which page a
    # caption is on -- which is how its figure is found.
    pages = [repair(t, ligs) for t in raw.split("\f")]
    lines, page_of = [], []
    for n, text in enumerate(pages):
        for line in text.split("\n"):
            lines.append(line.rstrip())
            page_of.append(n)
    return pymupdf.open(path), pages, lines, page_of, ligs, basis


def find_doi(pages, filename=""):
    """The article's own DOI: the one repeated on every page, not a reference.

    A reference list holds dozens of DOIs and each appears once. The running
    header carries the article's on every page, so frequency decides it -- and
    the file name is checked against the winner rather than parsed, because a
    downloaded name is a hint and the page furniture is evidence.
    """
    from collections import Counter
    hits = Counter(d.rstrip(".,;)") for d in DOI.findall("\n".join(pages)))
    if not hits:
        return None, "no DOI found in the text"
    doi, n = hits.most_common(1)[0]
    stem = re.sub(r"[^0-9a-z]", "", os.path.basename(filename).lower())
    tail = re.sub(r"[^0-9a-z]", "", doi.split("/")[-1].lower())
    agrees = tail and tail in stem
    return doi, (f"most frequent of {len(hits)} DOIs in the text, on {n} pages"
                 + ("; the file name agrees" if agrees else ""))


def find_title(lines):
    """The lines between the header DOI and the submission dates."""
    for i, line in enumerate(lines[:80]):
        if DOI.search(line):
            out = []
            for line2 in lines[i + 1:i + 12]:
                text = line2.strip()
                if re.match(r"^(Received|Accepted|Published|Check for)", text):
                    break
                # The by-line, which follows the title with nothing between
                # them: an author carries a superscript affiliation, so the
                # line holds a letter followed by a digit, or a digit pair.
                # A title does not. Without this the title ran on through
                # every author and into the abstract.
                if re.search(r"[A-Za-z]\d|\d,\s*\d", text):
                    break
                if text:
                    out.append(text)
                if len(out) >= 4:
                    break
            if out:
                return " ".join(out)
    return None


def find_journal(lines):
    for line in lines[:400]:
        m = re.match(r"^(.+?)\s*\|\s*\((\d{4})\)(\d+):(\d+)\s*$", line.strip())
        if m:
            return {"journal": m.group(1), "year": int(m.group(2)),
                    "volume": int(m.group(3)), "article": m.group(4)}
    return {}


def find_licence(text):
    """The article's licence, matched across line breaks.

    "licensed under a Creative Commons\nAttribution 4.0 International License"
    is one sentence set over two lines, and a pattern that does not collapse
    the whitespace first finds nothing -- which reads as "no licence" and is
    the difference between hosting a figure and refusing to.
    """
    flat = re.sub(r"\s+", " ", text)
    for pattern, spdx in OPEN_LICENCES:
        if pattern.search(flat):
            return spdx, True
    return None, False


def find_abstract(pages, title):
    """Page one's abstract, found by the sentence the Introduction repeats.

    Nature articles open the Introduction by restating the abstract's first
    words. So the abstract begins where a phrase that occurs twice in the
    paper occurs the FIRST time, and ends where that phrase occurs again.
    Where the trick does not apply, no abstract is written rather than a
    wrong one.
    """
    page1 = pages[0] if pages else ""
    body = "\n".join(pages[:3])
    for start in re.finditer(r"\n([A-Z][^\n]{40,120})\n", page1):
        probe = start.group(1)[:55]
        if body.count(probe) < 2:
            continue
        rest = page1[start.start(1):]
        second = rest.find(probe, 10)
        text = rest[:second] if second > 0 else rest
        text = re.sub(r"\s+", " ", text).strip()
        if 300 < len(text) < 3000:
            return text, "the paragraph the Introduction restates"
    return None, "not found"


def captions(lines, page_of):
    """One entry per `Fig. N |`, with where its text was cut and why.

    THE BOUNDARY IS THE HARD PART. A caption runs across column and page
    breaks, so a blank line cannot end it -- inside Fig. 1 there are three.
    Two signals do end one, and both are checked:

      * the sentence journals add at the end, "Source data are provided ...";
      * a blank line reached when the text so far ENDS A SENTENCE. Mid-caption
        breaks land mid-sentence, because a column break falls where it falls;
        the break after the caption lands after a full stop.

    Where neither fires the caption is left as collected and marked `run-on`,
    which is a thing to fix in paper.yaml, not a thing to hide.
    """
    starts = [i for i, line in enumerate(lines) if FIGSTART.match(line.strip())]
    out = []
    for k, i in enumerate(starts):
        stop = starts[k + 1] if k + 1 < len(starts) else len(lines)
        buf, basis = [], "run-on: no end marker found"
        for j in range(i, stop):
            line = lines[j].strip()
            if not line:
                joined = " ".join(buf).rstrip()
                if len(joined) > 200 and joined.endswith("."):
                    basis = "blank line after a complete sentence"
                    break
                continue
            if FURNITURE.match(line):
                continue
            buf.append(line)
        text = re.sub(r"\s+", " ", " ".join(buf)).strip()
        marker = re.search(r"Source data are provided[^.]*\.", text)
        if marker:
            text, basis = text[:marker.end()], "the journal's source-data sentence"
        if basis.startswith("run-on"):
            cut = text.rfind(". ")
            if cut > 200:
                text = text[:cut + 1]
        m = FIGSTART.match(lines[i].strip())
        head = re.sub(r"^Fig\.?\s*\d+\s*\|\s*", "", text).strip()
        # The caption's own title: up to the first full stop, which is what
        # journals set in bold and what a page wants as a heading.
        end = head.find(". ")
        title = head[:end + 1].strip() if 0 < end < 200 else ""
        out.append({"figure": int(m.group(1)), "page": page_of[i],
                    "title": title, "caption": head,
                    "caption_basis": basis})
    return out


def route(cap):
    for pattern, section in ROUTES:
        if pattern.search(cap.get("title") or cap.get("caption") or ""):
            return section
    return "paper"


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def extract_figures(doc, caps, out_dir, min_pixels=250000, write=True, dpi=200):
    """One image per caption, taken from the page that caption is on.

    Matched BY PAGE, never by order: an order match silently shifts every
    figure by one when a page carries a logo the filter did not exclude, and
    nothing about the result looks wrong. A page holding no image large enough
    yields no file, and the figure keeps its caption with `file: null`.

    THE PAGE REGION IS RENDERED, NOT THE STORED BITMAP. A PDF may hold an image
    at a different horizontal and vertical resolution and let the page geometry
    undo it: five of the six figures of the article this was written for are
    stored at 300 x 150 dpi, so their pixel arrays are squashed to half height
    and only the on-page rectangle says so. Saving the bitmap gives a figure
    twice as wide as it should be, which looks like a styling problem and is
    not one. So the figure's rectangle on the page is re-rendered at a single
    dpi in both directions. It costs a little more and it also keeps the vector
    labels drawn ON TOP of the bitmap, which the stored array does not contain
    at all. Where the rectangle cannot be found, the bitmap is saved and the
    figure records that it was.

    WHY 200 dpi. The page sets a figure across at most 980 CSS pixels, and a
    480 pt rectangle at 200 dpi is 1,333 -- sharp at full width, and 10 MB for
    six figures where 300 dpi is 22 MB on a page a visitor scrolls through.
    `--dpi 300` is there for a print-quality page.
    """
    import pymupdf
    if write:
        os.makedirs(out_dir, exist_ok=True)
    made, taken = [], set()
    # A page that carries another figure's caption is that figure's page, and
    # its image is not up for grabs. Without this, a figure whose own image is
    # missing quietly takes its neighbour's, and the page then shows Fig. 5's
    # picture under Fig. 4's caption -- which no amount of reading the output
    # would reveal.
    claimed = {c["page"] for c in caps}
    for cap in caps:
        # The caption's own page first, then the page before it: a figure at
        # the foot of a page pushes its caption over the break, which is where
        # Fig. 2 of this article lives. Never by order of appearance -- one
        # logo the size filter lets through shifts every figure by one, and
        # nothing about the result looks wrong.
        best, best_px, found_on = None, 0, None
        for page_no in (cap["page"], cap["page"] - 1, cap["page"] + 1):
            if not 0 <= page_no < doc.page_count:
                continue
            if page_no != cap["page"] and page_no in claimed:
                continue
            for info in doc[page_no].get_images(full=True):
                xref, w, h = info[0], info[2], info[3]
                if xref in taken or w * h <= best_px:
                    continue
                best, best_px, found_on = xref, w * h, page_no
            if best and best_px >= min_pixels:
                break
        if not best or best_px < min_pixels:
            cap["file"] = None
            cap["file_note"] = (f"no unclaimed image of at least {min_pixels:,} "
                                f"pixels on or beside page {cap['page'] + 1}")
            continue
        taken.add(best)
        if found_on != cap["page"]:
            cap["file_note"] = (f"image taken from page {found_on + 1}; the "
                                f"caption is on page {cap['page'] + 1}")
        name = f"fig{cap['figure']}.png"
        rects = doc[found_on].get_image_rects(best)
        rect = rects[0] if rects else None
        cap["render"] = (f"page {found_on + 1} region re-rendered at {dpi} dpi"
                         if rect else "the stored bitmap, at its own resolution")
        if write:
            if rect:
                pix = doc[found_on].get_pixmap(clip=rect, dpi=dpi)
            else:
                pix = pymupdf.Pixmap(doc, best)
            if pix.n - pix.alpha >= 4:                     # CMYK to RGB
                pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
            pix.save(os.path.join(out_dir, name))
            made.append(name)
            cap["pixels"] = f"{pix.width}x{pix.height}"
        cap["file"] = name
    return made


def retire(out_dir, keep):
    """Files a previous import left behind that this one did not write.

    Same rule as the slide importer: a second run must leave the directory
    holding what the paper holds, not the union of every run.
    """
    gone = []
    if not os.path.isdir(out_dir):
        return gone
    for name in sorted(os.listdir(out_dir)):
        if name.lower().endswith(".png") and name not in keep:
            os.remove(os.path.join(out_dir, name))
            gone.append(name)
    return gone


# --------------------------------------------------------------------------
# paper.yaml
# --------------------------------------------------------------------------
HEAD = """\
# paper.yaml -- the published article that goes with this dataset.
#
# Written by import_paper.py, edited by you. `page_paper_story.html` is built
# from this and from the measured rows in data/ together; page_story.html does
# not read it.
#
#   section:  where each figure goes on the page --
#             system | moves | touch | solvent | paper
#             A section that has no measured content on this page falls back to
#             `paper`, so a figure can never be routed into nothing.
#   caption:  as extracted. `caption_basis` says where the text was cut; a
#             `run-on` one is the one to check.
#   embed:    false keeps the figures out of the page and links to the article
#             instead. Set automatically from the licence.
"""


def write_yaml(path, doc):
    import yaml
    with open(path, "w") as fh:
        fh.write(HEAD)
        yaml.safe_dump(doc, fh, sort_keys=False, allow_unicode=True,
                       width=88, default_flow_style=False)
    return path


def bundle_dir(folder, record):
    return os.path.join(folder, RESULTS, ANALYSES, record)


def records_in(folder):
    d = os.path.join(folder, RESULTS, ANALYSES)
    return sorted(n for n in os.listdir(d)) if os.path.isdir(d) else []


def run(folder, pdf, record=None, write=True, dpi=200):
    records = records_in(folder)
    if not records:
        raise SystemExit(f"no analysis bundles under {folder}/{RESULTS}/{ANALYSES} "
                         f"-- run analyse.py first")
    if record and record not in records:
        raise SystemExit(f"no bundle '{record}'. Found: {', '.join(records)}")
    if not record:
        if len(records) > 1:
            raise SystemExit(f"{len(records)} bundles here; name one with "
                             f"--system: {', '.join(records)}")
        record = records[0]

    doc, pages, lines, page_of, ligs, lig_basis = read_pdf(pdf)
    whole = "\n".join(pages)
    doi, doi_basis = find_doi(pages, pdf)
    licence, may_embed = find_licence(whole)
    abstract, abs_basis = find_abstract(pages, None)
    caps = captions(lines, page_of)
    out = bundle_dir(folder, record)
    figdir = os.path.join(out, FIGDIR)
    made = extract_figures(doc, caps, figdir, write=write and may_embed,
                           dpi=dpi)
    gone = retire(figdir, set(made)) if (write and may_embed) else []

    for cap in caps:
        cap["section"] = route(cap)

    paper = {
        "title": find_title(lines),
        "doi": doi,
        "doi_basis": doi_basis,
        "url": f"https://doi.org/{doi}" if doi else None,
        # Left empty on purpose. The record beside this file already carries an
        # author line a person confirmed, and the page uses it. Parsing the
        # by-line off a two-column PDF produces a list that is subtly wrong --
        # a missing comma merges two people -- and a wrong author list is the
        # one error on a published page nobody forgives.
        "authors": None,
        "licence": licence,
        "licence_basis": ("the article's own licence statement" if licence
                          else "no open licence statement found in the text"),
        "embed": bool(may_embed),
        "source_pdf": os.path.abspath(pdf),
        "ligatures": {f"U+{ord(c):04X}": v for c, v in ligs.items()} or None,
        "ligature_basis": lig_basis,
        "abstract": abstract,
        "abstract_basis": abs_basis,
        "figures": [{k: c.get(k) for k in
                     ("figure", "file", "section", "title", "caption",
                      "caption_basis", "file_note", "render", "pixels")
                     if c.get(k) is not None} for c in caps],
    }
    paper.update(find_journal(lines))
    if write:
        os.makedirs(out, exist_ok=True)
        write_yaml(os.path.join(out, "paper.yaml"), paper)
    return record, paper, made, gone


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("pdf")
    ap.add_argument("--system", help="the record to attach the paper to")
    ap.add_argument("--dry-run", action="store_true",
                    help="say what would be written, write nothing")
    ap.add_argument("--dpi", type=int, default=200,
                    help="resolution the figure regions are re-rendered at "
                         "(default 200; 300 for print quality and twice the "
                         "page weight)")
    a = ap.parse_args(argv)

    record, paper, made, gone = run(a.folder, a.pdf, a.system,
                                    write=not a.dry_run, dpi=a.dpi)
    print(f"\n  {record}")
    print(f"    paper    : {paper.get('title')}")
    print(f"    doi      : {paper.get('doi')}  ({paper.get('doi_basis')})")
    print(f"    licence  : {paper.get('licence') or 'NOT ESTABLISHED'} "
          f"-- figures {'embedded' if paper['embed'] else 'linked, not copied'}")
    if paper.get("ligatures"):
        pretty = ", ".join(f"{k}={v}" for k, v in paper["ligatures"].items())
        print(f"    ligatures: {pretty}  ({paper['ligature_basis']})")
    print(f"    abstract : {len(paper['abstract'] or '')} chars "
          f"({paper['abstract_basis']})")
    print(f"    figures  : {len(paper['figures'])} caption(s), {len(made)} image(s)")
    for fig in paper["figures"]:
        flag = "" if "run-on" not in fig["caption_basis"] else "   <- CHECK"
        print(f"       Fig {fig['figure']:>2} -> {fig['section']:<7} "
              f"{(fig.get('file') or 'no image'):<10} "
              f"{len(fig['caption']):>5} chars  {fig['caption_basis']}{flag}")
    if gone:
        print(f"    retired  : {', '.join(gone)}")
    if a.dry_run:
        print("\n  --dry-run: nothing was written.")
    else:
        print(f"\n  Wrote {bundle_dir(a.folder, record)}/paper.yaml")
        print("  Edit `section:` and any caption, then build the page:")
        print(f"    analyse.py {a.folder} --system {record} --pages-only --page paper")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
