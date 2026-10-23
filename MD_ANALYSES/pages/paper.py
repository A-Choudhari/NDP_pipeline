#!/usr/bin/env python3
"""
paper.py -- the published article, as blocks the story page can place.

`page_paper_story.html` is `page_story.html` with the paper's figures set
among the measured ones. This module holds everything that knows about the
article, so the story page keeps one job -- arranging sections -- and gains
one optional argument rather than a second copy of itself.

WHAT MAKES THIS HONEST. The page carries two kinds of picture, and a reader
must never have to guess which is which. A figure from the article is tagged
with its number, its source and its licence, and links to the article. A plot
from the analysis is tagged with the rows it was drawn from. That distinction
is the whole reason the paper's material is allowed on the page at all: a
website page that mixes a published figure into a measurement without saying
so is claiming the paper's result as its own output.

THE LICENCE DECIDES, NOT THE AUTHOR. `import_paper.py` sets `embed` from the
article's own licence statement. Where it is false -- no open licence found --
the caption and the citation are still shown and the figure is a link. Nothing
here overrides that: a flag that says "publish it anyway" is a flag somebody
sets by accident.
"""

import html
import os

CSS = """
/* ---- the paper's own material, told apart from the measured kind -------- */
.paper-hd{display:flex;flex-wrap:wrap;gap:8px 14px;align-items:baseline;
          margin:0 0 10px}
.paper-hd .j{color:var(--muted);font-size:14px}
a.doi{color:var(--protein);text-decoration:none;font-size:13.5px;
      border-bottom:1px solid transparent;word-break:break-all}
a.doi:hover{border-bottom-color:var(--protein)}
.abstract{font-size:15px;line-height:1.66}
.abstract:first-letter{font-size:15px}
.srctag{display:inline-block;font-size:11px;font-weight:700;letter-spacing:.09em;
        text-transform:uppercase;color:var(--paper);border:1px solid var(--paper);
        border-radius:999px;padding:1px 9px;margin:0 7px 0 0;vertical-align:2px;
        white-space:nowrap}
.srctag.measured{color:var(--protein);border-color:var(--protein)}
figure.paperfig{margin:0 0 26px}
/* Wide figures fill the measure; a PORTRAIT one is held to most of the
   viewport instead, or Fig. 2 of a paper arrives 1,150 px tall and a reader
   meets one panel of it at a time. Both constraints, and the browser honours
   whichever binds. */
figure.paperfig img{display:block;margin:0 auto;border-radius:9px;
                    background:#fff;border:1px solid var(--line);
                    max-width:100%;max-height:82vh;width:auto;height:auto}
figure.paperfig figcaption{color:var(--ink);font-size:14px;margin-top:10px;
                           line-height:1.6}
figure.paperfig figcaption .t{font-weight:660}
.paperlink{background:var(--card);border:1px dashed var(--line);
           border-radius:11px;padding:13px 15px;margin:0 0 10px}
.licence{color:var(--muted);font-size:12.5px;margin:10px 0 0}
"""

#: The accent the paper's material is tagged in, so a reader learns one colour
#: for "this came from the article" and never has to read the tag twice.
ACCENT = {"light": "#7b5cd6", "dark": "#a48bff"}

#: Where a figure may be routed. `paper` is the section of its own, and is
#: also where anything routed at a section this page does not have ends up.
SECTIONS = ("system", "moves", "touch", "solvent", "paper")


def read(bundle_dir):
    """`paper.yaml` beside the bundle, or {} when there is none.

    A missing file is the ordinary case -- most datasets have no paper yet --
    so it is not an error and not a warning. A file that cannot be parsed IS
    reported, because silence there would drop the article from the page and
    look identical to not having one.
    """
    path = os.path.join(bundle_dir, "paper.yaml")
    if not os.path.exists(path):
        return {}
    try:
        import yaml
    except ImportError:                                    # pragma: no cover
        print("        paper.yaml found but PyYAML is not installed; "
              "the paper is left off the page")
        return {}
    try:
        with open(path) as fh:
            doc = yaml.safe_load(fh) or {}
    except Exception as exc:                               # noqa: BLE001
        print(f"        paper.yaml could not be read ({exc}); "
              f"the paper is left off the page")
        return {}
    if not isinstance(doc, dict):
        return {}
    doc["_dir"] = bundle_dir
    return doc


def figures(paper, section, live_sections=()):
    """The article's figures routed to one section, in figure order.

    A figure routed at a section this page does not draw would vanish -- the
    apo trimer has no interface section, so a figure sent to `touch` would be
    silently dropped. They fall back to the paper's own section instead, which
    is why `live_sections` is passed in rather than assumed.
    """
    out = []
    for fig in paper.get("figures") or []:
        want = str(fig.get("section") or "paper")
        if want not in live_sections:
            want = "paper"
        if want == section:
            out.append(fig)
    return sorted(out, key=lambda f: f.get("figure") or 0)


def _caption(fig, paper):
    """The caption, with its own title in bold and its provenance in front."""
    text = fig.get("caption") or ""
    title = fig.get("title") or ""
    rest = text[len(title):].strip() if title and text.startswith(title) else text
    n = fig.get("figure")
    url = paper.get("url")
    label = f"Fig. {n}, from the paper" if n else "From the paper"
    tag = (f'<a class="srctag" href="{html.escape(url)}" target="_blank" '
           f'rel="noopener">{html.escape(label)}</a>' if url else
           f'<span class="srctag">{html.escape(label)}</span>')
    body = (f'<span class="t">{html.escape(title)}</span> ' if title else "")
    return f'<figcaption>{tag}{body}{html.escape(rest)}</figcaption>'


def figure_html(fig, paper):
    """One figure: the image where the licence allows it, a link where not."""
    name = fig.get("file")
    if paper.get("embed") and name:
        alt = fig.get("title") or f"Figure {fig.get('figure')}"
        return (f'<figure class="paperfig">'
                f'<img src="figures_paper/{html.escape(name)}" '
                f'alt="{html.escape(alt)}" loading="lazy">'
                f'{_caption(fig, paper)}</figure>')
    # No licence to re-host it, or no image was found on the page.
    why = fig.get("file_note") or (
        "the article states no licence that permits a copy here"
        if not paper.get("embed") else "no image was found for this figure")
    link = ""
    if paper.get("url"):
        link = (f' <a class="doi" href="{html.escape(paper["url"])}" '
                f'target="_blank" rel="noopener">See it in the article</a>')
    return (f'<div class="paperlink"><figure class="paperfig">'
            f'{_caption(fig, paper)}</figure>'
            f'<p class="licence">Figure not reproduced here: '
            f'{html.escape(why)}.{link}</p></div>')


def section_html(paper, section, live_sections=()):
    figs = figures(paper, section, live_sections)
    return "".join(figure_html(f, paper) for f in figs)


def citation(paper, authors=None):
    """The card that says whose work the paper's figures are."""
    if not paper:
        return ""
    bits = []
    who = authors or paper.get("authors")
    if who:
        bits.append(f'<p class="j">{html.escape(str(who))}</p>')
    where = ", ".join(str(x) for x in [
        paper.get("journal"),
        (f"{paper['volume']}:{paper['article']}"
         if paper.get("volume") and paper.get("article") else None),
        paper.get("year")] if x)
    head = ""
    if where:
        head += f'<span class="j">{html.escape(where)}</span>'
    if paper.get("url"):
        head += (f'<a class="doi" href="{html.escape(paper["url"])}" '
                 f'target="_blank" rel="noopener">'
                 f'{html.escape(paper["url"])}</a>')
    lic = paper.get("licence")
    note = ("Figures marked <b>from the paper</b> are the article's own, "
            f"reproduced under {html.escape(lic)} with attribution. "
            "Everything else on this page was measured from the trajectories "
            "and is drawn from the rows in <code>data/</code>."
            if lic and paper.get("embed") else
            "The article's figures are not reproduced here. Everything on "
            "this page was measured from the trajectories.")
    body = (f'<div class="card"><div class="paper-hd">{head}</div>'
            + "".join(bits)
            + f'<p class="licence">{note}</p></div>')
    if paper.get("abstract"):
        body += (f'<div class="card"><h3>Abstract</h3>'
                 f'<p class="abstract">{html.escape(paper["abstract"])}</p>'
                 f'</div>')
    return body
