#!/usr/bin/env python3
"""
story.py -- the page that goes on the website: one column, read top to bottom.

The other page (`report.py`) is a working view: two dense columns, everything
visible at once, built for someone who already knows the system. This one has
a different job. It is published, and it is met by a reader who was not in the
room -- so it opens with a picture of the thing, then takes one analysis at a
time, each with a sentence saying what it shows and what the number is, and
ends with the record that will become the NDP entry.

SAME NUMBERS, SAME CLICKS, DIFFERENT ORDER. Both pages are built from the same
tidy rows in `data/` and both use the SAME drawing code -- `report.JS`,
imported, not copied. The div ids are kept identical for that reason: one
engine draws the plots and one provenance panel answers the click, so the two
pages can differ in arrangement and never in a value. Clicking a point still
names the frame it came from -- `merged.dcd frame 2016 of 2525` -- which is
the feature this whole thing was built around.

WHAT IS *NOT* HERE. A section whose rows do not exist is not drawn at all,
rather than drawn empty: the apo trimer has no antibody, so it has no
interface and no binding energy, and its page ends up shorter instead of
pocked with holes.
"""

import html
import json
import os

from . import captions as C
from . import layout as L
from . import paper as PA
from . import report as P
from . import series as S

CSS = """
:root{
  --ink:#16181d;--muted:#6a707c;--line:#e3e6ea;--bg:#fff;--card:#f7f8fa;
  --protein:#1596b8;--ligand:#d1495b;--paper:#7b5cd6;
  --sans:"Inter","Helvetica Neue",-apple-system,BlinkMacSystemFont,"Segoe UI",
         Roboto,"Noto Sans",Arial,sans-serif;
  --mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ink:#e9ebef;--muted:#9aa1ad;--line:#2b2f36;--bg:#0f1115;--card:#171a20;
  --protein:#35b8d6;--ligand:#e8743c;--paper:#a48bff}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);
     font-size:16px;line-height:1.62;-webkit-font-smoothing:antialiased}
/* ONE column, and a measure a person can actually read along: a caption set
   across 1800px of screen is a caption nobody finishes. */
.wrap{max-width:980px;margin:0 auto;padding:0 26px 90px}
header{padding:44px 0 20px;border-bottom:1px solid var(--line);margin-bottom:8px}
h1{margin:0 0 8px;font-size:30px;font-weight:660;letter-spacing:-.022em;
   line-height:1.22}
.sub{color:var(--muted);font-size:14.5px;letter-spacing:.01em}
/* The owner's own paragraph, when they wrote one. Set at the page's
   measure and a little larger than the body, because it is the one piece
   of prose here that a person chose rather than a stage produced. */
.blurb{margin:16px 0 0;font-size:16.5px;line-height:1.6;max-width:70ch}
section{padding:34px 0 6px;border-top:1px solid var(--line)}
section:first-of-type{border-top:none}
.eyebrow{font-size:11.5px;font-weight:700;letter-spacing:.14em;
         text-transform:uppercase;color:var(--protein);margin:0 0 4px}
.eyebrow.lig{color:var(--ligand)}
h2{margin:0 0 4px;font-size:20px;font-weight:640;letter-spacing:-.012em}
.lede{color:var(--muted);font-size:15px;margin:0 0 18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:11px;
      padding:13px 15px;margin:0 0 10px}
.card.plot{padding:7px 5px 3px}
figure{margin:0}
.figs{display:grid;gap:18px;grid-template-columns:repeat(auto-fit,minmax(min(100%,300px),1fr))}
.figs img{width:100%;height:auto;border-radius:9px;display:block;
          background:var(--card)}
figcaption{color:var(--muted);font-size:13px;margin-top:8px}
figcaption b{color:var(--ink);font-weight:640}
/* A caption belongs to the plot above it, so it sits tight under the card
   and is indented to read as one block with it rather than as loose text. */
.cap{margin:-2px 0 26px;padding-left:2px;font-size:14.5px;color:var(--ink)}
.cap b{font-weight:660}
.note{color:var(--muted);font-size:13px;margin:8px 0 22px}
.legend{display:inline-block;font-size:12px;color:var(--muted);
        border:1px solid var(--line);border-radius:999px;padding:2px 11px;
        margin-top:4px}
.pair{display:grid;gap:18px;grid-template-columns:repeat(auto-fit,minmax(min(100%,380px),1fr))}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:6px 12px 6px 0;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:600;font-size:11.5px;letter-spacing:.06em;
   text-transform:uppercase}
tr:last-child td{border-bottom:none}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums lining-nums}
td.sd{color:var(--muted);font-size:13px}
td.k{color:var(--muted);width:40%}
h3{font-size:14.5px;margin:0 0 9px;font-weight:650}
.twocol{display:grid;gap:0 26px;
        grid-template-columns:repeat(auto-fit,minmax(min(100%,230px),1fr))}
.twocol table tr:last-child td{border-bottom:none}
.tag{display:inline-block;font-size:11px;letter-spacing:.05em;color:var(--muted);
     border:1px solid var(--line);border-radius:999px;padding:1px 9px}
.pending{color:#a86b00}
@media (prefers-color-scheme:dark){.pending{color:#e0a72c}}
.warn{color:#b2182b}
@media (prefers-color-scheme:dark){.warn{color:#ef8a86}}
.hint{color:var(--muted);font-size:13px;margin:8px 0 0;line-height:1.5}
#detail{position:fixed;right:18px;bottom:18px;width:370px;max-width:92vw;
        background:var(--card);border:1px solid var(--line);border-radius:11px;
        padding:14px 16px;box-shadow:0 10px 34px rgba(0,0,0,.18);display:none;
        z-index:20;max-height:72vh;overflow:auto}
#detail.on{display:block}
#detail .k{color:var(--muted);font-size:10.5px;letter-spacing:.07em;
           text-transform:uppercase;margin-top:10px;width:auto}
#detail .v{font-size:14px;word-break:break-word}
#detail .v b{font-weight:650;font-size:16px}
#detail .x{float:right;cursor:pointer;color:var(--muted);border:0;
           background:none;font-size:18px;line-height:1}
pre{background:rgba(128,128,128,.13);padding:9px 10px;border-radius:7px;
    overflow-x:auto;font-family:var(--mono);font-size:11px;line-height:1.5;
    margin:7px 0;white-space:pre-wrap;word-break:break-all}
button.copy{font:inherit;font-size:12px;padding:4px 10px;border:1px solid var(--line);
            background:var(--bg);color:var(--ink);border-radius:6px;cursor:pointer}
button.copy:hover{border-color:var(--protein);color:var(--protein)}
.click{font-size:13px;color:var(--muted);margin:2px 0 20px}
.click b{color:var(--ink);font-weight:620}
"""


def _rich(text):
    """A caption's `**number**` markup, escaped and then emboldened.

    Escape first, mark up second: a residue name or a file path in a caption
    must never be able to open a tag, and doing it in this order means it
    cannot.
    """
    out, bold = html.escape(text or ""), False
    while "**" in out:
        out = out.replace("**", "<b>" if not bold else "</b>", 1)
        bold = not bold
    return out + ("</b>" if bold else "")


def _cap(text):
    return f'<p class="cap">{_rich(text)}</p>' if text else ""


def _section(eyebrow, heading, lede, body, family="protein", tag=""):
    """One section. `tag` marks a section that holds two kinds of thing.

    On the combined page a section can carry both the article's figures and
    this pipeline's plots. The figures are tagged individually; the tag here
    says what the rest of the section is, so neither kind is the one a reader
    has to infer.
    """
    if not body.strip():
        return ""
    cls = "eyebrow lig" if family == "ligand" else "eyebrow"
    return (f'<section><p class="{cls}">{html.escape(eyebrow)}</p>'
            f'<h2>{html.escape(heading)}</h2>'
            f'<p class="lede">{tag}{html.escape(lede)}</p>{body}</section>')


# --------------------------------------------------------------------------
# THE SYSTEM
# --------------------------------------------------------------------------
def _pictures(bundle_dir):
    """The supplied renders, side by side, each under its own caption.

    Referenced by relative path rather than inlined: an <img src> resolves
    over file:// where a fetch does not, and four 1 MB renders as base64
    would add 5 MB to a page that already carries plotly.
    """
    d = os.path.join(bundle_dir, "figures_user")
    if not os.path.isdir(d):
        return "", "", ""
    names = sorted(n for n in os.listdir(d)
                   if n.lower().endswith((".png", ".jpg", ".jpeg", ".webp",
                                          ".svg", ".gif")))
    if not names:
        return "", "", ""
    meta = L.read_credits(d)
    # A caption of its own comes from the slide; where two renders shared one
    # legend ("Heavy chain / Light chain") the importer recorded that instead,
    # and it is shown ONCE beneath both rather than repeated under each.
    order = {"top-view.png": 0, "side-view.png": 1,
             "view-left.png": 0, "view-right.png": 1}
    names.sort(key=lambda n: (order.get(n, 9), n))
    figs, legends, warn = "", [], ""
    for n in names:
        m = meta.get(n, {})
        cap = m.get("caption") or os.path.splitext(n)[0].replace("-", " ").title()
        if m.get("legend"):
            legends.append(m["legend"])
        state, rec = m.get("state_in_title"), m.get("records_state")
        if state and rec and state != rec and not warn:
            warn = (f'<p class="hint warn">The slide these came from is titled '
                    f'&lsquo;{html.escape(state)}&rsquo;; these records are '
                    f'&lsquo;{html.escape(rec)}&rsquo;.</p>')
        # A render the deck did not caption gets no caption line: "View Left"
        # under a picture tells a reader strictly less than the picture does,
        # and the shared legend beneath both says what there is to say.
        own = bool(m.get("caption"))
        figs += (f'<figure><img src="figures_user/{html.escape(n)}" '
                 f'alt="{html.escape(cap)}">'
                 + (f'<figcaption><b>{html.escape(cap)}</b></figcaption>'
                    if own else "")
                 + '</figure>')
    tail = ""
    uniq = sorted(set(legends))
    if uniq:
        tail = "".join(f'<span class="legend">{html.escape(g)}</span> '
                       for g in uniq)
    return f'<div class="figs">{figs}</div>', tail, warn


def _composition(rows):
    table = L.composition_single_column(rows)
    if not table:
        return ""

    def half(items):
        body = "".join(
            f"<tr><td class='k'>{html.escape(lbl)}</td>"
            f"<td class='num'>{P._fmt(v)}</td></tr>"
            for lbl, v, sd, n in items)
        return f"<table>{body}</table>"

    mid = (len(table) + 1) // 2
    return (f'<div class="card"><h3>What is in the box</h3>'
            f'<div class="twocol">{half(table[:mid])}{half(table[mid:])}</div>'
            f'</div>')


def _sd_like(value, sd):
    """A spread written to the same precision as the value it belongs to.

    `4.557e-06` beside `0.0312` is arithmetically right and unreadable, and
    a spread that vanishes at the value's own precision is better left out
    than printed as a row of zeros.
    """
    if not sd or value is None:
        return ""
    dp = max(0, 4 - len(str(int(abs(value))))) if abs(value) >= 1 else 4
    return "" if round(sd, dp) == 0 else f"\u00b1 {sd:,.{dp}f}"


def _binding_table(table):
    """The table alone: the explanation belongs to the caption under it.

    `report._binding_card` carries its own paragraph about pooling and SEM,
    which is right on a page with no captions. Here the caption says it, so
    printing both would tell a reader the same thing twice in a row -- the
    caveats a stage raised about its OWN numbers are kept, because those are
    not the caption's to carry.
    """
    if not table:
        return "", ""
    BOLD = " style='font-weight:650'"
    body = "".join(
        f"<tr{BOLD if r.get('total') else ''}>"
        f"<td>{html.escape(r['component'].replace('_', ' '))}</td>"
        f"<td class='num'>{P._fmt(r['value'])}</td>"
        f"<td class='sd num'>{('± ' + P._fmt(r['sd'])) if r['sd'] else ''}</td>"
        f"<td class='sd num'>{('± ' + P._fmt(r['sem'])) if r['sem'] else ''}</td>"
        f"</tr>" for r in table["rows"])
    caveats = ""
    lead, more = table.get("lead_caveat"), table.get("caveats") or []
    if lead:
        caveats += f"<p class='hint warn'>{html.escape(lead)}</p>"
    if more:
        caveats += ("<details><summary class='hint'>"
                    f"{len(more)} more caveat(s) about these numbers</summary>"
                    + "".join(f"<p class='hint warn'>{html.escape(c)}</p>"
                              for c in more) + "</details>")
    return (f'<div class="card"><h3>Binding free energy '
            f'<span class="tag">{html.escape(table["method"])}</span></h3>'
            f'<table><tr><th>component</th><th class="num">kcal/mol</th>'
            f'<th class="num">SD</th><th class="num">SEM</th></tr>{body}</table>'
            f'</div>'), caveats


# --------------------------------------------------------------------------
# THE RECORD
# --------------------------------------------------------------------------
PENDING = ("title", "related_publication", "doi", "authors", "license_holder")


def _record_section(doc, frames_doc, rows):
    """The metadata, with what is still missing shown rather than hidden.

    The owner tunes DOI and authors by hand at the end, so the page states
    which fields are waiting instead of quietly omitting them -- a reader can
    see the record is a draft, and whoever fills it in has the list.
    """
    if not doc:
        return ""
    samp = doc.get("sampling") or {}
    setup = doc.get("setup") or {}
    md = (frames_doc or {}).get("mds") or []
    analysed = ""
    if md:
        first = md[0]
        files = first.get("analysed_files") or []
        n_raw = first.get("n_analysed") or first.get("n_raw")
        analysed = " \u00b7 ".join(filter(None, [
            files[0] if files else "the production segments",
            f"{n_raw:,} frames" if isinstance(n_raw, int) else None,
            f"{first['framestep_ns']} ns per frame"
            if first.get("framestep_ns") else None]))

    def row(k, v, pending=False):
        if v in (None, "", "None"):
            v, pending = "to be added before registration", True
        cls = " class='pending'" if pending else ""
        return (f"<tr><td class='k'>{html.escape(k)}</td>"
                f"<td{cls}>{html.escape(str(v))}</td></tr>")

    body = "".join([
        row("identifier", doc.get("identifier")),
        row("title", doc.get("title")),
        row("topic", doc.get("topic")),
        row("lab", doc.get("lab")),
        row("contact", doc.get("contact")),
        row("licence", doc.get("license")),
        row("publication / DOI", doc.get("related_publication")),
        row("authors", doc.get("authors")),
        row("force field", setup.get("force_field")),
        row("software", setup.get("software")),
        row("ensemble", f"{setup.get('ensemble')}, {setup.get('temperature_K')} K"
            if setup.get("ensemble") else None),
        row("replicas", f"{samp.get('replicas')} × "
            f"{samp.get('length_per_replica_ns')} ns"
            if samp.get("replicas") else None),
        row("frame interval", f"{samp.get('frame_interval_ps')} ps"
            if samp.get("frame_interval_ps") else None),
        row("analysed", analysed or None),
        row("data", doc.get("data_availability")),
    ])
    n_rows = len(rows)
    return (f'<div class="card"><h3>Record '
            f'<span class="tag">{html.escape(str(doc.get("record") or ""))}'
            f'</span></h3><table>{body}</table>'
            f'<p class="hint">Every number on this page comes from '
            f'{n_rows:,} rows in <code>data/</code> (TSV and Parquet, one row '
            f'per measured value with the frame, replica and stride it came '
            f'from). Fields marked in amber are filled in by hand before the '
            f'record is registered.</p></div>')


# --------------------------------------------------------------------------
# the page
# --------------------------------------------------------------------------
def _blurb(doc):
    """The owner's description, under the title, when they wrote one.

    ONLY when they wrote one. Every record carries a `description`, and where
    nobody supplied it the metadata pipeline composes a serviceable sentence
    out of the atom count and the setup. That belongs in a catalogue listing.
    At the top of a published page it would restate the record table three
    screens further down and read as though a person had wanted it said.
    `description_basis` is what tells the two apart; the text alone cannot.

    The text takes the same `**bold**` markup the figure captions take, and is
    escaped before the markup is applied, so a residue name or a file path in
    it can never open a tag.
    """
    doc = doc or {}
    if doc.get("description_basis") != "declared":
        return ""
    text = str(doc.get("description") or "").strip()
    return f'<p class="blurb">{_rich(text)}</p>' if text else ""


def build(out_path, title, subtitle, rows, frames_doc, topology_name="merged.psf",
          doc=None, bundle=None, components=(), paper=None):
    """The story page. With `paper`, the article's figures are set among ours.

    ONE BUILDER, TWO PAGES. `page_paper_story.html` is this page with the
    paper's material placed into it, not a second page that happens to show
    the same curves. Anything else and the two would drift: a caption fixed
    here would go on saying the old number there, which is the failure this
    whole kit is built to make impossible.
    """
    from plotly.offline import get_plotlyjs
    paper = paper or {}

    two = L.has_ligand(rows)
    protein = S.panels_for(rows, "protein")
    ligand = S.panels_for(rows, "ligand") if two else []
    cmap = S.contact_map(rows) if two else None
    hist_r = S.histogram(rows, "receptor") if two else None
    hist_l = S.histogram(rows, "ligand") if two else None
    binding = S.binding_table(rows) if two else None

    panels = []
    for i, p in enumerate(protein):
        p["div"] = f"p{i}"
        panels.append(p)
    for i, p in enumerate(ligand):
        p["div"] = f"l{i}"
        panels.append(p)

    # --- THE SYSTEM
    figs, legend, figwarn = _pictures(os.path.dirname(out_path))
    system_body = (figs + legend
                   + _cap(C.system(rows, (doc or {}).get("setup"),
                                   bundle, components))
                   + _composition(rows) + figwarn)
    # Which sections this page actually draws. A figure routed at one that is
    # not here falls back to the paper's own section rather than vanishing.
    live = {"paper"}
    if system_body.strip():
        live.add("system")
    _MEASURED = '<span class="srctag measured">measured here</span>'

    def _with_figs(body, section):
        """(body + the article's figures for this section, the lede tag)."""
        figs = PA.section_html(paper, section, live) if paper else ""
        if not figs:
            return body, ""
        return body + figs, (_MEASURED if body.strip() else "")

    paper_html = _section(
        "The paper", paper.get("title") or "The published article",
        "This dataset is the one behind a published paper. Its figures appear "
        "on this page beside the measurements, each marked with its number.",
        PA.citation(paper, (doc or {}).get("authors"))) if paper else ""

    system_body, system_tag = _with_figs(system_body, "system")
    doc_html = _section(
        "The system", "What was simulated",
        "The assembly, as it was built and solvated.", system_body,
        tag=system_tag)

    # --- HOW IT MOVES. Protein and its ligand counterpart in sequence, so a
    # reader meets RMSD once and then RMSD of the other thing, rather than
    # meeting every protein analysis and then starting again.
    moves = ""
    by_obs = {}
    for p in ligand:
        by_obs.setdefault(p["observable"], []).append(p)
    used, ordered = set(), []
    for p in protein:
        ordered.append(p)
        for q in by_obs.get(p["observable"], []):
            if id(q) not in used:
                used.add(id(q))
                ordered.append(q)
                break
    ordered += [q for q in ligand if id(q) not in used]
    for p in ordered:
        moves += (f'<div class="card plot" id="{p["div"]}"></div>'
                  + _cap(C.for_panel(rows, p)))
    if moves:
        moves += ('<p class="click">Click any point on a curve to see which '
                  '<b>frame</b> it is, and a VMD command that opens exactly '
                  'that frame.</p>')
    if moves.strip():
        live.add("moves")
    moves, moves_tag = _with_figs(moves, "moves")
    moves_html = _section(
        "How it moves", "Stability and flexibility",
        "How far the structure travels from where it started, and which parts "
        "of it move most.", moves, tag=moves_tag)

    # --- WHAT IT TOUCHES
    touch = ""
    if cmap:
        touch += ('<div class="card plot" id="contactmap"></div>'
                  + _cap(C.contact_map(rows, cmap)))
    bars = ""
    if hist_r:
        bars += ('<div><div class="card plot" id="histrec"></div>'
                 + _cap(C.histogram(rows, "receptor")) + "</div>")
    if hist_l:
        bars += ('<div><div class="card plot" id="histlig"></div>'
                 + _cap(C.histogram(rows, "ligand")) + "</div>")
    if bars:
        touch += f'<div class="pair">{bars}</div>'
    if touch.strip():
        live.add("touch")
    touch, touch_tag = _with_figs(touch, "touch")
    touch_html = _section(
        "What it touches", "The interface",
        "Which residues of the two sides are in contact, and for what "
        "fraction of the simulation.", touch, family="ligand", tag=touch_tag)

    # --- HOW STRONGLY
    strong = ""
    if binding:
        card, caveats = _binding_table(binding)
        strong = card + _cap(C.binding(rows, binding)) + caveats
    strong_html = _section(
        "How strongly", "Binding free energy",
        "An end-state estimate of how favourable the bound state is.",
        strong, family="ligand")

    # --- THE SOLVENT
    curves = S.rdf_curves(rows)
    stable = S.solvent_table(rows)
    solvent = ""
    if curves:
        solvent += '<div class="card plot" id="rdfplot"></div>'
    if stable:
        trs = "".join(
            f"<tr><td class='k'>{html.escape(t['label'])}</td>"
            f"<td class='num'>{P._fmt(t['value'])}</td>"
            f"<td class='sd num'>{_sd_like(t['value'], t['sd'])}</td>"
            f"<td class='sd'>{html.escape(t['unit'])}</td></tr>" for t in stable)
        solvent += (f'<div class="card"><h3>Water and ions</h3><table>'
                    f'<tr><th>quantity</th><th class="num">value</th>'
                    f'<th class="num">SD</th><th></th></tr>{trs}</table></div>')
    if solvent:
        solvent += _cap(C.solvent(rows))
        grids = S.grid_files(rows)
        if grids:
            peak = max((g["peak"] or 0) for g in grids
                       if g["species"].startswith("water"))
            solvent += (f'<p class="note">{len(grids)} density volumes are in '
                        f'<code>grids/</code> (.dx, openable in VMD; an '
                        f'isovalue of 5 means five times bulk). Water peaks at '
                        f'{peak:.1f}× bulk.</p>')
    if solvent.strip():
        live.add("solvent")
    solvent, solvent_tag = _with_figs(solvent, "solvent")
    solvent_html = _section(
        "The solvent", "Water and ions",
        "The solvent was kept, not stripped, so it can be measured: where it "
        "sits, how it is structured, and what it bridges.", solvent,
        tag=solvent_tag)

    # Everything the paper carries that no measured section claimed, in the
    # article's own order.
    rest_html = _section(
        "From the paper", "The published figures",
        "The rest of the article's figures, with their captions as published.",
        PA.section_html(paper, "paper", live)) if paper else ""

    record_html = _section(
        "The record", "Metadata",
        "What accompanies this dataset when it is registered.",
        _record_section(doc, frames_doc, rows))

    data = {
        "panels": panels,
        "frames": frames_doc, "topology_name": topology_name,
        "contact_map": cmap, "hist_receptor": hist_r, "hist_ligand": hist_l,
        "color_protein": L.PROTEIN_GRADIENT[1], "color_ligand": L.LIGAND_GRADIENT[1],
        "grad_protein": L.PROTEIN_GRADIENT, "grad_ligand": L.LIGAND_GRADIENT,
        "rdf": curves,
    }
    blurb = _blurb(doc)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>{CSS}{PA.CSS if paper else ''}</style>
<script>{get_plotlyjs()}</script>
</head><body>
<div class="wrap">
<header>
  <h1>{html.escape(title)}</h1>
  <div class="sub">{html.escape(subtitle)}</div>
  {blurb}
</header>
{paper_html}
{doc_html}
{moves_html}
{touch_html}
{strong_html}
{solvent_html}
{rest_html}
{record_html}
</div>
<div id="detail"><button class="x" onclick="hideDetail()">×</button>
  <div id="detail-body"></div></div>
<script>window.__DATA__ = {json.dumps(data, separators=(",", ":"), default=str)};</script>
<script>{P.JS}</script>
</body></html>"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(page)
    return out_path
