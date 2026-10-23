#!/usr/bin/env python3
"""
report.py -- the page: protein on the left, ligand on the right, apo landscape.

Built entirely from the tidy rows in `data/`, so it can be rebuilt in seconds
(`analyse.py --pages-only`) and a change to how something is drawn can never
change what was measured. Everything is inlined -- plotly.js and the data --
because a page opened over `file://` cannot fetch its siblings and would
otherwise render empty.

WHAT DECIDES THE LAYOUT. Whether anything was measured on a ligand. An apo
trimer gets one landscape column; a complex gets two, protein left and ligand
right, with the contact map, the per-residue histograms and the binding-energy
table under the ligand. Nothing is configured: a system that has not had its
ligand stages run yet shows the landscape page and grows the second column when
the numbers arrive.
"""

import html
import json
import os

from . import captions as C
from . import layout as L
from . import series as S

CSS = """
:root{
  --ink:#16181d;--muted:#6a707c;--line:#e3e6ea;--bg:#fff;--card:#f7f8fa;
  --protein:#1596b8;--ligand:#d1495b;
  --sans:"Inter","Helvetica Neue",-apple-system,BlinkMacSystemFont,"Segoe UI",
         Roboto,"Noto Sans",Arial,sans-serif;
  --mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ink:#e9ebef;--muted:#9aa1ad;--line:#2b2f36;--bg:#0f1115;--card:#171a20;
  --protein:#35b8d6;--ligand:#e8743c}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);
     font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased}
header{padding:26px 34px 14px;border-bottom:1px solid var(--line)}
h1{margin:0 0 6px;font-size:23px;font-weight:650;letter-spacing:-.02em}
.sub{color:var(--muted);font-size:14px}
main{padding:10px 34px 60px;max-width:1800px;margin:0 auto}
.cols{display:grid;gap:20px 22px;align-items:stretch}
.cols.two{grid-template-columns:repeat(2,minmax(0,1fr))}
.cols.one{grid-template-columns:minmax(0,1fr)}
@media (max-width:1100px){.cols.two{grid-template-columns:minmax(0,1fr)}}
/* A panel with no counterpart takes the whole row rather than leaving a hole
   beside it: an empty half reads as something that failed to load. */
.span2{grid-column:1 / -1}
/* Card bands size themselves to their content and wrap, so a short table and
   a tall one never leave a column of white space under the shorter. */
.twocol{display:grid;gap:0 22px;
        grid-template-columns:repeat(auto-fit,minmax(min(100%,215px),1fr))}
.twocol table tr:last-child td{border-bottom:none}
.band{display:grid;gap:20px 22px;align-items:stretch;margin-top:4px;
      grid-template-columns:repeat(auto-fit,minmax(min(100%,430px),1fr))}
.heads{display:grid;gap:22px;grid-template-columns:repeat(2,minmax(0,1fr));
       margin:20px 0 2px}
@media (max-width:1100px){.heads{grid-template-columns:minmax(0,1fr)}}
.colhead{display:flex;align-items:baseline;gap:10px;margin:0 0 4px}
.colhead h2{margin:0;font-size:13px;font-weight:700;letter-spacing:.1em;
            text-transform:uppercase}
.colhead .what{color:var(--muted);font-size:12.5px}
.protein .colhead h2{color:var(--protein)}
.ligand  .colhead h2{color:var(--ligand)}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;
      padding:12px 14px;margin-bottom:16px}
.card.plot{padding:6px 4px 2px}
.tag{display:inline-block;font-size:11px;letter-spacing:.05em;color:var(--muted);
     border:1px solid var(--line);border-radius:999px;padding:1px 9px}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th,td{text-align:left;padding:5px 10px 5px 0;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:600;font-size:11.5px;letter-spacing:.05em;
   text-transform:uppercase}
tr:last-child td{border-bottom:none}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums lining-nums}
td.sd{color:var(--muted);font-size:12.5px}
h3{font-size:14px;margin:0 0 8px;font-weight:650}
.hint{color:var(--muted);font-size:12.5px;margin:8px 0 0;line-height:1.5}
.warn{color:#b2182b}
@media (prefers-color-scheme:dark){.warn{color:#ef8a86}}
#detail{position:fixed;right:18px;bottom:18px;width:370px;max-width:92vw;
        background:var(--card);border:1px solid var(--line);border-radius:11px;
        padding:14px 16px;box-shadow:0 10px 34px rgba(0,0,0,.18);display:none;
        z-index:20;max-height:72vh;overflow:auto}
#detail.on{display:block}
#detail .k{color:var(--muted);font-size:10.5px;letter-spacing:.07em;
           text-transform:uppercase;margin-top:10px}
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
"""

JS = r"""
const D = window.__DATA__;
const dark = window.matchMedia &&
             window.matchMedia("(prefers-color-scheme: dark)").matches;
const INK = dark ? "#e9ebef" : "#16181d";
const GRID = dark ? "rgba(255,255,255,.08)" : "rgba(0,0,0,.07)";

function rgba(hex, a){
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${(n>>16)&255},${(n>>8)&255},${n&255},${a})`;
}

function base(title, xlabel, ylabel, h){
  return {margin:{l:62,r:14,t:34,b:46}, height:h||320,
    title:{text:title, x:0, xanchor:"left",
           font:{size:13.5, color:INK, family:"Inter,sans-serif"}},
    xaxis:{title:{text:xlabel, font:{size:12}}, gridcolor:GRID, zeroline:false,
           ticks:"outside", ticklen:4, tickfont:{size:11.5}},
    yaxis:{title:{text:ylabel, font:{size:12}}, gridcolor:GRID, zeroline:false,
           ticks:"outside", ticklen:4, tickfont:{size:11.5}},
    legend:{orientation:"h", y:-0.22, font:{size:11}},
    paper_bgcolor:"rgba(0,0,0,0)", plot_bgcolor:"rgba(0,0,0,0)",
    font:{color:INK, family:"Inter,Helvetica Neue,Arial,sans-serif", size:11.5},
    hoverlabel:{font:{size:11.5}}};
}
const CFG = {displaylogo:false, responsive:true,
             modeBarButtonsToRemove:["select2d","lasso2d","autoScale2d"]};

function drawPanel(div, p){
  const traces = [];
  if(p.band){
    // The shade first, as two filled traces, so the mean sits on top of it.
    traces.push({x:p.band.x, y:p.band.hi, mode:"lines", type:"scatter",
                 line:{width:0}, hoverinfo:"skip", showlegend:false});
    traces.push({x:p.band.x, y:p.band.lo, mode:"lines", type:"scatter",
                 line:{width:0}, fill:"tonexty", fillcolor:rgba(p.color,.20),
                 hoverinfo:"skip", name:"± SD"});
    traces.push({x:p.band.x, y:p.band.mean, mode:"lines", type:"scattergl",
                 line:{width:1.7, color:p.color}, name:"mean",
                 customdata:p.band.x.map((_,i)=>[p.band.sd[i], p.band.labels[i]]),
                 hovertemplate:"<b>%{y:.2f} "+p.unit+"</b> ± %{customdata[0]:.2f}"+
                   "<br>%{x:.1f}<extra>"+p.short+"</extra>"});
  } else {
    for(const s of p.lines){
      traces.push({x:s.x, y:s.y, mode:"lines", type:"scattergl", name:s.name,
        line:{width:1.2, color:s.color},
        customdata:s.frames.map((f,i)=>[s.replica, f, s.labels[i]]),
        hovertemplate:"<b>%{y:.2f} "+p.unit+"</b><br>"+
          (p.observable==="rmsf" ? "%{customdata[2]}" : "%{x:.1f} ns · frame %{customdata[1]}")+
          "<extra>%{fullData.name}</extra>"});
    }
  }
  Plotly.newPlot(div, traces, base(p.title, p.xlabel, p.short+" ("+p.unit+")"), CFG);
  document.getElementById(div).on("plotly_click", ev => {
    const cd = ev.points[0].customdata;
    if(cd && cd.length > 1 && typeof cd[1] === "number") showDetail(cd[0], cd[1]);
  });
}

function drawMap(div, m){
  Plotly.newPlot(div, [{z:m.z, x:m.x, y:m.y, type:"heatmap",
    colorscale:[[0,"rgba(0,0,0,0)"],[0.15,"#f2d0a4"],[0.5,"#e8743c"],[1,"#a01526"]],
    zmin:0, zmax:1, colorbar:{title:{text:"occupancy",side:"right",font:{size:11}},
    thickness:10, len:0.8, tickfont:{size:10}},
    hovertemplate:"receptor %{x}<br>ligand %{y}<br><b>%{z:.2f}</b> of frames<extra></extra>"}],
    Object.assign(base("Interface contact map", "receptor residue", "ligand residue", 420),
                  {margin:{l:90,r:14,t:34,b:70}}), CFG);
}

function drawHist(div, h, palette, title){
  // One trace per chain so the legend says which chain each residue is on,
  // and horizontal so a residue name reads straight.
  const traces = h.chain_list.map((ch, i) => ({
    type:"bar", orientation:"h", name:ch,
    y:h.labels.filter((_,k)=>h.chains[k]===ch),
    x:h.values.filter((_,k)=>h.chains[k]===ch),
    error_x:{type:"data", visible:true, thickness:1, width:2,
             color:"rgba(128,128,128,.55)",
             array:h.sd.filter((_,k)=>h.chains[k]===ch)},
    marker:{color:palette[i % palette.length]},
    hovertemplate:"<b>%{y}</b><br>%{x:.2f} of frames<extra>"+ch+"</extra>"}));
  const lay = base(title, "fraction of frames in contact", "", null);
  lay.height = Math.max(240, 26 * h.labels.length + 96);
  lay.barmode = "stack";
  lay.margin = {l:86, r:16, t:34, b:44};
  lay.xaxis.range = [0, Math.min(1, Math.max(...h.values.map((v,k)=>v+h.sd[k])) * 1.12)];
  lay.yaxis.automargin = true;
  lay.legend = {orientation:"h", y:-0.14, font:{size:10.5}};
  Plotly.newPlot(div, traces, lay, CFG);
}

function provenance(replica, frame){
  const md = D.frames.mds.find(m => m.md_id === replica) || D.frames.mds[0];
  if(!md) return null;
  const row = md.frames[frame];
  if(!row) return null;
  const [step, segi, local, raw] = row;
  const file = (md.analysed_files||[])[0] || md.segments[segi].file;
  const path = md.analysed_path || md.segments[segi].path;
  // Which number names this frame INSIDE the file you are about to open.
  // A strip that kept the frames a restart re-simulated is the raw sequence;
  // one that dropped them is the de-duplicated sequence, and `frame` is
  // already that index. Guessing between them put every frame after the
  // first seam one off.
  let index = local, of = md.segments[segi].frames;
  if(md.analysed_path){
    if(md.analysed_indexing === "raw"){ index = raw; of = md.n_raw; }
    else if(md.analysed_indexing === "unique"){ index = frame; of = md.n_analysed; }
    else { index = null; of = md.n_analysed || null; }
  }
  // The step counter times a frame, except where there is no real one: a
  // post-processed copy carries a placeholder, so a declared interval times
  // the frame by its position instead.
  const time_ns = (md.framestep_basis === "declared" && md.framestep_ns)
      ? +(frame * md.framestep_ns).toFixed(2)
      : +(step * md.step_ps / 1000).toFixed(2);
  return {frame, raw, step, time_ns, replica:md.md_id, file, path, index, of,
          framestep_ns:md.framestep_ns, basis:md.framestep_basis};
}

function showDetail(replica, frame){
  const el = document.getElementById("detail");
  const body = document.getElementById("detail-body");
  const p = provenance(replica, frame);
  if(!p){ return; }
  body.innerHTML = `
    <div class="k">frame</div><div class="v"><b>${p.frame}</b></div>
    <div class="k">time</div><div class="v"><b>${p.time_ns} ns</b>
      <span class="hint">at ${p.framestep_ns} ns per frame${
        p.basis === "declared" ? ", declared in settings.yaml (the file's own"
          + " header does not state one)" : ""}</span></div>
    <div class="k">replica</div><div class="v">${p.replica}</div>
    <div class="k">lives in</div><div class="v"><b>${p.file}</b>
      <span class="hint">${p.index === null
        ? "this file's frame count matches neither the raw nor the de-duplicated"
          + " sequence, so its frame numbering cannot be stated"
        : `frame ${p.index} of ${p.of}`}</span></div>` +
    (p.index === null ? "" : `
    <div class="k">open exactly this frame</div>
    <pre id="cmd">mol new ${D.topology_name}
mol addfile ${p.path} type dcd first ${p.index} last ${p.index} waitfor all</pre>
    <button class="copy" onclick="navigator.clipboard.writeText(
        document.getElementById('cmd').textContent)">Copy VMD command</button>`);
  el.classList.add("on");
}
function hideDetail(){ document.getElementById("detail").classList.remove("on"); }

function drawRDF(div, curves){
  const traces = curves.map(c => ({x:c.x, y:c.y, mode:"lines", type:"scattergl",
    name:c.name, line:{width:1.4, color:c.color},
    hovertemplate:"r = %{x:.2f} \u00c5<br><b>g = %{y:.2f}</b><extra>%{fullData.name}</extra>"}));
  // Height chosen to finish level with the solvent table beside it: a plot
  // taller than its neighbour leaves a strip of empty card under the table.
  const lay = base("Radial distribution functions", "r (\u00c5)", "g(r)", 286);
  // g(r) -> 1 is what "no structure" looks like; drawing that line makes every
  // peak readable as a ratio against it rather than against nothing.
  lay.shapes = [{type:"line", xref:"paper", x0:0, x1:1, y0:1, y1:1,
                 line:{width:1, dash:"dot", color:"rgba(128,128,128,.6)"}}];
  Plotly.newPlot(div, traces, lay, CFG);
}

document.addEventListener("DOMContentLoaded", () => {
  for(const p of D.panels){ drawPanel(p.div, p); }
  if(D.contact_map) drawMap("contactmap", D.contact_map);
  if(D.hist_receptor) drawHist("histrec", D.hist_receptor, D.grad_protein,
                               "Receptor residues at the interface");
  if(D.hist_ligand) drawHist("histlig", D.hist_ligand, D.grad_ligand,
                             "Ligand residues at the interface");
  if(D.rdf) drawRDF("rdfplot", D.rdf);
});
"""


def _fmt(v, sig=4):
    if v is None:
        return "&mdash;"
    if isinstance(v, float):
        if abs(v) >= 10000:
            return f"{v:,.0f}"
        if abs(v) >= 1:
            return f"{v:,.{max(0, sig - len(str(int(abs(v)))))}f}"
        return f"{v:.4g}"
    return f"{v:,}" if isinstance(v, int) else html.escape(str(v))


def _composition_card(rows):
    """The composition, one value column, split into two side-by-side halves.

    Eleven rows in a single column made a tall thin card that set the height of
    the whole band; two halves of five or six read in one glance and sit level
    with whatever is beside them.
    """
    table = L.composition_single_column(rows)
    if not table:
        return ""

    def half(items):
        body = "".join(
            f"<tr><td>{html.escape(lbl)}</td><td class='num'>{_fmt(v)}</td></tr>"
            for lbl, v, sd, n in items)
        return f"<table>{body}</table>"

    mid = (len(table) + 1) // 2
    return f"""<div class="card">
  <h3>System composition</h3>
  <div class="twocol">{half(table[:mid])}{half(table[mid:])}</div>
</div>"""


def _binding_card(table):
    if not table:
        return ""
    # The total is the row a reader came for, so it is set apart. Built
    # without a backslash inside the f-string expression, which Python
    # refuses -- the reason this broke once already.
    BOLD = " style='font-weight:650'"
    body = "".join(
        f"<tr{BOLD if r.get('total') else ''}>"
        f"<td>{html.escape(r['component'].replace('_', ' '))}</td>"
        f"<td class='num'>{_fmt(r['value'])}</td>"
        f"<td class='sd num'>{('± ' + _fmt(r['sd'])) if r['sd'] else ''}</td>"
        f"<td class='sd num'>{('± ' + _fmt(r['sem'])) if r['sem'] else ''}</td>"
        f"</tr>"
        for r in table["rows"])
    zeros = table.get("dropped_zero") or []
    lead = table.get("lead_caveat")
    more = table.get("caveats") or []
    caveats = (f"<p class='hint warn'>{html.escape(lead)}</p>" if lead else "")
    if more:
        caveats += ("<details><summary class='hint'>"
                    f"{len(more)} more caveat(s) about these numbers</summary>"
                    + "".join(f"<p class='hint warn'>{html.escape(c)}</p>"
                              for c in more) + "</details>")
    return f"""<div class="card">
  <h3>Binding free energy <span class="tag">{html.escape(table['method'])}</span></h3>
  <table><tr><th>component</th><th class="num">kcal/mol</th>
  <th class="num">SD</th><th class="num">SEM</th></tr>{body}</table>
  <p class="hint">Pooled over {table['per_run']} protomer × replica runs. SD is the
  spread across them; SEM is over replicas only, because protomers of one trimer
  share a box and a thermostat.{(" The " + str(len(zeros)) + " bonded terms are identically zero and are not"
    " listed: complex, receptor and ligand come from the same frames, so they"
    " cancel exactly.") if zeros else ""}</p>{caveats}
</div>"""


def _user_figures(bundle_dir):
    """Pictures dropped into `figures_user/`, with what they were made for.

    Referenced by relative path rather than inlined: an <img src> resolves
    over file:// (a fetch does not), and four 1 MB renders inlined as base64
    would add 5 MB to a page that already carries plotly.
    """
    d = os.path.join(bundle_dir, "figures_user")
    if not os.path.isdir(d):
        return ""
    names = sorted(n for n in os.listdir(d)
                   if n.lower().endswith((".png", ".jpg", ".jpeg", ".gif",
                                          ".svg", ".webp")))
    if not names:
        return ""
    meta = L.read_credits(d)
    cards = ""
    for n in names:
        m = meta.get(n, {})
        cap = m.get("caption") or os.path.splitext(n)[0].replace("-", " ")
        src = m.get("title")
        state = m.get("state_in_title")
        rec = m.get("records_state")
        warn = ""
        if state and rec and state != rec:
            warn = (f"<p class='hint warn'>The slide it came from is titled "
                    f"'{html.escape(state)}'; these records are "
                    f"'{html.escape(rec)}'.</p>")
        cards += (f"<div class='card'>"
                  f"<img src='figures_user/{html.escape(n)}' alt='{html.escape(cap)}' "
                  f"style='width:100%;height:auto;border-radius:7px'>"
                  f"<p class='hint' style='margin-top:9px'><b>{html.escape(cap)}</b>"
                  + (f" &middot; {html.escape(src)}" if src else "") + "</p>"
                  + warn + "</div>")
    return f"""<section class="figures">
  <div class="heads"><div class="protein"><div class="colhead">
    <h2>Figures</h2><span class="what">supplied with the system</span>
  </div></div></div>
  <div class="band">{cards}</div>
</section>"""


def _solvent_section(rows):
    """Water and ions, full width under the two columns.

    The solvent belongs to neither half -- it surrounds both -- so it gets the
    landscape strip beneath them rather than being forced into a column where
    it would look like a property of the protein or of the ligand.
    """
    table = S.solvent_table(rows)
    grids = S.grid_files(rows)
    curves = S.rdf_curves(rows)
    if not (table or grids or curves):
        return ""
    body = ""
    if table:
        trs = "".join(
            f"<tr><td>{html.escape(t['label'])}</td>"
            f"<td class='num'>{_fmt(t['value'])}</td>"
            f"<td class='sd num'>{('± ' + _fmt(t['sd'])) if t['sd'] else ''}</td>"
            f"<td class='sd'>{html.escape(t['unit'])}</td></tr>" for t in table)
        body += (f"<div class='card'><h3>Solvent</h3><table>"
                 f"<tr><th>quantity</th><th class='num'>value</th>"
                 f"<th class='num'>SD</th><th></th></tr>{trs}</table>"
                 f"<p class='hint'>Averaged over the replicas; SD is their "
                 f"spread. " + html.escape(C._residence_note(rows)) + "</p>")
    if grids:
        # The per-replica grid table was a wall of filenames at the foot of the
        # page. The volumes themselves are what matter and they live in
        # `grids/`; one line says so, and every number that table held is still
        # in data/*density_grids*.tsv.
        unres = sum(1 for g in grids if not g["resolved"])
        peak = max((g["peak"] or 0) for g in grids if g["species"].startswith("water"))
        body += (f"<p class='hint'>{len(grids)} density volumes are in "
                 f"<code>grids/</code> (.dx, openable in VMD; an isovalue of 5 "
                 f"means five times bulk). Water peaks at {peak:.1f}× bulk. "
                 f"{unres} of them are marked unresolved in the tables — for the "
                 f"ions a bulk voxel expects about 0.01 counts, so those map "
                 f"where an ion was seen, not how dense it is.</p>")
    if body and not body.endswith("</div>"):
        body += "</div>"
    plot = ("<div class='card plot' id='rdfplot'></div>" if curves else "")
    # The two tables stack inside their own half so the RDF plot beside them
    # keeps a sensible aspect ratio instead of being squeezed to a strip.
    return f"""<section class="solvent">
  <div class="colhead"><h2 style="color:#3d9dc4">Solvent</h2>
    <span class="what">water and ions, measured around the whole system</span></div>
  <div class="band">{plot}<div>{body}</div></div>
</section>"""


def build(out_path, title, subtitle, rows, frames_doc, topology_name="merged.psf",
          notes=()):
    from plotly.offline import get_plotlyjs

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

    def head(family, heading, what):
        return (f'<div class="{family}"><div class="colhead">'
                f'<h2>{heading}</h2><span class="what">{html.escape(what)}</span>'
                f'</div></div>')

    def cell(p):
        return (f'<div><div class="card plot" id="{p["div"]}"></div>'
                f'<p class="hint" style="margin:-10px 0 0 4px">'
                f'{html.escape(p["note"])}</p></div>')

    # PAIR THE ROWS. The owner asked for protein left and ligand right; if the
    # two columns are simply stacked independently, the shorter one ends in a
    # column of white space. Pairing by observable puts each RMSD beside the
    # other RMSD, and a panel with no counterpart spans the full width instead
    # of leaving a hole.
    grid = ""
    if two:
        grid += '<div class="heads">' + \
                head("protein", "Protein", "receptor, glycans, the assembly") + \
                head("ligand", "Ligand", "what is bound, and the interface") + \
                "</div>"
        by_obs = {}
        for p in ligand:
            by_obs.setdefault(p["observable"], []).append(p)
        used = set()
        cells = []
        for p in protein:
            mate = None
            for q in by_obs.get(p["observable"], []):
                if id(q) not in used:
                    mate, _ = q, used.add(id(q))
                    break
            if mate:
                cells.append(cell(p) + cell(mate))
            else:
                cells.append(f'<div class="span2">{cell(p)}</div>')
        for q in ligand:
            if id(q) not in used:
                cells.append(f'<div class="span2">{cell(q)}</div>')
        grid += '<div class="cols two">' + "".join(cells) + "</div>"
    else:
        grid += head("protein", "Protein", "receptor, glycans, the assembly")
        grid += ('<div class="cols one">'
                 + "".join(cell(p) for p in protein) + "</div>")

    # The interface: the map is wide, the two histograms sit beside each other
    # beneath it, and the band wraps on a narrow screen instead of clipping.
    interface = ""
    if cmap or hist_r or hist_l:
        interface = ('<div class="heads"><div class="ligand"><div class="colhead">'
                     '<h2>Interface</h2><span class="what">'
                     'which residues touch, and how often</span></div></div></div>')
        if cmap:
            interface += (f'<div class="card plot span2" id="contactmap"></div>'
                          f'<p class="hint" style="margin:-10px 0 10px 4px">'
                          f'{html.escape(cmap["note"])}</p>')
        bars = ""
        if hist_r:
            bars += '<div class="card plot" id="histrec"></div>'
        if hist_l:
            bars += '<div class="card plot" id="histlig"></div>'
        if bars:
            interface += f'<div class="band">{bars}</div>'

    # Tables last, side by side, wrapping rather than leaving a gap.
    tables = _binding_card(binding) + _composition_card(rows)
    if tables:
        tables = f'<div class="band">{tables}</div>'

    data = {
        "panels": [{k: v for k, v in p.items() if k != "lines"} | {"lines": p["lines"]}
                   for p in panels],
        "frames": frames_doc, "topology_name": topology_name,
        "contact_map": cmap, "hist_receptor": hist_r, "hist_ligand": hist_l,
        "color_protein": L.PROTEIN_GRADIENT[1], "color_ligand": L.LIGAND_GRADIENT[1],
        "grad_protein": L.PROTEIN_GRADIENT, "grad_ligand": L.LIGAND_GRADIENT,
        "rdf": S.rdf_curves(rows),
    }
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script>
</head><body>
<header>
  <h1>{html.escape(title)}</h1>
  <div class="sub">{html.escape(subtitle)}</div>
</header>
<main>
  {grid}
  {interface}
  {tables}
  {_solvent_section(rows)}
  {_user_figures(os.path.dirname(out_path))}
  {"".join(f'<p class="hint">{html.escape(n)}</p>' for n in notes)}
</main>
<div id="detail"><button class="x" onclick="hideDetail()">×</button>
  <div id="detail-body"></div></div>
<script>window.__DATA__ = {json.dumps(data, separators=(",", ":"), default=str)};</script>
<script>{JS}</script>
</body></html>"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(doc)
    return out_path
