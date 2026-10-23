#!/usr/bin/env python3
"""
plotly_page.py -- variant (a): one self-contained page, plotly, click-to-frame.

WHAT MAKES IT WORTH BUILDING. A published RMSD curve is a picture; this is the
same curve where clicking a point answers "which frame is that?" -- the frame
number, the time, the replica, the exact `.dcd` it lives in and the local index
inside that file, plus a VMD command that opens it. A number nobody can trace
back to a frame cannot be checked.

SELF-CONTAINED, DELIBERATELY. Everything is inlined -- the data as JSON and
plotly.js itself. The existing viewer in GLOBAL_3STATE/80_products.py records
the reason: "a page opened over file:// cannot fetch sibling files and would
otherwise render empty". The tidy `.tsv`/`.parquet` files stay beside the page
for download and for scripts, but the page never depends on them.

NO IMPORTS FROM THE MD STACK. This module may not import MDAnalysis: pages are
built from the tables, never from trajectories, so a figure can always be
rebuilt in seconds without touching 23 GB of data.
"""

import html
import json
import os

TITLE_CSS = """
:root{
  --ink:#16181d;--muted:#6a707c;--line:#e3e6ea;--bg:#fff;--card:#f7f8fa;
  --accent:#1b6ca8;--warn:#b2182b;
  --sans:"Inter","Helvetica Neue",-apple-system,BlinkMacSystemFont,"Segoe UI",
         Roboto,"Noto Sans",Arial,sans-serif;
  --mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,
         "Liberation Mono",monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ink:#e9ebef;--muted:#9aa1ad;--line:#2b2f36;--bg:#0f1115;--card:#171a20;
  --accent:#6fb8ef;--warn:#ef8a86}}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);
     font-size:15px;line-height:1.55;letter-spacing:-0.005em;
     font-feature-settings:"kern" 1,"liga" 1;-webkit-font-smoothing:antialiased}
header{padding:26px 32px 16px;border-bottom:1px solid var(--line)}
h1{margin:0 0 6px;font-size:23px;font-weight:650;letter-spacing:-0.02em}
h2{font-size:12px;font-weight:650;letter-spacing:0.08em;text-transform:uppercase;
   color:var(--muted);margin:26px 0 10px}
.sub{color:var(--muted);font-size:14px}
main{display:grid;grid-template-columns:minmax(0,1fr) 360px;gap:26px;
     padding:8px 32px 56px;max-width:1500px}
@media (max-width:980px){main{grid-template-columns:minmax(0,1fr)}}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;
      padding:14px 16px;margin-bottom:16px}
.plot{padding:6px 4px 2px}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th,td{text-align:left;padding:5px 10px 5px 0;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:600;font-size:12px;letter-spacing:0.04em;
   text-transform:uppercase}
tr:last-child td{border-bottom:none}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums lining-nums}
#detail{position:sticky;top:18px}
#detail .k{color:var(--muted);font-size:11px;letter-spacing:0.07em;
           text-transform:uppercase;margin-top:11px}
#detail .v{font-size:14px;word-break:break-word}
#detail .v b{font-weight:650;font-size:16px}
code{font-family:var(--mono);font-size:12.5px;
     background:rgba(128,128,128,.13);padding:1px 5px;border-radius:4px}
pre{background:rgba(128,128,128,.13);padding:10px 11px;border-radius:7px;
    overflow-x:auto;font-family:var(--mono);font-size:11.5px;line-height:1.55;
    margin:8px 0;white-space:pre-wrap;word-break:break-all}
button{font:inherit;font-size:12.5px;padding:5px 11px;border:1px solid var(--line);
       background:var(--bg);color:var(--ink);border-radius:6px;cursor:pointer;
       transition:border-color .12s,color .12s}
button:hover{border-color:var(--accent);color:var(--accent)}
.hint{color:var(--muted);font-size:12.5px;margin-top:8px;line-height:1.5}
.tag{display:inline-block;font-size:11px;letter-spacing:.06em;color:var(--muted);
     border:1px solid var(--line);border-radius:999px;padding:1px 9px;
     margin-right:6px}
"""

JS = """
const D = window.__DATA__;
const COLOR = {receptor:"#1b6ca8", ligand:"#d1495b", glycans:"#3f8f4a",
               solute:"#7a5195", protein:"#8a8f98", epitope:"#e07b39"};
const DASH  = ["solid","dash","dot","longdash"];

function label(s){ return D.labels[s] || (s.charAt(0).toUpperCase()+s.slice(1)); }

function provenance(replica, frame){
  const md = D.frames.mds.find(m => m.md_id === replica);
  if(!md) return null;
  const row = md.frames[frame];
  if(!row) return null;
  const [step, segi, local, raw] = row;
  // The file the analysis read -- the merged trajectory -- and the frame's
  // index INSIDE it. The pr.part* segment it originally came from is an
  // implementation detail nobody opens.
  const file = (md.analysed_files||[])[0] || md.segments[segi].file;
  const path = md.analysed_path || md.segments[segi].path;
  return {frame, raw, step, time_ns:+(step*md.step_ps/1000).toFixed(2),
          replica, file, path, index:(md.analysed_path ? raw : local),
          of:(md.analysed_path ? md.n_raw : md.segments[segi].frames),
          framestep_ns:md.framestep_ns};
}

function valuesAt(replica, frame){
  const out = [];
  for(const s of D.series){
    if(s.replica !== replica || s.kind !== "time") continue;
    const i = s.frames.indexOf(frame);
    if(i >= 0) out.push([label(s.observable)+" · "+label(s.subset),
                         s.values[i].toFixed(2)+" "+s.unit]);
  }
  return out;
}

function showDetail(replica, frame){
  const el = document.getElementById("detail-body");
  if(frame < 0){
    el.innerHTML = "<p class='hint'>That point is a residue, not a frame — " +
      "per-residue values are averages over the whole run.</p>";
    return;
  }
  const p = provenance(replica, frame);
  if(!p){ el.innerHTML = "<p class='hint'>No frame record for that point.</p>"; return; }
  const vals = valuesAt(replica, frame)
      .map(([k,v]) => `<tr><td>${k}</td><td class="num">${v}</td></tr>`).join("");
  el.innerHTML = `
    <div class="k">frame</div><div class="v"><b>${p.frame}</b>
      <span class="hint">raw ${p.raw} in the concatenation</span></div>
    <div class="k">time</div><div class="v"><b>${p.time_ns} ns</b>
      <span class="hint">at ${p.framestep_ns} ns per frame</span></div>
    <div class="k">replica</div><div class="v">${p.replica}</div>
    <div class="k">step</div><div class="v">${p.step.toLocaleString()}</div>
    <div class="k">lives in</div><div class="v"><b>${p.file}</b>
      <span class="hint">frame ${p.index} of ${p.of}</span></div>
    <table style="margin-top:12px">${vals}</table>
    <div class="k">open exactly this frame</div>
    <pre id="cmd">mol new ${D.topology_name}
mol addfile ${p.path} type dcd first ${p.index} last ${p.index} waitfor all</pre>
    <button onclick="navigator.clipboard.writeText(
        document.getElementById('cmd').textContent)">Copy VMD command</button>`;
}

function draw(){
  const dark = window.matchMedia &&
               window.matchMedia("(prefers-color-scheme: dark)").matches;
  const ink = dark ? "#e9ebef" : "#16181d";
  const grid = dark ? "rgba(255,255,255,.08)" : "rgba(0,0,0,.07)";
  const by = {};
  for(const s of D.series){ (by[s.observable] = by[s.observable] || []).push(s); }
  const order = ["rmsd","rg","rmsf"];
  const keys = Object.keys(by).sort((a,b) =>
      (order.indexOf(a)+99)%99 - (order.indexOf(b)+99)%99);
  let n = 0;
  for(const obs of keys){
    const series = by[obs];
    const wrap = document.createElement("div");
    wrap.className = "card plot";
    wrap.id = "plot" + (n++);
    document.getElementById("plots").appendChild(wrap);
    const reps = [...new Set(series.map(s => s.replica))];
    const traces = series.map(s => ({
      x:s.x, y:s.values, name:label(s.subset)+" · "+s.replica,
      legendgroup:s.subset, type:"scattergl", mode:"lines",
      line:{width:1.3, color:COLOR[s.subset]||"#8a8f98",
            dash:DASH[reps.indexOf(s.replica) % DASH.length]},
      customdata:s.frames.map((f,i) => [s.replica, f, s.labels ? s.labels[i] : ""]),
      hovertemplate:"<b>%{y:.2f} "+s.unit+"</b><br>"+
        (s.kind==="time" ? "%{x:.1f} ns · frame %{customdata[1]}"
                         : "%{customdata[2]}")+
        "<extra>%{fullData.name}</extra>"
    }));
    Plotly.newPlot(wrap.id, traces, {
      margin:{l:64,r:16,t:42,b:52}, height:340,
      title:{text:D.titles[obs] || label(obs), x:0, xanchor:"left",
             font:{size:15, color:ink, family:"Inter,Helvetica Neue,sans-serif"}},
      xaxis:{title:{text:series[0].xlabel, font:{size:13}}, gridcolor:grid,
             zeroline:false, ticks:"outside", ticklen:4, tickfont:{size:12}},
      yaxis:{title:{text:label(obs)+(series[0].unit ? " ("+series[0].unit+")" : ""),
                    font:{size:13}},
             gridcolor:grid, zeroline:false, ticks:"outside", ticklen:4,
             tickfont:{size:12}},
      legend:{orientation:"h", y:-0.24, font:{size:11.5}, itemsizing:"constant"},
      hoverlabel:{font:{size:12, family:"Inter,sans-serif"}},
      paper_bgcolor:"rgba(0,0,0,0)", plot_bgcolor:"rgba(0,0,0,0)",
      font:{color:ink, family:"Inter,Helvetica Neue,Arial,sans-serif", size:12}
    }, {displaylogo:false, responsive:true,
        modeBarButtonsToRemove:["select2d","lasso2d","autoScale2d"]});
    document.getElementById(wrap.id).on("plotly_click", ev => {
      const cd = ev.points[0].customdata;
      if(cd) showDetail(cd[0], cd[1]);
    });
  }
}
document.addEventListener("DOMContentLoaded", draw);
"""


def _plotly_js():
    """plotly.js, inlined, so the page works from a file:// URL."""
    from plotly.offline import get_plotlyjs
    return get_plotlyjs()


def series_from(rows, observables=("rmsd", "rg"), subsets=None, scope="frame"):
    """Per-frame rows -> one series per (replica, subset, observable)."""
    out = {}
    for r in rows:
        if r.get("scope") != scope or r.get("observable") not in observables:
            continue
        if subsets and r.get("subset") not in subsets:
            continue
        key = (r["replica"], r["subset"], r["observable"])
        s = out.setdefault(key, {"replica": r["replica"], "subset": r["subset"],
                                 "observable": r["observable"], "kind": "time",
                                 "unit": r.get("unit") or "",
                                 "xlabel": "time (ns)", "xunit": "ns",
                                 "x": [], "values": [], "frames": []})
        s["x"].append(r.get("time_ns"))
        s["values"].append(float(r["value"]))
        s["frames"].append(int(r.get("frame") if r.get("frame") is not None
                               else r.get("scope_id")))
    return list(out.values())


def residue_series(rows, observable="rmsf", subsets=None):
    """Per-residue rows -> one series per (replica, subset), x = residue."""
    out = {}
    for r in rows:
        if r.get("scope") != "residue" or r.get("observable") != observable:
            continue
        if subsets and r.get("subset") not in subsets:
            continue
        key = (r["replica"], r["subset"])
        s = out.setdefault(key, {"replica": r["replica"], "subset": r["subset"],
                                 "observable": observable, "kind": "residue",
                                 "unit": r.get("unit") or "",
                                 "xlabel": "residue (in selection order)",
                                 "xunit": "", "x": [], "values": [],
                                 "frames": [], "labels": []})
        s["x"].append(len(s["x"]))
        s["values"].append(float(r["value"]))
        # A residue is not a frame: clicking one cannot open a frame, so the
        # provenance panel is told there is none rather than shown a wrong one.
        s["frames"].append(-1)
        s["labels"].append(f"{r.get('segid') or ''} {r.get('resname') or ''}"
                           f"{r.get('resid') or ''}".strip())
    return list(out.values())


def composition_table(rows):
    """The composition numbers, per replica, as rows for the table card."""
    order = ["water_molecules", "ions_SOD", "ions_CLA", "ions_total",
             "salt_concentration", "box_volume", "box_edge_a", "box_edge_b",
             "box_edge_c", "solute_atoms", "glycan_residues", "atoms_total"]
    pretty = {"water_molecules": "water molecules (TIP3P)",
              "ions_SOD": "Na+ ions", "ions_CLA": "Cl- ions",
              "ions_total": "ions, total",
              "salt_concentration": "salt concentration (mol/L)",
              "box_volume": "box volume (A^3)", "box_edge_a": "box a (A)",
              "box_edge_b": "box b (A)", "box_edge_c": "box c (A)",
              "solute_atoms": "solute atoms", "glycan_residues": "glycan residues",
              "atoms_total": "atoms, total"}
    got = {}
    reps = []
    for r in rows:
        if r.get("scope") != "run" or r.get("stage") != "composition":
            continue
        key = r["observable"] + (f"_{r['subset']}" if r.get("subset") else "")
        got.setdefault(key, {})[r["replica"]] = r["value"]
        if r["replica"] not in reps:
            reps.append(r["replica"])
    table = []
    for k in order:
        if k in got:
            table.append([pretty.get(k, k),
                          [got[k].get(rep) for rep in reps]])
    return reps, table


def _fmt(v):
    if v is None:
        return "&mdash;"
    if isinstance(v, float):
        return f"{v:,.4g}" if abs(v) < 1000 else f"{v:,.0f}"
    return f"{v:,}" if isinstance(v, int) else html.escape(str(v))


def build(out_path, title, subtitle, series, frames_doc, comp_rows,
          topology_name="merged.psf", notes=()):
    reps, table = composition_table(comp_rows)
    head = "".join(f"<th class='num'>{html.escape(str(r))}</th>" for r in reps)
    body = "".join(
        "<tr><td>" + html.escape(name) + "</td>" +
        "".join(f"<td class='num'>{_fmt(v)}</td>" for v in vals) + "</tr>"
        for name, vals in table)

    data = {"series": series, "frames": frames_doc,
            "topology_name": topology_name,
            # Observables keep their capitals -- RMSD and RMSF are acronyms,
            # not words -- and each plot gets a title that says what was
            # measured against what.
            "labels": {"rmsd": "RMSD", "rmsf": "RMSF", "rg": "Rg",
                       "sasa": "SASA", "receptor": "Receptor",
                       "ligand": "Ligand", "glycans": "Glycans",
                       "epitope": "Epitope", "solute": "Solute",
                       "protein": "Protein"},
            "titles": {
                "rmsd": "RMSD \u2014 heavy atoms, fitted on the receptor C\u03b1",
                "rg": "Rg \u2014 radius of gyration",
                "rmsf": "RMSF \u2014 per residue, about the mean structure"}}
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>{TITLE_CSS}</style>
<script>{_plotly_js()}</script>
</head><body>
<header>
  <h1>{html.escape(title)}</h1>
  <div class="sub">{html.escape(subtitle)}</div>
</header>
<main>
  <div>
    <h2>Analyses</h2>
    <div id="plots"></div>
    <p class="hint">Click any point on a time series to see which frame it is.</p>
  </div>
  <div>
    <div class="card" id="detail">
      <h2 style="margin-top:0">The frame you clicked</h2>
      <div id="detail-body"><p class="hint">Nothing selected yet. Click a
      point on a curve and this panel will name the frame, its time, the file
      it lives in and the command that opens it.</p></div>
    </div>
    <div class="card">
      <h2 style="margin-top:0">System composition</h2>
      <table><tr><th>quantity</th>{head}</tr>{body}</table>
      <p class="hint">Counted from the topology; box from the NAMD
      <code>.xst</code> records, salt computed from the ions and the mean box.</p>
    </div>
    {"".join(f'<div class="card"><p class="hint">{html.escape(n)}</p></div>' for n in notes)}
  </div>
</main>
<script>window.__DATA__ = {json.dumps(data, separators=(",", ":"))};</script>
<script>{JS}</script>
</body></html>"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(doc)
    return out_path
