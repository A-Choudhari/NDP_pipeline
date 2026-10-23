#!/usr/bin/env python3
"""
series.py -- turn tidy rows into the traces a page draws.

ONE RULE ABOVE ALL: prefer the pooled mean, fall back to the replicas, and say
which you drew. `scope="aggregate"` rows carry a mean over the three protomers
and three replicas with the spread beside it (`sd`) -- that is what the owner
asked to see, with the spread as a shaded band. Where a stage has not written
those rows yet, the page still draws the per-replica curves rather than showing
an empty panel, and labels them so nobody mistakes three replicas for a mean.
"""

from collections import defaultdict

from core import io as _io

from . import layout as L


def _num(v):
    try:
        return None if v is None else float(v)
    except (TypeError, ValueError):
        return None


def band(rows, observable, subset, x="time_ns"):
    """The pooled mean and its +/- SD band, or None when not pooled yet."""
    pts = []
    for r in rows:
        if r.get("scope") != "aggregate" or r.get("observable") != observable:
            continue
        if r.get("subset") != subset:
            continue
        xv, v, sd = _num(r.get(x)), _num(r.get("value")), _num(r.get("sd"))
        if xv is None or v is None:
            continue
        pts.append((xv, v, sd or 0.0, r.get("n"), r.get("resid"),
                    r.get("segid"), r.get("resname")))
    if not pts:
        return None
    pts.sort(key=lambda p: p[0])
    return {
        "kind": "band",
        "x": [p[0] for p in pts],
        "mean": [p[1] for p in pts],
        "lo": [p[1] - p[2] for p in pts],
        "hi": [p[1] + p[2] for p in pts],
        "sd": [p[2] for p in pts],
        "n": pts[0][3],
        "labels": [f"{p[5] or ''} {p[6] or ''}{p[4] or ''}".strip() for p in pts],
    }


def per_replica(rows, observable, subset, scope="frame", x="time_ns"):
    """One trace per replica (and per protomer where the rows carry one)."""
    out = defaultdict(lambda: {"x": [], "y": [], "frames": [], "labels": []})
    for r in rows:
        if r.get("scope") != scope or r.get("observable") != observable:
            continue
        if r.get("subset") != subset:
            continue
        v = _num(r.get("value"))
        if v is None:
            continue
        # A per-protomer row is a DIFFERENT observable from the whole-assembly
        # one, not another sample of it. Concatenating them made one replica's
        # curve out of three protomers' -- 504 points where 126 were measured.
        if not _io.absent(r.get("protomer")):
            continue
        key = (str(r.get("replica") or ""), "")
        s = out[key]
        if scope == "residue":
            s["x"].append(len(s["x"]))
            s["labels"].append(
                f"{r.get('segid') or ''} {r.get('resname') or ''}"
                f"{r.get('resid') or ''}".strip())
            s["frames"].append(-1)
        else:
            s["x"].append(_num(r.get(x)))
            s["labels"].append("")
            f = r.get("frame")
            s["frames"].append(int(f) if f is not None else -1)
        s["y"].append(v)
    series = []
    for (rep, prot), s in sorted(out.items()):
        name = rep + (f" · {prot}" if prot else "")
        series.append({"kind": "line", "name": name, "replica": rep,
                       "protomer": prot, **s})
    return series


def panel(rows, observable, subset, title, unit, family, xlabel):
    """Everything one plot needs: a band if there is one, traces if not."""
    b = band(rows, observable, subset,
             x="time_ns" if observable != "rmsf" else "index")
    if b is None and observable == "rmsf":
        b = band(rows, observable, subset, x="resid")
    lines = [] if b else per_replica(
        rows, observable, subset,
        scope="residue" if observable == "rmsf" else "frame")
    if b is None and not lines:
        return None
    for i, s in enumerate(lines):
        s["color"] = L.gradient(family, i, len(lines))
    return {
        "observable": observable, "subset": subset, "title": title,
        "unit": unit, "family": family, "xlabel": xlabel,
        "color": L.gradient(family, 0, 1),
        "band": b, "lines": lines,
        "pooled": b is not None,
        "note": (f"mean ± SD over {b['n']} protomer×replica units"
                 if b else "one curve per replica — not yet pooled"),
    }


def panels_for(rows, family):
    """Every panel a column can draw, skipping the ones with no rows."""
    out = []
    for observable, subset, short, unit, title in L.PANELS[family]:
        xlabel = "residue" if observable == "rmsf" else "time (ns)"
        p = panel(rows, observable, subset, title, unit, family, xlabel)
        if p:
            p["short"] = short
            out.append(p)
    return out


def contact_map(rows, top=35):
    """The interface contact map: receptor residues x ligand residues.

    Drawn from the POOLED rows (`scope="aggregate"`) and only the COGNATE
    pairs -- a receptor protomer with the ligand copy bound to it. The stage
    also measures cross-protomer contacts, where an Fv brushes the NEXT
    protomer; those are real but they are a different observable, and mixing
    them into one map would show a footprint no single Fv makes.

    Only the `top` most-occupied residues a side are drawn: 1,105 pooled pairs
    is not a figure, and its tail is near zero. The count is stated so nobody
    reads the picture as the whole interface.
    """
    pairs = [r for r in rows
             if r.get("stage") == "interface_contacts"
             and r.get("scope") == "aggregate"
             and r.get("observable") == "contact_occupancy"
             and r.get("subset") == "cognate"]
    if not pairs:
        return None
    rec_score, lig_score, cell = defaultdict(float), defaultdict(float), {}
    for r in pairs:
        v = _num(r.get("value")) or 0.0
        a = f"{r.get('resname') or ''}{r.get('resid') or ''}"
        b = f"{r.get('partner_resname') or ''}{r.get('partner_resid') or ''}"
        if not a.strip() or not b.strip():
            continue
        rec_score[a] += v
        lig_score[b] += v
        cell[(a, b)] = max(cell.get((a, b), 0.0), v)
    if not cell:
        return None
    rec = [k for k, _ in sorted(rec_score.items(), key=lambda kv: -kv[1])[:top]]
    lig = [k for k, _ in sorted(lig_score.items(), key=lambda kv: -kv[1])[:top]]
    num = lambda s: int("".join(c for c in s if c.isdigit()) or 0)
    rec.sort(key=num)
    lig.sort(key=num)
    z = [[cell.get((a, b), 0.0) for a in rec] for b in lig]
    return {"x": rec, "y": lig, "z": z, "n_pairs": len(cell),
            "note": f"{len(cell):,} pooled cognate residue pairs; the "
                    f"{len(rec)}×{len(lig)} most occupied are shown. Pooled over "
                    f"protomers and replicas; cross-protomer contacts are "
                    f"measured separately and not drawn here."}


def histogram(rows, subset, top=15):
    """Per-residue interface occupancy, pooled: the bars beside the map.

    Horizontal, fewest that say something (15), highest at the top, and
    coloured by the chain the residue belongs to. Twenty-five vertical bars
    with the labels turned 45 degrees was unreadable -- a residue name is a
    word, and words belong on the axis you can read along.
    """
    got = []
    for r in rows:
        if r.get("stage") != "interface_contacts":
            continue
        if r.get("scope") != "aggregate" or r.get("subset") != subset:
            continue
        if r.get("observable") != "residue_occupancy":
            continue
        v, sd = _num(r.get("value")), _num(r.get("sd"))
        if v is None:
            continue
        got.append((f"{r.get('resname') or ''}{r.get('resid') or ''}", v,
                    sd or 0.0, str(r.get("entity") or r.get("segid") or "")))
    if not got:
        return None
    got.sort(key=lambda g: -g[1])
    got = got[:top]
    got.reverse()                      # plotly draws the first bar at the bottom
    chains = sorted({g[3] for g in got})
    return {"labels": [g[0] for g in got], "values": [g[1] for g in got],
            "sd": [g[2] for g in got], "chains": [g[3] for g in got],
            "chain_list": chains,
            "note": f"the {len(got)} most-occupied {subset} residues, pooled over "
                    f"protomers and replicas; the bar is the fraction of frames "
                    f"in contact and the whisker the spread across the nine units"}


def _method_of(got):
    """Which end-state method produced these numbers, from its own settings.

    No row is labelled `method` -- the protocol is recorded as the parameters
    it ran with -- so the name is read from those: an implicit-solvent `igb`
    is Generalized Born, a `pb`/`istrng` setting is Poisson-Boltzmann. It
    matters on the page: the owner asked for "MM/PBSA or MM/GBSA or even FEP,
    depending on the method", and a table headed "binding energy" does not
    tell a reader which approximation the number carries.
    """
    params = {str(r.get("param_name")): str(r.get("param_value")) for r in got
              if r.get("param_name")}
    blob = " ".join(params) + " " + " ".join(params.values())
    igb = params.get("gb.igb")
    if igb or "gbsa" in blob.lower():
        return f"MM/GBSA (igb={igb})" if igb else "MM/GBSA"
    if "pb." in blob or "pbsa" in blob.lower():
        return "MM/PBSA"
    if "fep" in blob.lower() or "ti" in params:
        return "FEP"
    return params.get("method", "binding energy")


def binding_table(rows):
    """The binding free-energy table, with its method and its caveat."""
    got = [r for r in rows if r.get("stage") == "binding_energy"]
    if not got:
        return None
    method = _method_of(got)
    agg = [r for r in got if r.get("scope") == "aggregate"]
    # DISTINCT runs, not rows: there are fifteen components per run, and
    # counting rows said "pooled over 135 runs" where there are nine.
    per = {(r.get("protomer"), r.get("replica")) for r in got
           if r.get("scope") == "run"}
    caveats = [str(r.get("floor_note")) for r in got
               if r.get("floor_note") and r.get("resolved") in ("False", False)]
    # The total first, then the terms that make it, then the rest. A table
    # sorted alphabetically opens on `1-4_EEL` and buries the one number a
    # reader came for.
    order = ["DELTA_TOTAL", "VDWAALS", "EEL", "EGB", "ESURF",
             "DELTA_G_GAS", "DELTA_G_SOLV"]
    rank = {name: i for i, name in enumerate(order)}
    rowset, dropped = [], []
    for r in agg:
        name = str(r.get("subset"))
        v, sd = _num(r.get("value")), _num(r.get("sd"))
        # A bonded term is identically zero in a single-trajectory protocol --
        # complex, receptor and ligand come from the SAME frames, so those
        # terms cancel exactly. Printing eight zeros says only that the
        # protocol was what it says it was; the count is kept instead.
        if (v == 0 or v is None) and not sd:
            dropped.append(name)
            continue
        rowset.append({"component": name, "value": v, "sd": sd,
                       "sem": _num(r.get("sem")), "n": r.get("n"),
                       "unit": r.get("unit") or "kcal/mol",
                       "total": name == "DELTA_TOTAL"})
    rowset.sort(key=lambda r: (rank.get(r["component"], 99), r["component"]))
    # One caveat is read; four are skipped. The shortest honest one leads and
    # the rest fold away -- a wall of warning text is not more honest than a
    # sentence, it is just unread.
    caveats = sorted(set(caveats), key=len)
    return {"method": method, "rows": rowset, "per_run": len(per),
            "dropped_zero": sorted(dropped),
            "lead_caveat": (caveats[0] if caveats else None),
            "caveats": caveats[1:]}


# --------------------------------------------------------------------------
# the solvent section: measured for the whole system, so it belongs to neither
# column and sits full width beneath them
# --------------------------------------------------------------------------
SOLVENT_COLORS = ["#2f6f8f", "#3d9dc4", "#59b7d6", "#8fd3e8",
                  "#c99700", "#e0a72c", "#7a5195", "#a06cb8"]


def rdf_curves(rows, max_pairs=8):
    """g(r) for every pair the RDF stage measured, as one plot."""
    by = defaultdict(lambda: {"r": [], "g": []})
    for r in rows:
        if r.get("stage") != "rdf" or r.get("scope") != "bin":
            continue
        if r.get("observable") not in ("g_r", "g", "rdf"):
            continue
        x = _num(r.get("r_A")) or _num(r.get("scope_id")) or _num(r.get("label"))
        v = _num(r.get("value"))
        if x is None or v is None:
            continue
        by[str(r.get("subset"))]["r"].append(x)
        by[str(r.get("subset"))]["g"].append(v)
    if not by:
        return None
    out = []
    for i, (name, d) in enumerate(sorted(by.items())[:max_pairs]):
        order = sorted(range(len(d["r"])), key=lambda k: d["r"][k])
        out.append({"name": name.replace("_around_", " around "),
                    "x": [d["r"][k] for k in order],
                    "y": [d["g"][k] for k in order],
                    "color": SOLVENT_COLORS[i % len(SOLVENT_COLORS)]})
    return out


def solvent_table(rows):
    """The solvent numbers worth a table, each with what it was measured over."""
    want = [
        ("water_shell", "n_interfacial_mean", None, "Interfacial waters", ""),
        ("water_bridges", "wbridge_count_mean", None, "Water bridges", ""),
        ("ions", "coordination_number", "SOD", "Na⁺ coordination", "waters"),
        ("ions", "coordination_number", "CLA", "Cl⁻ coordination", "waters"),
        ("rdf", "bulk_density", "waterO", "Bulk water density", "Å⁻³"),
    ]
    out = []
    for stage, obs, sub, label, unit in want:
        vals = [_num(r.get("value")) for r in rows
                if str(r.get("stage", "")).startswith(stage)
                and r.get("observable") == obs
                and (sub is None or r.get("subset") == sub)
                and r.get("scope") in ("run", "aggregate")
                and (obs != "coordination_number" or r.get("label") == "first_shell")]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        mean = sum(vals) / len(vals)
        sd = (sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5 \
            if len(vals) > 1 else None
        out.append({"label": label, "value": mean, "sd": sd,
                    "n": len(vals), "unit": unit})
    return out or None


def grid_files(rows):
    """The .dx volumes, with the peak each one holds and whether it resolves."""
    out = []
    for r in rows:
        if r.get("stage") != "density_grids":
            continue
        if r.get("observable") != "grid_peak_density":
            continue
        out.append({"species": str(r.get("subset")),
                    "replica": str(r.get("replica")),
                    "peak": _num(r.get("value")),
                    "file": str(r.get("label") or ""),
                    "resolved": r.get("resolved") in ("True", True),
                    "note": str(r.get("floor_note") or "")})
    out.sort(key=lambda d: (d["species"], d["replica"]))
    return out or None
