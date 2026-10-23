#!/usr/bin/env python3
"""
captions.py -- one sentence per analysis, carrying its measured number.

The working page (`report.py`) is read by people who already know what they
are looking at. The story page is not: it goes on the lab website, where a
reader meets the system for the first time and has to be told what a curve
means before the curve is worth showing. So every plot gets a sentence.

TWO RULES, and they are the whole design.

1. THE NUMBER IN THE SENTENCE IS READ FROM THE SAME ROWS THE PLOT DRAWS.
   Not retyped, not rounded from memory, not passed in by the caller. A
   caption that says "settles at 4.2 Å" and a plot that shows 5.3 is worse
   than no caption at all, and the only way to make that impossible is for
   both to come from one table. This module therefore takes `rows` and
   nothing else. (It has been needed once already: a receptor RMSD was
   quoted from a 60-frame slice and was off by more than a third.)

2. WHERE A STAGE SAID IT COULD NOT RESOLVE SOMETHING, THE CAPTION SAYS SO IN
   WORDS AND DOES NOT QUOTE THE NUMBER. At the spacing these runs sample,
   a water residence
   time is not a short residence time, it is an unmeasured one; the tables
   carry `resolved=False` for exactly this, and a caption is where a reader
   who will never open a TSV finds out.

Captions are plain text with one piece of markup: `**...**` marks the
headline number, which the page renders bold. Keeping them free of HTML means
they can be tested as strings and reused in a record, a figure legend or an
email without carrying a page's markup along with them.
"""

from . import series as S
from core import io as _io

_num = S._num


# --------------------------------------------------------------------------
# reading the rows
# --------------------------------------------------------------------------
def _runs(rows, stage, observable, subset=None, scope="run"):
    """Every per-run value of one observable, as floats."""
    out = []
    for r in rows:
        if r.get("stage") != stage or r.get("observable") != observable:
            continue
        if r.get("scope") != scope:
            continue
        if subset is not None and r.get("subset") != subset:
            continue
        v = _num(r.get("value"))
        if v is not None:
            out.append(v)
    return out


def _agg(rows, stage, observable, subset=None):
    """One pooled row: (value, sd, n), or (None, None, None)."""
    for r in rows:
        if r.get("stage") != stage or r.get("observable") != observable:
            continue
        if r.get("scope") != "aggregate":
            continue
        if subset is not None and r.get("subset") != subset:
            continue
        return _num(r.get("value")), _num(r.get("sd")), r.get("n")
    return None, None, None


def _mean_sd(vals):
    if not vals:
        return None, None
    m = sum(vals) / len(vals)
    if len(vals) < 2:
        return m, None
    return m, (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5


def _unresolved(rows, stage):
    """The floor note a stage left, when it marked its own numbers unusable."""
    for r in rows:
        if r.get("stage") == stage and r.get("resolved") in ("False", False):
            note = r.get("floor_note")
            if note:
                return str(note)
    return None


# --------------------------------------------------------------------------
# the shape of a curve, described honestly
# --------------------------------------------------------------------------
def _plateau(band, tail=0.25):
    """Where a curve ends up, and when it got there.

    The last quarter is taken as the plateau and the settling time is the
    first point after which the mean never leaves that plateau's own SD
    again. Two things this deliberately does NOT do: it does not fit an
    exponential (nine units of 500 ns do not justify one), and it does not
    claim a plateau when the curve is still climbing -- `rising` is returned
    instead, and the caption then says the trimer is still moving, which is a
    result rather than a failure to report one.
    """
    y, x, sd = band["mean"], band["x"], band["sd"]
    n = len(y)
    if n < 8:
        return None
    k = max(2, int(n * tail))
    end, end_sd = _mean_sd(y[-k:])
    spread = max(sum(sd[-k:]) / k, end_sd or 0.0)
    mid, _ = _mean_sd(y[-2 * k:-k]) if n >= 2 * k else (end, None)
    # "Still rising" means the previous window sits more than its own spread
    # below the last one -- a drift larger than the disagreement between the
    # nine units is the only kind this sampling can honestly call a trend.
    rising = mid is not None and (end - mid) > max(spread, 1e-9)
    settled = None
    tol = max(spread, 0.15)
    for i in range(n):
        if all(abs(v - end) <= tol for v in y[i:]):
            settled = x[i]
            break
    return {"value": end, "spread": spread, "settled_ns": settled,
            "rising": rising, "start": y[0], "n": band.get("n")}


def _about(ns):
    """A settling time as a round number, because "about" means about.

    The curve is sampled every 0.2 ns, so the first point that happens to
    stay inside the band is not accurate to a tenth of a nanosecond; writing
    "about 173.8 ns" claims a precision the measurement does not have.
    """
    if ns is None:
        return None
    step = 100 if ns >= 250 else (50 if ns >= 100 else 10)
    return _io.human_ns(max(step, round(ns / step) * step))


def _fmt(v, dp=2):
    return "—" if v is None else f"{v:,.{dp}f}"


def _pm(v, sd, unit, dp=2):
    core = f"{_fmt(v, dp)}"
    # A spread that rounds to zero at the precision being printed is not a
    # spread worth printing: "0.0312 ± 0.0000" says the replicas differ by
    # something too small to show, which is better said by not saying it.
    if sd and round(sd, dp) != 0:
        core += f" ± {_fmt(sd, dp)}"
    return f"**{core} {unit}**" if unit else f"**{core}**"


# --------------------------------------------------------------------------
# one caption per analysis
# --------------------------------------------------------------------------
def for_panel(rows, panel):
    """The caption under an RMSD, RMSF or Rg plot."""
    obs, subset = panel["observable"], panel["subset"]
    what = "antibody" if panel["family"] == "ligand" else "trimer"
    band = panel.get("band")
    units = int(band["n"]) if band and band.get("n") else None
    over = (f"the line is the mean over {units} protomer × replica units and "
            f"the band is their spread — how much those {units} disagree, not "
            f"a confidence interval on the mean"
            if units else "one curve per replica; these are not pooled yet")

    if obs == "rmsf" and band:
        peaks = sorted(zip(band["mean"], band["labels"]), reverse=True)
        top = ", ".join(lbl for _, lbl in peaks[:3] if lbl)
        med = sorted(band["mean"])[len(band["mean"]) // 2]
        return (f"Per-residue fluctuation about the mean structure. Most of "
                f"the {what} moves by {_pm(med, None, 'Å')} or less; the "
                f"loosest residues are {top or 'at the termini'}, reaching "
                f"{_pm(peaks[0][0], None, 'Å')}. {over[0].upper() + over[1:]}.")

    if band:
        p = _plateau(band)
        if p:
            if p["rising"]:
                shape = (f"has not stopped moving by the end of the runs: it "
                         f"reaches {_pm(p['value'], p['spread'], 'Å')} and is "
                         f"still climbing")
            elif p["settled_ns"] is not None:
                shape = (f"settles at {_pm(p['value'], p['spread'], 'Å')} after "
                         f"about {_about(p['settled_ns'])} and stays there")
            else:
                shape = f"sits at {_pm(p['value'], p['spread'], 'Å')}"
            name = ("Backbone RMSD" if obs == "rmsd" else "Radius of gyration")
            fit = (" The fit is on the receptor protomer the antibody is bound "
                   "to, so this is the antibody's motion relative to its "
                   "epitope rather than its motion in the box."
                   if obs == "rmsd" and panel["family"] == "ligand" else "")
            return (f"{name} against the first frame. The {what} {shape}. "
                    f"{over[0].upper() + over[1:]}.{fit}")

    lines = panel.get("lines") or []
    tail = [y for s_ in lines for y in s_["y"][int(len(s_["y"]) * 0.75):]]
    if tail:
        m, sd = _mean_sd(tail)
        name = ("Radius of gyration" if obs == "rg" else panel["short"])
        return (f"{name} of the {what}, against time. Over the last quarter of "
                f"each run it is {_pm(m, sd, panel['unit'])}, averaged over "
                f"{len(lines)} replicas — a value that does not drift is the "
                f"sign of an assembly that stays folded rather than swelling "
                f"or collapsing. {over[0].upper() + over[1:]}.")
    return f"{panel['title']}. {over[0].upper() + over[1:]}."


def contact_map(rows, cmap=None):
    """The caption under the interface heat map."""
    pairs, _, _ = _agg(rows, "interface_contacts", "contact_pairs")
    per_frame, pf_sd, _ = _agg(rows, "interface_contacts", "contacts_per_frame")
    rec, _, _ = _agg(rows, "interface_contacts", "interface_residues", "receptor")
    lig, _, _ = _agg(rows, "interface_contacts", "interface_residues", "ligand")
    if per_frame is None and pairs is None:
        return ""
    bits = []
    if rec and lig:
        bits.append(f"{int(rec)} receptor residues meet {int(lig)} antibody "
                    f"residues")
    if per_frame:
        bits.append(f"with {_pm(per_frame, pf_sd, 'contacts', 1)} closer than "
                    f"the cutoff in a typical frame")
    lead = ("Which residues touch, and how often. "
            + ", ".join(bits) + ". ") if bits else ""
    tail = ("Colour is the fraction of frames a pair spends in contact, pooled "
            "over protomers and replicas. Only pairs on the SAME protomer and "
            "its own antibody copy are drawn: an Fv also brushes the next "
            "protomer, which is real but is a different measurement.")
    dyn = _unresolved(rows, "interface_contacts")
    if dyn:
        frac, _, _ = _agg(rows, "interface_contacts",
                          "contact_lifetime_single_frame_fraction")
        tail += (" How LONG a contact lasts is not measurable here: "
                 + (f"{frac * 100:.0f}% of contact episodes last a single "
                    f"frame at this sampling, " if frac else "")
                 + "so the map says which contacts are occupied, not how fast "
                   "they exchange.")
    return lead + tail


def histogram(rows, subset):
    """The caption under one of the two per-residue occupancy plots."""
    side = "receptor" if subset == "receptor" else "antibody"
    n, _, _ = _agg(rows, "interface_contacts", "interface_residues", subset)
    tops = [(_num(r.get("value")), f"{r.get('resname') or ''}{r.get('resid') or ''}")
            for r in rows
            if r.get("stage") == "interface_contacts"
            and r.get("scope") == "aggregate"
            and r.get("observable") == "residue_occupancy"
            and r.get("subset") == subset
            and _num(r.get("value")) is not None]
    if not tops:
        return ""
    tops.sort(reverse=True)
    v, name = tops[0]
    return (f"The {side} residues that spend the most time in contact, of "
            f"{int(n) if n else len(tops)} that touch at all. {name} leads at "
            f"{_pm(v, None, 'of frames')}; the whisker is the spread across "
            f"the nine protomer × replica units, so a short bar with a long "
            f"whisker is a contact one protomer makes and another does not.")


def binding(rows, table=None):
    """The caption under the binding free-energy table."""
    table = table or S.binding_table(rows)
    if not table:
        return ""
    total = next((r for r in table["rows"] if r.get("total")), None)
    if not total:
        return ""
    lead = (f"The antibody binds with {_pm(total['value'], total['sd'], 'kcal/mol', 1)}"
            f" by {table['method']}, pooled over {table['per_run']} protomer × "
            f"replica runs.")
    return (lead + " SD is how much those runs disagree; SEM is taken over "
            "replicas only, because the three protomers of one trimer share a "
            "box and a thermostat and are not independent samples. These are "
            "end-state energies from a single trajectory: they rank binders "
            "and they do not predict an affinity.")


def solvent(rows):
    """The caption under the RDF curves and the solvent table."""
    bulk, bulk_sd = _mean_sd(_runs(rows, "rdf", "bulk_density", "waterO"))
    pk = _runs(rows, "rdf", "first_peak_r", "waterO_around_protein_polar") \
        or _runs(rows, "rdf", "first_peak_r", "waterO_around_protein")
    peak_r, _ = _mean_sd(pk)
    wat, wat_sd = _mean_sd(_runs(rows, "water_bridges",
                                 "n_interfacial_water_mean", "all"))
    br, br_sd = _mean_sd(_runs(rows, "water_bridges", "wbridge_count_mean", "all"))
    bits = []
    if peak_r:
        bits.append(f"water sits {_pm(peak_r, None, 'Å')} from the protein's "
                    f"polar atoms — a hydrogen bond, measured rather than "
                    f"assumed")
    if bulk:
        bits.append(f"and the box reaches bulk density, "
                    f"{_pm(bulk, bulk_sd, 'Å⁻³', 4)}, away from the solute")
    if wat:
        bits.append(f"{_pm(wat, wat_sd, 'waters', 0)} sit in the interface at "
                    f"any moment")
    if br:
        bits.append(f"{_pm(br, br_sd, 'of them', 0)} bridging the two sides at "
                    f"once")
    if not bits:
        return ""
    body = "g(r) around each part of the system: " + ", ".join(bits[:2]) + ". "
    if len(bits) > 2:
        rest = ", ".join(bits[2:])
        body += rest[0].upper() + rest[1:] + ". "
    return body + _residence_note(rows)


def _spacing_ns(rows, stages):
    """Nanoseconds between the frames a stage actually sampled.

    `ns_per_frame` times `stride`, from the rows themselves. The sentence
    below used to say "at 0.2 ns per frame" whatever the system was -- true
    of the CDR campaign, and wrong by a factor of two for the occluded-open
    runs, which write a frame every 100 ps. A hardcoded number in a caption
    is the exact thing this module exists to prevent.
    """
    got = set()
    for r in rows:
        if r.get("stage") not in stages:
            continue
        ns, stride = _num(r.get("ns_per_frame")), _num(r.get("stride"))
        if ns:
            got.add(round(ns * (stride or 1), 4))
    return max(got) if got else None


def _residence_note(rows):
    """Why no residence time is quoted -- in this system's own numbers.

    The claim is only as strong as what was measured here: if a stage
    reported what fraction of episodes lasted a single frame, that figure is
    quoted; if none did, the sentence says the quantity was not measured and
    stops, rather than borrowing another system's evidence for it.
    """
    stages = ("water_residence", "water_shell", "water_bridges", "ions")
    fracs = [_num(r.get("value")) for r in rows
             if r.get("stage") in stages
             and str(r.get("observable", "")).endswith("single_frame_fraction")
             and _num(r.get("value")) is not None]
    spacing = _spacing_ns(rows, stages)
    at = f"at {spacing:g} ns between sampled frames " if spacing else ""
    if fracs:
        worst = max(fracs)
        return (f"Residence and exchange times are NOT reported: {at}"
                f"{worst * 100:.0f}% of episodes last a single frame, so this "
                f"sampling cannot resolve how long a water or an ion stays.")
    return ("Residence and exchange times are NOT reported: nothing here "
            "measured them at a spacing fine enough to resolve them"
            + (f", {at.strip()}" if at else "") + ".")


def system(rows, setup=None, bundle=None, components=()):
    """The opening paragraph: what this is, what is in the box, how it was run.

    Composed from the detected components, the counted composition and the
    record's `setup` block -- not from the record's prose `description`, which
    opens by restating the atom count the second sentence gives properly.
    Three short sentences, in the order someone meets a system: what it is,
    what surrounds it, how it was simulated.
    """
    from . import layout as L
    table = dict((lbl, v) for lbl, v, sd, n in
                 (L.composition_single_column(rows)))
    # The detector's vocabulary: `protein, glycans, water, ions, antibody,
    # membrane`. What surrounds the solute is NOT assumed: the full-length
    # gp120 trajectories are membrane-embedded and carry no water at all, and
    # an earlier version called them "solvated and neutralised" anyway --
    # a sentence asserting something nothing had measured, on a page whose
    # whole point is that it never does that.
    subject = ("a glycosylated protein assembly" if "glycans" in components
               else "a protein assembly")
    bound = [phrase for key, phrase in (("antibody", "an antibody bound to it"),
                                        ("ligand", "a bound ligand"))
             if key in components]
    first = ""
    if "protein" in components or bound:
        first = f"The system is {subject}"
        if bound:
            first += ", with " + " and ".join(bound)
        where = []
        if "membrane" in components:
            where.append("embedded in a lipid bilayer")
        # "Solvated" only when the water is there to be counted -- the next
        # sentence is about to count it.
        if table.get("Water molecules (TIP3P)"):
            where.append("solvated and neutralised")
        if where:
            first += ", " + " and ".join(where)
        first += ". "

    atoms = table.get("Atoms, total")
    solute = table.get("Solute atoms")
    waters = table.get("Water molecules (TIP3P)")
    box = table.get("Box X \u00d7 Y \u00d7 Z (\u00c5)")
    salt = table.get("Salt (mol/L)")
    bits = []
    if atoms:
        bits.append(f"**{atoms:,} atoms**")
    if solute:
        bits.append(f"of which {solute:,} are solute")
    if waters:
        bits.append(f"in {waters:,} waters")
    if salt:
        bits.append(f"at {salt:.2f} mol/L salt")
    second = ("The box holds " + ", ".join(bits) + ". ") if bits else ""
    if box:
        second += f"It measures {box} \u00c5. "

    setup = setup or {}
    run = []
    if setup.get("force_field"):
        run.append(str(setup["force_field"]))
    if setup.get("ensemble"):
        run.append(f"{setup['ensemble']} at {setup.get('temperature_K')} K"
                   if setup.get("temperature_K") else str(setup["ensemble"]))
    if setup.get("timestep_fs"):
        run.append(f"{setup['timestep_fs']:g} fs steps")
    third = ""
    if run:
        third = ("Run with " + ", ".join(run)
                 + (f" in {setup['software']}" if setup.get("software") else "")
                 + ". ")
    return first + second + third
