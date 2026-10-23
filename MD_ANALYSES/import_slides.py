#!/usr/bin/env python3
"""
import_slides.py -- take the pictures out of a deck and file them by system.

    python3 import_slides.py <folder> <deck.pptx> [--dry-run]

A .pptx is a zip: `ppt/slides/slideN.xml` holds the text, its `_rels` file
names the images, and `ppt/media/` holds them. Each slide is matched to a
record by the system name in its TITLE -- which is what the owner asked for --
and its images are copied into that record's `figures_user/` with the slide's
own labels as filenames.

WHAT IT REFUSES TO DO QUIETLY. A slide whose title names a different state
from the record it would land in is reported and, unless `--force`, skipped:
the deck this was written for says "tri-bound/Occluded-Open" on its antibody
slides while the records are the OPEN state, and a picture of one state filed
silently against another is worse than no picture. `credits.yaml` records for
every image which deck, which slide and which title it came from, so a reader
can always get back to the source.
"""

import argparse
import os
import re
import shutil
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

#: Tokens that identify a system in a slide title, in the order they are
#: tried. Deliberately explicit: a fuzzy match would put the wrong render on
#: the wrong record, and nothing downstream would catch it.
SYSTEM_TOKENS = [
    ("1F4", ("1f4", "if4")),
    ("JL427", ("jl427", "jl-427")),
    ("YZ23", ("yz23", "yz-23")),
    ("apo", ("apo", "unbound", "open state", "env only")),
]

#: States a title may name, so a mismatch can be seen rather than assumed away.
STATE_TOKENS = {
    "open": ("open state", "openstate"),
    "occluded-open": ("occluded-open", "occluded open", "tri-bound"),
    "closed": ("closed",),
}


EMU = 9525          # English Metric Units per pixel at 96 dpi


def _shapes(xml):
    """Every text box and picture of a slide, with where it sits.

    Position is the whole point: a deck says which picture is which by putting
    the words under the picture, not by the order the XML lists them. Matching
    by order gave slide 1's second render the label of the first.
    """
    texts, pics = [], []
    for m in re.finditer(r"<p:sp>.*?</p:sp>", xml, re.S):
        blk = m.group(0)
        t = " ".join(
            "".join(re.findall(r"<a:t>([^<]*)</a:t>", para)).strip()
            for para in re.findall(r"<a:p>(.*?)</a:p>", blk, re.S)).strip()
        off = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"', blk)
        if not t:
            continue
        texts.append({"text": t,
                      "x": int(off.group(1)) // EMU if off else 0,
                      "y": int(off.group(2)) // EMU if off else 0,
                      "is_title": "type=\"title\"" in blk or "type=\"ctrTitle\"" in blk,
                      "is_number": "type=\"sldNum\"" in blk})
    for m in re.finditer(r"<p:pic>.*?</p:pic>", xml, re.S):
        blk = m.group(0)
        rid = re.search(r'r:embed="([^"]+)"', blk)
        off = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"', blk)
        ext = re.search(r'<a:ext cx="(\d+)" cy="(\d+)"', blk)
        if not rid:
            continue
        x = int(off.group(1)) // EMU if off else 0
        y = int(off.group(2)) // EMU if off else 0
        w = int(ext.group(1)) // EMU if ext else 0
        h = int(ext.group(2)) // EMU if ext else 0
        pics.append({"rid": rid.group(1), "x": x, "y": y, "w": w, "h": h,
                     "cx": x + w / 2, "cy": y + h / 2})
    pics.sort(key=lambda p: (p["x"], p["y"]))
    return texts, pics


def _legend(texts):
    """Labels stacked at one x are a legend for the figure, not captions.

    On the antibody slides "Heavy chain" and "Light chain" sit at the same x,
    one above the other: they name the COLOURS, and giving one to each picture
    would caption two renders with halves of a key.
    """
    out = []
    by_x = {}
    for t in texts:
        by_x.setdefault(round(t["x"] / 40), []).append(t)
    for group in by_x.values():
        if len(group) >= 2:
            group.sort(key=lambda t: t["y"])
            out.extend(group)
    return out


def slide_size(z):
    """(width, height) of a slide in px, from the deck itself."""
    try:
        xml = z.read("ppt/presentation.xml").decode("utf8", "ignore")
        m = re.search(r'<p:sldSz cx="(\d+)" cy="(\d+)"', xml)
        if m:
            return int(m.group(1)) // EMU, int(m.group(2)) // EMU
    except KeyError:
        pass
    return 1280, 720


def slide_parts(z, n, size=None):
    """{title, captions per picture, legend, all text} for one slide."""
    width, height = size or slide_size(z)
    xml = z.read(f"ppt/slides/slide{n}.xml").decode("utf8", "ignore")
    texts, pics = _shapes(xml)
    rels = z.read(f"ppt/slides/_rels/slide{n}.xml.rels").decode("utf8", "ignore")
    media = dict(re.findall(r'Id="([^"]+)"[^>]*Target="\.\./media/([^"]+)"', rels))

    # The title is the title PLACEHOLDER, never "the first paragraph": slide 1's
    # first paragraph is its slide number, which is how every entry it wrote
    # ended up titled "1".
    title = next((t["text"] for t in texts if t["is_title"]), "")
    legend = _legend([t for t in texts if not t["is_title"] and not t["is_number"]])
    legend_ids = {id(t) for t in legend}
    labels = [t for t in texts
              if not t["is_title"] and not t["is_number"]
              and id(t) not in legend_ids and len(t["text"]) < 40]

    # A HEADING lives in the top band of the SLIDE, not merely above the
    # pictures: on the antibody slides a render starts at y = -2, so "above
    # every picture" caught nothing and the slide's own title was captioning
    # the left-hand render. The band is a fraction of the slide's real height.
    band = height * 0.12
    heading = [t for t in labels if t["y"] < band]
    labels = [t for t in labels if t["y"] >= band]
    if not title and heading:
        title = sorted(heading, key=lambda t: t["y"])[0]["text"]

    # WHICH PICTURE A LABEL BELONGS TO. Not the order it is listed in -- a
    # deck says which picture is which by where the words sit. The obvious
    # rule, "the label inside this picture's x-span", is not enough on its
    # own: slide 1's renders OVERLAP (33..674 and 599..1356), so "Side view"
    # at x = 636 falls inside both and the first picture claimed it.
    #
    # A caption is anchored at the same left edge as the picture it annotates,
    # so the label goes to the picture whose LEFT EDGE it is nearest, with
    # containment preferred over a picture it does not overlap at all. Every
    # label is assigned once, and each picture then keeps whichever of its
    # labels sits closest to its top or bottom edge.
    def contains(p, t):
        return p["x"] - 50 <= t["x"] <= p["x"] + max(p["w"], 1) + 50

    mine = {p["rid"]: [] for p in pics}
    for t in labels:
        if not pics:
            break
        owner = min(pics, key=lambda p: (0 if contains(p, t) else 1,
                                         abs(t["x"] - p["x"])))
        mine[owner["rid"]].append(t)
    captions = {}
    for p in pics:
        edge = lambda t: min(abs(t["y"] - p["y"]),
                             abs(t["y"] - (p["y"] + p["h"])))
        best = sorted(mine[p["rid"]], key=lambda t: (edge(t), t["x"]))
        captions[p["rid"]] = best[0]["text"] if best else ""

    return {"title": title, "heading": [t["text"] for t in heading],
            "pics": [dict(p, media=media.get(p["rid"], ""),
                          caption=captions.get(p["rid"], "")) for p in pics],
            "legend": [t["text"] for t in legend],
            "text": [t["text"] for t in texts]}


def state_of(text):
    low = " ".join(text).lower()
    # Longest token first: "occluded-open" contains "open".
    for state, tokens in sorted(STATE_TOKENS.items(),
                                key=lambda kv: -max(len(t) for t in kv[1])):
        if any(t in low for t in tokens):
            return state
    return None


def system_of(text):
    low = " ".join(text).lower()
    for name, tokens in SYSTEM_TOKENS:
        if any(t in low for t in tokens):
            return name
    return None


def read_deck(path):
    """Every slide as {n, title, pics, legend, text, system, state}."""
    out = []
    with zipfile.ZipFile(path) as z:
        ns = sorted(int(re.search(r"\d+", os.path.basename(m)).group())
                    for m in z.namelist()
                    if re.match(r"ppt/slides/slide\d+\.xml$", m))
        for n in ns:
            parts = slide_parts(z, n)
            text = parts["text"]
            out.append({"n": n, **parts, "system": system_of(text),
                        "state": state_of(text)})
    return out


def _label(caption, i, n=1):
    """A filename from the picture's own caption, or from where it sits.

    A render with no caption of its own still needs a name a person can read:
    the antibody slides put two views side by side under one shared legend,
    and `figure-2.png` says nothing that `view-right.png` does not say better.
    """
    slug = re.sub(r"[^a-z0-9]+", "-", (caption or "").lower()).strip("-")
    if slug:
        return slug
    if n == 2:
        return ("view-left", "view-right")[i]
    return f"view-{i + 1}"


def record_for(system, records):
    """The bundle a system's pictures belong in."""
    for stem in records:
        low = stem.lower()
        if system.lower() in low.replace("-", "").replace("_", ""):
            return stem
        if system == "apo" and low.endswith("_md_3x500ns"):
            return stem
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("deck")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="import even when the slide names a different state")
    a = ap.parse_args(argv)

    from analyse import bundle_dir, load_records
    folder = os.path.abspath(a.folder)
    records = sorted(load_records(folder))
    if not records:
        print(f"\n  No records under {folder}/MD_RESULTS.\n")
        return 1

    # What state are these records? Read it from the folder, not assumed.
    folder_state = state_of([os.path.basename(folder),
                             os.path.basename(os.path.dirname(folder))]) or "open"
    slides = read_deck(a.deck)
    print(f"\n  deck   : {os.path.basename(a.deck)}  ({len(slides)} slides)")
    print(f"  folder : {folder}  (records read as state '{folder_state}')\n")

    imported, skipped = 0, 0
    with zipfile.ZipFile(a.deck) as z:
        for s in slides:
            stem = record_for(s["system"], records) if s["system"] else None
            note = ""
            if not stem:
                note = "no record matches this title"
            elif s["state"] and s["state"] != folder_state:
                note = (f"slide says '{s['state']}', these records are "
                        f"'{folder_state}'")
            print(f"  slide {s['n']}  {s['title'][:44]:<44} -> "
                  f"{stem or '(unmatched)'}")
            if note:
                print(f"           {note}"
                      f"{'' if (a.force and stem) else ' -- SKIPPED'}")
            if not stem or (note and not a.force):
                skipped += len(s["pics"])
                continue

            out_dir = os.path.join(bundle_dir(folder, stem), "figures_user")
            os.makedirs(out_dir, exist_ok=True)
            lines = []
            for i, pic in enumerate(s["pics"]):
                if not pic["media"]:
                    continue
                ext = os.path.splitext(pic["media"])[1] or ".png"
                name = f"{_label(pic['caption'], i, len(s['pics']))}{ext}"
                dest = os.path.join(out_dir, name)
                if not a.dry_run:
                    with z.open(f"ppt/media/{pic['media']}") as src, \
                            open(dest, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                size = z.getinfo(f"ppt/media/{pic['media']}").file_size
                lines.append({"file": name, "caption": pic["caption"],
                              "order": i, "bytes": size})
                print(f"           {name:<22} {pic['caption'] or '(no caption)':<16}"
                      f" {size / 1e6:.1f} MB")
                imported += 1
            if s["legend"]:
                print(f"           legend: {' / '.join(s['legend'])}")
            if not a.dry_run:
                _write_credits(out_dir, a.deck, folder_state, lines, s)

    print(f"\n  {imported} image(s) imported, {skipped} skipped.\n")
    return 0


def _credited_files(path, deck):
    """Filenames a previous run imported from THIS deck, per credits.yaml."""
    got, current, deck_name = [], None, os.path.basename(deck)
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("- file:"):
                    current = line.split(":", 1)[1].strip()
                elif line.startswith("from:") and current:
                    if line.split(":", 1)[1].strip() == deck_name:
                        got.append(current)
                    current = None
    except OSError:
        pass
    return got


def _write_credits(out_dir, deck, folder_state, lines, slide):
    """Where every picture came from, beside the pictures.

    REWRITTEN, not appended: running the importer twice used to leave two
    entries per image (the apo record ended with four for its two pictures),
    and a credits file that grows on every run is one nobody can trust.
    """
    path = os.path.join(out_dir, "credits.yaml")
    # A rewritten credits file must also retire the pictures the PREVIOUS run
    # wrote. An early version of the caption matching produced
    # `heavy-chain.png` and `416-433-epitope.png`; correcting it left those
    # files sitting beside the right ones, unlisted and indistinguishable
    # from a figure a person had put there by hand. Only files this deck
    # imported before are removed -- never anything the importer did not write.
    keep = {ln["file"] for ln in lines}
    for old in _credited_files(path, deck):
        if old not in keep:
            try:
                os.remove(os.path.join(out_dir, old))
            except OSError:
                pass
    head = ("# Pictures supplied by the person who ran the system, imported\n"
            "# from a deck. Each entry says which slide and which caption it\n"
            "# came from, so a figure can always be traced back.\n"
            "images:\n")
    body = ""
    for ln in lines:
        body += (f"  - file: {ln['file']}\n"
                 f"    from: {os.path.basename(deck)}\n"
                 f"    slide: {slide['n']}\n"
                 f"    title: \"{slide['title']}\"\n"
                 f"    caption: \"{ln['caption']}\"\n"
                 f"    legend: \"{' / '.join(slide['legend'])}\"\n"
                 f"    order: {ln['order']}\n"
                 f"    state_in_title: {slide['state'] or 'unstated'}\n"
                 f"    records_state: {folder_state}\n")
    with open(path, "w") as fh:
        fh.write(head + body)
    return path


if __name__ == "__main__":
    sys.exit(main())
