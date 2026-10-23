#!/usr/bin/env python3
"""
selections.py -- selections written the way this lab writes them.

Everyone here types VMD/TCL selection strings, and every convention worth
keeping is already expressed in that language. So that is the input format, and
this module translates it to MDAnalysis and -- more importantly -- CHECKS it
against the actual topology before a single trajectory frame is read.

WHY THE CHECK IS THE POINT. The lab's own TCL scripts carry the warning "pick
the line matching your complex; the wrong one silently gives zero contacts". A
selection that matches nothing does not raise: it produces an empty analysis,
plausible-looking zeros, and a conclusion. Here, an empty selection is an error,
and every resolved selection reports its atom count so a person can see that
`receptor` really is 29,532 atoms and not 0.

TRANSLATION NOTES (each one is a real difference, not a cosmetic one):

  noh, not hydrogen  ->  not name H*        heavy atoms by NAME, never by mass:
                                            these topologies have broken masses
  within X of SEL    ->  (around X (SEL)) or (SEL)
                                            VMD's `within` INCLUDES the
                                            selection itself; MDAnalysis's
                                            `around` excludes it. Dropping this
                                            silently loses every atom of the
                                            reference group.
  exwithin X of SEL  ->  around X (SEL)     the exclusive form, which is what
                                            `around` already means
  chain              ->  chainID            and see the warning below
  segid A to C       ->  segid A B C        numeric/alphabetic suffix ranges are
                                            expanded, because MDAnalysis has no
                                            range syntax for segid
  resid 1 to 10      ->  resid 1:10
  water / ions /     ->  explicit resname lists from core.topology, so a VMD
  glycan               macro and this module pick the SAME atoms

`chainID` is accepted but warned about: in these systems it collapses all three
antibody copies onto one letter and every glycan onto `G`, which is exactly the
mistake `ANALYSIS/config.py` documents. Prefer `segid`.
"""

import re

from . import topology as T

HEAVY = "not name H*"

#: VMD macros that must mean the same atoms here as they do in a TCL script.
MACROS = {
    "water": "resname " + " ".join(sorted(T.WATER)),
    "waters": "resname " + " ".join(sorted(T.WATER)),
    "ions": "resname " + " ".join(sorted(T.ION)),
    "ion": "resname " + " ".join(sorted(T.ION)),
    "glycan": "resname " + " ".join(sorted(T.GLYCAN)),
    "glycans": "resname " + " ".join(sorted(T.GLYCAN)),
    "lipid": "resname " + " ".join(sorted(T.LIPID)),
    "lipids": "resname " + " ".join(sorted(T.LIPID)),
    "noh": HEAVY,
    "hydrogen": "name H*",
}

_WITHIN = re.compile(r"\b(ex)?within\s+([\d.]+)\s+of\s+", re.I)
#: One alternation over every macro name, longest first, matched case-sensitively
#: so an uppercase residue name (ION, MAN, CL) can never be mistaken for a macro.
_MACRO_RE = re.compile(
    r"\b(" + "|".join(sorted(MACROS, key=len, reverse=True)) + r")\b")
#: Keywords whose ARGUMENT is a list of values, which may contain ranges.
_KEYWORDS = ("segid", "segname", "resid", "resnum", "residue", "chain",
             "chainID", "resname", "name", "type", "index")
_OPERATORS = {"and", "or", "not", "(", ")", "byres", "around", "same",
              "residue", "as", "within", "of", "exwithin"}


def _expand_range(lo, hi):
    """`GA1 to GA23` -> the 23 segids; `384 to 401` -> `384:401`."""
    if lo.isdigit() and hi.isdigit():
        return f"{lo}:{hi}"
    m_lo = re.match(r"^([A-Za-z]*)(\d+)$", lo)
    m_hi = re.match(r"^([A-Za-z]*)(\d+)$", hi)
    if m_lo and m_hi and m_lo.group(1).upper() == m_hi.group(1).upper():
        stem = m_lo.group(1)
        a, b = int(m_lo.group(2)), int(m_hi.group(2))
        # `G001 to G020` keeps its padding; `GA1 to GA23` has none to keep.
        width = len(m_lo.group(2)) if m_lo.group(2).startswith("0") else 0
        if a <= b:
            return " ".join(f"{stem}{i:0{width}d}" for i in range(a, b + 1))
    if len(lo) == 1 and len(hi) == 1 and lo.isalpha() and hi.isalpha():
        return " ".join(chr(c) for c in range(ord(lo), ord(hi) + 1))
    raise ValueError(f"cannot expand the range '{lo} to {hi}'")


def _expand_ranges(text):
    """Expand EVERY `A to B` inside a keyword's value list, not just the first.

    `segid GA1 to GA23 GB1 to GB23 GC1 to GC23 GAE GAH` is one keyword with
    three ranges. A single regex substitution expanded the first and left the
    rest as literal segids called `to`, `GB23`, `GC1` -- which MDAnalysis
    accepted, matching 33 of the 75 glycan segments. Not empty, so no error;
    just a quietly wrong answer, which is the failure mode this whole module
    exists to prevent.
    """
    toks = text.replace("(", " ( ").replace(")", " ) ").split()
    out, i = [], 0
    while i < len(toks):
        tok = toks[i]
        low = tok.lower()
        if low in [k.lower() for k in _KEYWORDS]:
            key = ("resid" if low in ("resid", "resnum", "residue")
                   else "chainID" if low in ("chain", "chainid")
                   else "segid" if low in ("segid", "segname") else tok)
            out.append(key)
            i += 1
            values = []
            while i < len(toks) and toks[i].lower() not in _OPERATORS:
                values.append(toks[i])
                i += 1
            j = 0
            while j < len(values):
                if (j + 2 < len(values) and values[j + 1].lower() == "to"):
                    out.append(_expand_range(values[j], values[j + 2]))
                    j += 3
                else:
                    out.append(values[j])
                    j += 1
            continue
        out.append(tok)
        i += 1
    s = " ".join(out).replace("( ", "(").replace(" )", ")")
    if re.search(r"\bto\b", s):
        raise ValueError(
            f"a range was not expanded in {text!r} -> {s!r}: `to` survives "
            f"translation, and MDAnalysis would read it as a literal value")
    return s


def translate(text):
    """A VMD/CHARMM selection string as MDAnalysis understands it.

    Pure text in, pure text out: testable without a topology, and readable in a
    diff when someone asks what a stage actually selected.
    """
    s = " ".join(str(text or "").split())
    if not s:
        return ""

    # `within X of` before anything else: its argument may contain macros.
    def _w(m):
        return ("around " + m.group(2) + " ") if m.group(1) else \
            ("AROUNDPLUS " + m.group(2) + " ")
    s = _WITHIN.sub(_w, s)

    s = _expand_ranges(s)

    # Macros in ONE pass, case-sensitively, longest name first: a token loop
    # missed `noh)` for want of the bracket, and a second pass would happily
    # expand the `CL` it had just written for `ions`.
    s = _MACRO_RE.sub(lambda m: "(" + MACROS[m.group(1)] + ")", s)
    s = re.sub(r"\bsegname\b", "segid", s)
    s = re.sub(r"\bchain(ID)?\b", "chainID", s, flags=re.I)

    # VMD's inclusive `within`: the reference group is part of the answer.
    while "AROUNDPLUS" in s:
        i = s.index("AROUNDPLUS")
        head, tail = s[:i], s[i + len("AROUNDPLUS"):].strip()
        dist, rest = tail.split(" ", 1)
        ref, after = _first_group(rest)
        s = f"{head}((around {dist} ({ref})) or ({ref})) {after}"
    return " ".join(s.split())


def _first_group(text):
    """Split `(a and b) and c` into the reference group and what follows.

    Taking the whole tail as the reference -- which is what a naive split does
    -- silently swallowed any condition written after a `within`, turning
    `water within 3.5 of (X) and name OH2` into a selection of everything near
    X that is also named OH2 *inside the distance test*. Same atoms by luck
    here, wrong in general.
    """
    text = text.strip()
    if not text.startswith("("):
        parts = text.split(" ", 1)
        return parts[0], (parts[1] if len(parts) > 1 else "")
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return text[1:i], text[i + 1:].strip()
    raise ValueError(f"unbalanced parentheses in {text!r}")


class Resolved:
    """A selection that has been checked against a real topology."""

    __slots__ = ("name", "source", "mdanalysis", "n_atoms", "n_heavy",
                 "segids", "n_segments", "note")

    def __init__(self, name, source, mdanalysis, n_atoms, n_heavy, segids,
                 note, n_segments=None):
        self.name = name
        self.source = source
        self.mdanalysis = mdanalysis
        self.n_atoms = n_atoms
        self.n_heavy = n_heavy
        self.segids = segids
        self.n_segments = n_segments if n_segments is not None else len(segids)
        self.note = note

    def as_dict(self):
        return {"name": self.name, "selection": self.source,
                "mdanalysis": self.mdanalysis, "atoms": self.n_atoms,
                "heavy_atoms": self.n_heavy, "segments": self.n_segments,
                "segids": self.segids, "note": self.note}

    def __repr__(self):                                    # pragma: no cover
        return f"<{self.name}: {self.n_atoms} atoms ({self.n_heavy} heavy)>"


class EmptySelection(ValueError):
    """A selection matched nothing. Always a mistake, never a result."""


#: Names that can only be a water oxygen: CHARMM's OH2, and OW as GROMACS,
#: TIP4P and SPC files write it. The rdf, grid and ion stages accepted OH2
#: alone, so a GROMACS system passed the water stage and failed those three.
#: A bare `O` is NOT here: it is also every backbone carbonyl, and only the
#: water stage, which checks one oxygen per water residue, may fall back to it.
WATER_O_NAMES = ("OH2", "OW")


def water_oxygen_name(water_atoms):
    """The name of the water oxygen in `water_atoms`, or None.

    A name counts only when it picks exactly one atom per water residue, so a
    name that matches some waters cannot quietly give a smaller shell.
    """
    n_res = water_atoms.n_residues
    for nm in WATER_O_NAMES:
        if len(water_atoms.select_atoms(f"name {nm}")) == n_res:
            return nm
    return None


def resolve(universe, name, text, allow_empty=False):
    """Check one selection against a Universe and report what it holds."""
    sel = translate(text)
    if not sel:
        raise EmptySelection(f"{name}: no selection given")
    try:
        ag = universe.select_atoms(sel)
    except Exception as exc:                               # noqa: BLE001
        raise ValueError(f"{name}: {text!r} -> {sel!r}: {exc}") from exc
    if not len(ag) and not allow_empty:
        raise EmptySelection(
            f"{name}: {text!r} matched no atoms in this topology "
            f"(translated to {sel!r}). Check segid vs chainID, and that the "
            f"names match this system.")
    heavy = len(ag.select_atoms(HEAVY)) if len(ag) else 0
    segids = sorted({str(s) for s in getattr(ag, "segids", [])})
    note = None
    if "chainID" in sel:
        note = ("uses chainID: in these systems that collapses the antibody "
                "copies and all glycans -- prefer segid")
    # The list is capped for display; the COUNT is not, or a truncated list
    # reads as a truncated selection (it did, once, for 75 glycan segments).
    return Resolved(name, str(text), sel, len(ag), heavy, segids[:40], note,
                    n_segments=len(segids))


def resolve_all(universe, selections, allow_empty=()):
    """Resolve a whole spec's selections, reporting every failure at once.

    One failure at a time would mean one long trajectory-free round trip per
    mistake; a person fixing selections wants the whole list.
    """
    out, problems = {}, []
    for name, text in selections.items():
        if not str(text or "").strip():
            continue
        try:
            out[name] = resolve(universe, name, text,
                                allow_empty=name in allow_empty)
        except (EmptySelection, ValueError) as exc:
            problems.append(str(exc))
    return out, problems
