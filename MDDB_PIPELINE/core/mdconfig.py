#!/usr/bin/env python3
"""
mdconfig.py -- recover the six simulation-metadata fields from MD config files.

The NDP pilot document lists these as "the fields that currently live only in
people's heads":

    timestep (fs) . temp (K) . ensemble (NVT/NPT) . ff . wat . boxtype

None of them is in a trajectory header, so the substrate cannot derive them.
They ARE on disk, in the config file the job was launched with.

THE TRAP THAT MAKES THIS NON-TRIVIAL. NAMD configs are Tcl scripts, and the
interesting values are almost always variables:

    set temp        310
    ...
    langevinTemp    $temp

A scanner that greps for `langevinTemp` records the string "$temp" and writes
it into a metadata record that looks perfectly well-formed. So `set` is
resolved first, before any field is read. Verified against
ABZYME/4-OPEN/CDR/3-YZ23/0-MD/REPLICA/1R/PR.conf, where the six fields appear
as: `set temp 310`, `timestep 2.0`, `langevin on` + `langevinPiston on`,
`set dir_ff ./CHARMM36m` with `parameters .../toppar_all36_*.str`,
`parameters .../toppar_water_ions_namd.update.str`, and cell vectors in a
sibling `.xsc`.

WHAT IS NOT GUESSED. `wat` is left unset: the water model is implied by a
parameter FILENAME (`toppar_water_ions...`), and mapping that to "TIP3P" is an
inference about the science, not a reading of the file. The filename is
recorded as evidence so a person can decide. Same rule the substrate already
applies to `system: unknown` -- absence is recorded, never filled in.

Every value carries the file and line it came from, so any number in a
published record can be traced back.
"""

import os
import re

# NAMD writes fs; AMBER writes ps. Both are normalised to fs on the way out.
NAMD_CONF = ("*.conf", "*.namd", "*.inp")
AMBER_IN = ("*.mdin", "mdin", "*.in")
GROMACS_MDP = ("*.mdp",)

# Production first: a run's metadata should describe its production stage, not
# the minimisation that happens to sort first alphabetically.
PROD_HINTS = ("pr", "prod", "production", "md")
SKIP_HINTS = ("min", "heat", "npt", "nvt", "eq", "equil", "backup", "bak")

SET_RE = re.compile(r"^\s*set\s+([A-Za-z_]\w*)\s+(\S+)", re.M)
VAR_RE = re.compile(r"\$\{?([A-Za-z_]\w*)\}?")


def resolve_tcl(text):
    """Substitute `set` variables. Returns (resolved_text, {var: value}).

    Two passes, because a `set` may itself reference an earlier one
    (`set outname pr.part$seg`). Two is enough for every config seen here and
    it cannot loop; a variable that stays unresolved is left as-is so the
    caller can see that it failed rather than receiving a plausible wrong value.
    """
    variables = {m.group(1): m.group(2).strip('"') for m in SET_RE.finditer(text)}
    for _ in range(2):
        for key, val in list(variables.items()):
            variables[key] = VAR_RE.sub(
                lambda m: variables.get(m.group(1), m.group(0)), val)
    resolved = VAR_RE.sub(
        lambda m: variables.get(m.group(1), m.group(0)), text)
    return resolved, variables


IF_DISABLED_RE = re.compile(r"^\s*if\s*\{\s*0\s*\}", re.I)


def strip_disabled_blocks(text):
    """Blank out lines inside a Tcl `if {0} { ... }` block.

    NAMD templates ship every option and switch behaviour with `if {0}` /
    `if {1}`, so the file CONTAINS `langevinPiston on` whether or not pressure
    coupling is active. A scanner that greps for the directive reports NPT for
    the whole corpus. Measured across this workspace: naive grep gives 1,366
    NPT, brace-tracking gives 758 NPT and 608 NVT -- about 42% would carry the
    wrong ensemble into a deposited record, and nothing downstream could tell.

    Lines are blanked rather than deleted so that any line number quoted as
    evidence still refers to the real line in the file.
    """
    lines = text.split("\n")
    out, i = list(lines), 0
    while i < len(lines):
        if IF_DISABLED_RE.match(lines[i]):
            depth = 0
            started = False
            j = i
            while j < len(lines):
                depth += lines[j].count("{") - lines[j].count("}")
                if lines[j].count("{"):
                    started = True
                out[j] = ""
                j += 1
                if started and depth <= 0:
                    break
            i = j
            continue
        i += 1
    return "\n".join(out)


def _tleap_sidecar(run_dir):
    """AMBER records its force field and water model in the leap script, not in
    mdin. Those files sit in the run directory beside the trajectory (207 of
    207 here), so `wat` is genuinely DECLARED for AMBER -- unlike NAMD, where it
    can only ever be inferred from a parameter filename.
    """
    try:
        names = [n for n in sorted(os.listdir(run_dir))
                 if n.lower().startswith("tleap") and n.lower().endswith(".in")]
    except OSError:
        return None, None, None
    for name in names:
        path = os.path.join(run_dir, name)
        try:
            with open(path, errors="replace") as fh:
                txt = fh.read()
        except OSError:
            continue
        ff, wat = [], None
        for m in re.finditer(r"^\s*source\s+leaprc\.(\S+)", txt, re.I | re.M):
            token = m.group(1)
            if token.lower().startswith("water."):
                wat = token.split(".", 1)[1].upper().replace("TIP3P", "TIP3P")
            else:
                ff.append(token.split(".")[-1])
        if ff or wat:
            return ff or None, wat, name
    return None, None, None


def _directive(text, name):
    """Last value of a NAMD directive. Last, not first: configs commonly set a
    value, then override it inside a conditional later in the file."""
    hits = re.findall(rf"^\s*{name}\s+([^\s;#]+)", text, re.I | re.M)
    return hits[-1] if hits else None


def _all_directives(text, name):
    return re.findall(rf"^\s*{name}\s+([^\s;#]+)", text, re.I | re.M)


def _amber_nml(text, key):
    """AMBER mdin is Fortran namelist: `key = value`, comma or newline
    separated, and case-insensitive."""
    m = re.search(rf"\b{key}\s*=\s*([-\d.eE+]+)", text, re.I)
    return m.group(1) if m else None


def pick_config(run_dir):
    """The config that describes the production stage of this run, or None."""
    try:
        names = sorted(os.listdir(run_dir))
    except OSError:
        return None, None

    # The stems of every trajectory sitting in this directory. A config that
    # produced output beats one that did not, whatever it is called.
    produced = {n.lower().rsplit(".", 1)[0] for n in names
                if n.lower().endswith((".dcd", ".nc", ".netcdf", ".xtc",
                                       ".trr", ".mdcrd", ".crd.nc"))}

    def graded(cands, engine):
        scored = []
        for n in cands:
            low = n.lower()
            if any(low.endswith(s) or s in low for s in ("bak", "backup", ".old")):
                continue
            score = 0
            base = low.rsplit(".", 1)[0]
            if any(base.startswith(h) for h in PROD_HINTS):
                score += 2
            if any(h in base for h in SKIP_HINTS):
                score -= 2
            # OUTPUT IS THE EVIDENCE. A replica can hold two generations of
            # production config, and the shorter name is not the later one:
            # X-MODELS/MD_5REP holds `pr01..pr10.conf` at 2 fs beside
            # `prh01..prh20.conf` at 4 fs with hydrogen mass repartitioning,
            # and only the `prh` set ever ran -- there is no `pr01.dcd`. Name
            # length picked the abandoned one, so three records published a
            # 2 fs timestep for a 4 fs run while their lengths, read from the
            # trajectory headers, were right. A config whose own output is in
            # the directory is the one that describes this run.
            if base in produced:
                score += 3
            scored.append((score, len(n), n))
        if not scored:
            return None
        scored.sort(key=lambda t: (-t[0], t[1]))
        return scored[0][2] if scored[0][0] >= 0 else None

    import fnmatch
    for pats, engine in ((NAMD_CONF, "NAMD"), (AMBER_IN, "AMBER"),
                         (GROMACS_MDP, "GROMACS")):
        cands = [n for n in names
                 if any(fnmatch.fnmatch(n.lower(), p) for p in pats)]
        # AN EXTENSION IS NOT A PROMISE. `.in` is used by cpptraj
        # (`cpptraj-mergeddcd.in`), by MMPBSA (`mmpbsa.in`) and by tleap just as
        # much as by sander. Picking one of those yielded a "config" with no
        # timestep and no temperature for 123 runs -- a silent hole in the
        # metadata rather than an error. Confirm from the CONTENT that the file
        # is what the extension claims.
        cands = [n for n in cands
                 if _is_md_config(os.path.join(run_dir, n), engine)]
        pick = graded(cands, engine)
        if pick:
            return os.path.join(run_dir, pick), engine
    return None, None


# Signatures of files that wear an MD-config extension but are not one.
NOT_MD_CONFIG = ("trajin", "trajout", "parm ", "cpptraj", "&general",
                 "&gb", "&pb", "startframe", "loadpdb", "solvatebox",
                 "source leaprc", "saveamberparm")


def _is_md_config(path, engine):
    """True when the file actually configures a dynamics run."""
    try:
        with open(path, errors="replace") as fh:
            head = fh.read(4096).lower()
    except OSError:
        return False
    if any(sig in head for sig in NOT_MD_CONFIG):
        return False
    if engine == "AMBER":
        # sander/pmemd control namelist, with an actual integrator setting.
        return "&cntrl" in head and ("nstlim" in head or "dt" in head)
    if engine == "NAMD":
        return "timestep" in head or "numsteps" in head or "run " in head
    if engine == "GROMACS":
        return "integrator" in head
    return False


def _boxtype_from_xsc(run_dir):
    """NAMD .xsc carries the cell vectors: a_x a_y a_z b_x .. c_z.

    Zero off-diagonals mean an orthogonal box; anything else is triclinic. The
    values are read, not assumed -- an octahedral box is a real possibility and
    reporting "Cubic" for one would be a silent error in a published record.
    """
    try:
        xscs = [n for n in sorted(os.listdir(run_dir)) if n.endswith(".xsc")]
    except OSError:
        return None, None
    for name in xscs:
        path = os.path.join(run_dir, name)
        try:
            with open(path, errors="replace") as fh:
                rows = [ln for ln in fh if ln.strip()
                        and not ln.lstrip().startswith("#")]
        except OSError:
            continue
        if not rows:
            continue
        parts = rows[-1].split()
        if len(parts) < 10:
            continue
        try:
            v = [float(x) for x in parts[1:10]]
        except ValueError:
            continue
        off = [v[1], v[2], v[3], v[5], v[6], v[7]]
        kind = "Cubic" if all(abs(x) < 1e-6 for x in off) else "Triclinic"
        return kind, f"{name}: off-diagonals {'all zero' if kind=='Cubic' else 'non-zero'}"
    return None, None


def pick_config_inherited(run_dir, root, max_up=4):
    """Fall back to a config in an ancestor directory, when the run dir has none.

    880 of the 1,447 NAMD configs in this workspace live in template folders
    (`namd_inputfiles/`, `SCRIPT_AUDIT/MD_SETUP/`) that are copied into a run
    directory at launch; the copy is often not kept. Without a fallback, more
    than half of all runs report no simulation metadata at all.

    An inherited config is WEAKER EVIDENCE than one sitting next to the
    trajectory: it is the template the run was probably launched from, not
    provably the one it used. So it is labelled `inherited`, and the caller
    cross-checks it against the trajectory header before trusting it -- a
    template whose framestep disagrees with what was actually written plainly
    does not describe this run, and its values are dropped rather than
    published.
    """
    d = os.path.dirname(os.path.abspath(run_dir))
    root = os.path.abspath(root)
    for _ in range(max_up):
        # THE ROOT ITSELF IS A CANDIDATE. `d == root` used to break BEFORE the
        # directory was scanned, so a config at the top of the folder was never
        # seen -- and that is exactly where these runs keep it: the aptamer
        # complexes hold `PR.conf`, `NPT1/2.conf`, `heating.conf` and
        # `min1/2.conf` beside the replica directories `pr/`, `rep/`, `rep2/`,
        # one copy for the three. Every one of those replicas reported
        # `no_config` and published no force field, ensemble, temperature,
        # pressure or timestep at all.
        if not (d == root or d.startswith(root + os.sep)):
            break
        try:
            subs = [d] + [os.path.join(d, n) for n in sorted(os.listdir(d))
                          if os.path.isdir(os.path.join(d, n))
                          and ("input" in n.lower() or "script" in n.lower()
                               or "setup" in n.lower() or "conf" in n.lower())]
        except OSError:
            subs = [d]
        for cand_dir in subs:
            path, engine = pick_config(cand_dir)
            if path:
                return path, engine
        if d == root:
            break          # scanned it; never walk above the folder
        d = os.path.dirname(d)
    return None, None


def scan_run(run_dir, root=None, allow_inherited=True):
    """Return {field: value} plus `_evidence` mapping each field to its source.

    A field with no evidence is ABSENT from the dict -- never present-but-null,
    because a null in a metadata record reads as "measured and found to be
    nothing" rather than "not established".
    """
    out, ev = {}, {}
    path, engine = pick_config(run_dir)
    scope = "run"
    if not path and allow_inherited and root:
        path, engine = pick_config_inherited(run_dir, root)
        scope = "inherited"
    if not path:
        box, boxev = _boxtype_from_xsc(run_dir)
        if box:
            out["boxtype"], ev["boxtype"] = box, boxev
        return out, ev, None, None, None

    try:
        with open(path, errors="replace") as fh:
            raw = fh.read()
    except OSError:
        return out, ev, None, None, None

    base = os.path.basename(path)
    if engine == "NAMD":
        # Strip disabled blocks BEFORE resolving variables, so a value that only
        # exists inside `if {0}` never reaches the output.
        text, variables = resolve_tcl(strip_disabled_blocks(raw))
    else:
        text, variables = raw, {}

    if engine == "NAMD":
        ts = _directive(text, "timestep")
        if ts:
            try:
                out["timestep"] = float(ts)
                ev["timestep"] = f"{base}: timestep {ts}"
            except ValueError:
                pass
        tmp = _directive(text, "langevinTemp") or _directive(text, "temperature")
        if tmp and not tmp.startswith("$"):
            try:
                out["temp"] = int(round(float(tmp)))
                ev["temp"] = (f"{base}: langevinTemp {tmp}"
                              + (" (via set temp)" if "temp" in variables else ""))
            except ValueError:
                pass
        piston = _directive(text, "langevinPiston")
        if piston is not None:
            npt = piston.lower() in ("on", "yes", "true", "1")
            out["ensemble"] = "NPT" if npt else "NVT"
            ev["ensemble"] = f"{base}: langevinPiston {piston}"
        params = [os.path.basename(p) for p in _all_directives(text, "parameters")]
        if params:
            ff, water = [], None
            for p in params:
                low = p.lower()
                if "water" in low or "ions" in low:
                    water = p
                    continue
                ff.append(p)
            fam = None
            dir_ff = variables.get("dir_ff", "")
            if "charmm36m" in dir_ff.lower():
                fam = "CHARMM36m"
            elif any("all36" in p.lower() for p in params):
                fam = "CHARMM36"
            if fam:
                out["ff"] = [fam]
                ev["ff"] = (f"{base}: dir_ff {dir_ff}" if dir_ff
                            else f"{base}: parameters {ff[:2]}")
            elif ff:
                ev["ff"] = f"{base}: parameters {ff[:3]} (family not named)"
            if water:
                # Deliberately NOT translated into a model name. See docstring.
                ev["wat"] = (f"{base}: parameters {water} "
                             f"-- water model implied, not stated; set `wat` "
                             f"in datasets.yaml if you know it")
        dcd = _directive(text, "dcdfreq")
        if dcd and "timestep" in out:
            try:
                # ns per frame -- the same quantity the substrate derives from
                # the DCD header, which makes it a free cross-check.
                out["framestep_from_config"] = (
                    int(dcd) * out["timestep"] / 1e6)
                ev["framestep_from_config"] = f"{base}: dcdfreq {dcd} x timestep"
            except ValueError:
                pass

    elif engine == "AMBER":
        dt = _amber_nml(text, "dt")
        if dt:
            try:
                out["timestep"] = float(dt) * 1000.0      # ps -> fs
                ev["timestep"] = f"{base}: dt {dt} ps"
            except ValueError:
                pass
        t0 = _amber_nml(text, "temp0")
        if t0:
            try:
                out["temp"] = int(round(float(t0)))
                ev["temp"] = f"{base}: temp0 {t0}"
            except ValueError:
                pass
        ntp = _amber_nml(text, "ntp")
        if ntp is not None:
            try:
                out["ensemble"] = "NPT" if int(float(ntp)) > 0 else "NVT"
                ev["ensemble"] = f"{base}: ntp {ntp}"
            except ValueError:
                pass
        ff, wat, leapfile = _tleap_sidecar(run_dir)
        if ff:
            out["ff"] = ff
            ev["ff"] = f"{leapfile}: source leaprc.{ff[0]} (+{len(ff)-1} more)" \
                if len(ff) > 1 else f"{leapfile}: source leaprc.{ff[0]}"
        if wat:
            # AMBER STATES the water model; NAMD only implies it. Recording it
            # here is reading the file, not inferring from a filename.
            out["wat"] = wat
            ev["wat"] = f"{leapfile}: source leaprc.water.{wat.lower()}"
        ntb = _amber_nml(text, "ntb")
        if ntb is not None:
            try:
                if int(float(ntb)) == 1:
                    out["boxtype"] = "Cubic"
                    ev["boxtype"] = f"{base}: ntb 1 (constant volume, orthogonal)"
            except ValueError:
                pass

    elif engine == "GROMACS":
        dt = re.search(r"^\s*dt\s*=\s*([\d.eE+-]+)", text, re.I | re.M)
        if dt:
            out["timestep"] = float(dt.group(1)) * 1000.0
            ev["timestep"] = f"{base}: dt {dt.group(1)} ps"
        ref = re.search(r"^\s*ref[-_]?t\s*=\s*([\d.]+)", text, re.I | re.M)
        if ref:
            out["temp"] = int(round(float(ref.group(1))))
            ev["temp"] = f"{base}: ref_t {ref.group(1)}"
        pc = re.search(r"^\s*pcoupl\s*=\s*(\S+)", text, re.I | re.M)
        if pc:
            out["ensemble"] = "NVT" if pc.group(1).lower() == "no" else "NPT"
            ev["ensemble"] = f"{base}: pcoupl {pc.group(1)}"

    if "boxtype" not in out:
        box, boxev = _boxtype_from_xsc(run_dir)
        if box:
            out["boxtype"], ev["boxtype"] = box, boxev
    if scope == "inherited":
        ev["_scope"] = f"config inherited from {os.path.basename(os.path.dirname(path))}/"
    return out, ev, path, engine, scope


FIELDS = ("timestep", "temp", "ensemble", "ff", "wat", "boxtype")


def scan(root, rel_paths, verbose=True, header_dt_ps=None):
    """One row per run directory. Returns a DataFrame.

    `header_dt_ps` maps rel_path -> the frame spacing the SUBSTRATE read from
    the trajectory header. When a config was INHERITED from a template rather
    than found beside the trajectory, it is checked against that number and
    discarded if it disagrees: a template that describes a different frame
    spacing is not the config this run used, and publishing its temperature
    would be a guess wearing the appearance of a measurement.
    """
    import pandas as pd
    header_dt_ps = header_dt_ps or {}
    rows, dropped = [], 0
    for i, rel in enumerate(rel_paths):
        run_dir = os.path.join(root, str(rel))
        vals, ev, cfg_path, engine, scope = scan_run(run_dir, root=root)
        if scope == "inherited" and vals.get("framestep_from_config") is not None:
            hdr = header_dt_ps.get(str(rel))
            if hdr is not None and abs(
                    vals["framestep_from_config"] * 1000.0 - float(hdr)) > 0.01:
                ev = {"_rejected": f"inherited config {os.path.basename(cfg_path)} "
                                   f"implies {vals['framestep_from_config']*1000:.0f} ps/frame "
                                   f"but the header says {float(hdr):.0f} ps -- "
                                   f"not this run's config"}
                vals, cfg_path, engine, scope = {}, None, None, None
                dropped += 1
        row = {"rel_path": str(rel),
               "config_path": (os.path.relpath(cfg_path, root)
                               if cfg_path else None),
               "config_scope": scope,
               "config_engine": engine}
        for f in FIELDS:
            row[f] = vals.get(f)
        row["framestep_from_config"] = vals.get("framestep_from_config")
        row["evidence"] = "; ".join(f"{k}={v}" for k, v in sorted(ev.items())) or None
        rows.append(row)
        if verbose and (i + 1) % 200 == 0:
            print(f"      {i+1} run dirs scanned", flush=True)
    if verbose and dropped:
        print(f"      {dropped} inherited config(s) rejected: their frame "
              f"spacing disagrees with the trajectory header")
    spec = {"rel_path": "object", "config_path": "object",
            "config_scope": "object", "config_engine": "object", "timestep": "float64", "temp": "float64",
            "ensemble": "object", "ff": "object", "wat": "object",
            "boxtype": "object", "framestep_from_config": "float64",
            "evidence": "object"}
    if not rows:
        return pd.DataFrame({c: pd.Series(dtype=t) for c, t in spec.items()})
    return pd.DataFrame(rows)
