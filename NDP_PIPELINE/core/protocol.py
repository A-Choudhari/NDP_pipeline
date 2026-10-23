#!/usr/bin/env python3
"""
protocol.py -- how a run got to production: every stage, from its own files.

A deposited MD is its production trajectory, but the record is not complete
without the path to it: minimisation, heating, equilibration, and any restart
of production itself. ABZYME/3-OCCLUDED-OPEN/MDs/2-REPLICA/1R, for example:

    min1 -> min2 -> heating -> npt1 -> npt2 -> pr -> pr_r1

Each stage is read from three files, each for what only it can say:

    config      what was asked for: minimise or dynamics, timestep, temperature
                (or a heating ramp), ensemble, pressure, positional restraints,
                random seed. Disabled `if {0}` blocks are ignored and `set`
                variables resolved, exactly as mdconfig.py reads production.
    log         what actually ran: steps from FIRST TIMESTEP to the last ENERGY
                line, and whether it reached "End of program" or was cut off.
    trajectory  what is on disk: frames, atoms, frame spacing, size.

The order is the configs' own chain: each stage names the restart files it
starts from (binCoordinates npt2.coor -> after npt2). A config whose output
never appeared (PR2.conf, prepared and not run) is not a stage.

NAMD only, for now. Other engines return no stages rather than a guess.
Reads files; writes nothing.
"""

import os
import re

from . import mdconfig

_RAMP = re.compile(
    r"for\s*\{\s*set\s+(\w+)\s+([-\d.]+)\s*\}\s*\{\s*\$\1\s*<=?\s*([-\d.]+)\s*\}"
    r"\s*\{\s*incr\s+\1\s+([-\d.]+)\s*\}", re.I)
_ENERGY = re.compile(r"^ENERGY:\s+(\d+)", re.M)
_FIRST = re.compile(r"^Info: FIRST TIMESTEP\s+(\d+)", re.M)


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return int(f) if f == int(f) else f


def _is_namd_config(path):
    try:
        with open(path, errors="replace") as fh:
            head = fh.read(200_000)
    except OSError:
        return False
    return bool(re.search(r"^\s*outputName\s", head, re.I | re.M)
                and re.search(r"^\s*(structure|coordinates)\s", head, re.I | re.M))


def _read_log(path):
    """(first_step, last_step, completed) from a NAMD log, or Nones.

    The head holds FIRST TIMESTEP, the tail the last ENERGY line and "End of
    program". Only those two ends are read: production logs here are 130 MB.
    """
    if not os.path.isfile(path):
        return None, None, None
    try:
        with open(path, "rb") as fh:
            head = fh.read(300_000).decode("utf-8", "replace")
            size = fh.seek(0, os.SEEK_END)
            fh.seek(max(0, size - 200_000))
            tail = fh.read().decode("utf-8", "replace")
    except OSError:
        return None, None, None
    first = _FIRST.search(head)
    last = _ENERGY.findall(tail) or _ENERGY.findall(head)
    return (int(first.group(1)) if first else None,
            int(last[-1]) if last else None,
            "End of program" in tail)


def _trajectory(path):
    if not os.path.isfile(path):
        return None
    from core.trajscan import frame_info
    nf, natom, dt_ps, _ = frame_info(path)
    return {"file": os.path.basename(path), "frames": nf, "atoms": natom,
            "frame_interval_ps": _num(round(dt_ps, 4)) if dt_ps else None,
            "size_gb": round(os.path.getsize(path) / 1e9, 2)}


def read_stage(config_path, logs=True):
    """One stage from one NAMD config and the files it produced, or None.

    `logs=False` skips the log (steps actually run, completion) and the
    trajectory header, for callers that only need what the config asks for.
    """
    try:
        with open(config_path, errors="replace") as fh:
            raw = fh.read()
    except OSError:
        return None
    live = mdconfig.strip_disabled_blocks(raw)
    text, _ = mdconfig.resolve_tcl(live)
    d = lambda name: mdconfig._directive(text, name)        # noqa: E731
    run_dir = os.path.dirname(config_path)
    out = d("outputName")
    if not out:
        return None
    out = os.path.basename(out)
    dcd = os.path.join(run_dir, f"{out}.dcd")
    log = os.path.join(run_dir, f"{out}.log")
    if not (os.path.isfile(dcd) or os.path.isfile(log)):
        return None                                      # prepared, never run

    started_from = d("binCoordinates")
    if started_from:
        started_from = re.sub(r"(\.restart)?\.coor$", "", os.path.basename(started_from))

    minimise = d("minimize")
    ramp = _RAMP.search(live)
    piston = (d("langevinPiston") or "").lower() in ("on", "yes", "true", "1")
    restrained = (d("constraints") or "").lower() in ("on", "yes", "true", "1")
    timestep = _num(d("timestep"))

    if minimise:
        kind = "minimisation"
    elif ramp:
        kind = "heating"
    else:
        from core.trajscan import stage_of
        kind = {"equilibration": "equilibration", "production": "production",
                "heating": "heating", "minimisation": "minimisation"
                }.get(stage_of(f"{out}.dcd"), "dynamics")

    # What the config ASKS for: the step it starts from and the steps it runs.
    # `run` values are summed only outside a Tcl loop (a heating ramp runs its
    # `run` once per temperature, which a sum of the text would miss).
    first_cfg = int(_num(d("firsttimestep")) or 0)
    planned = None
    if minimise:
        planned = int(_num(minimise) or 0) or None
    elif not ramp:
        runs = [_num(v) for v in mdconfig._all_directives(text, "run")]
        if runs and all(r is not None for r in runs):
            planned = int(sum(runs))
        else:
            # `run [expr {$num_steps_pr - $ts}]`: the lab's production configs
            # compute how much is left from the restart's own step, so the
            # number in the file is the FINAL step, not a count. Taking it as a
            # count would report a 504 ns run as twice its length.
            target = _final_step(live, text)
            planned = None if target is None else max(0, target - first_cfg)

    first, last, completed = _read_log(log) if logs else (None, None, None)
    if first is None and last is not None:
        # NAMD prints FIRST TIMESTEP only when it is not zero; the config's own
        # `firsttimestep` (default 0) is then the start.
        first = first_cfg
    steps = (last - first) if (first is not None and last is not None) else None
    stage = {
        "stage": kind,
        "name": out,
        "config": os.path.basename(config_path),
        "log": os.path.basename(log) if os.path.isfile(log) else None,
        "started_from": started_from,
        "first_step": first_cfg,
        "planned_steps": planned,
        "steps": steps,
        "status": (None if completed is None else
                   "completed" if completed else f"interrupted at step {last}"),
        "timestep_fs": None if minimise else timestep,
        "duration_ns": (round(steps * timestep / 1e6, 3)
                        if steps is not None and timestep and not minimise else None),
        # A minimiser has no temperature; the value in its config only seeds
        # the velocities handed to the next stage.
        "temperature_K": (None if ramp or minimise
                          else _num(d("langevinTemp") or d("temperature"))),
        "temperature_ramp_K": ({"from": _num(ramp.group(2)), "to": _num(ramp.group(3)),
                                "step": _num(ramp.group(4))} if ramp else None),
        "ensemble": None if minimise else ("NPT" if piston else "NVT"),
        "pressure_bar": _num(d("langevinPistonTarget")) if piston and not minimise else None,
        "restraints": ({"reference": os.path.basename(d("consref") or "") or None,
                        "force_constants": os.path.basename(d("conskfile") or "") or None,
                        "column": d("conskcol")} if restrained else None),
        "seed": _num(d("seed")),
        "trajectory": _trajectory(dcd) if logs else None,
    }
    return stage


_RUN_EXPR = re.compile(r"^\s*run\s+\[\s*expr\s*\{?\s*\$(\w+)", re.I | re.M)


def _final_step(live_text, resolved_text):
    """The absolute step a `run [expr {$num_steps_pr - $ts}]` line aims at."""
    m = _RUN_EXPR.search(live_text)
    if not m:
        return None
    var = m.group(1)
    hit = re.search(rf"^\s*set\s+{var}\s+([\d.]+)", live_text, re.I | re.M)
    return int(float(hit.group(1))) if hit else None


def stages_of(run_dir, logs=True):
    """Every stage that ran in `run_dir`, in the order its configs chain.

    [] when no NAMD config there produced anything -- including for other
    engines, which this does not read yet.
    """
    try:
        names = sorted(os.listdir(run_dir))
    except OSError:
        return []
    stages = []
    for n in names:
        p = os.path.join(run_dir, n)
        if n.endswith((".conf", ".namd", ".inp")) and os.path.isfile(p) \
                and _is_namd_config(p):
            st = read_stage(p, logs=logs)
            if st:
                stages.append(st)
    # Two configs writing one output (a stage re-run under another name) keep
    # the newer config only: the files on disk are its.
    by_out = {}
    for st in stages:
        prev = by_out.get(st["name"])
        cur_m = os.path.getmtime(os.path.join(run_dir, st["config"]))
        if prev is None or cur_m >= os.path.getmtime(os.path.join(run_dir, prev["config"])):
            by_out[st["name"]] = st
    stages = list(by_out.values())

    def mtime(st):
        for f in (st.get("log"), (st.get("trajectory") or {}).get("file")):
            if f and os.path.exists(os.path.join(run_dir, f)):
                return os.path.getmtime(os.path.join(run_dir, f))
        return 0.0

    # Follow the chain: roots (nothing in this directory feeds them) first,
    # then each stage's children; ties by when their output was written.
    known = {st["name"] for st in stages}
    children = {}
    for st in stages:
        parent = st["started_from"] if st["started_from"] in known else None
        children.setdefault(parent, []).append(st)
    ordered, seen = [], set()

    def walk(parent):
        for st in sorted(children.get(parent, []), key=mtime):
            if st["name"] in seen:
                continue
            seen.add(st["name"])
            ordered.append(st)
            walk(st["name"])
    walk(None)
    ordered += [st for st in sorted(stages, key=mtime) if st["name"] not in seen]
    return ordered


def planned_production_ns(run_dir):
    """How long production was MEANT to run, in ns, from the configs; or None.

    The lab writes the target into the production config: PR.conf runs
    504,000,000 steps of 2 fs from step 0, and its restart pr_r1.conf starts at
    264,499,998 and runs to the same 504,000,000 -- one target, 1008 ns. So the
    target is from the first production stage's start to the furthest end any
    production config asks for. Configs only; no log or trajectory is opened.
    """
    prod = [s for s in stages_of(run_dir, logs=False)
            if s["stage"] == "production" and s.get("planned_steps")
            and s.get("timestep_fs")]
    if not prod:
        return None
    if len({s["timestep_fs"] for s in prod}) > 1:
        return None                          # mixed timesteps: no single target
    start = min(s["first_step"] for s in prod)
    end = max(s["first_step"] + s["planned_steps"] for s in prod)
    return round((end - start) * prod[0]["timestep_fs"] / 1e6, 3) or None


def summary(stages):
    """One sentence for a catalogue description, from a merged protocol."""
    def one(v):
        return None if isinstance(v, dict) and "per_replica" in v else v

    parts = []
    for st in stages:
        kind, bits = st.get("stage"), []
        if kind == "minimisation":
            steps = one(st.get("steps"))
            bits.append(f"{steps} steps" if steps else "")
        else:
            ramp = one(st.get("temperature_ramp_K"))
            if ramp:
                bits.append(f"{ramp['from']} to {ramp['to']} K")
            elif one(st.get("temperature_K")) is not None:
                bits.append(f"{st['temperature_K']} K")
            if one(st.get("ensemble")):
                bits.append(st["ensemble"])
            dur = st.get("duration_ns")
            if isinstance(dur, dict):
                vals = [v for v in dur["per_replica"].values() if v is not None]
                if vals:
                    bits.append(f"{min(vals):g}-{max(vals):g} ns")
            elif dur is not None:
                bits.append(f"{dur:g} ns")
        if st.get("restraints"):
            bits.append("restrained")
        if "interrupted" in repr(st.get("status")):
            bits.append("interrupted, continued by the next stage")
        detail = ", ".join(b for b in bits if b)
        parts.append(f"{st.get('name')} ({kind}{': ' + detail if detail else ''})")
    return "; ".join(parts)


def merge(per_replica):
    """One protocol for replicas of one system.

    `per_replica` is {replica label: [stages]}. Stages are matched by name, in
    the first replica's order. A field every replica agrees on is written once;
    one they differ on (the seed, how far an interrupted run got) becomes
    {"per_replica": {replica: value}}, so nothing is averaged or picked and a
    reader can tell it from a value that is itself a mapping (a heating ramp).
    """
    labels = list(per_replica)
    if not labels:
        return []
    order = []
    for lab in labels:
        for st in per_replica[lab]:
            if st["name"] not in order:
                order.append(st["name"])
    merged = []
    for name in order:
        have = {lab: next((s for s in per_replica[lab] if s["name"] == name), None)
                for lab in labels}
        present = {lab: s for lab, s in have.items() if s is not None}
        entry = {}
        keys = []
        for s in present.values():
            keys += [k for k in s if k not in keys]
        for k in keys:
            vals = {lab: s.get(k) for lab, s in present.items()}
            distinct = {repr(v) for v in vals.values()}
            entry[k] = (next(iter(vals.values())) if len(distinct) == 1
                        else {"per_replica": vals})
        if len(present) < len(labels):
            entry["replicas"] = sorted(present)
        merged.append(entry)
    return merged
