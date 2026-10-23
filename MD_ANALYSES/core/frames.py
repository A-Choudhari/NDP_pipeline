#!/usr/bin/env python3
"""
frames.py -- which frame is this, and where does it live?

The point of the whole interactive-graph exercise: a reader clicks a peak on an
RMSD curve and has to be told *frame 2016, t = 403.2 ns, replica 1-R1, file
pr.part6.dcd, local frame 44*. Without that, a number on a plot cannot be
checked, reproduced or reused.

Everything here comes from the DCD headers themselves -- `istart` (the step of
the first frame), `nsavc` (steps between frames), `delta` (the timestep) and the
frame count the file SIZE implies. Nothing is inferred from file names, and a
file that holds fewer frames than its header claims (a copy still in flight) is
counted at what it really holds.

TWO INDEX SPACES, and the difference matters.

    raw index   position in the plain concatenation of the segments, which is
                what MDAnalysis gives you when it chains files and what the
                existing ANALYSIS/agg/ALL.tsv calls `raw_frame`.
    frame index position after restart seams are dropped -- the time axis with
                every simulated step appearing exactly once.

A restart re-emits a step the previous segment already wrote, with DIFFERENT
coordinates: the window resumed from an earlier checkpoint and re-simulated up
to it, so the pair is a seam, not a duplicate. The later-written copy is the one
dropped, because the earlier segment's frame is the one the trajectory actually
continued from.
"""

import json
import os

from .dcdheader import AKMA_PS, read_dcd_header


class Segment:
    """One trajectory file's place on the time axis."""

    __slots__ = ("path", "name", "first_step", "n_frames", "nsavc", "step_ps",
                 "natom", "truncated", "raw_start")

    def __init__(self, path, head, raw_start):
        self.path = path
        self.name = os.path.basename(path)
        self.first_step = head["istart"]
        self.n_frames = head["n_frames"]
        self.nsavc = head["nsavc"]
        self.step_ps = (head["delta_akma"] or 0.0) * AKMA_PS
        self.natom = head["natom"]
        self.truncated = bool(head.get("truncated"))
        self.raw_start = raw_start

    @property
    def last_step(self):
        return self.first_step + (self.n_frames - 1) * self.nsavc

    def step_of(self, local):
        return self.first_step + local * self.nsavc

    def as_dict(self):
        return {"file": self.name, "path": self.path, "frames": self.n_frames,
                "first_step": self.first_step, "last_step": self.last_step,
                "nsavc": self.nsavc, "truncated": self.truncated}


#: Below this, a derived frame spacing is a post-processing placeholder
#: rather than a measurement. Production here writes 100-200 ps frames; a
#: picosecond is already two orders of magnitude finer than anything real.
PLACEHOLDER_FLOOR_NS = 0.001


class Timeline:
    """The frames of one MD, in time order, each knowing where it came from."""

    analysed = ()
    analysed_path = None
    n_analysed = None
    framestep_basis = None

    @property
    def analysed_indexing(self):
        """"raw", "unique", or None when the file's length says neither.

        Unknown is reported as unknown: a provenance panel that guesses is
        worse than one that says it cannot map the frame, because the whole
        point of it is that you can open exactly what you clicked.
        """
        if not self.n_analysed:
            return None
        if self.n_analysed == len(self.frames) + len(self.seams):
            return "raw"
        if self.n_analysed == len(self.frames):
            return "unique"
        return None

    def __init__(self, md_id, segments, frames, seams, framestep_ns):
        self.md_id = md_id
        self.segments = segments
        self.frames = frames            # [(step, seg_i, local, raw_index)]
        self.seams = seams              # raw indices dropped as re-simulated
        self.framestep_ns = framestep_ns

    # ---------------------------------------------------------------- lookups
    def __len__(self):
        return len(self.frames)

    def time_ns(self, i):
        """When frame `i` happened, from whichever clock is trustworthy.

        Normally the frame's own step counter times it, which survives
        restarts and gaps: frame 500 of a segment that resumed at step
        1,000,000 is not at 100 ns just because it is the 500th read. But a
        post-processed copy carries neither a real step counter nor a real
        delta, and there the only clock is the interval a person declared --
        so the frame's POSITION times it instead. Mixing the two is what put
        a declared 100 ps run at 4.9e-05 ns a frame.
        """
        if self.framestep_basis == "declared" and self.framestep_ns:
            return round(i * self.framestep_ns, 4)
        return round(self.frames[i][0] * self.step_ps / 1000.0, 4)

    @property
    def step_ps(self):
        for s in self.segments:
            if s.step_ps:
                return s.step_ps
        return 0.0

    def where(self, i):
        """Everything known about frame `i`, ready to show a person."""
        step, seg_i, local, raw = self.frames[i]
        seg = self.segments[seg_i]
        return {"frame": i, "raw_frame": raw, "step": step,
                "time_ns": self.time_ns(i), "md_id": self.md_id,
                "file": seg.name, "path": seg.path,
                "local_frame": local, "of": seg.n_frames}

    def vmd_command(self, i):
        """What to paste into VMD to look at exactly this frame."""
        w = self.where(i)
        return (f"mol addfile {w['path']} type dcd "
                f"first {w['local_frame']} last {w['local_frame']} waitfor all")

    def duration_ns(self):
        if self.framestep_basis == "declared" and self.framestep_ns:
            return round(len(self) * self.framestep_ns, 3)
        return round(len(self) * self.step_ps * self.nsavc_common / 1000.0, 3) \
            if self.nsavc_common else None

    @property
    def nsavc_common(self):
        vals = {s.nsavc for s in self.segments}
        return vals.pop() if len(vals) == 1 else None

    # ------------------------------------------------------------------ files
    def as_dict(self):
        return {
            "md_id": self.md_id,
            "n_frames": len(self),
            "framestep_ns": self.framestep_ns,
            # Where that number came from: "declared" (a person wrote it in
            # settings.yaml), "headers" (the trajectory states it), or null
            # (nobody knows, and no time axis may be drawn).
            "framestep_basis": self.framestep_basis,
            "duration_ns": self.duration_ns(),
            # ps per INTEGRATION STEP: the page multiplies a frame's step
            # counter by this to get its time, so leaving it out made every
            # provenance panel say NaN ns.
            "step_ps": self.step_ps,
            "segments": [s.as_dict() for s in self.segments],
            "seam_raw_indices": list(self.seams),
            # THE FILE A PERSON OPENS. The analysis runs on the merged
            # trajectory, so that -- not the pr.part* segment behind it -- is
            # what a click should name: "merged.dcd, frame 2019 of 2525".
            "analysed_files": list(self.analysed),
            "analysed_path": self.analysed_path,
            "n_raw": len(self.frames) + len(self.seams),
            # HOW THE ANALYSED FILE IS INDEXED. A strip either keeps the
            # frames a restart re-simulated or drops them, and the two give
            # different frame numbers for the same picosecond. The CDR strips
            # kept them, so `merged.dcd` is the raw sequence; the
            # occluded-open ones used `trajin ... 2 last` and dropped one per
            # seam, so their file is the de-duplicated sequence and is a
            # frame shorter. Saying which one this is turns "frame 2019 of
            # 2525" from a guess into a fact.
            "n_analysed": self.n_analysed,
            "analysed_indexing": self.analysed_indexing,
            # One row per frame: [step, segment, local frame, raw index].
            # Compact on purpose -- 2,525 frames must not bloat a page.
            "frames": [[f[0], f[1], f[2], f[3]] for f in self.frames],
        }


def build(md_id, paths, framestep_ns=None):
    """A Timeline for one MD from its trajectory segments, in any order.

    Segments are ordered by their first step, ties broken by frame count: a
    window that died after one frame shares its start with the window that
    resumed there, and either order gives the same time axis (the stub's frames
    are inside the other's span) -- the same rule the metadata pipeline uses.
    """
    segs = []
    for p in paths:
        h = read_dcd_header(p)
        if "error" in h:
            raise ValueError(f"{p}: {h['error']}")
        segs.append((h["istart"], h["n_frames"], p, h))
    segs.sort(key=lambda t: (t[0], t[1]))

    segments, frames, seams, seen = [], [], [], {}
    raw = 0
    for i, (_, _, path, head) in enumerate(segs):
        seg = Segment(path, head, raw)
        segments.append(seg)
        for local in range(seg.n_frames):
            step = seg.step_of(local)
            if step in seen:
                # Same step, re-simulated after a restart: keep the first.
                seams.append(raw)
            else:
                seen[step] = True
                frames.append((step, i, local, raw))
            raw += 1
    frames.sort(key=lambda f: f[0])

    basis = "declared" if framestep_ns is not None else None
    if framestep_ns is None and segments:
        nsavc = {s.nsavc for s in segments}
        step_ps = next((s.step_ps for s in segments if s.step_ps), 0.0)
        framestep_ns = (round(step_ps * nsavc.pop() / 1000.0, 6)
                        if len(nsavc) == 1 and step_ps else None)
        basis = "headers" if framestep_ns is not None else None
        if framestep_ns is not None and framestep_ns < PLACEHOLDER_FLOOR_NS:
            # NOT A SPACING, A PLACEHOLDER. cpptraj and VMD write a merged or
            # stripped DCD with `nsavc = 1` and a default delta, because the
            # frame interval of the file they read is not theirs to state.
            # Multiplied out that is 4.9e-05 ns per frame: a 1,050 ns replica
            # reported as 0.513 ns, wrong by a factor of two thousand, with
            # nothing raised. The lab's production runs write a frame every
            # 100-200 ps and nothing here is finer than a picosecond, so a
            # derived spacing below that is the placeholder, not the run.
            framestep_ns, basis = None, None
    t = Timeline(md_id, segments, frames, seams, framestep_ns)
    t.framestep_basis = basis
    return t


def gaps(timeline):
    """Steps the timeline does NOT hold, as (after_step, before_step, ns).

    A hole means a segment was never copied or the run restarted from an
    earlier checkpoint; an analysis that averages across it is averaging across
    time that was never simulated.
    """
    out = []
    nsavc = timeline.nsavc_common
    if not nsavc:
        return out
    steps = [f[0] for f in timeline.frames]
    for a, b in zip(steps, steps[1:]):
        if b - a > nsavc:
            out.append((a, b, round((b - a) * timeline.step_ps / 1000.0, 3)))
    return out


def write(path, timelines):
    """`frames.json` -- every replica's time axis, for the pages and for scripts."""
    doc = {"generated_by": "MD_ANALYSES/core/frames.py",
           "mds": [t.as_dict() for t in timelines]}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    text = json.dumps(doc, separators=(",", ":"))
    with open(path, "w") as fh:
        fh.write(text)
    return path


if __name__ == "__main__":                                 # pragma: no cover
    import sys
    t = build(sys.argv[1], sys.argv[2:])
    print(f"{t.md_id}: {len(t)} frames, {t.duration_ns()} ns, "
          f"{len(t.segments)} segments, {len(t.seams)} seam frame(s)")
    for g in gaps(t):
        print(f"  GAP {g[2]} ns between steps {g[0]:,} and {g[1]:,}")
    for i in (0, len(t) // 2, len(t) - 1):
        print(" ", t.where(i))
