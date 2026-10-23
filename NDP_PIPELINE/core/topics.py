#!/usr/bin/env python3
"""
topics.py -- propose research topics from directory names.

PROPOSES ONLY. Nothing here assigns a topic silently. substrate_init.py shows
the proposal, the human confirms or corrects, and the answer is written to
substrate.yaml. This mirrors the rule used for document status: a machine that
guesses a label at scale manufactures errors at scale.

The keyword table is the lab's own, in lab.yaml (`index_topics:`), so the
vocabulary matches what that lab actually names things.
"""

import re

try:
    from . import lab
except ImportError:              # run as a script: python3 core/topics.py <dir>
    import lab

# topic -> regexes matched case-insensitively against path components, from
# lab.yaml `index_topics:`. Ordered most-specific first; a path can match
# several topics. Bare `env` is deliberately NOT a pattern there: it matches
# python environment directories (af3_env, haddock3_env) far more often than
# HIV envelope.
TOPIC_PATTERNS = lab.index_topics()

# Directory names that carry no topic signal -- excluded before matching so
# they cannot drag in a false positive.
# Force-field parameter files carry chemistry words (carb, lipid, na_rna) that
# are not evidence of what anyone studies.
FORCEFIELD = re.compile(r"^(toppar|par_all|top_all|charmm\d*|amber\d*|ff\d+|"
                        r"gaff|tip3p|opls)", re.I)

NOISE = re.compile(
    r"^(data|analysis|scripts?|results?|tmp|temp|test|tests|old|new|backup|"
    r"figures?|figs?|plots?|docs?|input|inputs|output|outputs|run|runs|"
    r"work|misc|archive|shared|common|util|utils|bin|src|notebooks?)$", re.I)


def normalise(token):
    """Split lab-style names into words before matching.

    `\b` treats underscore as a word character, so `\bmmvt\b` does NOT match
    inside `hostguest_mmvt_input`. This lab names almost everything with
    underscores and hyphens, so without this every multi-word directory was
    silently missed. Verified: openmmvt returned only `hsp90` before this.
    """
    return re.sub(r"[_\-.]+", " ", str(token))


def propose(names, min_hits=1):
    """Given directory names, return {topic: [evidence strings]}.

    `names` should be depth-1 and depth-2 directory names -- deeper paths add
    noise without adding signal.
    """
    hits = {}
    for name in names:
        token = str(name).strip("/").split("/")[-1]
        if not token or NOISE.match(token) or FORCEFIELD.match(token):
            continue
        hay = normalise(token)
        for topic, pats in TOPIC_PATTERNS.items():
            for p in pats:
                if re.search(p, hay, re.I):
                    hits.setdefault(topic, [])
                    if token not in hits[topic]:
                        hits[topic].append(token)
                    break
    return {t: ev for t, ev in hits.items() if len(ev) >= min_hits}


def format_proposal(proposals):
    """Human-readable proposal block for substrate_init.py."""
    if not proposals:
        return ("  No topic could be inferred from directory names.\n"
                "  That is a normal outcome -- you can set topics by hand.")
    lines = []
    for topic, ev in sorted(proposals.items(), key=lambda kv: -len(kv[1])):
        shown = ", ".join(ev[:6]) + ("..." if len(ev) > 6 else "")
        lines.append(f"  {topic:<20} {len(ev):>3} match(es)   e.g. {shown}")
    return "\n".join(lines)


SKIP_EXT = {".log", ".out", ".err", ".png", ".jpg", ".svg", ".mp4", ".tmp",
            ".pyc", ".o", ".so", ".lock"}


def collect_names(base, max_depth=2, cap=4000):
    """Gather directory AND file names down to `max_depth`, for topic matching.

    Files matter: openmmvt stores its SEEKR/MMVT work as `.tgz` archives, so a
    directory-only scan reported no kinetics work at all.
    """
    import os
    names, base_depth = [], base.rstrip("/").count("/")
    for dirpath, dirnames, filenames in os.walk(base, topdown=True):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        if dirpath.count("/") - base_depth >= max_depth:
            dirnames[:] = []
        names.extend(dirnames)
        for f in filenames:
            stem, ext = os.path.splitext(f)
            if ext.lower() in SKIP_EXT:
                continue
            names.append(stem)
        if len(names) > cap:
            break
    return names


if __name__ == "__main__":
    import os, sys
    base = sys.argv[1] if len(sys.argv) > 1 else "."
    print(format_proposal(propose(collect_names(base))))
