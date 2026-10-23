#!/usr/bin/env python3
"""
safety.py -- what this kit refuses to do on a filesystem forty people share.

Every check here answers "could this run write somewhere it was not asked to",
and every one of them runs BEFORE the first byte moves.

THE SHAPE THAT MUST STILL BE ACCEPTED. The lab-wide build lives at
`/data/<user>/MD_RESULTS`, whose *folder* is the parent of
`MD_NDP_DATA`. Containment is therefore tested against the MD_RESULTS
directory itself, never against its folder -- test the folder and the collector
refuses every run it was built for.
"""

import json
import os
import socket
import time


class Refused(Exception):
    """Something that would write outside where it was asked to."""


def _real(p):
    return os.path.realpath(os.path.abspath(p))


def inside(child, parent):
    c, p = _real(child), _real(parent)
    return c == p or c.startswith(p + os.sep)


def check_source(path):
    """A source must be an MD_RESULTS directory that a build actually wrote."""
    if not os.path.isdir(path):
        raise Refused(f"{path}: not a directory")
    if os.path.basename(path.rstrip("/")) != "MD_RESULTS":
        raise Refused(f"{path}: not an MD_RESULTS directory")
    if not os.path.exists(os.path.join(path, "settings.yaml")):
        raise Refused(f"{path}: no settings.yaml -- the metadata pipeline has "
                      f"not run here, so there is nothing to collect")
    return _real(path)


def check_destination(dest, sources):
    """The share must not sit inside a source, nor a source inside the share."""
    d = _real(dest)
    for s in sources:
        if inside(d, s):
            raise Refused(f"the destination is inside a source:\n"
                          f"      destination {d}\n      source      {s}")
        if inside(s, d):
            raise Refused(f"a source is inside the destination:\n"
                          f"      source      {s}\n      destination {d}\n"
                          f"  Collecting it would copy the share into itself.")
    return d


def check_write_target(dest, target):
    """Every path written must resolve inside the share. Catches a `../..`
    link target and any path assembled from one."""
    if not inside(target, dest):
        raise Refused(f"refusing to write outside the share:\n"
                      f"      {target}\n  is not inside\n      {_real(dest)}")
    return target


def _alive(pid):
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes
        h = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)
        if not h:
            return False
        ctypes.windll.kernel32.CloseHandle(h)
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class Lock:
    """One collector at a time. Two concurrent swaps on a shared export is a
    corrupt collection, and the window is exactly as long as a rename.

    A file created with O_CREAT|O_EXCL, not flock: flock is unreliable on many
    NFS exports -- exactly the shared storage this guards -- and does not exist
    on Windows. The file names its holder. One left by a process that died on
    this host is taken over; one held from another host is never broken.
    """

    def __init__(self, path):
        self.path = path
        self.held = False

    def __enter__(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        me = {"pid": os.getpid(), "host": socket.gethostname(),
              "since": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
        for _ in range(2):
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            except FileExistsError:
                try:
                    with open(self.path) as fh:
                        held = json.load(fh)
                except (OSError, ValueError):
                    held = {}
                if (held.get("host") == me["host"]
                        and not _alive(int(held.get("pid") or 0))):
                    try:
                        os.unlink(self.path)       # left by a collector that died
                    except FileNotFoundError:
                        pass
                    continue
                raise Refused(
                    f"another collector holds {self.path} (pid "
                    f"{held.get('pid', '?')} on {held.get('host', '?')} since "
                    f"{held.get('since', '?')}). Wait for it, or remove the lock "
                    f"if you are sure nothing is running.")
            with os.fdopen(fd, "w") as fh:
                json.dump(me, fh)
            self.held = True
            return self
        raise Refused(f"could not take {self.path}")

    def __exit__(self, *exc):
        if self.held:
            try:
                os.unlink(self.path)
            except FileNotFoundError:
                pass
            self.held = False
        return False
