#!/usr/bin/env python3
"""
atomic.py -- whole-file writes and a per-folder lock.

Every output here is read by something else -- the next step, the collector, a
colleague's script, a harvester -- and a file opened with "w" is empty for as
long as the writer takes. A build that dies, or a reader that arrives early,
then sees a truncated record, catalogue or parquet. Everything is written to a
temporary name beside the target and renamed into place, which is atomic on
POSIX and on NTFS: a reader sees the old file or the new one, never half.

The lock is a file created with O_CREAT|O_EXCL. Unlike flock it is atomic on
NFS (v3 and later) and works on Windows, and the file says who holds it. A lock
left by a process that died on this host is taken over; one held from another
host is reported, never broken, because nothing here can tell whether that
process is alive.
"""

import contextlib
import json
import os
import socket
import time
import uuid


def atomic(path, writer, fsync=False):
    """Call `writer(tmp)` on a temporary name beside `path`, then rename.

    The temporary name carries the pid and a random suffix, so two writers
    never share one: a fixed `.tmp` name let concurrent runs clobber each
    other's half-written file.
    """
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    tmp = f"{path}.tmp-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    try:
        writer(tmp)
        if fsync:
            with open(tmp, "rb") as fh:
                os.fsync(fh.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return path


def write_text(path, text, newline=None, encoding=None):
    def _w(tmp):
        with open(tmp, "w", newline=newline, encoding=encoding) as fh:
            fh.write(text)
    return atomic(path, _w)


def write_json(path, obj, **dump_kw):
    dump_kw.setdefault("indent", 2)
    return write_text(path, json.dumps(obj, **dump_kw) + "\n")


def to_parquet(df, path, **kw):
    return atomic(path, lambda tmp: df.to_parquet(tmp, **kw))


class Locked(RuntimeError):
    """Another live process holds the folder's lock."""


def _alive(pid):
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes                                        # noqa: PLC0415
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


@contextlib.contextmanager
def lock(path, what="this folder"):
    """Hold `path` as a lock file for the duration of the block.

    Raises Locked, naming the holder, if another live process has it.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    me = {"pid": os.getpid(), "host": socket.gethostname(),
          "since": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    for _ in range(2):
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            try:
                with open(path) as fh:
                    held = json.load(fh)
            except (OSError, ValueError):
                held = {}
            same_host = held.get("host") == me["host"]
            if same_host and not _alive(int(held.get("pid") or 0)):
                with contextlib.suppress(FileNotFoundError):
                    os.unlink(path)                # left by a process that died
                continue
            raise Locked(
                f"{what} is being built by pid {held.get('pid', '?')} on "
                f"{held.get('host', '?')} since {held.get('since', '?')}.\n"
                f"  Wait for it to finish. If that process is gone, delete {path}")
        with os.fdopen(fd, "w") as fh:
            json.dump(me, fh)
        break
    else:
        raise Locked(f"could not take {path}")
    try:
        yield path
    finally:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(path)


_HELD = {}


def hold(path, what="this folder"):
    """Take `path` for the rest of this process; released at exit.

    For the command-line steps, which each own a folder from start to finish.
    Re-entrant within one process, so a test (or a step) that calls another
    step's main() does not lock itself out.
    """
    path = os.path.abspath(path)
    if path in _HELD:
        return path
    import atexit                                            # noqa: PLC0415
    cm = lock(path, what)
    cm.__enter__()
    _HELD[path] = cm
    atexit.register(release, path)
    return path


def release(path):
    cm = _HELD.pop(os.path.abspath(path), None)
    if cm is not None:
        cm.__exit__(None, None, None)
