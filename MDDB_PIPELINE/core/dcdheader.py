#!/usr/bin/env python3
"""
dcdheader.py -- read a DCD header without loading coordinates.

A 3.35 GB stripped DCD takes minutes to load with mdtraj and ~200 bytes to
interrogate this way. The frame count, atom count and timestep all live in the
first block, so an 8 TB campaign can be inventoried in minutes.

This follows the project's own rule -- read the value out of the data, not out
of a document. Frame counts quoted in prose have been wrong here before.

Layout (CHARMM/NAMD DCD):
    int32  84                     block size
    char4  'CORD'
    int32  icntrl[20]             [0]=NSET frames  [1]=ISTART  [2]=NSAVC
                                  [3]=NSTEP  [9]=DELTA(float32)  [10]=unitcell
    int32  84
    int32  titlesize              title block
    ...    titlesize bytes
    int32  titlesize
    int32  4
    int32  NATOM
    int32  4
"""

import os
import struct

# 1 AKMA time unit in picoseconds (NAMD/CHARMM internal time unit).
AKMA_PS = 0.04888821


def read_dcd_header(path):
    """Return a dict of header fields, or {'error': ...} if unreadable.

    Never raises -- a corrupt or truncated trajectory must not stop an
    inventory of thousands of them.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(92)
            if len(head) < 92:
                # An EMPTY file is a fact about the DATA, not a failure of this
                # reader: a job that died before writing a frame leaves a 0-byte
                # trajectory, and in an active workspace that happens often.
                # Reporting it as a reader error made the acceptance suite go
                # red for something the toolkit got exactly right -- and buried
                # the useful signal, which is "one of your runs produced
                # nothing". Distinguish them so each can be acted on.
                if not head:
                    return {"error": "empty file"}
                return {"error": "truncated header"}

            endian = "<"
            (blk,) = struct.unpack(endian + "i", head[0:4])
            if blk != 84:
                endian = ">"
                (blk,) = struct.unpack(endian + "i", head[0:4])
                if blk != 84:
                    return {"error": f"bad first block size {blk}"}

            if head[4:8] != b"CORD":
                return {"error": "missing CORD magic"}

            icntrl = struct.unpack(endian + "20i", head[8:88])
            nset, istart, nsavc, nstep = icntrl[0], icntrl[1], icntrl[2], icntrl[3]
            charmm_version = icntrl[19]
            # CHARMM-flavour flags only: in the X-PLOR form icntrl[9:11] hold
            # the float64 DELTA, so these slots are not flags there.
            namnf = icntrl[8] if charmm_version else 0     # fixed atoms
            has_unitcell = bool(icntrl[10]) if charmm_version else False
            has_4d = bool(icntrl[11]) if charmm_version else False

            # DELTA sits in icntrl[9]: float32 when CHARMM-flavoured (version
            # non-zero, which NAMD also sets), float64 in the old X-PLOR form.
            if charmm_version:
                (delta,) = struct.unpack(endian + "f", head[8 + 9 * 4: 8 + 10 * 4])
            else:
                (delta,) = struct.unpack(endian + "d", head[8 + 9 * 4: 8 + 11 * 4])

            # Title block, then the atom-count block.
            fh.seek(92)
            raw = fh.read(4)
            if len(raw) < 4:
                return {"error": "truncated title block"}
            (titlesize,) = struct.unpack(endian + "i", raw)
            # A title block is 4 bytes (NTITLE) plus 80 per line. Anything else
            # is a corrupt header, and seeking by it would read garbage as NATOM.
            if titlesize < 4 or (titlesize - 4) % 80 or titlesize > 4 + 80 * 1000:
                return {"error": f"bad title block size {titlesize}"}
            fh.seek(titlesize, 1)
            raw = fh.read(4)
            if len(raw) < 4 or struct.unpack(endian + "i", raw)[0] != titlesize:
                return {"error": "title block markers disagree"}
            raw = fh.read(8)
            if len(raw) < 8:
                return {"error": "truncated natom block"}
            _, natom = struct.unpack(endian + "ii", raw)
            fh.seek(4, 1)                      # natom block's trailing marker
            if natom <= 0 or not 0 <= namnf < natom:
                return {"error": f"bad atom counts natom={natom} fixed={namnf}"}
            # With fixed atoms the header carries the free-atom index list, and
            # only the FIRST frame stores every atom; later frames store the
            # free ones. Sizing every frame at NATOM undercounted such files and
            # flagged them truncated.
            if namnf:
                fh.seek((natom - namnf) * 4 + 8, 1)

            # HOW MANY FRAMES ARE REALLY THERE. `nset` is what the writer
            # INTENDED; the file size says what arrived. A trajectory still
            # being copied from the cluster (CDR 2-JL427/3R, 2026-09-22:
            # pr.part8.dcd, header 78 frames, 3 on disk) would otherwise be
            # published at its intended length. The smaller of the two is the
            # honest answer; a file with trailing bytes keeps `nset`.
            after_header = fh.tell()
            ncoord = 4 if has_4d else 3
            cell = 56 if has_unitcell else 0
            first_bytes = cell + ncoord * (natom * 4 + 8)
            frame_bytes = cell + ncoord * ((natom - namnf) * 4 + 8)
            on_disk = None
            try:
                body = os.path.getsize(path) - after_header
                on_disk = (0 if body < first_bytes
                           else 1 + (body - first_bytes) // frame_bytes)
            except OSError:
                on_disk = None
            truncated = on_disk is not None and on_disk < nset
            # NSET 0 with frames on disk: the writer died before it rewrote the
            # header (NAMD fills NSET in as it goes). The frames are real.
            unset = nset == 0 and bool(on_disk)
            dt_ps = delta * AKMA_PS * nsavc if (delta and nsavc) else None
            return {
                "n_frames": on_disk if (truncated or unset) else nset,
                "n_frames_header": nset,
                "n_frames_on_disk": on_disk,
                "truncated": truncated,
                "nset_unset": unset,
                "n_fixed": namnf,
                "has_4d": has_4d,
                "natom": natom,
                "istart": istart,
                "nsavc": nsavc,
                "nstep": nstep,
                "delta_akma": delta,
                "dt_ps_header": round(dt_ps, 6) if dt_ps else None,
                "has_unitcell": has_unitcell,
                "charmm_version": charmm_version,
                "endian": "little" if endian == "<" else "big",
            }
    except Exception as exc:                       # noqa: BLE001
        return {"error": f"{type(exc).__name__}: {exc}"}


if __name__ == "__main__":
    import json, sys
    for p in sys.argv[1:]:
        print(p)
        print(json.dumps(read_dcd_header(p), indent=2))
