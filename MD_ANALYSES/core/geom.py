#!/usr/bin/env python3
"""
geom.py -- the geometry and statistics kernels, lifted from
ABZYME/3-OCCLUDED-OPEN/CDR/ANALYSIS/geom.py (2026-09), which is generic and
correct. Only its `import config` (there purely to pin BLAS threads) is
replaced, by core.threads.

Kept verbatim from that module, because the reasoning is not obvious and was
paid for once already: frames are built by Gram-Schmidt from two physically
defined vectors, never from inertia eigenvectors -- inertia axes carry an
arbitrary sign and swap whenever two moments approach degeneracy, which WILL
happen for a near-cylindrical Fv, injecting spurious 180 degree jumps into
every orientation timeseries.
"""

from . import threads          # noqa: F401  (side effect: pins BLAS threads)

import numpy as np


def kabsch(P, Q, w=None):
    """Rotation and translation taking P onto Q; returns (R, t, rmsd)."""
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    if w is None:
        w = np.ones(len(P))
    w = np.asarray(w, dtype=np.float64)
    wn = w / w.sum()
    pc = (P * wn[:, None]).sum(0)
    qc = (Q * wn[:, None]).sum(0)
    P0, Q0 = P - pc, Q - qc
    H = (P0 * wn[:, None]).T @ Q0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    rmsd = float(np.sqrt(max(0.0, ((P0 @ R.T - Q0) ** 2 * wn[:, None]).sum())))
    return R, qc - R @ pc, rmsd


def rmsd_after_fit(X, ref, fit_idx, meas_idx=None):
    """RMSD of `meas_idx` after superposing `fit_idx` onto the reference.

    Fitting and measuring on different sets is the whole point: a ligand's RMSD
    means nothing unless the receptor it sits on was aligned first.
    """
    R, t, _ = kabsch(X[fit_idx], ref[fit_idx])
    Y = X @ R.T + t
    idx = meas_idx if meas_idx is not None else fit_idx
    d = Y[idx] - ref[idx]
    return float(np.sqrt((d * d).sum() / len(idx))), Y


def radius_of_gyration(X, masses=None):
    if masses is None:
        c = X.mean(0)
        return float(np.sqrt(((X - c) ** 2).sum(1).mean()))
    m = np.asarray(masses, dtype=np.float64)
    c = (X * m[:, None]).sum(0) / m.sum()
    return float(np.sqrt((m * ((X - c) ** 2).sum(1)).sum() / m.sum()))


class RunningRMSF:
    """Per-atom fluctuation without holding every frame in memory.

    A 2,520-frame x 55,818-atom trajectory is 1.7 GB of float32; nine of them
    do not fit in a comfortable process. Welford over the fitted coordinates
    costs two arrays and is exact.
    """

    def __init__(self, n_atoms):
        self.n = 0
        self.mean = np.zeros((n_atoms, 3))
        self.m2 = np.zeros((n_atoms, 3))

    def add(self, X):
        self.n += 1
        d = X - self.mean
        self.mean += d / self.n
        self.m2 += d * (X - self.mean)

    def rmsf(self):
        if self.n < 2:
            return np.zeros(len(self.mean))
        return np.sqrt((self.m2 / (self.n - 1)).sum(1))


def run_lengths(mask):
    """Start and length of every True run: the kernel behind every lifetime.

    Used for contact lifetimes, water residence and ion dwell alike, so those
    three cannot disagree about what "a run" is.
    """
    m = np.asarray(mask, dtype=bool)
    if not m.any():
        return np.empty((0, 2), dtype=int)
    d = np.diff(m.astype(np.int8))
    starts = np.flatnonzero(d == 1) + 1
    ends = np.flatnonzero(d == -1) + 1
    if m[0]:
        starts = np.r_[0, starts]
    if m[-1]:
        ends = np.r_[ends, len(m)]
    return np.c_[starts, ends - starts]


def blocking_curve(x, max_levels=14):
    """Flyvbjerg-Petersen blocking: sem estimate vs block level."""
    x = np.asarray(x, dtype=np.float64)
    out = []
    for _ in range(max_levels):
        n = len(x)
        if n < 4:
            break
        out.append(float(np.std(x, ddof=1) / np.sqrt(n)))
        x = 0.5 * (x[0:2 * (n // 2):2] + x[1:2 * (n // 2):2])
    return np.array(out)


def plateau_sem(x, max_levels=14):
    """The plateau of the blocking curve: the correlation-corrected s.e.m."""
    c = blocking_curve(x, max_levels)
    return float(c.max()) if len(c) else float("nan")


def autocorr_time(x, cut=None):
    """Integrated autocorrelation time, in samples."""
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    if not np.any(x):
        return 1.0
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conjugate(f))[:n].real
    ac /= ac[0]
    stop = cut or n
    tau, i = 1.0, 1
    while i < stop and ac[i] > 0:
        tau += 2.0 * ac[i] * (1.0 - i / n)
        if i > 6 * tau:
            break
        i += 1
    return float(tau)


def n_effective(x):
    """Independent samples in a correlated series."""
    t = autocorr_time(x)
    return float(len(x) / max(1.0, t))
