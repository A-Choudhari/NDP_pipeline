#!/usr/bin/env python3
"""
threads.py -- import this FIRST, before numpy, in every entry point.

Lifted from ANALYSIS/config.py's header. BLAS spawns one thread per core by
default; with several analysis processes on a 48-thread box that is 48 x 48
threads fighting over the same caches, and it makes a pooled stage SLOWER than
a serial one. Every worker here is single-threaded by contract, and
parallelism comes from the process pool -- which is the right granularity when
the bottleneck is a 94 MB/s NFS mount, not arithmetic.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
