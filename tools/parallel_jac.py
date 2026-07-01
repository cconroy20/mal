#!/usr/bin/env python3
"""Parallel finite-difference Jacobian for the RCG-forward-model fit.

The cost of a fit is (RCG ~0.9s + parse ~0.08s) x (n_params) x (n_LM_steps). The
n_params finite-difference evals per LM step are INDEPENDENT, but scipy's
least_squares computes them SERIALLY, leaving 7 of 8 cores idle. This computes the
Jacobian columns CONCURRENTLY across a pool of PROCESSES (each with its own RCG
run-dir so they don't clobber each other's ING11/OUTG11), then hands scipy an
analytic `jac` so it stops doing serial finite differences.

WHY PROCESSES, NOT THREADS: an earlier ThreadPool version only reached ~1.9x
concurrency on 8 cores -- GIL-bound. Each worker's resid_vector does substantial
Python work (write ING11, parse the 9 MB OUTG11, identity-match levels) that holds
the GIL, serializing the threads despite RCG itself being a subprocess. Separate
PROCESSES have no shared GIL, so all N run truly concurrently. Verified bit-
identical to scipy's serial 2-point forward differences (same step, same residual).

Each worker process builds its OWN Forward once (in the pool initializer), bound to
its own worker dir, from the SAME construction kwargs as the master -- so no Forward
object needs to be pickled across the process boundary (only the small (j, x) task).
"""
import os
import shutil
import concurrent.futures as cf

import numpy as np


# RCG input files copied into each worker dir (besides ING11, written per-eval).
_DECK = ["IN2", "FOR072", "FOR073", "FOR074", "SENIOR", "rwfn.dat",
         "tape2n", "TAPE2N", "in36", "ING11.abinitio", "OUTGINE.abinitio"]


def build_worker_dirs(master_dir, n):
    """Create n sibling worker dirs <master>_w{i} with the full RCG deck copied
    from master_dir. Returns the list of dir paths."""
    master_dir = os.path.abspath(master_dir)
    out = []
    for i in range(n):
        w = f"{master_dir}_w{i}"
        os.makedirs(w, exist_ok=True)
        for f in _DECK:
            src = os.path.join(master_dir, f)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(w, f))
        # seed the worker's ING11 so its first run has a valid deck
        seed = os.path.join(master_dir, "ING11")
        if os.path.exists(seed):
            shutil.copy(seed, os.path.join(w, "ING11"))
        with open(os.path.join(w, "cowan.cfg"), "w") as fh:
            fh.write(w + "/\n")
        out.append(w)
    return out


def cleanup_worker_dirs(worker_dirs):
    for w in worker_dirs:
        shutil.rmtree(w, ignore_errors=True)


# ----- per-process worker state (one Forward per process, built once) -----
_W = {}


def _worker_init(worker_dirs, fwd_kwargs, lam, ridge, scale, h):
    """Pool initializer: each process picks the NEXT unused worker dir (by a file
    lock-free claim using the process's position) and builds its own Forward there.
    We claim dirs round-robin via an OS-level atomic mkdir of a .claim marker."""
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gf_fit import Forward
    # claim a worker dir atomically (first process to mkdir <dir>/.claim wins it)
    mydir = None
    for d in worker_dirs:
        try:
            os.mkdir(os.path.join(d, ".claim"))
            mydir = d
            break
        except FileExistsError:
            continue
    if mydir is None:
        mydir = worker_dirs[0]                 # fallback (more procs than dirs)
    kw = dict(fwd_kwargs)
    kw["run_dir"] = mydir
    # the worker only needs to EVALUATE resid_vector, not re-pin targets from its
    # own seed; we rebuild Forward (which re-runs setup_targets in mydir) -- the
    # targets are identity-based and deterministic, so they match the master's.
    fwd = Forward(**kw)
    _W["fwd"] = fwd
    _W["lam"], _W["ridge"], _W["scale"], _W["h"] = lam, ridge, scale, h


def _worker_col(args):
    """Compute one Jacobian column in this process: (resid(x+h e_j) - f0)/h."""
    j, x, f0 = args
    fwd = _W["fwd"]
    xp = np.array(x, dtype=float)
    xp[j] += _W["h"]
    fp = fwd.resid_vector(xp, _W["scale"], _W["lam"], _W["ridge"])
    return j, (fp - f0) / _W["h"]


class ParallelJac:
    """Process-pool parallel finite-difference Jacobian for fit_lm.

    fwd_master: the Forward for the serial residual eval at x.
    fwd_kwargs: the exact kwargs used to build fwd_master (so each worker process
        can rebuild an identical Forward in its own dir). run_dir is overridden.
    """

    def __init__(self, fwd_master, worker_dirs, fwd_kwargs, scale, lam, ridge,
                 diff_step=1e-3, progress=None):
        self.fwd = fwd_master
        self.scale = scale
        self.lam = lam
        self.ridge = ridge
        self.h = diff_step
        self.progress = progress
        self.njac = 0
        self.nproc = len(worker_dirs)
        # clear any stale .claim markers, then start the process pool
        for d in worker_dirs:
            shutil.rmtree(os.path.join(d, ".claim"), ignore_errors=True)
        self.pool = cf.ProcessPoolExecutor(
            max_workers=self.nproc,
            initializer=_worker_init,
            initargs=(worker_dirs, fwd_kwargs, lam, ridge, scale, diff_step))

    def jac(self, x, f0=None):
        if f0 is None:
            f0 = self.fwd.resid_vector(x, self.scale, self.lam, self.ridge)
        n = len(x)
        m = len(f0)
        J = np.zeros((m, n))
        self.njac += 1
        chi2 = float(f0 @ f0)
        x = np.asarray(x, dtype=float)
        tasks = [(j, x, f0) for j in range(n)]
        done = 0
        for j, col in self.pool.map(_worker_col, tasks, chunksize=max(1, n // (4 * self.nproc))):
            J[:, j] = col
            done += 1
            if self.progress is not None and done % max(1, n // 4) == 0:
                self.progress(f"  jac {self.njac}: chi2={chi2:11.3f}  "
                              f"col {done}/{n}")
        return J

    def close(self):
        self.pool.shutdown(wait=False, cancel_futures=True)
