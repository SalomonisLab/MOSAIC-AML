"""Read only the rows you need out of an .h5ad CSR matrix.

The atlas pseudobulk file is ~1 GB and holds 12,255 (sample x cell-state) rows. Scoring one patient
needs ~30 of them. Materialising the whole CSR to slice 30 rows costs a gigabyte of RAM and most of a
minute, which is the difference between a usable per-patient entry point and an unusable one.

CSR makes the subset cheap: `indptr` alone locates every row's slice of `data`/`indices`. Contiguous
runs of wanted rows are read in one call each, so a sample whose rows sit together (the normal case)
costs a single read.
"""
from __future__ import annotations
import numpy as np
import h5py


def obs_column(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
        return cats[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) if v.dtype.kind == "S" else v


def var_index(f):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])


def read_rows(path, rows):
    """Dense (len(rows) x n_genes) array for the given row indices, read straight from the CSR."""
    rows = np.asarray(sorted(set(int(r) for r in rows)))
    with h5py.File(path, "r") as f:
        Xg = f["X"]
        shape = tuple(Xg.attrs["shape"])
        indptr = Xg["indptr"][:]
        out = np.zeros((len(rows), shape[1]), dtype=np.float32)
        # split the wanted rows into contiguous runs so each becomes one h5py read
        runs, start = [], 0
        for i in range(1, len(rows) + 1):
            if i == len(rows) or rows[i] != rows[i - 1] + 1:
                runs.append((start, i)); start = i
        for s, e in runs:
            lo, hi = int(indptr[rows[s]]), int(indptr[rows[e - 1] + 1])
            if hi <= lo:
                continue
            data = Xg["data"][lo:hi]
            idx = Xg["indices"][lo:hi]
            for k in range(s, e):
                a, b = int(indptr[rows[k]]) - lo, int(indptr[rows[k] + 1]) - lo
                out[k, idx[a:b]] = data[a:b]
    return rows, out
