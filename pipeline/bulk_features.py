#!/usr/bin/env python3
"""Whole-sample 'bulk RNA' feature block, synthesized from the scRNA atlas.

The cohort has NO paired bulk RNA-seq — every RNA asset is per-(sample x cell-state) pseudobulk. To compare
single-cell-resolved prediction against a plain bulk RNA classifier on the SAME samples/labels, we synthesize
a bulk-equivalent: sum ALL of a sample's cell-state pseudobulks into one whole-sample raw-count vector
(discarding cell-state structure, exactly what bulk RNA-seq gives), then CP10k + log1p.

  bulk_rna_matrix(ctx) -> DataFrame  (rows = sample_key 'Dataset::Sample', cols = genes)

Correct order of operations: combine RAW counts first (sum across cell-states), THEN normalize (CP10k) and
log (log1p) on the aggregate — never sum pre-normalized/logged per-cell-state values.
"""
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix


def bulk_rna_matrix(ctx):
    path = ctx._modality_paths["RNA"]
    with h5py.File(path, "r") as f:
        def col(n):
            g = f["obs"][n]
            if isinstance(g, h5py.Group):
                cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
                return cats[g["codes"][:]]
            v = g[:]
            return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) if v.dtype.kind == "S" else v
        ds, sm = col("Dataset"), col("Sample")
        genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]])
        Xg = f["X"]; shape = tuple(Xg.attrs["shape"])
        X = csr_matrix((Xg["data"][:], Xg["indices"][:], Xg["indptr"][:]), shape=shape)  # rows=(sample x state) raw counts
    skey = np.array(["%s::%s" % (d, s) for d, s in zip(ds, sm)])
    samples = sorted(set(skey))
    sidx = {k: j for j, k in enumerate(samples)}
    rows = np.array([sidx[k] for k in skey])
    # indicator (n_samples x n_pseudobulk_rows) @ X  ->  per-sample SUM of raw counts across ALL cell-states
    M = csr_matrix((np.ones(len(skey)), (rows, np.arange(len(skey)))), shape=(len(samples), X.shape[0]))
    bulk = np.asarray((M @ X).todense(), dtype=np.float64)                     # sample x gene, raw summed counts
    tot = bulk.sum(1); tot[tot == 0] = 1.0
    bulk = np.log1p(bulk / tot[:, None] * 1e4)                                 # CP10k + log1p AFTER aggregation
    df = pd.DataFrame(bulk, index=samples, columns=genes)
    df = df.loc[:, (df.values > 0).any(0)]                                     # drop genes zero in every sample
    return df


if __name__ == "__main__":
    import os, sys, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from amlmm.context import build_context, Config
    t = time.time(); m = bulk_rna_matrix(build_context(Config(run_id="bulkprobe")))
    print("bulk RNA matrix:", m.shape, "in %.1fs" % (time.time() - t), "| e.g.", list(m.index[:3]), list(m.columns[:3]))
