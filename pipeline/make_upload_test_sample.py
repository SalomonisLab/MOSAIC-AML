#!/usr/bin/env python3
"""Build a realistic single-cell upload from an atlas sample, so the upload path can be tested for real.

The point of the test is not "does the code run" but "does an uploaded sample come back with the right
answer". So the query is constructed with a KNOWN ground truth: one atlas sample's per-(cell-state)
pseudobulk is expanded into individual cells by Poisson resampling, in proportion to how many cells
that state actually had. The pipeline is then asked to recover the composition it was built from --
which it has to do through cosine assignment against the cellHarmony reference, never seeing the
labels.

  python make_upload_test_sample.py --sample "CCHMC::1009_AfInv16_29M" --cells 3000
      -> inbox/UPLOAD_TEST_<sample>.h5ad          the "uploaded" file
      -> inbox/UPLOAD_TEST_<sample>.truth.json    the composition it was built from
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, h5py

from amlmm.drug.h5rows import obs_column, var_index, read_rows

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
H5 = os.path.join(ROOT, "data", "RNA", "pseudobulk_counts_hashed.h5ad")
STATE_COL = "Hs-BM-titrated-reference-centroid"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="CCHMC::1009_AfInv16_29M")
    ap.add_argument("--cells", type=int, default=3000)
    ap.add_argument("--depth", type=int, default=4000, help="counts per synthetic cell")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "inbox"))
    a = ap.parse_args()

    with h5py.File(H5, "r") as f:
        key = np.array(["%s::%s" % (x, y) for x, y in zip(obs_column(f, "Dataset"), obs_column(f, "Sample"))])
        rows = np.where(key == a.sample)[0]
        if not len(rows):
            raise SystemExit("sample %r not in the atlas" % a.sample)
        states = obs_column(f, STATE_COL)[rows]
        ncell = np.asarray(obs_column(f, "n_cells"), dtype=float)[rows]
        genes = var_index(f)
    _, X = read_rows(H5, rows)

    # how many synthetic cells each state gets, in proportion to its real abundance
    w = ncell / ncell.sum()
    alloc = np.maximum(1, np.round(w * a.cells).astype(int))
    rng = np.random.RandomState(a.seed)

    cells, labels = [], []
    for i, st in enumerate(states):
        p = X[i].astype(np.float64)
        s = p.sum()
        if s <= 0:
            continue
        p = p / s
        for _ in range(int(alloc[i])):
            cells.append(rng.poisson(p * a.depth))          # a cell = a Poisson draw at that depth
            labels.append(str(st))
    M = np.vstack(cells).astype(np.float32)
    keep = M.sum(0) > 0                                     # drop genes zero in every synthetic cell
    M, gk = M[:, keep], genes[keep]
    truth = pd.Series(labels).value_counts(normalize=True).sort_values(ascending=False)

    os.makedirs(a.out_dir, exist_ok=True)
    tag = a.sample.replace("::", "_").replace("/", "_")
    h5out = os.path.join(a.out_dir, "UPLOAD_TEST_%s.h5ad" % tag)
    import anndata as ad
    ad.AnnData(X=M, obs=pd.DataFrame(index=["cell%05d" % i for i in range(M.shape[0])]),
               var=pd.DataFrame(index=[str(g) for g in gk])).write_h5ad(h5out)
    json.dump({"source_sample": a.sample, "n_cells": int(M.shape[0]), "n_genes": int(M.shape[1]),
               "composition": {k: round(float(v), 4) for k, v in truth.items()}},
              open(h5out.replace(".h5ad", ".truth.json"), "w"), indent=1)
    print("wrote %s  (%d cells x %d genes, %d states)" % (h5out, M.shape[0], M.shape[1], len(truth)))
    print("top states built in:", ", ".join("%s %.0f%%" % (k, 100 * v) for k, v in truth.head(6).items()))


if __name__ == "__main__":
    main()
