#!/usr/bin/env python3
"""Give the trained drug model the cohort-matched references it needs to score single-cell data.

Model A is fitted on BeatAML bulk. Handed a single-cell bulk-equivalent, or a single cell state's
pseudobulk, its raw decision scores land in a different band of the same scale -- so a probability
calibrated on BeatAML raw scores returns ~0.99 for essentially every inhibitor. The fix is the one the
mutation caller already uses: percentile each sample against samples OF ITS OWN KIND, then map that
percentile (which is cohort-invariant) through the calibration curve.

This builds and attaches three things:

  expression reference `sc`         gene-wise mean/sd over the 387 single-cell bulk-equivalents
  expression reference `sc_state`   the same over per-(sample x cell-state) pseudobulks
  score references `sc_sample`,     the model's own predicted-score distribution per inhibitor for
                  `sc_state`        whole single-cell samples and for individual cell states

Cell states get their own reference because a single state's profile is systematically more extreme
than any whole sample: scored against the sample-level distribution every state looks like an outlier.

  python build_drug_score_refs.py  ->  rewrites pipeline/drug_response_model.pkl in place
"""
import os, sys, time, pickle, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, h5py

from amlmm.drug import statemodel as SM
from amlmm.drug.h5rows import obs_column, var_index, read_rows

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MODEL = os.path.join(HERE, "drug_response_model.pkl")
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
H5 = os.path.join(ROOT, "data", "RNA", "pseudobulk_counts_hashed.h5ad")
STATE_COL = "Hs-BM-titrated-reference-centroid"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-state-rows", type=int, default=6000,
                    help="cap on per-(sample x cell-state) rows used for the state reference")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()

    with open(MODEL, "rb") as f:
        mod = pickle.load(f)
    fs = mod.fs
    d = np.load(BUNDLE, allow_pickle=True)

    # ---------- 1. whole single-cell samples (the bundle's bulk-equivalents) ----------
    sc_X = d["sc_X"].astype(float)
    Zsc = fs.add_reference("sc", sc_X)                    # gene-wise mean/sd over single-cell samples
    P = mod.predict_matrix(Zsc)
    mod.add_score_reference("sc_sample", {c: P[c].values for c in P.columns})
    print("sc_sample reference: %d samples x %d inhibitors (%.0fs)" % (P.shape[0], P.shape[1], time.time() - t0))

    # ---------- 2. individual cell states ----------
    with h5py.File(H5, "r") as f:
        ds, sm = obs_column(f, "Dataset"), obs_column(f, "Sample")
        st = obs_column(f, STATE_COL)
        nc = np.asarray(obs_column(f, "n_cells"), dtype=float)
        genes = var_index(f)
        n_rows = tuple(f["X"].attrs["shape"])[0]
    keep = np.where(nc >= SM.MIN_CELLS)[0]
    rng = np.random.RandomState(a.seed)
    if len(keep) > a.max_state_rows:
        keep = np.sort(rng.choice(keep, a.max_state_rows, replace=False))
    _, X = read_rows(H5, keep)
    print("cell-state rows read: %d of %d (>= %d cells) (%.0fs)"
          % (len(keep), n_rows, SM.MIN_CELLS, time.time() - t0))

    lin = SM.cp10k(X)
    aligned = pd.DataFrame(0.0, index=range(len(keep)), columns=fs.genes)
    gset = set(fs.genes)
    cols = {}
    for j, g in enumerate(genes):
        g = str(g)
        k = g if g in gset else mod.sym2ens.get(g)
        if k in gset and k not in cols:
            cols[k] = j
    take = [g for g in fs.genes if g in cols]
    aligned[take] = lin[:, [cols[g] for g in take]]
    Zst = fs.add_reference("sc_state", aligned.values)
    Ps = mod.predict_matrix(Zst)
    mod.add_score_reference("sc_state", {c: Ps[c].values for c in Ps.columns})
    print("sc_state reference: %d cell-state pseudobulks x %d inhibitors (%.0fs)"
          % (Ps.shape[0], Ps.shape[1], time.time() - t0))

    # ---------- 3. an OOD reference for single-cell samples ----------
    # Every single-cell sample is "far" from BeatAML because the assays differ, so penalising a patient
    # for that would just subtract a constant from everyone. The useful question is whether this
    # patient is unusual AMONG single-cell AML samples, which needs its own distance distribution.
    if getattr(mod, "nn", None):
        Psc = fs.pca.transform(Zsc[:, fs.sel])
        mu, sd = Psc.mean(0), Psc.std(0); sd[sd == 0] = 1.0
        mod.nn["sc_mu"], mod.nn["sc_sd"] = mu, sd
        mod.nn["ood_ref_sc"] = np.sqrt((((Psc - mu) / sd) ** 2).mean(1))
        print("sc OOD reference: median distance %.2f (BeatAML %.2f)"
              % (float(np.median(mod.nn["ood_ref_sc"])), float(np.median(mod.nn["ood_ref"]))))

    mod.score_ref_meta = {
        "built": time.strftime("%Y-%m-%d %H:%M"),
        "sc_sample_n": int(P.shape[0]), "sc_state_n": int(Ps.shape[0]),
        "expression_references": sorted(fs.ref),
        "score_references": sorted(mod.score_refs),
        "genes_matched_state_matrix": len(take),
    }
    with open(MODEL, "wb") as f:
        pickle.dump(mod, f)

    # a quick sanity read-out: are the three score scales actually different?
    print("\n%-26s %10s %10s %10s" % ("inhibitor", "BA median", "scSample", "scState"))
    for drug in ["Venetoclax", "Gilteritinib", "Trametinib (GSK1120212)", "Azacytidine", "Cytarabine"]:
        r = mod.score_refs
        def med(k):
            v = (r.get(k) or {}).get(drug)
            return "n/a" if v is None else "%.3f" % float(np.median(v))
        print("%-26s %10s %10s %10s" % (drug[:26], med("beataml"), med("sc_sample"), med("sc_state")))
    print("\nrewrote %s (%.0fs)" % (MODEL, time.time() - t0))


if __name__ == "__main__":
    main()
