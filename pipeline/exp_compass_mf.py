#!/usr/bin/env python3
"""A3 — multi-task matrix factorisation across inhibitors.

The deployed model fits a ridge per target-pathway family plus a shrunk per-drug residual. That shares
information only through hand-drawn family boundaries. The response matrix is 520 x 118 and ~87%
observed, which is exactly the regime where a low-rank factorisation shares strength across ALL drugs
and estimates a usable profile for the low-n ones.

    Y_ij  ~  u_i . v_j          u_i = patient latent factors,  v_j = drug latent factors

To generalise to a NEW patient (who has no response row), the patient factors are regressed on the
patient's molecular features, so at prediction time  y_hat_ij = f(x_i) . v_j.

Everything is fit inside the fold: the SVD, the feature->factor regression, and the within-drug
normalisation. Scored on BOTH targets — the current within-drug z, and the double-centred interaction
target that removes the patient main effect and is the honest measure of drug-specific skill.

  python exp_compass_mf.py [--rank 12] -> deliverables/exp_compass_mf.json
"""
import os, sys, json, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold

from amlmm.drug import data as D
from eval_drug_model import load_all
from train_drug_model import build_space, fit_predict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_compass_mf.json")


def soft_impute_svd(Y, rank, n_iter=60, tol=1e-4):
    """Rank-`rank` factorisation of a matrix with missing entries, by iterative SVD imputation.

    Missing entries start at zero (the matrix is already centred, so zero = "no information") and are
    refilled from the current low-rank reconstruction each pass. Simple, deterministic, and adequate at
    this size — no need for an optimiser.
    """
    obs = np.isfinite(Y)
    X = np.where(obs, Y, 0.0)
    prev = None
    for _ in range(n_iter):
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        Xl = (U[:, :rank] * s[:rank]) @ Vt[:rank]
        X = np.where(obs, Y, Xl)
        if prev is not None and np.linalg.norm(Xl - prev) / (np.linalg.norm(prev) + 1e-9) < tol:
            break
        prev = Xl
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :rank] * s[:rank], Vt[:rank]              # (patient factors, drug factors)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", type=int, default=12)
    ap.add_argument("--folds", type=int, default=5)
    a = ap.parse_args()
    t0 = time.time()

    A = load_all()
    long, row_of, drugs = A["long"], A["row_of"], A["drugs"]
    specs = sorted(long["specimen"].unique())
    sidx = {s: i for i, s in enumerate(specs)}
    didx = {d: j for j, d in enumerate(drugs)}
    subj = long.drop_duplicates("specimen").set_index("specimen")["subject"].to_dict()
    groups = np.array([subj[s] for s in specs])

    # raw AUC matrix; normalisation happens INSIDE each fold
    AUC = np.full((len(specs), len(drugs)), np.nan)
    for s, dr, v in zip(long["specimen"], long["inhibitor"], long["auc"]):
        if s in sidx and dr in didx:
            AUC[sidx[s], didx[dr]] = v

    pred_mf = np.full_like(AUC, np.nan)
    pred_dep = np.full_like(AUC, np.nan)
    ytar = np.full_like(AUC, np.nan)          # within-drug z (current target)
    yint = np.full_like(AUC, np.nan)          # double-centred (interaction target)

    gk = GroupKFold(n_splits=a.folds)
    for k, (tri, tei) in enumerate(gk.split(np.zeros(len(specs)), groups=groups)):
        # --- normalise on TRAIN specimens only -----------------------------------------
        med = np.nanmedian(AUC[tri], axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(AUC[tri] - med), axis=0)
        mad[~np.isfinite(mad) | (mad < 1e-9)] = np.nanstd(AUC[tri], axis=0)[~np.isfinite(mad) | (mad < 1e-9)]
        Z = -(AUC - med) / mad                                   # higher = more sensitive
        pm = np.nanmedian(Z[tri], axis=1)                        # patient main effect, train-derived
        Zi = Z - np.nanmedian(Z, axis=1, keepdims=True)          # interaction target
        ytar[tei] = Z[tei]; yint[tei] = Zi[tei]

        # --- A3: factorise the TRAIN block, learn features -> patient factors ----------
        Utr, V = soft_impute_svd(Z[tri], a.rank)
        rows_tr = sorted({row_of[specs[i]] for i in tri})
        fs, Zx = build_space(A["X"], rows_tr, A["meta"], A["mut"], 100, 4000, A["genes"])
        Xp, _, _ = fs.transform(Zx, meta=A["meta"], mut=A["mut"], blocks=["rna", "state", "mut", "clin"])
        Xs = np.vstack([Xp[row_of[specs[i]]] for i in range(len(specs))])
        reg = RidgeCV(alphas=(10.0, 100.0, 1000.0, 1e4)).fit(Xs[tri], Utr)
        pred_mf[tei] = reg.predict(Xs[tei]) @ V

        # --- deployed model, same fold, for a like-for-like comparison ------------------
        f_tr = long[long["specimen"].isin({specs[i] for i in tri})]
        f_te = long[long["specimen"].isin({specs[i] for i in tei})]
        mod, ltr, lte, p, _ = fit_predict(f_tr, f_te, Zx, row_of, drugs, fs, A["sym2ens"],
                                          ["rna", "state", "mut", "clin"], A["meta"], A["mut"], 0.20)
        for s, dr, v in zip(lte["specimen"], lte["inhibitor"], p.loc[lte.index]):
            if s in sidx and dr in didx and v == v:
                pred_dep[sidx[s], didx[dr]] = v
        print("   fold %d/%d (%.0fs)" % (k + 1, a.folds, time.time() - t0))

    # --- score ----------------------------------------------------------------------
    def per_drug(pred, target):
        sp, au = [], []
        for j, dr in enumerate(drugs):
            m = np.isfinite(pred[:, j]) & np.isfinite(target[:, j])
            if m.sum() < 50:
                continue
            y, p = target[m, j], pred[m, j]
            sp.append(float(spearmanr(y, p).statistic))
            lo, hi = np.quantile(y, [0.2, 0.8])
            yy = np.where(y >= hi, 1, np.where(y <= lo, 0, -1)); kk = yy >= 0
            if kk.sum() >= 20 and len(set(yy[kk])) == 2:
                au.append(float(roc_auc_score(yy[kk], p[kk])))
        return float(np.mean(sp)), float(np.mean(au)), len(sp)

    blend = np.where(np.isfinite(pred_mf) & np.isfinite(pred_dep),
                     0.5 * np.nan_to_num(pred_mf) + 0.5 * np.nan_to_num(pred_dep), np.nan)
    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "rank": a.rank, "arms": {}}
    print("\n== A3 · multi-task matrix factorisation (rank %d) ==" % a.rank)
    print("  %-22s %-24s %-24s" % ("", "current target (within-drug z)", "interaction target (honest)"))
    print("  %-22s %10s %10s   %10s %10s" % ("arm", "Spearman", "AUROC", "Spearman", "AUROC"))
    for nm, P in (("deployed (per-family)", pred_dep), ("matrix factorisation", pred_mf),
                  ("blend 50/50", blend)):
        s1, a1, n1 = per_drug(P, ytar)
        s2, a2, _ = per_drug(P, yint)
        res["arms"][nm] = {"n_drugs": n1, "spearman_current": round(s1, 4), "auroc_current": round(a1, 4),
                           "spearman_interaction": round(s2, 4), "auroc_interaction": round(a2, 4)}
        print("  %-22s %10.3f %10.3f   %10.3f %10.3f" % (nm, s1, a1, s2, a2))

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
