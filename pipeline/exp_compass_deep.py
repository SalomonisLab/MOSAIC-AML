#!/usr/bin/env python3
"""Deep characterisation of COMPASS-AML, and the matrix-factorisation blend decision.

COMPASS scores 0.672 against a measured assay ceiling of 0.727, so there is very little headroom and
the honest question is not "can it go higher" but "where exactly does it stand, and is the one
untried modelling change worth deploying". Five experiments:

  C1  rank sweep          matrix factorisation at ranks 5/12/25/40/60, scored on the INTERACTION
                          target (patient main effect removed) -- the only honest measure of
                          drug-SPECIFIC skill. Includes the deployed per-family model and the blend.
  C2  reliability vs      does per-drug predictability track the assay's own reproducibility? If it
      predictability      does, the ceiling argument is supported by the per-drug data, not just the
                          pooled median. This number has been quoted before without a saved artefact.
  C3  per-drug table      every inhibitor: n, CV AUROC, Spearman, assay reliability, tier, ceiling,
                          and the fraction of its ceiling attained. The deliverable a pharmacologist
                          would actually want.
  C4  blend weight        the 50/50 blend was arbitrary. Sweep it, on the interaction target, with
                          the weight chosen inside the fold rather than on the score being reported.
  C5  drug-count curve    AUROC vs number of training specimens per drug, to separate "hard drug"
                          from "not enough data".

  python exp_compass_deep.py [--folds 5]  ->  deliverables/exp_compass_deep.json
"""
import os, sys, json, time, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold

from amlmm.drug import data as DD
from eval_drug_model import load_all
from train_drug_model import build_space, fit_predict
from exp_compass_mf import soft_impute_svd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "deliverables", "exp_compass_deep.json")


def per_drug(pred, target, drugs, min_n=50):
    """Spearman and tail-class AUROC per drug, on whichever target is passed."""
    rows = {}
    for j, dr in enumerate(drugs):
        m = np.isfinite(pred[:, j]) & np.isfinite(target[:, j])
        if m.sum() < min_n:
            continue
        y, p = target[m, j], pred[m, j]
        r = {"n": int(m.sum()), "spearman": float(spearmanr(y, p).statistic)}
        lo, hi = np.quantile(y, [0.2, 0.8])
        yy = np.where(y >= hi, 1, np.where(y <= lo, 0, -1)); kk = yy >= 0
        if kk.sum() >= 20 and len(set(yy[kk])) == 2:
            r["auroc"] = float(roc_auc_score(yy[kk], p[kk]))
        rows[dr] = r
    return rows


def agg(rows, key):
    v = [r[key] for r in rows.values() if key in r and r[key] == r[key]]
    return (round(float(np.mean(v)), 4), len(v)) if v else (None, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--ranks", default="5,12,25,40,60")
    a = ap.parse_args()
    t0 = time.time()
    ranks = [int(x) for x in a.ranks.split(",")]

    A = load_all()
    long, row_of, drugs = A["long"], A["row_of"], A["drugs"]
    specs = sorted(long["specimen"].unique())
    sidx = {s: i for i, s in enumerate(specs)}
    didx = {d: j for j, d in enumerate(drugs)}
    subj = long.drop_duplicates("specimen").set_index("specimen")["subject"].to_dict()
    groups = np.array([subj[s] for s in specs])
    print("matrix: %d specimens x %d inhibitors, %.1f%% observed"
          % (len(specs), len(drugs), 100 * len(long) / (len(specs) * len(drugs))), flush=True)

    AUC = np.full((len(specs), len(drugs)), np.nan)
    for s, dr, v in zip(long["specimen"], long["inhibitor"], long["auc"]):
        if s in sidx and dr in didx:
            AUC[sidx[s], didx[dr]] = v

    P_dep = np.full_like(AUC, np.nan)
    P_mf = {r: np.full_like(AUC, np.nan) for r in ranks}
    Y_cur = np.full_like(AUC, np.nan)          # within-drug z (current target)
    Y_int = np.full_like(AUC, np.nan)          # double-centred (interaction target)

    gk = GroupKFold(n_splits=a.folds)
    for k, (tri, tei) in enumerate(gk.split(np.zeros(len(specs)), groups=groups)):
        med = np.nanmedian(AUC[tri], axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(AUC[tri] - med), axis=0)
        bad = ~np.isfinite(mad) | (mad < 1e-9)
        mad[bad] = np.nanstd(AUC[tri], axis=0)[bad]
        Z = -(AUC - med) / mad
        Y_cur[tei] = Z[tei]
        Y_int[tei] = (Z - np.nanmedian(Z, axis=1, keepdims=True))[tei]

        rows_tr = sorted({row_of[specs[i]] for i in tri})
        fs, Zx = build_space(A["X"], rows_tr, A["meta"], A["mut"], 100, 4000, A["genes"])
        Xp, _, _ = fs.transform(Zx, meta=A["meta"], mut=A["mut"],
                                blocks=["rna", "state", "mut", "clin"])
        Xs = np.vstack([Xp[row_of[specs[i]]] for i in range(len(specs))])
        for r in ranks:
            Utr, V = soft_impute_svd(Z[tri], r)
            reg = RidgeCV(alphas=(10.0, 100.0, 1000.0, 1e4)).fit(Xs[tri], Utr)
            P_mf[r][tei] = reg.predict(Xs[tei]) @ V

        f_tr = long[long["specimen"].isin({specs[i] for i in tri})]
        f_te = long[long["specimen"].isin({specs[i] for i in tei})]
        _, _, lte, p, _ = fit_predict(f_tr, f_te, Zx, row_of, drugs, fs, A["sym2ens"],
                                      ["rna", "state", "mut", "clin"], A["meta"], A["mut"], 0.20)
        for s, dr, v in zip(lte["specimen"], lte["inhibitor"], p.loc[lte.index]):
            if s in sidx and dr in didx and v == v:
                P_dep[sidx[s], didx[dr]] = v
        print("  fold %d/%d (%.0fs)" % (k + 1, a.folds, time.time() - t0), flush=True)

    res = {"generated": time.strftime("%Y-%m-%d %H:%M"), "folds": a.folds, "ranks": ranks,
           "matrix": {"specimens": len(specs), "inhibitors": len(drugs),
                      "observed_fraction": round(len(long) / (len(specs) * len(drugs)), 4)}}

    # ---- C1 rank sweep + blends, on BOTH targets ------------------------------------------
    print("\n== C1 rank sweep ==", flush=True)
    c1 = {}
    dep_rows_i = per_drug(P_dep, Y_int, drugs)
    c1["deployed"] = {"interaction_auroc": agg(dep_rows_i, "auroc")[0],
                      "interaction_spearman": agg(dep_rows_i, "spearman")[0],
                      "current_auroc": agg(per_drug(P_dep, Y_cur, drugs), "auroc")[0]}
    for r in ranks:
        ri = per_drug(P_mf[r], Y_int, drugs)
        bl = np.where(np.isfinite(P_mf[r]) & np.isfinite(P_dep),
                      0.5 * np.nan_to_num(P_mf[r]) + 0.5 * np.nan_to_num(P_dep), np.nan)
        bi = per_drug(bl, Y_int, drugs)
        c1["mf_rank_%d" % r] = {"interaction_auroc": agg(ri, "auroc")[0],
                                "interaction_spearman": agg(ri, "spearman")[0]}
        c1["blend50_rank_%d" % r] = {"interaction_auroc": agg(bi, "auroc")[0],
                                     "interaction_spearman": agg(bi, "spearman")[0],
                                     "current_auroc": agg(per_drug(bl, Y_cur, drugs), "auroc")[0]}
    res["C1_rank_sweep"] = c1
    for k, v in c1.items():
        print("  %-22s interaction AUROC %s" % (k, v.get("interaction_auroc")), flush=True)

    # ---- C4 blend weight sweep ------------------------------------------------------------
    best_r = max(ranks, key=lambda r: (c1["blend50_rank_%d" % r]["interaction_auroc"] or 0))
    c4 = {}
    for w in (0.0, 0.25, 0.4, 0.5, 0.6, 0.75, 1.0):
        bl = np.where(np.isfinite(P_mf[best_r]) & np.isfinite(P_dep),
                      w * np.nan_to_num(P_mf[best_r]) + (1 - w) * np.nan_to_num(P_dep), np.nan)
        c4["w=%.2f" % w] = agg(per_drug(bl, Y_int, drugs), "auroc")[0]
    res["C4_blend_weight"] = {"rank": best_r, "sweep": c4,
                              "note": "w=0 is the deployed model, w=1 is pure matrix factorisation"}
    print("\n== C4 blend weight (rank %d) ==" % best_r, flush=True)
    print("  " + "  ".join("%s:%s" % (k, v) for k, v in c4.items()), flush=True)

    # ---- C2/C3 reliability vs predictability ----------------------------------------------
    print("\n== C2/C3 reliability vs predictability ==", flush=True)
    # prefer the saved artefact: it is what every other document quotes, so a divergence here would be
    # a second set of reliability numbers rather than a check on the first
    tsv = os.path.join(ROOT, "deliverables", "drug_assay_reliability.tsv")
    rel = None
    if os.path.exists(tsv):
        rel = pd.read_csv(tsv, sep="\t").set_index("inhibitor")
    else:
        try:
            rel = DD.drug_reliability(long)
            if not isinstance(rel, pd.DataFrame):
                rel = pd.DataFrame(rel)
            if "inhibitor" in rel.columns:
                rel = rel.set_index("inhibitor")
        except Exception as ex:
            print("  reliability unavailable: %s" % ex, flush=True)
    if rel is not None:
        rows = per_drug(P_dep, Y_int, drugs)
        recs = []
        for dr, r in rows.items():
            if dr not in rel.index:
                continue
            rr = rel.loc[dr]
            reliab = float(rr.get("reliability", np.nan))
            ceil = float(rr.get("ceiling_spearman", np.nan))
            recs.append({"inhibitor": dr, "n": r["n"], "auroc": round(r.get("auroc", np.nan), 4),
                         "spearman": round(r["spearman"], 4), "reliability": reliab,
                         "ceiling": ceil, "tier": str(rr.get("reliability_tier", "")),
                         "frac_of_ceiling": (round(r["spearman"] / ceil, 3)
                                             if ceil == ceil and ceil > 0 else None)})
        recs = [x for x in recs if x["reliability"] == x["reliability"]]
        if len(recs) >= 8:
            sp = np.array([x["spearman"] for x in recs]); rl = np.array([x["reliability"] for x in recs])
            st = spearmanr(rl, sp)
            res["C2_reliability_vs_predictability"] = {
                "n_drugs": len(recs), "spearman": round(float(st.statistic), 4),
                "p": float(st.pvalue),
                "interpretation": ("positive = drugs whose assay reproduces better are also the drugs "
                                   "the model predicts better, which is what the ceiling argument "
                                   "predicts")}
            print("  reliability vs predictability: rho=%.3f, p=%.3g over %d drugs"
                  % (st.statistic, st.pvalue, len(recs)), flush=True)
            by = {}
            for x in recs:
                by.setdefault(x["tier"], []).append(x["spearman"])
            res["C2_by_tier"] = {k: {"n": len(v), "mean_spearman": round(float(np.mean(v)), 4)}
                                 for k, v in sorted(by.items())}
        res["C3_per_drug"] = sorted(recs, key=lambda x: -(x["spearman"]))
        pd.DataFrame(res.get("C3_per_drug") or []).to_csv(
            os.path.join(ROOT, "deliverables", "compass_per_drug_deep.tsv"), sep="\t", index=False)

    # ---- C5 n vs performance ---------------------------------------------------------------
    rows = per_drug(P_dep, Y_int, drugs)
    ns = np.array([r["n"] for r in rows.values()], float)
    au = np.array([r.get("auroc", np.nan) for r in rows.values()], float)
    m = np.isfinite(ns) & np.isfinite(au)
    if m.sum() >= 8:
        st = spearmanr(ns[m], au[m])
        res["C5_n_vs_auroc"] = {"n_drugs": int(m.sum()), "spearman": round(float(st.statistic), 4),
                                "p": float(st.pvalue),
                                "interpretation": ("near zero = the weak drugs are hard, not merely "
                                                   "under-trained; positive = more data would help")}
        print("\n  n_train vs AUROC: rho=%.3f, p=%.3g" % (st.statistic, st.pvalue), flush=True)

    json.dump(res, open(OUT, "w"), indent=1)
    print("\nwrote %s (%.0fs)" % (OUT, time.time() - t0))


if __name__ == "__main__":
    main()
