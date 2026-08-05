#!/usr/bin/env python3
"""BULK RNA-ONLY model (BeatAML2-trained) — sensitivity + specificity per mutation, per cohort.

One model under test: pipeline/bulk_mutation_predictor.pkl (BeatAML2, n=707, bulk RNA only, NO imputation).
Evaluated on four cohorts the boss asked for:
  BeatAML   (n=707)  its OWN training cohort -> honest 5-fold CV-OOF (re-derived here, not in-sample)
  Leucegene (n=367)  external bulk held-out
  held-out scRNA (n=29 sealed)   single-cell, bulk-equivalent input
  all scRNA (n=387)  single-cell, bulk-equivalent input (entirely external to a bulk caller)

Reports per category: n, n_pos, prevalence, sensitivity(=recall), specificity, precision, F1, AUROC,
and the raw TP/FP/FN/TN. Plus a cohort OVERALL (pooled + mean). Writes deliverables/bulk_matrix.json + .tsv.

  bsub -q test -W 60 -M 24000 -R "rusage[mem=24000]" -o bm.log \
    /usr/local/anaconda3-2020/bin/python bench_matrix.py
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from amlmm.bulk_predictor import BulkMutationPredictor, _pct

OUTD = os.path.join(ROOT, "deliverables"); os.makedirs(OUTD, exist_ok=True)
BP = BulkMutationPredictor.load(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
GENES = [str(g) for g in d["genes"]]; CATS = [str(c) for c in d["drivers"]]
HOLD = set(str(x) for x in d["holdout"]); MINP = 3
VARCAP = 2500


def confusion(y, call):
    y = np.asarray(y).astype(int); call = np.asarray(call).astype(int)
    tp = int((call & (y == 1)).sum()); fp = int((call & (y == 0)).sum())
    fn = int(((1 - call) & (y == 1)).sum()); tn = int(((1 - call) & (y == 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    f1 = 2 * prec * sens / (prec + sens) if (prec + sens) else 0.0
    return dict(tp=tp, fp=fp, fn=fn, tn=tn, sensitivity=round(sens, 4), specificity=round(spec, 4),
                precision=round(prec, 4), f1=round(f1, 4))


def full(y, p, thr):
    y = np.asarray(y).astype(int); p = np.asarray(p, float)
    m = confusion(y, (p >= thr).astype(int))
    try:
        m["auroc"] = round(float(roc_auc_score(y, p)), 4)
    except Exception:
        m["auroc"] = None
    m.update(n=int(len(y)), n_pos=int(y.sum()), prevalence=round(float(y.mean()), 4), threshold=round(float(thr), 4))
    return m


def external(Xk, Lk, Sk, ref, restrict=None):
    """trained bulk model applied to an external cohort (Leucegene / scRNA)."""
    X = pd.DataFrame(d[Xk].astype(float), index=[str(s) for s in d[Sk]], columns=GENES)
    L = pd.DataFrame(d[Lk].astype(float), index=X.index, columns=CATS)
    if restrict is not None:
        k = [s for s in X.index if s in restrict]; X, L = X.loc[k], L.loc[k]
    Z = {s: BP._z(BP._clog(BP._align(X.loc[s])), ref) for s in X.index}
    res = {}
    for c in BP.categories:
        if c not in L.columns:
            continue
        y = L[c].dropna()
        if int((y == 1).sum()) < MINP or int((y == 0).sum()) < MINP:
            continue
        ids = list(y.index)
        p = [BP.predict_one(c, Z[s], ref)["probability"] for s in ids]
        res[c] = full(y.loc[ids].values, p, BP.models[c]["threshold"])
    return res


def beataml_cv():
    """honest 5-fold CV-OOF on BeatAML — the model's OWN cohort, so in-sample would be optimistic.
    Replicates train_from_bundle's recipe: clog -> z(BeatAML) -> topvar(2500) -> logL2 -> OOF -> F1-max thr."""
    baX = d["ba_X"].astype(float); baL = d["ba_L"].astype(float)
    cl = np.log2(np.clip(baX, 0, None) + 1.0); mu = cl.mean(0); sd = cl.std(0); sd[sd == 0] = 1.0
    Zba = (cl - mu) / sd
    res = {}
    for j, cat in enumerate(CATS):
        if cat not in BP.categories:
            continue
        y = baL[:, j]; ok = ~np.isnan(y); yv = y[ok].astype(int); Z = Zba[ok]
        if yv.sum() < MINP or (yv == 0).sum() < MINP:
            continue
        oof = np.full(len(yv), np.nan)
        ns = min(5, int(yv.sum()), int((yv == 0).sum()))
        for tri, tei in StratifiedKFold(ns, shuffle=True, random_state=0).split(Z, yv):
            sel = np.argsort(Z[tri].var(0))[::-1][:VARCAP]
            est = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000).fit(Z[tri][:, sel], yv[tri])
            oof[tei] = est.decision_function(Z[tei][:, sel])
        om = ~np.isnan(oof)
        try:
            a = roc_auc_score(yv[om], oof[om])
        except Exception:
            a = np.nan
        sign = 1.0 if (a != a or a >= 0.5) else -1.0
        p = _pct(np.sort(oof[om]), oof[om]); p = p if sign > 0 else 1 - p
        yv_oof = yv[om]
        # F1-max threshold on the OOF (same rule the deployed model used)
        bt, bf1 = 0.5, -1.0
        for t in np.unique(p):
            m = confusion(yv_oof, (p >= t).astype(int))
            if m["f1"] > bf1:
                bf1, bt = m["f1"], float(t)
        res[cat] = full(yv_oof, p, bt)
    return res


def summarize(res):
    ok = {k: v for k, v in res.items() if "tp" in v}
    if not ok:
        return {}
    TP = sum(v["tp"] for v in ok.values()); FP = sum(v["fp"] for v in ok.values())
    FN = sum(v["fn"] for v in ok.values()); TN = sum(v["tn"] for v in ok.values())
    psens = TP / (TP + FN) if (TP + FN) else 0.0; pspec = TN / (TN + FP) if (TN + FP) else 0.0
    return dict(n_categories=len(ok),
                mean_sensitivity=round(float(np.mean([v["sensitivity"] for v in ok.values()])), 4),
                mean_specificity=round(float(np.mean([v["specificity"] for v in ok.values()])), 4),
                mean_f1=round(float(np.mean([v["f1"] for v in ok.values()])), 4),
                mean_auroc=round(float(np.mean([v["auroc"] for v in ok.values() if v["auroc"] is not None])), 4),
                pooled_sensitivity=round(psens, 4), pooled_specificity=round(pspec, 4),
                pooled_TP=TP, pooled_FP=FP, pooled_FN=FN, pooled_TN=TN)


COH = {"BeatAML_CV": beataml_cv(),
       "Leucegene": external("lg_X", "lg_L", "lg_samples", "leucegene"),
       "heldout_scRNA": external("sc_X", "sc_L", "sc_samples", "sc", restrict=HOLD),
       "all_scRNA": external("sc_X", "sc_L", "sc_samples", "sc")}

out = {"model": "bulk_mutation_predictor.pkl (BeatAML2-trained, bulk RNA only, no imputation)",
       "mean_beataml_cv_auroc_stored": BP.summary().get("mean_cv_auroc"),
       "cohorts": {k: {"per_category": v, "overall": summarize(v)} for k, v in COH.items()}}
json.dump(out, open(os.path.join(OUTD, "bulk_matrix.json"), "w"), indent=1)

# TSV: one row per (cohort, category)
with open(os.path.join(OUTD, "bulk_matrix.tsv"), "w") as fh:
    fh.write("cohort\tcategory\tn\tn_pos\tprevalence\tsensitivity\tspecificity\tprecision\tf1\tauroc\tTP\tFP\tFN\tTN\n")
    for coh, res in COH.items():
        for c, m in sorted(res.items()):
            fh.write("\t".join(str(x) for x in [coh, c, m["n"], m["n_pos"], m["prevalence"], m["sensitivity"],
                     m["specificity"], m["precision"], m["f1"], m["auroc"], m["tp"], m["fp"], m["fn"], m["tn"]]) + "\n")
    fh.write("\nCOHORT OVERALL (mean over categories | pooled over all calls)\n")
    fh.write("cohort\tn_cat\tmean_sens\tmean_spec\tmean_f1\tmean_auroc\tpooled_sens\tpooled_spec\n")
    for coh, res in COH.items():
        s = summarize(res)
        if s:
            fh.write("%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (coh, s["n_categories"],
                     s["mean_sensitivity"], s["mean_specificity"], s["mean_f1"], s["mean_auroc"],
                     s["pooled_sensitivity"], s["pooled_specificity"]))

print("=== BULK RNA-only model (BeatAML-trained) — overall by cohort ===")
for coh, res in COH.items():
    s = summarize(res)
    print("  %-16s  cats=%2d  mean sens=%.3f spec=%.3f F1=%.3f AUROC=%.3f | pooled sens=%.3f spec=%.3f"
          % (coh, s["n_categories"], s["mean_sensitivity"], s["mean_specificity"], s["mean_f1"],
             s["mean_auroc"], s["pooled_sensitivity"], s["pooled_specificity"]))
print("wrote deliverables/bulk_matrix.json + .tsv")
print("BULK MATRIX OK")
