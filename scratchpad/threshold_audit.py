#!/usr/bin/env python3
"""Is the deployed (BeatAML F1-max) threshold systematically wrong when applied to sc samples?

The thresholds are F1-max'd on BeatAML CV-OOF percentiles. At predict time an sc sample is z-scored
against the *sc* reference but percentile-mapped against the *BeatAML* training score distribution.
If the sc scores sit systematically lower/higher, every threshold mis-fires the same way.

Score all 387 sc bulk-equivalents, then for every category where the sc cohort actually has labels:
  * does the ranking transfer at all (sc AUROC)?
  * at the DEPLOYED threshold: sensitivity / precision / F1 on sc
  * what threshold would be F1-optimal ON SC, and how far is that from the deployed one?
  * prevalence in BeatAML (which set the threshold) vs prevalence in sc (where it's applied)
"""
import os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from amlmm.bulk_predictor import BulkMutationPredictor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BP = BulkMutationPredictor.load(os.path.join(ROOT, "pipeline", "bulk_mutation_predictor.pkl"))
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
genes = [str(g) for g in d["genes"]]; cats = [str(c) for c in d["drivers"]]
scX = pd.DataFrame(d["sc_X"].astype(float), index=[str(s) for s in d["sc_samples"]], columns=genes)
scL = pd.DataFrame(d["sc_L"].astype(float), index=scX.index, columns=cats)
baL = pd.DataFrame(d["ba_L"].astype(float), columns=cats)

# score every sc sample once (percentile per category, exactly as deployed)
Z = {s: BP._z(BP._clog(BP._align(scX.loc[s])), "sc") for s in scX.index}
print("scored %d sc samples\n" % len(Z))


def f1_at(y, p, t):
    pr = p >= t
    tp = int((pr & (y == 1)).sum()); fp = int((pr & (y == 0)).sum()); fn = int(((~pr) & (y == 1)).sum())
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    return (2 * prec * rec / (prec + rec) if (prec + rec) else 0.0), prec, rec, tp, fp, fn


rows = []
for cat in BP.categories:
    if cat not in scL.columns:
        continue
    y = scL[cat].dropna()
    if int((y == 1).sum()) < 5 or int((y == 0).sum()) < 5:
        continue
    ids = list(y.index)
    p = np.array([BP.predict_one(cat, Z[s])["probability"] for s in ids])
    yv = y.loc[ids].astype(int).values
    thr = BP.models[cat]["threshold"]
    try:
        auc = roc_auc_score(yv, p)
    except Exception:
        continue
    f1_dep, prec_d, rec_d, tp, fp, fn = f1_at(yv, p, thr)
    # best achievable threshold ON SC
    best_t, best_f1 = thr, f1_dep
    for t in np.unique(p):
        f1, _, _, _, _, _ = f1_at(yv, p, t)
        if f1 > best_f1:
            best_f1, best_t = f1, float(t)
    ba_prev = float((baL[cat] == 1).sum() / max((baL[cat].notna()).sum(), 1))
    sc_prev = float((yv == 1).mean())
    rows.append(dict(drv=cat, n_pos=int(yv.sum()), n=len(yv), auc=auc, thr=thr, sc_best_thr=best_t,
                     delta=thr - best_t, f1_dep=f1_dep, f1_best=best_f1, sens=rec_d, prec=prec_d,
                     tp=tp, fp=fp, fn=fn, ba_prev=ba_prev, sc_prev=sc_prev,
                     med_pct=float(np.median(p))))

df = pd.DataFrame(rows).sort_values("auc", ascending=False)
print("=" * 118)
print("%-12s %4s %5s %6s | %6s %8s %7s | %6s %6s | %6s %6s | %5s %5s %5s | %6s %6s" %
      ("category", "n+", "n", "scAUC", "thr", "scBestT", "delta", "F1@dep", "F1best", "sens", "prec", "TP", "FP", "FN", "BAprev", "SCprev"))
print("=" * 118)
for _, r in df.iterrows():
    print("%-12s %4d %5d %6.3f | %6.3f %8.3f %+7.3f | %6.3f %6.3f | %6.2f %6.2f | %5d %5d %5d | %6.3f %6.3f" %
          (r["drv"], r["n_pos"], r["n"], r["auc"], r["thr"], r["sc_best_thr"], r["delta"], r["f1_dep"], r["f1_best"],
           r["sens"], r["prec"], r["tp"], r["fp"], r["fn"], r["ba_prev"], r["sc_prev"]))

print("\n" + "=" * 60)
print("SYSTEMATIC CHECK")
print("=" * 60)
print("median(deployed_thr - sc_optimal_thr) = %+.3f" % df["delta"].median())
print("  categories where deployed is TOO TIGHT (delta>0, under-calls): %d/%d" % (int((df["delta"] > 0).sum()), len(df)))
print("  categories where deployed is TOO LOOSE (delta<0, over-calls) : %d/%d" % (int((df["delta"] < 0).sum()), len(df)))
print("mean sensitivity at the deployed threshold: %.2f   (mean precision %.2f)" % (df["sens"].mean(), df["prec"].mean()))
print("mean F1 @deployed %.3f  vs  mean F1 @sc-optimal %.3f   (headroom %+.3f)"
      % (df["f1_dep"].mean(), df["f1_best"].mean(), df["f1_best"].mean() - df["f1_dep"].mean()))
print("\nprevalence shift (threshold is set by BAprev, applied where SCprev holds):")
print("  mean BeatAML prevalence %.3f  vs  mean sc prevalence %.3f" % (df["ba_prev"].mean(), df["sc_prev"].mean()))
print("  corr(threshold, BeatAML prevalence) = %.2f   [a calibrated pctile thr should track 1-prev]"
      % df[["thr", "ba_prev"]].corr().iloc[0, 1])
print("\nmedian percentile assigned to sc samples (should be ~0.5 if the score distribution matches BeatAML):")
print("  %.3f" % df["med_pct"].median())
