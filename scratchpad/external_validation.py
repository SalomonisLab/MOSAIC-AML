#!/usr/bin/env python3
"""Expand external validation of the bulk mutation caller beyond the 26-sample sc held-out.

The sealed sc held-out is tiny (26 samples, 1-5 positives per driver), so per-driver AUROCs there are
noisy. But the bulk caller was trained on BeatAML ONLY -- which makes THREE cohorts fully external to it:

    Leucegene  367 samples   never trained on, dbGAP variant calls
    sc atlas   387 samples   never trained on (the bulk caller never saw sc at all -- not just the 26)
    Trumpp      16 samples   never trained on, Table S4 known drivers

AUROC is invariant to the percentile transform (it's monotonic in the decision score), so the choice of
score reference does not affect these numbers -- only calls/thresholds would be affected.

Trumpp truth is GENE-level (Table S4), while the caller is variant-level, so Trumpp is scored gene-level:
score = max prob over that gene's categories, y = gene in known_drivers. Leucegene/sc are scored at the
category level directly.
"""
import os, json, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from amlmm.bulk_predictor import BulkMutationPredictor

BP = BulkMutationPredictor.load(os.path.join(ROOT, "pipeline", "bulk_mutation_predictor.pkl"))
d = np.load(os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz"), allow_pickle=True)
genes = [str(g) for g in d["genes"]]
cats = [str(c) for c in d["drivers"]]
hold = set(str(x) for x in d["holdout"])


def score_cohort(Xkey, Lkey, skey, ref, restrict=None):
    """-> {category: (auroc, n_pos, n)} for every category this cohort actually labels."""
    X = pd.DataFrame(d[Xkey].astype(float), index=[str(s) for s in d[skey]], columns=genes)
    L = pd.DataFrame(d[Lkey].astype(float), index=X.index, columns=cats)
    if restrict is not None:
        keep = [s for s in X.index if s in restrict]
        X, L = X.loc[keep], L.loc[keep]
    Z = {s: BP._z(BP._clog(BP._align(X.loc[s])), ref) for s in X.index}
    out = {}
    for c in BP.categories:
        if c not in L.columns:
            continue
        y = L[c].dropna()
        if int((y == 1).sum()) < 3 or int((y == 0).sum()) < 3:
            continue
        ids = list(y.index)
        p = np.array([BP.predict_one(c, Z[s], ref)["probability"] for s in ids])
        yv = y.loc[ids].astype(int).values
        try:
            out[c] = (roc_auc_score(yv, p), int(yv.sum()), len(yv))
        except Exception:
            pass
    return out


print("scoring cohorts (bulk caller trained on BeatAML only -> all of these are external)...")
lg = score_cohort("lg_X", "lg_L", "lg_samples", "leucegene")
sc_all = score_cohort("sc_X", "sc_L", "sc_samples", "sc")
sc_ho = score_cohort("sc_X", "sc_L", "sc_samples", "sc", restrict=hold)
print("  leucegene: %d categories | sc(all 387): %d | sc(held-out 26): %d" % (len(lg), len(sc_all), len(sc_ho)))

# ---- Trumpp: gene-level (Table S4) ----
gene_of = lambda c: str(c).split("_")[0].split("-")[0].upper()
rows = []
for f in sorted(glob.glob(os.path.join(ROOT, "runs", "trumpp_*", "patient_report.json"))):
    r = json.load(open(f))
    kn = set(str(x).upper() for x in (r.get("known_drivers") or []))
    best = {}
    for p in (r.get("mutation_predictions") or []):
        g = gene_of(p["mutation"])
        pr = p.get("probability")
        if pr is None:
            continue
        best[g] = max(best.get(g, 0.0), float(pr))
    rows.append({"sample": r.get("sample_key"), "known": kn, "best": best})
tr = {}
allg = sorted(set().union(*[set(x["best"]) for x in rows])) if rows else []
for g in allg:
    y = [1 if g in x["known"] else 0 for x in rows if g in x["best"]]
    p = [x["best"][g] for x in rows if g in x["best"]]
    if sum(y) >= 3 and (len(y) - sum(y)) >= 3:
        try:
            tr[g] = (roc_auc_score(y, p), sum(y), len(y))
        except Exception:
            pass
print("  trumpp (gene-level, 16 samples): %d genes with >=3 pos" % len(tr))

# ---- BeatAML CV for reference ----
ba_cv = {c: BP.models[c]["cv_auroc"] for c in BP.categories}

print()
print("=" * 96)
print("%-24s %14s %16s %14s %12s" % ("category", "BeatAML CV", "Leucegene(ext)", "sc all(ext)", "sc held-out"))
print("=" * 96)
keys = sorted(set(lg) | set(sc_all), key=lambda c: -(ba_cv.get(c) or 0))
def fmt(d_, c):
    if c not in d_: return "        -"
    a, npos, n = d_[c]
    return "%.3f (%d/%d)" % (a, npos, n)
for c in keys:
    print("%-24s %14s %16s %14s %12s" % (c, ("%.3f" % ba_cv[c]) if ba_cv.get(c) else "-",
                                          fmt(lg, c), fmt(sc_all, c), fmt(sc_ho, c)))

print()
print("=" * 60)
print("SUMMARY - external validation size and mean AUROC")
print("=" * 60)
for name, dd, n in [("Leucegene", lg, 367), ("sc atlas (all)", sc_all, 387),
                    ("sc held-out only", sc_ho, 26), ("Trumpp (gene-level)", tr, 16)]:
    if dd:
        aucs = [v[0] for v in dd.values()]
        pos = sum(v[1] for v in dd.values())
        print("  %-20s %2d categories | %4d total positives | mean AUROC %.3f"
              % (name, len(dd), pos, float(np.mean(aucs))))
print()
print("  Trumpp per-gene:", ", ".join("%s %.2f(%d/%d)" % (g, v[0], v[1], v[2]) for g, v in sorted(tr.items())))

# ---- the payoff: how much does adding Leucegene+Trumpp grow the evidence base? ----
old_pos = sum(v[1] for v in sc_ho.values())
new_pos = old_pos + sum(v[1] for v in lg.values()) + sum(v[1] for v in tr.values())
print()
print("  EVIDENCE BASE: sc held-out alone = %d positives across %d categories" % (old_pos, len(sc_ho)))
print("                 + Leucegene + Trumpp = %d positives across %d categories"
      % (new_pos, len(set(sc_ho) | set(lg) | set(tr))))
