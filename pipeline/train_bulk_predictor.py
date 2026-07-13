#!/usr/bin/env python3
"""Train the deployable PRIMARY bulk mutation caller (variant-level, BeatAML-trained) and validate it on
our sealed single-cell held-out via each sample's bulk-equivalent.  -> bulk_mutation_predictor.pkl + card."""
import os, sys, json, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from amlmm.bulk_predictor import BulkMutationPredictor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BUNDLE = os.path.join(os.path.dirname(ROOT), "aml-bakeoff", "bundle_data.npz")
BA = os.path.join(ROOT, "data", "external", "beataml")

# gene symbol -> ENSG from BeatAML norm_exp header (cheap: two columns)
hdr = pd.read_csv(os.path.join(BA, "norm_exp.txt"), sep="\t", usecols=["stable_id", "display_label"])
sym2ens = {}
for e, s in zip(hdr["stable_id"].astype(str), hdr["display_label"].astype(str)):
    sym2ens.setdefault(s, e)

print("training bulk mutation predictor from", os.path.basename(BUNDLE))
P = BulkMutationPredictor.train_from_bundle(BUNDLE, sym2ens, varcap=2500, min_pos=6, base="logL2")
P.save(os.path.join(HERE, "bulk_mutation_predictor.pkl"))
print("trained %d categories | mean CV AUROC %.3f" % (len(P.categories), P.summary()["mean_cv_auroc"]))

# ---- validate on the sealed single-cell held-out (bulk-equivalent via ref='sc') ----
d = np.load(BUNDLE, allow_pickle=True)
genes = [str(g) for g in d["genes"]]; cats = [str(c) for c in d["drivers"]]
scX = pd.DataFrame(d["sc_X"].astype(float), index=[str(s) for s in d["sc_samples"]], columns=genes)
scL = pd.DataFrame(d["sc_L"].astype(float), index=scX.index, columns=cats)
hold = [s for s in scX.index if s in set(str(x) for x in d["holdout"])]

card = {"trained": P.meta, "summary": P.summary(), "categories": {}}
val_rows = []
for cat in P.categories:
    m = P.models[cat]
    entry = {"cv_auroc": m["cv_auroc"], "n_pos_train": m["n_pos"], "threshold": m["threshold"]}
    # sc held-out AUROC where we have sc labels
    yl = scL.loc[hold, cat] if cat in scL.columns else pd.Series(dtype=float)
    yl = yl.dropna()
    if yl.notna().sum() >= 4 and yl.nunique() == 2:
        probs = [P.predict_one(cat, P._z(P._clog(P._align(scX.loc[s])), "sc"))["probability"] for s in yl.index]
        try:
            ho = roc_auc_score(yl.astype(int).values, probs)
            entry["sc_heldout_auroc"] = round(float(ho), 3); entry["sc_heldout_npos"] = int(yl.sum())
            val_rows.append((cat, m["cv_auroc"], round(float(ho), 3), int(yl.sum())))
        except Exception:
            pass
    card["categories"][cat] = entry
json.dump(card, open(os.path.join(HERE, "bulk_model_card.json"), "w"), indent=1)

print("\n== validation on sealed single-cell held-out (bulk-equivalent, ref='sc') ==")
print("%-24s %8s %10s %6s" % ("category", "BA_CV", "sc_heldout", "n+"))
for cat, cv, ho, n in sorted(val_rows, key=lambda x: -x[2]):
    print("%-24s %8s %10s %6d" % (cat, cv, ho, n))
print("\ncategories with a bulk-CV AUROC (all): %d | mean %.3f" %
      (len(P.categories), np.mean([P.models[c]["cv_auroc"] for c in P.categories if P.models[c]["cv_auroc"]])))
print("wrote bulk_mutation_predictor.pkl + bulk_model_card.json")

# quick smoke: predict one held-out sample end-to-end (symbol- and ENSG-robust)
s0 = hold[0]
calls = P.predict(scX.loc[s0], ref="sc")
top = list(calls.items())[:6]
print("\nsmoke predict(%s): top calls =" % s0, [(c, v["probability"], v["call"]) for c, v in top])
print("BULK PREDICTOR TRAIN OK")
